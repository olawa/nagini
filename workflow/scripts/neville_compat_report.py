#!/usr/bin/env python3
"""Report amplicon variants in the neville_mx_amplicon layout.

nagini calls this panel with Rindels; neville_mx_amplicon calls the same 15
amplicons with a VarDict/ClairS-TO/DeepSomatic ensemble. This module renders
nagini's calls into neville's exact table shape so the two can be compared
sample by sample.

Mirrored from neville workflow/scripts/results_report_create_tables.py:
  * one row per variant, columns in neville's fixed order
  * variants split at 50 bp into SNV_Indels and SVs, SVs gaining "SV length"
  * AF reported as a FRACTION (0-1), as neville does, never a percentage
  * min AF 0.005, matching neville's VariantTable default
  * FILTER shown as a flag in its own column; nothing is dropped for being
    soft-filtered, again as neville does

Emits TSV (the comparable artefact) and optionally XLSX (the familiar shape).
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:
    import xlsxwriter
except ImportError:  # pragma: no cover - only needed for --output-xlsx
    xlsxwriter = None

# neville: VariantTable(min_sv_length=50, min_af=0.005)
MIN_AF = 0.005
MIN_SV_LENGTH = 50

# neville results_report_create_tables.py, create_snv_table()
SNV_COLUMNS = [
    "Filter", "RunID", "DNAnr", "Gene", "Chr", "Pos", "Ref", "Alt", "AF", "DP",
    "Transcript", "Mutation cds", "ENSP", "Consequence", "COSMIC ids on pos",
    "Clinical Significance", "dbSNP", "Max Pop AF", "Max Pop", "Callers",
]

# neville create_sv_table(): identical, with "SV length" inserted after Alt
SV_COLUMNS = SNV_COLUMNS[:8] + ["SV length"] + SNV_COLUMNS[8:]

CSQ_FORMAT_RE = re.compile(r"Format:\s*([^\"']+)")


def _open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") \
        else open(path, "r", encoding="utf-8", errors="replace")


def _to_int(value: Any, default: int = 0) -> int:
    try:
        return int(float(str(value).strip()))
    except (ValueError, TypeError, AttributeError):
        return default


def _to_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(str(value).strip())
    except (ValueError, TypeError, AttributeError):
        return default


def _split_field(value: str, index: int) -> str:
    """Element ``index`` of a comma-separated per-allele value, else the value."""
    if "," not in value:
        return value
    parts = value.split(",")
    return parts[index] if 0 <= index < len(parts) else ""


def parse_vcf_header(path: str) -> Tuple[List[str], str]:
    """Return (CSQ field names, the ##VEP header value).

    neville reads the CSQ layout out of the header the same way (index_vep);
    doing it here keeps the column meanings tied to the file rather than to an
    assumed VEP version.
    """
    csq_fields: List[str] = []
    vep_line = ""
    with _open_text(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            if line.startswith("##INFO=<ID=CSQ"):
                match = CSQ_FORMAT_RE.search(line)
                if match:
                    csq_fields = [f.strip() for f in match.group(1).strip().strip('">').split("|")]
            elif line.startswith("##VEP="):
                vep_line = line[len("##VEP="):].strip()
    return csq_fields, vep_line


def _csq_value(csq_parts: Sequence[str], csq_fields: Sequence[str], name: str) -> str:
    if name not in csq_fields:
        return ""
    idx = csq_fields.index(name)
    return csq_parts[idx] if idx < len(csq_parts) else ""


def _existing_variation(value: str) -> Tuple[str, str]:
    """Split VEP Existing_variation into (COSMIC ids, dbSNP ids), as neville does."""
    existing = value.split("&") if value else []
    cosmic = [c for c in existing if c.startswith("CO")]
    rs = [r for r in existing if r.startswith("rs")]
    return ", ".join(cosmic), ", ".join(rs)


def build_neville_rows(
    vcf_path: str,
    run_id: str = "",
    sample: str = "",
    callers: str = "rindels",
    min_af: float = MIN_AF,
    min_sv_length: int = MIN_SV_LENGTH,
) -> Tuple[List[List[Any]], List[List[Any]]]:
    """Render a VEP-annotated VCF into neville's SNV_Indels and SVs rows."""
    if not vcf_path or not os.path.isfile(vcf_path):
        return [], []

    csq_fields, _ = parse_vcf_header(vcf_path)
    snv_rows: List[List[Any]] = []
    sv_rows: List[List[Any]] = []

    with _open_text(vcf_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _vid, ref, alt, _qual, filt, info = fields[:8]
            fmt = fields[8] if len(fields) > 8 else ""
            sample_data = fields[9] if len(fields) > 9 else ""

            fmt_map = dict(zip(fmt.split(":"), sample_data.split(":"))) if fmt and sample_data else {}
            info_map: Dict[str, str] = {}
            for item in info.split(";"):
                if "=" in item:
                    k, v = item.split("=", 1)
                    info_map[k] = v
                elif item:
                    info_map[item] = "True"

            dp = _to_int(fmt_map.get("DP"), 0) or _to_int(info_map.get("DP"), 0)
            ad_values = [p for p in fmt_map.get("AD", "").split(",") if p != ""]
            csq_entries = info_map.get("CSQ", "").split(",") if info_map.get("CSQ") else []
            record_callers = ";".join(info_map["CALLERS"].split(",")) if info_map.get("CALLERS") else callers

            # One row per ALT. With upstream bcftools norm -m -any this is a
            # single allele, matching neville's decomposed input.
            for alt_index, alt_allele in enumerate(alt.split(",")):
                if not alt_allele or alt_allele == ".":
                    continue

                ad_alt = _to_int(ad_values[alt_index + 1], 0) if len(ad_values) > alt_index + 1 else 0

                af_raw = _split_field(fmt_map.get("AF", fmt_map.get("VAF", "")), alt_index)
                if af_raw in ("", "."):
                    af_raw = _split_field(info_map.get("AF", ""), alt_index)
                if af_raw not in ("", "."):
                    af = _to_float(af_raw)
                    if af > 1.0:  # a caller reporting percent
                        af /= 100.0
                elif dp > 0:
                    af = ad_alt / dp
                else:
                    af = 0.0

                if af <= min_af:
                    continue

                csq_raw = csq_entries[alt_index] if alt_index < len(csq_entries) else (
                    csq_entries[0] if csq_entries else ""
                )
                csq_parts = csq_raw.split("|") if csq_raw else []

                hgvsc = _csq_value(csq_parts, csq_fields, "HGVSc")
                transcript = hgvsc.split(":")[0] if hgvsc else ""
                coding_name = hgvsc.split(":")[1] if len(hgvsc.split(":")) > 1 else ""
                cosmic, rs = _existing_variation(_csq_value(csq_parts, csq_fields, "Existing_variation"))

                var_length = max(len(ref), len(alt_allele))
                head = [
                    filt,
                    run_id,
                    sample,
                    _csq_value(csq_parts, csq_fields, "SYMBOL"),
                    chrom,
                    _to_int(pos),
                    ref,
                    alt_allele,
                ]
                tail = [
                    af,
                    dp,
                    transcript,
                    coding_name,
                    _csq_value(csq_parts, csq_fields, "HGVSp"),
                    _csq_value(csq_parts, csq_fields, "Consequence"),
                    cosmic,
                    _csq_value(csq_parts, csq_fields, "CLIN_SIG"),
                    rs,
                    _csq_value(csq_parts, csq_fields, "MAX_AF"),
                    _csq_value(csq_parts, csq_fields, "MAX_AF_POPS"),
                    record_callers,
                ]

                if var_length <= min_sv_length:
                    snv_rows.append(head + tail)
                else:
                    sv_rows.append(head + [var_length] + tail)

    return snv_rows, sv_rows


def write_tsv(rows: List[List[Any]], columns: Sequence[str], path: str) -> None:
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join("" if v is None else str(v) for v in row) + "\n")


def write_xlsx(
    path: str,
    snv_rows: List[List[Any]],
    sv_rows: List[List[Any]],
    sample: str = "",
    run_id: str = "",
    reference: str = "",
    vep_line: str = "",
    caller_note: str = "",
) -> None:
    """Write the two variant sheets in neville's layout."""
    if xlsxwriter is None:
        raise RuntimeError("xlsxwriter is required for --output-xlsx")

    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    workbook = xlsxwriter.Workbook(path, {"nan_inf_to_errors": True})
    heading = workbook.add_format({"bold": True, "font_size": 18})
    bold = workbook.add_format({"bold": True})
    orange = workbook.add_format({"bg_color": "#ffd280"})

    def sheet(name, rows, columns, title):
        ws = workbook.add_worksheet(name)
        ws.write("A1", title, heading)
        ws.write("A3", "Sample: " + str(sample))
        ws.write("A4", "Reference used: " + str(reference))
        ws.write("A5", "RunID: " + str(run_id))
        ws.write("A6", "Databases used: " + str(vep_line))
        ws.write("A8", "Caller:", orange)
        ws.write("B8", caller_note, orange)
        ws.write_rich_string(
            "A10",
            "Layout follows neville_mx_amplicon. AF is a ", bold, "fraction", " (0-1). "
            "Variants are split at 50 bp; nothing is removed for being soft-filtered.",
        )
        top = 12
        end_col = len(columns) - 1
        height = max(len(rows), 1)
        ws.add_table(
            top - 1, 0, top - 1 + height, end_col,
            {
                "data": rows if rows else [[""] * len(columns)],
                "columns": [{"header": c} for c in columns],
                "style": "Table Style Light 1",
                "autofilter": True,
            },
        )

    sheet("SNV_Indels", snv_rows, SNV_COLUMNS, "Variants found")
    sheet("SVs", sv_rows, SV_COLUMNS, "Variants found")
    workbook.close()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vcf", required=True, help="VEP-annotated, normalized amplicon VCF")
    parser.add_argument("--sample", default="", help="Sample name (neville column DNAnr)")
    parser.add_argument("--run-id", default="", help="Run identifier (neville column RunID)")
    parser.add_argument("--reference", default="", help="Reference FASTA path, for the sheet header")
    parser.add_argument("--callers", default="rindels", help="Value for the Callers column")
    parser.add_argument("--min-af", type=float, default=MIN_AF)
    parser.add_argument("--min-sv-length", type=int, default=MIN_SV_LENGTH)
    parser.add_argument("--output-snv-tsv", help="Output TSV in neville SNV_Indels column order")
    parser.add_argument("--output-sv-tsv", help="Output TSV in neville SVs column order")
    parser.add_argument("--output-xlsx", help="Output XLSX with neville's two variant sheets")
    args = parser.parse_args()

    if not os.path.isfile(args.vcf):
        sys.stderr.write(f"Error: VCF not found: {args.vcf}\n")
        sys.exit(1)

    _, vep_line = parse_vcf_header(args.vcf)
    snv_rows, sv_rows = build_neville_rows(
        args.vcf, run_id=args.run_id, sample=args.sample, callers=args.callers,
        min_af=args.min_af, min_sv_length=args.min_sv_length,
    )

    if args.output_snv_tsv:
        write_tsv(snv_rows, SNV_COLUMNS, args.output_snv_tsv)
    if args.output_sv_tsv:
        write_tsv(sv_rows, SV_COLUMNS, args.output_sv_tsv)
    if args.output_xlsx:
        write_xlsx(
            args.output_xlsx, snv_rows, sv_rows, sample=args.sample, run_id=args.run_id,
            reference=args.reference, vep_line=vep_line,
            caller_note=f"{args.callers} (nagini). Compare against the neville ensemble "
                        "(clairs_to + deepsomatic + vardict).",
        )


if __name__ == "__main__":
    main()
