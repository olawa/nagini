"""Create the XLSX summary for a nagini analysis.

Sheets are built as plain lists of rows and handed to xlsxwriter, the same
writer the other Clinical Genomics Uppsala pipelines use.  A cell is either a
value or a ``(value, style_id)`` pair; the style vocabulary is defined next to
the writer below.
"""

import argparse
import bisect
import collections
import csv
import gzip
import os
import re
from datetime import datetime, timezone

import xlsxwriter


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r")


def parse_info(value):
    result = {}
    for item in value.split(";"):
        if "=" in item:
            key, item_value = item.split("=", 1)
            result[key] = item_value
        elif item:
            result[item] = True
    return result


class CytobandLookup:
    """Lookup UCSC cytobands from plain or gzip-compressed BED-like files."""

    def __init__(self, path=""):
        self.bands = collections.defaultdict(list)
        if path and os.path.isfile(path):
            self._load(path)

    def _load(self, path):
        with open_text(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    continue
                self.bands[fields[0]].append((end, fields[3], start))
        for chrom in self.bands:
            self.bands[chrom].sort()

    def get_band(self, chrom, position):
        if not self.bands or not chrom or not position:
            return ""
        keys = (chrom, chrom.removeprefix("chr"), f"chr{chrom.removeprefix('chr')}")
        entries = next((self.bands[key] for key in keys if key in self.bands), None)
        if not entries:
            return ""
        position = int(position)
        index = bisect.bisect_left(entries, (position, "", 0))
        if index >= len(entries):
            return ""
        end, band, start = entries[index]
        return band if start <= position <= end else ""


def _first_annotation_value(annotation, *names):
    for name in names:
        value = annotation.get(name, "")
        if value not in ("", "."):
            return value
    return ""


def _first_annotation(raw_value, field_count):
    """Extract one VEP annotation without splitting commas inside a field."""

    if not raw_value:
        return ""
    separators = 0
    for index, character in enumerate(str(raw_value)):
        if character == "|":
            separators += 1
        elif character == "," and separators >= field_count - 1:
            return str(raw_value)[:index]
    return str(raw_value)


def _annotation_fields(header_line, annotation_id):
    if f"ID={annotation_id}" not in header_line:
        return []
    match = re.search(r"(?:Format:|Functional annotations:)[ ']*([^\"'>]+)", header_line)
    return [field.strip() for field in match.group(1).rstrip(" '").split("|")] if match else []


def parse_small_variant_vcf(path, default_caller=""):
    columns = [
        "Caller", "Chromosome", "Position", "Reference", "Alternate", "Filter", "Quality",
        "VAF", "Depth", "Allele depth", "Gene", "Consequence", "Impact",
        "HGVSc", "HGVSp", "Existing variation", "ClinVar_Significance",
        "ClinVar_Disease", "ClinVar_Review",
    ]
    records = []
    csq_fields = []
    ann_fields = []
    detected_caller = default_caller

    path_str = str(path).lower()
    if "clairs" in path_str:
        detected_caller = "ClairS-TO"
    elif "mutect" in path_str:
        detected_caller = "Mutect2"

    if not path or not os.path.exists(path):
        return columns, records, (detected_caller or "ClairS-TO")

    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##"):
                line_lower = line.lower()
                if "clairs-to" in line_lower or "clairs_to" in line_lower or "source=clairs-to" in line_lower:
                    detected_caller = "ClairS-TO"
                elif "gatkcommandline=mutect2" in line_lower or "mutect2" in line_lower:
                    detected_caller = "Mutect2"
                if line.startswith("##INFO=<ID=CSQ"):
                    csq_fields = _annotation_fields(line, "CSQ")
                elif line.startswith("##INFO=<ID=ANN"):
                    ann_fields = _annotation_fields(line, "ANN")
                continue
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _record_id, ref, alt, qual, filter_value, info_text = fields[:8]
            info = parse_info(info_text)
            sample_values = {}
            if len(fields) > 9:
                sample_values = dict(zip(fields[8].split(":"), fields[9].split(":")))

            annotation = {}
            if csq_fields and info.get("CSQ"):
                annotation = dict(
                    zip(csq_fields, _first_annotation(info["CSQ"], len(csq_fields)).split("|"))
                )
            elif ann_fields and info.get("ANN"):
                annotation = dict(
                    zip(ann_fields, _first_annotation(info["ANN"], len(ann_fields)).split("|"))
                )

            clinvar_significance = _first_annotation_value(
                annotation, "ClinVar_CLNSIG", "CLIN_SIG", "CLNSIG"
            ) or str(info.get("CLNSIG", ""))
            clinvar_disease = _first_annotation_value(
                annotation, "ClinVar_CLNDN", "CLNDN"
            ) or str(info.get("CLNDN", ""))
            clinvar_review = _first_annotation_value(
                annotation, "ClinVar_CLNREVSTAT", "CLNREVSTAT"
            ) or str(info.get("CLNREVSTAT", ""))

            raw_vaf = sample_values.get("AF", sample_values.get("VAF", ""))
            try:
                vaf_float = float(str(raw_vaf).split(",")[0])
                vaf_display = f"{vaf_float:.3f}"
            except (ValueError, TypeError, AttributeError):
                vaf_display = raw_vaf

            records.append({
                "Caller": detected_caller or "ClairS-TO",
                "Chromosome": chrom,
                "Position": int(pos),
                "Reference": ref,
                "Alternate": alt,
                "Filter": filter_value,
                "Quality": qual,
                "VAF": vaf_display,
                "Depth": sample_values.get("DP", info.get("DP", "")),
                "Allele depth": sample_values.get("AD", ""),
                "Gene": annotation.get("SYMBOL", annotation.get("Gene_Name", "")),
                "Consequence": annotation.get("Consequence", annotation.get("Annotation", "")),
                "Impact": annotation.get("IMPACT", annotation.get("Annotation_Impact", "")),
                "HGVSc": annotation.get("HGVSc", annotation.get("HGVS.c", "")),
                "HGVSp": annotation.get("HGVSp", annotation.get("HGVS.p", "")),
                "Existing variation": annotation.get("Existing_variation", ""),
                "ClinVar_Significance": clinvar_significance,
                "ClinVar_Disease": clinvar_disease,
                "ClinVar_Review": clinvar_review,
            })
    return columns, records, (detected_caller or "ClairS-TO")


def parse_tsv(path):
    if not path or not os.path.exists(path):
        return [], []
    with open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
        return list(reader.fieldnames or []), rows


def _sv_position(value, default=0):
    try:
        return int(str(value).split(",", 1)[0])
    except (TypeError, ValueError):
        return default


def annotate_sv_rows(columns, rows, cytobands):
    """Add cytoband columns to SV summaries when a UCSC file is configured."""

    if not cytobands.bands or not rows:
        return columns, rows
    columns = list(columns)
    for name in ("Start_Band", "End_Band", "Cytogenetic_Notation"):
        if name not in columns:
            columns.append(name)
    for row in rows:
        chrom = row.get("chrom", row.get("Chr", ""))
        chrom2 = row.get("chrom2", "")
        pos = _sv_position(row.get("pos", row.get("Start", 0)))
        end = _sv_position(row.get("end", row.get("End", pos)), pos)
        band1 = row.get("Start_Band") or row.get("cytoband") or cytobands.get_band(chrom, pos)
        band2 = row.get("End_Band") or row.get("cytoband2") or cytobands.get_band(chrom2 or chrom, end)
        row["Start_Band"] = band1
        row["End_Band"] = band2
        svtype = row.get("svtype", row.get("SV_Type", "SV"))
        if row.get("Cytogenetic_Notation") or row.get("cytogenetic_notation"):
            row["Cytogenetic_Notation"] = row.get("Cytogenetic_Notation") or row.get("cytogenetic_notation")
            continue
        if chrom2 and chrom2 != chrom:
            c1 = str(chrom).removeprefix("chr")
            c2 = str(chrom2).removeprefix("chr")
            row["Cytogenetic_Notation"] = (
                f"t({c1};{c2})({band1};{band2})" if band1 and band2 else f"t({c1};{c2})"
            )
        elif str(svtype).upper() == "DEL" and band1 and band2:
            row["Cytogenetic_Notation"] = (
                f"del({str(chrom).removeprefix('chr')})({band1}{band2})"
            )
        else:
            row["Cytogenetic_Notation"] = ""
    return columns, rows


def _is_pathogenic_clinvar(value):
    normalized = re.sub(r"[\s-]+", "_", str(value or "").lower())
    if "benign" in normalized or "not_pathogenic" in normalized:
        return False
    return bool(re.search(r"(?:^|[|&,;/_])(?:likely_)?pathogenic(?:$|[|&,;/_])", normalized))


def parse_cramino(path):
    rows = []
    if not path or not os.path.exists(path):
        return ["Metric", "Value"], rows
    with open_text(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if "\t" in line:
                key, value = line.split("\t", 1)
            elif ":" in line:
                key, value = line.split(":", 1)
            else:
                continue
            rows.append({"Metric": key.strip(), "Value": value.strip()})
    return ["Metric", "Value"], rows


def _find_metric(data_dict, *keywords):
    for k, v in data_dict.items():
        k_lower = k.lower()
        if any(kw.lower() in k_lower for kw in keywords):
            return str(v).strip()
    return ""


def parse_bcr_abl_tsv(path):
    rows = []
    if not path or not os.path.exists(path):
        return rows
    with open_text(path) as handle:
        header = None
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields)))
    return rows


def _style_tki_response(val):
    norm = str(val or "").strip().lower()
    if norm == "resistant":
        return 4
    if norm in ("intermediate", "reduced"):
        return 6
    if norm == "sensitive":
        return 5
    return 0


def table(columns, records):
    return [columns] + [[record.get(column, "") for column in columns] for record in records]



# ---------------------------------------------------------------------------
# Workbook writing
# ---------------------------------------------------------------------------
#
# Style ids are referenced directly by the sheet builders below, so they are a
# small stable vocabulary rather than an implementation detail:
#
#   0 plain      1 table header   2 title banner   3 section header
#   4 alert      5 good           6 warning        7 field label
#   8 highlight  9 secondary alert               10 bold plain

_BASE = {"font_name": "Calibri", "font_size": 11}
_BORDER = {"border": 1, "border_color": "#D9D9D9"}

STYLE_SPECS = {
    0: {**_BASE, **_BORDER},
    1: {**_BASE, **_BORDER, "bold": True, "font_color": "#FFFFFF", "bg_color": "#1F4E78"},
    2: {**_BASE, "bold": True, "font_color": "#FFFFFF", "bg_color": "#14324F", "font_size": 13},
    3: {**_BASE, **_BORDER, "bold": True, "font_color": "#1F4E78", "bg_color": "#D9E1F2"},
    4: {**_BASE, **_BORDER, "bold": True, "font_color": "#9C0006", "bg_color": "#FFC7CE"},
    5: {**_BASE, **_BORDER, "bold": True, "font_color": "#006100", "bg_color": "#C6EFCE"},
    6: {**_BASE, **_BORDER, "bold": True, "font_color": "#9C6500", "bg_color": "#FFEB9C"},
    7: {**_BASE, **_BORDER, "bold": True, "font_color": "#1F4E78", "bg_color": "#F2F2F2"},
    8: {**_BASE, **_BORDER, "bold": True, "font_color": "#371485", "bg_color": "#EDEDF9"},
    9: {**_BASE, **_BORDER, "bold": True, "font_color": "#9C0006", "bg_color": "#FCE4D6"},
    10: {**_BASE, **_BORDER, "bold": True, "font_color": "#1F4E78"},
}

# Sheets that lay out a report rather than tabulate rows get no autofilter and
# no frozen header; a filter dropdown on a title banner is just noise.
LAYOUT = "layout"
TABLE = "table"

_ILLEGAL_SHEET_CHARS = re.compile(r"[\[\]:*?/\\]")
_CONTROL_CHARS = re.compile(r"[\x00-\x08\x0B\x0C\x0E-\x1F]")


def _cell_value_and_style(cell, row_index, column_index, default_header_style=1):
    if isinstance(cell, tuple) and len(cell) == 2:
        return cell[0], cell[1]
    if row_index == 1:
        return cell, default_header_style
    val_str = str(cell).strip()
    if val_str in ("PASS", "OK"):
        return cell, 5
    if _is_pathogenic_clinvar(val_str):
        return cell, 4
    if val_str.upper() in ("FAIL", "HIGH", "ACTION REQUIRED"):
        return cell, 4
    if val_str.upper() in ("BORDERLINE", "REVIEW", "WARNING"):
        return cell, 6
    return cell, 0


def safe_sheet_name(name, used=None):
    """Excel rejects some characters and anything past 31 chars."""

    cleaned = _ILLEGAL_SHEET_CHARS.sub("-", str(name)).strip() or "Sheet"
    cleaned = cleaned[:31]
    if used is None:
        return cleaned
    candidate, suffix = cleaned, 2
    while candidate in used:
        tail = f"_{suffix}"
        candidate = cleaned[: 31 - len(tail)] + tail
        suffix += 1
    used.add(candidate)
    return candidate


def _column_widths(rows):
    widths = []
    max_cols = max((len(row) for row in rows), default=0)
    for column_index in range(max_cols):
        lengths = []
        for row in rows:
            if column_index < len(row):
                value = row[column_index]
                if isinstance(value, tuple) and len(value) == 2:
                    value = value[0]
                lengths.append(len(str("" if value is None else value)))
        widths.append(min(max(max(lengths, default=8) + 2, 10), 45))
    return widths


def write_xlsx(path, sheets):
    """Write sheets to an XLSX workbook.

    ``sheets`` is a list of ``(name, rows)`` or ``(name, rows, kind)``.  A row
    is a list of plain values or ``(value, style_id)`` pairs.
    """

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    workbook = xlsxwriter.Workbook(
        path, {"constant_memory": False, "default_date_format": "yyyy-mm-dd"}
    )
    try:
        workbook.set_properties(
            {
                "title": "nagini analysis report",
                "author": "nagini",
                "company": "Clinical Genomics Uppsala",
                "created": datetime.now(timezone.utc).replace(tzinfo=None),
            }
        )
        formats = {style: workbook.add_format(spec) for style, spec in STYLE_SPECS.items()}
        used_names = set()

        for sheet in sheets:
            name, rows = sheet[0], sheet[1]
            kind = sheet[2] if len(sheet) > 2 else TABLE
            worksheet = workbook.add_worksheet(safe_sheet_name(name, used_names))

            for column_index, width in enumerate(_column_widths(rows)):
                worksheet.set_column(column_index, column_index, width)

            for row_index, row in enumerate(rows, 1):
                for column_index, cell in enumerate(row, 1):
                    value, style = _cell_value_and_style(cell, row_index, column_index)
                    cell_format = formats.get(style, formats[0])
                    if isinstance(value, bool) or value is None:
                        worksheet.write_string(row_index - 1, column_index - 1,
                                               "" if value is None else str(value), cell_format)
                    elif isinstance(value, (int, float)):
                        worksheet.write_number(row_index - 1, column_index - 1, value, cell_format)
                    else:
                        worksheet.write_string(
                            row_index - 1,
                            column_index - 1,
                            _CONTROL_CHARS.sub("", str(value)),
                            cell_format,
                        )

            if kind == TABLE and len(rows) > 1:
                worksheet.freeze_panes(1, 0)
                worksheet.autofilter(0, 0, len(rows) - 1, max(len(rows[0]) - 1, 0))
    finally:
        workbook.close()


def build_bcr_abl_sheet(sample, fusion_rows, mutation_rows, phasing_rows, exon_rows=None):
    layout = [
        [("BCR-ABL1 FUSION & RESISTANCE MUTATION REPORT", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2)],
        [("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0)],
        [("1. BCR-ABL1 FUSION ISOFORM BREAKDOWN", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3)],
        [
            ("Fusion Isoform", 1), ("Alt Name", 1), ("Protein Form", 1), ("Cluster", 1),
            ("BCR Exon", 1), ("ABL1 Exon", 1), ("Read Count", 1), ("Percent (%)", 1), ("Status", 1), ("Description", 1)
        ],
    ]
    if not fusion_rows:
        layout.append([("No BCR-ABL1 fusion events detected", 0)] + [("", 0)] * 9)
    else:
        for r in fusion_rows:
            st = str(r.get("status", "")).upper()
            st_style = 5 if st == "PRIMARY" else (6 if st == "SECONDARY" else 0)
            layout.append([
                (r.get("fusion_isoform", ""), 8 if st == "PRIMARY" else 0),
                (r.get("alt_name", ""), 0),
                (r.get("protein_form", ""), 0),
                (r.get("cluster", ""), 0),
                (r.get("bcr_exon", ""), 0),
                (r.get("abl_exon", ""), 0),
                (int(r.get("read_count", 0)), 0),
                (float(r.get("fraction_percent", 0.0)), 0),
                (st, st_style),
                (r.get("description", ""), 0),
            ])

    layout.append([("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0)])
    layout.append([("2. ABL1 KINASE DOMAIN RESISTANCE MUTATIONS (BRANFORD PANEL)", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3)])
    layout.append([
        ("Mutation", 1), ("AA Change", 1), ("Position", 1), ("Ref", 1), ("Alt", 1),
        ("Depth", 1), ("Alt Reads", 1), ("VAF (%)", 1), ("Status", 1),
        ("Imatinib", 1), ("Dasatinib", 1), ("Nilotinib", 1), ("Bosutinib", 1), ("Ponatinib", 1), ("Asciminib", 1)
    ])
    if not mutation_rows:
        layout.append([("No candidate resistance mutations called", 0)] + [("", 0)] * 14)
    else:
        for r in mutation_rows:
            mut_status = str(r.get("status", "")).upper()
            mut_style = 4 if mut_status == "DETECTED" else (6 if mut_status == "SUBCLONAL" else 0)
            layout.append([
                (r.get("mutation", ""), 8 if mut_status == "DETECTED" else 0),
                (r.get("aa_change", ""), 0),
                (int(r.get("pos", 0)), 0),
                (r.get("ref", ""), 0),
                (r.get("alt", ""), 0),
                (int(r.get("depth", 0)), 0),
                (int(r.get("ad_alt", 0)), 0),
                (float(r.get("vaf_percent", 0.0)), mut_style),
                (mut_status, mut_style),
                (r.get("imatinib", "Unknown"), _style_tki_response(r.get("imatinib", ""))),
                (r.get("dasatinib", "Unknown"), _style_tki_response(r.get("dasatinib", ""))),
                (r.get("nilotinib", "Unknown"), _style_tki_response(r.get("nilotinib", ""))),
                (r.get("bosutinib", "Unknown"), _style_tki_response(r.get("bosutinib", ""))),
                (r.get("ponatinib", "Unknown"), _style_tki_response(r.get("ponatinib", ""))),
                (r.get("asciminib", "Unknown"), _style_tki_response(r.get("asciminib", ""))),
            ])

    layout.append([("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0)])
    layout.append([("3. COMPOUND MUTATION PHASING (SINGLE-MOLECULE LONG-READ)", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3)])
    layout.append([
        ("Mutation 1", 1), ("Mutation 2", 1), ("Phase State", 1), ("Cis (%)", 1),
        ("Cis Reads (11)", 1), ("Trans Reads (10)", 1), ("Trans Reads (01)", 1), ("WT Reads (00)", 1),
        ("Spanning Depth", 1), ("Clinical Interpretation", 1)
    ])
    if not phasing_rows:
        layout.append([("No compound mutation pairs to phase (< 2 mutations detected)", 0)] + [("", 0)] * 9)
    else:
        for r in phasing_rows:
            ps = str(r.get("phase_state", ""))
            ps_style = 4 if "CIS" in ps.upper() else (6 if "TRANS" in ps.upper() else 0)
            layout.append([
                (r.get("mutation_1", ""), 0),
                (r.get("mutation_2", ""), 0),
                (ps, ps_style),
                (float(r.get("cis_fraction_percent", 0.0)), ps_style),
                (int(r.get("cis_reads_11", 0)), 0),
                (int(r.get("trans_reads_10", 0)), 0),
                (int(r.get("trans_reads_01", 0)), 0),
                (int(r.get("wt_reads_00", 0)), 0),
                (int(r.get("total_spanning_reads", 0)), 0),
                (r.get("clinical_interpretation", ""), 0),
            ])

    if exon_rows:
        layout.append([("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0)])
        layout.append([("4. BCR & ABL1 EXON COVERAGE STATISTICS", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3)])
        layout.append([
            ("Exon", 1), ("Chrom", 1), ("Start", 1), ("End", 1), ("Length (bp)", 1),
            ("Mean Depth", 1), ("Median Depth", 1), ("Min Depth", 1), ("Max Depth", 1),
            (">=100x (%)", 1), (">=500x (%)", 1), (">=1000x (%)", 1), ("Spanning Reads", 1)
        ])
        for r in exon_rows:
            mean_d = float(r.get("mean_depth", 0.0))
            d_style = 5 if mean_d >= 500 else (6 if mean_d >= 100 else 0)
            layout.append([
                (r.get("exon", ""), 0),
                (r.get("chrom", ""), 0),
                (int(r.get("start", 0)), 0),
                (int(r.get("end", 0)), 0),
                (int(r.get("length_bp", 0)), 0),
                (mean_d, d_style),
                (float(r.get("median_depth", 0.0)), 0),
                (int(r.get("min_depth", 0)), 0),
                (int(r.get("max_depth", 0)), 0),
                (float(r.get("pct_gte_100x", 0.0)), 0),
                (float(r.get("pct_gte_500x", 0.0)), 0),
                (float(r.get("pct_gte_1000x", 0.0)), 0),
                (int(r.get("spanning_reads", 0)), 0),
            ])

    return layout


def generate_report(output, sample, sample_type, reference_build, annotation, snv_vcf,
                    sv_tsv, translocations_tsv, wgs_coverage, target_coverage, cramino,
                    cytoband="", caller="", design_bed="", platform="",
                    bcr_abl_fusion="", bcr_abl_mutations="", bcr_abl_phasing="",
                    bcr_abl_exon_coverage=""):
    snv_columns, snvs, detected_caller = parse_small_variant_vcf(snv_vcf, default_caller=caller)
    sv_columns, svs = parse_tsv(sv_tsv)
    translocation_columns, translocations = parse_tsv(translocations_tsv)
    cytobands = CytobandLookup(cytoband)
    sv_columns, svs = annotate_sv_rows(sv_columns, svs, cytobands)
    translocation_columns, translocations = annotate_sv_rows(
        translocation_columns, translocations, cytobands
    )
    wgs_columns, wgs_rows = parse_tsv(wgs_coverage)
    target_columns, target_rows = parse_tsv(target_coverage)
    cramino_columns, cramino_rows = parse_cramino(cramino)

    bcr_abl_fusion_rows = parse_bcr_abl_tsv(bcr_abl_fusion)
    bcr_abl_mutation_rows = parse_bcr_abl_tsv(bcr_abl_mutations)
    bcr_abl_phasing_rows = parse_bcr_abl_tsv(bcr_abl_phasing)
    bcr_abl_exon_rows = parse_bcr_abl_tsv(bcr_abl_exon_coverage)
    has_bcr_abl = bool(bcr_abl_fusion_rows or bcr_abl_mutation_rows or bcr_abl_phasing_rows or bcr_abl_exon_rows)


    pass_count = sum(record["Filter"] in {"PASS", "."} for record in snvs)
    pathogenic_snvs = [
        r for r in snvs if _is_pathogenic_clinvar(r.get("ClinVar_Significance", ""))
    ]
    high_impact_snvs = [
        r for r in snvs if str(r.get("Impact", "")).upper() == "HIGH"
    ]
    target_svs = [
        r for r in svs if r.get("genes") and r.get("genes") != "."
    ]

    cramino_dict = {row["Metric"]: row["Value"] for row in cramino_rows}
    yield_gb = _find_metric(cramino_dict, "Yield [Gb]", "Yield", "Bases")
    total_reads = _find_metric(cramino_dict, "Number of reads", "Reads")
    mean_q = _find_metric(cramino_dict, "Mean quality", "Quality", "Q-score")
    n50_bp = _find_metric(cramino_dict, "Read length N50", "N50")
    median_len = _find_metric(cramino_dict, "Median read length", "Median length")

    wgs_cov = next(
        (str(row.get("mean", "")) for row in wgs_rows if str(row.get("chrom", "")).lower() == "total"),
        (str(wgs_rows[0].get("mean", "")) if wgs_rows else "")
    )
    target_cov = next(
        (str(row.get("mean", "")) for row in target_rows if "total" in str(row.get("chrom", "")).lower()),
        (str(target_rows[0].get("mean", "")) if target_rows else "")
    )

    try:
        q_val = float(mean_q.split()[0])
        q_status, q_style = ("PASS", 5) if q_val >= 12.0 else (("BORDERLINE", 6) if q_val >= 10.0 else ("LOW", 4))
    except (ValueError, TypeError, IndexError):
        q_status, q_style = ("INFO", 0)

    try:
        t_val = float(target_cov.split()[0])
        t_status, t_style = ("PASS", 5) if t_val >= 30.0 else (("BORDERLINE", 6) if t_val >= 15.0 else ("LOW", 4))
    except (ValueError, TypeError, IndexError):
        t_status, t_style = ("INFO", 0)

    caller_display = detected_caller
    if platform and "clairs" in caller_display.lower():
        caller_display = f"ClairS-TO ({platform})"
    elif not caller_display:
        caller_display = "ClairS-TO"

    target_display = os.path.basename(design_bed) if design_bed else "Whole Genome (WGS)"

    # Fluffy-Style Overview Layout with rich sections and status color coding
    overview = [
        [("NAGINI NANOPORE ANALYSIS REPORT", 2), ("", 2), ("", 2), ("", 2)],
        [("", 0), ("", 0), ("", 0), ("", 0)],
        [("SAMPLE & PIPELINE METADATA", 3), ("", 3), ("", 3), ("", 3)],
        [("Sample ID", 7), (sample, 0), ("Sequencing Platform", 7), ("Oxford Nanopore (ONT)", 0)],
        [("Sample Type", 7), (f"Tumor ({sample_type})" if sample_type == "T" else f"Normal ({sample_type})", 0), ("Reference Genome", 7), (reference_build, 0)],
        [("Small Variant Caller", 7), (caller_display, 0), ("Structural Variant Caller", 7), ("Sniffles2", 0)],
        [("Variant Annotation", 7), (annotation, 0), ("Target Region / Panel", 7), (target_display, 0)],
        [("Report Generated", 7), (datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"), 0), ("Pipeline Status", 7), ("COMPLETED", 5)],
        [("", 0), ("", 0), ("", 0), ("", 0)],
        [("QUALITY CONTROL & SEQUENCING METRICS", 3), ("", 3), ("VALUE", 3), ("STATUS", 3)],
        [("Total Yield [Gb]", 7), ("Total sequenced yield", 0), (yield_gb or "N/A", 0), ("PASS" if yield_gb else "INFO", 5 if yield_gb else 0)],
        [("Total Reads", 7), ("Total sequenced reads", 0), (total_reads or "N/A", 0), ("PASS" if total_reads else "INFO", 5 if total_reads else 0)],
        [("Mean Base Quality (Q)", 7), ("Average basecall Q-score", 0), (mean_q or "N/A", 0), (q_status, q_style)],
        [("Read Length N50 [bp]", 7), ("N50 read length", 0), (n50_bp or "N/A", 0), ("PASS" if n50_bp else "INFO", 5 if n50_bp else 0)],
        [("Median Read Length [bp]", 7), ("Median read length", 0), (median_len or "N/A", 0), ("PASS" if median_len else "INFO", 5 if median_len else 0)],
        [("Mean Target Coverage", 7), ("Mean depth in target BED", 0), (target_cov or "N/A", 0), (t_status, t_style)],
        [("Mean WGS Coverage", 7), ("Mean whole-genome depth", 0), (wgs_cov or "N/A", 0), ("INFO", 0)],
        [("", 0), ("", 0), ("", 0), ("", 0)],
        [("VARIANT CALLING SUMMARY", 3), ("", 3), ("COUNT", 3), ("FLAG / ACTION", 3)],
        [("Total Small Variants", 7), ("SNVs and short indels", 0), (len(snvs), 0), ("PASS", 5)],
        [("PASS Small Variants", 7), ("High-confidence calls", 0), (pass_count, 0), ("PASS", 5)],
        [("ClinVar Pathogenic Variants", 7), ("Pathogenic / Likely Pathogenic", 0), (len(pathogenic_snvs), 0), ("ACTION REQUIRED" if pathogenic_snvs else "NONE", 4 if pathogenic_snvs else 5)],
        [("High Impact Variants", 7), ("Stop gained / frameshift / splice", 0), (len(high_impact_snvs), 0), ("REVIEW" if high_impact_snvs else "NONE", 6 if high_impact_snvs else 5)],
        [("Total Structural Variants", 7), ("Sniffles2 SV calls", 0), (len(svs), 0), ("INFO", 0)],
        [("Target-Overlapping SVs", 7), ("SVs overlapping panel genes", 0), (len(target_svs), 0), ("REVIEW" if target_svs else "NONE", 6 if target_svs else 5)],
        [("Translocations / Breakends", 7), ("Inter-chromosomal rearrangements", 0), (len(translocations), 0), ("TRANSLOCATION DETECTED" if translocations else "NONE", 8 if translocations else 5)],
    ]

    # Add BCR-ABL summary block to overview if BCR-ABL data present
    if has_bcr_abl:
        prim_iso = next((r["fusion_isoform"] for r in bcr_abl_fusion_rows if r.get("status") == "PRIMARY"), "None detected")
        det_muts = [r["mutation"] for r in bcr_abl_mutation_rows if r.get("status") in ("DETECTED", "SUBCLONAL")]
        mut_summary = ", ".join(det_muts) if det_muts else "None detected"
        cis_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in bcr_abl_phasing_rows if "CIS" in r.get("phase_state", "").upper()]
        phase_summary = f"CIS ({', '.join(cis_pairs)})" if cis_pairs else ("TRANS" if bcr_abl_phasing_rows else "N/A (<2 muts)")

        overview.append([("", 0), ("", 0), ("", 0), ("", 0)])
        overview.append([("BCR-ABL1 TARGETED ASSAY SUMMARY", 3), ("", 3), ("", 3), ("", 3)])
        overview.append([("Primary Fusion Isoform", 7), ("Major/Minor BCR breakpoint", 0), (prim_iso, 8), ("PASS" if prim_iso != "None detected" else "INFO", 5 if prim_iso != "None detected" else 0)])
        overview.append([("Resistance Mutations", 7), ("Branford ABL1 panel", 0), (mut_summary, 4 if det_muts else 0), ("ACTION REQUIRED" if det_muts else "NONE", 4 if det_muts else 5)])
        overview.append([("Compound Phasing Status", 7), ("Single-molecule read phasing", 0), (phase_summary, 4 if cis_pairs else 0), ("COMPOUND MUTATION (CIS)" if cis_pairs else ("POLYCLONAL (TRANS)" if bcr_abl_phasing_rows else "PASS"), 4 if cis_pairs else 5)])

    # Add clinical findings preview tables directly on the Overview tab
    if pathogenic_snvs:
        overview.append([("", 0), ("", 0), ("", 0), ("", 0)])
        overview.append([("PATHOGENIC & LIKELY PATHOGENIC VARIANTS (CLINVAR)", 4), ("", 4), ("", 4), ("", 4), ("", 4), ("", 4), ("", 4)])
        overview.append([("Gene", 1), ("Chr", 1), ("Position", 1), ("Change", 1), ("VAF", 1), ("Consequence", 1), ("ClinVar Significance", 1)])
        for v in pathogenic_snvs[:15]:
            overview.append([
                (v.get("Gene", ""), 0),
                (v.get("Chromosome", ""), 0),
                (v.get("Position", 0), 0),
                (f"{v.get('Reference', '')}>{v.get('Alternate', '')}", 0),
                (v.get("VAF", ""), 0),
                (v.get("Consequence", ""), 0),
                (v.get("ClinVar_Significance", ""), 4),
            ])

    if translocations:
        overview.append([("", 0), ("", 0), ("", 0), ("", 0)])
        overview.append([("TRANSLOCATIONS & REARRANGEMENTS (SNIFFLES2)", 8), ("", 8), ("", 8), ("", 8), ("", 8), ("", 8), ("", 8)])
        overview.append([("Sample", 1), ("Type", 1), ("Chr1", 1), ("Pos1", 1), ("Chr2", 1), ("Notation", 1), ("Genes", 1)])
        for t in translocations[:15]:
            overview.append([
                (t.get("sample", sample), 0),
                (t.get("svtype", "BND"), 0),
                (t.get("chrom", ""), 0),
                (t.get("pos", 0), 0),
                (t.get("chrom2", ""), 0),
                (t.get("Cytogenetic_Notation", "") or t.get("svtype", ""), 8),
                (t.get("genes", ""), 0),
            ])

    sheets = [
        ("Overview", overview, LAYOUT),
    ]
    if has_bcr_abl:
        sheets.append(("BCR-ABL1", build_bcr_abl_sheet(sample, bcr_abl_fusion_rows, bcr_abl_mutation_rows, bcr_abl_phasing_rows, bcr_abl_exon_rows), LAYOUT))

    sheets.extend([
        ("SNVs and indels", table(snv_columns, snvs)),
        ("SVs", table(sv_columns, svs)),
        ("Translocations", table(translocation_columns, translocations)),
        ("WGS coverage", table(wgs_columns, wgs_rows)),
        ("Target coverage", table(target_columns, target_rows)),
        ("Cramino QC", table(cramino_columns, cramino_rows)),
    ])
    write_xlsx(output, sheets)


def _first_or_str(val):
    if isinstance(val, (list, tuple)):
        return str(val[0]) if val else ""
    return str(val) if val else ""


def run_from_snakemake(snakemake):
    generate_report(
        output=str(snakemake.output.xlsx),
        sample=str(snakemake.wildcards.sample),
        sample_type=str(snakemake.wildcards.type),
        reference_build=str(snakemake.params.reference_build),
        annotation=str(snakemake.params.snv_annotation),
        snv_vcf=_first_or_str(getattr(snakemake.input, "snv_vcf", "")),
        sv_tsv=_first_or_str(getattr(snakemake.input, "sv", "")),
        translocations_tsv=_first_or_str(getattr(snakemake.input, "translocations", "")),
        wgs_coverage=_first_or_str(getattr(snakemake.input, "wgs_coverage", "")),
        target_coverage=_first_or_str(getattr(snakemake.input, "target_coverage", "")),
        cramino=_first_or_str(getattr(snakemake.input, "cramino", "")),
        cytoband=str(getattr(snakemake.params, "cytoband", "") or ""),
        caller=str(getattr(snakemake.params, "small_variant_caller", "") or ""),
        design_bed=str(getattr(snakemake.params, "design_bed", "") or ""),
        platform=str(getattr(snakemake.params, "clairs_to_platform", "") or ""),
        bcr_abl_fusion=_first_or_str(getattr(snakemake.input, "bcr_abl_fusion", "")),
        bcr_abl_mutations=_first_or_str(getattr(snakemake.input, "bcr_abl_mutations", "")),
        bcr_abl_phasing=_first_or_str(getattr(snakemake.input, "bcr_abl_phasing", "")),
        bcr_abl_exon_coverage=_first_or_str(getattr(snakemake.input, "bcr_abl_exon_coverage", "")),
    )


def main():
    parser = argparse.ArgumentParser(description="Create a nagini XLSX report")
    parser.add_argument("--output", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--sample-type", required=True)
    parser.add_argument("--reference-build", required=True)
    parser.add_argument("--annotation", default="none")
    parser.add_argument("--snv-vcf", required=True)
    parser.add_argument("--sv-tsv", required=True)
    parser.add_argument("--translocations-tsv", required=True)
    parser.add_argument("--wgs-coverage", required=True)
    parser.add_argument("--target-coverage", required=True)
    parser.add_argument("--cramino", required=True)
    parser.add_argument("--cytoband", default="")
    parser.add_argument("--caller", default="")
    parser.add_argument("--design-bed", default="")
    parser.add_argument("--platform", default="")
    parser.add_argument("--bcr-abl-fusion", default="")
    parser.add_argument("--bcr-abl-mutations", default="")
    parser.add_argument("--bcr-abl-phasing", default="")
    parser.add_argument("--bcr-abl-exon-coverage", default="")
    args = parser.parse_args()
    generate_report(
        args.output, args.sample, args.sample_type, args.reference_build, args.annotation,
        args.snv_vcf, args.sv_tsv, args.translocations_tsv, args.wgs_coverage,
        args.target_coverage, args.cramino, args.cytoband, args.caller,
        args.design_bed, args.platform,
        bcr_abl_fusion=args.bcr_abl_fusion,
        bcr_abl_mutations=args.bcr_abl_mutations,
        bcr_abl_phasing=args.bcr_abl_phasing,
        bcr_abl_exon_coverage=args.bcr_abl_exon_coverage,
    )


if "snakemake" in globals():
    run_from_snakemake(snakemake)
elif __name__ == "__main__":
    main()

