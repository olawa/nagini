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

try:
    import xlsxwriter
except ImportError:
    xlsxwriter = None


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


def parse_mosdepth_dist(path):
    """Parse mosdepth global.dist.txt into a depth -> proportion mapping."""
    if not path or not os.path.exists(path):
        return {}
    dist = {}
    with open_text(path) as handle:
        for line in handle:
            parts = line.strip().split("\t")
            if len(parts) >= 3 and parts[0].lower() == "total":
                try:
                    cutoff = int(parts[1])
                    prop = float(parts[2])
                    dist[cutoff] = prop
                except ValueError:
                    continue
    return dist


def parse_demux_info(path, sample="", barcode=""):
    """Extract sample/barcode demux read count and % of run from demux TSV."""
    if not path or not os.path.exists(path):
        return {"barcode": barcode, "reads": "", "fraction": "", "status": ""}
    norm_sample = str(sample or "").strip().lower()
    norm_barcode = str(barcode or "").strip().lower()
    match_row = None
    with open_text(path) as handle:
        reader = csv.DictReader(
            (line for line in handle if not line.startswith("#")),
            delimiter="\t",
        )
        for row in reader:
            r_samp = str(row.get("Sample", row.get("sample", "")) or "").strip().lower()
            r_bar = str(row.get("Barcode", row.get("barcode", "")) or "").strip().lower()
            if (norm_sample and r_samp == norm_sample) or (norm_barcode and r_bar == norm_barcode):
                match_row = row
                break
    if match_row:
        reads = match_row.get("Reads", match_row.get("reads", ""))
        frac = match_row.get("Percent of run", match_row.get("percent_of_run", ""))
        bar = match_row.get("Barcode", match_row.get("barcode", barcode))
        return {
            "barcode": bar or barcode,
            "reads": reads,
            "fraction": f"{frac}%" if frac and "%" not in str(frac) else str(frac),
            "status": match_row.get("Status", match_row.get("status", "OK")),
        }
    return {"barcode": barcode, "reads": "", "fraction": "", "status": ""}


def parse_small_variant_vcf(path, default_caller=""):
    columns = [
        "Caller", "Chromosome", "Position", "Reference", "Alternate", "Filter", "Quality",
        "VAF", "Depth", "Allele depth", "Gene", "Consequence", "Impact",
        "HGVSc", "HGVSp", "Existing variation", "gnomAD_AF", "ClinVar_Significance",
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

            gnomad_af = _first_annotation_value(
                annotation, "gnomAD_AF", "gnomADg_AF", "gnomADe_AF", "MAX_AF", "AF"
            ) or str(info.get("gnomAD_AF", info.get("AF", "")))
            try:
                gnomad_float = float(str(gnomad_af).split(",")[0])
                gnomad_display = f"{gnomad_float:.5f}"
            except (ValueError, TypeError, AttributeError):
                gnomad_display = gnomad_af

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
                "gnomAD_AF": gnomad_display,
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


def merge_caller_rows(*sources):
    """Combine per-caller tables into one, with a leading Caller column.

    Column sets can differ between callers -- an older summary may lack a
    column a newer one writes -- so the union is taken and missing values are
    left blank rather than shifting a row's fields.
    """

    columns = ["caller"]
    for source_columns, _rows, _name in sources:
        for column in source_columns:
            if column != "caller" and column not in columns:
                columns.append(column)

    merged = []
    for _source_columns, rows, name in sources:
        for row in rows:
            new_row = dict(row)
            new_row["caller"] = str(row.get("caller", "") or name)
            merged.append({column: new_row.get(column, "") for column in columns})
    return columns, merged


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
    """Lookup a metric from Cramino or QC summary dict with priority matching."""
    # 1. Exact match (case-insensitive)
    for kw in keywords:
        kw_lower = kw.strip().lower()
        for k, v in data_dict.items():
            if k.strip().lower() == kw_lower:
                return str(v).strip()

    # 2. Match normalized names (without units/brackets)
    for kw in keywords:
        kw_norm = re.sub(r"\[.*?\]|\(.*?\)", "", kw).strip().lower()
        if not kw_norm:
            continue
        for k, v in data_dict.items():
            k_norm = re.sub(r"\[.*?\]|\(.*?\)", "", k).strip().lower()
            if k_norm == kw_norm:
                return str(v).strip()

    # 3. Substring match, skipping percentage/fraction fields when searching for counts
    for kw in keywords:
        kw_lower = kw.strip().lower()
        for k, v in data_dict.items():
            k_lower = k.strip().lower()
            if kw_lower in ("reads", "total reads", "number of reads"):
                if "aligned" in k_lower or "%" in k_lower or "fraction" in k_lower:
                    continue
            if kw_lower in k_lower:
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


class ReportLayout(list):
    """List subclass carrying optional embedded image metadata."""

    def __init__(self, items=None, image_info=None):
        super().__init__(items or [])
        self.image_info = image_info


def write_xlsx(path, sheets):
    """Write sheets to an XLSX workbook.

    ``sheets`` is a list of ``(name, rows)`` or ``(name, rows, kind)`` or
    ``(name, rows, kind, image_info)``. A row is a list of plain values or
    ``(value, style_id)`` pairs.
    """
    if xlsxwriter is None:
        raise RuntimeError("xlsxwriter is required to write XLSX workbooks")

    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    workbook = xlsxwriter.Workbook(
        path, {"constant_memory": False, "default_date_format": "yyyy-mm-dd"}
    )
    workbook.set_size(1600, 1200)
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
            image_info = getattr(rows, "image_info", None) or (sheet[3] if len(sheet) > 3 else None)
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

            if image_info and isinstance(image_info, dict) and os.path.isfile(image_info.get("path", "")) and os.path.getsize(image_info.get("path", "")) > 0:
                cell = image_info.get("cell", f"A{len(rows) + 3}")
                x_scale = float(image_info.get("x_scale", 0.68))
                y_scale = float(image_info.get("y_scale", 0.68))
                try:
                    worksheet.insert_image(cell, image_info["path"], {"x_scale": x_scale, "y_scale": y_scale})

                except Exception as exc:
                    sys.stderr.write(f"Warning: could not insert image {image_info['path']}: {exc}\n")

            if kind == TABLE and len(rows) > 1:
                worksheet.freeze_panes(1, 0)
                worksheet.autofilter(0, 0, len(rows) - 1, max(len(rows[0]) - 1, 0))
    finally:
        workbook.close()


def build_bcr_abl_sheet(sample, fusion_rows, mutation_rows, phasing_rows, exon_rows=None, plot_png=None):
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
    layout.append([("2. ABL1 KINASE DOMAIN RESISTANCE MUTATIONS (BRANFORD PANEL)", 3)] + [("", 3)] * 17)
    layout.append([
        ("Mutation", 1), ("AA Change", 1), ("Position", 1), ("Ref", 1), ("Alt", 1),
        ("Depth", 1), ("Alt Reads", 1), ("VAF (%)", 1), ("Status", 1),
        ("SNR (x)", 1), ("Strand Ratio", 1), ("Mean BaseQ", 1),
        ("Imatinib", 1), ("Dasatinib", 1), ("Nilotinib", 1), ("Bosutinib", 1), ("Ponatinib", 1), ("Asciminib", 1)
    ])
    if not mutation_rows:
        layout.append([("No candidate resistance mutations called", 0)] + [("", 0)] * 17)
    else:
        for r in mutation_rows:
            mut_status = str(r.get("status", "")).upper()
            mut_style = 4 if mut_status == "DETECTED" else (6 if mut_status == "SUBCLONAL" else 0)

            snr_val = r.get("snr", "")
            try:
                snr_display = f"{float(snr_val):.1f}x" if snr_val != "" and snr_val is not None else "N/A"
            except (ValueError, TypeError):
                snr_display = str(snr_val)

            strand_ratio_val = r.get("alt_minor_strand_ratio", "")
            try:
                strand_display = f"{float(strand_ratio_val):.2f}" if strand_ratio_val != "" and strand_ratio_val is not None else "N/A"
            except (ValueError, TypeError):
                strand_display = str(strand_ratio_val)

            baseq_val = r.get("mean_baseq_alt", "")
            try:
                baseq_display = f"{float(baseq_val):.1f}" if baseq_val != "" and baseq_val is not None else "N/A"
            except (ValueError, TypeError):
                baseq_display = str(baseq_val)

            layout.append([
                (r.get("mutation", ""), 8 if mut_status == "DETECTED" else (6 if mut_status == "SUBCLONAL" else 0)),
                (r.get("aa_change", ""), 0),
                (int(r.get("pos", 0)), 0),
                (r.get("ref", ""), 0),
                (r.get("alt", ""), 0),
                (int(r.get("depth", 0)), 0),
                (int(r.get("ad_alt", 0)), 0),
                (float(r.get("vaf_percent", 0.0)), mut_style),
                (mut_status, mut_style),
                (snr_display, 0),
                (strand_display, 0),
                (baseq_display, 0),
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

    if plot_png and os.path.isfile(plot_png) and os.path.getsize(plot_png) > 0:
        layout.append([("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0)])
        layout.append([("5. BCR-ABL1 FUSION ARCHITECTURE & JUNCTION PLOT (ARRIBA)", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3)])
        img_cell = f"B{len(layout) + 2}"
        return ReportLayout(layout, image_info={"path": plot_png, "cell": img_cell, "x_scale": 0.68, "y_scale": 0.68})


    return layout



def generate_report(output, sample, sample_type, reference_build, annotation, snv_vcf,
                    sv_tsv, translocations_tsv, wgs_coverage, target_coverage, cramino,
                    cytoband="", caller="", design_bed="", platform="",
                    target_dist="", demux_summary="", barcode="",
                    bcr_abl_fusion="", bcr_abl_mutations="", bcr_abl_phasing="",
                    bcr_abl_exon_coverage="", bcr_abl_plot_png="", cnv_tsv="",
                    rindels_sv_tsv="", rindels_translocations_tsv=""):

    snv_columns, snvs, detected_caller = parse_small_variant_vcf(snv_vcf, default_caller=caller)
    sv_columns, svs = parse_tsv(sv_tsv)
    translocation_columns, translocations = parse_tsv(translocations_tsv)
    cnv_columns, cnvs = parse_tsv(cnv_tsv)
    cytobands = CytobandLookup(cytoband)
    sv_columns, svs = annotate_sv_rows(sv_columns, svs, cytobands)
    translocation_columns, translocations = annotate_sv_rows(
        translocation_columns, translocations, cytobands
    )

    # Two callers reach the same breakpoints by different evidence: Sniffles2
    # from split and spanning reads, Rindels --sv-lite from pair and clip
    # evidence without assembly. One table with a Caller column lets a
    # dic(9;20) call be corroborated instead of read from one caller alone.
    rindels_sv_columns, rindels_svs = parse_tsv(rindels_sv_tsv)
    rindels_sv_columns, rindels_svs = annotate_sv_rows(
        rindels_sv_columns, rindels_svs, cytobands
    )
    rindels_trans_columns, rindels_translocations = parse_tsv(rindels_translocations_tsv)
    rindels_trans_columns, rindels_translocations = annotate_sv_rows(
        rindels_trans_columns, rindels_translocations, cytobands
    )
    sv_columns, svs = merge_caller_rows(
        (sv_columns, svs, "Sniffles2"), (rindels_sv_columns, rindels_svs, "Rindels")
    )
    translocation_columns, translocations = merge_caller_rows(
        (translocation_columns, translocations, "Sniffles2"),
        (rindels_trans_columns, rindels_translocations, "Rindels"),
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

    # Partition SVs into intra-chromosomal (DEL, DUP, INS, INV) vs inter-chromosomal / breakends
    svs_intra = [
        r for r in svs
        if str(r.get("translocation", "")).lower() not in ("true", "1")
        and str(r.get("svtype", "")).upper() not in ("BND", "TRA")
    ]
    target_svs = [
        r for r in svs_intra if r.get("genes") and r.get("genes") != "."
    ]

    cnv_non_diploid = [
        c for c in cnvs if str(c.get("event_type", "")).upper() not in ("DIPLOID", "")
    ]
    target_cnvs = [
        c for c in cnv_non_diploid
        if str(c.get("is_target_region", "")).lower() in ("true", "1")
    ]

    cramino_dict = {row["Metric"]: row["Value"] for row in cramino_rows}
    yield_gb = _find_metric(cramino_dict, "Yield [Gb]", "Yield", "Bases")
    total_reads = _find_metric(cramino_dict, "Number of reads", "Total reads", "Read count")
    mean_q = _find_metric(cramino_dict, "Mean read quality", "Mean quality", "Average read quality", "Mean base quality", "Median read quality", "Quality", "Q-score")
    n50_bp = _find_metric(cramino_dict, "Read length N50", "N50", "N50 read length")
    median_len = _find_metric(cramino_dict, "Median read length", "Median length", "Median")


    wgs_cov = next(
        (str(row.get("mean", "")) for row in wgs_rows if str(row.get("chrom", "")).lower() == "total"),
        (str(wgs_rows[0].get("mean", "")) if wgs_rows else "")
    )
    # mosdepth --by writes both a whole-reference "total" row and an in-region
    # "total_region" row, and "total" sorts first.  Substring-matching "total"
    # therefore picked the genome-wide mean as the target mean, which made the
    # target depth read like the background and fold enrichment like 1.0x.
    target_cov = next(
        (
            str(row.get("mean", ""))
            for row in target_rows
            if str(row.get("chrom", "")).strip().lower() == "total_region"
        ),
        next(
            (
                str(row.get("mean", ""))
                for row in target_rows
                if str(row.get("chrom", "")).strip().lower() == "total"
            ),
            (str(target_rows[0].get("mean", "")) if target_rows else ""),
        ),
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
        t_val = None
        t_status, t_style = ("INFO", 0)

    target_dist_map = parse_mosdepth_dist(target_dist)
    pct_20x = target_dist_map.get(20)
    pct_30x = target_dist_map.get(30)
    pct_20x_str = f"{pct_20x * 100:.1f}%" if pct_20x is not None else "N/A"
    pct_30x_str = f"{pct_30x * 100:.1f}%" if pct_30x is not None else "N/A"
    if pct_20x is not None:
        p20_val = pct_20x * 100.0
        p20_status, p20_style = (
            ("PASS", 5) if p20_val >= 90.0 else (("BORDERLINE", 6) if p20_val >= 80.0 else ("LOW", 4))
        )
    else:
        p20_status, p20_style = ("INFO", 0)

    total_fusion_reads = sum(int(r.get("read_count", 0)) for r in bcr_abl_fusion_rows)
    raw_reads_num = 0
    try:
        raw_reads_num = int(float(str(total_reads).replace(",", "")))
    except (ValueError, TypeError):
        raw_reads_num = 0

    if raw_reads_num > 0 and total_fusion_reads > 0:
        target_enrich_pct = (total_fusion_reads / raw_reads_num) * 100.0
        target_enrich_str = f"{target_enrich_pct:.2f}%"
        enrich_status, enrich_style = ("PASS", 5) if target_enrich_pct >= 0.5 else ("INFO", 0)
    else:
        target_enrich_str = "N/A"
        enrich_status, enrich_style = ("INFO", 0)

    if bcr_abl_exon_rows:
        total_exon_len = sum(int(r.get("length", 0)) for r in bcr_abl_exon_rows)
        if total_exon_len > 0:
            avg_exon_cov = sum(float(r.get("mean_depth", 0.0)) * int(r.get("length", 0)) for r in bcr_abl_exon_rows) / total_exon_len
            kd_depth_str = f"{avg_exon_cov:.1f}x"
            kd_val = avg_exon_cov
            kd_status, kd_style = ("PASS", 5) if kd_val >= 100.0 else (("BORDERLINE", 6) if kd_val >= 20.0 else ("LOW", 4))
        else:
            kd_depth_str = "N/A"
            kd_status, kd_style = ("INFO", 0)
    elif bcr_abl_mutation_rows:
        kd_val = sum(int(r.get("depth", 0)) for r in bcr_abl_mutation_rows) / len(bcr_abl_mutation_rows)
        kd_depth_str = f"{kd_val:.1f}x"
        kd_status, kd_style = ("PASS", 5) if kd_val >= 100.0 else (("BORDERLINE", 6) if kd_val >= 20.0 else ("LOW", 4))
    else:
        kd_depth_str = "N/A"
        kd_status, kd_style = ("INFO", 0)

    try:
        w_val = float(wgs_cov.split()[0])
    except (ValueError, TypeError, IndexError):
        w_val = None

    if t_val is not None and w_val is not None and w_val > 0:
        fold_enr = t_val / w_val
        fold_enr_str = f"{fold_enr:.1f}x"
        fold_enr_status, fold_enr_style = (
            ("PASS", 5) if fold_enr >= 5.0 else (("BORDERLINE", 6) if fold_enr >= 2.0 else ("LOW", 4))
        )
    else:
        fold_enr_str = "N/A"
        fold_enr_status, fold_enr_style = ("INFO", 0)

    low_cov_20x = []
    if not has_bcr_abl:
        for r in target_rows:
            chrom_name = str(r.get("chrom", "")).strip()
            if not chrom_name or chrom_name.lower() == "total":
                continue
            try:
                m_cov = float(str(r.get("mean", 0)).split()[0])
            except (ValueError, TypeError):
                continue
            if m_cov < 20.0:
                low_cov_20x.append((chrom_name, m_cov))
        low_cov_20x.sort(key=lambda x: x[1])
        low_20x_count = len(low_cov_20x)
        low_cov_status, low_cov_style = (
            ("PASS", 5) if low_20x_count == 0 else (("REVIEW", 6) if low_20x_count <= 5 else ("ACTION REQUIRED", 4))
        )

    caller_display = detected_caller
    if has_bcr_abl and not snvs:
        caller_display = "Branford Pileup Caller"
    elif platform and "clairs" in caller_display.lower():
        caller_display = f"ClairS-TO ({platform})"
    elif not caller_display:
        caller_display = "ClairS-TO"

    target_display = "BCR-ABL1 Target Exons" if has_bcr_abl else (os.path.basename(design_bed) if design_bed else "Whole Genome (WGS)")

    demux_meta = parse_demux_info(demux_summary, sample, barcode)
    barcode_display = demux_meta["barcode"] or barcode or "N/A"
    demux_yield_display = (
        f"{demux_meta['reads']} reads ({demux_meta['fraction']})"
        if demux_meta["reads"]
        else ("Assigned" if barcode else "N/A")
    )

    # Fluffy-Style Overview Layout with rich sections and status color coding
    overview = [
        [("NAGINI NANOPORE ANALYSIS REPORT", 2), ("", 2), ("", 2), ("", 2)],
        [("", 0), ("", 0), ("", 0), ("", 0)],
        [("SAMPLE & PIPELINE METADATA", 3), ("", 3), ("", 3), ("", 3)],
        [("Sample ID", 7), (sample, 0), ("Sequencing Platform", 7), ("Oxford Nanopore (ONT)", 0)],
        [("Sample Type", 7), (f"Tumor ({sample_type})" if sample_type == "T" else f"Normal ({sample_type})", 0), ("Reference Genome", 7), (reference_build, 0)],
        [("Barcode / Demux ID", 7), (barcode_display, 0), ("Demux Pool Fraction", 7), (demux_yield_display, 0)],
        [("Small Variant Caller", 7), (caller_display, 0), ("Structural Variant Caller", 7), ("Sniffles2", 0)],
        [("Variant Annotation", 7), (annotation, 0), ("Target Region / Panel", 7), (target_display, 0)],
        [("Report Generated", 7), (datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"), 0), ("Pipeline Status", 7), ("COMPLETED", 5)],
        [("", 0), ("", 0), ("", 0), ("", 0)],
    ]

    if has_bcr_abl:
        overview.extend([
            [("QUALITY CONTROL & SEQUENCING METRICS", 3), ("", 3), ("VALUE", 3), ("STATUS", 3)],
            [("Total Yield [Gb]", 7), ("Total sequenced yield", 0), (yield_gb or "N/A", 0), ("PASS" if yield_gb else "INFO", 5 if yield_gb else 0)],
            [("Total Reads", 7), ("Total sequenced reads", 0), (total_reads or "N/A", 0), ("PASS" if total_reads else "INFO", 5 if total_reads else 0)],
            [("Mean Base Quality (Q)", 7), ("Average basecall Q-score", 0), (mean_q or "N/A", 0), (q_status, q_style)],
            [("Read Length N50 [bp]", 7), ("N50 read length", 0), (n50_bp or "N/A", 0), ("PASS" if n50_bp else "INFO", 5 if n50_bp else 0)],
            [("Median Read Length [bp]", 7), ("Median read length", 0), (median_len or "N/A", 0), ("PASS" if median_len else "INFO", 5 if median_len else 0)],
            [("BCR-ABL1 Fusion Reads", 7), ("Total fusion-spanning reads (deduplicated)", 0), (total_fusion_reads, 0), ("PASS" if total_fusion_reads > 0 else "INFO", 5 if total_fusion_reads > 0 else 0)],
            [("BCR-ABL1 Target Enrichment", 7), ("Fusion reads / total sequenced reads", 0), (target_enrich_str, 0), (enrich_status, enrich_style)],
            [("Kinase Domain Mean Depth", 7), ("Mean depth across ABL1 resistance sites", 0), (kd_depth_str, 0), (kd_status, kd_style)],
        ])
    else:
        overview.extend([
            [("QUALITY CONTROL & SEQUENCING METRICS", 3), ("", 3), ("VALUE", 3), ("STATUS", 3)],
            [("Total Yield [Gb]", 7), ("Total sequenced yield", 0), (yield_gb or "N/A", 0), ("PASS" if yield_gb else "INFO", 5 if yield_gb else 0)],
            [("Total Reads", 7), ("Total sequenced reads", 0), (total_reads or "N/A", 0), ("PASS" if total_reads else "INFO", 5 if total_reads else 0)],
            [("Mean Base Quality (Q)", 7), ("Average basecall Q-score", 0), (mean_q or "N/A", 0), (q_status, q_style)],
            [("Read Length N50 [bp]", 7), ("N50 read length", 0), (n50_bp or "N/A", 0), ("PASS" if n50_bp else "INFO", 5 if n50_bp else 0)],
            [("Median Read Length [bp]", 7), ("Median read length", 0), (median_len or "N/A", 0), ("PASS" if median_len else "INFO", 5 if median_len else 0)],
            [("Mean Target Coverage", 7), ("Mean depth in target BED", 0), (target_cov or "N/A", 0), (t_status, t_style)],
            [("Target Fold Enrichment", 7), ("Target / WGS coverage ratio", 0), (fold_enr_str, 0), (fold_enr_status, fold_enr_style)],
            [("Target Breadth (>=20x)", 7), ("% target bases >= 20x depth", 0), (pct_20x_str, 0), (p20_status, p20_style)],
            [("Target Breadth (>=30x)", 7), ("% target bases >= 30x depth", 0), (pct_30x_str, 0), ("PASS" if (pct_30x and pct_30x >= 0.85) else ("INFO" if pct_30x is not None else "N/A"), 5 if (pct_30x and pct_30x >= 0.85) else 0)],
            [("Mean WGS Coverage", 7), ("Mean whole-genome depth", 0), (wgs_cov or "N/A", 0), ("INFO", 0)],
        ])

    overview.append([("", 0), ("", 0), ("", 0), ("", 0)])
    overview.append([("VARIANT CALLING SUMMARY", 3), ("", 3), ("COUNT", 3), ("FLAG / ACTION", 3)])
    overview.append([("Total Small Variants", 7), ("SNVs and short indels", 0), (len(snvs), 0), ("PASS", 5)])
    overview.append([("PASS Small Variants", 7), ("High-confidence calls", 0), (pass_count, 0), ("PASS", 5)])
    overview.append([("ClinVar Pathogenic Variants", 7), ("Pathogenic / Likely Pathogenic", 0), (len(pathogenic_snvs), 0), ("ACTION REQUIRED" if pathogenic_snvs else "NONE", 4 if pathogenic_snvs else 5)])
    overview.append([("High Impact Variants", 7), ("Stop gained / frameshift / splice", 0), (len(high_impact_snvs), 0), ("REVIEW" if high_impact_snvs else "NONE", 6 if high_impact_snvs else 5)])
    if cnvs:
        overview.append([("Total Copy Number Alterations", 7), ("Genome-wide CNV segments (Rindels)", 0), (len(cnv_non_diploid), 0), ("CNV DETECTED" if cnv_non_diploid else "NONE", 6 if cnv_non_diploid else 5)])
        if target_cnvs:
            overview.append([("Target-Overlapping CNVs", 7), ("CNVs overlapping design BED", 0), (len(target_cnvs), 0), ("REVIEW" if target_cnvs else "NONE", 6 if target_cnvs else 5)])
    overview.append([("Structural Variants (Intra-chrom)", 7), ("DEL / DUP / INS / INV, all callers", 0), (len(svs_intra), 0), ("INFO", 0)])
    overview.append([("Target-Overlapping SVs", 7), ("SVs overlapping panel genes", 0), (len(target_svs), 0), ("REVIEW" if target_svs else "NONE", 6 if target_svs else 5)])
    overview.append([("Translocations / Breakends", 7), ("Inter-chromosomal rearrangements", 0), (len(translocations), 0), ("TRANSLOCATION DETECTED" if translocations else "NONE", 8 if translocations else 5)])
    if not has_bcr_abl:
        overview.append([("Low-Coverage Regions (<20x)", 7), ("Target regions below 20x depth", 0), (low_20x_count, 0), (low_cov_status, low_cov_style)])


    # Add BCR-ABL summary block to overview if BCR-ABL data present
    if has_bcr_abl:
        prim_iso = next((r["fusion_isoform"] for r in bcr_abl_fusion_rows if r.get("status") == "PRIMARY"), "None detected")
        det_muts = [r["mutation"] for r in bcr_abl_mutation_rows if str(r.get("status", "")).upper() == "DETECTED"]
        subcl_muts = [r["mutation"] for r in bcr_abl_mutation_rows if str(r.get("status", "")).upper() == "SUBCLONAL"]

        if det_muts:
            mut_summary = ", ".join(det_muts)
            if subcl_muts:
                mut_summary += f" (Subclonal: {', '.join(subcl_muts)})"
            mut_action = "ACTION REQUIRED"
            mut_action_style = 4
        elif subcl_muts:
            mut_summary = f"Subclonal: {', '.join(subcl_muts)}"
            mut_action = "SUBCLONAL"
            mut_action_style = 6
        else:
            mut_summary = "None detected"
            mut_action = "NONE"
            mut_action_style = 5

        cis_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in bcr_abl_phasing_rows if "CIS" in r.get("phase_state", "").upper()]
        phase_summary = f"CIS ({', '.join(cis_pairs)})" if cis_pairs else ("TRANS" if bcr_abl_phasing_rows else "N/A (<2 muts)")

        overview.append([("", 0), ("", 0), ("", 0), ("", 0)])
        overview.append([("BCR-ABL1 TARGETED ASSAY SUMMARY", 3), ("", 3), ("", 3), ("", 3)])
        overview.append([("Primary Fusion Isoform", 7), ("Major/Minor BCR breakpoint", 0), (prim_iso, 8), ("PASS" if prim_iso != "None detected" else "INFO", 5 if prim_iso != "None detected" else 0)])
        overview.append([("Resistance Mutations", 7), ("Branford ABL1 panel", 0), (mut_summary, 4 if det_muts else (6 if subcl_muts else 0)), (mut_action, mut_action_style)])
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

    if cnv_non_diploid:
        overview.append([("", 0), ("", 0), ("", 0), ("", 0)])
        overview.append([("COPY NUMBER ALTERATIONS (RINDELS)", 6), ("", 6), ("", 6), ("", 6), ("", 6), ("", 6), ("", 6)])
        overview.append([("Chr", 1), ("Start", 1), ("End", 1), ("Length [kb]", 1), ("Copy Number", 1), ("Type", 1), ("Cytoband", 1)])
        for c in cnv_non_diploid[:15]:
            cn_style = 4 if str(c.get("event_type", "")).upper() == "LOSS" else 6
            overview.append([
                (c.get("chrom", ""), 0),
                (c.get("start", 0), 0),
                (c.get("end", 0), 0),
                (c.get("length_kb", 0), 0),
                (c.get("copy_number", "") or c.get("state", ""), cn_style),
                (c.get("event_type", ""), cn_style),
                (c.get("cytoband", "") or c.get("target_genes", ""), 0),
            ])

    if translocations:
        overview.append([("", 0), ("", 0), ("", 0), ("", 0)])
        overview.append([("TRANSLOCATIONS & REARRANGEMENTS", 8), ("", 8), ("", 8), ("", 8), ("", 8), ("", 8), ("", 8)])
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

    if low_cov_20x:
        overview.append([("", 0), ("", 0), ("", 0), ("", 0)])
        overview.append([("LOW-COVERAGE TARGET REGIONS (<20x)", 6), ("", 6), ("", 6), ("", 6)])
        overview.append([("Region / Target", 1), ("Mean Depth", 1), ("Status", 1), ("Description", 1)])
        for reg, d in low_cov_20x[:15]:
            st = "SEVERELY LOW (<10x)" if d < 10.0 else "SUBOPTIMAL (<20x)"
            st_style = 4 if d < 10.0 else 6
            overview.append([
                (reg, 0),
                (f"{d:.1f}x", st_style),
                (st, st_style),
                ("Target depth below 20x threshold", 0),
            ])

    sheets = [
        ("Overview", overview, LAYOUT),
    ]
    if has_bcr_abl:
        bcr_layout = build_bcr_abl_sheet(
            sample, bcr_abl_fusion_rows, bcr_abl_mutation_rows, bcr_abl_phasing_rows, bcr_abl_exon_rows, plot_png=bcr_abl_plot_png
        )
        sheets.append(("BCR-ABL1", bcr_layout, LAYOUT))

    sheets.extend([
        ("SNVs and indels", table(snv_columns, snvs)),
        ("CNVs", table(cnv_columns, cnvs)),
        ("SVs", table(sv_columns, svs_intra)),
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
        rindels_sv_tsv=_first_or_str(getattr(snakemake.input, "rindels_sv", "")),
        rindels_translocations_tsv=_first_or_str(
            getattr(snakemake.input, "rindels_translocations", "")
        ),
        wgs_coverage=_first_or_str(getattr(snakemake.input, "wgs_coverage", "")),
        target_coverage=_first_or_str(getattr(snakemake.input, "target_coverage", "")),
        cramino=_first_or_str(getattr(snakemake.input, "cramino", "")),
        cytoband=str(getattr(snakemake.params, "cytoband", "") or ""),
        caller=str(getattr(snakemake.params, "small_variant_caller", "") or ""),
        design_bed=str(getattr(snakemake.params, "design_bed", "") or ""),
        platform=str(getattr(snakemake.params, "clairs_to_platform", "") or ""),
        target_dist=_first_or_str(getattr(snakemake.input, "target_dist", "")),
        demux_summary=_first_or_str(getattr(snakemake.input, "demux_summary", "")),
        barcode=str(getattr(snakemake.params, "barcode", "") or ""),
        bcr_abl_fusion=_first_or_str(getattr(snakemake.input, "bcr_abl_fusion", "")),
        bcr_abl_mutations=_first_or_str(getattr(snakemake.input, "bcr_abl_mutations", "")),
        bcr_abl_phasing=_first_or_str(getattr(snakemake.input, "bcr_abl_phasing", "")),
        bcr_abl_exon_coverage=_first_or_str(getattr(snakemake.input, "bcr_abl_exon_coverage", "")),
        bcr_abl_plot_png=_first_or_str(getattr(snakemake.input, "bcr_abl_plot_png", "")),
        cnv_tsv=_first_or_str(getattr(snakemake.input, "cnv", "")),
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
    parser.add_argument("--target-dist", default="")
    parser.add_argument("--demux-summary", default="")
    parser.add_argument("--barcode", default="")
    parser.add_argument("--bcr-abl-fusion", default="")
    parser.add_argument("--bcr-abl-mutations", default="")
    parser.add_argument("--bcr-abl-phasing", default="")
    parser.add_argument("--bcr-abl-exon-coverage", default="")
    parser.add_argument("--bcr-abl-plot-png", default="")
    parser.add_argument("--cnv-tsv", default="")
    args = parser.parse_args()
    generate_report(
        args.output, args.sample, args.sample_type, args.reference_build, args.annotation,
        args.snv_vcf, args.sv_tsv, args.translocations_tsv, args.wgs_coverage,
        args.target_coverage, args.cramino, args.cytoband, args.caller,
        args.design_bed, args.platform,
        target_dist=args.target_dist,
        demux_summary=args.demux_summary,
        barcode=args.barcode,
        bcr_abl_fusion=args.bcr_abl_fusion,
        bcr_abl_mutations=args.bcr_abl_mutations,
        bcr_abl_phasing=args.bcr_abl_phasing,
        bcr_abl_exon_coverage=args.bcr_abl_exon_coverage,
        bcr_abl_plot_png=args.bcr_abl_plot_png,
        cnv_tsv=args.cnv_tsv,
    )


if "snakemake" in globals():
    run_from_snakemake(snakemake)
elif __name__ == "__main__":
    main()


