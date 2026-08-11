"""Create a dependency-free XLSX summary for a quicken analysis.

The writer uses only Python's standard library.  This keeps report generation
available in the same offline runtime as the other pipeline helper scripts.
"""

import argparse
import bisect
import collections
import csv
import gzip
import os
import re
import zipfile
from datetime import datetime, timezone
from xml.sax.saxutils import escape


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


def parse_small_variant_vcf(path):
    columns = [
        "Chromosome", "Position", "Reference", "Alternate", "Filter", "Quality",
        "VAF", "Depth", "Allele depth", "Gene", "Consequence", "Impact",
        "HGVSc", "HGVSp", "Existing variation", "ClinVar_Significance",
        "ClinVar_Disease", "ClinVar_Review",
    ]
    records = []
    csq_fields = []
    ann_fields = []
    with open_text(path) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                csq_fields = _annotation_fields(line, "CSQ")
                continue
            if line.startswith("##INFO=<ID=ANN"):
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

            records.append({
                "Chromosome": chrom,
                "Position": int(pos),
                "Reference": ref,
                "Alternate": alt,
                "Filter": filter_value,
                "Quality": qual,
                "VAF": sample_values.get("AF", sample_values.get("VAF", "")),
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
    return columns, records


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


def table(columns, records):
    return [columns] + [[record.get(column, "") for column in columns] for record in records]


def excel_column(index):
    name = ""
    while index:
        index, remainder = divmod(index - 1, 26)
        name = chr(65 + remainder) + name
    return name


def _cell_xml(reference, value, style=0):
    style_attr = f' s="{style}"' if style else ""
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return f'<c r="{reference}"{style_attr}><v>{value}</v></c>'
    text = re.sub(r"[\x00-\x08\x0B\x0C\x0E-\x1F]", "", "" if value is None else str(value))
    text = escape(text)
    preserve = ' xml:space="preserve"' if text != text.strip() else ""
    return f'<c r="{reference}" t="inlineStr"{style_attr}><is><t{preserve}>{text}</t></is></c>'


def worksheet_xml(rows):
    widths = []
    for column_index in range(max((len(row) for row in rows), default=0)):
        width = max((len(str(row[column_index])) for row in rows if column_index < len(row)), default=8)
        widths.append(min(max(width + 2, 10), 45))
    cols = "".join(
        f'<col min="{index}" max="{index}" width="{width}" customWidth="1"/>'
        for index, width in enumerate(widths, 1)
    )
    row_xml = []
    for row_index, row in enumerate(rows, 1):
        cells = "".join(
            _cell_xml(f"{excel_column(column_index)}{row_index}", value, 1 if row_index == 1 else 0)
            for column_index, value in enumerate(row, 1)
        )
        row_xml.append(f'<row r="{row_index}">{cells}</row>')
    max_col = excel_column(max((len(row) for row in rows), default=1))
    auto_filter = f'<autoFilter ref="A1:{max_col}{len(rows)}"/>' if len(rows) > 1 else ""
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        '<sheetViews><sheetView workbookViewId="0"><pane ySplit="1" topLeftCell="A2" '
        'activePane="bottomLeft" state="frozen"/></sheetView></sheetViews>'
        f'<cols>{cols}</cols><sheetData>{"".join(row_xml)}</sheetData>{auto_filter}</worksheet>'
    )


def write_xlsx(path, sheets):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    sheet_entries = "".join(
        f'<sheet name="{escape(name)}" sheetId="{index}" r:id="rId{index}"/>'
        for index, (name, _rows) in enumerate(sheets, 1)
    )
    relationships = "".join(
        f'<Relationship Id="rId{index}" '
        'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" '
        f'Target="worksheets/sheet{index}.xml"/>'
        for index in range(1, len(sheets) + 1)
    )
    overrides = "".join(
        f'<Override PartName="/xl/worksheets/sheet{index}.xml" '
        'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        for index in range(1, len(sheets) + 1)
    )
    timestamp = datetime.now(timezone.utc).replace(microsecond=0).isoformat().replace("+00:00", "Z")

    parts = {
        "[Content_Types].xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
            '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
            '<Default Extension="xml" ContentType="application/xml"/>'
            '<Override PartName="/xl/workbook.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
            '<Override PartName="/xl/styles.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
            '<Override PartName="/docProps/core.xml" '
            'ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>'
            '<Override PartName="/docProps/app.xml" '
            'ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>'
            f"{overrides}</Types>"
        ),
        "_rels/.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            '<Relationship Id="rId1" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" '
            'Target="xl/workbook.xml"/>'
            '<Relationship Id="rId2" '
            'Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" '
            'Target="docProps/core.xml"/>'
            '<Relationship Id="rId3" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/extended-properties" '
            'Target="docProps/app.xml"/>'
            '</Relationships>'
        ),
        "xl/workbook.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" '
            'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
            f"<sheets>{sheet_entries}</sheets></workbook>"
        ),
        "xl/_rels/workbook.xml.rels": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
            f'{relationships}<Relationship Id="rId{len(sheets) + 1}" '
            'Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" '
            'Target="styles.xml"/></Relationships>'
        ),
        "xl/styles.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
            '<fonts count="2"><font><sz val="11"/><name val="Calibri"/></font>'
            '<font><b/><color rgb="FFFFFFFF"/><sz val="11"/><name val="Calibri"/></font></fonts>'
            '<fills count="3"><fill><patternFill patternType="none"/></fill>'
            '<fill><patternFill patternType="gray125"/></fill>'
            '<fill><patternFill patternType="solid"><fgColor rgb="FF1F4E78"/>'
            '<bgColor indexed="64"/></patternFill></fill></fills>'
            '<borders count="1"><border><left/><right/><top/><bottom/><diagonal/></border></borders>'
            '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
            '<cellXfs count="2"><xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
            '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/></cellXfs>'
            '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
            '</styleSheet>'
        ),
        "docProps/core.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
            'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
            'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
            '<dc:title>quicken analysis report</dc:title><dc:creator>quicken</dc:creator>'
            f'<dcterms:created xsi:type="dcterms:W3CDTF">{timestamp}</dcterms:created></cp:coreProperties>'
        ),
        "docProps/app.xml": (
            '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
            '<Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties">'
            '<Application>quicken</Application></Properties>'
        ),
    }
    for index, (_name, rows) in enumerate(sheets, 1):
        parts[f"xl/worksheets/sheet{index}.xml"] = worksheet_xml(rows)

    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as workbook:
        for name, content in parts.items():
            workbook.writestr(name, content)


def generate_report(output, sample, sample_type, reference_build, annotation, snv_vcf,
                    sv_tsv, translocations_tsv, wgs_coverage, target_coverage, cramino,
                    cytoband=""):
    snv_columns, snvs = parse_small_variant_vcf(snv_vcf)
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

    pass_count = sum(record["Filter"] in {"PASS", "."} for record in snvs)
    clinvar_pathogenic_count = sum(
        _is_pathogenic_clinvar(record.get("ClinVar_Significance", ""))
        for record in snvs
    )
    overview = [
        ["Property", "Value"],
        ["Sample", sample],
        ["Type", sample_type],
        ["Platform", "Oxford Nanopore (ONT)"],
        ["Reference build", reference_build],
        ["Small-variant annotation", annotation],
        ["Small variants (all records)", len(snvs)],
        ["Small variants (PASS)", pass_count],
        ["ClinVar pathogenic / likely pathogenic SNVs", clinvar_pathogenic_count],
        ["Structural variants", len(svs)],
        ["Translocations / BNDs", len(translocations)],
        ["Interpretation", "Research-use summary; not clinically validated"],
    ]
    sheets = [
        ("Overview", overview),
        ("SNVs and indels", table(snv_columns, snvs)),
        ("SVs", table(sv_columns, svs)),
        ("Translocations", table(translocation_columns, translocations)),
        ("WGS coverage", table(wgs_columns, wgs_rows)),
        ("Target coverage", table(target_columns, target_rows)),
        ("Cramino QC", table(cramino_columns, cramino_rows)),
    ]
    write_xlsx(output, sheets)


def run_from_snakemake(snakemake):
    generate_report(
        output=str(snakemake.output.xlsx),
        sample=str(snakemake.wildcards.sample),
        sample_type=str(snakemake.wildcards.type),
        reference_build=str(snakemake.params.reference_build),
        annotation=str(snakemake.params.snv_annotation),
        snv_vcf=str(snakemake.input.snv_vcf),
        sv_tsv=str(snakemake.input.sv),
        translocations_tsv=str(snakemake.input.translocations),
        wgs_coverage=str(snakemake.input.wgs_coverage),
        target_coverage=str(snakemake.input.target_coverage),
        cramino=str(snakemake.input.cramino),
        cytoband=str(snakemake.params.cytoband or ""),
    )


def main():
    parser = argparse.ArgumentParser(description="Create a quicken XLSX report")
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
    args = parser.parse_args()
    generate_report(
        args.output, args.sample, args.sample_type, args.reference_build, args.annotation,
        args.snv_vcf, args.sv_tsv, args.translocations_tsv, args.wgs_coverage,
        args.target_coverage, args.cramino, args.cytoband,
    )


if "snakemake" in globals():
    run_from_snakemake(snakemake)
elif __name__ == "__main__":
    main()
