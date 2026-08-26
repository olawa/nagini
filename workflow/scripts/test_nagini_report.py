import gzip
import importlib.util
import tempfile
import unittest
import zipfile
from pathlib import Path
from xml.etree import ElementTree


SCRIPT = Path(__file__).with_name("nagini_report.py")
SPEC = importlib.util.spec_from_file_location("nagini_report", SCRIPT)
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)

NS = "{http://schemas.openxmlformats.org/spreadsheetml/2006/main}"
REL_NS = "{http://schemas.openxmlformats.org/officeDocument/2006/relationships}"


def read_workbook(path):
    """Read a workbook back as {sheet name: [[cell text, ...], ...]}.

    Asserting on rendered content rather than on the writer's XML keeps these
    tests about the report instead of about the library that emits it.
    """

    with zipfile.ZipFile(path) as workbook:
        names = workbook.namelist()
        shared = []
        if "xl/sharedStrings.xml" in names:
            root = ElementTree.fromstring(workbook.read("xl/sharedStrings.xml"))
            for item in root.findall(f"{NS}si"):
                shared.append("".join(node.text or "" for node in item.iter(f"{NS}t")))

        rels = {}
        root = ElementTree.fromstring(workbook.read("xl/_rels/workbook.xml.rels"))
        for relationship in root:
            rels[relationship.get("Id")] = relationship.get("Target").lstrip("/")

        sheets = {}
        root = ElementTree.fromstring(workbook.read("xl/workbook.xml"))
        for sheet in root.find(f"{NS}sheets"):
            target = rels[sheet.get(f"{REL_NS}id")]
            member = target if target.startswith("xl/") else f"xl/{target}"
            sheet_root = ElementTree.fromstring(workbook.read(member))
            rows = []
            for row in sheet_root.iter(f"{NS}row"):
                cells = []
                for cell in row.findall(f"{NS}c"):
                    value = cell.find(f"{NS}v")
                    text = "" if value is None else (value.text or "")
                    if cell.get("t") == "s":
                        text = shared[int(text)]
                    elif cell.get("t") == "inlineStr":
                        inline = cell.find(f"{NS}is/{NS}t")
                        text = "" if inline is None else (inline.text or "")
                    cells.append(text)
                rows.append(cells)
            sheets[sheet.get("name")] = rows
        return sheets


def sheet_xml(path, sheet_name):
    """Raw XML for one sheet, for the few structural assertions that need it."""

    with zipfile.ZipFile(path) as workbook:
        rels = {}
        root = ElementTree.fromstring(workbook.read("xl/_rels/workbook.xml.rels"))
        for relationship in root:
            rels[relationship.get("Id")] = relationship.get("Target").lstrip("/")
        root = ElementTree.fromstring(workbook.read("xl/workbook.xml"))
        for sheet in root.find(f"{NS}sheets"):
            if sheet.get("name") == sheet_name:
                target = rels[sheet.get(f"{REL_NS}id")]
                member = target if target.startswith("xl/") else f"xl/{target}"
                return workbook.read(member).decode()
    raise AssertionError(f"No sheet named {sheet_name}")


def cell_info(path, sheet_name, reference):
    """Return (type, style index, value) for one cell reference."""

    root = ElementTree.fromstring(sheet_xml(path, sheet_name))
    for cell in root.iter(f"{NS}c"):
        if cell.get("r") == reference:
            value = cell.find(f"{NS}v")
            return cell.get("t"), cell.get("s"), (None if value is None else value.text)
    raise AssertionError(f"No cell {reference} on {sheet_name}")


def write_test_inputs(root):
    vcf = root / "calls.vcf.gz"
    with gzip.open(vcf, "wt") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write(
            '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations. '
            'Format: Allele|Consequence|IMPACT|SYMBOL|HGVSc|HGVSp|Existing_variation|'
            'ClinVar_CLNSIG|ClinVar_CLNREVSTAT|ClinVar_CLNDN">\n'
        )
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\texample\n")
        handle.write(
            "chr1\t10\t.\tA\tT\t60\tPASS\t"
            "CSQ=T|missense_variant|MODERATE|GENE1|c.10A>T|p.Lys4Met|rs1|"
            "Pathogenic|criteria_provided,_multiple_submitters|Disease_one;DP=42\t"
            "GT:AF:DP:AD\t0/1:0.25:40:30,10\n"
        )
    sv = root / "sv.tsv"
    sv.write_text(
        "sample\ttype\tsvtype\tchrom\tpos\tend\tchrom2\ttranslocation\n"
        "example\tN\tDEL\tchr1\t20\t30\t\tFalse\n"
    )
    translocations = root / "translocations.tsv"
    translocations.write_text(
        "sample\ttype\tsvtype\tchrom\tpos\tend\tchrom2\ttranslocation\n"
        "example\tN\tBND\tchr1\t40\t40\tchr2\tTrue\n"
    )
    coverage = root / "coverage.tsv"
    coverage.write_text("chrom\tlength\tbases\tmean\nchr1\t100\t100\t12.5\ntotal\t100\t100\t12.5\n")
    cramino = root / "cramino.txt"
    cramino.write_text("Reads\t10\nYield [Gb]: 0.01\nMean quality: 14.2\n")
    cytoband = root / "cytoBand.txt"
    cytoband.write_text("chr1\t0\t100\tp36.33\tgneg\nchr2\t0\t100\tp25.3\tgneg\n")
    return vcf, sv, translocations, coverage, cramino, cytoband


class VcfParsingTest(unittest.TestCase):
    def test_vep_annotations_and_format_fields_are_extracted(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            vcf, *_ = write_test_inputs(Path(tmpdir))
            columns, records, caller = REPORT.parse_small_variant_vcf(vcf)
        self.assertIn("Gene", columns)
        self.assertIn("Caller", columns)
        self.assertEqual(records[0]["Gene"], "GENE1")
        self.assertEqual(records[0]["VAF"], "0.250")
        self.assertEqual(records[0]["ClinVar_Significance"], "Pathogenic")
        self.assertEqual(records[0]["ClinVar_Disease"], "Disease_one")


class WriterTest(unittest.TestCase):
    def test_values_survive_a_round_trip(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "report.xlsx"
            REPORT.write_xlsx(
                output, [("Overview", [["Property", "Value"], ["Sample", "example"]])]
            )
            self.assertTrue(zipfile.is_zipfile(output))
            sheets = read_workbook(output)
        self.assertEqual(sheets["Overview"], [["Property", "Value"], ["Sample", "example"]])

    def test_numbers_are_written_as_numbers(self):
        """A count must sort and sum in Excel, not sit there as text."""

        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "n.xlsx"
            REPORT.write_xlsx(output, [("Counts", [["Metric", "Value"], ["Variants", 42]])])
            kind, _style, value = cell_info(output, "Counts", "B2")
            text_kind, _, _ = cell_info(output, "Counts", "A2")
        self.assertIsNone(kind, "numeric cells carry no type attribute")
        self.assertEqual(value, "42")
        self.assertEqual(text_kind, "s")

    def test_different_style_ids_produce_different_formats(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "s.xlsx"
            REPORT.write_xlsx(
                output,
                [(
                    "S",
                    [[("Title", 2)], [("Alert", 4)], [("Good", 5)], [("Alert again", 4)]],
                    REPORT.LAYOUT,
                )],
            )
            title = cell_info(output, "S", "A1")[1]
            alert = cell_info(output, "S", "A2")[1]
            good = cell_info(output, "S", "A3")[1]
            alert_again = cell_info(output, "S", "A4")[1]
        self.assertEqual(len({title, alert, good}), 3, "each style id needs its own format")
        self.assertEqual(alert, alert_again, "the same style id must reuse one format")

    def test_status_words_are_coloured_without_an_explicit_style(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "auto.xlsx"
            REPORT.write_xlsx(
                output, [("S", [["Header"], ["PASS"], ["REVIEW"], ["something else"]])]
            )
            good = cell_info(output, "S", "A2")[1]
            review = cell_info(output, "S", "A3")[1]
            plain = cell_info(output, "S", "A4")[1]
        self.assertNotEqual(good, plain)
        self.assertNotEqual(review, plain)
        self.assertNotEqual(good, review)

    def test_layout_sheets_get_no_autofilter(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "l.xlsx"
            REPORT.write_xlsx(
                output,
                [
                    ("Layout", [["a", "b"], ["c", "d"]], REPORT.LAYOUT),
                    ("Table", [["a", "b"], ["c", "d"]]),
                ],
            )
            layout = sheet_xml(output, "Layout")
            table = sheet_xml(output, "Table")
        self.assertNotIn("autoFilter", layout)
        self.assertNotIn("frozen", layout)
        self.assertIn("autoFilter", table)
        self.assertIn("frozen", table)

    def test_sheet_names_are_made_legal_and_unique(self):
        self.assertEqual(REPORT.safe_sheet_name("SNVs/indels"), "SNVs-indels")
        self.assertEqual(len(REPORT.safe_sheet_name("x" * 40)), 31)
        used = set()
        first = REPORT.safe_sheet_name("Report", used)
        second = REPORT.safe_sheet_name("Report", used)
        self.assertEqual(first, "Report")
        self.assertNotEqual(first, second)

    def test_control_characters_are_stripped(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            output = Path(tmpdir) / "c.xlsx"
            REPORT.write_xlsx(output, [("S", [["Header"], ["a\x07b"]])])
            sheets = read_workbook(output)
        self.assertEqual(sheets["S"][1][0], "ab")


class FullReportTest(unittest.TestCase):
    def _build(self, tmpdir):
        root = Path(tmpdir)
        vcf, sv, translocations, coverage, cramino, cytoband = write_test_inputs(root)
        output = root / "full_report.xlsx"
        REPORT.generate_report(
            str(output), "example", "N", "hg38", "VEP", str(vcf),
            str(sv), str(translocations), str(coverage), str(coverage), str(cramino),
            str(cytoband), caller="ClairS-TO", design_bed="/data/panel.bed",
            platform="ont_r10_dorado_hac_4khz",
        )
        return output

    def test_every_sheet_is_present(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sheets = read_workbook(self._build(tmpdir))
        self.assertEqual(
            list(sheets),
            [
                "Overview", "SNVs and indels", "SVs", "Translocations",
                "WGS coverage", "Target coverage", "Cramino QC",
            ],
        )

    def test_overview_reports_metadata_and_findings(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sheets = read_workbook(self._build(tmpdir))
        flat = ["\t".join(row) for row in sheets["Overview"]]
        text = "\n".join(flat)
        self.assertIn("NAGINI NANOPORE ANALYSIS REPORT", text)
        self.assertIn("ClairS-TO (ont_r10_dorado_hac_4khz)", text)
        self.assertIn("panel.bed", text)
        self.assertIn("PATHOGENIC", text)
        self.assertIn("TRANSLOCATION DETECTED", text)

    def test_structural_variants_are_annotated_with_cytobands(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sheets = read_workbook(self._build(tmpdir))
        self.assertIn("Cytogenetic_Notation", sheets["SVs"][0])
        self.assertIn("Start_Band", sheets["SVs"][0])

    def test_variant_row_reaches_the_snv_sheet(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            sheets = read_workbook(self._build(tmpdir))
        header = sheets["SNVs and indels"][0]
        row = dict(zip(header, sheets["SNVs and indels"][1]))
        self.assertEqual(row["Gene"], "GENE1")
        self.assertEqual(row["Caller"], "ClairS-TO")
        self.assertEqual(row["ClinVar_Significance"], "Pathogenic")


if __name__ == "__main__":
    unittest.main()
