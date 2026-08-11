import gzip
import importlib.util
import tempfile
import unittest
import zipfile
from pathlib import Path


SCRIPT = Path(__file__).with_name("nagini_report.py")
SPEC = importlib.util.spec_from_file_location("nagini_report", SCRIPT)
REPORT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(REPORT)


class NaginiReportTest(unittest.TestCase):
    def test_parse_vep_vcf_and_write_workbook(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
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

            columns, records = REPORT.parse_small_variant_vcf(vcf)
            self.assertIn("Gene", columns)
            self.assertEqual(records[0]["Gene"], "GENE1")
            self.assertEqual(records[0]["VAF"], "0.25")
            self.assertEqual(records[0]["ClinVar_Significance"], "Pathogenic")
            self.assertEqual(records[0]["ClinVar_Disease"], "Disease_one")

            output = root / "report.xlsx"
            REPORT.write_xlsx(output, [("Overview", [["Property", "Value"], ["Sample", "example"]])])
            self.assertTrue(zipfile.is_zipfile(output))
            with zipfile.ZipFile(output) as workbook:
                self.assertIn("xl/worksheets/sheet1.xml", workbook.namelist())
                sheet = workbook.read("xl/worksheets/sheet1.xml").decode()
                self.assertIn("example", sheet)

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
            cramino.write_text("Reads\t10\nYield [Gb]: 0.01\n")

            cytoband = root / "cytoBand.txt"
            cytoband.write_text(
                "chr1\t0\t100\tp36.33\tgneg\n"
                "chr2\t0\t100\tp25.3\tgneg\n"
            )

            full_output = root / "full_report.xlsx"
            REPORT.generate_report(
                str(full_output), "example", "N", "hg38", "VEP", str(vcf),
                str(sv), str(translocations), str(coverage), str(coverage), str(cramino),
                str(cytoband),
            )
            with zipfile.ZipFile(full_output) as workbook:
                workbook_xml = workbook.read("xl/workbook.xml").decode()
                self.assertEqual(workbook_xml.count("<sheet "), 7)
                self.assertIn("Translocations", workbook_xml)
                overview = workbook.read("xl/worksheets/sheet1.xml").decode()
                self.assertIn("ClinVar pathogenic / likely pathogenic SNVs", overview)
                sv_sheet = workbook.read("xl/worksheets/sheet3.xml").decode()
                self.assertIn("Cytogenetic_Notation", sv_sheet)


if __name__ == "__main__":
    unittest.main()
