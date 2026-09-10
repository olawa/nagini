import tempfile
import unittest
from pathlib import Path

try:
    import xlsxwriter
    import nagini_summary_report as SUMMARY_MODULE
    from nagini_summary_report import build_summary_report
except ImportError:
    xlsxwriter = None
    SUMMARY_MODULE = None
    build_summary_report = None


@unittest.skipIf(xlsxwriter is None, "xlsxwriter not available")
class SummaryReportTest(unittest.TestCase):

    def test_build_summary_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            out_xlsx = str(tmp / "summary.xlsx")

            # Create test input files for sample1
            cramino = tmp / "sample1.cramino.txt"
            cramino.write_text("Yield [Gb]\t0.15\nNumber of reads\t50000\nMean quality\t18.5\nRead length N50\t1400\nMedian read length\t700\n")

            fusion_tsv = tmp / "sample1.fusion.tsv"
            fusion_tsv.write_text("sample\tfusion_isoform\tread_count\tfraction_percent\tstatus\nsample1\te14a2\t1200\t80.0\tPRIMARY\nsample1\te13a2\t300\t20.0\tSECONDARY\n")

            mut_tsv = tmp / "sample1.mut.tsv"
            mut_tsv.write_text("sample\tmutation\tchrom\tpos\tref\talt\tstatus\tdepth\talt_reads\tvaf_percent\ttki_resistance\nsample1\tT315I\tchr9\t133748283\tC\tT\tDETECTED\t5000\t2500\t50.0\tPan-resistant\n")

            phase_tsv = tmp / "sample1.phase.tsv"
            phase_tsv.write_text("sample\tmutation_1\tmutation_2\tphase_state\tcis_reads\ttrans_reads\tpct_cis\nsample1\tT315I\tY253H\tCIS\t1000\t10\t99.0\n")

            sample_records = [
                {
                    "sample": "sample1",
                    "type": "N",
                    "barcode": "barcode01",
                    "assay": "bcr_abl",
                    "cramino": str(cramino),
                    "wgs_coverage": "",
                    "snv_vcf": "",
                    "bcr_abl_fusion": str(fusion_tsv),
                    "bcr_abl_mutations": str(mut_tsv),
                    "bcr_abl_phasing": str(phase_tsv),
                    "bcr_abl_exon_coverage": "",
                }
            ]

            build_summary_report(
                output_path=out_xlsx,
                sequence_id="TestRun",
                reference_build="hg19",
                sample_records=sample_records,
            )

            self.assertTrue(Path(out_xlsx).exists())
            self.assertGreater(Path(out_xlsx).stat().st_size, 1000)

            import zipfile
            with zipfile.ZipFile(out_xlsx, "r") as z:
                wb_xml = z.read("xl/workbook.xml").decode("utf-8")
                self.assertIn("Overview", wb_xml)
                self.assertIn("BCR-ABL1 Overview", wb_xml)
                self.assertIn("BCR-ABL (sample1)", wb_xml)
            # A bcr_abl run keeps the columns an adaptive run drops.
            text = _shared_strings(out_xlsx)
            self.assertIn("BCR-ABL1 Fusion Reads", text)
            self.assertIn("Mean Target / KD Depth", text)


    def test_a_run_without_small_variant_calling_says_so(self):
        """"Not assessed" and "none found" are different clinical statements.

        An assay with snv_calling off produces no VCF, and reporting that as
        "No high-confidence PASS variants detected" reads as a negative
        finding rather than as a question the run did not ask.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            cramino = tmp / "s.cramino.txt"
            cramino.write_text("Yield [Gb]\t0.15\nNumber of reads\t50000\n")
            record = {
                "sample": "gustaf1", "type": "T", "barcode": "barcode01",
                "assay": "adaptive_cnv_sv", "cramino": str(cramino),
                "wgs_coverage": "", "snv_vcf": "", "bcr_abl_fusion": "",
                "bcr_abl_mutations": "", "bcr_abl_phasing": "",
                "bcr_abl_exon_coverage": "",
            }
            out = str(tmp / "summary.xlsx")
            build_summary_report(
                output_path=out, sequence_id="T2T", reference_build="hs1",
                sample_records=[record],
            )
            text = _shared_strings(out)

        self.assertIn("was not part of this run's assay", text)
        self.assertNotIn("No high-confidence PASS variants detected", text)
        # An adaptive run has no fusion, isoform, resistance or phasing data.
        for column in (
            "BCR-ABL1 Fusion Reads",
            "BCR-ABL1 Enrichment",
            "Primary Isoform",
            "Resistance Mutations",
            "Compound Phasing",
        ):
            self.assertNotIn(column, text)
        self.assertIn("Mean Target Depth", text)


class StackedTablesTest(unittest.TestCase):
    def test_the_sample_column_leads_and_columns_are_unioned(self):
        """Five samples' calls should be readable in one place.

        Callers and versions can write different column sets, so the union is
        taken and blanks are left rather than shifting a row's fields.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            first = tmp / "a.tsv"
            first.write_text(
                "sample\ttype\tcaller\tsvtype\tchrom\n"
                "gustaf1\tT\tSniffles2\tDEL\tchr9\n"
            )
            # A second sample whose table carries an extra column.
            second = tmp / "b.tsv"
            second.write_text(
                "sample\ttype\tcaller\tsvtype\tchrom\tgenes\n"
                "gustaf2\tT\tRindels\tBND\tchr21\tRUNX1\n"
            )
            records = [
                {"sample": "gustaf1", "type": "T", "sv_tsv": str(first)},
                {"sample": "gustaf2", "type": "T", "sv_tsv": str(second)},
            ]
            rows = SUMMARY_MODULE.stack_sample_tables(records, ("sv_tsv",))

        header = [cell[0] for cell in rows[0]]
        self.assertEqual(header[0], "sample")
        self.assertIn("genes", header)
        body = [[cell[0] for cell in row] for row in rows[1:]]
        self.assertEqual([row[0] for row in body], ["gustaf1", "gustaf2"])
        # The sample without a genes column gets a blank, not a shifted row.
        self.assertEqual(body[0][header.index("genes")], "")
        self.assertEqual(body[1][header.index("genes")], "RUNX1")

    def test_no_sheet_when_no_sample_produced_the_table(self):
        self.assertEqual(SUMMARY_MODULE.stack_sample_tables([{"sample": "x"}], ("cnv_tsv",)), [])


def _shared_strings(xlsx_path):
    """Every string in the workbook, which is where cell text lives."""

    import zipfile

    with zipfile.ZipFile(xlsx_path, "r") as archive:
        return archive.read("xl/sharedStrings.xml").decode("utf-8")


if __name__ == "__main__":
    unittest.main()

