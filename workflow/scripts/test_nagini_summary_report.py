import tempfile
import unittest
from pathlib import Path

try:
    import xlsxwriter
    from nagini_summary_report import build_summary_report
except ImportError:
    xlsxwriter = None
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


def _shared_strings(xlsx_path):
    """Every string in the workbook, which is where cell text lives."""

    import zipfile

    with zipfile.ZipFile(xlsx_path, "r") as archive:
        return archive.read("xl/sharedStrings.xml").decode("utf-8")


if __name__ == "__main__":
    unittest.main()

