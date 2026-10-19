import tempfile
import unittest
from pathlib import Path

try:
    import openpyxl
except ImportError:
    openpyxl = None

import nagini_report as report
import nagini_summary_report as summary


@unittest.skipIf(openpyxl is None, "openpyxl is required to inspect generated workbooks")
class BcrAblReportRegressionTests(unittest.TestCase):
    def _write(self, root, name, content):
        path = root / name
        path.write_text(content)
        return str(path)

    @staticmethod
    def _sheet_text(path, sheet):
        book = openpyxl.load_workbook(path, read_only=True, data_only=True)
        try:
            return "\n".join("\t".join("" if v is None else str(v) for v in row) for row in book[sheet].iter_rows(values_only=True))
        finally:
            book.close()

    def test_sample_report_preserves_uncertainty_and_audit_findings(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            muts = self._write(root, "mut.tsv", "mutation\taa_change\tconsequence\tresistance_class\tstatus\tchrom\tpos\tref\talt\tdepth\tad_alt\tvaf_percent\nL387L\tp.Leu387Leu\tSYNONYMOUS\tSYNONYMOUS\tDETECTED\tchr9\t100\tA\tG\t100\t20\t20\nT315I\tp.Thr315Ile\tMISSENSE\tCURATED_RESISTANCE\tSUBCLONAL\tchr9\t200\tC\tT\t1000\t5\t0.5\n")
            phase = self._write(root, "phase.tsv", "mutation_1\tmutation_2\tphase_state\ntestA\ttestB\tLOW_SPANNING_DEPTH\ntestC\ttestD\tCIS (Compound)\n")
            sites = self._write(root, "all.tsv", "chrom\tpos\tstatus\tdepth\tfilter_reason\tpanel_sites_expected\nchr9\t100\tLOW_DEPTH\t0\tLOW_DEPTH\t1\n")
            out = root / "sample.xlsx"
            report.generate_report(
                str(out), "S1", "T", "hg38", "", "", "", "", "", "", "",
                bcr_abl_mutations=muts, bcr_abl_phasing=phase, bcr_abl_all_sites=sites,
            )
            self.assertTrue(out.exists())
            text = self._sheet_text(out, "Overview")
            self.assertIn("T315I", text)
            self.assertIn("SUBCLONAL", text)
            self.assertIn("LOW_SPANNING_DEPTH", text)
            self.assertIn("CIS (testC+testD)", text)
            self.assertIn("panel incomplete (0/1 sites callable)", text)
            self.assertNotIn("POLYCLONAL (TRANS)", text)
            detail = self._sheet_text(out, "BCR-ABL1")
            self.assertIn("L387L", detail)
            self.assertIn("SYNONYMOUS", detail)

    def test_run_summary_does_not_call_empty_panel_sensitive(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            muts = self._write(root, "mut.tsv", "mutation\tstatus\tdepth\tvaf_percent\n")
            sites = self._write(root, "sites.tsv", "chrom\tpos\tstatus\tdepth\tfilter_reason\tpanel_sites_expected\nchr9\t100\tLOW_DEPTH\t0\tLOW_DEPTH\t1\n")
            phase = self._write(root, "phase.tsv", "mutation_1\tmutation_2\tphase_state\n")
            out = root / "run.xlsx"
            summary.build_summary_report(str(out), "RUN", "hg38", [{
                "sample": "S1", "type": "T", "assay": "bcr_abl",
                "bcr_abl_mutations": muts, "bcr_abl_phasing": phase,
                "bcr_abl_all_sites": sites,
            }])
            self.assertTrue(out.exists())
            text = self._sheet_text(out, "BCR-ABL1 Overview")
            self.assertIn("Panel incomplete (0/1 sites callable)", text)
            self.assertIn("Resistance panel incomplete", text)
            self.assertNotIn("Normal / Sensitive", text)
            self.assertNotIn("POLYCLONAL (TRANS)", text)


if __name__ == "__main__":
    unittest.main()
