import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

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

    def test_multi_caller_consolidation_in_individual_and_summary_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            muts = self._write(
                root, "mut.tsv",
                "mutation\taa_change\tconsequence\tresistance_class\tstatus\tchrom\tpos\tref\talt\tdepth\tad_alt\tvaf_percent\n"
                "T315I\tp.Thr315Ile\tMISSENSE\tCURATED_RESISTANCE\tDETECTED\tchr9\t133748283\tC\tT\t1000\t250\t25.0\n"
            )
            vardict = self._write(
                root, "vardict.vcf",
                "##fileformat=VCFv4.2\n"
                "##source=VarDictRust\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n"
                "##FORMAT=<ID=VD,Number=1,Type=Integer,Description=\"Variant Depth\">\n"
                "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr9\t133748283\t.\tC\tT\t100\tPASS\t.\tGT:DP:VD:AF\t0/1:1000:240:0.240\n"
                "chr9\t133750500\t.\tG\tA\t80\tPASS\t.\tGT:DP:VD:AF\t0/1:800:160:0.200\n"
            )
            bcftools = self._write(
                root, "bcftools.vcf",
                "##fileformat=VCFv4.2\n"
                "##source=bcftools\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr9\t133748283\t.\tC\tT\t100\tPASS\t.\tGT:DP:AD\t0/1:1000:760,240\n"
            )
            rindels = self._write(
                root, "rindels.vcf",
                "##fileformat=VCFv4.2\n"
                "##source=rindels\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
                "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr9\t133748283\t.\tC\tT\t100\tPASS\t.\tGT:DP:AD:AF\t0/1:1000:750,250:0.250\n"
            )

            # Test individual sample report
            out_sample = root / "sample.xlsx"
            report.generate_report(
                str(out_sample), "S1", "T", "hg19", "", "", "", "", "", "", "",
                bcr_abl_mutations=muts,
                bcr_abl_vardictrust_vcf=vardict,
                bcr_abl_rindels_vcf=rindels,
                bcr_abl_bcftools_vcf=bcftools,
            )
            self.assertTrue(out_sample.exists())

            overview_txt = self._sheet_text(out_sample, "Overview")
            self.assertIn("Small Variant Caller", overview_txt)
            self.assertIn("Branford", overview_txt)
            self.assertIn("VarDictRust", overview_txt)

            snv_txt = self._sheet_text(out_sample, "SNVs and indels")
            self.assertIn("Caller", snv_txt)
            self.assertIn("Branford", snv_txt)
            self.assertIn("VarDictRust", snv_txt)
            self.assertIn("bcftools", snv_txt)
            self.assertIn("Rindels", snv_txt)
            self.assertIn("133748283", snv_txt)
            self.assertIn("133750500", snv_txt)

            bcr_abl_txt = self._sheet_text(out_sample, "BCR-ABL1")
            self.assertIn("Supporting Callers", bcr_abl_txt)
            self.assertIn("VarDictRust", bcr_abl_txt)
            self.assertIn("bcftools", bcr_abl_txt)
            self.assertIn("Rindels", bcr_abl_txt)
            self.assertIn("ADDITIONAL ABL1 EXON VARIANTS (OTHER CALLERS)", bcr_abl_txt)
            self.assertIn("133750500", bcr_abl_txt)

            # Test summary report
            out_summary = root / "summary.xlsx"
            summary.build_summary_report(
                str(out_summary), "RUN_01", "hg19", [{
                    "sample": "S1", "type": "T", "assay": "bcr_abl",
                    "bcr_abl_mutations": muts,
                    "bcr_abl_vardictrust_vcf": vardict,
                    "bcr_abl_rindels_vcf": rindels,
                    "bcr_abl_bcftools_vcf": bcftools,
                }]
            )
            self.assertTrue(out_summary.exists())
            pass_var_txt = self._sheet_text(out_summary, "PASS Variants")
            self.assertIn("Caller", pass_var_txt)
            self.assertIn("Branford", pass_var_txt)
            self.assertIn("VarDictRust", pass_var_txt)
            self.assertIn("bcftools", pass_var_txt)
            self.assertIn("Rindels", pass_var_txt)

    def test_low_depth_caller_filtering_and_vaf_formatting(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            muts = self._write(
                root, "mut.tsv",
                "mutation\taa_change\tconsequence\tresistance_class\tstatus\tchrom\tpos\tref\talt\tdepth\tad_alt\tvaf_percent\n"
                "M351T\tp.Met351Thr\tMISSENSE\tCURATED_RESISTANCE\tDETECTED\tchr9\t133748391\tT\tC\t10000\t9500\t95.00\n"
            )
            sites = self._write(
                root, "sites.tsv",
                "chrom\tpos\tstatus\tdepth\tfilter_reason\tpanel_sites_expected\n"
                "chr9\t133748391\tPASS\t10000\tNONE\t1\n"
            )
            rindels = self._write(
                root, "rindels.vcf",
                "##fileformat=VCFv4.2\n"
                "##source=rindels\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
                "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr9\t133748391\t.\tT\tC\t100\tPASS\t.\tGT:DP:AD:AF\t0/1:10000:500,9500:0.950\n"
                "chr9\t133755936\t.\tC\tCT\t20\tPASS\t.\tGT:DP:AD:AF\t0/1:20:4,16:0.800\n"
            )

            # Test individual report filters out the 20x call (cutoff is 1000x at 10000x mean depth)
            out_sample = root / "sample.xlsx"
            report.generate_report(
                str(out_sample), "S1", "T", "hg19", "", "", "", "", "", "", "",
                bcr_abl_mutations=muts,
                bcr_abl_all_sites=sites,
                bcr_abl_rindels_vcf=rindels,
                bcr_abl_min_relative_depth=0.10,
                bcr_abl_min_depth=50.0,
            )
            self.assertTrue(out_sample.exists())
            snv_txt = self._sheet_text(out_sample, "SNVs and indels")
            self.assertIn("133748391", snv_txt)
            self.assertNotIn("133755936", snv_txt)
            self.assertIn("95.00%", snv_txt)

            # Test summary report also filters out 20x call and formats VAF as %
            out_summary = root / "summary.xlsx"
            summary.build_summary_report(
                str(out_summary), "RUN_01", "hg19", [{
                    "sample": "S1", "type": "T", "assay": "bcr_abl",
                    "bcr_abl_mutations": muts,
                    "bcr_abl_all_sites": sites,
                    "bcr_abl_rindels_vcf": rindels,
                }],
                bcr_abl_min_relative_depth=0.10,
                bcr_abl_min_depth=50.0,
            )
            self.assertTrue(out_summary.exists())
            pass_var_txt = self._sheet_text(out_summary, "PASS Variants")
            self.assertIn("133748391", pass_var_txt)
            self.assertNotIn("133755936", pass_var_txt)
            self.assertIn("95.00%", pass_var_txt)


    # Rows exactly as a real run reported them. The BCR-ABL1 sheet showed these
    # correctly; the SNVs and indels sheet of the same workbook showed T315A as
    # 65.20% and E499E as 41.40%, because format_vaf guessed that any value up
    # to 1 was a fraction -- and the Branford caller writes percentages.
    REAL_ROWS = (
        "mutation\taa_change\tconsequence\tresistance_class\tstatus\tchrom\tpos\tref\talt\tdepth\tad_alt\tvaf_percent\n"
        "Y253H\tp.Tyr253His\tMISSENSE\tCURATED_RESISTANCE\tDETECTED\tchr9\t133738357\tT\tC\t1966\t307\t15.615\n"
        "T315I\tp.Thr315Ile\tMISSENSE\tCURATED_RESISTANCE\tDETECTED\tchr9\t133748282\tACT\tATT\t1993\t309\t15.504\n"
        "T315A\tp.Thr315Ala\tMISSENSE\tCURATED_RESISTANCE\tSUBCLONAL\tchr9\t133748282\tACT\tGCT\t1993\t13\t0.652\n"
        "E499E\tp.Glu499Glu\tSYNONYMOUS\tSYNONYMOUS\tSUBCLONAL\tchr9\t133755528\tA\tG\t1933\t8\t0.414\n"
    )

    @staticmethod
    def _rows_by(path, sheet, key_column):
        book = openpyxl.load_workbook(path, read_only=True, data_only=True)
        try:
            rows = [list(r) for r in book[sheet].iter_rows(values_only=True)]
        finally:
            book.close()
        for i, row in enumerate(rows):
            if key_column in row:
                header = row
                return [dict(zip(header, r)) for r in rows[i + 1:] if any(v is not None for v in r)]
        raise AssertionError(f"{key_column} not found in {sheet}")

    def test_subclonal_branford_vaf_is_not_scaled_by_100(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            muts = self._write(root, "mut.tsv", self.REAL_ROWS)
            out = root / "sample.xlsx"
            report.generate_report(str(out), "S1", "T", "hg19", "", "", "", "", "", "", "", bcr_abl_mutations=muts)
            snvs = {r["HGVSp"]: r for r in self._rows_by(out, "SNVs and indels", "HGVSp") if r.get("HGVSp")}
            self.assertEqual(snvs["p.Tyr253His"]["VAF"], "15.62%")
            self.assertEqual(snvs["p.Thr315Ala"]["VAF"], "0.65%")
            self.assertEqual(snvs["p.Glu499Glu"]["VAF"], "0.41%")

    def test_branford_rows_carry_no_invented_annotation(self):
        """Nothing in a ClinVar or VEP column that ClinVar or VEP did not say."""
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            muts = self._write(root, "mut.tsv", self.REAL_ROWS)
            out = root / "sample.xlsx"
            report.generate_report(str(out), "S1", "T", "hg19", "", "", "", "", "", "", "", bcr_abl_mutations=muts)
            for row in self._rows_by(out, "SNVs and indels", "HGVSp"):
                if not row.get("HGVSp"):
                    continue
                for column in ("Impact", "ClinVar_Significance", "ClinVar_Disease"):
                    self.assertIn(row.get(column), (None, ""), f"{row['HGVSp']}: {column}={row.get(column)!r}")
            text = self._sheet_text(out, "SNVs and indels")
            self.assertNotIn("Resistance mutation", text)
            self.assertNotIn("Chronic Myeloid Leukemia", text)

            summ = root / "summary.xlsx"
            summary.build_summary_report(str(summ), "RUN", "hg19", [{"sample": "S1", "type": "T", "assay": "bcr_abl", "bcr_abl_mutations": muts}])
            self.assertNotIn("Resistance mutation", self._sheet_text(summ, "PASS Variants"))

    def test_detected_call_at_exactly_one_percent(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            muts = self._write(
                root, "mut.tsv",
                "mutation\taa_change\tconsequence\tresistance_class\tstatus\tchrom\tpos\tref\talt\tdepth\tad_alt\tvaf_percent\n"
                "T315I\tp.Thr315Ile\tMISSENSE\tCURATED_RESISTANCE\tDETECTED\tchr9\t133748282\tACT\tATT\t2000\t20\t1.0\n",
            )
            summ = root / "summary.xlsx"
            summary.build_summary_report(str(summ), "RUN", "hg19", [{"sample": "S1", "type": "T", "assay": "bcr_abl", "bcr_abl_mutations": muts}])
            text = self._sheet_text(summ, "PASS Variants")
            self.assertIn("1.00%", text)
            self.assertNotIn("100.00%", text)

    def test_supporting_caller_vaf_is_not_formatted_twice(self):
        """The SNV list is formatted in place before the BCR-ABL1 sheet reads it.

        Y253H rather than T315A: the caller writes codon calls as MNVs
        (133748282 ACT>GCT) and the VCF callers the normalised SNV
        (133748283 C>T), which the supporting-caller match does not join.
        """
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            muts = self._write(root, "mut.tsv", self.REAL_ROWS)
            rindels = self._write(
                root, "rindels.vcf",
                "##fileformat=VCFv4.2\n##source=rindels\n"
                "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
                "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">\n"
                "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depths\">\n"
                "##FORMAT=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr9\t133738357\t.\tT\tC\t50\tPASS\t.\tGT:DP:AD:AF\t0/1:1966:1952,14:0.007\n",
            )
            out = root / "sample.xlsx"
            report.generate_report(
                str(out), "S1", "T", "hg19", "", "", "", "", "", "", "",
                bcr_abl_mutations=muts, bcr_abl_rindels_vcf=rindels,
                bcr_abl_min_relative_depth=0.0, bcr_abl_min_depth=0.0,
            )
            text = self._sheet_text(out, "BCR-ABL1")
            self.assertIn("(0.70%)", text)
            self.assertNotIn("(70.00%)", text)


class FormatVafTests(unittest.TestCase):
    """The scale cannot be read off the number, so the caller must say it."""

    def setUp(self):
        from bcr_abl_report_utils import format_vaf
        self.format_vaf = format_vaf

    def test_percentages_stay_percentages(self):
        for value, shown in ((0.2, "0.20%"), (0.652, "0.65%"), (1.0, "1.00%"), (15.615, "15.62%"), (100, "100.00%")):
            self.assertEqual(self.format_vaf(value, scale="percent"), shown)

    def test_fractions_become_percentages(self):
        for value, shown in ((0.002, "0.20%"), ("0.007", "0.70%"), (0.95, "95.00%"), (1.0, "100.00%")):
            self.assertEqual(self.format_vaf(value, scale="fraction"), shown)

    def test_formatting_twice_changes_nothing(self):
        once = self.format_vaf(0.005, scale="fraction")
        self.assertEqual(once, "0.50%")
        self.assertEqual(self.format_vaf(once, scale="fraction"), "0.50%")
        self.assertEqual(self.format_vaf(once, scale="percent"), "0.50%")

    def test_missing_values(self):
        for value in (None, "", ".", "-"):
            self.assertEqual(self.format_vaf(value, scale="fraction"), "-")

    def test_the_scale_is_required(self):
        with self.assertRaises(TypeError):
            self.format_vaf(0.5)
        with self.assertRaises(ValueError):
            self.format_vaf(0.5, scale="ratio")


if __name__ == "__main__":
    unittest.main()
