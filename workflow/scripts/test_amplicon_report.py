import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from amplicon_report import (
    audit_hotspots,
    derive_flt3_itd_from_variants,
    parse_amplicon_pileup_counts,
    parse_amplicon_stats_json,
    parse_amplicon_vcf,
    parse_flt3_itd_tsv,
    write_amplicon_stats_tsv,
    write_amplicon_summary_tsv,
    write_flt3_itd_tsv,
)

PILEUP_HEADER = (
    "chrom\tpos\tref\talt\traw_depth\tfiltered_depth\tA\tC\tG\tT\tA_fwd\tA_rev\tC_fwd\tC_rev"
    "\tG_fwd\tG_rev\tT_fwd\tT_rev\tINS\tINS_fwd\tINS_rev\tDEL\tDEL_fwd\tDEL_rev\tDEL_span"
    "\tS\tinsertions\tdeletions\tis_capped\tcap_value\n"
)


def write_tmp(content, suffix=".tsv"):
    with tempfile.NamedTemporaryFile("w", suffix=suffix, delete=False) as tf:
        tf.write(content)
        return tf.name


class AmpliconReportTest(unittest.TestCase):
    def test_parse_amplicon_stats_json(self):
        dummy_json = {
            "amplicons": {
                "TP53_3kb_A1": {
                    "count": 1500,
                    "median_length": 2932,
                    "avg_qs": 18.2,
                    "chrom": "chr17",
                    "start": 7667467,
                    "end": 7670399,
                },
                "FLT3_ITD_3kb_1": {
                    "count": 80,
                    "median_length": 3167,
                    "avg_qs": 19.1,
                    "chrom": "chr13",
                    "start": 28032002,
                    "end": 28035169,
                }
            }
        }
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as tf:
            json.dump(dummy_json, tf)
            json_path = tf.name

        rows = parse_amplicon_stats_json(json_path)
        self.assertEqual(len(rows), 2)
        r_tp53 = next(r for r in rows if r["amplicon"] == "TP53_3kb_A1")
        self.assertEqual(r_tp53["gene"], "TP53")
        self.assertEqual(r_tp53["pool"], "1")
        self.assertEqual(r_tp53["read_count"], 1500)
        self.assertEqual(r_tp53["status"], "PASS")

        r_flt3 = next(r for r in rows if r["amplicon"] == "FLT3_ITD_3kb_1")
        self.assertEqual(r_flt3["gene"], "FLT3")
        self.assertEqual(r_flt3["pool"], "3")
        self.assertEqual(r_flt3["read_count"], 80)
        self.assertEqual(r_flt3["status"], "FAIL")

        Path(json_path).unlink()

    def test_parse_amplicon_vcf(self):
        vcf_content = """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
chr17\t7674220\t.\tC\tT\t60\tPASS\tDP=1200;ANN=T|missense_variant|MODERATE|TP53||transcript|exon7||7/11|c.743G>A|p.Arg248Gln\tGT:DP:AD:AF\t0/1:1200:240:0.20
"""
        with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as tf:
            tf.write(vcf_content)
            vcf_path = tf.name

        variants = parse_amplicon_vcf(vcf_path)
        self.assertEqual(len(variants), 1)
        v = variants[0]
        self.assertEqual(v["gene"], "TP53")
        self.assertEqual(v["pos"], 7674220)
        self.assertEqual(v["vaf"], 20.0)
        self.assertEqual(v["hgvs_p"], "p.Arg248Gln")

        Path(vcf_path).unlink()

    def test_parse_flt3_itd_tsv(self):
        """Round-trip of the TSV this module now writes itself."""
        itd_content = (
            "gene\tchrom\tdup_start\tdup_end\tdup_len\taf\taf_pct\tar\tin_frame\t"
            "reading_frame\tar_category\tstatus\n"
            "FLT3\tchr13\t28034000\t28034078\t78\t0.35\t35.0\t0.54\ttrue\tIn-frame\tAR >= 0.5\tDETECTED\n"
        )
        tsv_path = write_tmp(itd_content)

        itds = parse_flt3_itd_tsv(tsv_path)
        self.assertEqual(len(itds), 1)
        self.assertEqual(itds[0]["itd_length"], "78")
        self.assertEqual(itds[0]["af_pct"], 35.0)
        self.assertEqual(itds[0]["status"], "DETECTED")
        self.assertEqual(itds[0]["reading_frame"], "In-frame")
        self.assertEqual(itds[0]["ar"], 0.54)
        self.assertEqual(itds[0]["ar_category"], "AR >= 0.5")

        Path(tsv_path).unlink()

    def test_parse_flt3_itd_tsv_unique_reads_and_frameshift(self):
        itd_content = (
            "sample\tgene\tdup_len\talt_unique_reads\tref_unique_reads\taf\tin_frame\tstatus\n"
            "sample1\tFLT3\t79\t200\t800\t0.20\tfalse\tDETECTED\n"
        )
        tsv_path = write_tmp(itd_content)

        itds = parse_flt3_itd_tsv(tsv_path)
        self.assertEqual(len(itds), 1)
        self.assertEqual(itds[0]["itd_length"], "79")
        self.assertEqual(itds[0]["ar"], 0.25)  # 200 / 800
        self.assertEqual(itds[0]["reading_frame"], "Frameshift")
        self.assertEqual(itds[0]["ar_category"], "AR < 0.5")

        Path(tsv_path).unlink()

    def test_parse_flt3_itd_tsv_tolerates_non_numeric_af(self):
        """A '.' in af must not raise or drop the row (it used to crash)."""
        itd_content = (
            "gene\tchrom\tdup_start\tdup_end\tdup_len\taf\tstatus\n"
            "FLT3\tchr13\t28034000\t28034078\t78\t.\tDETECTED\n"
        )
        tsv_path = write_tmp(itd_content)

        itds = parse_flt3_itd_tsv(tsv_path)
        self.assertEqual(len(itds), 1)
        self.assertEqual(itds[0]["status"], "DETECTED")
        self.assertEqual(itds[0]["af_pct"], 0.0)

        Path(tsv_path).unlink()

    def test_explicit_status_is_never_overridden_by_af(self):
        """A caller-filtered row must not be promoted to DETECTED by its AF."""
        itd_content = "gene\tdup_len\taf_pct\tstatus\nFLT3\t78\t12.0\tFILTERED\n"
        tsv_path = write_tmp(itd_content)

        itds = parse_flt3_itd_tsv(tsv_path)
        self.assertEqual(itds[0]["status"], "NEGATIVE")

        Path(tsv_path).unlink()

    # -- FLT3-ITD derived from the long-read VCF ---------------------------

    def test_derive_flt3_itd_from_insertion(self):
        variants = [{
            "chrom": "chr13", "pos": 28034100, "ref": "A", "alt": "A" + "G" * 78,
            "filter": "PASS", "depth": 1000, "alt_depth": 300, "vaf": 30.0,
            "length_change": 78,
        }]
        rows = derive_flt3_itd_from_variants(variants)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["status"], "DETECTED")
        self.assertEqual(rows[0]["dup_len"], "78")
        self.assertEqual(rows[0]["af_pct"], 30.0)
        self.assertEqual(rows[0]["reading_frame"], "In-frame")
        self.assertEqual(rows[0]["ar"], 0.43)  # 300 / 700
        self.assertEqual(rows[0]["ar_category"], "AR < 0.5")

    def test_derive_flt3_itd_ignores_snvs_and_out_of_region(self):
        variants = [
            {"chrom": "chr13", "pos": 28034100, "ref": "A", "alt": "G", "filter": "PASS",
             "depth": 1000, "alt_depth": 300, "vaf": 30.0, "length_change": 0},
            {"chrom": "chr13", "pos": 28018500, "ref": "A", "alt": "A" + "G" * 30, "filter": "PASS",
             "depth": 1000, "alt_depth": 300, "vaf": 30.0, "length_change": 30},
            {"chrom": "chr13", "pos": 28034120, "ref": "A", "alt": "A" + "G" * 30, "filter": "LowQual",
             "depth": 1000, "alt_depth": 300, "vaf": 30.0, "length_change": 30},
        ]
        rows = derive_flt3_itd_from_variants(variants)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["status"], "NEGATIVE")

    def test_derive_flt3_itd_without_vcf_is_not_assessed(self):
        """A missing call must never be delivered as a negative result."""
        rows = derive_flt3_itd_from_variants([], vcf_available=False)
        self.assertEqual(rows[0]["status"], "NOT_ASSESSED")

    def test_write_and_reread_derived_itd(self):
        variants = [{
            "chrom": "chr13", "pos": 28034100, "ref": "A", "alt": "A" + "G" * 30,
            "filter": "PASS", "depth": 1000, "alt_depth": 600, "vaf": 60.0,
            "length_change": 30,
        }]
        out = write_tmp("")
        write_flt3_itd_tsv(derive_flt3_itd_from_variants(variants), out)
        rows = parse_flt3_itd_tsv(out)
        self.assertEqual(rows[0]["status"], "DETECTED")
        self.assertEqual(rows[0]["itd_length"], "30")
        self.assertEqual(rows[0]["ar_category"], "AR >= 0.5")  # 600 / 400 = 1.5
        Path(out).unlink()

    # -- VCF parsing -------------------------------------------------------

    def test_vaf_falls_back_to_ad_over_dp(self):
        """No AF/VAF in FORMAT must not report a real variant at 0 %."""
        vcf = (
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr17\t7674220\t.\tC\tT\t50\tPASS\tDP=800\tGT:DP:AD\t0/1:800:700,100\n"
        )
        path = write_tmp(vcf, ".vcf")
        variants = parse_amplicon_vcf(path)
        self.assertEqual(variants[0]["vaf"], 12.5)
        Path(path).unlink()

    def test_multiallelic_record_becomes_one_row_per_allele(self):
        vcf = (
            "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
            "chr17\t7674220\t.\tC\tT,A\t50\tPASS\t.\tGT:DP:AD:AF\t0/1:1000:700,200,100:0.2,0.1\n"
        )
        path = write_tmp(vcf, ".vcf")
        variants = parse_amplicon_vcf(path)
        self.assertEqual([v["alt"] for v in variants], ["T", "A"])
        self.assertEqual([v["alt_depth"] for v in variants], [200, 100])
        self.assertEqual([v["vaf"] for v in variants], [20.0, 10.0])
        Path(path).unlink()

    # -- Hotspot audit -----------------------------------------------------

    def test_parse_amplicon_pileup_counts_and_audit(self):
        pileup_content = PILEUP_HEADER + (
            "chr17\t7675088\tC\tT\t1000\t1000\t1\t995\t2\t2\t1\t0\t500\t495\t1\t1\t1\t1\t0\t0\t0\t0\t0\t0\t0\t0\t.\t.\t0\t0\n"
            "chr17\t7674220\tC\tT\t1000\t1000\t5\t975\t10\t10\t2\t3\t490\t485\t5\t5\t5\t5\t0\t0\t0\t0\t0\t0\t0\t0\t.\t.\t0\t0\n"
            "chr13\t28018505\tC\tT\t200\t200\t0\t200\t0\t0\t0\t0\t100\t100\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t.\t.\t0\t0\n"
        )
        pileup_path = write_tmp(pileup_content)

        pileup_map = parse_amplicon_pileup_counts(pileup_path)
        self.assertIn(("chr17", 7675088), pileup_map)

        audit_res = audit_hotspots(pileup_map, [])
        r175 = next(r for r in audit_res if r["codon"] == "p.R175")
        self.assertEqual(r175["depth"], 1000)
        self.assertEqual(r175["non_ref_reads"], 5)
        self.assertEqual(r175["noise_pct"], 0.5)
        self.assertIn("CLEAN", r175["status"])

        r248 = next(r for r in audit_res if r["codon"] == "p.R248")
        self.assertEqual(r248["noise_pct"], 2.5)
        self.assertIn("ELEVATED NOISE", r248["status"])

        d835 = next(r for r in audit_res if r["codon"] == "p.D835")
        self.assertIn("LOW DEPTH", d835["status"])

        w288 = next(r for r in audit_res if r["codon"] == "p.W288")
        self.assertEqual(w288["status"], "NOT ASSESSED")

        Path(pileup_path).unlink()

    def test_hotspot_coordinates_are_grch38(self):
        """Guards the published GRCh38 positions of the reported hotspots.

        R172K, D835Y and the NPM1c+ duplication agree with the coordinates used
        by neville_mx_amplicon, which calls the same panel.
        """
        from amplicon_report import CLINICAL_HOTSPOTS
        by_codon = {h["codon"]: h for h in CLINICAL_HOTSPOTS}
        self.assertEqual((by_codon["p.R175"]["chrom"], by_codon["p.R175"]["pos"]), ("chr17", 7675088))
        self.assertEqual((by_codon["p.R248"]["chrom"], by_codon["p.R248"]["pos"]), ("chr17", 7674220))
        self.assertEqual((by_codon["p.R273"]["chrom"], by_codon["p.R273"]["pos"]), ("chr17", 7673802))
        self.assertEqual((by_codon["p.R132"]["chrom"], by_codon["p.R132"]["pos"]), ("chr2", 208248388))
        self.assertEqual((by_codon["p.R140"]["chrom"], by_codon["p.R140"]["pos"]), ("chr15", 90088702))
        self.assertEqual((by_codon["p.R172"]["chrom"], by_codon["p.R172"]["pos"]), ("chr15", 90088606))
        self.assertEqual((by_codon["p.D835"]["chrom"], by_codon["p.D835"]["pos"]), ("chr13", 28018505))
        # R172 is 32 codons from R140 on the minus strand, and in exon 4.
        self.assertEqual(by_codon["p.R140"]["pos"] - by_codon["p.R172"]["pos"], (172 - 140) * 3)
        self.assertEqual(by_codon["p.R172"]["exon"], "Exon 4")

    def test_hotspot_matches_any_base_of_the_codon(self):
        """A call on the third base of R175 belongs to R175."""
        variants = [{
            "chrom": "chr17", "pos": 7675087, "ref": "G", "alt": "A", "filter": "PASS",
            "gene": "TP53", "hgvs_p": "p.R175H", "vaf": 22.0, "depth": 900,
            "alt_depth": 198, "length_change": 0,
        }]
        r175 = next(r for r in audit_hotspots({}, variants) if r["codon"] == "p.R175")
        self.assertIn("MUTATION DETECTED", r175["status"])
        self.assertIn("p.R175H", r175["status"])

    def test_hotspot_does_not_claim_a_variant_from_another_gene(self):
        variants = [{
            "chrom": "chr17", "pos": 7675088, "ref": "C", "alt": "T", "filter": "PASS",
            "gene": "SOMETHINGELSE", "hgvs_p": "p.X1Y", "vaf": 22.0, "depth": 900,
            "alt_depth": 198, "length_change": 0,
        }]
        r175 = next(r for r in audit_hotspots({}, variants) if r["codon"] == "p.R175")
        self.assertEqual(r175["status"], "NOT ASSESSED")

    def test_npm1_insertion_is_linked_to_w288(self):
        """The NPM1c+ duplication is left-aligned to chr5:171410539 (as in neville)."""
        variants = [{
            "chrom": "chr5", "pos": 171410539, "ref": "C", "alt": "CTCTG", "filter": "PASS",
            "gene": "NPM1", "hgvs_p": "p.W288Cfs*12", "vaf": 38.0, "depth": 1200,
            "alt_depth": 456, "length_change": 4,
        }]
        w288 = next(r for r in audit_hotspots({}, variants) if r["codon"] == "p.W288")
        self.assertIn("MUTATION DETECTED", w288["status"])

    def test_write_amplicon_summary_tsv_with_hotspots(self):
        out_path = write_tmp("")

        itd_rows = [{
            "gene": "FLT3", "status": "DETECTED", "itd_length": "54", "af_pct": 25.0,
            "ar": 0.33, "reading_frame": "In-frame", "ar_category": "AR < 0.5",
            "region": "chr13:28033864-28034426",
        }]
        hotspot_rows = [{
            "gene": "TP53", "codon": "p.R175", "exon": "Exon 5", "chrom": "chr17",
            "pos": 7675088, "ref": "C", "depth": 1500, "non_ref_reads": 6,
            "noise_pct": 0.40, "status": "CLEAN (Noise 0.40%)",
        }]
        variants = [
            {"chrom": "chr17", "pos": 7674220, "gene": "TP53", "ref": "C", "alt": "T",
             "vaf": 12.5, "depth": 800, "alt_depth": 100, "filter": "PASS"},
            {"chrom": "chr17", "pos": 7674999, "gene": "TP53", "ref": "A", "alt": "AT",
             "vaf": 3.0, "depth": 800, "alt_depth": 24, "filter": "LowQual"},
        ]
        write_amplicon_summary_tsv([], variants, itd_rows, out_path, sample="S1", hotspot_rows=hotspot_rows)

        content = Path(out_path).read_text()
        self.assertIn("# Section: FLT3_ITD", content)
        self.assertIn("ar\treading_frame\tar_category", content)
        self.assertIn("0.33\tIn-frame\tAR < 0.5", content)
        self.assertIn("# Section: CLINICAL_HOTSPOT_BACKGROUND_NOISE", content)
        self.assertIn("p.R175\tExon 5\tchr17\t7675088\tC\t1500\t6\t0.4\tCLEAN (Noise 0.40%)", content)

        # Filtered calls are delivered, but not inside the clinical variant list.
        clinical = content.split("# Section: FILTERED_VARIANTS")[0]
        self.assertIn("7674220", clinical)
        self.assertNotIn("7674999", clinical)
        self.assertIn("7674999", content)

        Path(out_path).unlink()


if __name__ == "__main__":
    unittest.main()
