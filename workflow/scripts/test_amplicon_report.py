import json
import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from amplicon_report import (
    parse_amplicon_stats_json,
    parse_amplicon_vcf,
    parse_flt3_itd_tsv,
    write_amplicon_stats_tsv,
    write_amplicon_summary_tsv,
)


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
        itd_content = "gene\tregion\tlength\taf\tstatus\nFLT3\tchr13:28032002-28035169\t78\t0.35\tDETECTED\n"
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
            tf.write(itd_content)
            tsv_path = tf.name

        itds = parse_flt3_itd_tsv(tsv_path)
        self.assertEqual(len(itds), 1)
        self.assertEqual(itds[0]["itd_length"], "78")
        self.assertEqual(itds[0]["af"], "0.35")
        self.assertEqual(itds[0]["status"], "DETECTED")
        self.assertEqual(itds[0]["reading_frame"], "In-frame")
        self.assertEqual(itds[0]["ar"], 0.54)  # 35 / (100 - 35) = 0.538 -> 0.54
        self.assertEqual(itds[0]["eln_risk"], "HIGH AR (>=0.5)")

        Path(tsv_path).unlink()

    def test_parse_flt3_itd_tsv_unique_reads_and_frameshift(self):
        itd_content = (
            "sample\tgene\tdup_len\talt_unique_reads\tref_unique_reads\taf\tin_frame\tstatus\n"
            "sample1\tFLT3\t79\t200\t800\t0.20\tfalse\tDETECTED\n"
        )
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
            tf.write(itd_content)
            tsv_path = tf.name

        itds = parse_flt3_itd_tsv(tsv_path)
        self.assertEqual(len(itds), 1)
        self.assertEqual(itds[0]["itd_length"], "79")
        self.assertEqual(itds[0]["ar"], 0.25)  # 200 / 800 = 0.25
        self.assertEqual(itds[0]["reading_frame"], "Frameshift")
        self.assertEqual(itds[0]["eln_risk"], "LOW AR (<0.5)")

        Path(tsv_path).unlink()

    def test_parse_amplicon_pileup_counts_and_audit(self):
        from amplicon_report import parse_amplicon_pileup_counts, audit_hotspots

        pileup_content = (
            "chrom\tpos\tref\talt\traw_depth\tfiltered_depth\tA\tC\tG\tT\tA_fwd\tA_rev\tC_fwd\tC_rev\tG_fwd\tG_rev\tT_fwd\tT_rev\tINS\tINS_fwd\tINS_rev\tDEL\tDEL_fwd\tDEL_rev\tDEL_span\tS\tinsertions\tdeletions\tis_capped\tcap_value\n"
            # TP53 p.R175 (chr17:7675107, ref C): 1000x depth, 5 non-ref -> 0.50% noise -> CLEAN
            "chr17\t7675107\tC\tT\t1000\t1000\t1\t995\t2\t2\t1\t0\t500\t495\t1\t1\t1\t1\t0\t0\t0\t0\t0\t0\t0\t0\t.\t.\t0\t0\n"
            # TP53 p.R248 (chr17:7674220, ref C): 1000x depth, 25 non-ref -> 2.50% noise -> ELEVATED NOISE
            "chr17\t7674220\tC\tT\t1000\t1000\t5\t975\t10\t10\t2\t3\t490\t485\t5\t5\t5\t5\t0\t0\t0\t0\t0\t0\t0\t0\t.\t.\t0\t0\n"
            # FLT3 p.D835 (chr13:28018504, ref G): 200x depth -> LOW DEPTH
            "chr13\t28018504\tG\tT\t200\t200\t0\t0\t200\t0\t0\t0\t0\t0\t100\t100\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t.\t.\t0\t0\n"
        )
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
            tf.write(pileup_content)
            pileup_path = tf.name

        pileup_map = parse_amplicon_pileup_counts(pileup_path)
        self.assertIn(("chr17", 7675107), pileup_map)
        self.assertIn(("chr17", 7674220), pileup_map)
        self.assertIn(("chr13", 28018504), pileup_map)

        # Audit with no called variants
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

        # Hotspot without pileup data
        w288 = next(r for r in audit_res if r["codon"] == "p.W288")
        self.assertEqual(w288["status"], "NOT ASSESSED")

        Path(pileup_path).unlink()

    def test_write_amplicon_summary_tsv_with_hotspots(self):
        from amplicon_report import write_amplicon_summary_tsv

        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
            out_path = tf.name

        itd_rows = [{
            "gene": "FLT3",
            "status": "DETECTED",
            "itd_length": "54",
            "af_pct": 25.0,
            "ar": 0.33,
            "reading_frame": "In-frame",
            "eln_risk": "LOW AR (<0.5)",
            "region": "chr13:28032002-28035169",
        }]
        hotspot_rows = [{
            "gene": "TP53",
            "codon": "p.R175",
            "exon": "Exon 5",
            "chrom": "chr17",
            "pos": 7675107,
            "ref": "C",
            "depth": 1500,
            "non_ref_reads": 6,
            "noise_pct": 0.40,
            "status": "CLEAN (Noise 0.40%)",
        }]
        write_amplicon_summary_tsv([], [], itd_rows, out_path, sample="S1", hotspot_rows=hotspot_rows)

        content = Path(out_path).read_text()
        self.assertIn("# Section: FLT3_ITD", content)
        self.assertIn("ar\treading_frame\teln_risk", content)
        self.assertIn("0.33\tIn-frame\tLOW AR (<0.5)", content)
        self.assertIn("# Section: CLINICAL_HOTSPOT_BACKGROUND_NOISE", content)
        self.assertIn("p.R175\tExon 5\tchr17\t7675107\tC\t1500\t6\t0.4\tCLEAN (Noise 0.40%)", content)

        Path(out_path).unlink()


if __name__ == "__main__":
    unittest.main()
