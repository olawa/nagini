"""Integration checks for the files exchanged by the caller and phaser."""

import json
import tempfile
import unittest
from pathlib import Path

import bcr_abl_phaser as phaser
from test_bcr_abl_caller_regressions import run_caller, write_bam


class CallerPhaserIntegrationTests(unittest.TestCase):
    def test_detected_single_site_variants_reach_phasing_as_a_cis_pair(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            reads = []
            for index in range(100):
                sequence = list("N" * 150)
                sequence[49] = "A" if index < 40 else "G"
                sequence[99] = "T" if index < 40 else "C"
                reads.append({"name": f"molecule{index}", "seq": "".join(sequence)})
            bam = tmp / "reads.bam"
            write_bam(bam, reads)
            audit, variants, _ = run_caller(
                tmp, bam,
                "chr9\t49\t50\tNM_005157.6:c.(763G>A)\tG\tA\tp.Glu255Lys\n"
                "chr9\t99\t100\tNM_005157.6:c.(944C>T)\tC\tT\tp.Thr315Ile\n",
            )
            self.assertEqual(len(variants), 2)
            self.assertTrue(all(row["record_type"] == "SITE_AUDIT" for row in audit))
            self.assertEqual({row["mutation"] for row in variants}, {"E255K", "T315I"})

            output_json = tmp / "phasing.json"
            phaser.main([
                "--bam", str(bam), "--mutations-tsv", str(tmp / "mutations.tsv"),
                "--sample", "S", "--output-tsv", str(tmp / "phasing.tsv"),
                "--output-json", str(output_json),
            ])
            result = json.loads(output_json.read_text())
            self.assertEqual(result["detected_mutations_count"], 2)
            self.assertEqual(result["phased_pairs_count"], 1)
            pair = result["phased_pairs"][0]
            self.assertEqual(pair["phase_state"], "CIS (Compound)")
            self.assertEqual(pair["cis_reads_11"], 40)
            self.assertEqual(pair["wt_reads_00"], 60)
            self.assertEqual(pair["usable_spanning_reads"], 100)


if __name__ == "__main__":
    unittest.main()
