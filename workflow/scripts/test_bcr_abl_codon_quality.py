"""Regression tests for final haplotype quality and incomplete codon evidence."""
import tempfile
import unittest
from pathlib import Path

import pysam
from bcr_abl_report_utils import panel_assessment
from test_bcr_abl_caller_regressions import caller, phaser, run_caller, write_bam

BED = ("chr9\t99\t100\tNM_005157.6:c.(943A>G)\tA\tG\tp.Thr315Ala\n"
       "chr9\t100\t101\tNM_005157.6:c.(944C>T)\tC\tT\tp.Thr315Ile\n")


def full(codon, flag=0, quality=30, dirty=False):
    sequence = list("N" * 150)
    sequence[99:102] = codon
    if dirty:
        sequence[98] = "A"
    return {"seq": "".join(sequence), "flag": flag, "quals": [quality] * 150}


class CodonQualityTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.tmp = Path(self.directory.name)
        reference = list("N" * 200)
        reference[99:102] = "ACT"
        self.fasta = self.tmp / "ref.fa"
        self.fasta.write_text(">chr9\n" + "".join(reference) + "\n")
        pysam.faidx(str(self.fasta))

    def call(self, reads, extra=()):
        bam = self.tmp / "reads.bam"
        write_bam(bam, reads)
        return run_caller(self.tmp, bam, BED, ("--reference-fasta", str(self.fasta), *extra))

    def test_codon_strand_bias_is_measured_and_filtered(self):
        reads = ([full("ACT")] * 30 + [full("ACT", flag=16)] * 30
                 + [full("GCT", flag=16)] * 20 + [full("ATT", flag=16)] * 20
                 + [full("GTT")] * 20)
        audit, variants, _ = self.call(reads, ("--include-artifacts",))
        self.assertTrue(all(row["status"] == "DETECTED" for row in audit))
        self.assertEqual(len(variants), 3)
        for row in variants:
            self.assertEqual(row["status"], "ARTIFACT")
            self.assertIn("STRAND_BIAS", row["filter_reason"])
            self.assertEqual(float(row["sor"]), 3.04)
            self.assertEqual(float(row["mean_baseq_ref"]), 30)
        self.assertEqual(panel_assessment(audit)["state"], "INCOMPLETE")

    def test_codon_base_quality_filter_cannot_be_bypassed(self):
        reads = [full("ACT")] * 60 + [full("GCT")] * 20 + [full("ATT")] * 20 + [full("GTT", quality=12)] * 10
        audit, variants, _ = self.call(reads, ("--include-artifacts",))
        self.assertTrue(all(row["status"] == "DETECTED" for row in audit))
        valine = next(row for row in variants if row["mutation"] == "T315V")
        self.assertEqual(valine["status"], "ARTIFACT")
        self.assertIn("LOW_BASEQ", valine["filter_reason"])
        self.assertEqual(valine["mean_baseq_alt"], "12.0")

    def test_codon_reads_with_a_dirty_context_are_evidence_for_no_haplotype(self):
        """As for single sites: the dirty GTT reads are not counted, so no T315V
        haplotype is formed from them -- not even as an artefact."""
        reads = [full("ACT")] * 60 + [full("GCT")] * 20 + [full("ATT")] * 20 + [full("GTT", dirty=True)] * 10
        audit, variants, _ = self.call(reads, ("--include-artifacts",))
        self.assertTrue(all(row["status"] == "DETECTED" for row in audit))
        self.assertFalse([row for row in variants if row["mutation"] == "T315V"])

    def test_insufficient_complete_codon_evidence_is_preserved_as_unresolved(self):
        left = [{"seq": "N" * 99 + ("G" if i < 20 else "A")} for i in range(100)]
        right = [{"seq": ("T" if i < 20 else "C") + "T" + "N" * 48, "start": 100} for i in range(100)]
        # Cover no codons, too few codons, or enough WT codons but no ALT codons.
        for full_reference_count in (0, 5, 50):
            with self.subTest(full_reference_count=full_reference_count):
                audit, variants, vcf = self.call(left + [full("ACT")] * full_reference_count + right)
                self.assertEqual(len(variants), 2)
                self.assertTrue(all(row["record_type"] == "CODON_CONTEXT_UNRESOLVED" for row in variants))
                self.assertTrue(all(row["codon_context_status"] == "CODON_CONTEXT_UNRESOLVED" for row in audit))
                self.assertEqual(panel_assessment(audit)["state"], "INCOMPLETE")
                self.assertEqual(phaser.load_detected_mutations(self.tmp / "mutations.tsv"), [])
                records = [line.split("\t") for line in vcf.read_text().splitlines() if not line.startswith("#")]
                self.assertEqual(len(records), 2)
                self.assertTrue(all(row[6] == "CODON_CONTEXT_UNRESOLVED" and row[9].startswith("./.:") for row in records))

    def test_complete_supplementary_codon_survives_partial_primary(self):
        primary = {"name": "molecule", "seq": "GT", "cigar": [(0, 2), (5, 148)], "start": 99}
        supplementary = {"name": "molecule", "seq": "GTT", "cigar": [(5, 99), (0, 3), (5, 48)], "start": 99, "flag": 2048}
        bam = self.tmp / "segments.bam"
        write_bam(bam, [primary, supplementary])
        with pysam.AlignmentFile(str(bam), "rb") as handle:
            counts = caller.collect_codon_haplotypes(handle, "chr9", [99, 100, 101], 12, 1)
        self.assertEqual(counts[("GTT", False)]["count"], 1)

    def test_failed_panel_competitor_remains_in_noise(self):
        rows = []
        for base, count, quality in (("T", 79, 30), ("A", 15, 30), ("C", 6, 12)):
            for _ in range(count):
                sequence = "N" * 49 + base + "N" * 100
                quals = [30] * 150
                quals[49] = quality
                rows.append({"seq": sequence, "quals": quals})
        bam = self.tmp / "competitors.bam"
        write_bam(bam, rows)
        audit, variants, _ = run_caller(self.tmp, bam,
            "chr9\t49\t50\tNM:c.(949T>A)\tT\tA\tp.Phe317Ile\n"
            "chr9\t49\t50\tNM:c.(949T>C)\tT\tC\tp.Phe317Leu\n")
        self.assertEqual(audit[0]["status"], "ARTIFACT")
        self.assertEqual(float(audit[0]["snr"]), 2.5)
        self.assertNotIn("C", audit[0]["credible_panel_alt_competitors"])
        self.assertEqual(variants, [])


if __name__ == "__main__":
    unittest.main()
