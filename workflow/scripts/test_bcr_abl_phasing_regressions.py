"""Regression tests for conservative BCR-ABL read-backed phasing."""

import importlib.util
import tempfile
import unittest
from unittest import mock
from pathlib import Path

import pysam


HERE = Path(__file__).parent
SPEC = importlib.util.spec_from_file_location("bcr_abl_phaser", HERE / "bcr_abl_phaser.py")
PHASER = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PHASER)


def write_bam(path: Path, records, contigs=("chr9",)):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": contig, "LN": 500} for contig in contigs]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for index, record in enumerate(records):
            alignment = pysam.AlignedSegment()
            alignment.query_name = record.get("name", f"read{index}")
            alignment.query_sequence = record["seq"]
            alignment.query_qualities = record.get("qualities", [30] * len(record["seq"]))
            alignment.reference_id = contigs.index(record.get("contig", contigs[0]))
            alignment.reference_start = record.get("start", 0)
            alignment.mapping_quality = record.get("mapq", 60)
            alignment.flag = record.get("flag", 0)
            alignment.cigartuples = record.get("cigar", [(0, len(record["seq"]))])
            out.write(alignment)
    pysam.index(str(path))


def sequence(base1="A", base2="T"):
    seq = list("A" * 150)
    seq[49] = base1
    seq[99] = base2
    return "".join(seq)


MUT1 = {"chrom": "chr9", "pos": 50, "ref": "G", "alt": "A", "mutation": "E255K"}
MUT2 = {"chrom": "chr9", "pos": 100, "ref": "C", "alt": "T", "mutation": "T315I"}


class TestBcrAblPhasingRegressions(unittest.TestCase):
    def phase(self, tmp, records, mut1=MUT1, mut2=MUT2, min_spanning=10, min_bq=12):
        path = tmp / "reads.bam"
        write_bam(path, records, contigs=("chr9", "chr8"))
        with pysam.AlignmentFile(str(path), "rb") as bam:
            return PHASER.phase_mutation_pair(bam, mut1, mut2, "S", min_spanning, 1, min_bq)

    def test_trans_requires_support_for_both_mutant_only_haplotypes(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            records = ([{"seq": sequence("A", "C")}] * 10 + [{"seq": sequence("G", "C")}] * 90)
            result = self.phase(tmp, records)
            self.assertNotEqual(result["phase_state"], "TRANS (Polyclonal)")
            self.assertEqual(result["trans_reads_01"], 0)

    def test_nested_compound_subclone_is_not_called_trans(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            records = ([{"seq": sequence("A", "T")}] * 10 + [{"seq": sequence("A", "C")}] * 90)
            result = self.phase(tmp, records)
            self.assertEqual(result["phase_state"], "AMBIGUOUS")
            self.assertEqual(result["cis_reads_11"], 10)
            self.assertEqual(result["trans_reads_10"], 90)
            self.assertEqual(result["cis_fraction_given_alt2_percent"], 100.0)

    def test_low_quality_and_refskip_do_not_satisfy_spanning_depth(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            low_q = [30] * 150
            low_q[99] = 0
            records = ([{"seq": sequence("A", "T"), "qualities": low_q}] * 20 +
                       [{"seq": sequence("G", "T")}] * 20 + [{"seq": sequence("G", "C")}] * 60)
            result = self.phase(tmp, records)
            self.assertEqual(result["cis_reads_11"], 0)
            self.assertEqual(result["quality_filtered_alignments"], 20)

            skip_records = ([{"seq": sequence("A", "T")}] * 3 +
                            [{"seq": "A" * 150, "cigar": [(0, 75), (3, 50), (0, 75)]}] * 7)
            result = self.phase(tmp, skip_records)
            self.assertEqual(result["total_spanning_reads"], 3)
            self.assertEqual(result["unusable_base_alignments"], 7)
            self.assertEqual(result["phase_state"], "LOW_SPANNING_DEPTH")

    def test_duplicate_alignments_do_not_inflate_support_and_secondary_is_skipped(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            records = ([{"name": "one", "seq": sequence(), "flag": 0}] +
                       [{"name": "one", "seq": sequence(), "flag": 2048}] * 9 +
                       [{"name": "secondary", "seq": sequence(), "flag": 256}] * 10)
            result = self.phase(tmp, records)
            self.assertEqual(result["total_spanning_reads"], 1)
            self.assertEqual(result["duplicate_alignments"], 9)
            self.assertEqual(result["phase_state"], "LOW_SPANNING_DEPTH")

    def test_conflicting_equal_priority_alignments_are_excluded(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            records = [
                {"name": "conflict", "seq": sequence("A", "T"), "flag": 0},
                {"name": "conflict", "seq": sequence("G", "C"), "flag": 0},
            ]
            result = self.phase(tmp, records)
            self.assertEqual(result["conflicting_molecules"], 1)
            self.assertEqual(result["total_spanning_reads"], 0)

    def test_supplementary_only_alignment_can_support_fusion_pair(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            records = [{"seq": sequence(), "flag": 2048} for _ in range(10)]
            result = self.phase(tmp, records)
            self.assertEqual(result["cis_reads_11"], 10)
            self.assertEqual(result["phase_state"], "CIS (Compound)")

    def test_non_distinct_positions_and_indels_are_explicitly_unphaseable(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            same_pos = {**MUT2, "pos": 50}
            result = self.phase(tmp, [], mut2=same_pos)
            self.assertEqual(result["phase_state"], "UNPHASEABLE_SAME_POSITION")
            diff_contig = {**MUT2, "chrom": "chr8"}
            result = self.phase(tmp, [], mut2=diff_contig)
            self.assertEqual(result["phase_state"], "UNPHASEABLE_DIFFERENT_CONTIG")
            indel = {**MUT1, "alt": "AT"}
            result = self.phase(tmp, [], mut1=indel)
            self.assertEqual(result["phase_state"], "UNPHASEABLE_INDEL")

    def test_same_length_codon_haplotype_is_read_as_one_allele(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            codon = {"chrom": "chr9", "pos": 50, "ref": "ACT", "alt": "GTT", "mutation": "T315V"}
            other = {"chrom": "chr9", "pos": 100, "ref": "C", "alt": "T", "mutation": "T315I"}
            def seq(codon_value, second):
                chars = list("A" * 150)
                chars[49:52] = codon_value
                chars[99] = second
                return "".join(chars)
            records = ([{"seq": seq("GTT", "T")}] * 10 + [{"seq": seq("ACT", "C")}] * 10)
            result = self.phase(tmp, records, codon, other)
            self.assertEqual(result["cis_reads_11"], 10)
            self.assertEqual(result["phase_state"], "CIS (Compound)")

    def test_insertion_inside_codon_is_not_joined_into_a_codon_allele(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            codon = {"chrom": "chr9", "pos": 50, "ref": "ACT", "alt": "GTT", "mutation": "T315V"}
            other = {"chrom": "chr9", "pos": 100, "ref": "C", "alt": "T", "mutation": "T315I"}
            seq = list("A" * 151)
            seq[49], seq[51], seq[52], seq[100] = "G", "T", "T", "T"
            records = [{"seq": "".join(seq), "cigar": [(0, 50), (1, 1), (0, 100)]} for _ in range(10)]
            result = self.phase(tmp, records, codon, other)
            self.assertEqual(result["total_spanning_reads"], 0)
            self.assertEqual(result["unusable_base_alignments"], 10)

    def test_loader_excludes_synonymous_and_unresolved_records(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mutations.tsv"
            path.write_text(
                "chrom\tpos\tref\talt\tstatus\tmutation\trecord_type\tconsequence\tresistance_class\n"
                "chr9\t10\tA\tG\tDETECTED\tL387L\tSITE_AUDIT\tMISSENSE\tUNCLASSIFIED\n"
                "chr9\t20\tA\tG\tDETECTED\tCODON_CONTEXT_UNRESOLVED\tCODON_CONTEXT_UNRESOLVED\tCODON_CONTEXT_UNRESOLVED\tUNCLASSIFIED\n"
                "chr9\t30\tA\tG\tDETECTED\tE255K\tVARIANT\tMISSENSE\tCURATED_RESISTANCE\n"
            )
            result = PHASER.load_detected_mutations(path)
            self.assertEqual([row["mutation"] for row in result], ["E255K"])

    def test_load_detected_mutations_handles_trailing_empty_fields(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "mutations.tsv"
            path.write_text(
                "chrom\tpos\tref\talt\tstatus\tmutation\textra1\textra2\n"
                "chr9\t30\tA\tG\tDETECTED\tE255K\t\t\n"
            )
            result = PHASER.load_detected_mutations(path)
            self.assertEqual(len(result), 1)
            self.assertEqual(result[0]["mutation"], "E255K")

    def test_missing_or_malformed_mutation_tsv_fails(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            with self.assertRaises(FileNotFoundError):
                PHASER.load_detected_mutations(tmp / "missing.tsv")
            path = tmp / "bad.tsv"
            path.write_text("chrom\tpos\nchr9\t50\n")
            with self.assertRaises(ValueError):
                PHASER.load_detected_mutations(path)
            aa3_path = tmp / "legacy_aa3.tsv"
            aa3_path.write_text(
                "chrom\tpos\tref\talt\tstatus\tmutation\taa_change\n"
                "chr9\t10\tA\tG\tDETECTED\tL387L\tp.Leu387Leu\n"
            )
            self.assertEqual(PHASER.load_detected_mutations(aa3_path), [])

            bad_bam = tmp / "bad.bam"
            bad_bam.write_text("not a BAM")
            empty_mutations = tmp / "empty.tsv"
            empty_mutations.write_text("chrom\tpos\tref\talt\tstatus\n")
            output = tmp / "out.tsv"
            with mock.patch.object(PHASER, "parse_args", return_value=type("Args", (), {
                "bam": bad_bam, "mutations_tsv": empty_mutations, "sample": "S", "output_tsv": output,
                "output_json": None, "min_spanning_reads": 10, "min_mapq": 1,
                "min_base_quality": 12, "threads": 1,
            })()):
                with self.assertRaises((ValueError, OSError)):
                    PHASER.main()
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
