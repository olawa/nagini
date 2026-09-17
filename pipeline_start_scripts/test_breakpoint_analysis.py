import importlib.util
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("breakpoint_analysis.py")
SPEC = importlib.util.spec_from_file_location("breakpoint_analysis", SCRIPT)
BREAKPOINTS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(BREAKPOINTS)


class SampleNameTest(unittest.TestCase):
    def test_a_nagini_target_bam_keeps_its_type_suffix(self):
        """A run can carry both N and T of one sample, so the type stays."""

        self.assertEqual(
            BREAKPOINTS.sample_name_for("alignment/target/gustaf1_T.target.bam", None),
            "gustaf1_T",
        )

    def test_the_old_sam_conventions_still_resolve(self):
        for path, expected in (
            ("sampleA_chr9_to_chr20.sam", "sampleA"),
            ("sampleA_chr20_to_chr9.sam", "sampleA"),
            ("/tmp/gustaf3.sam", "gustaf3"),
            ("x/gustaf4_T.target.cram", "gustaf4_T"),
        ):
            with self.subTest(path=path):
                self.assertEqual(BREAKPOINTS.sample_name_for(path, None), expected)

    def test_an_explicit_regex_wins(self):
        self.assertEqual(
            BREAKPOINTS.sample_name_for("run1/DNA2025-1001_T.target.bam", r"(DNA\d{4}-\d+)"),
            "DNA2025-1001",
        )


class ClusteringTest(unittest.TestCase):
    def test_both_coordinates_must_be_within_the_window(self):
        """A shared chr9 position is not a shared breakpoint on its own."""

        records = [
            (1000, 5000, "s1", "a"),
            (1003, 5002, "s1", "b"),   # within 5 bp on both -> same cluster
            (1002, 9000, "s1", "c"),   # close on chr9 only -> its own cluster
        ]
        groups = BREAKPOINTS.cluster_records(records, window=5)
        sizes = sorted(len(g) for g in groups)
        self.assertEqual(sizes, [1, 2])

    def test_a_read_seen_in_both_directions_is_one_observation(self):
        """The supplementary alignment carries the mirror of the same link.

        parse_alignment_file skips flag 0x800, so a fragment contributes once
        however many of its alignments the region filter returned.
        """

        import tempfile

        sam = (
            "@HD\tVN:1.6\tSO:coordinate\n"
            "@SQ\tSN:chr9\tLN:150617247\n"
            "@SQ\tSN:chr20\tLN:66210255\n"
            "r1\t0\tchr9\t40000100\t60\t60M\t*\t0\t0\t*\t*\tSA:Z:chr20,20000200,+,60M,60,0;\n"
            "r1\t2048\tchr20\t20000200\t60\t60M\t*\t0\t0\t*\t*\tSA:Z:chr9,40000100,+,60M,60,0;\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "one.sam"
            path.write_text(sam)
            links = list(
                BREAKPOINTS.parse_alignment_file(str(path), "chr9", "chr20", [])
            )
        self.assertEqual(links, [(40000100, 20000200, "r1")])


if __name__ == "__main__":
    unittest.main()
