import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("pad_bed.py")
SPEC = importlib.util.spec_from_file_location("pad_bed", SCRIPT)
PAD_BED = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PAD_BED)


class PadBedTest(unittest.TestCase):
    def test_pad_single_interval(self):
        intervals = [(100000, 200000, "GENE1")]
        padded = PAD_BED.pad_and_merge_intervals(intervals, padding_bp=25000)
        self.assertEqual(padded, [(75000, 225000, ["GENE1"])])

    def test_start_does_not_go_below_zero(self):
        intervals = [(10000, 50000, "GENE2")]
        padded = PAD_BED.pad_and_merge_intervals(intervals, padding_bp=25000)
        self.assertEqual(padded, [(0, 75000, ["GENE2"])])

    def test_merge_overlapping_intervals_after_padding(self):
        # Two intervals separated by 40 kb, with 25 kb padding on each side (+/- 25kb), they overlap by 10 kb
        intervals = [
            (100000, 200000, "REGION1"),
            (240000, 300000, "REGION2"),
        ]
        padded = PAD_BED.pad_and_merge_intervals(intervals, padding_bp=25000)
        self.assertEqual(len(padded), 1)
        self.assertEqual(padded[0][0], 75000)
        self.assertEqual(padded[0][1], 325000)
        self.assertEqual(padded[0][2], ["REGION1", "REGION2"])

    def test_pad_and_merge_bed_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            in_bed = Path(tmpdir) / "in.bed"
            out_bed = Path(tmpdir) / "out.bed"
            in_bed.write_text(
                "chr9\t38000000\t70000000\tREGION_A\n"
                "chr20\t100000\t500000\tREGION_B\n"
            )
            PAD_BED.pad_and_merge_bed(str(in_bed), str(out_bed), padding_bp=25000)
            lines = out_bed.read_text().strip().split("\n")
            self.assertEqual(len(lines), 2)
            self.assertEqual(lines[0], "chr20\t75000\t525000\tREGION_B")
            self.assertEqual(lines[1], "chr9\t37975000\t70025000\tREGION_A")


if __name__ == "__main__":
    unittest.main()
