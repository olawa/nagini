"""The same QS/base filters must select the same observations in both tools."""

import tempfile
import unittest
from pathlib import Path

import pysam

import bcr_abl_branford_caller as caller
from bcr_abl_read_quality import ReadQualityFilter
try:
    import matplotlib  # noqa: F401  (the sweep and the plotter need it)
except ImportError:
    matplotlib = None

if matplotlib is not None:
    import sweep_bcr_abl_qs as sweep
else:
    sweep = None


@unittest.skipIf(matplotlib is None, "matplotlib is required for the QS sweep and background plots")
class SweepCallerParityTests(unittest.TestCase):
    def compare_counts(self, records, min_baseq):
        with tempfile.TemporaryDirectory() as directory:
            bam_path = Path(directory) / "reads.bam"
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr9", "LN": 300}]}
            with pysam.AlignmentFile(str(bam_path), "wb", header=header) as output:
                for record in records:
                    alignment = pysam.AlignedSegment()
                    alignment.query_name = record.get("name", "read")
                    alignment.query_sequence = record.get("base", "A") * 150
                    alignment.query_qualities = record.get("quals", [35] * 150)
                    alignment.reference_id = 0
                    alignment.reference_start = 0
                    alignment.mapping_quality = 60
                    alignment.flag = record.get("flag", 0)
                    alignment.cigartuples = record.get("cigar", [(0, 150)])
                    if "qs" in record:
                        alignment.set_tag("qs", record["qs"], "f")
                    output.write(alignment)
            pysam.index(str(bam_path))
            site = {
                "chrom": "chr9", "pos0": 49, "pos1": 50, "ref": "G", "alt": "A",
                "name": "NM:c.(763G>A)", "aa_change": "p.Glu255Lys", "mutation": "E255K",
            }
            _, counts = sweep.collect_qs_pileup(
                bam_path, [site], min_base_quality=min_baseq, min_mapping_quality=1,
                cumulative_thresholds=[20],
            )
            bin_def = sweep.make_bins(20, 20, "cumulative")[0]
            row = sweep.rows_for_bin([site], counts, bin_def, "S")[0]
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                expected = caller.count_position(
                    bam, "chr9", 50, 1000, min_baseq, 1,
                    read_filter=ReadQualityFilter(20),
                )
            for base in "ACGT":
                self.assertEqual(
                    int(row[f"count_{base}"]), expected[f"{base}_fwd"] + expected[f"{base}_rev"],
                    f"QS >= 20 count for {base}",
                )

    def test_unknown_quality_hard_clipped_supplementary_read_is_kept(self):
        self.compare_counts([
            {"flag": 2048, "cigar": [(5, 50), (0, 150)], "quals": [10] * 150},
        ], min_baseq=0)

    def test_low_qs_tag_takes_priority_over_high_base_qualities(self):
        self.compare_counts([{"qs": 10.0}], min_baseq=12)

    def test_best_alignment_is_chosen_before_base_quality_filtering(self):
        qualities = [35] * 150
        qualities[49] = 5
        self.compare_counts([
            {"qs": 30.0, "base": "G", "quals": qualities},
            {"qs": 30.0, "base": "A", "flag": 2048},
        ], min_baseq=12)


if __name__ == "__main__":
    unittest.main()
