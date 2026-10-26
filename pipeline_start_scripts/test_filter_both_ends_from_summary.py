#!/usr/bin/env python3
"""Tests for filter_both_ends_from_summary.py."""

import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import filter_both_ends_from_summary as fmod


class TestFilterBothEndsFromSummary(unittest.TestCase):
    def test_normalize_barcode(self):
        self.assertEqual(fmod.normalize_barcode("barcode01"), "barcode01")
        self.assertEqual(fmod.normalize_barcode("barcode1"), "barcode01")
        self.assertEqual(fmod.normalize_barcode("BC05"), "barcode05")
        self.assertEqual(fmod.normalize_barcode("bc12"), "barcode12")

    def test_parse_summary_file(self):
        content = (
            "read_id\tbarcode_arrangement\tmean_qscore_template\tsequence_length_template\tbarcode_front_score\tbarcode_rear_score\n"
            "read1\tbarcode01\t15.5\t1500\t85.0\t80.0\n"   # Pass all
            "read2\tbarcode01\t12.0\t1500\t85.0\t80.0\n"   # Fail QV (<14)
            "read3\tbarcode01\t16.0\t500\t85.0\t80.0\n"    # Fail len (<800)
            "read4\tbarcode01\t16.0\t3000\t85.0\t80.0\n"   # Fail len (>2500)
            "read5\tbarcode01\t16.0\t1500\t85.0\t0.0\n"    # Fail rear (single end)
            "read6\tbarcode02\t18.0\t1200\t90.0\t95.0\n"   # Pass all barcode02
            "read7\tunclassified\t18.0\t1200\t90.0\t95.0\n" # Unclassified
        )
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as f:
            f.write(content)
            temp_path = Path(f.name)

        try:
            reads_by_bc, stats = fmod.parse_summary_file(
                temp_path, min_qv=14.0, min_len=800, max_len=2500, min_barcode_score=60.0
            )
            self.assertEqual(len(reads_by_bc["barcode01"]), 1)
            self.assertIn("read1", reads_by_bc["barcode01"])
            self.assertEqual(len(reads_by_bc["barcode02"]), 1)
            self.assertIn("read6", reads_by_bc["barcode02"])

            self.assertEqual(stats["barcode01"]["total"], 5)
            self.assertEqual(stats["barcode01"]["pass_all"], 1)
            self.assertEqual(stats["barcode01"]["pass_both_ends"], 4)
        finally:
            if temp_path.is_file():
                temp_path.unlink()


if __name__ == "__main__":
    unittest.main()
