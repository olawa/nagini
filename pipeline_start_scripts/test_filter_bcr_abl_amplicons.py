import csv
import gzip
import importlib.util
import os
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("filter_bcr_abl_amplicons.py")
SPEC = importlib.util.spec_from_file_location("filter_bcr_abl_amplicons", SCRIPT)
FILTER_MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(FILTER_MOD)


class FilterBcrAblAmpliconsTest(unittest.TestCase):
    def test_find_default_primers(self):
        primers_file = FILTER_MOD.find_default_primers()
        self.assertTrue(os.path.isfile(primers_file))
        with open(primers_file, "r") as f:
            content = f.read()
        self.assertIn("BCR-ABL1major_fwd", content)
        self.assertIn("BCR-ABL1_rev", content)

    def test_units_parse_and_write(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            units_path = os.path.join(tmpdir, "units.tsv")
            rows = [
                {"sample": "sample1", "type": "T", "fastq1": "fq1.fastq.gz"},
                {"sample": "sample2", "type": "T", "fastq1": "fq2.fastq.gz"},
            ]
            fieldnames = ["sample", "type", "fastq1"]
            FILTER_MOD.write_units_tsv(units_path, rows, fieldnames)
            parsed = FILTER_MOD.parse_units_tsv(units_path)
            self.assertEqual(len(parsed), 2)
            self.assertEqual(parsed[0]["sample"], "sample1")
            self.assertEqual(parsed[1]["fastq1"], "fq2.fastq.gz")


if __name__ == "__main__":
    unittest.main()
