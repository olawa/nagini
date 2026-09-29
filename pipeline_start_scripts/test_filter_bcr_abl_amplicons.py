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

    def test_bam_units_handling(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            units_path = os.path.join(tmpdir, "units.tsv")
            rows = [
                {"sample": "sample1", "type": "T", "bam": "sample1.bam", "input_type": "bam"},
            ]
            fieldnames = ["sample", "type", "bam", "input_type"]
            FILTER_MOD.write_units_tsv(units_path, rows, fieldnames)
            parsed = FILTER_MOD.parse_units_tsv(units_path)
            self.assertEqual(len(parsed), 1)
            self.assertEqual(parsed[0]["bam"], "sample1.bam")

    def test_run_nanostream_amplicons_args(self):
        import unittest.mock as mock
        with mock.patch("subprocess.run") as mock_run:
            mock_run.return_value = mock.Mock(returncode=0, stdout="", stderr="")
            FILTER_MOD.run_nanostream_amplicons(
                nanostream_bin="nanostream",
                primers_tsv="primers.tsv",
                input_file="input.bam",
                output_json="out.json",
                output_fastq_prefix="out.fastq.gz",
                threads=4,
                split_by_amplicon=True,
                split_chimeras=False,
                min_qs=14.0,
                len_filter="1000-2500",
                mode="coords",
            )
            cmd = mock_run.call_args[0][0]
            self.assertIn("--min-qs", cmd)
            self.assertIn("14.0", cmd)
            self.assertIn("--len", cmd)
            self.assertIn("1000-2500", cmd)
            self.assertIn("--mode", cmd)
            self.assertIn("coords", cmd)
            self.assertIn("input.bam", cmd)


if __name__ == "__main__":
    unittest.main()
