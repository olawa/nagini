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
    def test_handoff_config_repoints_samples_and_units(self):
        """The copied config must read the filtered tables, not the originals."""
        import yaml

        with tempfile.TemporaryDirectory() as tmpdir:
            original = Path(tmpdir) / "original"
            original.mkdir()
            out_dir = Path(tmpdir) / "amplicon_filtered"
            out_dir.mkdir()
            config_path = original / "config.nagini.yaml"
            # prepare_nagini_inputs.py writes these as resolved absolute paths.
            config_path.write_text(yaml.safe_dump({
                "samples": str(original / "samples.tsv"),
                "units": str(original / "units.tsv"),
                "reference_build": "hg19",
                "read_filter": {"enabled": True, "min_qv": 20},
            }))

            dest = FILTER_MOD.write_handoff_config(config_path, out_dir)
            written = yaml.safe_load(dest.read_text())

            self.assertEqual(written["samples"], str((out_dir / "samples.tsv").resolve()))
            self.assertEqual(written["units"], str((out_dir / "units.tsv").resolve()))
            # Nothing may still point at the pre-filtering directory.
            self.assertNotIn(str(original), written["samples"])
            self.assertNotIn(str(original), written["units"])

    def test_handoff_config_disables_read_filter_and_keeps_the_rest(self):
        import yaml

        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir) / "out"
            out_dir.mkdir()
            config_path = Path(tmpdir) / "config.nagini.yaml"
            config_path.write_text(yaml.safe_dump({
                "samples": "/elsewhere/samples.tsv",
                "units": "/elsewhere/units.tsv",
                "reference_build": "hg19",
                "read_filter": {"enabled": True, "min_qv": 20},
                "bcr_abl": {"enabled": True},
            }))

            written = yaml.safe_load(
                FILTER_MOD.write_handoff_config(config_path, out_dir).read_text()
            )

            # Already-filtered amplicons must not be filtered again.
            self.assertIs(written["read_filter"]["enabled"], False)
            # ... but unrelated settings survive the round trip.
            self.assertEqual(written["read_filter"]["min_qv"], 20)
            self.assertEqual(written["reference_build"], "hg19")
            self.assertEqual(written["bcr_abl"], {"enabled": True})

    def test_handoff_config_adds_read_filter_when_absent(self):
        import yaml

        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir) / "out"
            out_dir.mkdir()
            config_path = Path(tmpdir) / "config.nagini.yaml"
            config_path.write_text(yaml.safe_dump({"reference_build": "hg19"}))

            written = yaml.safe_load(
                FILTER_MOD.write_handoff_config(config_path, out_dir).read_text()
            )
            self.assertIs(written["read_filter"]["enabled"], False)

    def test_handoff_config_fails_loudly_on_a_broken_config(self):
        """A verbatim copy would be exactly the bug this replaces."""
        with tempfile.TemporaryDirectory() as tmpdir:
            out_dir = Path(tmpdir) / "out"
            out_dir.mkdir()
            config_path = Path(tmpdir) / "config.nagini.yaml"
            config_path.write_text("samples: [unclosed\n  - broken: :\n")

            with self.assertRaises(Exception):
                FILTER_MOD.write_handoff_config(config_path, out_dir)
            self.assertFalse((out_dir / "config.nagini.yaml").exists())

    def test_nagini_root_is_this_checkout(self):
        root = FILTER_MOD.nagini_root()
        self.assertTrue((root / "workflow" / "Snakefile").is_file())
        self.assertTrue((root / "pipeline_start_scripts").is_dir())
        # No one machine's layout may be baked into the script.
        offenders = [
            line.strip()
            for line in SCRIPT.read_text().splitlines()
            if "/home/cgu-lab" in line and not line.lstrip().startswith("#")
        ]
        self.assertEqual(offenders, [])

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

    def test_is_valid_input_path(self):
        self.assertFalse(FILTER_MOD.is_valid_input_path(None))
        self.assertFalse(FILTER_MOD.is_valid_input_path(""))
        self.assertFalse(FILTER_MOD.is_valid_input_path("   "))
        self.assertFalse(FILTER_MOD.is_valid_input_path("UNUSED"))
        self.assertFalse(FILTER_MOD.is_valid_input_path("unused"))
        self.assertFalse(FILTER_MOD.is_valid_input_path("NONE"))
        self.assertFalse(FILTER_MOD.is_valid_input_path("."))
        self.assertFalse(FILTER_MOD.is_valid_input_path("NA"))
        self.assertFalse(FILTER_MOD.is_valid_input_path("null"))
        self.assertTrue(FILTER_MOD.is_valid_input_path("sample.bam"))
        self.assertTrue(FILTER_MOD.is_valid_input_path("/path/to/reads.fastq.gz"))

    def test_resolve_sample_input_ignores_unused_fastq(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            bam_path = run_dir / "ratio0.1_T.bam"
            bam_path.touch()

            row = {
                "sample": "ratio0.1",
                "type": "T",
                "input_type": "aligned_bam",
                "fastq1": "UNUSED",
                "fastq2": "UNUSED",
                "bam": str(bam_path),
            }
            resolved = FILTER_MOD.resolve_sample_input(row, run_dir, fq_col="fastq1", bam_col="bam")
            self.assertEqual(resolved, str(bam_path))

    def test_resolve_sample_input_relative_bam(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            bams_subdir = run_dir / "bams"
            bams_subdir.mkdir()
            bam_path = bams_subdir / "sample.bam"
            bam_path.touch()

            row = {
                "sample": "sample",
                "type": "T",
                "input_type": "aligned_bam",
                "fastq1": "UNUSED",
                "bam": "bams/sample.bam",
            }
            resolved = FILTER_MOD.resolve_sample_input(row, run_dir, fq_col="fastq1", bam_col="bam")
            self.assertEqual(resolved, str(bam_path))

    def test_resolve_sample_input_fallback_merged(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            merged_dir = run_dir / "alignment" / "merged"
            merged_dir.mkdir(parents=True)
            merged_bam = merged_dir / "sample1_T.bam"
            merged_bam.touch()

            row = {
                "sample": "sample1",
                "type": "T",
                "input_type": "aligned_bam",
                "fastq1": "UNUSED",
                "bam": "UNUSED",
            }
            resolved = FILTER_MOD.resolve_sample_input(row, run_dir, fq_col="fastq1", bam_col="bam")
            self.assertEqual(resolved, str(merged_bam))


if __name__ == "__main__":
    unittest.main()
