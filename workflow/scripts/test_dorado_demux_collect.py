import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest import mock

import dorado_demux_collect as dc


def touch(root, *names):
    for name in names:
        path = Path(root) / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(name.encode())


class Recorder:
    def __init__(self, code=0):
        self.calls, self.code = [], code

    def __call__(self, command, **kwargs):
        self.calls.append((list(command), kwargs))
        return SimpleNamespace(returncode=self.code, stdout=b"@HD\tVN:1.6\n")


class FindBamsTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = self.tmp.name

    def tearDown(self):
        self.tmp.cleanup()

    def names(self, barcode):
        return [str(p.relative_to(self.root)) for p in dc.find_barcode_bams(self.root, barcode)]

    def test_flat_layout_named_by_kit_and_barcode(self):
        touch(self.root, "SQK-NBD114-24_barcode05.bam", "SQK-NBD114-24_barcode06.bam", "unclassified.bam")
        self.assertEqual(self.names("barcode05"), ["SQK-NBD114-24_barcode05.bam"])

    def test_directory_per_barcode_layout_and_bare_name(self):
        touch(self.root, "barcode05/reads_0.bam", "barcode05/reads_1.bam", "barcode06/reads_0.bam")
        self.assertEqual(self.names("barcode05"), ["barcode05/reads_0.bam", "barcode05/reads_1.bam"])
        touch(self.root, "barcode07.bam")
        self.assertEqual(self.names("barcode07"), ["barcode07.bam"])

    def test_a_barcode_is_not_a_prefix_or_suffix_of_another(self):
        touch(self.root, "kit_barcode10.bam", "kit_barcode100.bam", "kit_xbarcode1.bam", "kit_barcode1.bam")
        self.assertEqual(self.names("barcode1"), ["kit_barcode1.bam"])
        self.assertEqual(self.names("barcode10"), ["kit_barcode10.bam"])

    def test_zero_padding_matters(self):
        touch(self.root, "kit_barcode005.bam", "kit_barcode05.bam")
        self.assertEqual(self.names("barcode05"), ["kit_barcode05.bam"])

    def test_non_bam_files_are_ignored(self):
        touch(self.root, "kit_barcode05.bam.bai", "kit_barcode05.txt")
        self.assertEqual(self.names("barcode05"), [])


class CollectTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.raw = str(self.root / "raw")
        self.out = str(self.root / "out" / "barcode05.bam")

    def tearDown(self):
        self.tmp.cleanup()

    def test_single_bam_is_hard_linked_and_replaces_an_old_output(self):
        touch(self.raw, "kit_barcode05.bam")
        Path(self.out).parent.mkdir(parents=True)
        Path(self.out).write_text("stale")
        self.assertEqual(dc.collect_barcode(self.raw, "barcode05", "b.bam", self.out, False), "linked")
        self.assertEqual(os.stat(self.out).st_ino, os.stat(Path(self.raw) / "kit_barcode05.bam").st_ino)

    def test_copy_when_a_hard_link_is_not_possible(self):
        touch(self.raw, "kit_barcode05.bam")
        with mock.patch("os.link", side_effect=OSError("cross-device")):
            self.assertEqual(dc.collect_barcode(self.raw, "barcode05", "b.bam", self.out, False), "linked")
        self.assertEqual(Path(self.out).read_bytes(), b"kit_barcode05.bam")

    def test_several_bams_are_concatenated_in_name_order(self):
        touch(self.raw, "barcode05/b.bam", "barcode05/a.bam")
        rec = Recorder()
        self.assertEqual(dc.collect_barcode(self.raw, "barcode05", "b.bam", self.out, False, run=rec), "concatenated")
        command = rec.calls[0][0]
        self.assertEqual(command[:4], ["samtools", "cat", "-o", self.out])
        self.assertEqual([Path(p).name for p in command[4:]], ["a.bam", "b.bam"])

    def test_a_failing_cat_is_an_error(self):
        touch(self.raw, "barcode05/a.bam", "barcode05/b.bam")
        with self.assertRaises(dc.DemuxError):
            dc.collect_barcode(self.raw, "barcode05", "b.bam", self.out, False, run=Recorder(code=1))

    def test_no_reads_is_an_error_unless_empty_is_allowed(self):
        touch(self.raw, "kit_barcode06.bam")
        with self.assertRaises(dc.DemuxError) as raised:
            dc.collect_barcode(self.raw, "barcode05", "b.bam", self.out, False)
        self.assertIn("barcode05", str(raised.exception))
        rec = Recorder()
        self.assertEqual(dc.collect_barcode(self.raw, "barcode05", "b.bam", self.out, True, run=rec), "empty")
        self.assertEqual([c[0][:3] for c in rec.calls], [["samtools", "view", "-H"], ["samtools", "view", "-b"]])
        self.assertEqual(rec.calls[1][1]["input"], b"@HD\tVN:1.6\n")


class SummaryTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_summary_is_copied_and_the_output_directory_created(self):
        touch(self.root / "raw", "barcoding_summary.txt")
        out = self.root / "deep" / "demux_summary.tsv"
        dc.collect_summary(str(self.root / "raw"), str(out))
        self.assertEqual(out.read_bytes(), b"barcoding_summary.txt")

    def test_a_barcoding_summary_wins_over_other_summaries(self):
        touch(self.root / "raw", "a_summary.tsv", "z_barcoding_summary.txt", "reads.bam")
        self.assertEqual(dc.find_summary(str(self.root / "raw")).name, "z_barcoding_summary.txt")

    def test_the_choice_is_stable_when_nothing_says_barcoding(self):
        touch(self.root / "raw", "b_SUMMARY.TSV", "a_summary.txt")
        self.assertEqual(dc.find_summary(str(self.root / "raw")).name, "a_summary.txt")

    def test_no_summary_is_an_error(self):
        touch(self.root / "raw", "kit_barcode05.bam")
        with self.assertRaises(dc.DemuxError):
            dc.find_summary(str(self.root / "raw"))


class CliTest(unittest.TestCase):
    def test_exit_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            self.assertEqual(dc.main(["summary", "--raw-dir", tmp, "--output", f"{tmp}/s.tsv"]), 1)
            touch(tmp, "kit_barcode05.bam")
            self.assertEqual(dc.main(["barcode", "--raw-dir", tmp, "--barcode", "barcode05",
                                      "--basecalled", "b.bam", "--output", f"{tmp}/o/b.bam"]), 0)
            self.assertEqual(dc.main(["barcode", "--raw-dir", tmp, "--barcode", "barcode09",
                                      "--basecalled", "b.bam", "--output", f"{tmp}/o/c.bam"]), 1)


if __name__ == "__main__":
    unittest.main()
