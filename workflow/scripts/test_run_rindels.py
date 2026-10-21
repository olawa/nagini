import argparse
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import run_rindels as rr


def done(code=0):
    return SimpleNamespace(returncode=code)


class Recorder:
    """A subprocess.run stand-in that records commands and can create files."""

    def __init__(self, code=0, creates=(), codes=None):
        self.calls, self.code, self.creates, self.codes = [], code, creates, codes or {}

    def __call__(self, command):
        self.calls.append(list(command))
        for path in self.creates:
            Path(path).write_text("x")
        return done(self.codes.get(command[0], self.code))


def have(*names):
    return lambda name: f"/bin/{name}" if name in names else None


class BinaryTest(unittest.TestCase):
    def test_binary_on_path_or_executable_file_is_accepted(self):
        rr.check_binary("rindels", which=have("rindels"))
        with tempfile.TemporaryDirectory() as tmp:
            exe = Path(tmp) / "rindels"
            exe.write_text("#!/bin/sh\n")
            exe.chmod(0o755)
            rr.check_binary(str(exe), which=have())

    def test_missing_binary_is_a_clear_error(self):
        with self.assertRaises(rr.RindelsError) as raised:
            rr.check_binary("/nope/rindels", which=have())
        self.assertIn("/nope/rindels", str(raised.exception))


class CnvTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.args = argparse.Namespace(
            binary="rindels", bam="a.bam", ref="r.fa", bins=str(root / "bins.bed"), bin_size=50000,
            vcf=str(root / "s.vcf"), ascat=str(root / "s.ascat.bed"), pdf=str(root / "s.pdf"),
            report=True, annotation="--gtf g.gtf --annotate --annot-coding-only",
            threads=4, extra="--ont --sv-lite",
        )

    def tearDown(self):
        self.tmp.cleanup()

    def test_command_order_is_the_historical_one(self):
        command = rr.cnv_command(self.args)
        self.assertEqual(command[:11], ["rindels", "-b", "a.bam", "-r", "r.fa", "--cnv-output", self.args.bins,
                                         "--cnv-bin-size", "50000", "-o", self.args.vcf])
        self.assertEqual(command[11:14], ["--report", "--report-path", self.args.pdf])
        self.assertEqual(command[14:16], ["-t", "4"])
        self.assertEqual(command[16:], ["--gtf", "g.gtf", "--annotate", "--annot-coding-only", "--ont", "--sv-lite"])

    def test_report_can_be_switched_off(self):
        self.args.report = False
        self.assertNotIn("--report", rr.cnv_command(self.args))

    def test_placeholders_when_rindels_skips_segmentation_and_report(self):
        rr.run_cnv(self.args, run=Recorder(), which=have("rindels"))
        self.assertTrue(Path(self.args.ascat).read_text().startswith("# ASCAT-style Heuristic Copy Number Segments"))
        self.assertEqual(Path(self.args.pdf).stat().st_size, 0)

    def test_real_outputs_are_left_alone(self):
        rec = Recorder(creates=[self.args.ascat, self.args.pdf])
        rr.run_cnv(self.args, run=rec, which=have("rindels"))
        self.assertEqual(Path(self.args.ascat).read_text(), "x")
        self.assertEqual(Path(self.args.pdf).read_text(), "x")

    def test_a_failing_rindels_run_fails_and_leaves_no_placeholders(self):
        with self.assertRaises(rr.RindelsError):
            rr.run_cnv(self.args, run=Recorder(code=2), which=have("rindels"))
        self.assertFalse(Path(self.args.ascat).exists())
        self.assertFalse(Path(self.args.pdf).exists())

    def test_missing_binary_stops_before_running(self):
        rec = Recorder()
        with self.assertRaises(rr.RindelsError):
            rr.run_cnv(self.args, run=rec, which=have())
        self.assertEqual(rec.calls, [])


class IdeogramTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        self.args = argparse.Namespace(
            binary="rindels", bins=str(self.root / "bins.bed"), parquet=str(self.root / "s.vcf.parquet"),
            baf=str(self.root / "s.baf.bed"), centromeres=str(self.root / "cen.bed"), gaps="",
            bin_size=50000, width=2000, height=2400, extra="", output=str(self.root / "out" / "i.png"),
        )
        self.messages = []

    def tearDown(self):
        self.tmp.cleanup()

    def test_nothing_optional_available(self):
        command = rr.ideogram_command(self.args, say=self.messages.append)
        self.assertEqual(command[:4], ["rindels", "ideogram", "--bed", self.args.bins])
        self.assertNotIn("--baf-bed", command)
        self.assertNotIn("--centromeres", command)
        self.assertNotIn("--gaps", command)
        self.assertEqual(len(self.messages), 3)
        self.assertEqual(command[command.index("--res") + 1], "50")
        self.assertEqual(command[-2:], ["-o", self.args.output])

    def test_everything_optional_available(self):
        for name in ("parquet", "baf", "centromeres"):
            Path(getattr(self.args, name)).write_text("data")
        self.args.gaps = str(self.root / "gaps.bed")
        Path(self.args.gaps).write_text("g")
        command = rr.ideogram_command(self.args, say=self.messages.append)
        self.assertEqual(command[2:4], ["--cnv", self.args.parquet])
        self.assertIn("--baf-bed", command)
        self.assertIn("--centromeres", command)
        self.assertIn("--gaps", command)
        self.assertNotIn("--bed", command)
        self.assertEqual(self.messages, [])

    def test_an_empty_sidecar_counts_as_missing(self):
        Path(self.args.parquet).write_text("")
        self.assertIn("--bed", rr.ideogram_command(self.args, say=lambda m: None))

    def test_output_directory_is_created_and_extra_is_split(self):
        self.args.extra = "--foo 1"
        rec = Recorder()
        rr.run_ideogram(self.args, run=rec, which=have("rindels"))
        self.assertTrue(Path(self.args.output).parent.is_dir())
        self.assertIn("--foo", rec.calls[0])


class BcrAblTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = Path(self.tmp.name)
        self.args = argparse.Namespace(
            binary="rindels", bam="a.bam", ref="r.fa", vcf=str(root / "s.rindels.vcf.gz"),
            haplotypes=str(root / "s.haplotypes.tsv"), genome="hg19", threads=2,
            extra="--ont --amplicon --bcr-abl",
        )

    def tearDown(self):
        self.tmp.cleanup()

    def test_command(self):
        self.assertEqual(
            rr.bcr_abl_command(self.args),
            ["rindels", "-b", "a.bam", "-r", "r.fa", "-o", self.args.vcf, "--ont", "--amplicon", "--bcr-abl",
             "--bcr-abl-genome", "hg19", "--write-haplotypes", self.args.haplotypes, "-t", "2"],
        )

    def test_vcf_is_indexed_with_bcftools_then_tabix(self):
        rec = Recorder(creates=[self.args.vcf])
        rr.run_bcr_abl(self.args, run=rec, which=have("rindels", "bcftools", "tabix"))
        self.assertEqual(rec.calls[-1], ["bcftools", "index", "-t", self.args.vcf])
        rec = Recorder(creates=[self.args.vcf])
        rr.run_bcr_abl(self.args, run=rec, which=have("rindels", "tabix"))
        self.assertEqual(rec.calls[-1][:3], ["tabix", "-p", "vcf"])

    def test_no_index_tool_or_empty_vcf_is_not_an_error(self):
        rec = Recorder(creates=[self.args.vcf])
        rr.run_bcr_abl(self.args, run=rec, which=have("rindels"))
        self.assertEqual(len(rec.calls), 1)
        Path(self.args.vcf).unlink()
        rec = Recorder()
        rr.run_bcr_abl(self.args, run=rec, which=have("rindels", "bcftools"))
        self.assertEqual(len(rec.calls), 1)

    def test_a_failed_index_warns_but_a_failed_run_fails(self):
        rec = Recorder(creates=[self.args.vcf], codes={"bcftools": 1})
        rr.run_bcr_abl(self.args, run=rec, which=have("rindels", "bcftools"))
        with self.assertRaises(rr.RindelsError):
            rr.run_bcr_abl(self.args, run=Recorder(codes={"rindels": 3}), which=have("rindels", "bcftools"))

    def test_ensure_compressed_vcf_compresses_raw_vcf_when_gz_is_empty(self):
        # Simulate rindels creating 0-byte vcf.gz and populating raw vcf
        gz_vcf = Path(self.args.vcf)
        gz_vcf.write_bytes(b"")  # 0 bytes
        raw_vcf = gz_vcf.with_suffix("")  # strip .gz -> .vcf
        raw_vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\nchr9\t100\n")

        rec = Recorder(creates=[])
        rr.ensure_compressed_vcf(str(gz_vcf), run=rec, which=have("bcftools"))
        # Verify bcftools view -Oz was called or gzip populated the gz file
        self.assertTrue(gz_vcf.exists())
        self.assertGreater(gz_vcf.stat().st_size, 0)


class CliTest(unittest.TestCase):
    def test_main_reports_a_missing_binary_with_status_one(self):
        with tempfile.TemporaryDirectory() as tmp:
            code = rr.main(["bcr-abl", "--binary", "/nope/rindels", "--bam", "a", "--ref", "r",
                            "--vcf", f"{tmp}/v.vcf.gz", "--haplotypes", f"{tmp}/h.tsv", "--genome", "hg19"])
        self.assertEqual(code, 1)


if __name__ == "__main__":
    unittest.main()
