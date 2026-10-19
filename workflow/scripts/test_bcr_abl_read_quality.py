import csv
import importlib.util
import math
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pysam

import bcr_abl_read_quality as rq

HERE = Path(__file__).parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


caller = load("bcr_abl_branford_caller")
phaser = load("bcr_abl_phaser")


def fake_aln(name="r", qs=None, quals=None, read_length=None, supplementary=False, start=0):
    tags = {} if qs is None else {"qs": qs}
    return SimpleNamespace(
        query_name=name, is_supplementary=supplementary, reference_start=start, is_read1=False,
        has_tag=lambda tag: tag in tags, get_tag=lambda tag: tags[tag],
        query_qualities=quals, infer_read_length=lambda: read_length if read_length is not None else len(quals or []),
    )


class FilterTest(unittest.TestCase):
    def test_mean_quality_is_the_error_probability_mean_not_the_phred_mean(self):
        self.assertAlmostEqual(rq.mean_read_quality([20] * 10), 20.0)
        self.assertAlmostEqual(rq.mean_read_quality([30] * 5 + [10] * 5), -10 * math.log10((5 * 0.001 + 5 * 0.1) / 10))
        self.assertIsNone(rq.mean_read_quality([]))

    def test_off_by_default_and_when_zero(self):
        flt = rq.ReadQualityFilter()
        self.assertFalse(flt.enabled)
        self.assertTrue(flt.passes(fake_aln(qs=1.0)))
        self.assertEqual(flt.summary(), "read quality filter off")

    def test_the_qs_tag_wins_over_the_qualities(self):
        flt = rq.ReadQualityFilter(15)
        self.assertFalse(flt.passes(fake_aln("low", qs=12.0, quals=[40] * 50)))
        self.assertTrue(flt.passes(fake_aln("ok", qs=18.0, quals=[5] * 50)))

    def test_without_the_tag_the_value_is_recomputed(self):
        flt = rq.ReadQualityFilter(15)
        self.assertFalse(flt.passes(fake_aln("low", quals=[10] * 50)))
        self.assertTrue(flt.passes(fake_aln("ok", quals=[20] * 50)))

    def test_a_hard_clipped_alignment_without_a_tag_is_kept_as_unknown(self):
        flt = rq.ReadQualityFilter(15)
        clipped = fake_aln("clip", quals=[5] * 50, read_length=100, supplementary=True)
        self.assertTrue(flt.passes(clipped))
        self.assertEqual(len(flt.unknown), 1)
        self.assertEqual(len(flt.dropped), 0)

    def test_a_read_dropped_on_one_alignment_is_dropped_on_all(self):
        flt = rq.ReadQualityFilter(15)
        self.assertFalse(flt.passes(fake_aln("r1", qs=5.0)))
        clipped = fake_aln("r1", quals=[40] * 20, read_length=100, supplementary=True, start=500)
        self.assertFalse(flt.passes(clipped))
        self.assertEqual(flt.dropped, {"r1"})
        self.assertIn("1 reads dropped", flt.summary())


def write_bam(path, records, length=300):
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr9", "LN": length}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, rec in enumerate(records):
            aln = pysam.AlignedSegment()
            aln.query_name = rec.get("name", f"read{i}")
            aln.query_sequence = rec["seq"]
            aln.query_qualities = rec.get("quals", [30] * len(rec["seq"]))
            aln.reference_id = 0
            aln.reference_start = 0
            aln.cigartuples = [(0, len(rec["seq"]))]
            aln.mapping_quality = 60
            aln.flag = 0
            if "qs" in rec:
                aln.set_tag("qs", rec["qs"], "f")
            out.write(aln)
    pysam.index(str(path))


def seq_with(base1, base2=None):
    seq = list("A" * 150)
    seq[49] = base1
    if base2:
        seq[99] = base2
    return "".join(seq)


def run_caller(tmp, bam, extra=()):
    bed = tmp / "panel.bed"
    bed.write_text("chr9\t49\t50\tNM:c.(949T>A)\tT\tA\tp.Phe317Ile\n")
    paths = {n: tmp / f"{n}.tsv" for n in ("all", "mut")}
    argv = ["caller", "--bam", str(bam), "--bed", str(bed), "--sample", "S",
            "--output-all", str(paths["all"]), "--output-mutations", str(paths["mut"]),
            "--output-vcf", str(tmp / "m.vcf"), "--min-depth", "10", "--min-alt-reads", "5", *extra]
    with patch.object(caller.sys, "argv", argv):
        caller.main()
    text = paths["all"].read_text()
    rows = list(csv.DictReader([l for l in text.splitlines() if l and not l.startswith("#")], delimiter="\t"))
    return rows[0], text


class CallerTest(unittest.TestCase):
    def records(self):
        # 60 good reference reads and 20 alt reads, the alt reads of low read quality
        good = [{"seq": seq_with("T"), "name": f"g{i}", "qs": 25.0} for i in range(60)]
        bad = [{"seq": seq_with("A"), "name": f"b{i}", "qs": 8.0} for i in range(20)]
        return good + bad

    def test_by_default_every_read_counts(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d); write_bam(tmp / "r.bam", self.records())
            row, text = run_caller(tmp, tmp / "r.bam")
            self.assertEqual(row["ad_alt"], "20")
            self.assertIn("# min_read_qs=0 reads_dropped=0", text)

    def test_a_read_quality_cut_off_removes_the_low_quality_reads_from_the_counts(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d); write_bam(tmp / "r.bam", self.records())
            row, text = run_caller(tmp, tmp / "r.bam", ["--min-read-qs", "10"])
            self.assertEqual(row["ad_alt"], "0")
            self.assertEqual(row["depth"], "60")
            self.assertIn("# min_read_qs=10 reads_dropped=20", text)

    def test_without_tags_the_quality_is_recomputed_from_the_bases(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            records = ([{"seq": seq_with("T"), "name": f"g{i}", "quals": [30] * 150} for i in range(60)]
                       + [{"seq": seq_with("A"), "name": f"b{i}", "quals": [10] * 150} for i in range(20)])
            write_bam(tmp / "r.bam", records)
            common = ["--min-base-quality", "0", "--min-mean-baseq", "0"]
            self.assertEqual(run_caller(tmp, tmp / "r.bam", common)[0]["ad_alt"], "20")
            self.assertEqual(run_caller(tmp, tmp / "r.bam", common + ["--min-read-qs", "15"])[0]["ad_alt"], "0")


class PhaserTest(unittest.TestCase):
    MUT1 = {"chrom": "chr9", "pos": 50, "ref": "T", "alt": "A", "mutation": "E255K"}
    MUT2 = {"chrom": "chr9", "pos": 100, "ref": "A", "alt": "T", "mutation": "T315I"}

    def phase(self, tmp, records, read_filter=None):
        path = tmp / "r.bam"; write_bam(path, records)
        with pysam.AlignmentFile(str(path), "rb") as bam:
            return phaser.phase_mutation_pair(bam, self.MUT1, self.MUT2, "S", 10, 1, 12, read_filter=read_filter)

    def test_low_read_quality_reads_leave_the_phasing_counts(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            records = ([{"seq": seq_with("A", "T"), "name": f"hq{i}", "qs": 25.0} for i in range(30)]
                       + [{"seq": seq_with("A", "T"), "name": f"lq{i}", "qs": 6.0} for i in range(20)])
            self.assertEqual(self.phase(tmp, records)["cis_reads_11"], 50)
            flt = rq.ReadQualityFilter(10)
            result = self.phase(tmp, records, read_filter=flt)
            self.assertEqual(result["cis_reads_11"], 30)
            self.assertEqual(len(flt.dropped), 20)
            self.assertEqual(result["overlapping_alignments"], 30)


if __name__ == "__main__":
    unittest.main()
