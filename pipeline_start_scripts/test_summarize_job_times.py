#!/usr/bin/env python3
import sys
import tempfile
import unittest
from datetime import datetime
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

from summarize_job_times import (
    format_duration,
    parse_timestamp,
    parse_snakemake_log,
    parse_benchmark_files,
    aggregate_stats,
    main,
)


class TestSummarizeJobTimes(unittest.TestCase):
    def test_parse_timestamp(self):
        dt1 = parse_timestamp("[Wed Oct 21 13:44:14 2026]")
        self.assertIsNotNone(dt1)
        self.assertEqual(dt1.hour, 13)
        self.assertEqual(dt1.minute, 44)
        self.assertEqual(dt1.second, 14)

        dt2 = parse_timestamp("2026-10-21 13:44:14")
        self.assertIsNotNone(dt2)
        self.assertEqual(dt2.year, 2026)

    def test_format_duration(self):
        self.assertEqual(format_duration(0.45), "0.45s")
        self.assertEqual(format_duration(14.2), "14.2s")
        self.assertEqual(format_duration(105), "1m 45s")
        self.assertEqual(format_duration(3665), "1h 01m 05s")

    def test_parse_snakemake_log_interleaved_jobs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            log_path = tmp / "test.snakemake.log"
            log_path.write_text(
                "[Wed Oct 21 13:44:10 2026]\n"
                "Job 1: alignment_minimap2_fastq_align: start 1\n"
                "[Wed Oct 21 13:44:15 2026]\n"
                "Job 2: bcr_abl_branford_calling: start 2\n"
                "[Wed Oct 21 13:44:20 2026]\n"
                "Finished job 2.\n"
                "1 of 2 steps done\n"
                "[Wed Oct 21 13:46:10 2026]\n"
                "Finished job 1.\n"
                "2 of 2 steps done\n"
            )

            res = parse_snakemake_log(log_path)
            self.assertEqual(res["completed_jobs"], 2)
            self.assertEqual(res["wall_clock_seconds"], 120.0)

            rule_jobs = res["rule_jobs"]
            self.assertIn("alignment_minimap2_fastq_align", rule_jobs)
            self.assertIn("bcr_abl_branford_calling", rule_jobs)

            # Minimap: 13:44:10 to 13:46:10 = 120s
            self.assertEqual(rule_jobs["alignment_minimap2_fastq_align"][0]["duration"], 120.0)
            # Branford: 13:44:15 to 13:44:20 = 5s
            self.assertEqual(rule_jobs["bcr_abl_branford_calling"][0]["duration"], 5.0)

            stats = aggregate_stats(res)
            self.assertEqual(len(stats), 2)
            self.assertEqual(stats[0]["rule"], "alignment_minimap2_fastq_align")
            self.assertEqual(stats[0]["total_seconds"], 120.0)
            self.assertEqual(stats[1]["rule"], "bcr_abl_branford_calling")
            self.assertEqual(stats[1]["total_seconds"], 5.0)

    def test_parse_benchmark_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bench_dir = tmp / "benchmarks" / "bcr_abl_vardictrust"
            bench_dir.mkdir(parents=True)
            tsv = bench_dir / "sample1.tsv"
            tsv.write_text("s\th:m:s\tmax_rss\tmax_vms\tmax_uss\tmax_pss\tio_in\tio_out\tmean_load\tcpu_time\n45.2\t0:00:45\t512.4\t1000\t500\t500\t0\t0\t1.5\t60.1\n")

            benchmarks = parse_benchmark_files(tmp / "benchmarks")
            self.assertIn("bcr_abl_vardictrust", benchmarks)
            self.assertAlmostEqual(benchmarks["bcr_abl_vardictrust"]["max_rss_mb"], 512.4)

    def test_main_cli(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            log_dir = tmp / ".snakemake" / "log"
            log_dir.mkdir(parents=True)
            log_file = log_dir / "run.snakemake.log"
            log_file.write_text(
                "[Wed Oct 21 13:44:10 2026]\n"
                "Job 1: test_rule: run test\n"
                "[Wed Oct 21 13:44:20 2026]\n"
                "Finished job 1.\n"
            )
            out_tsv = tmp / "stats.tsv"
            ret = main([str(tmp), "--tsv", str(out_tsv)])
            self.assertEqual(ret, 0)
            self.assertTrue(out_tsv.exists())
            self.assertIn("test_rule", out_tsv.read_text())


if __name__ == "__main__":
    unittest.main()
