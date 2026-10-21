#!/usr/bin/env python3
"""Summarize job execution times, memory usage, and concurrency from Snakemake logs.

Parses Snakemake run logs (e.g. .snakemake/log/*.snakemake.log) and optional
benchmarks directories to report per-rule execution time, average durations,
peak memory, and workflow parallelization efficiency.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime
import os
from pathlib import Path
import re
import sys
from typing import Any, Dict, List, Optional, Tuple


TIMESTAMP_PATTERNS = [
    "%a %b %d %H:%M:%S %Y",      # e.g. Wed Oct 21 13:44:14 2026
    "%Y-%m-%d %H:%M:%S",         # e.g. 2026-10-21 13:44:14
    "%Y-%m-%dT%H:%M:%S.%f",      # ISO with microsecond
    "%Y-%m-%dT%H:%M:%S",         # ISO without microsecond
]


def parse_timestamp(text: str) -> Optional[datetime]:
    clean = text.strip("[] \t\r\n")
    for fmt in TIMESTAMP_PATTERNS:
        try:
            return datetime.strptime(clean, fmt)
        except ValueError:
            pass
    return None


def format_duration(seconds: float) -> str:
    if seconds < 0:
        return "0s"
    if seconds < 1.0:
        return f"{seconds:.2f}s"
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes = int(seconds // 60)
    rem_seconds = int(seconds % 60)
    if minutes < 60:
        return f"{minutes}m {rem_seconds:02d}s"
    hours = int(minutes // 60)
    rem_minutes = int(minutes % 60)
    return f"{hours}h {rem_minutes:02d}m {rem_seconds:02d}s"


def format_mb(mb: float) -> str:
    if mb <= 0:
        return "-"
    if mb >= 1024:
        return f"{mb / 1024:.2f} GB"
    return f"{mb:.1f} MB"


def find_latest_snakemake_log(search_dir: Path) -> Optional[Path]:
    log_dir = search_dir / ".snakemake" / "log"
    if not log_dir.is_dir():
        log_dir = search_dir / "log"
    if not log_dir.is_dir():
        return None
    logs = [p for p in log_dir.glob("*.snakemake.log") if p.is_file()]
    if not logs:
        logs = [p for p in log_dir.glob("*.log") if p.is_file()]
    if not logs:
        return None
    return max(logs, key=lambda p: p.stat().st_mtime)


def parse_snakemake_log(log_path: Path) -> Dict[str, Any]:
    if not log_path.is_file():
        raise FileNotFoundError(f"Log file not found: {log_path}")

    last_ts: Optional[datetime] = None
    first_ts: Optional[datetime] = None
    finish_ts: Optional[datetime] = None

    running_jobs: Dict[int, Dict[str, Any]] = {}
    rule_jobs: Dict[str, List[Dict[str, Any]]] = {}
    completed_jobs_count = 0
    failed_jobs_count = 0

    ts_regex = re.compile(r"^\[([A-Za-z]{3}\s+[A-Za-z]{3}\s+\d+\s+\d{2}:\d{2}:\d{2}\s+\d{4}|\d{4}-\d{2}-\d{2}[ T]\d{2}:\d{2}:\d{2})\]")
    job_start_regex = re.compile(r"^Job\s+(\d+):\s*([A-Za-z0-9_.-]+)")
    job_finish_regex = re.compile(r"^Finished job\s+(\d+)\.")
    job_error_regex = re.compile(r"Error in rule\s+([A-Za-z0-9_.-]+):")

    with log_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line_str = line.strip()

            # Check timestamp line
            m_ts = ts_regex.match(line_str)
            if m_ts:
                parsed_ts = parse_timestamp(m_ts.group(1))
                if parsed_ts:
                    last_ts = parsed_ts
                    if first_ts is None:
                        first_ts = parsed_ts
                    finish_ts = parsed_ts
                continue

            # Check job start line
            m_start = job_start_regex.match(line_str)
            if m_start:
                job_id = int(m_start.group(1))
                rule_name = m_start.group(2)
                running_jobs[job_id] = {
                    "job_id": job_id,
                    "rule": rule_name,
                    "start": last_ts,
                    "line": line_str,
                }
                continue

            # Check job finished line
            m_fin = job_finish_regex.match(line_str)
            if m_fin:
                job_id = int(m_fin.group(1))
                completed_jobs_count += 1
                if job_id in running_jobs:
                    job_info = running_jobs.pop(job_id)
                    rule_name = job_info["rule"]
                    start_time = job_info["start"]
                    end_time = last_ts or finish_ts
                    duration = (end_time - start_time).total_seconds() if (start_time and end_time) else 0.0
                    if duration < 0:
                        duration = 0.0
                    job_info["end"] = end_time
                    job_info["duration"] = duration
                    rule_jobs.setdefault(rule_name, []).append(job_info)
                continue

            # Check error line
            if job_error_regex.search(line_str):
                failed_jobs_count += 1

    wall_clock = (finish_ts - first_ts).total_seconds() if (first_ts and finish_ts) else 0.0

    return {
        "log_path": log_path,
        "first_ts": first_ts,
        "finish_ts": finish_ts,
        "wall_clock_seconds": wall_clock,
        "completed_jobs": completed_jobs_count,
        "failed_jobs": failed_jobs_count,
        "rule_jobs": rule_jobs,
    }


def parse_benchmark_files(benchmark_dir: Path) -> Dict[str, Dict[str, float]]:
    """Parse Snakemake benchmark TSVs for max_rss (MB) and cpu_time."""
    bench_data: Dict[str, List[Dict[str, float]]] = {}
    if not benchmark_dir.is_dir():
        return {}

    for tsv_path in benchmark_dir.rglob("*.tsv"):
        rule_name = tsv_path.parent.name
        try:
            with tsv_path.open("r", encoding="utf-8", errors="replace") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                for row in reader:
                    try:
                        s_time = float(row.get("s", row.get("cpu_time", 0.0)) or 0.0)
                        max_rss = float(row.get("max_rss", 0.0) or 0.0)
                        bench_data.setdefault(rule_name, []).append({
                            "s": s_time,
                            "max_rss": max_rss,
                        })
                    except (ValueError, TypeError):
                        pass
        except Exception:
            pass

    summary: Dict[str, Dict[str, float]] = {}
    for rname, rows in bench_data.items():
        if rows:
            max_rss_vals = [r["max_rss"] for r in rows if r["max_rss"] > 0]
            summary[rname] = {
                "max_rss_mb": max(max_rss_vals) if max_rss_vals else 0.0,
                "avg_rss_mb": sum(max_rss_vals) / len(max_rss_vals) if max_rss_vals else 0.0,
            }
    return summary


def aggregate_stats(
    parsed_log: Dict[str, Any],
    benchmarks: Optional[Dict[str, Dict[str, float]]] = None,
    sort_by: str = "total",
) -> List[Dict[str, Any]]:
    rule_jobs: Dict[str, List[Dict[str, Any]]] = parsed_log.get("rule_jobs", {})
    benchmarks = benchmarks or {}

    total_cpu_seconds = sum(
        sum(j["duration"] for j in jobs)
        for jobs in rule_jobs.values()
    )

    rows: List[Dict[str, Any]] = []
    for rule, jobs in rule_jobs.items():
        durations = [j["duration"] for j in jobs]
        count = len(durations)
        total_d = sum(durations)
        mean_d = total_d / count if count > 0 else 0.0
        min_d = min(durations) if durations else 0.0
        max_d = max(durations) if durations else 0.0
        pct = (total_d / total_cpu_seconds * 100.0) if total_cpu_seconds > 0 else 0.0

        b_info = benchmarks.get(rule, {})
        max_rss = b_info.get("max_rss_mb", 0.0)

        rows.append({
            "rule": rule,
            "count": count,
            "total_seconds": total_d,
            "mean_seconds": mean_d,
            "min_seconds": min_d,
            "max_seconds": max_d,
            "percent_cpu": pct,
            "max_rss_mb": max_rss,
        })

    sort_keys = {
        "total": lambda r: r["total_seconds"],
        "mean": lambda r: r["mean_seconds"],
        "max": lambda r: r["max_seconds"],
        "min": lambda r: r["min_seconds"],
        "count": lambda r: r["count"],
        "rule": lambda r: r["rule"].lower(),
        "name": lambda r: r["rule"].lower(),
        "memory": lambda r: r["max_rss_mb"],
    }
    key_func = sort_keys.get(sort_by, sort_keys["total"])
    reverse = sort_by not in ("rule", "name")
    rows.sort(key=key_func, reverse=reverse)
    return rows


def print_table(
    parsed_log: Dict[str, Any],
    rows: List[Dict[str, Any]],
    has_benchmarks: bool = False,
) -> None:
    wall_clock = parsed_log["wall_clock_seconds"]
    total_cpu = sum(r["total_seconds"] for r in rows)
    first_ts = parsed_log["first_ts"]
    finish_ts = parsed_log["finish_ts"]
    log_path = parsed_log["log_path"]
    speedup = (total_cpu / wall_clock) if wall_clock > 0 else 1.0

    print("=" * 96)
    print("NAGINI SNAKEMAKE WORKFLOW RUNTIME SUMMARY")
    print("=" * 96)
    print(f"Log file:        {log_path}")
    if first_ts and finish_ts:
        print(f"Start time:      {first_ts.strftime('%Y-%m-%d %H:%M:%S')}")
        print(f"Finish time:     {finish_ts.strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Wall-clock time: {format_duration(wall_clock)} ({wall_clock:.0f}s)")
    print(f"Total CPU time:  {format_duration(total_cpu)} (cumulative across jobs)")
    print(f"Speedup ratio:   {speedup:.2f}x (effective core concurrency)")
    print(f"Jobs completed:  {parsed_log['completed_jobs']} (Failed: {parsed_log['failed_jobs']})")
    print("=" * 96)

    rule_col_w = max(32, max((len(r["rule"]) for r in rows), default=32))
    header = (
        f"{'Rule Name':<{rule_col_w}} "
        f"{'Count':>6} "
        f"{'Total Time':>12} "
        f"{'Mean/Job':>10} "
        f"{'Min':>8} "
        f"{'Max':>8} "
        f"{'% Total':>8}"
    )
    if has_benchmarks:
        header += f" {'Max RSS':>10}"

    print(header)
    print("-" * len(header))

    for r in rows:
        line = (
            f"{r['rule']:<{rule_col_w}} "
            f"{r['count']:>6d} "
            f"{format_duration(r['total_seconds']):>12} "
            f"{format_duration(r['mean_seconds']):>10} "
            f"{format_duration(r['min_seconds']):>8} "
            f"{format_duration(r['max_seconds']):>8} "
            f"{r['percent_cpu']:>7.1f}%"
        )
        if has_benchmarks:
            line += f" {format_mb(r['max_rss_mb']):>10}"
        print(line)

    print("-" * len(header))
    summary_line = (
        f"{'TOTAL':<{rule_col_w}} "
        f"{sum(r['count'] for r in rows):>6d} "
        f"{format_duration(total_cpu):>12} "
        f"{format_duration(total_cpu / sum(r['count'] for r in rows) if rows else 0):>10} "
        f"{'':>8} "
        f"{'':>8} "
        f"{'100.0%':>8}"
    )
    print(summary_line)
    print("=" * len(header))


def export_tsv(rows: List[Dict[str, Any]], out_path: Path, has_benchmarks: bool = False) -> None:
    headers = [
        "rule", "count", "total_seconds", "total_duration",
        "mean_seconds", "mean_duration", "min_seconds", "max_seconds", "percent_cpu"
    ]
    if has_benchmarks:
        headers.append("max_rss_mb")

    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(headers)
        for r in rows:
            row_vals = [
                r["rule"],
                r["count"],
                f"{r['total_seconds']:.2f}",
                format_duration(r["total_seconds"]),
                f"{r['mean_seconds']:.2f}",
                format_duration(r["mean_seconds"]),
                f"{r['min_seconds']:.2f}",
                f"{r['max_seconds']:.2f}",
                f"{r['percent_cpu']:.2f}",
            ]
            if has_benchmarks:
                row_vals.append(f"{r['max_rss_mb']:.1f}")
            writer.writerow(row_vals)
    print(f"[TSV] Saved stats table to: {out_path}", file=sys.stderr)


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Summarize job runtimes and statistics from a Snakemake log file."
    )
    parser.add_argument(
        "log_or_dir",
        nargs="?",
        default=".",
        help="Path to Snakemake log file (*.snakemake.log) or analysis directory (default: current dir)",
    )
    parser.add_argument(
        "--benchmarks",
        type=Path,
        default=None,
        help="Optional path to benchmarks directory (default: auto-detect if benchmarks/ exists)",
    )
    parser.add_argument(
        "--tsv",
        type=Path,
        default=None,
        help="Optional path to output stats table as TSV",
    )
    parser.add_argument(
        "--sort-by",
        choices=["total", "mean", "max", "min", "count", "rule", "name", "memory"],
        default="total",
        help="Sort table by column (default: total)",
    )
    args = parser.parse_args(argv)

    target_path = Path(args.log_or_dir).resolve()
    if target_path.is_dir():
        log_file = find_latest_snakemake_log(target_path)
        if not log_file:
            sys.exit(
                f"Error: Could not find any .snakemake.log in '{target_path}' or '{target_path}/.snakemake/log'."
            )
        base_dir = target_path
    else:
        log_file = target_path
        base_dir = target_path.parent
        if ".snakemake" in base_dir.parts:
            # Climb up to analysis dir
            while base_dir.name in (".snakemake", "log") and len(base_dir.parts) > 1:
                base_dir = base_dir.parent

    parsed = parse_snakemake_log(log_file)

    # Check benchmarks
    bench_dir = args.benchmarks
    if bench_dir is None:
        cand_bench = base_dir / "benchmarks"
        if cand_bench.is_dir():
            bench_dir = cand_bench
    bench_data = parse_benchmark_files(bench_dir) if bench_dir else {}

    rows = aggregate_stats(parsed, bench_data, sort_by=args.sort_by)
    print_table(parsed, rows, has_benchmarks=bool(bench_data))

    if args.tsv:
        export_tsv(rows, args.tsv, has_benchmarks=bool(bench_data))

    return 0


if __name__ == "__main__":
    sys.exit(main())
