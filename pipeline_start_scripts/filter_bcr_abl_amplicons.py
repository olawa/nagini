#!/usr/bin/env python3
"""Run nanostream amplicons filtering on an existing Nagini run directory.

This script:
1. Reads `units.tsv` and `samples.tsv` from an existing Nagini run directory.
2. Runs `nanostream amplicons` using BCR-ABL1 primers with `--split-by-amplicon`.
3. Collects stats (counts, lengths, Q-scores, primer dimers, off-target reads).
4. Sets up a new analysis folder with filtered FASTQs, updated `units.tsv`, `samples.tsv`,
   and `config.nagini.yaml`.
5. Outputs the exact Snakemake command to run the whole pipeline on the amplicon-filtered data.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


def find_nanostream(custom_path: str | None = None) -> str:
    if custom_path and os.path.isfile(custom_path) and os.access(custom_path, os.X_OK):
        return custom_path

    env_path = os.environ.get("NANOSTREAM_BIN")
    if env_path and os.path.isfile(env_path) and os.access(env_path, os.X_OK):
        return env_path

    which_path = shutil.which("nanostream")
    if which_path:
        return which_path

    candidates = [
        os.path.expanduser("~/bin/nanostream"),
        "/home/cgu-lab/bin/nanostream",
        "/usr/local/bin/nanostream",
        os.path.expanduser("~/.cargo/bin/nanostream"),
        os.path.expanduser("~/tmp/tools/nanoMonitor/target/release/nanostream"),
        "/home/cgu-lab/tmp/tools/nanoMonitor/target/release/nanostream",
    ]

    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c

    raise FileNotFoundError(
        "Could not find 'nanostream' executable. Specify path with --nanostream /path/to/nanostream"
    )


def find_default_primers() -> str:
    script_dir = Path(__file__).resolve().parent
    repo_root = script_dir.parent
    primers_path = repo_root / "workflow" / "resources" / "bcr_abl" / "bcr_abl_primers.tsv"
    if primers_path.is_file():
        return str(primers_path)
    raise FileNotFoundError(f"Default primers file not found at {primers_path}")


def parse_units_tsv(units_file: str) -> list[dict[str, str]]:
    with open(units_file, "r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        return list(reader)


def write_units_tsv(units_file: str, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    with open(units_file, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def run_nanostream_amplicons(
    nanostream_bin: str,
    primers_tsv: str,
    input_fastq: str,
    output_json: str,
    output_fastq_prefix: str,
    threads: int = 8,
    split_by_amplicon: bool = True,
    max_edit_dist: int = 3,
    end_length: int = 150,
) -> dict:
    os.makedirs(os.path.dirname(output_json), exist_ok=True)
    os.makedirs(os.path.dirname(output_fastq_prefix), exist_ok=True)

    cmd = [
        nanostream_bin,
        "amplicons",
        "--primers", primers_tsv,
        "--output", output_json,
        "--output-fastq", output_fastq_prefix,
        "--threads", str(threads),
        "--max-edit-dist", str(max_edit_dist),
        "--end-length", str(end_length),
        "--summary",
    ]
    if split_by_amplicon:
        cmd.append("--split-by-amplicon")
    cmd.append(input_fastq)

    sys.stderr.write(f"\n[Nanostream] Running amplicon filter for {Path(input_fastq).name}...\n")
    sys.stderr.write(f"Command: {' '.join(cmd)}\n")
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.stderr.write(f"Error running nanostream: {proc.stderr}\n")
        raise RuntimeError(f"nanostream amplicons failed with exit code {proc.returncode}")

    if os.path.exists(output_json):
        with open(output_json, "r", encoding="utf-8") as jf:
            return json.load(jf)
    return {}


def main():
    parser = argparse.ArgumentParser(
        description="Filter BCR-ABL1 amplicons using Nanostream and prepare a fresh Nagini analysis folder."
    )
    parser.add_argument(
        "--run-dir", "-i",
        default=".",
        help="Path to existing Nagini analysis run directory (default: current directory)",
    )
    parser.add_argument(
        "--out-dir", "-o",
        required=True,
        help="Path to new analysis folder for the amplicon-filtered experiment",
    )
    parser.add_argument(
        "--primers", "-p",
        default=None,
        help="Path to primers TSV file (default: workflow/resources/bcr_abl/bcr_abl_primers.tsv)",
    )
    parser.add_argument(
        "--nanostream",
        default=None,
        help="Path to nanostream executable",
    )
    parser.add_argument(
        "--threads", "-t",
        type=int,
        default=8,
        help="Threads per sample for nanostream (default: 8)",
    )
    parser.add_argument(
        "--keep-all-amplicons",
        action="store_true",
        help="If set, include both major and minor amplicons in a combined FASTQ for each sample instead of picking the primary",
    )
    args = parser.parse_args()

    run_dir = Path(args.run_dir).resolve()
    out_dir = Path(args.out_dir).resolve()
    primers_tsv = args.primers or find_default_primers()
    nanostream_bin = find_nanostream(args.nanostream)

    units_path = run_dir / "units.tsv"
    if not units_path.is_file():
        sys.exit(f"Error: {units_path} not found. Specify valid --run-dir.")

    samples_path = run_dir / "samples.tsv"
    config_path = run_dir / "config.nagini.yaml"

    units = parse_units_tsv(str(units_path))
    if not units:
        sys.exit(f"Error: No samples found in {units_path}")

    fieldnames = list(units[0].keys())
    fq_col = "fastq1" if "fastq1" in fieldnames else ("fq1" if "fq1" in fieldnames else None)
    if not fq_col:
        sys.exit(f"Error: Neither 'fastq1' nor 'fq1' column found in {units_path}")

    out_dir.mkdir(parents=True, exist_ok=True)
    fastq_amplicon_dir = out_dir / "fastq_amplicon"
    fastq_amplicon_dir.mkdir(parents=True, exist_ok=True)
    stats_dir = out_dir / "amplicon_stats"
    stats_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("NAGINI & NANOSTREAM BCR-ABL1 AMPLICON FILTERING EXPERIMENT")
    print(f"Source run dir:  {run_dir}")
    print(f"Output run dir:  {out_dir}")
    print(f"Nanostream bin:  {nanostream_bin}")
    print(f"Primers TSV:     {primers_tsv}")
    print(f"Samples to process: {len(units)}")
    print("=" * 80)

    new_units = []
    summary_report_table = []

    for row in units:
        sample = row["sample"]
        stype = row.get("type", "N")
        orig_fq = row[fq_col]
        if not os.path.isabs(orig_fq):
            orig_fq = str(run_dir / orig_fq)

        if not os.path.exists(orig_fq):
            sys.stderr.write(f"Warning: FASTQ not found for {sample}: {orig_fq}, keeping original\n")
            new_units.append(row)
            continue

        sample_prefix = f"{sample}_{stype}"
        stats_json_path = str(stats_dir / f"{sample_prefix}.amplicon_stats.json")
        out_fq_prefix = str(fastq_amplicon_dir / f"{sample_prefix}.matched.fastq.gz")

        stats_data = run_nanostream_amplicons(
            nanostream_bin=nanostream_bin,
            primers_tsv=primers_tsv,
            input_fastq=orig_fq,
            output_json=stats_json_path,
            output_fastq_prefix=out_fq_prefix,
            threads=args.threads,
            split_by_amplicon=True,
        )

        amplicons_dict = stats_data.get("amplicons", {})
        total_raw = stats_data.get("total_raw_reads", 0)
        primer_dimers = stats_data.get("failed_pair_match", 0)
        unmatched = stats_data.get("unmatched_count", 0)
        total_matched = sum(info.get("count", 0) for info in amplicons_dict.values())

        # Find the split fastq files created by nanostream
        # They will be named: {out_fq_prefix}_<amplicon_name>.fastq.gz
        prefix_stem = out_fq_prefix[:-9] if out_fq_prefix.endswith(".fastq.gz") else out_fq_prefix
        matching_files = list(fastq_amplicon_dir.glob(f"{Path(prefix_stem).name}_*.fastq.gz"))

        final_fq_path = fastq_amplicon_dir / f"{sample_prefix}.amplicon.fastq.gz"

        if matching_files:
            if len(matching_files) == 1 or not args.keep_all_amplicons:
                # Pick the file with the largest size or major amplicon
                matching_files.sort(key=lambda p: p.stat().st_size, reverse=True)
                primary_split = matching_files[0]
                if final_fq_path.exists():
                    final_fq_path.unlink()
                # Create a relative symlink
                os.symlink(primary_split.name, str(final_fq_path))
            else:
                # Concatenate all split amplicons into one final fastq
                with open(final_fq_path, "wb") as outfile:
                    for mf in matching_files:
                        with open(mf, "rb") as infile:
                            shutil.copyfileobj(infile, outfile)
        else:
            sys.stderr.write(f"Warning: No amplicon FASTQs generated for {sample_prefix}\n")
            final_fq_path = Path(orig_fq)

        new_row = dict(row)
        # Explicitly set input_type to fastq and clear old BAM path so Nagini realigns from filtered FASTQ
        new_row["input_type"] = "fastq"
        if "bam" in new_row:
            new_row["bam"] = ""

        # Store path relative to out_dir if inside, else absolute
        try:
            rel_path = final_fq_path.relative_to(out_dir)
            new_row[fq_col] = str(rel_path)
        except ValueError:
            new_row[fq_col] = str(final_fq_path)

        new_units.append(new_row)


        summary_report_table.append({
            "sample": sample,
            "type": stype,
            "raw_reads": total_raw,
            "matched_amplicons": total_matched,
            "match_rate": f"{(total_matched / max(total_raw, 1)) * 100.0:.1f}%",
            "primer_dimers": primer_dimers,
            "unmatched": unmatched,
            "amplicons": ", ".join(f"{k} ({v['count']})" for k, v in amplicons_dict.items()) or "None",
        })

    # Write new units.tsv
    new_units_path = out_dir / "units.tsv"
    write_units_tsv(str(new_units_path), new_units, fieldnames)

    # Copy samples.tsv and config.nagini.yaml
    if samples_path.is_file():
        shutil.copy2(str(samples_path), str(out_dir / "samples.tsv"))
    if config_path.is_file():
        shutil.copy2(str(config_path), str(out_dir / "config.nagini.yaml"))

    # Print summary table
    print("\n" + "=" * 80)
    print("AMPLICON FILTERING RESULTS SUMMARY")
    print("=" * 80)
    print(f"{'Sample':<20} | {'Raw Reads':<10} | {'Matched':<10} | {'Match %':<8} | {'Dimers':<8} | {'Amplicons (count)'}")
    print("-" * 80)
    for r in summary_report_table:
        print(f"{r['sample']:<20} | {r['raw_reads']:<10} | {r['matched_amplicons']:<10} | {r['match_rate']:<8} | {r['primer_dimers']:<8} | {r['amplicons']}")
    print("=" * 80)

    print("\nSUCCESS! New amplicon-filtered analysis directory is ready at:")
    print(f"  {out_dir}\n")
    print("To run Nagini on the amplicon-filtered dataset, run:")
    print(f"  cd {out_dir}")
    print("  snakemake --profile /home/cgu-lab/tmp/tools/nagini/profiles/workstation \\")
    print("    --snakefile /home/cgu-lab/tmp/tools/nagini/workflow/Snakefile \\")
    print("    --configfile config.nagini.yaml\n")


if __name__ == "__main__":
    main()
