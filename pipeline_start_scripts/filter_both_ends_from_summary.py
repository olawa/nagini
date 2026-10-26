#!/usr/bin/env python3
"""Filter existing Nanopore BAM files for both-ends barcoding, length, and QV
using ONT sequencing_summary.txt.

MinKNOW/Guppy sequencing_summary.txt already records:
- barcode_arrangement (barcode01, barcode02, etc.)
- barcode_front_score and barcode_rear_score
- mean_qscore_template (QV)
- sequence_length_template (length in bp)

Reads with barcode_front_score >= 60 and barcode_rear_score >= 60 have verified
barcodes at BOTH ends (two-sided barcoding), avoiding chimeric / misassigned reads.

This script extracts valid read IDs meeting these criteria and filters existing
aligned or unmapped BAM files, producing a clean directory structure:
    output_dir/barcodeXX/<sample>.bam
ready for Nagini:
    start_nagini.sh --input-type aligned_bam --inbox-path output_dir ...
"""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


def normalize_barcode(val: str) -> str:
    """Normalize barcode token to lowercase barcodeXX."""
    if not val:
        return ""
    val_clean = str(val).strip().lower()
    m = re.search(r"(?:barcode|bc)(\d+)", val_clean)
    if m:
        num = int(m.group(1))
        return f"barcode{num:02d}"
    return val_clean


def parse_summary_file(
    summary_path: Path,
    min_qv: float = 14.0,
    min_len: int = 800,
    max_len: int = 2500,
    min_barcode_score: float = 60.0,
    target_barcodes: Optional[Set[str]] = None,
) -> Tuple[Dict[str, Set[str]], Dict[str, Dict[str, int]]]:
    """Parse sequencing_summary.txt and return matching read IDs per barcode.
    
    Returns:
        (reads_by_barcode, stats_by_barcode)
    """
    reads_by_barcode: Dict[str, Set[str]] = defaultdict(set)
    stats: Dict[str, Dict[str, int]] = defaultdict(lambda: {
        "total": 0,
        "pass_qv": 0,
        "pass_len": 0,
        "pass_both_ends": 0,
        "pass_all": 0,
    })

    open_fn = gzip.open if summary_path.suffix == ".gz" else open
    mode = "rt" if summary_path.suffix == ".gz" else "r"

    sys.stderr.write(f"[Summary] Scanning {summary_path.name}...\n")
    with open_fn(summary_path, mode, encoding="utf-8", errors="replace") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            arr = normalize_barcode(row.get("barcode_arrangement", ""))
            if not arr or arr == "unclassified":
                continue
            if target_barcodes and arr not in target_barcodes:
                continue

            read_id = row.get("read_id", "").strip()
            if not read_id:
                continue

            stats[arr]["total"] += 1

            try:
                q = float(row.get("mean_qscore_template", 0.0))
            except (ValueError, TypeError):
                q = 0.0

            try:
                length = int(row.get("sequence_length_template", 0))
            except (ValueError, TypeError):
                length = 0

            try:
                f_score = float(row.get("barcode_front_score", 0.0))
                r_score = float(row.get("barcode_rear_score", 0.0))
            except (ValueError, TypeError):
                f_score, r_score = 0.0, 0.0

            pass_q = (q >= min_qv)
            pass_l = (min_len <= length <= max_len)
            pass_ends = (f_score >= min_barcode_score and r_score >= min_barcode_score)

            if pass_q:
                stats[arr]["pass_qv"] += 1
            if pass_l:
                stats[arr]["pass_len"] += 1
            if pass_ends:
                stats[arr]["pass_both_ends"] += 1

            if pass_q and pass_l and pass_ends:
                stats[arr]["pass_all"] += 1
                reads_by_barcode[arr].add(read_id)

    return dict(reads_by_barcode), dict(stats)


def find_bams_for_barcodes(
    input_dir: Path, target_barcodes: Set[str]
) -> Dict[str, Path]:
    """Find the best matching BAM file for each barcode in input_dir."""
    barcode_to_bam: Dict[str, Path] = {}
    all_bams = sorted(input_dir.rglob("*.bam"))

    for bc in target_barcodes:
        bc_num = int(re.search(r"\d+", bc).group(0))
        candidates = []
        for p in all_bams:
            # Match folder name (e.g. /barcode01/) or filename token
            parts = [part.lower() for part in p.parts]
            stem = p.stem.lower()
            if bc in parts or f"bc{bc_num}" in parts or f"bc{bc_num:02d}" in parts:
                candidates.append(p)
            elif re.search(rf"(?:barcode|bc){bc_num:02d}(?:[^0-9]|$)", stem) or re.search(rf"(?:barcode|bc){bc_num}(?:[^0-9]|$)", stem):
                candidates.append(p)

        if candidates:
            # Prefer files with 'merged' or 'sorted' or shortest path
            candidates.sort(key=lambda x: (
                0 if "merge" in x.name.lower() or "sort" in x.name.lower() else 1,
                len(str(x))
            ))
            barcode_to_bam[bc] = candidates[0]

    return barcode_to_bam


def filter_bam_samtools(
    input_bam: Path,
    output_bam: Path,
    read_ids: Set[str],
    threads: int = 4
) -> None:
    """Filter input_bam to output_bam retaining only read_ids using samtools."""
    output_bam.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile("w", delete=False) as qfile:
        qfile_path = Path(qfile.name)
        for rid in sorted(read_ids):
            qfile.write(f"{rid}\n")

    try:
        cmd = [
            "samtools", "view",
            "-N", str(qfile_path),
            "-b",
            "-@", str(threads),
            "-o", str(output_bam),
            str(input_bam)
        ]
        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f"samtools view failed: {res.stderr}")

        # Index the output BAM
        index_cmd = ["samtools", "index", "-@", str(threads), str(output_bam)]
        subprocess.run(index_cmd, capture_output=True, check=True)
    finally:
        if qfile_path.is_file():
            qfile_path.unlink()


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Filter existing Nanopore BAM files by both-ends barcoding, QV, and length from sequencing_summary.txt."
    )
    parser.add_argument(
        "--summary", "-s", type=Path, required=True,
        help="Path to ONT sequencing_summary.txt or .txt.gz"
    )
    parser.add_argument(
        "--input-bams", "-i", type=Path, required=True,
        help="Directory containing existing BAM files to filter (e.g. from previous run)"
    )
    parser.add_argument(
        "--output-dir", "-o", type=Path, required=True,
        help="Directory to write filtered BAMs (structured as output_dir/barcodeXX/<sample>.bam)"
    )
    parser.add_argument(
        "--barcode-map", "-m", type=Path, default=None,
        help="Optional barcode map CSV/TSV to assign sample names to BAM files"
    )
    parser.add_argument(
        "--min-qv", type=float, default=14.0,
        help="Minimum mean Q-score (default: 14.0)"
    )
    parser.add_argument(
        "--min-len", type=int, default=800,
        help="Minimum read length in bp (default: 800)"
    )
    parser.add_argument(
        "--max-len", type=int, default=2500,
        help="Maximum read length in bp (default: 2500)"
    )
    parser.add_argument(
        "--min-barcode-score", type=float, default=60.0,
        help="Minimum score for BOTH front and rear barcode (default: 60.0)"
    )
    parser.add_argument(
        "--threads", "-t", type=int, default=8,
        help="Threads for samtools (default: 8)"
    )

    args = parser.parse_args(argv)

    if not args.summary.is_file():
        sys.stderr.write(f"Error: summary file not found: {args.summary}\n")
        return 1

    if not args.input_bams.is_dir():
        sys.stderr.write(f"Error: input BAM directory not found: {args.input_bams}\n")
        return 1

    # Load barcode map if provided
    sample_names: Dict[str, str] = {}
    if args.barcode_map and args.barcode_map.is_file():
        delimiter = "," if args.barcode_map.suffix == ".csv" else "\t"
        with open(args.barcode_map, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter=delimiter)
            for r in reader:
                bc = normalize_barcode(r.get("barcode", "") or r.get("barcode_name", ""))
                sname = r.get("sample", "") or r.get("sample_id", "") or r.get("sample_name", "")
                if bc and sname:
                    sample_names[bc] = sname

    # 1. Parse sequencing_summary
    reads_by_bc, stats = parse_summary_file(
        args.summary,
        min_qv=args.min_qv,
        min_len=args.min_len,
        max_len=args.max_len,
        min_barcode_score=args.min_barcode_score,
        target_barcodes=set(sample_names.keys()) if sample_names else None,
    )

    if not reads_by_bc:
        sys.stderr.write("No reads met the filtering criteria in the summary file.\n")
        return 1

    # 2. Find input BAMs
    bams = find_bams_for_barcodes(args.input_bams, set(reads_by_bc.keys()))
    if not bams:
        sys.stderr.write(f"Could not find matching BAM files for barcodes in {args.input_bams}\n")
        return 1

    print("\n" + "=" * 80)
    print("NAGINI SEQUENCING_SUMMARY BOTH-ENDS FILTER")
    print(f"Summary file:       {args.summary}")
    print(f"Input BAMs:         {args.input_bams}")
    print(f"Output directory:   {args.output_dir}")
    print(f"Filter settings:    QV >= {args.min_qv}, Length {args.min_len}–{args.max_len} bp, Both-ends score >= {args.min_barcode_score}")
    print("=" * 80 + "\n")

    print(f"{'Barcode':<12} {'Sample Name':<20} {'Summary Total':<14} {'Both-Ends %':<12} {'Passed Reads':<14} {'Status'}")
    print("-" * 80)

    # 3. Filter BAMs
    for bc in sorted(reads_by_bc.keys()):
        valid_ids = reads_by_bc[bc]
        bc_stats = stats[bc]
        sname = sample_names.get(bc, bc)
        in_bam = bams.get(bc)

        if not in_bam:
            print(f"{bc:<12} {sname:<20} {bc_stats['total']:<14} {bc_stats['pass_both_ends']/max(1, bc_stats['total'])*100:5.1f}%     {len(valid_ids):<14} [MISSING INPUT BAM]")
            continue

        pct = (len(valid_ids) / max(1, bc_stats["total"])) * 100
        print(f"{bc:<12} {sname:<20} {bc_stats['total']:<14} {bc_stats['pass_both_ends']/max(1, bc_stats['total'])*100:5.1f}%     {len(valid_ids):<14} Filtering...")

        out_bam = args.output_dir / bc / f"{sname}.bam"
        filter_bam_samtools(in_bam, out_bam, valid_ids, threads=args.threads)

    print("-" * 80)
    print(f"\n[Done] Filtered BAMs saved to: {args.output_dir}")
    print("You can now run Nagini with:")
    print(f"  pipeline_start_scripts/start_nagini.sh \\")
    print(f"    --site workstation \\")
    print(f"    --inbox-path {args.output_dir} \\")
    print(f"    --input-type aligned_bam \\")
    print(f"    --analysis-dir <YOUR_ANALYSIS_DIR> \\")
    print(f"    --assay bcr_abl --library-type cdna --reference-build hg19 ...\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
