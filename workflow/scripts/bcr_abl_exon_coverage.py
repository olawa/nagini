#!/usr/bin/env python3
"""Calculate per-exon coverage statistics for BCR-ABL1 amplicon sequencing.

Computes mean, median, min, max depth, read counts, and percentage of bases
meeting diagnostic thresholds (>=100x, >=500x, >=1000x) per BCR and ABL1 exon.
"""

from __future__ import annotations

import argparse
import statistics
import sys
from pathlib import Path
from typing import Any, Dict, List

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path, help="Input aligned BAM file")
    parser.add_argument("--bed", required=True, type=Path, help="Exon BED file (bcr_abl1_exons_*.bed)")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--output", required=True, type=Path, help="Output exon coverage TSV path")
    parser.add_argument("--min-mapq", type=int, default=1, help="Minimum mapping quality (default 1)")
    parser.add_argument("--max-depth", type=int, default=1_000_000, help="Maximum pileup depth")
    return parser.parse_args()


def read_exon_bed(bed_path: Path) -> List[Dict[str, Any]]:
    """Parse exon BED file (chrom, start_0based, end_1based, exon_name)."""
    exons = []
    with bed_path.open() as handle:
        for line in handle:
            text = line.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t")
            if len(fields) < 4:
                fields = text.split()
            if len(fields) < 4:
                continue
            chrom, start, end, name = fields[:4]
            exons.append({
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "name": name,
            })
    return exons


def calculate_exon_coverage(
    bam_path: Path,
    exons: List[Dict[str, Any]],
    sample: str,
    min_mapq: int = 1,
    max_depth: int = 1_000_000,
) -> List[Dict[str, Any]]:
    """Calculate depth metrics for each exon in the BED file."""
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    bam_contigs = set(bam.references)

    results = []
    try:
        for exon in exons:
            chrom = exon["chrom"]
            if chrom not in bam_contigs:
                alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
                if alt_chrom in bam_contigs:
                    chrom = alt_chrom

            start_0based = exon["start"]
            end_1based = exon["end"]
            length = end_1based - start_0based
            exon_name = exon["name"]

            base_depths: Dict[int, int] = {pos: 0 for pos in range(start_0based, end_1based)}
            covering_reads = set()

            if chrom in bam_contigs and length > 0:
                for column in bam.pileup(
                    chrom,
                    start_0based,
                    end_1based,
                    truncate=True,
                    stepper="nofilter",
                    max_depth=max_depth,
                    min_base_quality=0,
                    ignore_overlaps=False,
                ):
                    pos = column.reference_pos
                    if pos not in base_depths:
                        continue

                    valid_depth = 0
                    for read in column.pileups:
                        aln = read.alignment
                        if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail:
                            continue
                        if (aln.mapping_quality or 0) < min_mapq:
                            continue
                        if read.is_refskip:
                            continue
                        if read.is_del or read.query_position is None:
                            continue

                        valid_depth += 1
                        covering_reads.add(aln.query_name)

                    base_depths[pos] = valid_depth

            depth_values = list(base_depths.values())
            if not depth_values:
                depth_values = [0]

            mean_d = sum(depth_values) / len(depth_values)
            median_d = statistics.median(depth_values)
            min_d = min(depth_values)
            max_d = max(depth_values)

            gte_100 = sum(1 for d in depth_values if d >= 100) / len(depth_values) * 100.0
            gte_500 = sum(1 for d in depth_values if d >= 500) / len(depth_values) * 100.0
            gte_1000 = sum(1 for d in depth_values if d >= 1000) / len(depth_values) * 100.0

            results.append({
                "sample": sample,
                "exon": exon_name,
                "chrom": exon["chrom"],
                "start": start_0based + 1,  # 1-based start
                "end": end_1based,
                "length_bp": length,
                "mean_depth": round(mean_d, 1),
                "median_depth": round(float(median_d), 1),
                "min_depth": min_d,
                "max_depth": max_d,
                "pct_gte_100x": round(gte_100, 1),
                "pct_gte_500x": round(gte_500, 1),
                "pct_gte_1000x": round(gte_1000, 1),
                "spanning_reads": len(covering_reads),
            })
    finally:
        bam.close()

    return results


def write_exon_coverage_tsv(results: List[Dict[str, Any]], output_path: Path) -> None:
    """Write exon coverage results to TSV."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    columns = [
        "sample", "exon", "chrom", "start", "end", "length_bp",
        "mean_depth", "median_depth", "min_depth", "max_depth",
        "pct_gte_100x", "pct_gte_500x", "pct_gte_1000x", "spanning_reads",
    ]
    with output_path.open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in results:
            line = "\t".join(str(row.get(col, "")) for col in columns)
            handle.write(line + "\n")


def main() -> None:
    args = parse_args()
    if not args.bam.is_file():
        sys.exit(f"BAM file not found: {args.bam}")
    if not args.bed.is_file():
        sys.exit(f"BED file not found: {args.bed}")

    exons = read_exon_bed(args.bed)
    if not exons:
        sys.exit(f"No exons parsed from {args.bed}")

    results = calculate_exon_coverage(
        bam_path=args.bam,
        exons=exons,
        sample=args.sample,
        min_mapq=args.min_mapq,
        max_depth=args.max_depth,
    )
    write_exon_coverage_tsv(results, args.output)
    print(f"Wrote coverage statistics for {len(results)} exons to {args.output}")


if __name__ == "__main__":
    main()
