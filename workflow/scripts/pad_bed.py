#!/usr/bin/env python3
"""Pad BED intervals by a fixed number of base pairs and merge overlapping intervals."""

import argparse
import collections
import os
import sys


def pad_and_merge_intervals(intervals, padding_bp=25000):
    """
    Pad each interval by +/- padding_bp and merge overlapping/adjacent intervals.
    intervals: list of (start, end, name) tuples.
    Returns: list of (start, end, names_list) tuples.
    """
    padded = []
    for item in intervals:
        start = max(0, int(item[0]) - padding_bp)
        end = int(item[1]) + padding_bp
        name = str(item[2]) if len(item) > 2 else ""
        padded.append((start, end, name))

    padded.sort(key=lambda x: (x[0], x[1]))
    merged = []
    for start, end, name in padded:
        if not merged:
            merged.append([start, end, [name] if name else []])
        else:
            prev_start, prev_end, prev_names = merged[-1]
            if start <= prev_end:  # Overlapping or adjacent
                merged[-1][1] = max(prev_end, end)
                if name and name not in prev_names:
                    prev_names.append(name)
            else:
                merged.append([start, end, [name] if name else []])

    return [(start, end, names) for start, end, names in merged]


def pad_and_merge_bed(input_bed_path, output_bed_path, padding_bp=25000):
    """Read a BED file, apply padding, merge overlapping intervals per chromosome, and write output BED."""
    chrom_intervals = collections.defaultdict(list)
    if not os.path.exists(input_bed_path):
        raise FileNotFoundError(f"Input BED file not found: {input_bed_path}")

    with open(input_bed_path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            name = parts[3] if len(parts) > 3 else ""
            chrom_intervals[chrom].append((start, end, name))

    os.makedirs(os.path.dirname(os.path.abspath(output_bed_path)), exist_ok=True)
    with open(output_bed_path, "w") as out:
        for chrom in sorted(chrom_intervals.keys()):
            merged = pad_and_merge_intervals(chrom_intervals[chrom], padding_bp=padding_bp)
            for start, end, names in merged:
                name_str = ";".join(names) if names else f"{chrom}:{start}-{end}"
                out.write(f"{chrom}\t{start}\t{end}\t{name_str}\n")


def main():
    parser = argparse.ArgumentParser(description="Pad and merge BED intervals")
    parser.add_argument("--input-bed", required=True, help="Input BED file path")
    parser.add_argument("--output-bed", required=True, help="Output padded BED file path")
    parser.add_argument("--padding-bp", type=int, default=25000, help="Padding in base pairs (default: 25000)")
    args = parser.parse_args()
    pad_and_merge_bed(args.input_bed, args.output_bed, padding_bp=args.padding_bp)


if "snakemake" in globals():
    pad_and_merge_bed(
        input_bed_path=str(snakemake.input.bed),
        output_bed_path=str(snakemake.output.bed),
        padding_bp=int(getattr(snakemake.params, "padding_bp", 25000) or 25000),
    )
elif __name__ == "__main__":
    main()
