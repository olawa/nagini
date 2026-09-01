#!/usr/bin/env python3
"""Fallback read filter implementing nanofilter CLI syntax for FASTQ and BAM."""

import argparse
import gzip
import os
import sys
from pathlib import Path
from typing import Optional, Tuple


def parse_len_range(len_arg: Optional[str]) -> Tuple[Optional[int], Optional[int]]:
    if not len_arg:
        return None, None
    raw = str(len_arg).strip().replace("..", "-").replace(":", "-")
    if "-" in raw:
        parts = raw.split("-", 1)
        min_l = int(parts[0]) if parts[0].strip() else None
        max_l = int(parts[1]) if parts[1].strip() else None
        return min_l, max_l
    try:
        val = int(raw)
        return val, None
    except ValueError:
        return None, None


def mean_q_from_qual_str(qual_str: str) -> float:
    if not qual_str:
        return 0.0
    return sum(ord(c) - 33 for c in qual_str) / len(qual_str)


def filter_fastq(input_path: Path, output_path: Path, min_qv: float, min_len: Optional[int], max_len: Optional[int]) -> None:
    is_in_gz = str(input_path).endswith(".gz")
    is_out_gz = str(output_path).endswith(".gz")

    open_in = gzip.open(input_path, "rt") if is_in_gz else open(input_path, "r")
    open_out = gzip.open(output_path, "wt") if is_out_gz else open(output_path, "w")

    total_reads = 0
    passed_reads = 0

    with open_in as handle_in, open_out as handle_out:
        while True:
            header = handle_in.readline()
            if not header:
                break
            seq = handle_in.readline()
            plus = handle_in.readline()
            qual = handle_in.readline()

            if not qual:
                break

            total_reads += 1
            seq_str = seq.strip()
            qual_str = qual.strip()
            seq_len = len(seq_str)

            if min_len is not None and seq_len < min_len:
                continue
            if max_len is not None and seq_len > max_len:
                continue
            if min_qv > 0 and mean_q_from_qual_str(qual_str) < min_qv:
                continue

            passed_reads += 1
            handle_out.write(header)
            handle_out.write(seq)
            handle_out.write(plus)
            handle_out.write(qual)

    pct = (passed_reads / total_reads * 100.0) if total_reads > 0 else 0.0
    print(f"nanofilter [fallback]: filtered {input_path} -> {output_path}", file=sys.stderr)
    print(f"Total reads: {total_reads}, Passed: {passed_reads} ({pct:.1f}%), Filtered: {total_reads - passed_reads}", file=sys.stderr)


def filter_bam(input_path: Path, output_path: Path, min_qv: float, min_len: Optional[int], max_len: Optional[int], threads: int = 1) -> None:
    import pysam

    in_bam = pysam.AlignmentFile(str(input_path), "rb", threads=threads)
    out_bam = pysam.AlignmentFile(str(output_path), "wb", template=in_bam, threads=threads)

    total_reads = 0
    passed_reads = 0

    try:
        for read in in_bam.fetch(until_eof=True):
            total_reads += 1
            qlen = read.query_length
            if min_len is not None and qlen < min_len:
                continue
            if max_len is not None and qlen > max_len:
                continue
            if min_qv > 0:
                quals = read.query_qualities
                if quals is not None and len(quals) > 0:
                    mean_q = sum(quals) / len(quals)
                    if mean_q < min_qv:
                        continue
            passed_reads += 1
            out_bam.write(read)
    finally:
        in_bam.close()
        out_bam.close()

    pct = (passed_reads / total_reads * 100.0) if total_reads > 0 else 0.0
    print(f"nanofilter [fallback]: filtered {input_path} -> {output_path}", file=sys.stderr)
    print(f"Total reads: {total_reads}, Passed: {passed_reads} ({pct:.1f}%), Filtered: {total_reads - passed_reads}", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Fallback read filter")
    subparsers = parser.add_subparsers(dest="subcommand")

    filter_parser = subparsers.add_parser("filter", help="Filter reads")
    filter_parser.add_argument("-i", "--input", dest="input_file", type=Path, default=None, help="Input BAM/FASTQ(.gz) file")
    filter_parser.add_argument("positional_input", nargs="?", type=Path, default=None, help="Positional input BAM/FASTQ(.gz) file")
    filter_parser.add_argument("-o", "--output", type=Path, required=True, help="Output file")
    filter_parser.add_argument("--qv", type=float, default=20.0, help="Minimum mean QV [default: 20]")
    filter_parser.add_argument("--len", type=str, default="", help="Length range min-max")
    filter_parser.add_argument("-t", "--threads", type=int, default=1, help="Threads")

    args, unknown = parser.parse_known_args()
    if args.subcommand != "filter":
        # Check if called as: script -i input -o output --qv ...
        parser_direct = argparse.ArgumentParser()
        parser_direct.add_argument("-i", "--input", dest="input_file", type=Path, default=None)
        parser_direct.add_argument("positional_input", nargs="?", type=Path, default=None)
        parser_direct.add_argument("-o", "--output", type=Path, required=True)
        parser_direct.add_argument("--qv", type=float, default=20.0)
        parser_direct.add_argument("--len", type=str, default="")
        parser_direct.add_argument("-t", "--threads", type=int, default=1)
        args, _ = parser_direct.parse_known_args()

    input_target = args.input_file or getattr(args, "positional_input", None)
    if not input_target:
        sys.exit("nanofilter error: input file is required via -i/--input or positional argument.")

    min_l, max_l = parse_len_range(args.len)
    in_str = str(input_target).lower()
    if in_str.endswith(".bam"):
        filter_bam(input_target, args.output, args.qv, min_l, max_l, threads=args.threads)
    else:
        filter_fastq(input_target, args.output, args.qv, min_l, max_l)



if __name__ == "__main__":
    main()
