#!/usr/bin/env python3
"""Collect what ``dorado demux`` wrote into one file per barcode and one summary.

``summary``  copy Dorado's demultiplexing summary out of the raw output directory
``barcode``  make one unaligned BAM for a barcode from the BAMs Dorado wrote for it

Dorado normally writes one BAM per barcode, named ``<kit>_<barcode>.bam`` (older
versions used a directory per barcode).  Both layouts are understood.  A barcode
matches a file only as a whole token, so ``barcode1`` never picks up the reads of
``barcode10``.  These steps were shell blocks in dorado.smk.

Standard library only, plus ``samtools`` for the two cases that need it.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Callable, List, Optional, Sequence


class DemuxError(RuntimeError):
    pass


def barcode_pattern(barcode: str) -> "re.Pattern[str]":
    """The barcode as a whole token: not glued to letters or digits on either side."""

    return re.compile(rf"(?<![0-9A-Za-z]){re.escape(barcode)}(?![0-9A-Za-z])")


def find_barcode_bams(raw_dir: str, barcode: str) -> List[Path]:
    """BAMs under raw_dir that belong to the barcode, in a stable order."""

    root = Path(raw_dir)
    token = barcode_pattern(barcode)
    found = []
    for path in sorted(root.rglob("*.bam")):
        if not path.is_file():
            continue
        relative_dirs = path.relative_to(root).parts[:-1]
        if barcode in relative_dirs or path.stem == barcode or token.search(path.name):
            found.append(path)
    return found


def find_summary(raw_dir: str) -> Path:
    """Dorado's summary file.  If there are several, a barcoding summary wins."""

    candidates = sorted(
        path for path in Path(raw_dir).rglob("*")
        if path.is_file()
        and "summary" in path.name.lower()
        and path.suffix.lower() in (".txt", ".tsv")
    )
    if not candidates:
        raise DemuxError("Dorado demux did not emit the requested summary file")
    for path in candidates:
        if "barcod" in path.name.lower():
            return path
    return candidates[0]


def collect_summary(raw_dir: str, output: str) -> Path:
    source = find_summary(raw_dir)
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(source, output)
    print(f"Copied {source} to {output}", flush=True)
    return source


def header_only_bam(basecalled: str, output: str, run: Callable = subprocess.run) -> None:
    """A BAM with the header of the basecalled file and no reads."""

    header = run(["samtools", "view", "-H", basecalled], stdout=subprocess.PIPE)
    if header.returncode != 0:
        raise DemuxError(f"samtools view -H failed on {basecalled}")
    written = run(["samtools", "view", "-b", "-o", output, "-"], input=header.stdout)
    if written.returncode != 0:
        raise DemuxError(f"samtools could not write the header-only BAM {output}")


def collect_barcode(
    raw_dir: str,
    barcode: str,
    basecalled: str,
    output: str,
    allow_empty: bool,
    run: Callable = subprocess.run,
) -> str:
    """Write the barcode's BAM.  Returns "empty", "linked" or "concatenated"."""

    bams = find_barcode_bams(raw_dir, barcode)
    Path(output).parent.mkdir(parents=True, exist_ok=True)

    if not bams:
        if not allow_empty:
            raise DemuxError(f"No demultiplexed reads found for {barcode}")
        print(f"WARNING: no reads for {barcode}; writing a header-only BAM", flush=True)
        header_only_bam(basecalled, output, run)
        return "empty"

    if len(bams) == 1:
        # One BAM per barcode is the normal case.  A hard link avoids rewriting
        # the whole file while the raw demux directory is still there.
        if os.path.lexists(output):
            os.remove(output)
        try:
            os.link(bams[0], output)
            print(f"Linked {bams[0]} to {output}", flush=True)
        except OSError:
            shutil.copyfile(bams[0], output)
            print(f"Copied {bams[0]} to {output} (a hard link was not possible)", flush=True)
        return "linked"

    print(f"Concatenating {len(bams)} BAMs for {barcode}", flush=True)
    result = run(["samtools", "cat", "-o", output] + [str(path) for path in bams])
    if result.returncode != 0:
        raise DemuxError(f"samtools cat failed for {barcode}")
    return "concatenated"


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="mode", required=True)

    summary = sub.add_parser("summary")
    summary.add_argument("--raw-dir", required=True)
    summary.add_argument("--output", required=True)

    barcode = sub.add_parser("barcode")
    barcode.add_argument("--raw-dir", required=True)
    barcode.add_argument("--barcode", required=True)
    barcode.add_argument("--basecalled", required=True)
    barcode.add_argument("--output", required=True)
    barcode.add_argument("--allow-empty", action="store_true")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        if args.mode == "summary":
            collect_summary(args.raw_dir, args.output)
        else:
            collect_barcode(args.raw_dir, args.barcode, args.basecalled, args.output, args.allow_empty)
    except DemuxError as error:
        print(f"Error: {error}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
