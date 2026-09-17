#!/usr/bin/env python3
"""Derive a centromere BED from a cytoband table.

Rindels' ideogram falls back to a built-in hg38 centromere set when it is not
given one. On T2T that draws the centromeres in the wrong places, which for a
dic(9;20) ideogram is precisely the wrong kind of wrong, so the positions are
taken from whichever cytoband table the reference build actually registered.

The ``acen`` rows of a UCSC-style cytoband file are the centromere, one per
arm. They are merged per chromosome so each chromosome yields a single
interval, which is what an ideogram wants to draw.
"""

from __future__ import annotations

import argparse
import gzip
import sys
from pathlib import Path


CENTROMERE_STAINS = {"acen"}


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def centromeres(handle) -> list[tuple[str, int, int]]:
    """One merged interval per chromosome, in first-seen order."""

    spans: dict[str, list[int]] = {}
    order: list[str] = []
    for line in handle:
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 5:
            continue
        stain = fields[4].strip().lower()
        if stain not in CENTROMERE_STAINS:
            continue
        try:
            start, end = int(fields[1]), int(fields[2])
        except ValueError:
            continue
        if end <= start:
            continue
        chrom = fields[0].strip()
        if chrom not in spans:
            spans[chrom] = [start, end]
            order.append(chrom)
        else:
            spans[chrom][0] = min(spans[chrom][0], start)
            spans[chrom][1] = max(spans[chrom][1], end)
    return [(chrom, spans[chrom][0], spans[chrom][1]) for chrom in order]


def write_bed(rows, destination: Path) -> None:
    with destination.open("w") as handle:
        for chrom, start, end in rows:
            handle.write(f"{chrom}\t{start}\t{end}\tcentromere\n")


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cytoband", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args(argv)

    destination = Path(args.output)
    destination.parent.mkdir(parents=True, exist_ok=True)

    source = Path(args.cytoband)
    if not source.is_file():
        # An empty BED is a valid answer: the ideogram then draws none rather
        # than the ones from another reference build.
        write_bed([], destination)
        print(f"No cytoband table at {source}; wrote an empty centromere BED.", file=sys.stderr)
        return 0

    with open_text(source) as handle:
        rows = centromeres(handle)
    write_bed(rows, destination)
    if not rows:
        print(
            f"{source} has no acen rows; wrote an empty centromere BED.",
            file=sys.stderr,
        )
    else:
        print(f"Wrote {len(rows)} centromere interval(s) to {destination}")
    return 0


if __name__ == "__main__":
    if "snakemake" in globals():  # noqa: F821
        main(
            [
                "--cytoband",
                str(globals()["snakemake"].params.cytoband),  # noqa: F821
                "--output",
                str(globals()["snakemake"].output.bed),  # noqa: F821
            ]
        )
    else:
        raise SystemExit(main())
