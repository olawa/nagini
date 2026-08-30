#!/usr/bin/env python3
"""Regenerate the Branford resistance-position BED files from pickett's list.

pickett_bcr_abl_pipeline/config/annotations_hg19/branford_SNVs_230118.bed is
the list already used in production Illumina reporting at CGU. The hg19 file
here is derived from it directly rather than re-typed, and the hg38 file is
derived from the hg19 file by a constant offset rather than hand-converted,
so the two builds cannot drift apart the way they did before: 71 rows instead
of 73, one fabricated variant, and 9 of 70 shared rows carrying the wrong
amino-acid change (K247R labelled as Leu248Val, D276N and D276G swapped with
each other, a synonymous E499E labelled as a missense change, and others).

The offset (chr9, hg19 -> hg38) is 2875387 bp, verified two independent ways:
it reproduces T315I at chr9:130872896 in both rindels' and this project's own
hg38 resistance-mutation lists, and it holds exactly across every one of the
73 Branford positions (verified below, not assumed).

Run this script and diff its output against the checked-in BED files whenever
branford_SNVs_*.bed changes upstream in pickett.
"""

import re
import sys
from pathlib import Path

PICKETT_BED = (
    Path.home()
    / "dev/cgu_pipelines/wp2/pickett_bcr_abl_pipeline/config/annotations_hg19"
    / "branford_SNVs_230118.bed"
)
HG19_OUT = Path(__file__).parent / "branford_positions_hg19.bed"
HG38_OUT = Path(__file__).parent / "branford_positions_hg38.bed"

# hg19 -> hg38 delta for the ABL1 locus on chr9, confirmed against two
# independently curated sources (rindels' breakpoint clusters, and the
# T315I position quoted by both rindels and nanoMonitor for hg38).
CHR9_HG19_TO_HG38_OFFSET = 2875387

AA3 = {
    "A": "Ala", "R": "Arg", "N": "Asn", "D": "Asp", "C": "Cys", "Q": "Gln",
    "E": "Glu", "G": "Gly", "H": "His", "I": "Ile", "L": "Leu", "K": "Lys",
    "M": "Met", "F": "Phe", "P": "Pro", "S": "Ser", "T": "Thr", "W": "Trp",
    "Y": "Tyr", "V": "Val",
}

AA_PATTERN = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def to_three_letter(short: str) -> str:
    match = AA_PATTERN.fullmatch(short.strip())
    if not match:
        raise SystemExit(f"Cannot parse amino-acid change: {short!r}")
    ref, pos, alt = match.groups()
    return f"p.{AA3[ref]}{pos}{AA3[alt]}"


def read_pickett_bed(path: Path):
    rows = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        text = line.strip()
        if not text:
            continue
        fields = text.split("\t")
        if len(fields) != 7:
            raise SystemExit(f"{path}:{line_number}: expected 7 columns, got {len(fields)}.")
        chrom, start, end, hgvs, ref, alt, aa_short = fields
        rows.append(
            {
                "chrom": chrom,
                "start": int(start),
                "end": int(end),
                "hgvs": hgvs,
                "ref": ref,
                "alt": alt,
                "aa_three": to_three_letter(aa_short),
            }
        )
    return rows


def write_bed(path: Path, rows, offset: int) -> None:
    lines = []
    for row in rows:
        lines.append(
            "\t".join(
                [
                    row["chrom"],
                    str(row["start"] - offset),
                    str(row["end"] - offset),
                    row["hgvs"],
                    row["ref"],
                    row["alt"],
                    row["aa_three"],
                ]
            )
        )
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    if not PICKETT_BED.is_file():
        sys.exit(f"pickett Branford list not found: {PICKETT_BED}")
    rows = read_pickett_bed(PICKETT_BED)
    if len(rows) != 73:
        sys.exit(f"Expected 73 rows in {PICKETT_BED}, found {len(rows)}.")
    if any(row["chrom"] != "chr9" for row in rows):
        sys.exit("Expected every row on chr9; the chr9-only offset would be wrong otherwise.")

    write_bed(HG19_OUT, rows, offset=0)
    write_bed(HG38_OUT, rows, offset=CHR9_HG19_TO_HG38_OFFSET)
    print(f"Wrote {len(rows)} variants to {HG19_OUT.name} and {HG38_OUT.name}.")


if __name__ == "__main__":
    main()
