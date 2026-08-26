#!/usr/bin/env python3
"""Turn Dorado's per-read demux summary into a barcode yield table.

Dorado emits one row per read, which answers no question on its own.  What a
run needs to show is reads per barcode, which sample each barcode belongs to,
and how much of the flowcell went unclassified -- the three numbers that
explain a barcode with no data before anyone opens a BAM.

Output is MultiQC custom content, so the table lands in the run report next to
the rest of the QC rather than in a file nobody opens.
"""

import argparse
import csv
import gzip
import re
import sys
from collections import Counter
from pathlib import Path


BARCODE_COLUMNS = ("barcode", "barcode_arrangement", "alias", "barcode_name")
UNCLASSIFIED = "unclassified"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--units", type=Path, default=None)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--kit", default="")
    return parser.parse_args()


def open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path, "r", newline="")


def normalize_barcode(value: str) -> str:
    """barcode05, bc5 and 05 all name the same barcode."""

    text = str(value or "").strip().lower()
    if not text or text in ("na", "none", "-"):
        return ""
    if text.startswith(UNCLASSIFIED) or text == "unknown":
        return UNCLASSIFIED
    match = re.fullmatch(r"(?:barcode|bc)?0*(\d+)", text)
    return f"barcode{int(match.group(1)):02d}" if match else text


def count_barcodes(path: Path) -> Counter:
    with open_text(path) as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration:
            return Counter()
        lowered = [column.strip().lower() for column in header]
        index = next((lowered.index(name) for name in BARCODE_COLUMNS if name in lowered), None)
        if index is None:
            raise SystemExit(
                f"{path} has no barcode column; found: {', '.join(header)}."
            )
        counts = Counter()
        for row in reader:
            if index < len(row):
                barcode = normalize_barcode(row[index])
                if barcode:
                    counts[barcode] += 1
    return counts


def read_barcode_samples(path) -> dict:
    """Map barcode -> sample from units.tsv, when it is available."""

    if not path or not Path(path).is_file():
        return {}
    mapping = {}
    with open(path, newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            barcode = normalize_barcode(row.get("barcode_name", ""))
            sample = str(row.get("sample", "") or "").strip()
            if barcode and sample:
                mapping.setdefault(barcode, sample)
    return mapping


def build_rows(counts: Counter, barcode_samples: dict) -> list:
    classified = sum(value for key, value in counts.items() if key != UNCLASSIFIED)
    total = sum(counts.values())

    # Every expected barcode appears, including the ones with no reads: a
    # missing row and a zero row mean very different things here.
    barcodes = sorted(set(counts) | set(barcode_samples))
    rows = []
    for barcode in barcodes:
        if barcode == UNCLASSIFIED:
            continue
        reads = counts.get(barcode, 0)
        rows.append(
            {
                "Barcode": barcode,
                "Sample": barcode_samples.get(barcode, ""),
                "Reads": reads,
                "Percent of classified": round(100.0 * reads / classified, 2) if classified else 0.0,
                "Percent of run": round(100.0 * reads / total, 2) if total else 0.0,
                "Status": "OK" if reads else "NO READS",
            }
        )
    if UNCLASSIFIED in counts or total:
        reads = counts.get(UNCLASSIFIED, 0)
        rows.append(
            {
                "Barcode": UNCLASSIFIED,
                "Sample": "",
                "Reads": reads,
                "Percent of classified": "",
                "Percent of run": round(100.0 * reads / total, 2) if total else 0.0,
                "Status": "",
            }
        )
    return rows


COLUMNS = ["Barcode", "Sample", "Reads", "Percent of classified", "Percent of run", "Status"]


def write_mqc(path: Path, rows: list, kit: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    description = "Reads assigned per barcode by Dorado demux"
    if kit:
        description += f" (kit {kit})"
    with path.open("w") as handle:
        handle.write("# id: 'dorado_demux'\n")
        handle.write("# section_name: 'Dorado demultiplexing'\n")
        handle.write(f"# description: '{description}'\n")
        handle.write("# plot_type: 'table'\n")
        handle.write("# pconfig:\n")
        handle.write("#     id: 'dorado_demux_table'\n")
        handle.write("#     namespace: 'Dorado'\n")
        handle.write("#     col1_header: 'Barcode'\n")
        handle.write("\t".join(COLUMNS) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[column]) for column in COLUMNS) + "\n")


def main() -> None:
    args = parse_args()
    if not args.summary.is_file():
        raise SystemExit(f"Demux summary does not exist: {args.summary}")
    counts = count_barcodes(args.summary)
    barcode_samples = read_barcode_samples(args.units)
    rows = build_rows(counts, barcode_samples)
    write_mqc(args.output, rows, args.kit)

    empty = [row["Barcode"] for row in rows if row["Status"] == "NO READS"]
    if empty:
        print(
            "WARNING: no reads demultiplexed for: " + ", ".join(empty),
            file=sys.stderr,
        )
    unclassified = next((row for row in rows if row["Barcode"] == UNCLASSIFIED), None)
    if unclassified and isinstance(unclassified["Percent of run"], float):
        print(f"Unclassified reads: {unclassified['Percent of run']}% of the run.")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
