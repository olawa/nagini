#!/usr/bin/env python3
"""Summarise Rindels CNV / ASCAT copy number segments and annotate with cytobands and target BED."""

import argparse
import bisect
import collections
import json
import os
import sys


class CytobandLookup:
    """Fast lookup of cytobands by genomic coordinates using binary search."""

    def __init__(self, cytoband_file):
        self.bands_by_chr = collections.defaultdict(list)
        if cytoband_file and os.path.isfile(cytoband_file):
            with open(cytoband_file, "r") as handle:
                for line in handle:
                    if line.startswith("#") or not line.strip():
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) >= 4:
                        chrom = parts[0]
                        try:
                            start = int(parts[1])
                            end = int(parts[2])
                            band = parts[3]
                            self.bands_by_chr[chrom].append((start, end, band))
                        except ValueError:
                            continue
            for chrom in self.bands_by_chr:
                self.bands_by_chr[chrom].sort(key=lambda x: x[0])

    def get_band(self, chrom, pos):
        if not self.bands_by_chr or chrom not in self.bands_by_chr:
            return ""
        bands = self.bands_by_chr[chrom]
        idx = bisect.bisect_right(bands, (pos, float("inf"), "")) - 1
        if 0 <= idx < len(bands):
            start, end, name = bands[idx]
            if start <= pos <= end:
                return name
        return ""


def read_target_bed(bed_path):
    """Read a target BED file into a mapping of chrom -> list of (start, end, label)."""
    intervals = collections.defaultdict(list)
    if not bed_path or not os.path.isfile(bed_path):
        return intervals
    with open(bed_path, "r") as handle:
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
            label = parts[3] if len(parts) > 3 else f"{chrom}:{start}-{end}"
            intervals[chrom].append((start, end, label))
    for chrom in intervals:
        intervals[chrom].sort(key=lambda x: x[0])
    return intervals


def find_overlapping_targets(target_intervals, chrom, start, end):
    """Return list of target labels overlapping the given interval."""
    if not target_intervals or chrom not in target_intervals:
        return []
    overlapping = []
    for t_start, t_end, label in target_intervals[chrom]:
        if max(start, t_start) < min(end, t_end):
            if label and label not in overlapping:
                overlapping.append(label)
    return overlapping


def parse_ascat_bed(ascat_bed_path):
    """
    Parse Rindels .ascat.bed format:
    Header:
      # ASCAT-style Heuristic Copy Number Segments (exploratory, internal markers)
      # Heuristic Purity: 1.000
      # Heuristic Ploidy: 1.250
      # Columns: chrom, start, end, CN_A, CN_B, BAF, logR, n_sites, state, n_het, n_hom, n_total
    Data lines:
      chr1 10394 248946189 1 2 0.217 -1.229 7 Amp 39 9 48
    """
    metadata = {
        "purity": None,
        "ploidy": None,
    }
    segments = []
    if not ascat_bed_path or not os.path.isfile(ascat_bed_path):
        return metadata, segments

    with open(ascat_bed_path, "r") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "Heuristic Purity:" in line:
                    try:
                        metadata["purity"] = float(line.split("Heuristic Purity:")[1].split()[0])
                    except (ValueError, IndexError):
                        pass
                elif "Heuristic Ploidy:" in line:
                    try:
                        metadata["ploidy"] = float(line.split("Heuristic Ploidy:")[1].split()[0])
                    except (ValueError, IndexError):
                        pass
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

            cn_a = None
            cn_b = None
            baf = None
            log_r = None
            state = "Diploid"

            if len(parts) >= 5:
                try:
                    cn_a = int(parts[3])
                    cn_b = int(parts[4])
                except ValueError:
                    pass
            if len(parts) >= 6:
                try:
                    baf = float(parts[5])
                except ValueError:
                    pass
            if len(parts) >= 7:
                try:
                    log_r = float(parts[6])
                except ValueError:
                    pass
            if len(parts) >= 9:
                state = parts[8]

            segments.append({
                "chrom": chrom,
                "start": start,
                "end": end,
                "cn_a": cn_a,
                "cn_b": cn_b,
                "baf": baf,
                "log_r": log_r,
                "state": state,
            })

    return metadata, segments


def classify_event(cn_a, cn_b, state, log_r):
    """Classify CNV event into GAIN, LOSS, LOH, or DIPLOID."""
    state_lower = str(state).lower()
    if cn_a is not None and cn_b is not None:
        total_cn = cn_a + cn_b
        if total_cn > 2:
            return "AMPLIFICATION" if total_cn >= 5 else "GAIN"
        if total_cn < 2:
            return "LOSS"
        if total_cn == 2:
            if cn_a == 0 or cn_b == 0 or "loh" in state_lower:
                return "LOH"
            return "DIPLOID"

    if "amp" in state_lower or "gain" in state_lower:
        return "GAIN"
    if "del" in state_lower or "loss" in state_lower:
        return "LOSS"
    if "loh" in state_lower:
        return "LOH"
    if log_r is not None:
        if log_r > 0.35:
            return "GAIN"
        if log_r < -0.35:
            return "LOSS"
    return "DIPLOID"


def process_cnv_segments(ascat_bed, target_bed, cytoband_file, sample="sample", sample_type="T"):
    """Annotate segments and return structured rows and summary metadata."""
    metadata, raw_segments = parse_ascat_bed(ascat_bed)
    targets = read_target_bed(target_bed)
    cytobands = CytobandLookup(cytoband_file)

    rows = []
    gains = 0
    losses = 0
    target_cnvs = 0

    for seg in raw_segments:
        chrom = seg["chrom"]
        start = seg["start"]
        end = seg["end"]
        length_kb = max(0.0, (end - start) / 1000.0)
        cn_a = seg["cn_a"]
        cn_b = seg["cn_b"]
        total_cn = (cn_a + cn_b) if (cn_a is not None and cn_b is not None) else None
        event_type = classify_event(cn_a, cn_b, seg["state"], seg["log_r"])

        band_start = cytobands.get_band(chrom, start)
        band_end = cytobands.get_band(chrom, end)
        if band_start and band_end:
            cyto_str = band_start if band_start == band_end else f"{band_start}-{band_end}"
        else:
            cyto_str = band_start or band_end or ""

        overlapping = find_overlapping_targets(targets, chrom, start, end)
        is_target = bool(overlapping)

        if event_type in ("GAIN", "AMPLIFICATION"):
            gains += 1
        elif event_type == "LOSS":
            losses += 1
        if is_target and event_type != "DIPLOID":
            target_cnvs += 1

        rows.append({
            "sample": sample,
            "type": sample_type,
            "chrom": chrom,
            "start": start,
            "end": end,
            "length_kb": round(length_kb, 1),
            "cn_a": cn_a if cn_a is not None else "",
            "cn_b": cn_b if cn_b is not None else "",
            "copy_number": total_cn if total_cn is not None else "",
            "log2_ratio": round(seg["log_r"], 3) if seg["log_r"] is not None else "",
            "event_type": event_type,
            "state": seg["state"],
            "cytoband": cyto_str,
            "is_target_region": is_target,
            "target_genes": ";".join(overlapping) if overlapping else ".",
        })

    summary = {
        "sample": sample,
        "type": sample_type,
        "purity": metadata.get("purity"),
        "ploidy": metadata.get("ploidy"),
        "total_segments": len(rows),
        "gains": gains,
        "losses": losses,
        "target_overlapping_cnvs": target_cnvs,
    }
    return rows, summary


def write_outputs(rows, summary, tsv_path, json_path):
    columns = [
        "sample",
        "type",
        "chrom",
        "start",
        "end",
        "length_kb",
        "cn_a",
        "cn_b",
        "copy_number",
        "log2_ratio",
        "event_type",
        "state",
        "cytoband",
        "is_target_region",
        "target_genes",
    ]
    os.makedirs(os.path.dirname(os.path.abspath(tsv_path)), exist_ok=True)
    with open(tsv_path, "w") as out:
        out.write("\t".join(columns) + "\n")
        for r in rows:
            out.write("\t".join(str(r.get(c, "")) for c in columns) + "\n")

    if json_path:
        os.makedirs(os.path.dirname(os.path.abspath(json_path)), exist_ok=True)
        with open(json_path, "w") as out:
            data = dict(summary)
            data["segments"] = rows
            json.dump(data, out, indent=2)


def main():
    parser = argparse.ArgumentParser(description="Summarise Rindels CNV / ASCAT copy number segments.")
    parser.add_argument("--ascat-bed", required=True, help="Input .ascat.bed from Rindels")
    parser.add_argument("--target-bed", default="", help="Design BED file for target overlap annotation")
    parser.add_argument("--cytoband", default="", help="UCSC cytoband TSV")
    parser.add_argument("--sample", default="sample", help="Sample name")
    parser.add_argument("--type", default="T", help="Sample type (T/N)")
    parser.add_argument("--output-tsv", required=True, help="Output summary TSV")
    parser.add_argument("--output-json", default="", help="Output summary JSON")
    args = parser.parse_args()

    rows, summary = process_cnv_segments(
        args.ascat_bed,
        args.target_bed,
        args.cytoband,
        sample=args.sample,
        sample_type=args.type,
    )
    write_outputs(rows, summary, args.output_tsv, args.output_json)


if "snakemake" in globals():
    target_bed_input = getattr(snakemake.input, "target_bed", "")
    target_bed_val = str(target_bed_input) if target_bed_input else ""
    cyto_input = getattr(snakemake.input, "cytoband", "")
    cyto_val = str(cyto_input) if cyto_input else ""

    rows, summary = process_cnv_segments(
        ascat_bed=str(snakemake.input.ascat),
        target_bed=target_bed_val,
        cytoband_file=cyto_val,
        sample=str(snakemake.wildcards.sample),
        sample_type=str(snakemake.wildcards.type),
    )
    json_path = str(snakemake.output.json) if hasattr(snakemake.output, "json") else ""
    write_outputs(rows, summary, str(snakemake.output.tsv), json_path)
elif __name__ == "__main__":
    main()
