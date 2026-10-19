#!/usr/bin/env python3
"""Sweep read quality scores (QS) across BCR-ABL1 kinase domain Branford panel positions.

Reads the BAM file over the ABL1 kinase domain in a single pass, extracts the read
`qs` tag (integer Q-score), stratifies per-base counts across distinct QS values,
and outputs:
1. A multi-page PDF with:
   - Page 1: Titration / Sensitivity overview curve (Median Depth, Mean Noise, Mutation VAF vs QS)
   - Pages 2..N: Full 3-panel stacked bar plot for each distinct QS score / threshold
2. An optional tidy TSV with all counts, depth, and VAFs across all sites and QS bins.
"""

from __future__ import annotations

import argparse
import csv
import math
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import numpy as np
    import pysam
except ImportError as err:
    sys.exit(f"Error: missing required package ({err}). Ensure pysam, matplotlib, and numpy are installed.")

from plot_branford_background import build_plot_figure, BASE_COLORS, STATUS_COLORS

AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}


def short_aa(aa: str) -> str:
    text = aa.strip()
    if text.startswith("p."):
        text = text[2:]
    for three, one in AA3_TO_AA1.items():
        text = text.replace(three, one)
    return text


def load_branford_bed(bed_path: Path) -> List[Dict[str, Any]]:
    sites: List[Dict[str, Any]] = []
    with bed_path.open() as handle:
        for line_idx, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t")
            if len(fields) < 7:
                fields = text.split()
            if len(fields) < 6:
                continue
            chrom = fields[0]
            start_pos = int(fields[1])
            name = fields[3]
            ref = fields[4].upper()
            alt = fields[5].upper()
            aa_change = fields[6] if len(fields) >= 7 else name
            short_name = short_aa(aa_change)
            sites.append({
                "chrom": chrom,
                "pos0": start_pos,
                "pos1": start_pos + 1,
                "name": name,
                "ref": ref,
                "alt": alt,
                "aa_change": aa_change,
                "mutation": short_name,
            })
    return sites


def collect_qs_pileup(
    bam_path: Path,
    sites: List[Dict[str, Any]],
    min_base_quality: int = 0,
    min_mapping_quality: int = 0,
) -> Tuple[str, Dict[int, Dict[Optional[int], Dict[str, int]]]]:
    """Collect base counts for each genomic position stratified by integer read QS."""
    by_pos: Dict[int, Dict[Optional[int], Dict[str, int]]] = {
        s["pos0"]: defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0, "total": 0})
        for s in sites
    }
    target_pos = set(by_pos.keys())
    if not target_pos:
        return "", by_pos

    with pysam.AlignmentFile(str(bam_path), "rb") as bam:
        ref_chrom = sites[0]["chrom"]
        chrom = ref_chrom
        if chrom not in bam.references:
            alt_c = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            if alt_c in bam.references:
                chrom = alt_c
            else:
                sys.stderr.write(f"Warning: neither {chrom} nor {alt_c} in BAM references!\n")
                return ref_chrom, by_pos

        min_pos = min(target_pos)
        max_pos = max(target_pos) + 1

        print(f"[Pileup] Scanning BAM {bam_path.name} across {chrom}:{min_pos:,}-{max_pos:,}...", file=sys.stderr)
        for col in bam.pileup(
            chrom, min_pos, max_pos, truncate=True,
            stepper="nofilter", max_depth=2000000,
            min_base_quality=0, ignore_overlaps=False
        ):
            pos0 = col.reference_pos
            if pos0 not in target_pos:
                continue

            best_by_read: Dict[Tuple[str, bool, bool], Tuple[Tuple[int, int], Any, str]] = {}
            for read in col.pileups:
                aln = read.alignment
                if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail or aln.is_secondary:
                    continue
                if (aln.mapping_quality or 0) < min_mapping_quality:
                    continue
                if read.is_del or read.is_refskip:
                    continue
                qpos = read.query_position
                if qpos is None:
                    continue
                if aln.query_qualities is not None and aln.query_qualities[qpos] < min_base_quality:
                    continue
                base = aln.query_sequence[qpos].upper()
                if base not in ("A", "C", "G", "T"):
                    continue

                # Deduplicate supplementary / split alignments of the same read
                key = (aln.query_name, aln.is_read1, aln.is_read2)
                rank = (0 if not aln.is_supplementary else 1, -(aln.mapping_quality or 0))
                if key not in best_by_read or rank < best_by_read[key][0]:
                    best_by_read[key] = (rank, aln, base)

            for _, aln, base in best_by_read.values():
                qs_int: Optional[int] = None
                if aln.has_tag("qs"):
                    try:
                        qs_val = float(aln.get_tag("qs"))
                        qs_int = int(qs_val)
                    except (ValueError, TypeError):
                        qs_int = None

                counts = by_pos[pos0][qs_int]
                counts[base] += 1
                counts["total"] += 1

    return ref_chrom, by_pos


def make_bins(qs_min: int, qs_max: int, mode: str) -> List[Dict[str, Any]]:
    bins = []
    if mode in ("binned", "both"):
        # Distinct integer bins
        bins.append({
            "type": "binned",
            "label": f"QS < {qs_min}",
            "short_label": f"<{qs_min}",
            "filter_fn": lambda q: q is None or q < qs_min,
        })
        for q in range(qs_min, qs_max + 1):
            bins.append({
                "type": "binned",
                "label": f"QS = {q}",
                "short_label": f"{q}",
                "filter_fn": (lambda target_q: lambda q: q == target_q)(q),
            })
        bins.append({
            "type": "binned",
            "label": f"QS > {qs_max}",
            "short_label": f">{qs_max}",
            "filter_fn": lambda q: q is not None and q > qs_max,
        })

    if mode in ("cumulative", "both"):
        # Cumulative thresholds
        for q in range(qs_min, qs_max + 1):
            bins.append({
                "type": "cumulative",
                "label": f"QS >= {q}",
                "short_label": f">={q}",
                "filter_fn": (lambda target_q: lambda q: q is not None and q >= target_q)(q),
            })
    return bins


def rows_for_bin(
    sites: List[Dict[str, Any]],
    by_pos: Dict[int, Dict[Optional[int], Dict[str, int]]],
    bin_def: Dict[str, Any],
    sample: str,
) -> List[Dict[str, str]]:
    rows = []
    filter_fn = bin_def["filter_fn"]

    for s in sites:
        pos0 = s["pos0"]
        counts_for_pos = by_pos.get(pos0, {})

        c_a = sum(c["A"] for q, c in counts_for_pos.items() if filter_fn(q))
        c_c = sum(c["C"] for q, c in counts_for_pos.items() if filter_fn(q))
        c_g = sum(c["G"] for q, c in counts_for_pos.items() if filter_fn(q))
        c_t = sum(c["T"] for q, c in counts_for_pos.items() if filter_fn(q))
        depth = c_a + c_c + c_g + c_t

        ref = s["ref"]
        alt = s["alt"]
        base_counts = {"A": c_a, "C": c_c, "G": c_g, "T": c_t}
        ad_ref = base_counts.get(ref, 0)
        ad_alt = base_counts.get(alt, 0)

        vaf_pct = (ad_alt / depth * 100.0) if depth > 0 else 0.0

        r = {
            "record_type": "SITE_AUDIT",
            "sample": sample,
            "mutation": s["mutation"],
            "aa_change": s["aa_change"],
            "name": s["name"],
            "chrom": s["chrom"],
            "pos": str(s["pos1"]),
            "ref": ref,
            "alt": alt,
            "depth": str(depth),
            "count_A": str(c_a),
            "count_C": str(c_c),
            "count_G": str(c_g),
            "count_T": str(c_t),
            "pct_A": f"{(c_a / depth * 100.0) if depth > 0 else 0.0:.2f}",
            "pct_C": f"{(c_c / depth * 100.0) if depth > 0 else 0.0:.2f}",
            "pct_G": f"{(c_g / depth * 100.0) if depth > 0 else 0.0:.2f}",
            "pct_T": f"{(c_t / depth * 100.0) if depth > 0 else 0.0:.2f}",
            "ad_ref": str(ad_ref),
            "ad_alt": str(ad_alt),
            "vaf_percent": f"{vaf_pct:.2f}",
            "status": "DETECTED" if vaf_pct >= 1.0 else ("SUBCLONAL" if vaf_pct >= 0.2 else "NEGATIVE"),
            "filter_reason": "PASS" if vaf_pct >= 0.2 else "NOISE",
            "qs_bin": bin_def["label"],
            "qs_type": bin_def["type"],
        }
        rows.append(r)
    return rows


def build_summary_figure(
    bin_summaries: List[Dict[str, Any]],
    sample: str,
) -> plt.Figure:
    """Build multi-panel summary curve across QS bins."""
    labels = [b["short_label"] for b in bin_summaries]
    x = np.arange(len(labels))
    depths = [b["median_depth"] for b in bin_summaries]
    noises = [b["mean_background_noise"] for b in bin_summaries]

    mut_vafs: Dict[str, List[float]] = {}
    for b in bin_summaries:
        for mut_name, vaf in b["mutation_vafs"].items():
            mut_vafs.setdefault(mut_name, []).append(vaf)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(max(14, len(labels) * 0.45), 10), sharex=True)
    fig.suptitle(f"BCR-ABL1 Kinase Domain | Quality Score (QS) Sweep Summary ({sample})", fontsize=14, fontweight="bold", y=0.98)

    # 1. Median Depth vs QS
    ax1.plot(x, depths, marker="o", color="#3470a3", linewidth=2.5, label="Median Depth across Sites")
    ax1.fill_between(x, depths, color="#3470a3", alpha=0.15)
    ax1.set_ylabel("Median Depth (x)", fontsize=11, fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.5)
    import matplotlib.ticker as ticker
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{int(y):,}"))
    ax1.set_title("Sequencing Coverage Retention across QS Bins", fontsize=11, loc="left", fontweight="semibold")
    ax1.legend(loc="upper right", framealpha=0.9)

    # 2. Mean Background Noise (%) vs QS
    ax2.plot(x, noises, marker="s", color="#d95f02", linewidth=2.2, label="Mean Non-Ref Noise (%) on Negative Sites")
    ax2.axhline(1.0, color="gray", linestyle=":", linewidth=1, label="1.0% Threshold")
    ax2.axhline(0.2, color="lightgray", linestyle=":", linewidth=0.8, label="0.2% Subclonal")
    ax2.set_ylabel("Background Noise (%)", fontsize=11, fontweight="bold")
    max_noise = max(noises) if noises else 0.5
    ax2.set_ylim(0, max(1.5, max_noise * 1.25))
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.set_title("Background Substitution Noise across Negative Sites", fontsize=11, loc="left", fontweight="semibold")
    ax2.legend(loc="upper right", framealpha=0.9)

    # 3. Detected Mutations VAF (%) vs QS
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628"]
    if mut_vafs:
        for idx, (m_name, v_list) in enumerate(mut_vafs.items()):
            col = colors[idx % len(colors)]
            ax3.plot(x, v_list, marker="^", color=col, linewidth=2.2, label=f"{m_name} VAF (%)")
    else:
        ax3.text(0.5, 0.5, "No variants with VAF >= 1.0% detected in sample", ha="center", va="center", transform=ax3.transAxes)

    ax3.set_ylabel("Mutation VAF (%)", fontsize=11, fontweight="bold")
    all_vafs = [v for v_list in mut_vafs.values() for v in v_list]
    max_vaf = max(all_vafs) if all_vafs else 5.0
    ax3.set_ylim(0, max(5.0, max_vaf * 1.25))
    ax3.grid(True, linestyle="--", alpha=0.5)
    ax3.set_title("Detected Resistance Mutation Stability across QS Bins", fontsize=11, loc="left", fontweight="semibold")
    if mut_vafs:
        ax3.legend(loc="upper right", framealpha=0.9)

    ax3.set_xticks(x)
    ax3.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax3.set_xlabel("Read Quality Score (QS) Bin / Cutoff", fontsize=11, fontweight="bold")

    plt.tight_layout()
    return fig


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", "-b", required=True, type=Path, help="Input aligned BAM file")
    parser.add_argument(
        "--bed", type=Path,
        default=Path("workflow/resources/bcr_abl/branford_positions_hg38.bed"),
        help="Branford positions BED file (default: hg38 BED)"
    )
    parser.add_argument("--output-pdf", "-o", required=True, type=Path, help="Path to output multi-page PDF")
    parser.add_argument("--output-tsv", type=Path, default=None, help="Optional path to output tidy summary TSV")
    parser.add_argument("--sample", default=None, help="Sample name (defaults to BAM filename stem)")
    parser.add_argument(
        "--mode", choices=["binned", "cumulative", "both"], default="binned",
        help="Sweep mode: 'binned' (distinct integer bins: <20, 20, 21... 30, >30), "
             "'cumulative' (>=20, >=21... >=30), or 'both' (default: 'binned')"
    )
    parser.add_argument("--qs-min", type=int, default=20, help="Minimum QS value to evaluate (default: 20)")
    parser.add_argument("--qs-max", type=int, default=30, help="Maximum QS value to evaluate (default: 30)")
    parser.add_argument("--min-base-quality", type=int, default=0, help="Minimum base quality to count base (default: 0)")
    parser.add_argument("--min-mapping-quality", type=int, default=0, help="Minimum mapping quality (default: 0)")
    parser.add_argument("--zoom-max", type=float, default=5.0, help="Y-axis maximum for background zoom panel (default: 5.0)")
    parser.add_argument("--color-mode", choices=["bases", "relative"], default="bases", help="Coloring mode for bar plots")
    return parser.parse_args(argv)


def main():
    args = parse_args()
    if not args.bam.is_file():
        sys.exit(f"Error: BAM file not found: {args.bam}")
    if not args.bed.is_file():
        sys.exit(f"Error: BED file not found: {args.bed}")

    sample = args.sample or args.bam.stem.replace("_T", "").replace(".bam", "")
    sites = load_branford_bed(args.bed)
    if not sites:
        sys.exit(f"Error: no valid sites loaded from BED {args.bed}")

    print(f"[Init] Loaded {len(sites)} Branford resistance positions.", file=sys.stderr)
    chrom, by_pos = collect_qs_pileup(
        bam_path=args.bam,
        sites=sites,
        min_base_quality=args.min_base_quality,
        min_mapping_quality=args.min_mapping_quality,
    )

    bins = make_bins(args.qs_min, args.qs_max, args.mode)
    print(f"[Sweep] Evaluating {len(bins)} QS bins ({args.mode} mode)...", file=sys.stderr)

    bin_data: List[Tuple[Dict[str, Any], List[Dict[str, str]]]] = []
    bin_summaries: List[Dict[str, Any]] = []
    all_tsv_rows: List[Dict[str, str]] = []

    # Identify primary mutations across the whole sample (to track on summary curve)
    # First check an unconstrained pool
    all_reads_fn = lambda q: True
    baseline_rows = rows_for_bin(sites, by_pos, {"filter_fn": all_reads_fn, "label": "ALL", "type": "all"}, sample)
    candidate_mutations: Set[str] = {
        r["mutation"] for r in baseline_rows
        if float(r["vaf_percent"]) >= 1.0
    }
    if not candidate_mutations:
        candidate_mutations = {r["mutation"] for r in baseline_rows if float(r["vaf_percent"]) >= 0.5}

    for b in bins:
        rows = rows_for_bin(sites, by_pos, b, sample)
        depths = [float(r["depth"]) for r in rows]
        med_depth = float(np.median(depths)) if depths else 0.0

        # Background noise on non-candidate sites
        non_candidate_noises = [
            float(r["vaf_percent"]) for r in rows
            if r["mutation"] not in candidate_mutations
        ]
        mean_noise = float(np.mean(non_candidate_noises)) if non_candidate_noises else 0.0

        # VAF of tracked candidate mutations
        mut_vafs = {
            r["mutation"]: float(r["vaf_percent"])
            for r in rows
            if r["mutation"] in candidate_mutations
        }

        bin_summary = {
            "label": b["label"],
            "short_label": b["short_label"],
            "type": b["type"],
            "median_depth": med_depth,
            "mean_background_noise": mean_noise,
            "mutation_vafs": mut_vafs,
        }
        bin_summaries.append(bin_summary)
        bin_data.append((b, rows))
        all_tsv_rows.extend(rows)

    # Generate multi-page PDF
    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.output_pdf) as pdf:
        print(f"[PDF] Rendering summary curve page...", file=sys.stderr)
        summary_fig = build_summary_figure(bin_summaries, sample)
        pdf.savefig(summary_fig, dpi=300)
        plt.close(summary_fig)

        for b, rows in bin_data:
            med_depth = int(np.median([float(r["depth"]) for r in rows]))
            title = f"BCR-ABL1 Kinase Domain ({sample}) | {b['label']} | Median Depth: {med_depth:,}x"
            print(f"[PDF] Rendering page: {b['label']} (Median Depth {med_depth:,}x)...", file=sys.stderr)

            fig = build_plot_figure(
                rows=rows,
                sample=sample,
                title=title,
                zoom_max=args.zoom_max,
                color_mode=args.color_mode,
                include_depth=True,
                depth_position="top",
            )
            if fig is not None:
                pdf.savefig(fig, dpi=300)
                plt.close(fig)

    print(f"[PDF] Successfully created multi-page QS sweep report: {args.output_pdf}", file=sys.stderr)

    # Optional TSV output
    if args.output_tsv:
        headers = [
            "sample", "qs_bin", "qs_type", "mutation", "aa_change", "chrom", "pos", "ref", "alt",
            "depth", "count_A", "count_C", "count_G", "count_T",
            "pct_A", "pct_C", "pct_G", "pct_T",
            "ad_ref", "ad_alt", "vaf_percent", "status", "filter_reason"
        ]
        args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
        with args.output_tsv.open("w") as handle:
            handle.write("\t".join(headers) + "\n")
            for r in all_tsv_rows:
                line = [r.get(h, "") for h in headers]
                handle.write("\t".join(str(val) for val in line) + "\n")
        print(f"[TSV] Saved tidy QS sweep counts summary: {args.output_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
