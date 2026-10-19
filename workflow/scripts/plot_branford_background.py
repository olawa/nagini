#!/usr/bin/env python3
"""Generate stacked bar chart of per-base counts and background frequencies across Branford panel sites.

Reads `*.branford_all_sites.tsv` produced by `bcr_abl_branford_caller.py` and creates:
1. A stacked bar plot (PNG/PDF) showing base frequencies (A, C, G, T) or (Ref vs Alt vs Other).
   Features a dual-panel layout:
   - Upper panel: Overview (0 - 100%) showing overall base composition.
   - Lower panel: Background noise zoom (0 - 5%) highlighting low-frequency errors and subclonal mutations.
2. An optional clean TSV summary with per-base counts and percentages for each Branford site.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    plt = None
    np = None


BASE_COLORS = {
    "A": "#2ca02c",  # Green
    "C": "#1f77b4",  # Blue
    "G": "#ff7f0e",  # Orange / Yellow-orange
    "T": "#d62728",  # Red
}

STATUS_COLORS = {
    "REF": "#90caf9",    # Light Blue
    "ALT": "#e53935",    # Red
    "OTHER_1": "#ffa726",# Orange
    "OTHER_2": "#ab47bc",# Purple
}


def parse_args(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--all-sites", "-i", required=True, type=Path,
        help="Input reports/bcr_abl/{sample}_{type}.branford_all_sites.tsv"
    )
    parser.add_argument(
        "--output-png", "-o", required=True, type=Path,
        help="Path to output PNG image"
    )
    parser.add_argument(
        "--output-tsv", type=Path, default=None,
        help="Optional path to output simplified base counts TSV"
    )
    parser.add_argument(
        "--sample", default=None,
        help="Sample name (defaults to metadata in all_sites file)"
    )
    parser.add_argument(
        "--title", default=None,
        help="Plot title"
    )
    parser.add_argument(
        "--zoom-max", type=float, default=5.0,
        help="Y-axis maximum for the background zoom panel in percent (default: 5.0)"
    )
    parser.add_argument(
        "--color-mode", choices=["bases", "relative"], default="bases",
        help="Color bars by absolute nucleotide ('bases': A=green, C=blue, G=orange, T=red) "
             "or relative allele ('relative': Ref=light blue, Alt=red, Others=orange/purple)"
    )
    parser.add_argument(
        "--no-depth-panel", action="store_true", default=False,
        help="Disable the total used depth / coverage subplot panel"
    )
    parser.add_argument(
        "--depth-position", choices=["top", "bottom"], default="top",
        help="Position of the depth panel ('top' or 'bottom', default: 'top')"
    )
    return parser.parse_args(argv)


def load_all_sites(path: Path) -> Tuple[Optional[str], List[Dict[str, str]]]:
    sample_name = None
    rows = []
    with path.open() as handle:
        first_line = ""
        for line in handle:
            line_str = line.strip()
            if not line_str:
                continue
            if line_str.startswith("#"):
                if "sample=" in line_str and not sample_name:
                    sample_name = line_str.split("sample=")[1].split()[0]
                continue
            first_line = line_str
            break

        if not first_line:
            return sample_name, rows

        fieldnames = first_line.split("\t")
        reader = csv.DictReader(handle, fieldnames=fieldnames, delimiter="\t")
        is_rindels = "filtered_depth" in fieldnames or ("A" in fieldnames and "A_fwd" in fieldnames)

        for r in reader:
            if not is_rindels and r.get("record_type") != "SITE_AUDIT":
                continue
            rows.append(r)

    if not sample_name:
        stem = path.name
        for suffix in [".branford_all_sites.tsv", ".pileup_counts.tsv", ".target_raw_counts.tsv", ".tsv"]:
            if stem.endswith(suffix):
                stem = stem[:-len(suffix)]
                break
        sample_name = stem

    return sample_name, rows


def build_plot_figure(
    rows: List[Dict[str, str]],
    sample: str = "",
    title: Optional[str] = None,
    zoom_max: float = 5.0,
    color_mode: str = "bases",
    include_depth: bool = True,
    depth_position: str = "top",
) -> Optional[Any]:
    if plt is None:
        sys.exit("Error: matplotlib and numpy are required to generate the plot.")

    if not rows:
        sys.stderr.write("Warning: no SITE_AUDIT rows to plot.\n")
        return None

    # Extract variant names and per-base percentages
    labels = []
    depths = []
    data_a = []
    data_c = []
    data_g = []
    data_t = []

    # Relative mode data
    data_ref = []
    data_alt = []
    data_other = []

    for r in rows:
        mut = r.get("mutation") or r.get("name", "Unknown")
        ref = (r.get("ref") or r.get("ref_base") or "").upper()
        alt = (r.get("alt") or r.get("alt_base") or "").upper()
        depth = float(r.get("depth") or r.get("filtered_depth") or r.get("raw_depth") or 0)
        depths.append(depth)
        labels.append(f"{mut}\n({ref}>{alt})")

        # Counts or percentages
        pct_a = float(r.get("pct_A", 0) or 0)
        pct_c = float(r.get("pct_C", 0) or 0)
        pct_g = float(r.get("pct_G", 0) or 0)
        pct_t = float(r.get("pct_T", 0) or 0)

        c_a = float(r.get("count_A") or r.get("A") or 0)
        c_c = float(r.get("count_C") or r.get("C") or 0)
        c_g = float(r.get("count_G") or r.get("G") or 0)
        c_t = float(r.get("count_T") or r.get("T") or 0)

        # Fallback if pct_A wasn't populated
        if pct_a == 0 and pct_c == 0 and pct_g == 0 and pct_t == 0 and depth > 0:
            pct_a = c_a / depth * 100.0
            pct_c = c_c / depth * 100.0
            pct_g = c_g / depth * 100.0
            pct_t = c_t / depth * 100.0

        data_a.append(pct_a)
        data_c.append(pct_c)
        data_g.append(pct_g)
        data_t.append(pct_t)

        if "ad_ref" in r:
            c_ref = float(r.get("ad_ref") or 0)
        elif ref in ("A", "C", "G", "T"):
            c_ref = float(r.get(f"count_{ref}") or r.get(ref) or 0)
        else:
            c_ref = float(r.get("ref_count") or 0)

        if "ad_alt" in r:
            c_alt = float(r.get("ad_alt") or 0)
        elif alt in ("A", "C", "G", "T"):
            c_alt = float(r.get(f"count_{alt}") or r.get(alt) or 0)
        else:
            c_alt = float(r.get("alt_count") or 0)

        pct_ref = c_ref / depth * 100.0 if depth > 0 else 0.0
        pct_alt = c_alt / depth * 100.0 if depth > 0 else 0.0
        pct_other = max(0.0, 100.0 - pct_ref - pct_alt) if depth > 0 else 0.0

        data_ref.append(pct_ref)
        data_alt.append(pct_alt)
        data_other.append(pct_other)

    x = np.arange(len(labels))
    width = 0.75

    if include_depth:
        fig_height = 12
        if depth_position == "bottom":
            fig, (ax_top, ax_bot, ax_cov) = plt.subplots(
                3, 1, figsize=(max(14, len(labels) * 0.28), fig_height),
                gridspec_kw={"height_ratios": [1.0, 1.4, 0.8]},
                sharex=True
            )
            labels_ax = ax_cov
        else:
            fig, (ax_cov, ax_top, ax_bot) = plt.subplots(
                3, 1, figsize=(max(14, len(labels) * 0.28), fig_height),
                gridspec_kw={"height_ratios": [0.8, 1.0, 1.4]},
                sharex=True
            )
            labels_ax = ax_bot
    else:
        fig, (ax_top, ax_bot) = plt.subplots(
            2, 1, figsize=(max(14, len(labels) * 0.28), 9),
            gridspec_kw={"height_ratios": [1, 1.4]},
            sharex=True
        )
        ax_cov = None
        labels_ax = ax_bot

    plot_title = title or f"BCR-ABL1 Kinase Domain Base Composition & Background ({sample or 'Sample'})"
    fig.suptitle(plot_title, fontsize=14, fontweight="bold", y=0.98)

    # 1. Total Depth Plot (if enabled)
    if ax_cov is not None:
        ax_cov.bar(x, depths, width, color="#3470a3", edgecolor="none", alpha=0.9)
        ax_cov.set_ylabel("Used Depth", fontsize=11, fontweight="bold")
        ax_cov.grid(axis="y", linestyle="--", alpha=0.5)

        import matplotlib.ticker as ticker
        ax_cov.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{int(y):,}"))

        max_d = max(depths) if depths else 0
        if max_d > 0:
            ax_cov.set_ylim(0, max_d * 1.15)
            med_d = float(np.median(depths))
            ax_cov.axhline(
                med_d, color="#d95f02", linestyle="--", linewidth=1.2,
                label=f"Median Depth ({int(med_d):,}x)"
            )
            ax_cov.legend(loc="upper right", framealpha=0.9, fontsize=9)
        else:
            ax_cov.set_ylim(0, 1)

        ax_cov.set_title("Total Used Sequencing Depth per Site", fontsize=11, loc="left", fontweight="semibold")

    # 2. Composition and Zoom Plots
    if color_mode == "bases":
        # Stacked A, C, G, T
        bottom_top = np.zeros(len(labels))
        for base, values in [("A", data_a), ("C", data_c), ("G", data_g), ("T", data_t)]:
            v = np.array(values)
            ax_top.bar(x, v, width, bottom=bottom_top, label=f"Base {base}", color=BASE_COLORS[base], edgecolor="none")
            bottom_top += v

        # Bottom panel: zoom in on non-ref bases
        non_ref_a = [pa if (r.get("ref") or r.get("ref_base", "")).upper() != "A" else 0.0 for pa, r in zip(data_a, rows)]
        non_ref_c = [pc if (r.get("ref") or r.get("ref_base", "")).upper() != "C" else 0.0 for pc, r in zip(data_c, rows)]
        non_ref_g = [pg if (r.get("ref") or r.get("ref_base", "")).upper() != "G" else 0.0 for pg, r in zip(data_g, rows)]
        non_ref_t = [pt if (r.get("ref") or r.get("ref_base", "")).upper() != "T" else 0.0 for pt, r in zip(data_t, rows)]

        bottom_bot = np.zeros(len(labels))
        for base, values in [("A", non_ref_a), ("C", non_ref_c), ("G", non_ref_g), ("T", non_ref_t)]:
            v = np.array(values)
            ax_bot.bar(x, v, width, bottom=bottom_bot, label=f"Non-Ref {base}", color=BASE_COLORS[base], edgecolor="none")
            bottom_bot += v

    else:
        # Relative: Ref vs Alt vs Others
        ax_top.bar(x, data_ref, width, label="Reference Base", color=STATUS_COLORS["REF"])
        ax_top.bar(x, data_alt, width, bottom=data_ref, label="Candidate Alt Base", color=STATUS_COLORS["ALT"])
        ax_top.bar(x, data_other, width, bottom=np.array(data_ref) + np.array(data_alt), label="Other Non-Candidate Bases", color=STATUS_COLORS["OTHER_1"])

        # Bottom panel: zoom on Alt and Other
        ax_bot.bar(x, data_alt, width, label="Candidate Alt Base", color=STATUS_COLORS["ALT"])
        ax_bot.bar(x, data_other, width, bottom=data_alt, label="Other Noise Bases", color=STATUS_COLORS["OTHER_1"])

    # Configure top composition panel
    ax_top.set_ylabel("Allele Fraction (%)", fontsize=11, fontweight="bold")
    ax_top.set_ylim(0, 100)
    ax_top.grid(axis="y", linestyle="--", alpha=0.5)
    ax_top.legend(loc="upper right", framealpha=0.9, fontsize=9)
    ax_top.set_title("Full Composition (0 - 100%)", fontsize=11, loc="left", fontweight="semibold")

    # Configure bottom background zoom panel
    ax_bot.set_ylabel("Background / Alt AF (%)", fontsize=11, fontweight="bold")
    ax_bot.set_ylim(0, zoom_max)
    ax_bot.grid(axis="y", linestyle="--", alpha=0.5)
    ax_bot.axhline(1.0, color="gray", linestyle=":", linewidth=1, label="1.0% Threshold")
    ax_bot.axhline(0.2, color="lightgray", linestyle=":", linewidth=0.8, label="0.2% Subclonal")
    ax_bot.legend(loc="upper right", framealpha=0.9, fontsize=9)
    ax_bot.set_title(f"Background Error & Alternate Allele Zoom (0 - {zoom_max:g}%)", fontsize=11, loc="left", fontweight="semibold")

    # Set x-ticks and tick labels only on the bottom-most axis
    for ax in [ax_cov, ax_top, ax_bot]:
        if ax is not None and ax is not labels_ax:
            ax.tick_params(labelbottom=False)

    labels_ax.set_xticks(x)
    labels_ax.set_xticklabels(labels, rotation=90, fontsize=8, ha="center")
    labels_ax.set_xlim(-0.8, len(labels) - 0.2)

    plt.tight_layout()
    return fig


def generate_plot(
    rows: List[Dict[str, str]],
    output_png: Path,
    sample: str = "",
    title: Optional[str] = None,
    zoom_max: float = 5.0,
    color_mode: str = "bases",
    include_depth: bool = True,
    depth_position: str = "top",
):
    fig = build_plot_figure(
        rows=rows,
        sample=sample,
        title=title,
        zoom_max=zoom_max,
        color_mode=color_mode,
        include_depth=include_depth,
        depth_position=depth_position,
    )
    if fig is None:
        return
    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300)
    plt.close(fig)
    print(f"[Plot] Saved background base frequency plot: {output_png}", file=sys.stderr)


def write_summary_tsv(rows: List[Dict[str, str]], output_tsv: Path, sample: str):
    headers = [
        "sample", "mutation", "aa_change", "chrom", "pos", "ref", "alt",
        "depth", "count_A", "count_C", "count_G", "count_T",
        "pct_A", "pct_C", "pct_G", "pct_T",
        "ad_ref", "ad_alt", "vaf_percent", "status", "filter_reason"
    ]
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with output_tsv.open("w") as handle:
        handle.write("\t".join(headers) + "\n")
        for r in rows:
            depth = float(r.get("depth") or r.get("filtered_depth") or r.get("raw_depth") or 0)
            c_a = float(r.get("count_A") or r.get("A") or 0)
            c_c = float(r.get("count_C") or r.get("C") or 0)
            c_g = float(r.get("count_G") or r.get("G") or 0)
            c_t = float(r.get("count_T") or r.get("T") or 0)
            ref = (r.get("ref") or r.get("ref_base") or "").upper()
            alt = (r.get("alt") or r.get("alt_base") or "").upper()

            pct_a = float(r.get("pct_A") or (c_a / depth * 100.0 if depth > 0 else 0.0))
            pct_c = float(r.get("pct_C") or (c_c / depth * 100.0 if depth > 0 else 0.0))
            pct_g = float(r.get("pct_G") or (c_g / depth * 100.0 if depth > 0 else 0.0))
            pct_t = float(r.get("pct_T") or (c_t / depth * 100.0 if depth > 0 else 0.0))

            if "ad_ref" in r:
                ad_ref = r.get("ad_ref", "")
            elif ref in ("A", "C", "G", "T"):
                ad_ref = int(r.get(f"count_{ref}") or r.get(ref) or 0)
            else:
                ad_ref = int(r.get("ref_count") or 0)

            if "ad_alt" in r:
                ad_alt = r.get("ad_alt", "")
            elif alt in ("A", "C", "G", "T"):
                ad_alt = int(r.get(f"count_{alt}") or r.get(alt) or 0)
            else:
                ad_alt = int(r.get("alt_count") or 0)

            if "vaf_percent" in r:
                vaf_pct = r.get("vaf_percent", "")
            elif "vaf" in r:
                vaf_pct = f"{float(r.get('vaf', 0)) * 100.0:.2f}"
            elif depth > 0:
                vaf_pct = f"{float(ad_alt) / depth * 100.0:.2f}"
            else:
                vaf_pct = "0.0"

            line = [
                sample or r.get("sample", ""),
                r.get("mutation") or r.get("name", ""),
                r.get("aa_change", ""),
                r.get("chrom", ""),
                r.get("pos", ""),
                ref,
                alt,
                int(depth),
                int(c_a),
                int(c_c),
                int(c_g),
                int(c_t),
                f"{pct_a:.2f}",
                f"{pct_c:.2f}",
                f"{pct_g:.2f}",
                f"{pct_t:.2f}",
                ad_ref,
                ad_alt,
                vaf_pct,
                r.get("status") or r.get("vcf_filter", ""),
                r.get("filter_reason", ""),
            ]
            handle.write("\t".join(str(val) for val in line) + "\n")
    print(f"[TSV] Saved base counts summary: {output_tsv}", file=sys.stderr)


def main():
    args = parse_args()
    sample_meta, rows = load_all_sites(args.all_sites)
    sample = args.sample or sample_meta or args.all_sites.stem.split(".")[0]

    generate_plot(
        rows=rows,
        output_png=args.output_png,
        sample=sample,
        title=args.title,
        zoom_max=args.zoom_max,
        color_mode=args.color_mode,
        include_depth=not args.no_depth_panel,
        depth_position=args.depth_position,
    )

    if args.output_tsv:
        write_summary_tsv(rows, args.output_tsv, sample)


if __name__ == "__main__":
    main()
