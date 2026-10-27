#!/usr/bin/env python3
"""Generate publication-ready figures for BCR-ABL1 dilution series (Table 7).

Compares Oxford Nanopore (Nagini Branford caller) against Illumina MiSeq across
spiked-in resistance mutations (T315I, Y253H at nominal 1%, 5%, 10%) across multiple
BCR-ABL1 background ratios (10%, 5%, 0.5%, 0.1%).

Creates:
1. Multi-panel figure (fig4_dilution_series.png / .svg):
   - Panel A: Linearity & Limit of Detection (Nominal vs Measured VAF %)
   - Panel B: Cross-Platform Concordance (Illumina MiSeq vs Oxford Nanopore)
2. Single-panel hero figure (fig4_dilution_series_single.png / .svg):
   - Focused cross-platform scatter plot with regression and low-VAF inset.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Visual Palette
RATIO_COLORS = {
    "10%": "#1E3A8A",    # Dark Royal Blue
    "5%": "#0284C7",     # Sky Blue
    "0.5%": "#059669",   # Emerald Green
    "0.1%": "#D97706",   # Amber / Orange
}

MUTATION_MARKERS = {
    "T315I": "o",   # Circle
    "Y253H": "^",   # Triangle
}


def parse_vaf(val: Any) -> Optional[float]:
    """Parse VAF string or float to percentage (0.0 to 100.0)."""
    if val is None:
        return None
    s = str(val).strip()
    if not s or s.lower() in ("-", "n/a", "n d", "nd", ".", "none", "ses ej"):
        return None
    s = s.replace("%", "").replace(",", ".").strip()
    try:
        return round(float(s), 3)
    except ValueError:
        return None


def parse_table7(table_path: Path) -> Dict[Tuple[str, float, str], float]:
    """Parse Table 7 text file for MiSeq benchmark VAFs.
    Returns: {(ratio_str, nominal_andel, mutation): miseq_vaf}
    """
    res: Dict[Tuple[str, float, str], float] = {}
    if not table_path.is_file():
        return res

    with open(table_path, "r", encoding="utf-8") as f:
        lines = [line.strip() for line in f if line.strip()]

    if len(lines) < 3:
        return res

    # Line 1: Tabell 7 \t T315I \t Y253H \t T315I \t Y253H \t T315I \t Y253H
    # Line 2: ratio \t 10% \t 10% \t 5% \t 5% \t 1% \t 1%
    mut_headers = [c.strip() for c in lines[0].split("\t")]
    andel_headers = [c.strip() for c in lines[1].split("\t")]

    # Find columns for (mutation, andel)
    col_specs: List[Tuple[int, str, float]] = []
    for idx in range(1, min(len(mut_headers), len(andel_headers))):
        mut = mut_headers[idx].strip()
        andel_str = andel_headers[idx].replace("%", "").replace(",", ".").strip()
        try:
            andel = float(andel_str)
            if mut in ("T315I", "Y253H"):
                col_specs.append((idx, mut, andel))
        except ValueError:
            continue

    # Data rows starting at line 2 (0-indexed)
    for line in lines[2:]:
        parts = [p.strip() for p in line.split("\t")]
        if not parts:
            continue
        raw_ratio = parts[0].replace(",", ".").strip()
        # Normalize ratio string: e.g. "0.50%" -> "0.5%", "0.10%" -> "0.1%"
        try:
            ratio_val = float(raw_ratio.replace("%", ""))
            ratio_str = f"{ratio_val:g}%"
        except ValueError:
            ratio_str = raw_ratio

        for col_idx, mut, andel in col_specs:
            if col_idx < len(parts):
                val = parse_vaf(parts[col_idx])
                if val is not None:
                    res[(ratio_str, andel, mut)] = val

    return res


def parse_miseq_results_dir(results_dir: Path) -> Dict[Tuple[str, float, str], float]:
    """Parse MiSeq *_summary.xlsx files in the results directory."""
    res: Dict[Tuple[str, float, str], float] = {}
    if not results_dir.is_dir():
        return res

    try:
        import openpyxl
    except ImportError:
        return res

    sample_id_map = {
        "01-1": ("0.1%", 1.0),
        "01-5": ("0.1%", 5.0),
        "01-10": ("0.1%", 10.0),
        "05-1": ("0.5%", 1.0),
        "05-5": ("0.5%", 5.0),
        "05-10": ("0.5%", 10.0),
        "5-1": ("5%", 1.0),
        "5-5": ("5%", 5.0),
        "5-10": ("5%", 10.0),
        "10-1": ("10%", 1.0),
        "10-5": ("10%", 5.0),
        "10-10": ("10%", 10.0),
    }

    for s_id, (ratio_str, andel) in sample_id_map.items():
        summary_file = results_dir / f"{s_id}_R_summary.xlsx"
        if not summary_file.is_file():
            continue
        wb = openpyxl.load_workbook(summary_file, data_only=True)
        if "Branford variants" not in wb.sheetnames:
            continue
        ws = wb["Branford variants"]
        for row in ws.iter_rows(values_only=True):
            if not row or not row[0]:
                continue
            mut = str(row[0]).strip()
            if mut in ("T315I", "Y253H") and len(row) > 6 and row[6] is not None:
                vaf = parse_vaf(row[6])
                if vaf is not None:
                    # In MiSeq xlsx, Allelic Freq is fraction (0.094 -> 9.4%)
                    if vaf <= 1.0 and "0." in str(row[6]):
                        vaf = round(vaf * 100, 3)
                    res[(ratio_str, andel, mut)] = vaf

    return res


def parse_ont_report(report_path: Path) -> Dict[Tuple[str, float, str], float]:
    """Parse Oxford Nanopore Nagini summary report Excel for Branford VAFs.
    Returns: {(ratio_str, nominal_andel, mutation): ont_vaf}
    """
    res: Dict[Tuple[str, float, str], float] = {}
    if not report_path.is_file():
        return res

    try:
        import openpyxl
    except ImportError:
        sys.stderr.write("Error: openpyxl is required to read Excel summary report.\n")
        return res

    wb = openpyxl.load_workbook(report_path, data_only=True)

    # Sample sheets follow pattern BCR-ABL (ratio<ratio>-andel<andel>)
    sample_regex = re.compile(r"^BCR-ABL \(ratio(?P<ratio>0\.1|0\.5|5|10)-andel(?P<andel>1|5|10)\)$")

    for sheet_name in wb.sheetnames:
        m = sample_regex.match(sheet_name)
        if not m:
            continue
        raw_ratio = m.group("ratio")
        raw_andel = m.group("andel")
        ratio_str = f"{raw_ratio}%"
        andel = float(raw_andel)

        ws = wb[sheet_name]
        in_branford_section = False
        for row in ws.iter_rows(values_only=True):
            if not row or not row[0]:
                continue
            first_cell = str(row[0]).strip()
            if "2. ABL1 KINASE DOMAIN RESISTANCE" in first_cell:
                in_branford_section = True
                continue
            if in_branford_section:
                if first_cell == "Mutation":
                    continue
                if "3." in first_cell or "COMPOUND" in first_cell:
                    break
                mut = first_cell
                if mut in ("T315I", "Y253H") and len(row) > 7 and row[7] is not None:
                    vaf = parse_vaf(row[7])
                    if vaf is not None:
                        res[(ratio_str, andel, mut)] = vaf

    return res


def build_combined_dataset(
    miseq_data: Dict[Tuple[str, float, str], float],
    ont_data: Dict[Tuple[str, float, str], float]
) -> List[Dict[str, Any]]:
    """Combine MiSeq and ONT datasets."""
    all_keys = sorted(set(miseq_data.keys()) | set(ont_data.keys()))
    dataset: List[Dict[str, Any]] = []

    for ratio_str, andel, mut in all_keys:
        m_vaf = miseq_data.get((ratio_str, andel, mut))
        o_vaf = ont_data.get((ratio_str, andel, mut))
        is_background = (ratio_str == "0.1%" and andel == 5.0)

        dataset.append({
            "ratio": ratio_str,
            "nominal_vaf": andel,
            "mutation": mut,
            "miseq_vaf": m_vaf,
            "ont_vaf": o_vaf,
            "is_background": is_background,
        })

    return dataset


def setup_plot_style() -> None:
    """Set publication-grade font and styling."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica Neue", "Arial"],
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "xtick.labelsize": 10.5,
        "ytick.labelsize": 10.5,
        "legend.fontsize": 10,
        "figure.titlesize": 15,
        "axes.linewidth": 1.2,
        "axes.edgecolor": "#334155",
        "grid.color": "#CBD5E1",
        "grid.linestyle": "--",
        "grid.alpha": 0.5,
    })


def create_two_panel_figure(
    dataset: List[Dict[str, Any]],
    output_path: Path,
    dpi: int = 300
) -> None:
    """Generate 2-panel publication figure:
    Panel A: Linearity & Limit of Detection (Nominal vs Measured VAF %)
    Panel B: Cross-Platform Concordance (Illumina MiSeq vs Oxford Nanopore)
    """
    setup_plot_style()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15.5, 6.8), constrained_layout=True)

    # -------------------------------------------------------------
    # PANEL A: Linearity & Limit of Detection
    # -------------------------------------------------------------
    ax1.set_title("A. Linearity & Limit of Detection (LOD)", weight="bold", pad=12)
    ax1.set_xlabel("Expected Nominal Mutation Fraction (%)", weight="bold")
    ax1.set_ylabel("Measured VAF (%) [Oxford Nanopore]", weight="bold")

    # Plot ideal line
    x_ideal = np.linspace(0, 12, 100)
    ax1.plot(x_ideal, x_ideal, color="#64748B", linestyle="--", linewidth=1.5,
             alpha=0.8, label="Ideal linearity ($y = x$)")

    # Group points by ratio and mutation
    ratio_order = ["10%", "5%", "0.5%", "0.1%"]

    for d in dataset:
        if d["ont_vaf"] is None:
            continue
        ratio = d["ratio"]
        mut = d["mutation"]
        nom = d["nominal_vaf"]
        vaf = d["ont_vaf"]
        color = RATIO_COLORS.get(ratio, "#334155")
        marker = MUTATION_MARKERS.get(mut, "o")

        # Jitter nominal slightly for visual separation between variants
        x_plot = nom - 0.15 if mut == "T315I" else nom + 0.15

        if d["is_background"]:
            # Highlight known pre-analytical background sample
            ax1.scatter(
                x_plot, vaf, color=color, marker=marker, s=110,
                edgecolors="#DC2626", linewidths=2.2, alpha=0.9, zorder=5
            )
            ax1.annotate(
                f"{mut}\n(Sample background\npre-analytical)",
                (x_plot, vaf),
                xytext=(x_plot + 0.8, vaf - 0.5),
                fontsize=8.5,
                color="#991B1B",
                arrowprops=dict(arrowstyle="->", color="#DC2626", lw=1.2)
            )
        else:
            ax1.scatter(
                x_plot, vaf, color=color, marker=marker, s=95,
                edgecolors="#0F172A", linewidths=1.2, alpha=0.9, zorder=4
            )

    # Inset for low-VAF range (0 - 2.5%)
    ax_ins = ax1.inset_axes([0.18, 0.54, 0.40, 0.40])
    ax_ins.plot([0, 2.5], [0, 2.5], color="#64748B", linestyle="--", lw=1.2, alpha=0.8)

    for d in dataset:
        if d["ont_vaf"] is None or d["is_background"]:
            continue
        if d["nominal_vaf"] <= 2.0:
            ratio = d["ratio"]
            mut = d["mutation"]
            nom = d["nominal_vaf"]
            vaf = d["ont_vaf"]
            color = RATIO_COLORS.get(ratio, "#334155")
            marker = MUTATION_MARKERS.get(mut, "o")
            x_plot = nom - 0.04 if mut == "T315I" else nom + 0.04
            ax_ins.scatter(
                x_plot, vaf, color=color, marker=marker, s=80,
                edgecolors="#0F172A", linewidths=1.0, alpha=0.95
            )

    ax_ins.set_xlim(0.5, 1.5)
    ax_ins.set_ylim(0.5, 1.8)
    ax_ins.set_xticks([1.0])
    ax_ins.set_yticks([1.0, 1.5])
    ax_ins.set_title("1.0% Low VAF Dynamic Range", fontsize=9.5, weight="bold")
    ax_ins.tick_params(labelsize=8.5)
    ax_ins.grid(True, linestyle=":", alpha=0.6)
    ax1.indicate_inset_zoom(ax_ins, edgecolor="#475569", lw=1.2)

    ax1.set_xlim(-0.5, 12.5)
    ax1.set_ylim(-0.5, 20.0)
    ax1.set_xticks([1, 5, 10])
    ax1.set_xticklabels(["1%", "5%", "10%"])
    ax1.grid(True, linestyle="--", alpha=0.5)

    # Custom legends for Panel A
    legend_ratios = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=RATIO_COLORS[r],
                   markeredgecolor="#0F172A", markersize=9, label=f"BCR-ABL1 ratio {r}")
        for r in ratio_order
    ]
    legend_muts = [
        plt.Line2D([0], [0], marker=MUTATION_MARKERS[m], color="#334155", markerfacecolor="#94A3B8",
                   markeredgecolor="#0F172A", markersize=9, label=f"Mutation: {m}")
        for m in ["T315I", "Y253H"]
    ]
    first_legend = ax1.legend(handles=legend_ratios, loc="lower right", framealpha=0.92, title="Background Ratio")
    ax1.add_artist(first_legend)
    ax1.legend(handles=legend_muts, loc="center left", bbox_to_anchor=(0.02, 0.45), framealpha=0.92)

    # -------------------------------------------------------------
    # PANEL B: Cross-Platform Concordance (MiSeq vs Nanopore)
    # -------------------------------------------------------------
    ax2.set_title("B. Cross-Platform Concordance: MiSeq vs. Nanopore", weight="bold", pad=12)
    ax2.set_xlabel("Illumina MiSeq VAF (%)", weight="bold")
    ax2.set_ylabel("Oxford Nanopore Branford VAF (%)", weight="bold")

    # Filter pairs with both measurements
    pairs = [d for d in dataset if d["miseq_vaf"] is not None and d["ont_vaf"] is not None]

    x_vals = np.array([p["miseq_vaf"] for p in pairs])
    y_vals = np.array([p["ont_vaf"] for p in pairs])

    # Statistics
    slope, intercept = np.polyfit(x_vals, y_vals, 1)
    r_matrix = np.corrcoef(x_vals, y_vals)
    r_val = r_matrix[0, 1]
    r2_val = r_val ** 2
    mae_val = np.mean(np.abs(x_vals - y_vals))

    # Plot identity line
    max_val = max(np.max(x_vals), np.max(y_vals)) + 2.0
    ax2.plot([0, max_val], [0, max_val], color="#DC2626", linestyle="--", linewidth=1.5,
             alpha=0.85, label="Identity ($y = x$)", zorder=2)

    # Plot regression line
    x_fit = np.linspace(0, max_val, 100)
    y_fit = slope * x_fit + intercept
    ax2.plot(x_fit, y_fit, color="#1E3A8A", linestyle="-", linewidth=2.0,
             alpha=0.9, label=f"Fit: $y = {slope:.2f}x {'+' if intercept >= 0 else '-'} {abs(intercept):.2f}$", zorder=3)

    # Plot points
    for p in pairs:
        ratio = p["ratio"]
        mut = p["mutation"]
        color = RATIO_COLORS.get(ratio, "#334155")
        marker = MUTATION_MARKERS.get(mut, "o")

        if p["is_background"]:
            ax2.scatter(
                p["miseq_vaf"], p["ont_vaf"], color=color, marker=marker, s=120,
                edgecolors="#DC2626", linewidths=2.2, alpha=0.95, zorder=5
            )
            ax2.annotate(
                f"Sample background\nMiSeq: {p['miseq_vaf']:.1f}%\nONT: {p['ont_vaf']:.1f}%",
                (p["miseq_vaf"], p["ont_vaf"]),
                xytext=(p["miseq_vaf"] - 6.5, p["ont_vaf"] + 1.2),
                fontsize=8.5,
                color="#991B1B",
                arrowprops=dict(arrowstyle="->", color="#DC2626", lw=1.2)
            )
        else:
            ax2.scatter(
                p["miseq_vaf"], p["ont_vaf"], color=color, marker=marker, s=95,
                edgecolors="#0F172A", linewidths=1.2, alpha=0.9, zorder=4
            )

    # Stats box
    stats_text = (
        f"Concordance Metrics (N = {len(pairs)}):\n"
        f"• Pearson r = {r_val:.4f}\n"
        f"• $R^2$ = {r2_val:.4f}\n"
        f"• MAE = {mae_val:.2f}%\n"
        f"• Slope = {slope:.3f}"
    )
    ax2.text(
        0.05, 0.95, stats_text, transform=ax2.transAxes,
        fontsize=10, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#F8FAFC", edgecolor="#CBD5E1", alpha=0.95)
    )

    ax2.set_xlim(0, 20.0)
    ax2.set_ylim(0, 20.0)
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.legend(loc="lower right", framealpha=0.92)

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    sys.stdout.write(f"Generated 2-panel dilution figure: {output_path}\n")


def create_single_panel_figure(
    dataset: List[Dict[str, Any]],
    output_path: Path,
    dpi: int = 300
) -> None:
    """Generate single-panel hero scatter plot of MiSeq vs Oxford Nanopore."""
    setup_plot_style()
    fig, ax = plt.subplots(figsize=(8.5, 7.5), constrained_layout=True)

    pairs = [d for d in dataset if d["miseq_vaf"] is not None and d["ont_vaf"] is not None]
    if not pairs:
        sys.stderr.write("No matching pairs found for single panel figure.\n")
        return

    x_vals = np.array([p["miseq_vaf"] for p in pairs])
    y_vals = np.array([p["ont_vaf"] for p in pairs])

    slope, intercept = np.polyfit(x_vals, y_vals, 1)
    r_matrix = np.corrcoef(x_vals, y_vals)
    r_val = r_matrix[0, 1]
    r2_val = r_val ** 2
    mae_val = np.mean(np.abs(x_vals - y_vals))

    ax.set_title("Cross-Platform Dilution Concordance: MiSeq vs. Nanopore", weight="bold", pad=14, fontsize=14)
    ax.set_xlabel("Illumina MiSeq Allele Frequency (%)", weight="bold", fontsize=12)
    ax.set_ylabel("Oxford Nanopore Allele Frequency (%)", weight="bold", fontsize=12)

    max_val = 20.0
    ax.plot([0, max_val], [0, max_val], color="#DC2626", linestyle="--", linewidth=1.6,
            alpha=0.85, label="Identity line ($y = x$)", zorder=2)

    x_fit = np.linspace(0, max_val, 100)
    y_fit = slope * x_fit + intercept
    ax.plot(x_fit, y_fit, color="#1E3A8A", linestyle="-", linewidth=2.2,
            alpha=0.9, label=f"Linear Fit ($y = {slope:.2f}x {'+' if intercept >= 0 else '-'} {abs(intercept):.2f}$)", zorder=3)

    for p in pairs:
        ratio = p["ratio"]
        mut = p["mutation"]
        color = RATIO_COLORS.get(ratio, "#334155")
        marker = MUTATION_MARKERS.get(mut, "o")

        if p["is_background"]:
            ax.scatter(
                p["miseq_vaf"], p["ont_vaf"], color=color, marker=marker, s=130,
                edgecolors="#DC2626", linewidths=2.2, alpha=0.95, zorder=5
            )
            ax.annotate(
                f"Sample background\nMiSeq: {p['miseq_vaf']:.1f}%\nONT: {p['ont_vaf']:.1f}%",
                (p["miseq_vaf"], p["ont_vaf"]),
                xytext=(p["miseq_vaf"] - 6.8, p["ont_vaf"] + 1.2),
                fontsize=9.0,
                color="#991B1B",
                arrowprops=dict(arrowstyle="->", color="#DC2626", lw=1.2)
            )
        else:
            ax.scatter(
                p["miseq_vaf"], p["ont_vaf"], color=color, marker=marker, s=100,
                edgecolors="#0F172A", linewidths=1.2, alpha=0.9, zorder=4
            )

    # Inset for low-VAF range (0.5% - 2.0%)
    ax_ins = ax.inset_axes([0.16, 0.52, 0.40, 0.40])
    ax_ins.plot([0, 2.5], [0, 2.5], color="#DC2626", linestyle="--", lw=1.2, alpha=0.8)
    x_ins_fit = np.linspace(0, 2.5, 50)
    ax_ins.plot(x_ins_fit, slope * x_ins_fit + intercept, color="#1E3A8A", lw=1.5, alpha=0.9)

    for p in pairs:
        if p["is_background"] or p["nominal_vaf"] > 2.0:
            continue
        ratio = p["ratio"]
        mut = p["mutation"]
        color = RATIO_COLORS.get(ratio, "#334155")
        marker = MUTATION_MARKERS.get(mut, "o")
        ax_ins.scatter(
            p["miseq_vaf"], p["ont_vaf"], color=color, marker=marker, s=85,
            edgecolors="#0F172A", linewidths=1.0, alpha=0.95
        )

    ax_ins.set_xlim(0.8, 1.8)
    ax_ins.set_ylim(0.7, 1.5)
    ax_ins.set_title("Low VAF Range (1.0% Spike-in)", fontsize=9.5, weight="bold")
    ax_ins.tick_params(labelsize=8.5)
    ax_ins.grid(True, linestyle=":", alpha=0.6)
    ax.indicate_inset_zoom(ax_ins, edgecolor="#475569", lw=1.2)

    stats_text = (
        f"Concordance Metrics (N = {len(pairs)}):\n"
        f"• Pearson r = {r_val:.4f}\n"
        f"• $R^2$ = {r2_val:.4f}\n"
        f"• MAE = {mae_val:.2f}%\n"
        f"• Slope = {slope:.3f}"
    )
    ax.text(
        0.05, 0.44, stats_text, transform=ax.transAxes,
        fontsize=10, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.6", facecolor="#F8FAFC", edgecolor="#CBD5E1", alpha=0.95)
    )

    ax.set_xlim(0, 20.0)
    ax.set_ylim(0, 20.0)
    ax.grid(True, linestyle="--", alpha=0.5)

    ratio_order = ["10%", "5%", "0.5%", "0.1%"]
    legend_ratios = [
        plt.Line2D([0], [0], marker="o", color="w", markerfacecolor=RATIO_COLORS[r],
                   markeredgecolor="#0F172A", markersize=9, label=f"Ratio {r}")
        for r in ratio_order
    ]
    legend_muts = [
        plt.Line2D([0], [0], marker=MUTATION_MARKERS[m], color="#334155", markerfacecolor="#94A3B8",
                   markeredgecolor="#0F172A", markersize=9, label=m)
        for m in ["T315I", "Y253H"]
    ]
    leg1 = ax.legend(handles=legend_ratios + legend_muts, loc="lower right", framealpha=0.92, title="Sample Details")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=dpi)
    plt.close()
    sys.stdout.write(f"Generated single-panel dilution figure: {output_path}\n")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Generate publication-ready figures for BCR-ABL1 dilution series."
    )
    parser.add_argument(
        "--report-xlsx",
        type=Path,
        default=Path("/Users/olwal516/data/bcr-abl/bcrabl_dilution_all_q14_800-2500_nagini_summary_report.xlsx"),
        help="Path to Nagini summary report Excel file for dilution run."
    )
    parser.add_argument(
        "--truth-table",
        type=Path,
        default=Path("/Users/olwal516/data/bcr-abl/table_7-dilution.txt"),
        help="Path to Table 7 dilution truth table."
    )
    parser.add_argument(
        "--miseq-dir",
        type=Path,
        default=Path("/Users/olwal516/data/bcr-abl/CGU_2022_12_BCR-ABL1/Workarea/230126_M01652_0666_000000000-GFYK4/Results"),
        help="Path to MiSeq results directory with *_summary.xlsx."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/poster/fig4_dilution_series.png"),
        help="Output path for multi-panel figure."
    )
    parser.add_argument(
        "--layout",
        choices=["all", "multi", "single"],
        default="all",
        help="Figure layout to generate."
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Figure resolution in DPI."
    )

    args = parser.parse_args(argv)

    # 1. Parse MiSeq data
    miseq_data = {}
    if args.miseq_dir and args.miseq_dir.is_dir():
        miseq_data = parse_miseq_results_dir(args.miseq_dir)
    if not miseq_data and args.truth_table and args.truth_table.is_file():
        miseq_data = parse_table7(args.truth_table)

    # 2. Parse ONT report
    ont_data = {}
    if args.report_xlsx and args.report_xlsx.is_file():
        ont_data = parse_ont_report(args.report_xlsx)

    if not ont_data:
        sys.stderr.write("Warning: No Oxford Nanopore dilution data found.\n")

    dataset = build_combined_dataset(miseq_data, ont_data)
    if not dataset:
        sys.stderr.write("Error: No combined dataset could be assembled.\n")
        return 1

    # 3. Generate outputs
    if args.layout in ("all", "multi"):
        create_two_panel_figure(dataset, args.output, dpi=args.dpi)
        if args.output.suffix.lower() == ".png":
            create_two_panel_figure(dataset, args.output.with_suffix(".svg"), dpi=args.dpi)

    if args.layout in ("all", "single"):
        if args.layout == "single":
            single_png = args.output
        else:
            single_png = args.output.parent / f"{args.output.stem}_single{args.output.suffix}"
        create_single_panel_figure(dataset, single_png, dpi=args.dpi)
        if single_png.suffix.lower() == ".png":
            create_single_panel_figure(dataset, single_png.with_suffix(".svg"), dpi=args.dpi)

    return 0


if __name__ == "__main__":
    sys.exit(main())
