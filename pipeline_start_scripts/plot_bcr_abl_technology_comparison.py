#!/usr/bin/env python3
"""Generate publication-ready comparison figures comparing Nagini BCR-ABL1 callers
against Illumina MiSeq and PacBio Sequel.

Creates:
1. Panel A: Cross-Platform Linear Correlation & Inset (0-3% low VAF range)
2. Panel B: Grouped VAF Concordance across Dynamic Range (Low, Intermediate, High)
3. Panel C: Bland-Altman Residuals / Delta VAF (%)
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


# Standard Branford amino acid 3-letter to 1-letter mapping
THREE_TO_ONE = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
    "Glu": "E", "Gln": "Q", "Gly": "G", "His": "H", "Ile": "I",
    "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
    "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"
}

# Caller styling
CALLER_COLORS = {
    "MiSeq": "#334155",        # Charcoal / Slate
    "Sequel": "#0284C7",       # Sky Blue / Cyan
    "Branford": "#1E3A8A",     # Royal Navy Blue
    "bcftools": "#D97706",     # Amber / Orange
    "Rindels": "#7C3AED",      # Violet / Purple
    "VarDictRust": "#DC2626",  # Crimson / Red
}

CALLER_MARKERS = {
    "Branford": "o",
    "bcftools": "s",
    "Rindels": "^",
    "VarDictRust": "D",
}


def normalize_aa(s: str) -> str:
    """Normalize amino acid change representation (e.g. p.Met351Thr -> M351T)."""
    if not s:
        return ""
    clean = str(s).strip()
    m = re.match(r"^p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})$", clean)
    if m:
        ref, pos, alt = m.groups()
        return f"{THREE_TO_ONE.get(ref, ref)}{pos}{THREE_TO_ONE.get(alt, alt)}"
    m2 = re.match(r"^([A-Z])(\d+)([A-Z])$", clean)
    if m2:
        return clean
    return clean


def parse_vaf(val: Any) -> Optional[float]:
    """Parse VAF string or float to percentage (0.0 to 100.0)."""
    if val is None:
        return None
    s = str(val).strip()
    if not s or s in ("-", "N/A", "n d", "nd", ".", "None", "ses ej"):
        return None
    s = s.replace("%", "").replace(",", ".").strip()
    try:
        num = float(s)
        return round(num, 2)
    except ValueError:
        return None


def load_branford_bed(bed_path: Path) -> Dict[Tuple[str, int, str], str]:
    """Load Branford BED file to map (chrom, pos, alt) to AA change."""
    coords_to_aa: Dict[Tuple[str, int, str], str] = {}
    if not bed_path.is_file():
        return coords_to_aa

    with open(bed_path, "r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 7:
                chrom = parts[0]
                pos = int(parts[2])  # 1-based end position
                alt = parts[5]
                aa_raw = parts[6]
                aa_short = normalize_aa(aa_raw)
                coords_to_aa[(chrom, pos, alt)] = aa_short
                coords_to_aa[(chrom.removeprefix("chr"), pos, alt)] = aa_short

    return coords_to_aa


def load_caller_results(
    tsv_path: Path, coords_to_aa: Optional[Dict[Tuple[str, int, str], str]] = None
) -> Dict[Tuple[str, str], Dict[str, float]]:
    """Load caller VAFs from summary TSV.
    Returns: {(dna_number, aa_change): {caller: vaf}}
    """
    coords_to_aa = coords_to_aa or {}
    results: Dict[Tuple[str, str], Dict[str, float]] = {}
    if not tsv_path.is_file():
        return results

    with open(tsv_path, "r", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for r in reader:
            dna = r.get("DNA Number", "").strip()
            caller = r.get("Caller", "").strip()
            chrom = r.get("Chr", "").strip()
            alt = r.get("Alt", "").strip()
            try:
                pos = int(r.get("Position", 0))
            except ValueError:
                pos = 0

            aa = normalize_aa(r.get("Consequence / AA Change", "").strip())
            if not aa and pos:
                aa = coords_to_aa.get((chrom, pos, alt), "") or coords_to_aa.get((chrom.removeprefix("chr"), pos, alt), "")

            vaf_str = r.get("VAF", "").strip()
            vaf = parse_vaf(vaf_str)
            if not dna or not caller or not aa or vaf is None:
                continue

            key = (dna, aa)
            if key not in results:
                results[key] = {}
            results[key][caller] = vaf

    return results


def load_truth_table(tsv_path: Path) -> List[Dict[str, Any]]:
    """Parse Table 6 truth table for MiSeq & Sequel AFs."""
    records: List[Dict[str, Any]] = []
    if not tsv_path.is_file():
        return records

    with open(tsv_path, "r", encoding="utf-8") as f:
        lines = [l.strip() for l in f if l.strip()]

    start_idx = 0
    header_cols = []
    for idx, l in enumerate(lines):
        if l.startswith("Provnr") or "aa-variant" in l or l.startswith("Sample ID"):
            start_idx = idx + 1
            header_cols = [c.strip().lower() for c in l.split("\t")]

    def split_vars_and_afs(var_s: str, af_s: str) -> List[Tuple[str, Optional[float]]]:
        v_clean = var_s.strip()
        if not v_clean or v_clean.lower() in ("ingen variant", "ingen variant*", "ses ej", "låg cov", "snv i bcr"):
            return []
        vars_list = [v.strip() for v in re.split(r"[/,\n]+", v_clean) if v.strip()]
        afs_list = [a.strip() for a in re.split(r"[/;\n\s]+", af_s) if a.strip()]
        res = []
        for i, v in enumerate(vars_list):
            af_val = parse_vaf(afs_list[i]) if i < len(afs_list) else None
            res.append((normalize_aa(v), af_val))
        return res

    has_sample_col = False
    sample_col_idx = None
    for c_idx, h in enumerate(header_cols):
        if h in ("sample id", "sample_id", "sample", "bcrabl-bc", "bcrabl_bc", "barcode"):
            has_sample_col = True
            sample_col_idx = c_idx
            break

    for l in lines[start_idx:]:
        parts = [p.strip() for p in l.split("\t")]
        if len(parts) < 8:
            continue

        provnr = parts[0]
        sample_id = ""

        if has_sample_col and sample_col_idx is not None and sample_col_idx < len(parts):
            sample_id = parts[sample_col_idx]
            clean_parts = [p for i, p in enumerate(parts) if i != sample_col_idx]
            seq_var_str = clean_parts[2]
            seq_af_str = clean_parts[3]
            mis_var_str = clean_parts[5] if len(clean_parts) > 5 else ""
            mis_af_str = clean_parts[6] if len(clean_parts) > 6 else ""
        elif len(parts) >= 9 and (parts[1].lower().startswith("bcrabl-") or parts[1].lower().startswith("barcode")):
            sample_id = parts[1]
            seq_var_str = parts[3]
            seq_af_str = parts[4]
            mis_var_str = parts[6]
            mis_af_str = parts[7]
        else:
            seq_var_str = parts[2]
            seq_af_str = parts[3]
            mis_var_str = parts[5]
            mis_af_str = parts[6]
            if provnr.lower().startswith("bcrabl-") or provnr.lower().startswith("barcode"):
                sample_id = provnr

        seq_items = dict(split_vars_and_afs(seq_var_str, seq_af_str))
        mis_items = dict(split_vars_and_afs(mis_var_str, mis_af_str))

        all_vars = list(dict.fromkeys(list(seq_items.keys()) + list(mis_items.keys())))
        for v in all_vars:
            records.append({
                "dna_number": provnr,
                "sample_id": sample_id,
                "variant": v,
                "sequel_af": seq_items.get(v),
                "miseq_af": mis_items.get(v),
            })

    return records


def build_plot_dataset(
    truth_records: List[Dict[str, Any]],
    caller_results: Dict[Tuple[str, str], Dict[str, float]],
    sample_map: Optional[Dict[str, str]] = None
) -> List[Dict[str, Any]]:
    """Merge truth table and caller results into unified plot records."""
    dataset = []
    sample_map = sample_map or {}

    # 1. Add truth records
    seen_keys = set()
    for tr in truth_records:
        dna = tr.get("dna_number") or tr.get("provnr", "")
        sid = tr.get("sample_id", "") or sample_map.get(dna, "")
        if not sid and (dna.lower().startswith("bcrabl-") or dna.lower().startswith("barcode")):
            sid = dna
        var = tr.get("variant", "")
        key = (dna, var)
        seen_keys.add(key)
        if sid:
            seen_keys.add((sid, var))
        calls = caller_results.get(key) or (caller_results.get((sid, var), {}) if sid else {})

        dataset.append({
            "dna_number": dna,
            "sample_id": sid,
            "variant": var,
            "miseq": tr.get("miseq_af"),
            "sequel": tr.get("sequel_af"),
            "branford": calls.get("Branford"),
            "bcftools": calls.get("bcftools"),
            "rindels": calls.get("Rindels"),
            "vardict": calls.get("VarDictRust"),
        })

    # 2. Add caller results not in truth table (e.g. F359V, E255K, K375E, subclonal E499E, F359I)
    for (dna_or_sid, var), calls in caller_results.items():
        if not var or var.lower().startswith("snv") or "chr22" in var:
            continue
        key = (dna_or_sid, var)
        if key not in seen_keys:
            sid = sample_map.get(dna_or_sid, "") or (dna_or_sid if dna_or_sid.lower().startswith("bcrabl-") else "")
            dataset.append({
                "dna_number": dna_or_sid,
                "sample_id": sid,
                "variant": var,
                "miseq": None,
                "sequel": None,
                "branford": calls.get("Branford"),
                "bcftools": calls.get("bcftools"),
                "rindels": calls.get("Rindels"),
                "vardict": calls.get("VarDictRust"),
            })
            seen_keys.add(key)

    return dataset


def create_comparison_figure(
    dataset: List[Dict[str, Any]],
    out_path: Path,
    dpi: int = 300,
) -> None:
    """Create publication-ready multi-panel figure."""
    # Setup styling
    plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.edgecolor"] = "#475569"
    plt.rcParams["axes.linewidth"] = 0.9

    fig = plt.figure(figsize=(18.2, 5.8), dpi=dpi)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.15, 1.45, 0.95], wspace=0.28, left=0.06, right=0.97, top=0.90, bottom=0.12)

    ax_corr = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_bland = fig.add_subplot(gs[0, 2])

    # -------------------------------------------------------------
    # Panel A: Correlation & Linearity (Nanopore vs Orthogonal Truth)
    # -------------------------------------------------------------
    x_miseq_branford = []
    y_branford_for_mis = []

    # Gather paired points
    paired_pts = []
    seen_labels = set()

    for d in dataset:
        mis = d["miseq"]
        seq = d["sequel"]
        bf = d["branford"]
        bcf = d["bcftools"]
        rin = d["rindels"]
        vd = d["vardict"]

        if mis is not None and bf is not None:
            x_miseq_branford.append(mis)
            y_branford_for_mis.append(bf)

        # Plot all callers against MiSeq
        for caller_name, val, color, marker in [
            ("Branford", bf, CALLER_COLORS["Branford"], CALLER_MARKERS["Branford"]),
            ("bcftools", bcf, CALLER_COLORS["bcftools"], CALLER_MARKERS["bcftools"]),
            ("Rindels", rin, CALLER_COLORS["Rindels"], CALLER_MARKERS["Rindels"]),
            ("VarDictRust", vd, CALLER_COLORS["VarDictRust"], CALLER_MARKERS["VarDictRust"]),
        ]:
            if mis is not None and val is not None:
                lbl = f"ONT {caller_name} vs MiSeq"
                ax_corr.scatter(
                    mis, val, color=color, marker=marker, s=55, alpha=0.88,
                    edgecolors="#0F172A", linewidth=0.6,
                    label=lbl if lbl not in seen_labels else ""
                )
                seen_labels.add(lbl)
                paired_pts.append((mis, val, caller_name, "MiSeq"))
            if seq is not None and val is not None:
                lbl = f"ONT {caller_name} vs Sequel"
                ax_corr.scatter(
                    seq, val, facecolor="none", edgecolor=color, marker=marker, s=65,
                    linewidth=1.3, alpha=0.9,
                    label=lbl if lbl not in seen_labels else ""
                )
                seen_labels.add(lbl)
                paired_pts.append((seq, val, caller_name, "Sequel"))

    # Identity diagonal line
    ax_corr.plot([0, 105], [0, 105], linestyle="--", color="#94A3B8", linewidth=1.2, zorder=1, label="Identity (y = x)")

    # Linear regression fit for Branford vs MiSeq
    if len(x_miseq_branford) >= 2:
        coeffs = np.polyfit(x_miseq_branford, y_branford_for_mis, deg=1)
        poly_fn = np.poly1d(coeffs)
        x_line = np.linspace(0, 102, 100)
        ax_corr.plot(x_line, poly_fn(x_line), color="#1E3A8A", linewidth=1.6, alpha=0.95, zorder=2)

        # Compute R2 and Pearson r
        r_mat = np.corrcoef(x_miseq_branford, y_branford_for_mis)
        r_val = r_mat[0, 1]
        r2_val = r_val ** 2
        mae = np.mean(np.abs(np.array(x_miseq_branford) - np.array(y_branford_for_mis)))

        reg_text = (
            f"$\\mathbf{{Linear\\ Fit\\ (Branford\\ vs\\ MiSeq):}}$\n"
            f"$y = {coeffs[0]:.3f}x + {coeffs[1]:.2f}$\n"
            f"$R^2 = {r2_val:.4f}$  $(r = {r_val:.4f})$\n"
            f"Mean Abs Error $= {mae:.2f}\\%$\n"
            f"$n = {len(x_miseq_branford)}$ orthogonal sites"
        )
        ax_corr.text(
            0.04, 0.96, reg_text, transform=ax_corr.transAxes,
            fontsize=8.5, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#F8FAFC", edgecolor="#CBD5E1", alpha=0.95)
        )

    ax_corr.set_xlim(-2, 105)
    ax_corr.set_ylim(-2, 105)
    ax_corr.set_xlabel("Orthogonal Validation VAF (%) [MiSeq / Sequel]", fontsize=10, fontweight="bold", labelpad=6)
    ax_corr.set_ylabel("Oxford Nanopore VAF (%) [Nagini callers]", fontsize=10, fontweight="bold", labelpad=6)
    ax_corr.set_title("A. Cross-Platform Variant Allele Frequency Correlation", fontsize=11, fontweight="bold", loc="left", pad=10)
    ax_corr.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")

    # Inset for low-frequency region (0 - 3%)
    ax_inset = ax_corr.inset_axes([0.54, 0.08, 0.42, 0.38])
    for mis, val, caller, truth_src in paired_pts:
        if mis <= 3.5 and val <= 3.5:
            c = CALLER_COLORS[caller]
            m = CALLER_MARKERS[caller]
            if truth_src == "MiSeq":
                ax_inset.scatter(mis, val, color=c, marker=m, s=60, edgecolors="#0F172A", linewidth=0.6)
            else:
                ax_inset.scatter(mis, val, facecolor="none", edgecolor=c, marker=m, s=65, linewidth=1.2)

    ax_inset.plot([0, 3.5], [0, 3.5], linestyle="--", color="#94A3B8", linewidth=1.0)
    ax_inset.set_xlim(0, 3.0)
    ax_inset.set_ylim(0, 3.0)
    ax_inset.set_title("Subclonal Range (0–3%)", fontsize=7.5, fontweight="bold", pad=2)
    ax_inset.tick_params(labelsize=7)
    ax_inset.grid(True, linestyle=":", alpha=0.4)
    # Highlight E459K point annotation in inset
    ax_inset.annotate("E459K (1.7% / 2.0%)", xy=(1.7, 2.01), xytext=(0.2, 2.4),
                      arrowprops=dict(arrowstyle="->", color="#0F172A", lw=0.7),
                      fontsize=7, fontweight="bold", color="#1E3A8A")

    # -------------------------------------------------------------
    # Panel B: Grouped VAF Concordance across Dynamic Range
    # -------------------------------------------------------------
    # Select variants for bar chart
    bar_variants = [
        d for d in dataset
        if d.get("branford") is not None or d.get("miseq") is not None
    ]
    # Filter out entries without variant name
    bar_variants = [d for d in bar_variants if d.get("variant") and not d["variant"].lower().startswith("ins")]

    # Sort order: by VAF ascending
    def get_sort_key(d):
        return d.get("miseq") or d.get("sequel") or d.get("branford") or 0.0

    bar_variants.sort(key=get_sort_key)

    labels = []
    miseq_vals = []
    sequel_vals = []
    branford_vals = []
    bcftools_vals = []
    rindels_vals = []
    vardict_vals = []

    for d in bar_variants:
        bc_tag = d['sample_id'].replace("BCRABL-", "") if d['sample_id'] else d['dna_number']
        lbl = f"{d['variant']}\n({bc_tag})"
        labels.append(lbl)
        miseq_vals.append(d.get("miseq") or 0.0)
        sequel_vals.append(d.get("sequel") or 0.0)
        branford_vals.append(d.get("branford") or 0.0)
        bcftools_vals.append(d.get("bcftools") or 0.0)
        rindels_vals.append(d.get("rindels") or 0.0)
        vardict_vals.append(d.get("vardict") or 0.0)

    x = np.arange(len(labels))
    width = 0.13

    ax_bar.bar(x - 2.5 * width, miseq_vals, width, label="Illumina MiSeq", color=CALLER_COLORS["MiSeq"], edgecolor="#1E293B", linewidth=0.5)
    ax_bar.bar(x - 1.5 * width, sequel_vals, width, label="PacBio Sequel", color=CALLER_COLORS["Sequel"], edgecolor="#0369A1", linewidth=0.5)
    ax_bar.bar(x - 0.5 * width, branford_vals, width, label="ONT Branford", color=CALLER_COLORS["Branford"], edgecolor="#172554", linewidth=0.5)
    ax_bar.bar(x + 0.5 * width, bcftools_vals, width, label="ONT bcftools", color=CALLER_COLORS["bcftools"], edgecolor="#78350F", linewidth=0.5)
    ax_bar.bar(x + 1.5 * width, rindels_vals, width, label="ONT Rindels", color=CALLER_COLORS["Rindels"], edgecolor="#4C1D95", linewidth=0.5)
    ax_bar.bar(x + 2.5 * width, vardict_vals, width, label="ONT VarDictRust", color=CALLER_COLORS["VarDictRust"], edgecolor="#7F1D1D", linewidth=0.5)

    ax_bar.set_ylabel("Variant Allele Frequency (%)", fontsize=10, fontweight="bold", labelpad=6)
    ax_bar.set_title("B. Allele Frequency Concordance across Callers & Dynamic Range", fontsize=11, fontweight="bold", loc="left", pad=10)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels, fontsize=7.8, rotation=0, ha="center")
    ax_bar.set_ylim(0, 118)
    ax_bar.grid(True, axis="y", linestyle=":", alpha=0.5, color="#CBD5E1")
    ax_bar.legend(loc="upper left", frameon=True, fontsize=8.0, ncol=3,
                  facecolor="#F8FAFC", edgecolor="#CBD5E1", framealpha=0.95)

    # -------------------------------------------------------------
    # Panel C: Residuals / Bland-Altman (Delta VAF % = ONT - MiSeq)
    # -------------------------------------------------------------
    deltas_branford = []
    means_branford = []
    deltas_bcftools = []
    means_bcftools = []
    deltas_rindels = []
    means_rindels = []
    deltas_vardict = []
    means_vardict = []

    for d in dataset:
        mis = d["miseq"]
        if mis is None:
            continue
        if d["branford"] is not None:
            deltas_branford.append(d["branford"] - mis)
            means_branford.append((d["branford"] + mis) / 2.0)
        if d["bcftools"] is not None:
            deltas_bcftools.append(d["bcftools"] - mis)
            means_bcftools.append((d["bcftools"] + mis) / 2.0)
        if d["rindels"] is not None:
            deltas_rindels.append(d["rindels"] - mis)
            means_rindels.append((d["rindels"] + mis) / 2.0)
        if d["vardict"] is not None:
            deltas_vardict.append(d["vardict"] - mis)
            means_vardict.append((d["vardict"] + mis) / 2.0)

    ax_bland.scatter(means_branford, deltas_branford, color=CALLER_COLORS["Branford"], marker=CALLER_MARKERS["Branford"], s=55, label="Branford - MiSeq", alpha=0.9, edgecolors="#0F172A", linewidth=0.6)
    ax_bland.scatter(means_bcftools, deltas_bcftools, color=CALLER_COLORS["bcftools"], marker=CALLER_MARKERS["bcftools"], s=55, label="bcftools - MiSeq", alpha=0.9, edgecolors="#0F172A", linewidth=0.6)
    ax_bland.scatter(means_rindels, deltas_rindels, color=CALLER_COLORS["Rindels"], marker=CALLER_MARKERS["Rindels"], s=55, label="Rindels - MiSeq", alpha=0.9, edgecolors="#0F172A", linewidth=0.6)
    ax_bland.scatter(means_vardict, deltas_vardict, color=CALLER_COLORS["VarDictRust"], marker=CALLER_MARKERS["VarDictRust"], s=55, label="VarDict - MiSeq", alpha=0.9, edgecolors="#0F172A", linewidth=0.6)

    # Reference line at zero bias
    ax_bland.axhline(0, color="#0F172A", linestyle="-", linewidth=1.1, alpha=0.7)
    ax_bland.axhline(1.0, color="#94A3B8", linestyle="--", linewidth=0.9, label="±1.0% tolerance")
    ax_bland.axhline(-1.0, color="#94A3B8", linestyle="--", linewidth=0.9)

    # Mean difference and bounds for Branford
    if deltas_branford:
        mean_diff = np.mean(deltas_branford)
        sd_diff = np.std(deltas_branford)
        ax_bland.axhline(mean_diff, color="#1E3A8A", linestyle=":", linewidth=1.4, alpha=0.85)

        stat_box = (
            f"$\\mathbf{{Mean\\ Bias\\ (Branford):}}$\n"
            f"${mean_diff:+.2f}\\% \\pm {sd_diff:.2f}\\%$\n"
            f"Max $|\\Delta| = {np.max(np.abs(deltas_branford)):.2f}\\%$"
        )
        ax_bland.text(
            0.05, 0.15, stat_box, transform=ax_bland.transAxes,
            fontsize=8.2, verticalalignment="bottom",
            bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8FAFC", edgecolor="#CBD5E1", alpha=0.95)
        )

    ax_bland.set_xlim(-2, 105)
    ax_bland.set_ylim(-4.5, 4.5)
    ax_bland.set_xlabel("Mean VAF (%) [(ONT + MiSeq) / 2]", fontsize=10, fontweight="bold", labelpad=6)
    ax_bland.set_ylabel("$\\Delta$ VAF (%) [Nanopore − MiSeq]", fontsize=10, fontweight="bold", labelpad=6)
    ax_bland.set_title("C. VAF Residuals & Precision", fontsize=11, fontweight="bold", loc="left", pad=10)
    ax_bland.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
    ax_bland.legend(loc="upper right", frameon=True, fontsize=7.8,
                    facecolor="#F8FAFC", edgecolor="#CBD5E1", framealpha=0.95)

    # Clean styling on all axes
    for ax in (ax_corr, ax_bar, ax_bland):
        ax.tick_params(axis="both", which="major", labelsize=8.5)
        ax.set_facecolor("#FFFFFF")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[Nagini] Successfully saved comparison figure to: {out_path}")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plot publication-ready comparison figure across MiSeq, Sequel, and Nagini callers."
    )
    parser.add_argument(
        "--truth-table", "-t", type=Path, default=Path("references/table6_precision.txt"),
        help="Path to Table 6 text file (default: references/table6_precision.txt)"
    )
    parser.add_argument(
        "--variants-tsv", "-v", type=Path, default=Path("references/bcr_abl_caller_results.tsv"),
        help="Path to caller results TSV (default: references/bcr_abl_caller_results.tsv)"
    )
    parser.add_argument(
        "--bed", "-b", type=Path, default=Path("workflow/resources/bcr_abl/branford_positions_hg19.bed"),
        help="Path to Branford positions BED (default: workflow/resources/bcr_abl/branford_positions_hg19.bed)"
    )
    parser.add_argument(
        "--sample-map", "-m", type=Path, default=Path("references/bcr_abl_sample_map.tsv"),
        help="Path to sample mapping TSV (default: references/bcr_abl_sample_map.tsv)"
    )
    parser.add_argument(
        "--output", "-o", type=Path, default=Path("docs/poster/fig3_technology_comparison.png"),
        help="Output plot image path (PNG or SVG, default: docs/poster/fig3_technology_comparison.png)"
    )
    parser.add_argument(
        "--dpi", type=int, default=300,
        help="DPI for raster outputs (default: 300)"
    )

    args = parser.parse_args(argv)

    # 1. Parse Branford BED
    coords_to_aa = {}
    if args.bed and args.bed.is_file():
        coords_to_aa = load_branford_bed(args.bed)

    # 2. Parse truth table
    truth_records = load_truth_table(args.truth_table)
    if not truth_records:
        sys.stderr.write(f"Warning: no truth records loaded from {args.truth_table}\n")

    # 3. Parse caller results
    caller_results = load_caller_results(args.variants_tsv, coords_to_aa=coords_to_aa)
    if not caller_results:
        sys.stderr.write(f"Warning: no caller results loaded from {args.variants_tsv}\n")

    # 4. Parse sample map if available
    s_map = {}
    if args.sample_map and args.sample_map.is_file():
        with open(args.sample_map, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for r in reader:
                dna = r.get("dna_number", "").strip()
                sid = r.get("sample_id", "").strip()
                if dna and sid:
                    s_map[dna] = sid

    # 5. Assemble dataset
    dataset = build_plot_dataset(truth_records, caller_results, sample_map=s_map)

    # 6. Generate figure
    create_comparison_figure(dataset, args.output, dpi=args.dpi)

    # Also generate vector SVG alongside PNG if saving to docs/poster
    if args.output.suffix.lower() == ".png":
        svg_path = args.output.with_suffix(".svg")
        try:
            create_comparison_figure(dataset, svg_path, dpi=args.dpi)
        except Exception as exc:
            sys.stderr.write(f"Notice: could not save SVG copy: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
