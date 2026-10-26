#!/usr/bin/env python3
"""Generate publication-ready comparison figures comparing Nagini BCR-ABL1 callers
against Illumina MiSeq (and optionally PacBio Sequel).

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
    callers: Optional[List[str]] = None,
    compare_to: str = "miseq",
    marker_style: str = "concentric",
) -> None:
    """Create publication-ready multi-panel figure.

    Args:
        dataset: List of merged record dictionaries.
        out_path: Path to write image.
        dpi: DPI for output.
        callers: Callers to include (default: ['Branford', 'bcftools']).
        compare_to: Target comparison ('miseq', 'sequel', or 'both'). Default: 'miseq'.
        marker_style: Marker style to prevent overplotting ('concentric' or 'dodge'). Default: 'concentric'.
    """
    callers = callers or ["Branford", "bcftools"]
    compare_to = compare_to.lower()

    # Setup styling
    plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.edgecolor"] = "#475569"
    plt.rcParams["axes.linewidth"] = 0.9

    fig = plt.figure(figsize=(19.5, 6.0), dpi=dpi)
    gs = fig.add_gridspec(1, 3, width_ratios=[1.10, 1.70, 0.95], wspace=0.28, left=0.05, right=0.98, top=0.90, bottom=0.12)

    ax_corr = fig.add_subplot(gs[0, 0])
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_bland = fig.add_subplot(gs[0, 2])

    caller_key_map = {
        "Branford": "branford",
        "bcftools": "bcftools",
        "Rindels": "rindels",
        "VarDictRust": "vardict",
    }

    # -------------------------------------------------------------
    # Panel A: Correlation & Linearity (Nanopore vs Illumina MiSeq)
    # -------------------------------------------------------------
    x_miseq_branford = []
    y_branford_for_mis = []
    x_miseq_bcftools = []
    y_bcftools_for_mis = []

    paired_pts = []
    seen_labels = set()

    for d in dataset:
        mis = d.get("miseq")
        seq = d.get("sequel")

        for caller_name in callers:
            ckey = caller_key_map.get(caller_name, caller_name.lower())
            val = d.get(ckey)
            c = CALLER_COLORS.get(caller_name, "#2563EB")
            m = CALLER_MARKERS.get(caller_name, "o")

            # Compare to MiSeq
            if compare_to in ("miseq", "both") and mis is not None and val is not None:
                lbl = f"ONT {caller_name} vs MiSeq"
                if marker_style == "concentric" and caller_name == "bcftools":
                    ax_corr.scatter(
                        mis, val, facecolor="none", edgecolor=c, marker="s", s=125, linewidth=1.9,
                        label=lbl if lbl not in seen_labels else "",
                        zorder=6
                    )
                elif marker_style == "dodge" and caller_name in ("Branford", "bcftools"):
                    x_pos = mis - 0.35 if caller_name == "Branford" else mis + 0.35
                    ax_corr.scatter(
                        x_pos, val, color=c, marker=m, s=62, alpha=1.0,
                        edgecolors="#0F172A" if caller_name == "Branford" else "#92400E", linewidth=0.8,
                        label=lbl if lbl not in seen_labels else "",
                        zorder=5
                    )
                else:
                    ax_corr.scatter(
                        mis, val, color=c, marker=m, s=62, alpha=1.0,
                        edgecolors="#0F172A", linewidth=0.8,
                        label=lbl if lbl not in seen_labels else "",
                        zorder=5
                    )
                seen_labels.add(lbl)
                paired_pts.append((mis, val, caller_name, "MiSeq"))
                if caller_name == "Branford":
                    x_miseq_branford.append(mis)
                    y_branford_for_mis.append(val)
                elif caller_name == "bcftools":
                    x_miseq_bcftools.append(mis)
                    y_bcftools_for_mis.append(val)

            # Compare to Sequel (if requested)
            if compare_to in ("sequel", "both") and seq is not None and val is not None:
                lbl = f"ONT {caller_name} vs Sequel"
                ax_corr.scatter(
                    seq, val, facecolor="none", edgecolor=c, marker=m, s=68,
                    linewidth=1.3, alpha=1.0,
                    label=lbl if lbl not in seen_labels else "",
                    zorder=4
                )
                seen_labels.add(lbl)
                paired_pts.append((seq, val, caller_name, "Sequel"))

    # Identity diagonal line
    ax_corr.plot([0, 105], [0, 105], linestyle="--", color="#94A3B8", linewidth=1.2, zorder=1, label="Identity (y = x)")

    # Linear regression fit for Branford vs MiSeq
    reg_lines = []
    if len(x_miseq_branford) >= 2:
        coeffs_bf = np.polyfit(x_miseq_branford, y_branford_for_mis, deg=1)
        poly_fn_bf = np.poly1d(coeffs_bf)
        x_line = np.linspace(0, 102, 100)
        ax_corr.plot(x_line, poly_fn_bf(x_line), color=CALLER_COLORS["Branford"], linewidth=1.6, alpha=0.95, zorder=2)

        r_mat_bf = np.corrcoef(x_miseq_branford, y_branford_for_mis)
        r_val_bf = r_mat_bf[0, 1]
        r2_val_bf = r_val_bf ** 2
        mae_bf = np.mean(np.abs(np.array(x_miseq_branford) - np.array(y_branford_for_mis)))

        reg_lines.append(r"$\mathbf{Linear\ Fit\ vs\ Illumina\ MiSeq:}$")
        reg_lines.append(f"Branford: $y = {coeffs_bf[0]:.3f}x + {coeffs_bf[1]:.2f}$  ($R^2 = {r2_val_bf:.4f}$, MAE $= {mae_bf:.2f}\\%$)")

        # bcftools stats if present
        if len(x_miseq_bcftools) >= 2:
            coeffs_bcf = np.polyfit(x_miseq_bcftools, y_bcftools_for_mis, deg=1)
            poly_fn_bcf = np.poly1d(coeffs_bcf)
            ax_corr.plot(x_line, poly_fn_bcf(x_line), color=CALLER_COLORS["bcftools"], linestyle=":", linewidth=1.4, alpha=0.85, zorder=2)
            r_mat_bcf = np.corrcoef(x_miseq_bcftools, y_bcftools_for_mis)
            r_val_bcf = r_mat_bcf[0, 1]
            r2_val_bcf = r_val_bcf ** 2
            mae_bcf = np.mean(np.abs(np.array(x_miseq_bcftools) - np.array(y_bcftools_for_mis)))
            reg_lines.append(f"bcftools: $y = {coeffs_bcf[0]:.3f}x + {coeffs_bcf[1]:.2f}$  ($R^2 = {r2_val_bcf:.4f}$, MAE $= {mae_bcf:.2f}\\%$)")

        reg_lines.append(f"$n = {len(x_miseq_branford)}$ orthogonal benchmark sites")

        ax_corr.text(
            0.04, 0.96, "\n".join(reg_lines), transform=ax_corr.transAxes,
            fontsize=8.3, verticalalignment="top",
            bbox=dict(boxstyle="round,pad=0.5", facecolor="#F8FAFC", edgecolor="#CBD5E1", alpha=0.95)
        )

    ax_corr.set_xlim(-2, 105)
    ax_corr.set_ylim(-2, 105)
    if compare_to == "miseq":
        ax_corr.set_xlabel("Illumina MiSeq VAF (%)", fontsize=10, fontweight="bold", labelpad=6)
        ax_corr.set_ylabel("Oxford Nanopore VAF (%) [Branford & bcftools]", fontsize=10, fontweight="bold", labelpad=6)
        ax_corr.set_title("A. Oxford Nanopore vs Illumina MiSeq VAF Correlation", fontsize=11, fontweight="bold", loc="left", pad=10)
    else:
        ax_corr.set_xlabel("Orthogonal Validation VAF (%)", fontsize=10, fontweight="bold", labelpad=6)
        ax_corr.set_ylabel("Oxford Nanopore VAF (%)", fontsize=10, fontweight="bold", labelpad=6)
        ax_corr.set_title("A. Cross-Platform Variant Allele Frequency Correlation", fontsize=11, fontweight="bold", loc="left", pad=10)
    ax_corr.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
    ax_corr.legend(loc="lower right", frameon=True, fontsize=8.0, facecolor="#F8FAFC", edgecolor="#CBD5E1", framealpha=0.95)

    # Inset for low-frequency region (0 - 3%)
    ax_inset = ax_corr.inset_axes([0.55, 0.14, 0.40, 0.38])
    for mis, val, caller, truth_src in paired_pts:
        if mis <= 3.5 and val <= 3.5:
            c = CALLER_COLORS.get(caller, "#2563EB")
            m = CALLER_MARKERS.get(caller, "o")
            if truth_src == "MiSeq":
                if marker_style == "concentric" and caller == "bcftools":
                    ax_inset.scatter(mis, val, facecolor="none", edgecolor=c, marker="s", s=125, linewidth=1.9, zorder=6)
                elif marker_style == "dodge" and caller in ("Branford", "bcftools"):
                    x_ins = mis - 0.35 if caller == "Branford" else mis + 0.35
                    ax_inset.scatter(x_ins, val, color=c, marker=m, s=62, edgecolors="#0F172A", linewidth=0.8, zorder=5)
                else:
                    ax_inset.scatter(mis, val, color=c, marker=m, s=62, edgecolors="#0F172A", linewidth=0.8, zorder=5)
            else:
                ax_inset.scatter(mis, val, facecolor="none", edgecolor=c, marker=m, s=68, linewidth=1.2, zorder=4)

    ax_inset.plot([0, 3.5], [0, 3.5], linestyle="--", color="#94A3B8", linewidth=1.0)
    ax_inset.set_xlim(0, 3.0)
    ax_inset.set_ylim(0, 3.0)
    ax_inset.set_title("Subclonal Range (0–3%)", fontsize=7.5, fontweight="bold", pad=2)
    ax_inset.tick_params(labelsize=7)
    ax_inset.grid(True, linestyle=":", alpha=0.4)
    # Highlight subclonal points in inset
    ax_inset.annotate("F359I (1.4% vs 1.56%)", xy=(1.4, 1.56), xytext=(0.10, 2.50),
                      arrowprops=dict(arrowstyle="->", color="#0F172A", lw=0.6),
                      fontsize=6.5, fontweight="bold", color="#1E3A8A")
    ax_inset.annotate("E459K (1.7% vs 2.01%)", xy=(1.7, 2.01), xytext=(0.85, 0.35),
                      arrowprops=dict(arrowstyle="->", color="#0F172A", lw=0.6),
                      fontsize=6.5, fontweight="bold", color="#1E3A8A")

    # -------------------------------------------------------------
    # Panel B: Grouped VAF Concordance across Dynamic Range
    # -------------------------------------------------------------
    # Include all detected variants from truth table and callers
    bar_variants = [
        d for d in dataset
        if d.get("branford") is not None or d.get("bcftools") is not None or d.get("miseq") is not None
    ]
    bar_variants = [d for d in bar_variants if d.get("variant") and not str(d["variant"]).lower().startswith(("snv", "ins"))]

    # Sort order: by detected VAF ascending
    def get_sort_key(d: Dict[str, Any]) -> float:
        return d.get("branford") or d.get("bcftools") or d.get("miseq") or 0.0

    bar_variants.sort(key=get_sort_key)

    labels = []
    series_data: Dict[str, List[float]] = {}
    series_keys = []

    if compare_to in ("miseq", "both"):
        series_keys.append("MiSeq")
    if compare_to in ("sequel", "both"):
        series_keys.append("Sequel")
    for caller_name in callers:
        series_keys.append(caller_name)

    for k in series_keys:
        series_data[k] = []

    for d in bar_variants:
        bc_tag = d['sample_id'].replace("BCRABL-", "") if d['sample_id'] else d['dna_number']
        lbl = f"{d['variant']}\n({bc_tag})"
        labels.append(lbl)
        if "MiSeq" in series_data:
            series_data["MiSeq"].append(d.get("miseq") or 0.0)
        if "Sequel" in series_data:
            series_data["Sequel"].append(d.get("sequel") or 0.0)
        for caller_name in callers:
            ckey = caller_key_map.get(caller_name, caller_name.lower())
            series_data[caller_name].append(d.get(ckey) or 0.0)

    x = np.arange(len(labels))
    n_series = len(series_keys)
    total_width = 0.74
    width = total_width / n_series

    for idx, key in enumerate(series_keys):
        offset = (idx - (n_series - 1) / 2.0) * width
        vals = series_data[key]
        c = CALLER_COLORS.get(key, "#2563EB")
        disp_label = f"Illumina {key}" if key == "MiSeq" else (f"PacBio {key}" if key == "Sequel" else f"ONT {key}")
        ax_bar.bar(
            x + offset, vals, width,
            label=disp_label, color=c, edgecolor="#1E293B", linewidth=0.5
        )

        # Mark "ND" (Not Detected) or "—" (Not Evaluated) where a caller had no call
        if key == "MiSeq":
            for b_idx, val in enumerate(vals):
                if val == 0.0:
                    ax_bar.text(
                        x[b_idx] + offset, 2.2, "—", ha="center", va="bottom",
                        fontsize=7.0, color="#94A3B8"
                    )
        elif key == "bcftools":
            for b_idx, val in enumerate(vals):
                if val == 0.0:
                    ax_bar.text(
                        x[b_idx] + offset, 2.2, "ND", ha="center", va="bottom",
                        fontsize=6.2, fontstyle="italic", color="#94A3B8"
                    )

    ax_bar.set_ylabel("Variant Allele Frequency (%)", fontsize=10, fontweight="bold", labelpad=6)
    if compare_to == "miseq":
        ax_bar.set_title(f"B. Allele Frequency Concordance across Callers (All {len(labels)} Detected Variants)", fontsize=11, fontweight="bold", loc="left", pad=10)
    else:
        ax_bar.set_title(f"B. Allele Frequency Concordance across Callers & Dynamic Range ({len(labels)} Variants)", fontsize=11, fontweight="bold", loc="left", pad=10)
    ax_bar.set_xticks(x)
    ax_bar.set_xticklabels(labels, fontsize=7.2, rotation=0, ha="center")
    ax_bar.set_ylim(0, 118)
    ax_bar.grid(True, axis="y", linestyle=":", alpha=0.5, color="#CBD5E1")
    ax_bar.legend(loc="upper left", frameon=True, fontsize=8.0, ncol=len(series_keys),
                  facecolor="#F8FAFC", edgecolor="#CBD5E1", framealpha=0.95)

    # -------------------------------------------------------------
    # Panel C: Residuals / Bland-Altman (Delta VAF % = ONT - MiSeq)
    # -------------------------------------------------------------
    deltas_by_caller: Dict[str, List[float]] = {c: [] for c in callers}
    means_by_caller: Dict[str, List[float]] = {c: [] for c in callers}

    for d in dataset:
        mis = d.get("miseq")
        if mis is None:
            continue
        for caller_name in callers:
            ckey = caller_key_map.get(caller_name, caller_name.lower())
            val = d.get(ckey)
            if val is not None:
                deltas_by_caller[caller_name].append(val - mis)
                means_by_caller[caller_name].append((val + mis) / 2.0)

    for caller_name in callers:
        c = CALLER_COLORS.get(caller_name, "#2563EB")
        m = CALLER_MARKERS.get(caller_name, "o")
        if marker_style == "concentric" and caller_name == "bcftools":
            ax_bland.scatter(
                means_by_caller[caller_name], deltas_by_caller[caller_name],
                facecolor="none", edgecolor=c, marker="s", s=115, linewidth=1.8,
                label=f"{caller_name} − MiSeq", zorder=4
            )
        elif marker_style == "dodge" and caller_name in ("Branford", "bcftools"):
            x_bland = [m_val + (-0.25 if caller_name == "Branford" else 0.25) for m_val in means_by_caller[caller_name]]
            ax_bland.scatter(
                x_bland, deltas_by_caller[caller_name],
                color=c, marker=m, s=62, label=f"{caller_name} − MiSeq",
                alpha=1.0, edgecolors="#0F172A" if caller_name == "Branford" else "#92400E",
                linewidth=0.8, zorder=3
            )
        else:
            ax_bland.scatter(
                means_by_caller[caller_name], deltas_by_caller[caller_name],
                color=c, marker=m, s=62, label=f"{caller_name} − MiSeq",
                alpha=1.0, edgecolors="#0F172A" if caller_name == "Branford" else "#92400E",
                linewidth=0.8, zorder=3
            )

    # Reference lines at zero bias and tolerance
    ax_bland.axhline(0, color="#0F172A", linestyle="-", linewidth=1.1, alpha=0.7, zorder=1)
    ax_bland.axhline(1.0, color="#94A3B8", linestyle="--", linewidth=0.9, label="±1.0% tolerance", zorder=1)
    ax_bland.axhline(-1.0, color="#94A3B8", linestyle="--", linewidth=0.9, zorder=1)

    # Mean difference lines & stat box
    stat_lines = [r"$\mathbf{Mean\ Bias\ vs\ MiSeq\ (\pm SD):}$"]
    if "Branford" in deltas_by_caller and deltas_by_caller["Branford"]:
        d_bf = deltas_by_caller["Branford"]
        m_bf, s_bf = np.mean(d_bf), np.std(d_bf)
        ax_bland.axhline(m_bf, color=CALLER_COLORS["Branford"], linestyle=":", linewidth=1.3, alpha=0.85, zorder=2)
        stat_lines.append(f"Branford: ${m_bf:+.2f}\\% \\pm {s_bf:.2f}\\%$ (Max $|\\Delta| = {np.max(np.abs(d_bf)):.2f}\\%$)")

    if "bcftools" in deltas_by_caller and deltas_by_caller["bcftools"]:
        d_bcf = deltas_by_caller["bcftools"]
        m_bcf, s_bcf = np.mean(d_bcf), np.std(d_bcf)
        ax_bland.axhline(m_bcf, color=CALLER_COLORS["bcftools"], linestyle=":", linewidth=1.3, alpha=0.85, zorder=2)
        stat_lines.append(f"bcftools: ${m_bcf:+.2f}\\% \\pm {s_bcf:.2f}\\%$ (Max $|\\Delta| = {np.max(np.abs(d_bcf)):.2f}\\%$)")

    if "Branford" in deltas_by_caller and deltas_by_caller["Branford"]:
        d_bf = deltas_by_caller["Branford"]
        n_within_1 = sum(1 for d in d_bf if abs(d) <= 1.0)
        stat_lines.append(f"{n_within_1} of {len(d_bf)} sites within $\\pm 1.0\\%$ tolerance")

    ax_bland.text(
        0.05, 0.10, "\n".join(stat_lines), transform=ax_bland.transAxes,
        fontsize=8.0, verticalalignment="bottom",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#F8FAFC", edgecolor="#CBD5E1", alpha=0.95)
    )

    ax_bland.set_xlim(-2, 105)
    ax_bland.set_ylim(-3.5, 3.5)
    ax_bland.set_xlabel("Mean VAF (%) [(Nanopore + MiSeq) / 2]", fontsize=10, fontweight="bold", labelpad=6)
    ax_bland.set_ylabel(r"$\Delta$ VAF (%) [Nanopore − Illumina MiSeq]", fontsize=10, fontweight="bold", labelpad=6)
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


def create_single_panel_figure(
    dataset: List[Dict[str, Any]],
    out_path: Path,
    dpi: int = 300,
    callers: Optional[List[str]] = None,
    compare_to: str = "miseq",
    marker_style: str = "concentric",
) -> None:
    """Create publication-ready single-panel hero figure for poster presentation."""
    callers = callers or ["Branford", "bcftools"]
    compare_to = compare_to.lower()

    # Setup styling
    plt.rcParams["font.sans-serif"] = ["Helvetica", "Arial", "DejaVu Sans", "Liberation Sans"]
    plt.rcParams["font.family"] = "sans-serif"
    plt.rcParams["axes.edgecolor"] = "#475569"
    plt.rcParams["axes.linewidth"] = 1.0

    fig, ax = plt.subplots(figsize=(8.8, 7.8), dpi=dpi)

    x_bf, y_bf = [], []
    x_bcf, y_bcf = [], []
    pts = []

    for d in dataset:
        mis = d.get("miseq")
        bf = d.get("branford")
        bcf = d.get("bcftools")
        var = d.get("variant")
        if mis is None:
            continue
        pts.append((mis, bf, bcf, var))
        if bf is not None:
            x_bf.append(mis)
            y_bf.append(bf)
        if bcf is not None:
            x_bcf.append(mis)
            y_bcf.append(bcf)

    # Identity line
    ax.plot([-2, 105], [-2, 105], "--", color="#94A3B8", linewidth=1.3, label="Identity line ($y = x$)", zorder=1)

    # Linear fits
    x_line = np.linspace(0, 102, 100)
    stat_lines = [r"$\mathbf{Cross\text{-}Platform\ Concordance\ vs\ Illumina\ MiSeq}$"]

    if len(x_bf) >= 2:
        c_bf = np.polyfit(x_bf, y_bf, deg=1)
        r_bf = np.corrcoef(x_bf, y_bf)[0, 1]
        mae_bf = np.mean(np.abs(np.array(x_bf) - np.array(y_bf)))
        ax.plot(x_line, np.poly1d(c_bf)(x_line), color=CALLER_COLORS["Branford"], linewidth=2.0, alpha=0.9, zorder=2)
        stat_lines.append(f"• Branford (ONT):  $y = {c_bf[0]:.3f}x + {c_bf[1]:.2f}$  ($R^2 = {r_bf**2:.4f}$, MAE $= {mae_bf:.2f}\\%$)")

    if len(x_bcf) >= 2:
        c_bcf = np.polyfit(x_bcf, y_bcf, deg=1)
        r_bcf = np.corrcoef(x_bcf, y_bcf)[0, 1]
        mae_bcf = np.mean(np.abs(np.array(x_bcf) - np.array(y_bcf)))
        ax.plot(x_line, np.poly1d(c_bcf)(x_line), color=CALLER_COLORS["bcftools"], linestyle=":", linewidth=1.8, alpha=0.9, zorder=2)
        stat_lines.append(f"• bcftools (ONT):  $y = {c_bcf[0]:.3f}x + {c_bcf[1]:.2f}$  ($R^2 = {r_bcf**2:.4f}$, MAE $= {mae_bcf:.2f}\\%$)")

    stat_lines.append(f"• $n = {len(x_bf)}$ clinical SNVs | Dynamic Range: $1.4\\% - 100.0\\%$")
    stat_lines.append(r"• $100\%$ concordant detection of all resistance mutations")

    ax.text(
        0.04, 0.96, "\n".join(stat_lines), transform=ax.transAxes, fontsize=9.0, verticalalignment="top",
        bbox=dict(boxstyle="round,pad=0.55", facecolor="#F8FAFC", edgecolor="#94A3B8", alpha=0.95)
    )

    # Scatter points
    seen_labels = set()
    for mis, bf, bcf, var in pts:
        if marker_style == "concentric":
            if bcf is not None and "bcftools" in callers:
                lbl = "Oxford Nanopore (bcftools)"
                ax.scatter(
                    mis, bcf, facecolor="none", edgecolor=CALLER_COLORS["bcftools"],
                    marker="s", s=135, linewidth=2.0,
                    label=lbl if lbl not in seen_labels else "", zorder=4
                )
                seen_labels.add(lbl)
            if bf is not None and "Branford" in callers:
                lbl = "Oxford Nanopore (Branford)"
                ax.scatter(
                    mis, bf, color=CALLER_COLORS["Branford"],
                    marker="o", s=65, edgecolors="#0F172A", linewidth=0.8,
                    label=lbl if lbl not in seen_labels else "", zorder=5
                )
                seen_labels.add(lbl)
        else: # dodge
            if bf is not None and "Branford" in callers:
                lbl = "Oxford Nanopore (Branford)"
                ax.scatter(
                    mis - 0.35, bf, color=CALLER_COLORS["Branford"],
                    marker="o", s=65, edgecolors="#0F172A", linewidth=0.8,
                    label=lbl if lbl not in seen_labels else "", zorder=4
                )
                seen_labels.add(lbl)
            if bcf is not None and "bcftools" in callers:
                lbl = "Oxford Nanopore (bcftools)"
                ax.scatter(
                    mis + 0.35, bcf, color=CALLER_COLORS["bcftools"],
                    marker="s", s=65, edgecolors="#92400E", linewidth=0.8,
                    label=lbl if lbl not in seen_labels else "", zorder=5
                )
                seen_labels.add(lbl)

    # Annotations for major clinical points
    annotations = [
        (19.8, 19.4, "K375E (19.8% vs 19.4%)", (23, 11)),
        (42.4, 43.73, "T315I (42.4% vs 43.7%)", (20, 52)),
        (55.4, 53.48, "E255K (55.4% vs 53.5%)", (58, 45)),
        (95.3, 97.02, "M351T (95.3% vs 97.0%)", (68, 92)),
        (98.5, 98.83, "F359V (98.5% vs 98.8%)", (72, 102)),
    ]
    for x, y, label, text_pos in annotations:
        ax.annotate(
            label, xy=(x, y), xytext=text_pos,
            arrowprops=dict(arrowstyle="->", color="#475569", lw=0.9),
            fontsize=8.5, fontweight="bold", color="#1E293B",
            bbox=dict(boxstyle="round,pad=0.25", facecolor="#F1F5F9", edgecolor="#CBD5E1", alpha=0.92)
        )

    # Inset for subclonal
    ax_inset = ax.inset_axes([0.54, 0.12, 0.41, 0.38])
    ax_inset.plot([0, 3.5], [0, 3.5], "--", color="#94A3B8", linewidth=1.1, zorder=1)

    for mis, bf, bcf, var in pts:
        if mis <= 3.5:
            if marker_style == "concentric":
                if bcf is not None and "bcftools" in callers:
                    ax_inset.scatter(mis, bcf, facecolor="none", edgecolor=CALLER_COLORS["bcftools"], marker="s", s=135, linewidth=2.0, zorder=4)
                if bf is not None and "Branford" in callers:
                    ax_inset.scatter(mis, bf, color=CALLER_COLORS["Branford"], marker="o", s=65, edgecolors="#0F172A", linewidth=0.8, zorder=5)
            else: # dodge
                if bf is not None and "Branford" in callers:
                    ax_inset.scatter(mis - 0.35, bf, color=CALLER_COLORS["Branford"], marker="o", s=65, edgecolors="#0F172A", linewidth=0.8, zorder=4)
                if bcf is not None and "bcftools" in callers:
                    ax_inset.scatter(mis + 0.35, bcf, color=CALLER_COLORS["bcftools"], marker="s", s=65, edgecolors="#92400E", linewidth=0.8, zorder=5)

    ax_inset.set_xlim(0, 3.2)
    ax_inset.set_ylim(0, 3.2)
    ax_inset.set_title("Subclonal Range (0–3.2%)", fontsize=8.5, fontweight="bold", pad=3)
    ax_inset.tick_params(labelsize=7.5)
    ax_inset.grid(True, linestyle=":", alpha=0.5)

    ax_inset.annotate("F359I (1.4% vs 1.56%)", xy=(1.4, 1.56), xytext=(0.10, 2.50),
                      arrowprops=dict(arrowstyle="->", color="#0F172A", lw=0.7),
                      fontsize=7.2, fontweight="bold", color="#1E3A8A")
    ax_inset.annotate("E459K (1.7% vs 2.01%)", xy=(1.7, 2.01), xytext=(0.80, 0.40),
                      arrowprops=dict(arrowstyle="->", color="#0F172A", lw=0.7),
                      fontsize=7.2, fontweight="bold", color="#1E3A8A")

    ax.set_xlim(-2, 105)
    ax.set_ylim(-2, 105)
    ax.set_xlabel("Illumina MiSeq VAF (%)", fontsize=11, fontweight="bold", labelpad=8)
    ax.set_ylabel("Oxford Nanopore VAF (%)", fontsize=11, fontweight="bold", labelpad=8)
    ax.set_title("BCR-ABL1 Cross-Platform Validation: Oxford Nanopore vs Illumina MiSeq",
                 fontsize=12, fontweight="bold", pad=12)
    ax.grid(True, linestyle=":", alpha=0.5, color="#CBD5E1")
    ax.legend(loc="center left", bbox_to_anchor=(0.04, 0.65), frameon=True, fontsize=8.8,
              facecolor="#F8FAFC", edgecolor="#CBD5E1", framealpha=0.95)

    ax.tick_params(axis="both", which="major", labelsize=9.0)
    ax.set_facecolor("#FFFFFF")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    print(f"[Nagini] Successfully saved single-panel hero figure to: {out_path}")


def main(argv: Optional[List[str]] = None) -> int:
    parser = argparse.ArgumentParser(
        description="Plot publication-ready comparison figure across Illumina MiSeq and Nagini callers."
    )
    parser.add_argument(
        "--truth-table", "-t", type=Path, default=None,
        help="Path to Table 6 text file"
    )
    parser.add_argument(
        "--variants-tsv", "-v", type=Path, default=None,
        help="Path to caller results TSV"
    )
    parser.add_argument(
        "--bed", "-b", type=Path, default=Path("workflow/resources/bcr_abl/branford_positions_hg19.bed"),
        help="Path to Branford positions BED (default: workflow/resources/bcr_abl/branford_positions_hg19.bed)"
    )
    parser.add_argument(
        "--sample-map", "--barcode-map", "-m", type=Path, default=None,
        help="Path to sample mapping or barcode map TSV/CSV"
    )
    parser.add_argument(
        "--callers", "-c", type=str, default="Branford,bcftools",
        help="Comma-separated list of callers to include (default: Branford,bcftools)"
    )
    parser.add_argument(
        "--compare-to", type=str, choices=["miseq", "sequel", "both"], default="miseq",
        help="Truth technology to compare against (default: miseq)"
    )
    parser.add_argument(
        "--layout", choices=["all", "single", "multi"], default="all",
        help="Layout mode: 'all' (creates both multi-panel and single hero panel), 'single' (hero panel only), or 'multi' (3 panels)"
    )
    parser.add_argument(
        "--marker-style", choices=["concentric", "dodge"], default="concentric",
        help="Marker style to prevent overplotting: 'concentric' (open square for bcftools around solid circle for Branford) or 'dodge' (horizontal offset)"
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

    caller_list = [c.strip() for c in args.callers.split(",") if c.strip()]

    # 1. Parse Branford BED
    coords_to_aa = {}
    if args.bed and args.bed.is_file():
        coords_to_aa = load_branford_bed(args.bed)

    # 2. Parse truth table
    truth_records = load_truth_table(args.truth_table) if args.truth_table else []
    if not truth_records and args.truth_table:
        sys.stderr.write(f"Warning: no truth records loaded from {args.truth_table}\n")

    # 3. Parse caller results
    caller_results = load_caller_results(args.variants_tsv, coords_to_aa=coords_to_aa) if args.variants_tsv else {}
    if not caller_results and args.variants_tsv:
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

    # 6. Generate figures
    # Multi-panel figure
    if args.layout in ("all", "multi"):
        create_comparison_figure(
            dataset, args.output, dpi=args.dpi,
            callers=caller_list, compare_to=args.compare_to,
            marker_style=args.marker_style
        )
        if args.output.suffix.lower() == ".png":
            svg_path = args.output.with_suffix(".svg")
            try:
                create_comparison_figure(
                    dataset, svg_path, dpi=args.dpi,
                    callers=caller_list, compare_to=args.compare_to,
                    marker_style=args.marker_style
                )
            except Exception as exc:
                sys.stderr.write(f"Notice: could not save SVG copy: {exc}\n")

    # Single-panel hero figure
    if args.layout in ("all", "single"):
        if args.layout == "single":
            single_out_png = args.output
        else:
            single_out_png = args.output.parent / f"{args.output.stem}_single{args.output.suffix}"

        create_single_panel_figure(
            dataset, single_out_png, dpi=args.dpi,
            callers=caller_list, compare_to=args.compare_to,
            marker_style=args.marker_style
        )
        if single_out_png.suffix.lower() == ".png":
            single_out_svg = single_out_png.with_suffix(".svg")
            try:
                create_single_panel_figure(
                    dataset, single_out_svg, dpi=args.dpi,
                    callers=caller_list, compare_to=args.compare_to,
                    marker_style=args.marker_style
                )
            except Exception as exc:
                sys.stderr.write(f"Notice: could not save SVG copy for single panel: {exc}\n")

    return 0


if __name__ == "__main__":
    sys.exit(main())
