"""Generate publication-ready PDF diagrams for BCR-ABL1 gene fusions."""

import argparse
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

# Set writable MPLCONFIGDIR before importing matplotlib
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_nagini")

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.patches as patches
    import matplotlib.pyplot as plt
    from matplotlib.path import Path as MplPath
    HAS_MATPLOTLIB = True
except Exception:
    HAS_MATPLOTLIB = False


def parse_tsv(path: Optional[str]) -> List[Dict[str, str]]:
    if not path or not os.path.isfile(path):
        return []
    rows = []
    with open(path, "r", encoding="utf-8") as handle:
        header = None
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields)))
    return rows


def draw_matplotlib_diagram(
    sample: str,
    output_pdf: str,
    fusion_rows: List[Dict[str, str]],
    mutation_rows: List[Dict[str, str]],
    phasing_rows: List[Dict[str, str]],
    exon_cov_rows: List[Dict[str, str]],
    build: str = "hg19",
) -> None:
    fig = plt.figure(figsize=(12, 8), dpi=300)
    gs = fig.add_gridspec(3, 2, height_ratios=[1.2, 4.0, 2.8], width_ratios=[1.2, 1.0], hspace=0.35, wspace=0.25)

    # -----------------------------------------------------------------------
    # 1. Header & Summary Banner
    # -----------------------------------------------------------------------
    ax_header = fig.add_subplot(gs[0, :])
    ax_header.axis("off")

    total_fusions = sum(int(r.get("read_count", 0)) for r in fusion_rows)
    prim_row = next((r for r in fusion_rows if r.get("status") == "PRIMARY"), None)
    prim_str = f"{prim_row['fusion_isoform']} ({prim_row.get('alt_name', '')}) — {prim_row.get('fraction_percent', '0')}% ({prim_row.get('read_count', '0')} reads)" if prim_row else "No primary fusion detected"

    det_muts = [r for r in mutation_rows if str(r.get("status", "")).upper() == "DETECTED"]
    mut_str = ", ".join(f"{r['mutation']} (VAF {float(r.get('vaf_percent', r.get('vaf', 0))):.1f}%)" for r in det_muts) if det_muts else "None detected"

    cis_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in phasing_rows if "CIS" in r.get("phase_state", "").upper()]
    phase_str = f"CIS ({', '.join(cis_pairs)})" if cis_pairs else ("TRANS" if phasing_rows else "N/A")

    ax_header.text(0.0, 0.85, "NAGINI BCR-ABL1 FUSION & ARCHITECTURE REPORT", fontsize=14, fontweight="bold", color="#14324F")
    ax_header.text(0.0, 0.45, f"Sample: {sample}   |   Genome Build: {build}   |   Total Fusion Reads: {total_fusions:,}", fontsize=11, fontweight="bold", color="#1F4E78")
    
    status_box_color = "#C6EFCE" if (prim_row and not det_muts) else ("#FFC7CE" if det_muts or cis_pairs else "#D9E1F2")
    status_text_color = "#006100" if (prim_row and not det_muts) else ("#9C0006" if det_muts or cis_pairs else "#1F4E78")
    
    summary_text = f"Primary Isoform: {prim_str}\nResistance Mutations: {mut_str}   |   Phasing: {phase_str}"
    bbox_props = dict(boxstyle="round,pad=0.5", facecolor=status_box_color, edgecolor=status_text_color, alpha=0.9, lw=1.2)
    ax_header.text(0.0, -0.05, summary_text, fontsize=9.5, color=status_text_color, fontweight="semibold", bbox=bbox_props)

    # -----------------------------------------------------------------------
    # 2. Gene Architecture & Fusion Chimeras Track
    # -----------------------------------------------------------------------
    ax_gene = fig.add_subplot(gs[1, :])
    ax_gene.set_xlim(-5, 105)
    ax_gene.set_ylim(-10, 110)
    ax_gene.axis("off")

    # Draw BCR Track (Top, Y=80)
    bcr_y = 80
    ax_gene.plot([0, 95], [bcr_y, bcr_y], color="#7F7F7F", lw=3, zorder=1)
    ax_gene.text(-4, bcr_y, "BCR\n(chr22)", fontsize=10, fontweight="bold", color="#1F4E78", va="center", ha="right")

    # BCR Exons
    bcr_exons = [
        ("e1", 5, 8, "#41719C"),
        ("e2", 18, 5, "#41719C"),
        ("e6", 30, 5, "#41719C"),
        ("e12", 42, 5, "#41719C"),
        ("e13", 52, 6, "#2E75B6"),  # b2
        ("e14", 62, 6, "#1F4E78"),  # b3
        ("e19", 78, 6, "#41719C"),  # micro
        ("e23", 90, 5, "#41719C"),
    ]
    bcr_exon_centers = {}
    for name, x, w, col in bcr_exons:
        rect = patches.FancyBboxPatch((x, bcr_y - 6), w, 12, boxstyle="round,pad=0.2", facecolor=col, edgecolor="#14324F", lw=1.2, zorder=3)
        ax_gene.add_patch(rect)
        ax_gene.text(x + w/2, bcr_y, name, color="white", fontsize=8, fontweight="bold", ha="center", va="center", zorder=4)
        bcr_exon_centers[name] = x + w/2

    # Breakpoint clusters
    ax_gene.text(9, bcr_y + 10, "m-bcr (Minor)", fontsize=7.5, color="#595959", ha="center", style="italic")
    ax_gene.text(58, bcr_y + 10, "M-bcr (Major)", fontsize=8, color="#C00000", ha="center", fontweight="bold")
    ax_gene.text(81, bcr_y + 10, "μ-bcr (Micro)", fontsize=7.5, color="#595959", ha="center", style="italic")

    # Draw ABL1 Track (Bottom, Y=20)
    abl_y = 20
    ax_gene.plot([0, 95], [abl_y, abl_y], color="#7F7F7F", lw=3, zorder=1)
    ax_gene.text(-4, abl_y, "ABL1\n(chr9)", fontsize=10, fontweight="bold", color="#A61C1C", va="center", ha="right")

    # ABL1 Exons
    abl_exons = [
        ("a2", 8, 7, "#C00000", False),
        ("a3", 20, 6, "#C00000", False),
        ("a4", 32, 6, "#ED7D31", True),  # Kinase domain start
        ("a5", 42, 6, "#ED7D31", True),
        ("a6", 52, 6, "#ED7D31", True),
        ("a7", 62, 6, "#ED7D31", True),
        ("a8", 72, 6, "#ED7D31", True),
        ("a9", 82, 6, "#ED7D31", True),  # Kinase domain end
        ("a10", 91, 5, "#C00000", False),
    ]
    abl_exon_centers = {}
    for name, x, w, col, is_kd in abl_exons:
        edge = "#70AD47" if is_kd else "#800000"
        rect = patches.FancyBboxPatch((x, abl_y - 6), w, 12, boxstyle="round,pad=0.2", facecolor=col, edgecolor=edge, lw=1.5 if is_kd else 1.0, zorder=3)
        ax_gene.add_patch(rect)
        ax_gene.text(x + w/2, abl_y, name, color="white", fontsize=8, fontweight="bold", ha="center", va="center", zorder=4)
        abl_exon_centers[name] = x + w/2

    # Kinase domain label
    ax_gene.text(57, abl_y - 12, "Kinase Domain (Exons 4–9 / Codons 242–493)", fontsize=8, color="#ED7D31", ha="center", fontweight="bold")

    # Highlight detected mutations on ABL1 exons
    for mut in det_muts:
        mname = mut.get("mutation", "")
        # Map codon position to exon
        target_center = abl_exon_centers.get("a6", 55) # default around T315I (exon 6)
        if any(c in mname for c in ("250", "252", "253", "255")):
            target_center = abl_exon_centers.get("a4", 35) # P-loop
        elif any(c in mname for c in ("315", "317")):
            target_center = abl_exon_centers.get("a6", 55) # Gatekeeper
        elif any(c in mname for c in ("359", "351", "355")):
            target_center = abl_exon_centers.get("a7", 65)
        elif any(c in mname for c in ("382", "384", "387", "396")):
            target_center = abl_exon_centers.get("a8", 75)

        ax_gene.annotate(
            f"★ {mname}",
            xy=(target_center, abl_y + 7),
            xytext=(target_center, abl_y + 20),
            arrowprops=dict(facecolor="#C00000", edgecolor="#800000", arrowstyle="->", lw=1.5),
            fontsize=8.5,
            fontweight="bold",
            color="#C00000",
            ha="center",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#FFEB9C", edgecolor="#C00000", lw=1),
            zorder=6,
        )

    # Draw Fusion Ribbons connecting BCR to ABL1
    active_fusions = [r for r in fusion_rows if int(r.get("read_count", 0)) > 0]
    if not active_fusions:
        ax_gene.text(50, 50, "No chimeric fusion junctions detected", fontsize=12, color="#595959", ha="center", style="italic")
    else:
        for r in active_fusions:
            iso = r["fusion_isoform"]
            cnt = int(r.get("read_count", 0))
            pct = float(r.get("fraction_percent", 0.0))
            status = r.get("status", "SECONDARY")

            bcr_ex = iso.split("a")[0] if "a" in iso else "e14"
            abl_ex = "a" + iso.split("a")[1] if "a" in iso else "a2"

            x_bcr = bcr_exon_centers.get(bcr_ex, 65)
            x_abl = abl_exon_centers.get(abl_ex, 11)

            # Draw Bezier Curve
            verts = [
                (x_bcr, bcr_y - 6),
                (x_bcr, (bcr_y + abl_y) / 2),
                (x_abl, (bcr_y + abl_y) / 2),
                (x_abl, abl_y + 6),
            ]
            codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4]
            path = MplPath(verts, codes)

            lw = min(max(cnt / max(total_fusions, 1) * 6.0, 1.5), 8.0)
            color = "#1F4E78" if status == "PRIMARY" else "#8EA9DB"
            alpha = 0.85 if status == "PRIMARY" else 0.55

            patch = patches.PathPatch(path, facecolor="none", edgecolor=color, lw=lw, alpha=alpha, zorder=2)
            ax_gene.add_patch(patch)

            # Label on primary fusion junction
            if status == "PRIMARY" or pct >= 10.0:
                mid_x = (x_bcr + x_abl) / 2
                mid_y = (bcr_y + abl_y) / 2
                ax_gene.text(
                    mid_x, mid_y,
                    f"{iso} ({pct:.1f}%)\n{cnt:,} reads",
                    fontsize=8.5,
                    fontweight="bold",
                    color=color,
                    ha="center",
                    va="center",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="#FFFFFF", edgecolor=color, lw=1.2, alpha=0.9),
                    zorder=5,
                )

    # -----------------------------------------------------------------------
    # 3. Lower Left: Isoform Abundance Bar Chart
    # -----------------------------------------------------------------------
    ax_bar = fig.add_subplot(gs[2, 0])
    iso_names = [r["fusion_isoform"] for r in fusion_rows if int(r.get("read_count", 0)) > 0]
    iso_counts = [int(r.get("read_count", 0)) for r in fusion_rows if int(r.get("read_count", 0)) > 0]
    iso_pcts = [float(r.get("fraction_percent", 0.0)) for r in fusion_rows if int(r.get("read_count", 0)) > 0]
    colors = ["#1F4E78" if r.get("status") == "PRIMARY" else "#8EA9DB" for r in fusion_rows if int(r.get("read_count", 0)) > 0]

    if iso_names:
        bars = ax_bar.barh(iso_names, iso_counts, color=colors, edgecolor="#14324F", height=0.55)
        ax_bar.set_xlabel("Fusion Read Count", fontsize=9, fontweight="bold")
        ax_bar.set_title("Isoform Quantification", fontsize=10, fontweight="bold", color="#1F4E78")
        ax_bar.grid(axis="x", linestyle="--", alpha=0.5)
        for bar, pct, cnt in zip(bars, iso_pcts, iso_counts):
            ax_bar.text(bar.get_width() + max(iso_counts)*0.02, bar.get_y() + bar.get_height()/2, f"{pct:.1f}% ({cnt:,})", va="center", fontsize=8, fontweight="bold", color="#14324F")
        ax_bar.set_xlim(0, max(iso_counts) * 1.3 if iso_counts else 10)
    else:
        ax_bar.text(0.5, 0.5, "No fusion reads detected", ha="center", va="center", fontsize=10, color="#7F7F7F")
        ax_bar.axis("off")

    # -----------------------------------------------------------------------
    # 4. Lower Right: Exon Coverage Profile
    # -----------------------------------------------------------------------
    ax_cov = fig.add_subplot(gs[2, 1])
    if exon_cov_rows:
        import re
        ex_names = [r.get("exon", "") for r in exon_cov_rows]
        ex_depths = [float(r.get("mean_depth", 0.0)) for r in exon_cov_rows]
        
        def _is_kd_exon(ex_str):
            if "ABL1" in ex_str.upper():
                m = re.search(r"(?:exon|e|a)(\d+)", ex_str, re.IGNORECASE)
                if m:
                    return int(m.group(1)) in range(4, 10)
            return False

        cov_colors = ["#ED7D31" if _is_kd_exon(ex) else "#41719C" for ex in ex_names]
        short_names = [re.sub(r"(?i)(BCR|ABL1)[_ ]+Exon(\d+).*", r"\1 e\2", ex) for ex in ex_names]
        short_names = [re.sub(r"(?i)ABL1 e(\d+)", r"ABL1 a\1", s) for s in short_names]
        
        ax_cov.bar(range(len(short_names)), ex_depths, color=cov_colors, edgecolor="#14324F", width=0.6)
        ax_cov.set_xticks(range(len(short_names)))
        ax_cov.set_xticklabels(short_names, rotation=45, ha="right", fontsize=7.5)
        ax_cov.set_ylabel("Mean Depth (x)", fontsize=9, fontweight="bold")
        ax_cov.set_title("Exon Mean Coverage Depth", fontsize=10, fontweight="bold", color="#1F4E78")

        ax_cov.grid(axis="y", linestyle="--", alpha=0.5)
    else:
        ax_cov.text(0.5, 0.5, "Exon coverage data not available", ha="center", va="center", fontsize=10, color="#7F7F7F")
        ax_cov.axis("off")

    # Save PDF
    plt.savefig(output_pdf, bbox_inches="tight")
    plt.close()


def generate_fallback_pdf(sample: str, output_pdf: str, message: str) -> None:
    """Generate a clean 1-page PDF using basic R or Python if matplotlib is missing."""
    import subprocess
    r_script = f"""
    pdf("{output_pdf}", width=11, height=7)
    plot(0:10, 0:10, type="n", axes=FALSE, xlab="", ylab="", main="NAGINI BCR-ABL1 FUSION REPORT")
    text(5, 7.5, "{sample}", cex=1.8, font=2, col="#14324F")
    text(5, 5.0, "{message}", cex=1.3, col="#1F4E78")
    text(5, 3.0, "Detailed variant and isoform metrics are available in the Excel report", cex=1.0, col="#595959")
    dev.off()
    """
    try:
        subprocess.run(["Rscript", "-e", r_script], check=True, capture_output=True)
    except Exception:
        # If Rscript also fails, write minimal valid PDF header
        with open(output_pdf, "wb") as f:
            f.write(b"%PDF-1.4\n1 0 obj<</Type/Catalog/Pages 2 0 R>>endobj\n2 0 obj<</Type/Pages/Kids[3 0 R]/Count 1>>endobj\n3 0 obj<</Type/Page/MediaBox[0 0 612 792]/Parent 2 0 R/Resources<<>>>>endobj\nxref\n0 4\n0000000000 65535 f\n0000000009 00000 n\n0000000052 00000 n\n0000000101 00000 n\ntrailer<</Size 4/Root 1 0 R>>\nstartxref\n178\n%%EOF\n")


def main():
    parser = argparse.ArgumentParser(description="Draw BCR-ABL1 fusion diagram")
    parser.add_argument("--sample", required=True, help="Sample ID")
    parser.add_argument("--fusions", required=True, help="Path to fusion_summary.tsv")
    parser.add_argument("--mutations", default="", help="Path to branford_mutations.tsv")
    parser.add_argument("--phasing", default="", help="Path to phasing.tsv")
    parser.add_argument("--exon-coverage", default="", help="Path to exon_coverage.tsv")
    parser.add_argument("--output", required=True, help="Output PDF path")
    parser.add_argument("--build", default="hg19", help="Genome build")

    args = parser.parse_args()

    fusion_rows = parse_tsv(args.fusions)
    mutation_rows = parse_tsv(args.mutations)
    phasing_rows = parse_tsv(args.phasing)
    exon_cov_rows = parse_tsv(args.exon_coverage)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)

    if HAS_MATPLOTLIB:
        draw_matplotlib_diagram(
            sample=args.sample,
            output_pdf=args.output,
            fusion_rows=fusion_rows,
            mutation_rows=mutation_rows,
            phasing_rows=phasing_rows,
            exon_cov_rows=exon_cov_rows,
            build=args.build,
        )
    else:
        generate_fallback_pdf(
            sample=args.sample,
            output_pdf=args.output,
            message="BCR-ABL1 Analysis Completed",
        )

    print(f"Generated BCR-ABL1 fusion plot: {args.output} ({os.path.getsize(args.output):,} bytes)")


if __name__ == "__main__":
    main()
