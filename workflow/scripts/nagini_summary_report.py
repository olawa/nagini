"""Create a multi-sample XLSX summary report for a Nagini run."""

import argparse
import os
import re
import shutil
import sys
from datetime import datetime, timezone

from pathlib import Path
from typing import Any, Dict, List, Optional

import xlsxwriter

# Import shared functions and styles from nagini_report
try:
    from nagini_report import (
        LAYOUT,
        STYLE_SPECS,
        TABLE,
        _find_metric,
        _is_pathogenic_clinvar,
        build_bcr_abl_sheet,
        parse_bcr_abl_tsv,
        parse_cramino,
        parse_small_variant_vcf,
        parse_tsv,
        safe_sheet_name,
        write_xlsx,
    )
except ImportError:
    # Direct script execution fallback
    sys.path.insert(0, os.path.dirname(__file__))
    from nagini_report import (
        LAYOUT,
        STYLE_SPECS,
        TABLE,
        _find_metric,
        _is_pathogenic_clinvar,
        build_bcr_abl_sheet,
        parse_bcr_abl_tsv,
        parse_cramino,
        parse_small_variant_vcf,
        parse_tsv,
        safe_sheet_name,
        write_xlsx,
    )


def build_summary_report(
    output_path: str,
    sequence_id: str,
    reference_build: str,
    sample_records: List[Dict[str, Any]],
) -> None:
    """Build the multi-sample workbook."""
    sheets = []
    used_sheet_names = set()

    out_dir = os.path.dirname(output_path) or "."
    parent_dir = os.path.basename(os.path.abspath(os.path.join(out_dir, "..")))
    seq_clean = str(sequence_id or "").strip()
    display_id = seq_clean if (seq_clean and seq_clean.lower() not in {"run", "default"}) else parent_dir
    date_str = datetime.now().strftime("%Y-%m-%d")

    # -----------------------------------------------------------------------
    # 1. OVERVIEW SHEET
    # -----------------------------------------------------------------------
    overview_rows = [
        [(f"NAGINI RUN SUMMARY - {display_id} ({date_str})", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2)],
        [("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0)],
        [("RUN & PIPELINE METADATA", 3), ("", 3), ("", 3), ("", 3)],
        [("Sequence ID / Run", 7), (display_id, 0), ("Sequencing Platform", 7), ("Oxford Nanopore (ONT)", 0)],
        [("Reference Genome", 7), (reference_build, 0), ("Total Samples", 7), (len(sample_records), 0)],
        [("Report Generated", 7), (datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"), 0), ("Pipeline Status", 7), ("COMPLETED", 5)],

        [("", 0), ("", 0), ("", 0), ("", 0)],
        [("SAMPLE QUALITY CONTROL & SEQUENCING SUMMARY", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3)],
        [
            ("Sample ID", 1),
            ("Type", 1),
            ("Barcode", 1),
            ("Assay", 1),
            ("Yield [Gb]", 1),
            ("Total Reads", 1),
            ("Mean Q", 1),
            ("N50 [bp]", 1),
            ("Median [bp]", 1),
            ("Mean Target / KD Depth", 1),
            ("BCR-ABL1 Fusion Reads", 1),
            ("BCR-ABL1 Enrichment", 1),
            ("Primary Isoform", 1),
            ("Resistance Mutations", 1),
            ("Compound Phasing", 1),
            ("Status", 1),
        ],
    ]

    has_any_bcr_abl = False
    all_pass_variants = []
    bcr_abl_samples_data = []

    for s in sample_records:
        sample_name = s["sample"]
        sample_type = s["type"]
        barcode = s.get("barcode", "")
        assay = s.get("assay", "wgs")
        cramino_path = s.get("cramino", "")
        wgs_cov_path = s.get("wgs_coverage", "")
        snv_vcf_path = s.get("snv_vcf", "")
        fusion_tsv_path = s.get("bcr_abl_fusion", "")
        mut_tsv_path = s.get("bcr_abl_mutations", "")
        phase_tsv_path = s.get("bcr_abl_phasing", "")
        exon_cov_path = s.get("bcr_abl_exon_coverage", "")

        # Parse cramino metrics
        _, cramino_rows = parse_cramino(cramino_path)
        cramino_dict = {row["Metric"]: row["Value"] for row in cramino_rows}
        yield_gb = _find_metric(cramino_dict, "Yield [Gb]", "Yield", "Bases")
        total_reads_str = _find_metric(cramino_dict, "Number of reads", "Total reads", "Read count")
        mean_q = _find_metric(cramino_dict, "Mean read quality", "Mean quality", "Average read quality", "Mean base quality", "Median read quality", "Quality", "Q-score")
        n50_bp = _find_metric(cramino_dict, "Read length N50", "N50", "N50 read length")
        median_len = _find_metric(cramino_dict, "Median read length", "Median length", "Median")


        # Parse BCR-ABL1 metrics if present
        fusion_rows = parse_bcr_abl_tsv(fusion_tsv_path)
        mutation_rows = parse_bcr_abl_tsv(mut_tsv_path)
        phasing_rows = parse_bcr_abl_tsv(phase_tsv_path)
        exon_rows = parse_bcr_abl_tsv(exon_cov_path)
        is_bcr_abl = bool(fusion_rows or mutation_rows or phasing_rows or exon_rows or assay == "bcr_abl")

        if is_bcr_abl:
            has_any_bcr_abl = True
            bcr_abl_samples_data.append({
                "sample": sample_name,
                "type": sample_type,
                "barcode": barcode,
                "fusion_rows": fusion_rows,
                "mutation_rows": mutation_rows,
                "phasing_rows": phasing_rows,
                "exon_rows": exon_rows,
                "plot_png": s.get("bcr_abl_plot_png", ""),
            })


        total_fusion_reads = sum(int(r.get("read_count", 0)) for r in fusion_rows)
        raw_reads_num = 0
        try:
            raw_reads_num = int(float(str(total_reads_str).replace(",", "")))
        except (ValueError, TypeError):
            raw_reads_num = 0

        target_enrich_str = "N/A"
        if raw_reads_num > 0 and total_fusion_reads > 0:
            pct = (total_fusion_reads / raw_reads_num) * 100.0
            target_enrich_str = f"{pct:.2f}%"

        # Target / Kinase coverage
        kd_depth_str = "N/A"
        if exon_rows:
            tot_len = sum(int(r.get("length", 0)) for r in exon_rows)
            if tot_len > 0:
                avg_d = sum(float(r.get("mean_depth", 0.0)) * int(r.get("length", 0)) for r in exon_rows) / tot_len
                kd_depth_str = f"{avg_d:.1f}x"
        elif mutation_rows:
            avg_d = sum(int(r.get("depth", 0)) for r in mutation_rows) / len(mutation_rows)
            kd_depth_str = f"{avg_d:.1f}x"
        elif wgs_cov_path and os.path.exists(wgs_cov_path):
            _, wgs_rows = parse_tsv(wgs_cov_path)
            for wr in wgs_rows:
                if str(wr.get("chrom", "")).lower() == "total":
                    kd_depth_str = f"{float(wr.get('mean', 0.0)):.2f}x"
                    break

        prim_iso = next((r["fusion_isoform"] for r in fusion_rows if r.get("status") == "PRIMARY"), ("None detected" if is_bcr_abl else "-"))
        det_muts = [r["mutation"] for r in mutation_rows if str(r.get("status", "")).upper() == "DETECTED"]
        mut_str = ", ".join(det_muts) if det_muts else ("None detected" if is_bcr_abl else "-")

        cis_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in phasing_rows if "CIS" in r.get("phase_state", "").upper()]
        trans_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in phasing_rows if "TRANS" in r.get("phase_state", "").upper()]
        phase_str = "-"
        if cis_pairs:
            phase_str = f"CIS ({', '.join(cis_pairs)})"
        elif trans_pairs:
            phase_str = f"TRANS ({', '.join(trans_pairs)})"
        elif is_bcr_abl:
            phase_str = "N/A (<2 muts)"

        sample_status, sample_style = ("PASS", 5)
        if det_muts or cis_pairs:
            sample_status, sample_style = ("ACTION REQUIRED", 4)
        elif not yield_gb or (raw_reads_num < 100):
            sample_status, sample_style = ("REVIEW", 6)

        overview_rows.append([
            (sample_name, 10),
            (sample_type, 0),
            (barcode or "-", 0),
            (assay, 0),
            (yield_gb or "-", 0),
            (total_reads_str or "-", 0),
            (mean_q or "-", 0),
            (n50_bp or "-", 0),
            (median_len or "-", 0),
            (kd_depth_str, 0),
            (total_fusion_reads if is_bcr_abl else "-", 0),
            (target_enrich_str if is_bcr_abl else "-", 0),
            (prim_iso, 8 if prim_iso not in ("-", "None detected") else 0),
            (mut_str, 4 if det_muts else 0),
            (phase_str, 4 if cis_pairs else (6 if trans_pairs else 0)),
            (sample_status, sample_style),
        ])

        # Collect PASS variants for PASS Variants tab
        if mutation_rows:
            for r in mutation_rows:
                if str(r.get("status", "")).upper() == "DETECTED":
                    all_pass_variants.append({
                        "sample": sample_name,
                        "type": sample_type,
                        "gene": "ABL1",
                        "chrom": r.get("chrom", "chr9"),
                        "pos": r.get("pos", ""),
                        "ref": r.get("ref", ""),
                        "alt": r.get("alt", ""),
                        "vaf": f"{float(r.get('vaf_percent', r.get('vaf', 0.0))):.2f}%",
                        "depth": r.get("depth", ""),
                        "alt_reads": r.get("ad_alt", r.get("alt_reads", "")),
                        "consequence": r.get("mutation", ""),

                        "clinical_note": r.get("tki_resistance", "Resistance mutation"),
                        "status": "DETECTED",
                    })

        if snv_vcf_path and os.path.exists(snv_vcf_path):
            _, snvs, _ = parse_small_variant_vcf(snv_vcf_path)
            for v in snvs:
                if v.get("Filter") in {"PASS", "."}:
                    all_pass_variants.append({
                        "sample": sample_name,
                        "type": sample_type,
                        "gene": v.get("Gene", ""),
                        "chrom": v.get("Chromosome", ""),
                        "pos": v.get("Position", ""),
                        "ref": v.get("Reference", ""),
                        "alt": v.get("Alternate", ""),
                        "vaf": v.get("VAF", ""),
                        "depth": v.get("Depth", ""),
                        "alt_reads": v.get("Alt_Reads", ""),
                        "consequence": v.get("Consequence", ""),
                        "clinical_note": v.get("ClinVar_Significance", ""),
                        "status": "PASS",
                    })

    sheets.append(("Overview", overview_rows, LAYOUT))

    # -----------------------------------------------------------------------
    # 2. PASS VARIANTS SHEET
    # -----------------------------------------------------------------------
    pass_var_rows = [
        [
            ("Sample ID", 1),
            ("Type", 1),
            ("Gene", 1),
            ("Chr", 1),
            ("Position", 1),
            ("Ref", 1),
            ("Alt", 1),
            ("VAF", 1),
            ("Depth", 1),
            ("Alt Reads", 1),
            ("Consequence / AA Change", 1),
            ("Clinical Significance / Resistance Profile", 1),
            ("Status", 1),
        ]
    ]
    for pv in all_pass_variants:
        pass_var_rows.append([
            (pv["sample"], 10),
            (pv["type"], 0),
            (pv["gene"], 0),
            (pv["chrom"], 0),
            (pv["pos"], 0),
            (pv["ref"], 0),
            (pv["alt"], 0),
            (pv["vaf"], 0),
            (pv["depth"], 0),
            (pv["alt_reads"], 0),
            (pv["consequence"], 0),
            (pv["clinical_note"], 4 if _is_pathogenic_clinvar(pv["clinical_note"]) or "Resistant" in pv["clinical_note"] else 0),
            (pv["status"], 5),
        ])
    if len(pass_var_rows) == 1:
        # An assay with snv_calling off produces no VCF at all. Reporting that
        # as "none detected" reads as a negative finding rather than as a
        # question this run did not ask.
        assessed = any(
            record.get("snv_vcf") or record.get("bcr_abl_mutations")
            for record in sample_records
        )
        note = (
            "No high-confidence PASS variants detected across samples"
            if assessed
            else "Small variant calling was not part of this run's assay; no variants were assessed"
        )
        pass_var_rows.append([(note, 0)] + [("", 0)] * 12)
    sheets.append(("PASS Variants", pass_var_rows, TABLE))

    # -----------------------------------------------------------------------
    # 3. MULTI-SAMPLE BCR-ABL1 OVERVIEW (if BCR-ABL1 assay data present)
    # -----------------------------------------------------------------------
    if has_any_bcr_abl:
        bcr_summary_rows = [
            [("MULTI-SAMPLE BCR-ABL1 ASSAY COMPARISON", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2), ("", 2)],
            [("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0)],
            [("1. FUSION ISOFORM COMPARISON ACROSS SAMPLES", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3)],
            [
                ("Sample ID", 1),
                ("Barcode", 1),
                ("Primary Isoform", 1),
                ("e14a2 (b3a2) Reads", 1),
                ("e14a2 %", 1),
                ("e13a2 (b2a2) Reads", 1),
                ("e13a2 %", 1),
                ("e1a2 Reads", 1),
                ("e1a2 %", 1),
                ("Total Fusion Reads", 1),
            ],
        ]

        for bd in bcr_abl_samples_data:
            f_rows = bd["fusion_rows"]
            iso_map = {r["fusion_isoform"]: (int(r.get("read_count", 0)), float(r.get("fraction_percent", 0.0))) for r in f_rows}
            tot_reads = sum(int(r.get("read_count", 0)) for r in f_rows)
            prim = next((r["fusion_isoform"] for r in f_rows if r.get("status") == "PRIMARY"), "None")

            e14_cnt, e14_pct = iso_map.get("e14a2", (0, 0.0))
            e13_cnt, e13_pct = iso_map.get("e13a2", (0, 0.0))
            e1a_cnt, e1a_pct = iso_map.get("e1a2", (0, 0.0))

            bcr_summary_rows.append([
                (bd["sample"], 10),
                (bd["barcode"], 0),
                (prim, 8 if prim != "None" else 0),
                (e14_cnt, 0),
                (f"{e14_pct:.1f}%", 0),
                (e13_cnt, 0),
                (f"{e13_pct:.1f}%", 0),
                (e1a_cnt, 0),
                (f"{e1a_pct:.1f}%", 0),
                (tot_reads, 0),
            ])

        bcr_summary_rows.append([("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0), ("", 0)])
        bcr_summary_rows.append([("2. RESISTANCE MUTATIONS & COMPOUND PHASING SUMMARY", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3)])
        bcr_summary_rows.append([
            ("Sample ID", 1),
            ("Detected Mutations", 1),
            ("VAF Details", 1),
            ("Compound Phasing", 1),
            ("Cis Pair(s)", 1),
            ("Trans Pair(s)", 1),
            ("Clinical Assessment", 1),
            ("", 1), ("", 1), ("", 1)
        ])

        for bd in bcr_abl_samples_data:
            m_rows = bd["mutation_rows"]
            p_rows = bd["phasing_rows"]
            det_m = [r for r in m_rows if str(r.get("status", "")).upper() == "DETECTED"]
            subcl_m = [r for r in m_rows if str(r.get("status", "")).upper() == "SUBCLONAL"]

            if det_m:
                mut_names = ", ".join(r["mutation"] for r in det_m)
                if subcl_m:
                    mut_names += f" (Subcl: {', '.join(r['mutation'] for r in subcl_m)})"
            elif subcl_m:
                mut_names = f"Subclonal: {', '.join(r['mutation'] for r in subcl_m)}"
            else:
                mut_names = "None detected"

            all_active = det_m + subcl_m
            vaf_details = ", ".join(f"{r['mutation']}: {float(r.get('vaf_percent', r.get('vaf', 0.0))):.1f}% ({r.get('ad_alt', r.get('alt_reads', 0))}/{r.get('depth')})" for r in all_active) if all_active else "-"

            cis_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in p_rows if "CIS" in r.get("phase_state", "").upper()]
            trans_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in p_rows if "TRANS" in r.get("phase_state", "").upper()]
            phase_state = "PASS (No resistance muts)" if not all_active else ("SINGLE MUTATION" if len(all_active) == 1 else ("COMPOUND (CIS)" if cis_pairs else "POLYCLONAL (TRANS)"))

            clinical_note = "Normal / Sensitive" if not all_active else ("Compound mutation (high TKI resistance)" if cis_pairs else ("Resistance mutation detected" if det_m else "Subclonal resistance mutation"))

            bcr_summary_rows.append([
                (bd["sample"], 10),
                (mut_names, 4 if det_m else (6 if subcl_m else 0)),
                (vaf_details, 0),
                (phase_state, 4 if cis_pairs else (6 if trans_pairs else (4 if det_m else (6 if subcl_m else 5)))),
                (", ".join(cis_pairs) if cis_pairs else "-", 4 if cis_pairs else 0),
                (", ".join(trans_pairs) if trans_pairs else "-", 6 if trans_pairs else 0),
                (clinical_note, 4 if det_m else (6 if subcl_m else 5)),
                ("", 0), ("", 0), ("", 0),
            ])

        sheets.append(("BCR-ABL1 Overview", bcr_summary_rows, LAYOUT))

        # -------------------------------------------------------------------
        # 4. PER-SAMPLE BCR-ABL1 DETAIL SHEETS
        # -------------------------------------------------------------------
        for bd in bcr_abl_samples_data:
            sheet_title = safe_sheet_name(f"BCR-ABL ({bd['sample']})", used_sheet_names)
            sample_bcr_rows = build_bcr_abl_sheet(
                bd["sample"],
                bd["fusion_rows"],
                bd["mutation_rows"],
                bd["phasing_rows"],
                bd["exon_rows"],
                plot_png=bd.get("plot_png", ""),
            )
            sheets.append((sheet_title, sample_bcr_rows, LAYOUT))

    write_xlsx(output_path, sheets)


    # Also create unique, run-specific named copies so multiple summaries can be opened in Excel simultaneously
    clean_id = re.sub(r"[^\w\.-]", "_", display_id)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")

    # 1. Folder/run-named copy: e.g. reports/dilution_full_run_nagini_summary_report.xlsx
    named_file = os.path.join(out_dir, f"{clean_id}_nagini_summary_report.xlsx")
    if os.path.abspath(named_file) != os.path.abspath(output_path):
        try:
            shutil.copyfile(output_path, named_file)
            sys.stderr.write(f"[Nagini] Created run-named summary workbook: {named_file}\n")
        except Exception as exc:
            sys.stderr.write(f"Warning: could not create run-named summary copy: {exc}\n")

    # 2. Timestamped copy: e.g. reports/dilution_full_run_summary_20260902_1128.xlsx
    time_file = os.path.join(out_dir, f"{clean_id}_summary_{timestamp}.xlsx")
    if os.path.abspath(time_file) != os.path.abspath(output_path) and os.path.abspath(time_file) != os.path.abspath(named_file):
        try:
            shutil.copyfile(output_path, time_file)
            sys.stderr.write(f"[Nagini] Created timestamped summary workbook: {time_file}\n")
        except Exception as exc:
            sys.stderr.write(f"Warning: could not create timestamped summary copy: {exc}\n")



def run_from_snakemake(snakemake):
    sample_names = snakemake.params.sample_names
    sample_types = snakemake.params.sample_types
    barcodes = getattr(snakemake.params, "barcodes", []) or []
    assays = getattr(snakemake.params, "assays", []) or []
    ref_build = getattr(snakemake.params, "reference_build", "hg38")
    seq_id = getattr(snakemake.params, "sequence_id", "Run")

    craminos = list(getattr(snakemake.input, "craminos", []))
    wgs_covs = list(getattr(snakemake.input, "wgs_coverages", []))

    # The optional per-sample files come from params, not input: a Snakemake
    # input list cannot hold an empty entry, so it is compacted to the samples
    # that have the file, and indexing that by sample position attributed one
    # sample's variants to another. The params lists carry one entry per
    # sample with "" where the file does not exist.
    snv_vcfs = list(getattr(snakemake.params, "snv_vcf_paths", []))
    bcr_fusions = list(getattr(snakemake.params, "bcr_abl_fusion_paths", []))
    bcr_mutations = list(getattr(snakemake.params, "bcr_abl_mutation_paths", []))
    bcr_phasings = list(getattr(snakemake.params, "bcr_abl_phasing_paths", []))
    bcr_exons = list(getattr(snakemake.params, "bcr_abl_exon_paths", []))
    bcr_plots = list(getattr(snakemake.params, "bcr_abl_plot_paths", []))

    sample_records = []
    for i, name in enumerate(sample_names):
        stype = sample_types[i] if i < len(sample_types) else "N"
        barcode = barcodes[i] if i < len(barcodes) else ""
        assay = assays[i] if i < len(assays) else "wgs"
        sample_records.append({
            "sample": name,
            "type": stype,
            "barcode": barcode,
            "assay": assay,
            "cramino": craminos[i] if i < len(craminos) else "",
            "wgs_coverage": wgs_covs[i] if i < len(wgs_covs) else "",
            "snv_vcf": snv_vcfs[i] if i < len(snv_vcfs) else "",
            "bcr_abl_fusion": bcr_fusions[i] if i < len(bcr_fusions) else "",
            "bcr_abl_mutations": bcr_mutations[i] if i < len(bcr_mutations) else "",
            "bcr_abl_phasing": bcr_phasings[i] if i < len(bcr_phasings) else "",
            "bcr_abl_exon_coverage": bcr_exons[i] if i < len(bcr_exons) else "",
            "bcr_abl_plot_png": bcr_plots[i] if i < len(bcr_plots) else "",
        })


    build_summary_report(
        output_path=snakemake.output.xlsx,
        sequence_id=seq_id,
        reference_build=ref_build,
        sample_records=sample_records,
    )


if "snakemake" in globals():
    run_from_snakemake(snakemake)
