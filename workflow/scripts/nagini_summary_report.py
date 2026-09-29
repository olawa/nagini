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
        _amplicon_float,
        _amplicon_int,
        _find_metric,
        _is_pathogenic_clinvar,
        build_amplicon_sheet,
        build_bcr_abl_sheet,
        parse_amplicon_stats_tsv,
        parse_amplicon_summary_tsv,
        parse_amplicon_vcf,
        parse_bcr_abl_tsv,
        parse_cramino,
        parse_flt3_itd_tsv,
        derive_flt3_itd_from_variants,
        is_pass_filter,
        itd_is_assessed,
        itd_is_detected,
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
        _amplicon_float,
        _amplicon_int,
        _find_metric,
        _is_pathogenic_clinvar,
        build_amplicon_sheet,
        build_bcr_abl_sheet,
        parse_amplicon_stats_tsv,
        parse_amplicon_summary_tsv,
        parse_amplicon_vcf,
        parse_bcr_abl_tsv,
        parse_cramino,
        parse_flt3_itd_tsv,
        derive_flt3_itd_from_variants,
        is_pass_filter,
        itd_is_assessed,
        itd_is_detected,
        parse_small_variant_vcf,
        parse_tsv,
        safe_sheet_name,
        write_xlsx,
    )


def stack_sample_tables(sample_records, keys):
    """Concatenate one TSV per sample into styled rows, sample column first.

    Column sets can differ between samples and between callers, so the union
    is taken in first-seen order and missing values are left blank rather than
    shifting a row's fields. Returns [] when no sample produced the table, so
    the caller can leave the sheet out entirely.
    """

    tables = []
    for record in sample_records:
        for key in keys:
            path = record.get(key, "")
            if path and os.path.isfile(path):
                columns, rows = parse_tsv(path)
                if rows:
                    tables.append((columns, rows, record))

    if not tables:
        return []

    ordered = []
    for columns, _rows, _record in tables:
        for column in columns:
            if column not in ordered:
                ordered.append(column)
    # The sample is what makes a stacked table readable, so it leads.
    for leading in ("type", "sample"):
        if leading in ordered:
            ordered.insert(0, ordered.pop(ordered.index(leading)))

    # The lab identifier belongs in front of the pipeline's sample id here
    # too: it is what a reader files a finding under.
    has_dna = any(str(record.get("dna_number", "") or "").strip() for _c, _r, record in tables)
    if has_dna:
        ordered.insert(0, "dna_number")

    out = [[(column, 1) for column in ordered]]
    for _columns, rows, record in tables:
        for row in rows:
            values = dict(row)
            values.setdefault("sample", record.get("sample", ""))
            values.setdefault("type", record.get("type", ""))
            if has_dna:
                values["dna_number"] = str(record.get("dna_number", "") or "").strip() or "-"
            out.append([(values.get(column, ""), 0) for column in ordered])
    return out


def build_summary_report(
    output_path: str,
    sequence_id: str,
    reference_build: str,
    sample_records: List[Dict[str, Any]],
    design_bed: str = "",
) -> None:
    """Build the multi-sample workbook."""
    sheets = []
    used_sheet_names = set()

    out_dir = os.path.dirname(output_path) or "."
    analysis_dir_abs = (
        os.path.abspath(os.path.join(out_dir, ".."))
        if out_dir not in (".", "")
        else os.path.abspath(".")
    )
    analysis_folder = os.path.basename(analysis_dir_abs)
    if not analysis_folder or analysis_folder in ("/", "\\", "."):
        analysis_folder = os.path.basename(os.getcwd())
    clean_analysis_name = re.sub(r"[^\w\.-]", "_", analysis_folder)

    seq_clean = str(sequence_id or "").strip()
    display_id = (
        clean_analysis_name
        if (clean_analysis_name and clean_analysis_name.lower() not in {"run", "default", "reports"})
        else (seq_clean or "Run")
    )
    date_str = datetime.now().strftime("%Y-%m-%d")

    # Whether any sample in this run is a BCR-ABL1 workup, decided before the
    # header is built. An adaptive-sampling run carries no fusion, isoform,
    # resistance or phasing data, and five columns of "-" per sample invite
    # the reader to look for a result that was never part of the assay.
    run_has_bcr_abl = any(
        str(record.get("assay", "")).strip().lower() == "bcr_abl"
        or record.get("bcr_abl_fusion")
        or record.get("bcr_abl_mutations")
        or record.get("bcr_abl_phasing")
        or record.get("bcr_abl_exon_coverage")
        for record in sample_records
    )
    run_has_amplicon = any(
        str(record.get("assay", "")).strip().lower() in ("amplicon", "flt3_tp53")
        or record.get("amplicon_stats")
        or record.get("amplicon_summary")
        or record.get("amplicon_itd")
        or record.get("amplicon_vcf")
        for record in sample_records
    )

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
        [
            ("Target Region / Panel", 7),
            (os.path.basename(design_bed) if design_bed else "Whole Genome (WGS)", 0),
            ("Assays", 7),
            (
                ", ".join(
                    sorted(
                        {
                            str(record.get("assay", "") or "").strip()
                            for record in sample_records
                            if str(record.get("assay", "") or "").strip()
                        }
                    )
                )
                or "-",
                0,
            ),
        ],

        [("", 0), ("", 0), ("", 0), ("", 0)],
        [("SAMPLE QUALITY CONTROL & SEQUENCING SUMMARY", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3), ("", 3)],
        [
            ("DNA Number", 1),
            ("Sample ID", 1),
            ("Type", 1),
            ("Barcode", 1),
            ("Assay", 1),
            ("Yield [Gb]", 1),
            ("Total Reads", 1),
            ("Mean Q", 1),
            ("N50 [bp]", 1),
            ("Median [bp]", 1),
            ("Mean Target / KD Depth" if run_has_bcr_abl else "Mean Target Depth", 1),
            ("Mean WGS Depth", 1),
            ("Fold Enrichment", 1),
            ("On-Target Accuracy", 1),
            ("On-Target N50 [bp]", 1),
            ("On-Target Median [bp]", 1),
        ]
        + (
            [
                ("BCR-ABL1 Fusion Reads", 1),
                ("BCR-ABL1 Enrichment", 1),
                ("Primary Isoform", 1),
                ("Resistance Mutations", 1),
                ("Compound Phasing", 1),
            ]
            if run_has_bcr_abl
            else []
        )
        + (
            [
                ("Total Amplicon Reads", 1),
                ("FLT3-ITD Status", 1),
                ("Target Mutations", 1),
            ]
            if run_has_amplicon
            else []
        )
        + [
            ("Status", 1),
        ],
    ]

    has_any_bcr_abl = False
    has_any_amplicon = False
    all_pass_variants = []
    bcr_abl_samples_data = []
    amplicon_samples_data = []

    for s in sample_records:
        sample_name = s["sample"]
        sample_type = s["type"]
        barcode = s.get("barcode", "")
        assay = s.get("assay", "wgs")
        cramino_path = s.get("cramino", "")
        wgs_cov_path = s.get("wgs_coverage", "")
        target_cov_path = s.get("target_coverage", "")
        target_cramino_path = s.get("target_cramino", "")
        dna_number = str(s.get("dna_number", "") or "").strip()
        snv_vcf_path = s.get("snv_vcf", "")
        fusion_tsv_path = s.get("bcr_abl_fusion", "")
        mut_tsv_path = s.get("bcr_abl_mutations", "")
        phase_tsv_path = s.get("bcr_abl_phasing", "")
        exon_cov_path = s.get("bcr_abl_exon_coverage", "")

        amp_stats_path = s.get("amplicon_stats", "")
        amp_sum_path = s.get("amplicon_summary", "")
        amp_itd_path = s.get("amplicon_itd", "")
        amp_vcf_path = s.get("amplicon_vcf", "")
        amp_pileup_path = s.get("amplicon_pileup", "")

        amp_stats_rows = parse_amplicon_stats_tsv(amp_stats_path)
        amp_itd_rows = parse_flt3_itd_tsv(amp_itd_path)
        amp_var_rows = parse_amplicon_vcf(amp_vcf_path)
        if amp_sum_path and (not amp_itd_rows or not amp_var_rows or not amp_stats_rows):
            sec = parse_amplicon_summary_tsv(amp_sum_path)
            if not amp_itd_rows and sec.get("itd"):
                amp_itd_rows = sec["itd"]
            if not amp_var_rows and sec.get("variants"):
                amp_var_rows = sec["variants"]
            if not amp_stats_rows and sec.get("coverage"):
                amp_stats_rows = sec["coverage"]

        is_amplicon = bool(amp_stats_rows or amp_itd_rows or amp_var_rows or assay in ("amplicon", "flt3_tp53"))
        if is_amplicon and not amp_itd_rows:
            amp_itd_rows = derive_flt3_itd_from_variants(
                amp_var_rows, vcf_available=bool(amp_vcf_path and os.path.exists(amp_vcf_path))
            )
        if is_amplicon:
            has_any_amplicon = True
            amplicon_samples_data.append({
                "sample": sample_name,
                "type": sample_type,
                "barcode": barcode,
                "stats_rows": amp_stats_rows,
                "itd_rows": amp_itd_rows,
                "variant_rows": amp_var_rows,
                "pileup_counts_path": amp_pileup_path,
            })

        # Parse cramino metrics
        _, cramino_rows = parse_cramino(cramino_path)
        cramino_dict = {row["Metric"]: row["Value"] for row in cramino_rows}
        yield_gb = _find_metric(cramino_dict, "Yield [Gb]", "Yield", "Bases")
        total_reads_str = _find_metric(
            cramino_dict,
            "Number of reads",
            # cramino labels this "Number of alignments".
            "Number of alignments",
            "Total reads",
            "Read count",
        )
        # cramino reports alignment identity rather than a basecall Q in the
        # versions in use, which is why this column read "-".
        quality_keys = (
            "Mean read quality",
            "Mean quality",
            "Average read quality",
            "Mean base quality",
            "Median read quality",
            "Quality",
            "Q-score",
        )
        # Median before mean: on an adaptive run most reads are the few
        # hundred bases the pore rejected, so the mean describes neither
        # population.
        identity_keys = ("Median identity", "Mean identity", "Identity")
        mean_q = _find_metric(cramino_dict, *quality_keys) or _find_metric(
            cramino_dict, *identity_keys
        )
        n50_bp = _find_metric(cramino_dict, "Read length N50", "N50", "N50 read length")
        median_len = _find_metric(cramino_dict, "Median read length", "Median length", "Median")

        # The same statistics over the target BAM. On an adaptive run the
        # genome-wide numbers describe mostly rejected reads, so accuracy and
        # length on target are the ones that say anything about the calls.
        target_cramino_dict = {}
        if target_cramino_path and os.path.exists(target_cramino_path):
            target_cramino_dict = {
                row["Metric"]: row["Value"] for row in parse_cramino(target_cramino_path)[1]
            }
        on_target_accuracy = _find_metric(
            target_cramino_dict, *identity_keys
        ) or _find_metric(target_cramino_dict, *quality_keys)
        on_target_n50 = _find_metric(
            target_cramino_dict, "Read length N50", "N50", "N50 read length"
        )
        on_target_median = _find_metric(
            target_cramino_dict, "Median read length", "Median length", "Median"
        )


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

        wgs_depth_str = "N/A"
        if wgs_cov_path and os.path.exists(wgs_cov_path):
            _, rows = parse_tsv(wgs_cov_path)
            for row in rows:
                if str(row.get("chrom", "")).strip().lower() == "total":
                    wgs_depth_str = f"{float(row.get('mean', 0.0)):.2f}x"
                    break

        # Target / Kinase coverage
        kd_depth_str = "N/A"
        if exon_rows:
            tot_len = sum(int(r.get("length_bp", r.get("length", 0))) for r in exon_rows)
            if tot_len > 0:
                avg_d = sum(float(r.get("mean_depth", 0.0)) * int(r.get("length_bp", r.get("length", 0))) for r in exon_rows) / tot_len
                kd_depth_str = f"{avg_d:.1f}x"
        elif mutation_rows:
            avg_d = sum(int(r.get("depth", 0)) for r in mutation_rows) / len(mutation_rows)
            kd_depth_str = f"{avg_d:.1f}x"
        elif target_cov_path and os.path.exists(target_cov_path):
            # mosdepth --by writes both "total" (whole reference) and
            # "total_region" (inside the BED). Reading the wrong one reported
            # background depth in a column labelled target depth.
            _, target_rows = parse_tsv(target_cov_path)
            for row in target_rows:
                if str(row.get("chrom", "")).strip().lower() == "total_region":
                    kd_depth_str = f"{float(row.get('mean', 0.0)):.2f}x"
                    break
        else:
            kd_depth_str = wgs_depth_str

        prim_iso = next((r["fusion_isoform"] for r in fusion_rows if r.get("status") == "PRIMARY"), ("None detected" if is_bcr_abl else "-"))
        det_muts = [r["mutation"] for r in mutation_rows if str(r.get("status", "")).upper() == "DETECTED"]
        mut_str = ", ".join(det_muts) if det_muts else ("None detected" if is_bcr_abl else "-")

        cis_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in phasing_rows if "CIS" in r.get("phase_state", "").upper() and "COEXISTING" not in r.get("phase_state", "").upper()]
        coexist_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in phasing_rows if "COEXISTING" in r.get("phase_state", "").upper()]
        trans_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in phasing_rows if "TRANS" in r.get("phase_state", "").upper() and "COEXISTING" not in r.get("phase_state", "").upper()]

        cis_details = []
        for r in phasing_rows:
            ps = r.get("phase_state", "").upper()
            if "CIS" in ps:
                c_rds = int(r.get("cis_reads_11", 0))
                c_span = int(r.get("total_spanning_reads", 0))
                c_af = float(r.get("cis_vaf_percent", (c_rds / c_span * 100.0) if c_span > 0 else 0.0))
                p_name = f"{r.get('mutation_1')}+{r.get('mutation_2')}"
                if "COEXISTING" in ps:
                    cis_details.append(f"{p_name} (coexist, {c_af:.1f}% AF, {c_rds:,} rds)")
                else:
                    cis_details.append(f"{p_name} ({c_af:.1f}% AF, {c_rds:,} rds)")

        phase_str = "-"
        if cis_pairs:
            phase_str = f"CIS ({', '.join(cis_details if cis_details else cis_pairs)})"
        elif coexist_pairs:
            phase_str = f"COEXISTING ({', '.join(cis_details if cis_details else coexist_pairs)})"
        elif trans_pairs:
            phase_str = f"TRANS ({', '.join(trans_pairs)})"
        elif is_bcr_abl:
            phase_str = "N/A (<2 muts)"

        # Target depth alone says nothing about whether the enrichment worked;
        # the ratio to background does.
        fold_enrich_str = "-"
        if target_cov_path or (is_bcr_abl and kd_depth_str != "N/A"):
            try:
                fold = float(kd_depth_str.rstrip("x")) / float(wgs_depth_str.rstrip("x"))
                fold_enrich_str = f"{fold:.1f}x"
            except (ValueError, ZeroDivisionError, AttributeError):
                fold_enrich_str = "-"

        assessed_itd = [r for r in amp_itd_rows if itd_is_assessed(r)]
        det_itd = [r for r in amp_itd_rows if itd_is_detected(r)]
        if det_itd:
            itd_entry = det_itd[0]
            itd_af = _amplicon_float(itd_entry.get("af_pct"))
            itd_len = _amplicon_int(itd_entry.get("dup_len", itd_entry.get("itd_length")))
            ar_val = _amplicon_float(itd_entry.get("ar"))
            ar_part = ", AR >100" if (ar_val >= 999.0 or str(itd_entry.get("ar", "")).startswith(">")) else (f", AR {ar_val:.2f}" if ar_val > 0 else "")
            itd_display = f"DETECTED ({itd_len} bp, {itd_af:.1f}%{ar_part})"
            itd_style = 4
        elif assessed_itd:
            itd_display = "NEGATIVE"
            itd_style = 5
        elif is_amplicon:
            # No ITD result for an amplicon sample is a gap, not a negative.
            itd_display = "NOT ASSESSED"
            itd_style = 6
        else:
            itd_display = "-"
            itd_style = 0

        det_amp_vars = [v for v in amp_var_rows if is_pass_filter(v.get("filter"))]
        if det_amp_vars:
            v_labels = []
            for v in det_amp_vars[:3]:
                g = v.get("gene", "")
                p = v.get("hgvs_p", "") or v.get("hgvs_c", "")
                try:
                    vaf = float(v.get("vaf", 0.0))
                except (ValueError, TypeError):
                    vaf = 0.0
                v_labels.append(f"{g} {p} ({vaf:.1f}%)" if p else f"{g} ({vaf:.1f}%)")
            amp_mut_str = ", ".join(v_labels)
            if len(det_amp_vars) > 3:
                amp_mut_str += f" (+{len(det_amp_vars) - 3})"
            amp_mut_style = 4 if any(_amplicon_float(v.get("vaf")) >= 5.0 for v in det_amp_vars) else 6
        elif is_amplicon:
            amp_mut_str = "None detected"
            amp_mut_style = 5
        else:
            amp_mut_str = "-"
            amp_mut_style = 0

        total_amp_reads_val = sum(_amplicon_int(r.get("read_count")) for r in amp_stats_rows) if amp_stats_rows else ("-" if not is_amplicon else 0)

        sample_status, sample_style = ("PASS", 5)
        if det_muts or cis_pairs or (det_itd and itd_style == 4) or (det_amp_vars and amp_mut_style == 4):
            sample_status, sample_style = ("ACTION REQUIRED", 4)
        elif det_amp_vars or not yield_gb or (raw_reads_num < 100):
            sample_status, sample_style = ("REVIEW", 6)

        overview_rows.append(
            [
                (dna_number or "-", 10),
                (sample_name, 0),
                (sample_type, 0),
                (barcode or "-", 0),
                (assay, 0),
                (yield_gb or "-", 0),
                (total_reads_str or "-", 0),
                (mean_q or "-", 0),
                (n50_bp or "-", 0),
                (median_len or "-", 0),
                (kd_depth_str, 0),
                (wgs_depth_str, 0),
                (fold_enrich_str, 0),
                (on_target_accuracy or "-", 0),
                (on_target_n50 or "-", 0),
                (on_target_median or "-", 0),
            ]
            + (
                [
                    (total_fusion_reads if is_bcr_abl else "-", 0),
                    (target_enrich_str if is_bcr_abl else "-", 0),
                    (prim_iso, 8 if prim_iso not in ("-", "None detected") else 0),
                    (mut_str, 4 if det_muts else 0),
                    (phase_str, 4 if cis_pairs else (6 if trans_pairs else 0)),
                ]
                if run_has_bcr_abl
                else []
            )
            + (
                [
                    (total_amp_reads_val, 0),
                    (itd_display, itd_style),
                    (amp_mut_str, amp_mut_style),
                ]
                if run_has_amplicon
                else []
            )
            + [
                (sample_status, sample_style),
            ]
        )

        # Collect PASS variants for PASS Variants tab
        if mutation_rows:
            for r in mutation_rows:
                if str(r.get("status", "")).upper() == "DETECTED":
                    all_pass_variants.append({
                        "sample": sample_name,
                        "dna_number": dna_number,
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

        if is_amplicon and amp_var_rows:
            for v in amp_var_rows:
                if str(v.get("filter", "PASS")).upper() in ("PASS", "."):
                    try:
                        vaf_val = float(v.get("vaf", 0.0))
                    except (ValueError, TypeError):
                        vaf_val = 0.0
                    all_pass_variants.append({
                        "sample": sample_name,
                        "dna_number": dna_number,
                        "type": sample_type,
                        "gene": v.get("gene", ""),
                        "chrom": v.get("chrom", ""),
                        "pos": v.get("pos", ""),
                        "ref": v.get("ref", ""),
                        "alt": v.get("alt", ""),
                        "vaf": f"{vaf_val:.2f}%",
                        "depth": v.get("depth", ""),
                        "alt_reads": v.get("alt_depth", ""),
                        "consequence": v.get("consequence", ""),
                        "clinical_note": v.get("hgvs_p", "") or v.get("hgvs_c", ""),
                        "status": "DETECTED",
                    })

        if snv_vcf_path and os.path.exists(snv_vcf_path):
            _, snvs, _ = parse_small_variant_vcf(snv_vcf_path)
            for v in snvs:
                if v.get("Filter") in {"PASS", "."}:
                    all_pass_variants.append({
                        "sample": sample_name,
                        "dna_number": dna_number,
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
            ("DNA Number", 1),
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
            (pv.get("dna_number", "") or "-", 10),
            (pv["sample"], 0),
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
    # 2b. THE PER-SAMPLE CALL TABLES, STACKED
    # -----------------------------------------------------------------------
    # The same sheets the per-sample workbooks carry, concatenated across the
    # run. The TSVs already lead with a sample column, so a finding can be read
    # across five samples in one place instead of by opening five files.
    for sheet_name, keys in (
        ("SVs", ("sv_tsv", "rindels_sv_tsv", "savana_sv_tsv")),
        ("Translocations", ("translocations_tsv", "rindels_translocations_tsv", "savana_translocations_tsv")),
        ("CNVs", ("cnv_tsv",)),
        ("Amplicon Stats", ("amplicon_stats",)),
    ):
        stacked = stack_sample_tables(sample_records, keys)
        if stacked:
            sheets.append((sheet_name, stacked, TABLE))

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

            cis_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in p_rows if "CIS" in r.get("phase_state", "").upper() and "COEXISTING" not in r.get("phase_state", "").upper()]
            coexist_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in p_rows if "COEXISTING" in r.get("phase_state", "").upper()]
            trans_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in p_rows if "TRANS" in r.get("phase_state", "").upper() and "COEXISTING" not in r.get("phase_state", "").upper()]
            phase_state = "PASS (No resistance muts)" if not all_active else ("SINGLE MUTATION" if len(all_active) == 1 else ("COMPOUND (CIS)" if cis_pairs else ("COEXISTING (CIS+TRANS)" if coexist_pairs else "POLYCLONAL (TRANS)")))

            clinical_note = "Normal / Sensitive" if not all_active else ("Compound mutation (high TKI resistance)" if (cis_pairs or coexist_pairs) else ("Resistance mutation detected" if det_m else "Subclonal resistance mutation"))

            all_cis_display = []
            for r in p_rows:
                ps = r.get("phase_state", "").upper()
                if "CIS" in ps:
                    c_rds = int(r.get("cis_reads_11", 0))
                    c_span = int(r.get("total_spanning_reads", 0))
                    c_af = float(r.get("cis_vaf_percent", (c_rds / c_span * 100.0) if c_span > 0 else 0.0))
                    p_name = f"{r.get('mutation_1')}+{r.get('mutation_2')}"
                    coexist_tag = " (coexist)" if "COEXISTING" in ps else ""
                    all_cis_display.append(f"{p_name}{coexist_tag}: {c_af:.1f}% ({c_rds:,} rds)")
            bcr_summary_rows.append([
                (bd["sample"], 10),
                (mut_names, 4 if det_m else (6 if subcl_m else 0)),
                (vaf_details, 0),
                (phase_state, 4 if (cis_pairs or coexist_pairs) else (6 if trans_pairs else (4 if det_m else (6 if subcl_m else 5)))),
                (", ".join(all_cis_display) if all_cis_display else "-", 4 if all_cis_display else 0),
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

    # -------------------------------------------------------------------
    # 5. PER-SAMPLE AMPLICON DETAIL SHEETS (if amplicon assay data present)
    # -------------------------------------------------------------------
    if has_any_amplicon:
        for ad in amplicon_samples_data:
            sheet_title = safe_sheet_name(f"Amplicon ({ad['sample']})", used_sheet_names)
            sample_amp_rows = build_amplicon_sheet(
                ad["sample"],
                ad["stats_rows"],
                ad["variant_rows"],
                ad["itd_rows"],
                pileup_counts_path=ad.get("pileup_counts_path", ""),
            )
            sheets.append((sheet_title, sample_amp_rows, LAYOUT))

    write_xlsx(output_path, sheets)


    # Also create unique, analysis-folder named copies so multiple summaries can be opened in Excel simultaneously
    clean_id = clean_analysis_name if (clean_analysis_name and clean_analysis_name.lower() not in {"run", "default", "reports"}) else re.sub(r"[^\w\.-]", "_", display_id)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")

    # Target directories to place named copies in: both the local reports dir and results/reports
    target_dirs = [out_dir]
    results_reports = os.path.join(analysis_dir_abs, "results", "reports")
    if os.path.abspath(results_reports) != os.path.abspath(out_dir):
        target_dirs.append(results_reports)

    for target_dir in target_dirs:
        try:
            os.makedirs(target_dir, exist_ok=True)
            # Create both {analysis_folder}_summary_report.xlsx and {analysis_folder}_nagini_summary_report.xlsx
            for name_pattern in (
                f"{clean_id}_summary_report.xlsx",
                f"{clean_id}_nagini_summary_report.xlsx",
            ):
                named_file = os.path.join(target_dir, name_pattern)
                if os.path.abspath(named_file) != os.path.abspath(output_path):
                    shutil.copyfile(output_path, named_file)
                    sys.stderr.write(f"[Nagini] Created analysis-named summary workbook: {named_file}\n")

            # Timestamped copy: e.g. {analysis_folder}_summary_20260902_1128.xlsx
            time_file = os.path.join(target_dir, f"{clean_id}_summary_{timestamp}.xlsx")
            if os.path.abspath(time_file) != os.path.abspath(output_path):
                shutil.copyfile(output_path, time_file)
        except Exception as exc:
            sys.stderr.write(f"Warning: could not create analysis-named summary copy in {target_dir}: {exc}\n")



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
    target_covs = list(getattr(snakemake.params, "target_coverage_paths", []))
    target_craminos = list(getattr(snakemake.params, "target_cramino_paths", []))
    dna_numbers = list(getattr(snakemake.params, "dna_numbers", []))
    sv_tsvs = list(getattr(snakemake.params, "sv_tsv_paths", []))
    trans_tsvs = list(getattr(snakemake.params, "translocation_tsv_paths", []))
    rindels_sv_tsvs = list(getattr(snakemake.params, "rindels_sv_tsv_paths", []))
    rindels_trans_tsvs = list(getattr(snakemake.params, "rindels_translocation_tsv_paths", []))
    savana_sv_tsvs = list(getattr(snakemake.params, "savana_sv_tsv_paths", []))
    savana_trans_tsvs = list(getattr(snakemake.params, "savana_translocation_tsv_paths", []))
    cnv_tsvs = list(getattr(snakemake.params, "cnv_tsv_paths", []))
    snv_vcfs = list(getattr(snakemake.params, "snv_vcf_paths", []))
    bcr_fusions = list(getattr(snakemake.params, "bcr_abl_fusion_paths", []))
    bcr_mutations = list(getattr(snakemake.params, "bcr_abl_mutation_paths", []))
    bcr_phasings = list(getattr(snakemake.params, "bcr_abl_phasing_paths", []))
    bcr_exons = list(getattr(snakemake.params, "bcr_abl_exon_paths", []))
    bcr_plots = list(getattr(snakemake.params, "bcr_abl_plot_paths", []))
    amplicon_stats = list(getattr(snakemake.params, "amplicon_stats_paths", []))
    amplicon_summaries = list(getattr(snakemake.params, "amplicon_summary_paths", []))
    amplicon_itds = list(getattr(snakemake.params, "amplicon_itd_paths", []))
    amplicon_vcfs = list(getattr(snakemake.params, "amplicon_vcf_paths", []))
    amplicon_pileups = list(getattr(snakemake.params, "amplicon_pileup_paths", []))

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
            "target_coverage": target_covs[i] if i < len(target_covs) else "",
            "target_cramino": target_craminos[i] if i < len(target_craminos) else "",
            "dna_number": dna_numbers[i] if i < len(dna_numbers) else "",
            "sv_tsv": sv_tsvs[i] if i < len(sv_tsvs) else "",
            "translocations_tsv": trans_tsvs[i] if i < len(trans_tsvs) else "",
            "rindels_sv_tsv": rindels_sv_tsvs[i] if i < len(rindels_sv_tsvs) else "",
            "rindels_translocations_tsv": (
                rindels_trans_tsvs[i] if i < len(rindels_trans_tsvs) else ""
            ),
            "savana_sv_tsv": savana_sv_tsvs[i] if i < len(savana_sv_tsvs) else "",
            "savana_translocations_tsv": (
                savana_trans_tsvs[i] if i < len(savana_trans_tsvs) else ""
            ),
            "cnv_tsv": cnv_tsvs[i] if i < len(cnv_tsvs) else "",
            "snv_vcf": snv_vcfs[i] if i < len(snv_vcfs) else "",
            "bcr_abl_fusion": bcr_fusions[i] if i < len(bcr_fusions) else "",
            "bcr_abl_mutations": bcr_mutations[i] if i < len(bcr_mutations) else "",
            "bcr_abl_phasing": bcr_phasings[i] if i < len(bcr_phasings) else "",
            "bcr_abl_exon_coverage": bcr_exons[i] if i < len(bcr_exons) else "",
            "bcr_abl_plot_png": bcr_plots[i] if i < len(bcr_plots) else "",
            "amplicon_stats": amplicon_stats[i] if i < len(amplicon_stats) else "",
            "amplicon_summary": amplicon_summaries[i] if i < len(amplicon_summaries) else "",
            "amplicon_itd": amplicon_itds[i] if i < len(amplicon_itds) else "",
            "amplicon_vcf": amplicon_vcfs[i] if i < len(amplicon_vcfs) else "",
            "amplicon_pileup": amplicon_pileups[i] if i < len(amplicon_pileups) else "",
        })


    build_summary_report(
        output_path=snakemake.output.xlsx,
        sequence_id=seq_id,
        reference_build=ref_build,
        sample_records=sample_records,
        design_bed=str(getattr(snakemake.params, "design_bed", "") or ""),
    )


if "snakemake" in globals():
    run_from_snakemake(snakemake)
