"""Create a multi-sample XLSX summary report for a Nagini run."""

import argparse
import os
import re
import shutil
import sys

try:
    from bcr_abl_report_utils import (
        active_resistance_mutations,
        calculate_bcr_abl_depth_cutoff,
        format_vaf,
        is_resistance_associated,
        panel_assessment,
        summarize_phasing,
    )
except ImportError:
    from .bcr_abl_report_utils import (
        active_resistance_mutations,
        calculate_bcr_abl_depth_cutoff,
        format_vaf,
        is_resistance_associated,
        panel_assessment,
        summarize_phasing,
    )
from datetime import datetime, timezone

from pathlib import Path
from typing import Any, Dict, List, Optional

try:
    import xlsxwriter
except ImportError:
    xlsxwriter = None

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
    bcr_abl_min_relative_depth: float = 0.10,
    bcr_abl_min_depth: float = 50.0,
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
        or record.get("bcr_abl_all_sites")
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
                ("Mean ABL1 Panel-site Depth" if run_has_bcr_abl else "Mean Target Depth", 1),
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
        + [
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
        target_cov_path = s.get("target_coverage", "")
        target_cramino_path = s.get("target_cramino", "")
        dna_number = str(s.get("dna_number", "") or "").strip()
        snv_vcf_path = s.get("snv_vcf", "")
        fusion_tsv_path = s.get("bcr_abl_fusion", "")
        mut_tsv_path = s.get("bcr_abl_mutations", "")
        phase_tsv_path = s.get("bcr_abl_phasing", "")
        exon_cov_path = s.get("bcr_abl_exon_coverage", "")
        all_sites_path = s.get("bcr_abl_all_sites", "")

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
        all_sites_rows = parse_bcr_abl_tsv(all_sites_path)
        is_bcr_abl = bool(fusion_rows or mutation_rows or phasing_rows or exon_rows or all_sites_rows or assay == "bcr_abl")

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
                "all_sites_rows": all_sites_rows,
                "plot_png": s.get("bcr_abl_plot_png", ""),
            })

            # Check if read_stats exist from BCR-ABL1 fusion summary JSON or TSV

            bcr_read_stats = {}
            if fusion_tsv_path:
                fusion_json_path = fusion_tsv_path.rsplit(".", 1)[0] + ".json"
                if os.path.isfile(fusion_json_path):
                    try:
                        with open(fusion_json_path, "r", encoding="utf-8") as f:
                            bcr_read_stats = json.load(f).get("read_stats", {})
                    except Exception:
                        bcr_read_stats = {}
                if not bcr_read_stats and os.path.isfile(fusion_tsv_path):
                    try:
                        with open(fusion_tsv_path, "r", encoding="utf-8") as f:
                            for line in f:
                                if not line.startswith("#"):
                                    break
                                if "=" in line:
                                    k, v = line[1:].strip().split("=", 1)
                                    bcr_read_stats[k.strip()] = v.strip()
                    except Exception:
                        pass

            if bcr_read_stats:
                if "yield_gb" in bcr_read_stats and bcr_read_stats["yield_gb"] is not None:
                    yield_gb = f"{float(bcr_read_stats['yield_gb']):.2f}"
                if "total_reads" in bcr_read_stats and bcr_read_stats["total_reads"] is not None:
                    total_reads_str = f"{int(float(bcr_read_stats['total_reads'])):,}"
                if "mean_q" in bcr_read_stats and bcr_read_stats["mean_q"] is not None:
                    mean_q = f"{float(bcr_read_stats['mean_q']):.2f}"
                if "n50_bp" in bcr_read_stats and bcr_read_stats["n50_bp"] is not None:
                    n50_bp = str(int(float(bcr_read_stats["n50_bp"])))
                if "median_bp" in bcr_read_stats and bcr_read_stats["median_bp"] is not None:
                    median_len = f"{float(bcr_read_stats['median_bp']):.2f}"
                if "on_target_accuracy" in bcr_read_stats and bcr_read_stats["on_target_accuracy"] is not None:
                    on_target_accuracy = f"{float(bcr_read_stats['on_target_accuracy']):.2f}"
                if "on_target_n50_bp" in bcr_read_stats and bcr_read_stats["on_target_n50_bp"] is not None:
                    on_target_n50 = str(int(float(bcr_read_stats["on_target_n50_bp"])))
                if "on_target_median_bp" in bcr_read_stats and bcr_read_stats["on_target_median_bp"] is not None:
                    on_target_median = f"{float(bcr_read_stats['on_target_median_bp']):.2f}"


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
        mean_target_dp = 0.0
        if all_sites_rows:
            site_depths = {}
            for row in all_sites_rows:
                pos = str(row.get("pos", row.get("position", "")))
                if pos:
                    try:
                        site_depths.setdefault(f"{row.get('chrom', '')}:{pos}", float(row.get("depth", row.get("dp", 0)) or 0))
                    except (ValueError, TypeError):
                        pass
            if site_depths:
                mean_target_dp = sum(site_depths.values()) / len(site_depths)
                kd_depth_str = f"{mean_target_dp:.1f}x"
        elif exon_rows:
            abl_exons = [r for r in exon_rows if "ABL1" in str(r.get("exon", "")).upper()]
            tot_len = sum(int(r.get("length_bp", r.get("length", 0))) for r in abl_exons)
            if tot_len > 0:
                mean_target_dp = sum(float(r.get("mean_depth", 0.0)) * int(r.get("length_bp", r.get("length", 0))) for r in abl_exons) / tot_len
                kd_depth_str = f"{mean_target_dp:.1f}x"
        elif not is_bcr_abl and target_cov_path and os.path.exists(target_cov_path):
            # mosdepth --by writes both "total" (whole reference) and
            # "total_region" (inside the BED). Reading the wrong one reported
            # background depth in a column labelled target depth.
            _, target_rows = parse_tsv(target_cov_path)
            for row in target_rows:
                if str(row.get("chrom", "")).strip().lower() == "total_region":
                    kd_depth_str = f"{float(row.get('mean', 0.0)):.2f}x"
                    break
        elif is_bcr_abl:
            kd_depth_str = "N/A"
        else:
            kd_depth_str = wgs_depth_str

        prim_iso = next((r["fusion_isoform"] for r in fusion_rows if r.get("status") == "PRIMARY"), ("None detected" if is_bcr_abl else "-"))
        active_muts = active_resistance_mutations(mutation_rows)
        det_muts = [r["mutation"] for r in active_muts if str(r.get("status", "")).upper() == "DETECTED"]
        subcl_muts = [r["mutation"] for r in active_muts if str(r.get("status", "")).upper() == "SUBCLONAL"]
        unclassified = [r for r in mutation_rows if str(r.get("status", "")).upper() in {"DETECTED", "SUBCLONAL"} and not is_resistance_associated(r)]
        assessment = panel_assessment(all_sites_rows)
        mut_str = ", ".join(det_muts + [f"{m} (SUBCLONAL)" for m in subcl_muts]) if det_muts or subcl_muts else ("Resistance panel not assessed" if assessment["state"] == "NOT_ASSESSED" else "No resistance-associated variants detected")
        if unclassified:
            mut_str += f"; unclassified for review: {', '.join(r.get('mutation', '?') for r in unclassified)}"

        phase_str = summarize_phasing(phasing_rows) if phasing_rows else "NOT_ASSESSED"
        cis_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in phasing_rows if str(r.get("phase_state", "")).strip().upper() == "CIS (COMPOUND)"]
        coexist_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in phasing_rows if str(r.get("phase_state", "")).strip().upper() == "COEXISTING (COMPOUND + POLYCLONAL)"]
        trans_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in phasing_rows if str(r.get("phase_state", "")).strip().upper() == "TRANS (POLYCLONAL)"]

        # Target depth alone says nothing about whether the enrichment worked;
        # the ratio to background does.
        fold_enrich_str = "-"
        if target_cov_path or (is_bcr_abl and kd_depth_str != "N/A"):
            try:
                fold = float(kd_depth_str.rstrip("x")) / float(wgs_depth_str.rstrip("x"))
                fold_enrich_str = f"{fold:.1f}x"
            except (ValueError, ZeroDivisionError, AttributeError):
                fold_enrich_str = "-"

        sample_status, sample_style = ("PASS", 5)
        if det_muts or cis_pairs:
            sample_status, sample_style = ("ACTION REQUIRED", 4)
        elif not yield_gb or (raw_reads_num < 100):
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
                        "vaf": format_vaf(r.get("vaf_percent", ""), scale="percent"),
                        "depth": r.get("depth", ""),
                        "alt_reads": r.get("ad_alt", r.get("alt_reads", "")),
                        "consequence": r.get("mutation", ""),
                        "clinical_note": "",  # no ClinVar evidence for a Branford-caller row
                        "status": "DETECTED",
                        "caller": "Branford",
                    })

        min_depth_cutoff = 0.0
        if is_bcr_abl and mean_target_dp > 0:
            min_depth_cutoff = calculate_bcr_abl_depth_cutoff(
                mean_target_dp, bcr_abl_min_relative_depth, bcr_abl_min_depth
            )

        if snv_vcf_path and os.path.exists(snv_vcf_path):
            _, snvs, caller_name = parse_small_variant_vcf(snv_vcf_path)
            for v in snvs:
                if v.get("Filter") in {"PASS", "."}:
                    try:
                        v_dp = float(str(v.get("Depth", 0)).split(",")[0])
                    except (ValueError, TypeError):
                        v_dp = 0.0
                    if is_bcr_abl and min_depth_cutoff > 0 and v_dp < min_depth_cutoff:
                        continue
                    all_pass_variants.append({
                        "sample": sample_name,
                        "dna_number": dna_number,
                        "type": sample_type,
                        "gene": v.get("Gene", ""),
                        "chrom": v.get("Chromosome", ""),
                        "pos": v.get("Position", ""),
                        "ref": v.get("Reference", ""),
                        "alt": v.get("Alternate", ""),
                        "vaf": format_vaf(v.get("VAF", ""), scale="fraction"),
                        "depth": v.get("Depth", ""),
                        "alt_reads": v.get("Allele depth", v.get("Alt_Reads", "")),
                        "consequence": v.get("Consequence", ""),
                        "clinical_note": v.get("ClinVar_Significance", ""),
                        "status": "PASS",
                        "caller": v.get("Caller") or caller_name or "ClairS-TO",
                    })

        for caller_key, default_caller_name in (
            ("bcr_abl_vardictrust_vcf", "VarDictRust"),
            ("bcr_abl_rindels_vcf", "Rindels"),
            ("bcr_abl_bcftools_vcf", "bcftools"),
        ):
            c_vcf = s.get(caller_key, "")
            if c_vcf and c_vcf.endswith(".gz") and (not os.path.isfile(c_vcf) or os.path.getsize(c_vcf) == 0):
                uncompressed = c_vcf[:-3]
                if os.path.isfile(uncompressed) and os.path.getsize(uncompressed) > 0:
                    c_vcf = uncompressed
            if c_vcf and os.path.exists(c_vcf):
                _, c_snvs, c_name = parse_small_variant_vcf(c_vcf, default_caller=default_caller_name)
                for v in c_snvs:
                    if v.get("Filter") in {"PASS", "."}:
                        try:
                            v_dp = float(str(v.get("Depth", 0)).split(",")[0])
                        except (ValueError, TypeError):
                            v_dp = 0.0
                        if is_bcr_abl and min_depth_cutoff > 0 and v_dp < min_depth_cutoff:
                            continue
                        all_pass_variants.append({
                            "sample": sample_name,
                            "dna_number": dna_number,
                            "type": sample_type,
                            "gene": v.get("Gene", "") or ("ABL1" if str(v.get("Chromosome", "")).removeprefix("chr").lower() in ("9", "abl1") else ("BCR" if str(v.get("Chromosome", "")).removeprefix("chr").lower() in ("22", "bcr") else "")),
                            "chrom": v.get("Chromosome", ""),
                            "pos": v.get("Position", ""),
                            "ref": v.get("Reference", ""),
                            "alt": v.get("Alternate", ""),
                            "vaf": format_vaf(v.get("VAF", ""), scale="fraction"),
                            "depth": v.get("Depth", ""),
                            "alt_reads": v.get("Allele depth", v.get("Alt_Reads", "")),
                            "consequence": v.get("Consequence", ""),
                            "clinical_note": v.get("ClinVar_Significance", "") or f"{default_caller_name} call",
                            "status": "PASS",
                            "caller": v.get("Caller") or c_name or default_caller_name,
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
            ("Caller", 1),
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
            (pv.get("caller", "ClairS-TO"), 0),
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
            record.get("snv_vcf")
            or record.get("bcr_abl_mutations")
            or record.get("bcr_abl_vardictrust_vcf")
            or record.get("bcr_abl_rindels_vcf")
            or record.get("bcr_abl_bcftools_vcf")
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
            active = active_resistance_mutations(m_rows)
            det_m = [r for r in active if str(r.get("status", "")).upper() == "DETECTED"]
            subcl_m = [r for r in active if str(r.get("status", "")).upper() == "SUBCLONAL"]
            unclassified = [r for r in m_rows if str(r.get("status", "")).upper() in {"DETECTED", "SUBCLONAL"} and not is_resistance_associated(r)]
            assessment = panel_assessment(bd.get("all_sites_rows", []))
            if det_m:
                mut_names = ", ".join(r["mutation"] for r in det_m)
                if subcl_m:
                    mut_names += f" (Subclonal: {', '.join(r['mutation'] for r in subcl_m)})"
            elif subcl_m:
                mut_names = f"Subclonal: {', '.join(r['mutation'] for r in subcl_m)}"
            elif assessment["state"] == "ASSESSED":
                mut_names = "No resistance-associated variants detected (panel assessed)"
            elif assessment["state"] == "INCOMPLETE":
                mut_names = f"Panel incomplete ({assessment['callable_positions']}/{assessment['positions']} sites callable)"
            elif assessment["state"] == "COMPLETENESS_UNVERIFIED":
                mut_names = "No resistance-associated variants detected; completeness unverified"
            else:
                mut_names = "Resistance panel not assessed"
            if unclassified:
                mut_names += f"; unclassified for review: {', '.join(r.get('mutation', '?') for r in unclassified)}"

            all_active = det_m + subcl_m
            vaf_details = ", ".join(f"{r['mutation']}: {float(r.get('vaf_percent', r.get('vaf', 0.0))):.1f}% ({r.get('ad_alt', r.get('alt_reads', 0))}/{r.get('depth')})" for r in all_active) if all_active else "-"
            phase_state = summarize_phasing(p_rows)
            phase_classes = {str(r.get("phase_state", "")).strip().upper() for r in p_rows}
            cis_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in p_rows if str(r.get("phase_state", "")).strip().upper() == "CIS (COMPOUND)"]
            coexist_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in p_rows if str(r.get("phase_state", "")).strip().upper() == "COEXISTING (COMPOUND + POLYCLONAL)"]
            trans_pairs = [f"{r['mutation_1']}+{r['mutation_2']}" for r in p_rows if str(r.get("phase_state", "")).strip().upper() == "TRANS (POLYCLONAL)"]
            clinical_note = "Resistance-associated variant detected" if det_m else ("Subclonal resistance-associated variant detected" if subcl_m else ("Resistance panel incomplete" if assessment["state"] == "INCOMPLETE" else ("Resistance panel not assessed" if assessment["state"] == "NOT_ASSESSED" else "No resistance-associated variant detected; sensitivity not inferred")))
            if unclassified:
                clinical_note += "; unclassified variants require review"

            all_cis_display = []
            for r in p_rows:
                ps = r.get("phase_state", "").upper()
                if ps in {"CIS (COMPOUND)", "COEXISTING (COMPOUND + POLYCLONAL)"}:
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
                (phase_state, 4 if (cis_pairs or coexist_pairs) else (6 if p_rows else 0)),
                (", ".join(all_cis_display) if all_cis_display else "-", 4 if all_cis_display else 0),
                (", ".join(trans_pairs) if trans_pairs else "-", 6 if trans_pairs else 0),
                (clinical_note, 4 if det_m else (6 if subcl_m or unclassified or assessment["state"] != "ASSESSED" else 0)),
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
                all_sites_rows=bd.get("all_sites_rows", []),
            )
            sheets.append((sheet_title, sample_bcr_rows, LAYOUT))

    write_xlsx(output_path, sheets)


    # Also create unique, analysis-folder named copy so multiple summaries can be opened in Excel simultaneously
    clean_id = clean_analysis_name if (clean_analysis_name and clean_analysis_name.lower() not in {"run", "default", "reports"}) else re.sub(r"[^\w\.-]", "_", display_id)

    try:
        named_file = os.path.join(out_dir, f"{clean_id}_summary_report.xlsx")
        if os.path.abspath(named_file) != os.path.abspath(output_path):
            shutil.copyfile(output_path, named_file)
            sys.stderr.write(f"[Nagini] Created analysis-named summary workbook: {named_file}\n")
    except Exception as exc:
        sys.stderr.write(f"Warning: could not create analysis-named summary copy in {out_dir}: {exc}\n")



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
    bcr_all_sites = list(getattr(snakemake.params, "bcr_abl_all_sites_paths", []))
    bcr_plots = list(getattr(snakemake.params, "bcr_abl_plot_paths", []))
    bcr_vardict_vcfs = list(getattr(snakemake.params, "bcr_abl_vardictrust_vcf_paths", []))
    bcr_rindels_vcfs = list(getattr(snakemake.params, "bcr_abl_rindels_vcf_paths", []))
    bcr_bcftools_vcfs = list(getattr(snakemake.params, "bcr_abl_bcftools_vcf_paths", []))

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
            "bcr_abl_all_sites": bcr_all_sites[i] if i < len(bcr_all_sites) else "",
            "bcr_abl_plot_png": bcr_plots[i] if i < len(bcr_plots) else "",
            "bcr_abl_vardictrust_vcf": bcr_vardict_vcfs[i] if i < len(bcr_vardict_vcfs) else "",
            "bcr_abl_rindels_vcf": bcr_rindels_vcfs[i] if i < len(bcr_rindels_vcfs) else "",
            "bcr_abl_bcftools_vcf": bcr_bcftools_vcfs[i] if i < len(bcr_bcftools_vcfs) else "",
        })


    bcr_abl_min_relative_depth = float(
        getattr(snakemake.params, "bcr_abl_min_relative_depth", 0.10) or 0.10
    )
    bcr_abl_min_depth = float(
        getattr(snakemake.params, "bcr_abl_min_depth", 50.0) or 50.0
    )

    build_summary_report(
        output_path=snakemake.output.xlsx,
        sequence_id=seq_id,
        reference_build=ref_build,
        sample_records=sample_records,
        design_bed=str(getattr(snakemake.params, "design_bed", "") or ""),
        bcr_abl_min_relative_depth=bcr_abl_min_relative_depth,
        bcr_abl_min_depth=bcr_abl_min_depth,
    )


if "snakemake" in globals():
    run_from_snakemake(snakemake)
