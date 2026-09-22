__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule nagini_excel_report:
    """Build a portable multi-sheet workbook from final QC and variant outputs."""
    input:
        snv_vcf=get_report_snv_vcf,
        snv_tbi=get_report_snv_tbi,
        sv=get_report_sv_tsv,
        translocations=get_report_translocations_tsv,
        cnv_ideogram=get_report_cnv_ideogram,
        rindels_sv=get_report_rindels_sv_tsv,
        rindels_translocations=get_report_rindels_translocations_tsv,
        savana_sv=get_report_savana_sv_tsv,
        savana_translocations=get_report_savana_translocations_tsv,
        cnv=get_report_cnv_tsv,
        wgs_coverage="qc/mosdepth/{sample}_{type}.mosdepth.summary.txt",
        target_coverage=get_report_target_coverage,
        target_dist=get_report_target_dist,
        demux_summary=get_report_demux_summary,
        cramino="qc/cramino/{sample}_{type}.txt",
        target_cramino=get_report_target_cramino,
        bcr_abl_fusion=get_report_bcr_abl_fusion,
        bcr_abl_mutations=get_report_bcr_abl_mutations,
        bcr_abl_phasing=get_report_bcr_abl_phasing,
        bcr_abl_exon_coverage=get_report_bcr_abl_exon_coverage,
        bcr_abl_plot_png=get_report_bcr_abl_plot_png,
        amplicon_stats=get_report_amplicon_stats,
        amplicon_summary=get_report_amplicon_summary,
        amplicon_itd=get_report_amplicon_itd,
        amplicon_vcf=get_report_amplicon_vcf,
        amplicon_pileup=get_report_amplicon_pileup,
    output:
        xlsx="reports/{sample}_{type}.nagini_report.xlsx",
    params:
        reference_build=lambda wildcards: config.get("reference", {}).get("build", config.get("reference_build", "hg38")),
        snv_annotation=(
            "VEP" if vep_enabled() else "none"
        ),
        cytoband=get_cytoband_file,
        small_variant_caller=get_small_variant_caller(),
        design_bed=lambda wildcards: get_sample_design_bed(wildcards.sample),
        clairs_to_platform=lambda wildcards: config.get("clairs_to_call", {}).get("platform", ""),
        barcode=get_sample_barcode_name,
        dna_number=lambda wildcards: get_sample_dna_number(wildcards.sample),
    log:
        "logs/nagini_excel_report/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/nagini_excel_report/{sample}_{type}.tsv",
            config.get("nagini_excel_report", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("nagini_excel_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("nagini_excel_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("nagini_excel_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("nagini_excel_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("nagini_excel_report", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: compile Excel report for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/nagini_report.py"


# Keys whose file exists only for some assays. Everything that reads them has
# to align them to the sample list by position, so they are built once, per
# sample, with an empty string where the file does not exist.
_SUMMARY_OPTIONAL_KEYS = (
    "target_coverages",
    "target_craminos",
    "sv_tsvs",
    "translocation_tsvs",
    "rindels_sv_tsvs",
    "rindels_translocation_tsvs",
    "savana_sv_tsvs",
    "savana_translocation_tsvs",
    "cnv_tsvs",
    "snv_vcfs",
    "bcr_abl_fusions",
    "bcr_abl_mutations",
    "bcr_abl_phasings",
    "bcr_abl_exons",
    "bcr_abl_plots",
    "amplicon_stats",
    "amplicon_summaries",
    "amplicon_itds",
    "amplicon_vcfs",
    "amplicon_pileups",
)


def summary_report_files():
    """One row per sample/type with every summary input, "" when absent.

    The optional files exist only for some assays. Appending them only for
    qualifying samples and then indexing the compacted list by sample position
    attributed one sample's variants to another as soon as a run mixed assays,
    so the alignment is built here and both the rule's input and its params
    are derived from it.
    """

    rows = []
    for sample in get_samples(samples):
        s = str(sample).strip()
        for unit_type in get_unit_types_for_sample(sample):
            t = str(unit_type).strip()
            snv_vcf = ""
            if sample_needs_dna_variant_calling(s):
                caller = get_small_variant_caller(s, t)
                if caller == "deepvariant":
                    caller_dir = "snv_indels/deepvariant_postprocess"
                elif caller == "mutect2":
                    caller_dir = "snv_indels/mutect2_postprocess"
                else:
                    caller_dir = "snv_indels/clairs_to_postprocess"
                suffix = ".filtered.norm.vep.vcf.gz" if vep_enabled() else ".filtered.norm.vcf.gz"
                snv_vcf = f"{caller_dir}/{s}_{t}{suffix}"
            bcr = sample_has_assay_module(s, "bcr_abl")
            amp = sample_has_assay_module(s, "amplicon")
            rows.append(
                {
                    "reports": f"reports/{s}_{t}.nagini_report.xlsx",
                    "craminos": f"qc/cramino/{s}_{t}.txt",
                    "wgs_coverages": f"qc/mosdepth/{s}_{t}.mosdepth.summary.txt",
                    "target_craminos": (
                        f"qc/cramino_target/{s}_{t}.txt"
                        if sample_has_target_metrics(s)
                        else ""
                    ),
                    "target_coverages": (
                        f"qc/mosdepth_bed/{s}_{t}.mosdepth.summary.txt"
                        if sample_has_target_metrics(s)
                        else ""
                    ),
                    "snv_vcfs": snv_vcf,
                    "sv_tsvs": (
                        f"reports/sv_summary/{s}_{t}.tsv"
                        if sniffles2_enabled(s)
                        else ""
                    ),
                    "translocation_tsvs": (
                        f"reports/sv_summary/{s}_{t}.translocations.tsv"
                        if sniffles2_enabled(s)
                        else ""
                    ),
                    "rindels_sv_tsvs": (
                        f"reports/sv_summary_rindels/{s}_{t}.tsv"
                        if sample_has_rindels_sv(s)
                        else ""
                    ),
                    "rindels_translocation_tsvs": (
                        f"reports/sv_summary_rindels/{s}_{t}.translocations.tsv"
                        if sample_has_rindels_sv(s)
                        else ""
                    ),
                    "savana_sv_tsvs": (
                        f"reports/sv_summary_savana/{s}_{t}.tsv"
                        if savana_enabled(s)
                        else ""
                    ),
                    "savana_translocation_tsvs": (
                        f"reports/sv_summary_savana/{s}_{t}.translocations.tsv"
                        if savana_enabled(s)
                        else ""
                    ),
                    "cnv_tsvs": (
                        f"reports/cnv_summary/{s}_{t}.cnv.tsv"
                        if sample_needs_cnv_calling(s)
                        else ""
                    ),
                    "bcr_abl_fusions": f"reports/bcr_abl/{s}_{t}.fusion_summary.tsv" if bcr else "",
                    "bcr_abl_mutations": f"reports/bcr_abl/{s}_{t}.branford_mutations.tsv" if bcr else "",
                    "bcr_abl_phasings": f"reports/bcr_abl/{s}_{t}.phasing.tsv" if bcr else "",
                    "bcr_abl_exons": f"reports/bcr_abl/{s}_{t}.exon_coverage.tsv" if bcr else "",
                    "bcr_abl_plots": f"reports/bcr_abl/{s}_{t}.fusion_plot.png" if bcr else "",
                    "amplicon_stats": f"reports/amplicon/{s}_{t}.amplicon_stats.tsv" if amp else "",
                    "amplicon_summaries": f"reports/amplicon/{s}_{t}.amplicon_summary.tsv" if amp else "",
                    "amplicon_itds": f"reports/amplicon/{s}_{t}.flt3_itd.tsv" if amp else "",
                    "amplicon_vcfs": f"snv_indels/rindels/{s}_{t}.amplicon.vcf.gz" if amp else "",
                    "amplicon_pileups": f"reports/amplicon/{s}_{t}.pileup_counts.tsv" if amp else "",
                }
            )
    return rows


def get_summary_report_inputs(wildcards):
    """The files the summary depends on. Snakemake inputs cannot be empty."""

    rows = summary_report_files()
    return {
        key: [row[key] for row in rows if row[key]]
        for key in rows[0]
    } if rows else {}


def get_summary_report_paths(key):
    """One entry per sample for an optional input, "" where it does not exist."""

    return [row[key] for row in summary_report_files()]


rule nagini_summary_excel_report:
    """Build a multi-sample summary Excel workbook across all samples in the run."""
    input:
        unpack(get_summary_report_inputs),
    output:
        xlsx="reports/nagini_summary_report.xlsx",
    params:
        sample_names=lambda wildcards: [
            str(s).strip()
            for s in get_samples(samples)
            for _ in get_unit_types_for_sample(s)
        ],
        sample_types=lambda wildcards: [
            str(t).strip()
            for s in get_samples(samples)
            for t in get_unit_types_for_sample(s)
        ],
        dna_numbers=lambda wildcards: [
            get_sample_dna_number(str(s).strip())
            for s in get_samples(samples)
            for _ in get_unit_types_for_sample(s)
        ],
        barcodes=lambda wildcards: [
            get_sample_barcode_name_str(str(s).strip(), str(t).strip())
            for s in get_samples(samples)
            for t in get_unit_types_for_sample(s)
        ],
        assays=lambda wildcards: [
            get_assay(str(s).strip())
            for s in get_samples(samples)
            for _ in get_unit_types_for_sample(s)
        ],
        reference_build=lambda wildcards: config.get("reference", {}).get("build", config.get("reference_build", "hg38")),
        sequence_id=lambda wildcards: config.get("sequence_id", os.path.basename(os.getcwd())),
        # Aligned to sample_names, unlike input, which cannot carry a gap.
        target_coverage_paths=lambda wildcards: get_summary_report_paths("target_coverages"),
        target_cramino_paths=lambda wildcards: get_summary_report_paths("target_craminos"),
        sv_tsv_paths=lambda wildcards: get_summary_report_paths("sv_tsvs"),
        translocation_tsv_paths=lambda wildcards: get_summary_report_paths("translocation_tsvs"),
        rindels_sv_tsv_paths=lambda wildcards: get_summary_report_paths("rindels_sv_tsvs"),
        rindels_translocation_tsv_paths=lambda wildcards: get_summary_report_paths(
            "rindels_translocation_tsvs"
        ),
        savana_sv_tsv_paths=lambda wildcards: get_summary_report_paths("savana_sv_tsvs"),
        savana_translocation_tsv_paths=lambda wildcards: get_summary_report_paths(
            "savana_translocation_tsvs"
        ),
        cnv_tsv_paths=lambda wildcards: get_summary_report_paths("cnv_tsvs"),
        design_bed=lambda wildcards: get_resolved_design_bed(),
        snv_vcf_paths=lambda wildcards: get_summary_report_paths("snv_vcfs"),
        bcr_abl_fusion_paths=lambda wildcards: get_summary_report_paths("bcr_abl_fusions"),
        bcr_abl_mutation_paths=lambda wildcards: get_summary_report_paths("bcr_abl_mutations"),
        bcr_abl_phasing_paths=lambda wildcards: get_summary_report_paths("bcr_abl_phasings"),
        bcr_abl_exon_paths=lambda wildcards: get_summary_report_paths("bcr_abl_exons"),
        bcr_abl_plot_paths=lambda wildcards: get_summary_report_paths("bcr_abl_plots"),
        amplicon_stats_paths=lambda wildcards: get_summary_report_paths("amplicon_stats"),
        amplicon_summary_paths=lambda wildcards: get_summary_report_paths("amplicon_summaries"),
        amplicon_itd_paths=lambda wildcards: get_summary_report_paths("amplicon_itds"),
        amplicon_vcf_paths=lambda wildcards: get_summary_report_paths("amplicon_vcfs"),
        amplicon_pileup_paths=lambda wildcards: get_summary_report_paths("amplicon_pileups"),
    log:
        "logs/nagini_summary_excel_report/summary.log",
    threads: 1
    resources:
        mem_mb=config.get("nagini_excel_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        threads=1,
    container:
        config.get("nagini_excel_report", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: compile multi-sample Excel summary report"
    script:
        "../scripts/nagini_summary_report.py"


