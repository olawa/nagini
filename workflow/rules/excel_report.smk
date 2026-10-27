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
        bcr_abl_all_sites=get_report_bcr_abl_all_sites,
        bcr_abl_phasing=get_report_bcr_abl_phasing,
        bcr_abl_exon_coverage=get_report_bcr_abl_exon_coverage,
        bcr_abl_plot_png=get_report_bcr_abl_plot_png,
        bcr_abl_vardictrust_vcf=get_report_bcr_abl_vardictrust_vcf,
        bcr_abl_rindels_vcf=get_report_bcr_abl_rindels_vcf,
        bcr_abl_bcftools_vcf=get_report_bcr_abl_bcftools_vcf,
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
        bcr_abl_min_relative_depth=lambda wildcards: config.get("bcr_abl", {}).get("min_relative_depth", 0.10),
        bcr_abl_min_depth=lambda wildcards: config.get("bcr_abl", {}).get("min_depth", 50.0),
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


rule nagini_summary_excel_report:
    """Build a multi-sample summary Excel workbook across all samples in the run."""
    input:
        unpack(get_summary_report_inputs),
    output:
        xlsx="reports/nagini_summary_report.xlsx",
    params:
        sample_names=lambda wildcards: [
            clean_report_path(s)
            for s in get_samples(samples)
            for _ in get_unit_types_for_sample(s)
        ],
        sample_types=lambda wildcards: [
            clean_report_path(t)
            for s in get_samples(samples)
            for t in get_unit_types_for_sample(s)
        ],
        dna_numbers=lambda wildcards: [
            clean_report_path(get_sample_dna_number(clean_report_path(s)))
            for s in get_samples(samples)
            for _ in get_unit_types_for_sample(s)
        ],
        barcodes=lambda wildcards: [
            clean_report_path(get_sample_barcode_name_str(clean_report_path(s), clean_report_path(t)))
            for s in get_samples(samples)
            for t in get_unit_types_for_sample(s)
        ],
        assays=lambda wildcards: [
            clean_report_path(get_assay(clean_report_path(s)))
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
        bcr_abl_all_sites_paths=lambda wildcards: get_summary_report_paths("bcr_abl_all_sites"),
        bcr_abl_phasing_paths=lambda wildcards: get_summary_report_paths("bcr_abl_phasings"),
        bcr_abl_exon_paths=lambda wildcards: get_summary_report_paths("bcr_abl_exons"),
        bcr_abl_plot_paths=lambda wildcards: get_summary_report_paths("bcr_abl_plots"),
        bcr_abl_vardictrust_vcf_paths=lambda wildcards: get_summary_report_paths("bcr_abl_vardictrust_vcfs"),
        bcr_abl_rindels_vcf_paths=lambda wildcards: get_summary_report_paths("bcr_abl_rindels_vcfs"),
        bcr_abl_bcftools_vcf_paths=lambda wildcards: get_summary_report_paths("bcr_abl_bcftools_vcfs"),
        bcr_abl_min_relative_depth=lambda wildcards: config.get("bcr_abl", {}).get("min_relative_depth", 0.10),
        bcr_abl_min_depth=lambda wildcards: config.get("bcr_abl", {}).get("min_depth", 50.0),
    log:
        "logs/nagini_summary_excel_report/summary.log",
    benchmark:
        repeat(
            "benchmarks/nagini_summary_excel_report/nagini_summary_excel_report.tsv",
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
        "{rule}: compile multi-sample Excel summary report"
    script:
        "../scripts/nagini_summary_report.py"


