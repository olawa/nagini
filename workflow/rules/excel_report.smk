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
        wgs_coverage="qc/mosdepth/{sample}_{type}.mosdepth.summary.txt",
        target_coverage=get_report_target_coverage,
        target_dist=get_report_target_dist,
        demux_summary=get_report_demux_summary,
        cramino="qc/cramino/{sample}_{type}.txt",
        bcr_abl_fusion=get_report_bcr_abl_fusion,
        bcr_abl_mutations=get_report_bcr_abl_mutations,
        bcr_abl_phasing=get_report_bcr_abl_phasing,
        bcr_abl_exon_coverage=get_report_bcr_abl_exon_coverage,
    output:
        xlsx="reports/{sample}_{type}.nagini_report.xlsx",

    params:
        reference_build=lambda wildcards: config.get("reference", {}).get("build", config.get("reference_build", "hg38")),
        snv_annotation=(

            "VEP" if vep_enabled() else "none"
        ),
        cytoband=get_cytoband_file,
        small_variant_caller=get_small_variant_caller(),
        design_bed=lambda wildcards: config.get("reference", {}).get("design_bed", ""),
        clairs_to_platform=lambda wildcards: config.get("clairs_to_call", {}).get("platform", ""),
        barcode=get_sample_barcode_name,
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
