__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule quicken_excel_report:
    """Build a portable multi-sheet workbook from final QC and variant outputs."""
    input:
        snv_vcf=get_small_variant_report_vcf,
        snv_tbi=get_small_variant_report_vcf_index,
        sv="reports/sv_summary/{sample}_{type}.tsv",
        translocations="reports/sv_summary/{sample}_{type}.translocations.tsv",
        wgs_coverage="qc/mosdepth/{sample}_{type}.mosdepth.summary.txt",
        target_coverage="qc/mosdepth_bed/{sample}_{type}.mosdepth.summary.txt",
        cramino="qc/cramino/{sample}_{type}.txt",
    output:
        xlsx="reports/{sample}_{type}.quicken_report.xlsx",
    params:
        reference_build=config["reference"]["build"],
        snv_annotation=(
            "VEP" if vep_enabled() else "none"
        ),
        cytoband=get_cytoband_file,
    log:
        "logs/quicken_excel_report/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/quicken_excel_report/{sample}_{type}.tsv",
            config.get("quicken_excel_report", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("quicken_excel_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("quicken_excel_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("quicken_excel_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("quicken_excel_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("quicken_excel_report", {}).get("container", config["default_container"])
    conda:
        "../envs/quicken.yaml"
    message:
        "{rule}: compile Excel report for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/quicken_report.py"
