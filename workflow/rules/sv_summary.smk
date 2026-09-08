__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule sv_summary:
    """Summarise Sniffles2 SVs and extract translocations/BNDs."""
    input:
        vcf=lambda wildcards: (
            "cnv_sv/sniffles2_call/{sample}_{type}.snpeff_annotated.ss_annotated.vcf"
            if config.get("annotate_svs", False)
            else "cnv_sv/sniffles2_call/{sample}_{type}.vcf.gz"
        ),
        tbi=lambda wildcards: (
            []
            if config.get("annotate_svs", False)
            else "cnv_sv/sniffles2_call/{sample}_{type}.vcf.gz.tbi"
        ),
        annotation_bed=lambda wildcards: (
            config.get("reference", {}).get("sv_annotation_bed", "")
            or get_sample_design_bed(wildcards.sample)
            or []
        ),
        cytoband=lambda wildcards: get_cytoband_file(wildcards) or [],
    output:
        tsv="reports/sv_summary/{sample}_{type}.tsv",
        translocations="reports/sv_summary/{sample}_{type}.translocations.tsv",
        json="reports/sv_summary/{sample}_{type}.json",
    params:
        extra=config.get("sv_summary", {}).get("extra", ""),
    log:
        "logs/sv_summary/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/sv_summary/{sample}_{type}.tsv",
            config.get("sv_summary", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("sv_summary", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("sv_summary", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("sv_summary", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("sv_summary", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("sv_summary", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("sv_summary", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("sv_summary", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: summarise structural variants in {input.vcf}"
    script:
        "../scripts/sv_summary.py"
