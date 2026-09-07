__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnv_rindels:
    """Run Rindels in CNV-only mode to produce whole-genome binned read counts and ASCAT segments."""
    input:
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
        ref=config["reference"]["fasta"],
    output:
        bins="cnv/rindels/{sample}_{type}.cnv_bins.bed",
        ascat="cnv/rindels/{sample}_{type}.ascat.bed",
        vcf=temp("cnv/rindels/{sample}_{type}.rindels.vcf"),
    params:
        binary=lambda wildcards: _find_rindels_binary(),
        bin_size=config.get("rindels", {}).get("cnv_bin_size", 50000),
        extra=config.get("rindels", {}).get("extra_cnv", "--ont --only-cnv"),
    log:
        "logs/cnv_rindels/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/cnv_rindels/{sample}_{type}.tsv",
            config.get("cnv_rindels", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("rindels", {}).get("threads", 4)
    resources:
        mem_mb=config.get("rindels", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("rindels", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("rindels", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("rindels", {}).get("threads", 4),
        time=config.get("rindels", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: call whole-genome copy number variations with Rindels for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        {params.binary} \
            -b {input.bam} \
            -r {input.ref} \
            --cnv-output {output.bins} \
            --cnv-bin-size {params.bin_size} \
            -o {output.vcf} \
            {params.extra} &> {log}
        """


rule cnv_summary:
    """Summarise Rindels ASCAT CNV segments and annotate with cytobands and target regions."""
    input:
        ascat="cnv/rindels/{sample}_{type}.ascat.bed",
        target_bed=lambda wildcards: get_resolved_design_bed() or [],
        cytoband=lambda wildcards: get_cytoband_file(wildcards) or [],
    output:
        tsv="reports/cnv_summary/{sample}_{type}.cnv.tsv",
        json="reports/cnv_summary/{sample}_{type}.cnv.json",
    params:
        extra=config.get("cnv_summary", {}).get("extra", ""),
    log:
        "logs/cnv_summary/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/cnv_summary/{sample}_{type}.tsv",
            config.get("cnv_summary", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("cnv_summary", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnv_summary", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnv_summary", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("cnv_summary", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnv_summary", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: summarise copy number alterations for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/cnv_summary.py"
