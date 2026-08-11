__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule stage_aligned_bam:
    """Prepare an existing aligned BAM without invoking an aligner."""
    input:
        bam=lambda wildcards: get_unit_input(wildcards, "bam"),
    output:
        bam="alignment/quicken_align/{sample}_{type}_{processing_unit}_{barcode}.aligned_bam.bam",
        bai="alignment/quicken_align/{sample}_{type}_{processing_unit}_{barcode}.aligned_bam.bam.bai",
    params:
        extra=config.get("stage_aligned_bam", {}).get("extra", ""),
    log:
        "logs/stage_aligned_bam/{sample}_{type}_{processing_unit}_{barcode}.log",
    benchmark:
        repeat(
            "benchmarks/stage_aligned_bam/{sample}_{type}_{processing_unit}_{barcode}.tsv",
            config.get("stage_aligned_bam", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("stage_aligned_bam", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("stage_aligned_bam", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("stage_aligned_bam", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("stage_aligned_bam", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("stage_aligned_bam", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("stage_aligned_bam", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("stage_aligned_bam", {}).get("container", config["default_container"])
    conda:
        "../envs/quicken.yaml"
    message:
        "{rule}: sort/index an existing aligned BAM; no realignment is performed"
    shell:
        """
        samtools sort -@ {threads} {params.extra} -o {output.bam} {input.bam} &> {log}
        samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1
        """
