__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule stage_aligned_bam:
    """Prepare an existing aligned BAM without invoking an aligner."""
    input:
        bam=lambda wildcards: get_unit_input(wildcards, "bam"),
    output:
        bam="alignment/nagini_align/{sample}_{type}_{processing_unit}_{barcode}.from_aligned.bam",
        bai="alignment/nagini_align/{sample}_{type}_{processing_unit}_{barcode}.from_aligned.bam.bai",
    params:
        extra=config.get("stage_aligned_bam", {}).get("extra", ""),
        # GATK identifies the sample from the BAM read-group SM tag.  Existing
        # aligned BAMs from Dorado commonly have no @RG header, so make the
        # sample identity explicit while staging the input.
        read_group=lambda wildcards: get_minimap2_read_group(wildcards),
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
        "../../envs/nagini.yaml"
    message:
        "{rule}: add sample metadata and sort/index an existing aligned BAM; no realignment is performed"
    shell:
        """
        samtools addreplacerg \
            -w \
            -r '{params.read_group}' \
            -m overwrite_all \
            -o {output.bam}.with_rg.bam \
            {input.bam} &> {log}
        samtools sort -@ {threads} {params.extra} -o {output.bam} {output.bam}.with_rg.bam >> {log} 2>&1
        rm -f {output.bam}.with_rg.bam
        samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1
        """
