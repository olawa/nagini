__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule samtools_merge_bam:
    """Merge the per-unit BAMs into one coordinate-sorted BAM per sample/type."""
    input:
        bams=get_unit_alignment_bams,
    output:
        bam="alignment/quicken_bam/{sample}_{type}.bam",
        bai="alignment/quicken_bam/{sample}_{type}.bam.bai",
    params:
        extra=config.get("samtools_merge_bam", {}).get("extra", ""),
    log:
        "logs/samtools_merge_bam/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/samtools_merge_bam/{sample}_{type}.tsv",
            config.get("samtools_merge_bam", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("samtools_merge_bam", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("samtools_merge_bam", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("samtools_merge_bam", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("samtools_merge_bam", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("samtools_merge_bam", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("samtools_merge_bam", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("samtools_merge_bam", {}).get("container", config["default_container"])
    conda:
        "../envs/quicken.yaml"
    message:
        "{rule}: merge per-unit BAMs into {output.bam}"
    shell:
        """
        samtools merge -f -@ {threads} {params.extra} {output.bam} {input.bams} &> {log}
        samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1
        """
