__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule merged_fastq:
    """One FASTQ per sample from its unaligned source BAMs, mod tags kept.

    Each BAM is converted separately and the BGZF streams are concatenated,
    which is a valid bgzip file.  Merging the BAMs first would mean reconciling
    headers from different flowcells for no gain, since the read groups are
    already in the per-unit alignments.
    """
    input:
        bams=get_sample_source_bams,
    output:
        fastq="reads/merged/{sample}_{type}.fastq.gz",
    params:
        extra=config.get("samtools_fastq", {}).get("extra", "-n -T MM,ML"),
    log:
        "logs/merged_fastq/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/merged_fastq/{sample}_{type}.tsv",
            config.get("samtools_fastq", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("samtools_fastq", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("samtools_fastq", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("samtools_fastq", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("samtools_fastq", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("samtools_fastq", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("samtools_fastq", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("samtools_fastq", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: extract one FASTQ for {wildcards.sample}_{wildcards.type} from its source BAMs"
    shell:
        """
        : > {output.fastq}
        for bam in {input.bams}; do
            samtools fastq -@ {threads} {params.extra} "$bam" 2>> {log} | \
            bgzip -@ {threads} -c >> {output.fastq} 2>> {log}
        done
        """
