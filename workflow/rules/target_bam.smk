__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule pad_target_bed:
    """Pad the design BED file by +/- 25 kb (or configurable bp) and merge overlapping intervals."""
    input:
        bed=lambda wildcards: get_resolved_design_bed(),
    output:
        bed="alignment/target/targets_padded_25kb.bed",
    params:
        padding_bp=lambda wildcards: config.get("target_bam", {}).get("padding_bp", 25000),
    log:
        "logs/target_bam/pad_target_bed.log",
    benchmark:
        repeat(
            "benchmarks/target_bam/pad_target_bed.tsv",
            config.get("pad_target_bed", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("pad_target_bed", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pad_target_bed", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pad_target_bed", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("pad_target_bed", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("pad_target_bed", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: pad design BED intervals by +/-25kb for target region BAM extraction"
    script:
        "../scripts/pad_bed.py"


rule extract_target_bam:
    """Extract reads overlapping padded target regions into a target-only BAM."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        bed="alignment/target/targets_padded_25kb.bed",
    output:
        bam="alignment/target/{sample}_{type}.target.bam",
        bai="alignment/target/{sample}_{type}.target.bam.bai",
    params:
        extra=config.get("extract_target_bam", {}).get("extra", ""),
    log:
        "logs/target_bam/{sample}_{type}.extract_target_bam.log",
    benchmark:
        repeat(
            "benchmarks/target_bam/{sample}_{type}.extract_target_bam.tsv",
            config.get("extract_target_bam", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("extract_target_bam", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("extract_target_bam", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("extract_target_bam", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("extract_target_bam", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("extract_target_bam", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("extract_target_bam", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("extract_target_bam", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: extract target region BAM (+/-25kb) for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        samtools view -b -L {input.bed} -@ {threads} {params.extra} {input.bam} -o {output.bam} &> {log}
        samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1
        """
