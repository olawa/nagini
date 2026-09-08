__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule mosdepth_genome:
    """Calculate windowed whole-genome coverage with the Hydra wrapper."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
    output:
        bed="qc/mosdepth/{sample}_{type}.regions.bed.gz",
        dist="qc/mosdepth/{sample}_{type}.mosdepth.global.dist.txt",
        bed_dist="qc/mosdepth/{sample}_{type}.mosdepth.region.dist.txt",
        summary="qc/mosdepth/{sample}_{type}.mosdepth.summary.txt",
    params:
        by=config.get("mosdepth", {}).get("by", 1000),
        extra=config.get("mosdepth", {}).get("extra", ""),
    log:
        "logs/mosdepth_genome/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/mosdepth_genome/{sample}_{type}.tsv",
            config.get("mosdepth", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mosdepth", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosdepth", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mosdepth", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mosdepth", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mosdepth", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mosdepth", {}).get("container", config["default_container"])
    message:
        "{rule}: calculate whole-genome coverage for {input.bam}"
    wrapper:
        "0.80.2/bio/mosdepth"
