__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule mosdepth_targets:
    """Calculate target-region and per-base coverage with the Hydra wrapper."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        bed=lambda wildcards: get_resolved_design_bed(),
    output:
        bed="qc/mosdepth_bed/{sample}_{type}.regions.bed.gz",
        per_base="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz",
        dist="qc/mosdepth_bed/{sample}_{type}.mosdepth.global.dist.txt",
        bed_dist="qc/mosdepth_bed/{sample}_{type}.mosdepth.region.dist.txt",
        summary="qc/mosdepth_bed/{sample}_{type}.mosdepth.summary.txt",
    params:
        extra=config.get("mosdepth_bed", {}).get("extra", ""),
    log:
        "logs/mosdepth_targets/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/mosdepth_targets/{sample}_{type}.tsv",
            config.get("mosdepth_bed", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mosdepth_bed", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mosdepth_bed", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosdepth_bed", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mosdepth_bed", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mosdepth_bed", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mosdepth_bed", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mosdepth_bed", {}).get("container", config["default_container"])
    message:
        "{rule}: calculate target coverage for {input.bam}"
    wrapper:
        "0.80.2/bio/mosdepth"
