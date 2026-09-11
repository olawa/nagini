__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"

import os


def _get_ref(key, config):
    """Return a reference file path resolved via reference_build."""
    build = config.get("reference_build", "hg38")
    return config["references"][build][key]


rule cnvkit_batch:
    """Run CNVkit in WGS mode on a tumor-only ONT BAM.

    No Panel of Normals is needed: -n builds an internal flat reference.
    --no-gc skips GC-bias correction (ONT is PCR-free, so correction
    introduces more noise than it removes).
    --count-reads counts read midpoints instead of computing per-base depth,
    which is orders of magnitude faster for long-read data.
    """
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        fasta=lambda wildcards: _get_ref("fasta", config),
    output:
        cnr="cnv/cnvkit/{sample}_{type}.cnr",
        cns="cnv/cnvkit/{sample}_{type}.cns",
    params:
        tmpdir="cnv/cnvkit/tmp_{sample}_{type}",
        method=config.get("cnvkit_batch", {}).get("method", "wgs"),
        target_avg_size=lambda wildcards: f"--target-avg-size {config.get('cnvkit_batch', {}).get('target_avg_size', 50000)}" if config.get("cnvkit_batch", {}).get("target_avg_size", 50000) else "",
        extra=config.get("cnvkit_batch", {}).get("extra", "--count-reads"),
    log:
        "cnv/cnvkit/{sample}_{type}.cnvkit_batch.log",
    benchmark:
        repeat(
            "cnv/cnvkit/{sample}_{type}.cnvkit_batch.benchmark.tsv",
            config.get("cnvkit_batch", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_batch", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_batch", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_batch", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_batch", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_batch", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_batch", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_batch", {}).get("container", config["default_container"])
    message:
        "{rule}: CNV calling with CNVkit (WGS, tumor-only) for {wildcards.sample}_{wildcards.type}"
    shell:
        "mkdir -p {params.tmpdir} && "
        "cnvkit.py batch {input.bam} -n -f {input.fasta} "
        "--output-reference {params.tmpdir}/reference.cnn "
        "-m {params.method} {params.target_avg_size} -d {params.tmpdir} -p {threads} {params.extra} "
        "> {log} 2>&1 && "
        "bam_base=$(basename {input.bam} .bam) && "
        "mv {params.tmpdir}/${{bam_base}}.cnr {output.cnr} && "
        "mv {params.tmpdir}/${{bam_base}}.cns {output.cns} && "
        "rm -rf {params.tmpdir}"
