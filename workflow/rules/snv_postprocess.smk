__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"

import os


def get_vep_custom_annotations(wildcards):
    """Add ClinVar as an optional, indexed VEP custom annotation."""

    clinvar = get_clinvar_vcf(wildcards)
    if not clinvar or not os.path.isfile(f"{clinvar}.tbi"):
        return ""
    return (
        f"--custom {clinvar},short_name=ClinVar,format=vcf,type=exact,"
        "fields=CLNSIG%CLNREVSTAT%CLNDN"
    )


rule create_genome_contigs_bed:
    """Create a whole-genome BED from FASTA index when no design BED is provided."""
    input:
        fai=config["reference"]["fai"],
    output:
        bed="references/genome_contigs.bed",
    log:
        "logs/create_genome_contigs_bed/genome_contigs.log",
    benchmark:
        repeat(
            "benchmarks/create_genome_contigs_bed/genome_contigs.tsv",
            config.get("create_genome_contigs_bed", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("create_genome_contigs_bed", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_genome_contigs_bed", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_genome_contigs_bed", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("create_genome_contigs_bed", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_genome_contigs_bed", {}).get("container", config["default_container"])
    message:
        "{rule}: write a whole-genome BED from {input.fai}"
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} {{print $1, 0, $2}}' {input.fai} > {output.bed} 2> {log}
        """
