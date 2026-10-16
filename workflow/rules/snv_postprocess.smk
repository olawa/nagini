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


rule vep_annotate_mutect2:
    """Annotate normalized small variants with an offline VEP cache."""
    input:
        vcf="snv_indels/mutect2_postprocess/{sample}_{type}.filtered.norm.vcf.gz",
        tbi="snv_indels/mutect2_postprocess/{sample}_{type}.filtered.norm.vcf.gz.tbi",
        fasta=config["reference"]["fasta"],
        cache=lambda wildcards: config["vep"]["vep_cache"],
    output:
        vcf=temp("snv_indels/mutect2_postprocess/{sample}_{type}.filtered.norm.vep.vcf"),
        stats="snv_indels/mutect2_postprocess/{sample}_{type}.filtered.norm.vep.summary.html",
    params:
        mode=config.get("vep", {}).get("mode", "--offline --cache"),
        custom=get_vep_custom_annotations,
        extra=config.get("vep", {}).get("extra", "--assembly GRCh38 --check_existing --pick --symbol --hgvs"),
    log:
        "logs/vep_annotate_mutect2/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/vep_annotate_mutect2/{sample}_{type}.tsv",
            config.get("vep", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("vep", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("vep", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vep", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vep", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vep", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vep", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vep", {}).get("container", config["default_container"])
    conda:
        "../../envs/vep.yaml"
    message:
        "{rule}: VEP annotate filtered Mutect2 calls for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        vep \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --vcf \
            --force_overwrite \
            --no_progress \
            --stats_file {output.stats} \
            --dir_cache {input.cache} \
            --fasta {input.fasta} \
            --fork {threads} \
            {params.mode} \
            {params.custom} \
            {params.extra} &> {log}
        """


rule vep_annotate_clairs_to:
    """Annotate normalized ClairS-TO small variants with an offline VEP cache."""
    input:
        vcf="snv_indels/clairs_to_postprocess/{sample}_{type}.filtered.norm.vcf.gz",
        tbi="snv_indels/clairs_to_postprocess/{sample}_{type}.filtered.norm.vcf.gz.tbi",
        fasta=config["reference"]["fasta"],
        cache=lambda wildcards: config["vep"]["vep_cache"],
    output:
        vcf=temp("snv_indels/clairs_to_postprocess/{sample}_{type}.filtered.norm.vep.vcf"),
        stats="snv_indels/clairs_to_postprocess/{sample}_{type}.filtered.norm.vep.summary.html",
    params:
        mode=config.get("vep", {}).get("mode", "--offline --cache"),
        custom=get_vep_custom_annotations,
        extra=config.get("vep", {}).get("extra", "--assembly GRCh38 --check_existing --pick --symbol --hgvs"),
    log:
        "logs/vep_annotate_clairs_to/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/vep_annotate_clairs_to/{sample}_{type}.tsv",
            config.get("vep", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("vep", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("vep", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vep", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vep", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vep", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vep", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vep", {}).get("container", config["default_container"])
    conda:
        "../../envs/vep.yaml"
    message:
        "{rule}: VEP annotate ClairS-TO calls for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        vep \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --vcf \
            --force_overwrite \
            --no_progress \
            --stats_file {output.stats} \
            --dir_cache {input.cache} \
            --fasta {input.fasta} \
            --fork {threads} \
            {params.mode} \
            {params.custom} \
            {params.extra} &> {log}
        """


rule create_genome_contigs_bed:
    """Create a whole-genome BED from FASTA index when no design BED is provided."""
    input:
        fai=config["reference"]["fai"],
    output:
        bed="references/genome_contigs.bed",
    threads: 1
    container:
        config["default_container"]
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} {{print $1, 0, $2}}' {input.fai} > {output.bed}
        """


rule vep_annotate_deepvariant:
    """Annotate normalized DeepVariant small variants with an offline VEP cache."""
    input:
        vcf="snv_indels/deepvariant_postprocess/{sample}_{type}.filtered.norm.vcf.gz",
        tbi="snv_indels/deepvariant_postprocess/{sample}_{type}.filtered.norm.vcf.gz.tbi",
        fasta=config["reference"]["fasta"],
        cache=lambda wildcards: config["vep"]["vep_cache"],
    output:
        vcf=temp("snv_indels/deepvariant_postprocess/{sample}_{type}.filtered.norm.vep.vcf"),
        stats="snv_indels/deepvariant_postprocess/{sample}_{type}.filtered.norm.vep.summary.html",
    params:
        mode=config.get("vep", {}).get("mode", "--offline --cache"),
        custom=get_vep_custom_annotations,
        extra=config.get("vep", {}).get("extra", "--assembly GRCh38 --check_existing --pick --symbol --hgvs"),
    log:
        "logs/vep_annotate_deepvariant/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/vep_annotate_deepvariant/{sample}_{type}.tsv",
            config.get("vep", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("vep", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("vep", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vep", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vep", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vep", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vep", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vep", {}).get("container", config["default_container"])
    conda:
        "../../envs/vep.yaml"
    message:
        "{rule}: VEP annotate DeepVariant calls for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        vep \
            --input_file {input.vcf} \
            --output_file {output.vcf} \
            --vcf \
            --force_overwrite \
            --no_progress \
            --stats_file {output.stats} \
            --dir_cache {input.cache} \
            --fasta {input.fasta} \
            --fork {threads} \
            {params.mode} \
            {params.custom} \
            {params.extra} &> {log}
        """

