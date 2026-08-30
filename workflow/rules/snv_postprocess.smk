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


rule gatk_merge_mutect_stats:
    """Merge the per-chromosome Mutect2 statistics required for filtering."""
    input:
        stats=expand(
            "snv_indels/gatk_mutect2/{{sample}}_{{type}}_{chr}.unfiltered.vcf.gz.stats",
            chr=reference_chromosomes,
        ),
    output:
        stats=temp("snv_indels/mutect2_postprocess/{sample}_{type}.unfiltered.merged.stats"),
    params:
        stats=lambda wildcards, input: " ".join(
            f"--stats {shlex.quote(str(path))}" for path in input.stats
        ),
        extra=config.get("gatk_merge_mutect_stats", {}).get("extra", ""),
    log:
        "logs/gatk_merge_mutect_stats/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/gatk_merge_mutect_stats/{sample}_{type}.tsv",
            config.get("gatk_merge_mutect_stats", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_merge_mutect_stats", {}).get("threads", 1)
    resources:
        mem_mb=config.get("gatk_merge_mutect_stats", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_merge_mutect_stats", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_merge_mutect_stats", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_merge_mutect_stats", {}).get("threads", 1),
        time=config.get("gatk_merge_mutect_stats", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_merge_mutect_stats", {}).get("container", config.get("gatk_mutect2", {}).get("container", config["default_container"]))
    conda:
        "../../envs/gatk.yaml"
    message:
        "{rule}: merge Mutect2 statistics for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        gatk MergeMutectStats {params.stats} {params.extra} \
            --output {output.stats} &> {log}
        """


rule gatk_filter_mutect_calls:
    """Apply GATK's Mutect2 filters to the concatenated callset."""
    input:
        vcf="snv_indels/bcftools_concat/{sample}_{type}.unfiltered.vcf.gz",
        tbi="snv_indels/bcftools_concat/{sample}_{type}.unfiltered.vcf.gz.tbi",
        stats="snv_indels/mutect2_postprocess/{sample}_{type}.unfiltered.merged.stats",
        ref=config["reference"]["fasta"],
        fai=config["reference"]["fai"],
        dictionary=config["reference"]["dict"],
    output:
        vcf="snv_indels/mutect2_postprocess/{sample}_{type}.filtered.vcf.gz",
        stats="snv_indels/mutect2_postprocess/{sample}_{type}.filteringStats.tsv",
    params:
        extra=config.get("gatk_filter_mutect_calls", {}).get("extra", ""),
    log:
        "logs/gatk_filter_mutect_calls/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/gatk_filter_mutect_calls/{sample}_{type}.tsv",
            config.get("gatk_filter_mutect_calls", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("gatk_filter_mutect_calls", {}).get("threads", 1)
    resources:
        mem_mb=config.get("gatk_filter_mutect_calls", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("gatk_filter_mutect_calls", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("gatk_filter_mutect_calls", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("gatk_filter_mutect_calls", {}).get("threads", 1),
        time=config.get("gatk_filter_mutect_calls", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("gatk_filter_mutect_calls", {}).get("container", config.get("gatk_mutect2", {}).get("container", config["default_container"]))
    conda:
        "../../envs/gatk.yaml"
    message:
        "{rule}: filter Mutect2 calls for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        gatk FilterMutectCalls \
            --reference {input.ref} \
            --variant {input.vcf} \
            --stats {input.stats} \
            --filtering-stats {output.stats} \
            --create-output-variant-index false \
            {params.extra} \
            --output {output.vcf} &> {log}
        """


rule bcftools_norm_mutect2:
    """Left-align, normalize and split multiallelic filtered Mutect2 calls."""
    input:
        vcf="snv_indels/mutect2_postprocess/{sample}_{type}.filtered.vcf.gz",
        tbi="snv_indels/mutect2_postprocess/{sample}_{type}.filtered.vcf.gz.tbi",
        ref=config["reference"]["fasta"],
        fai=config["reference"]["fai"],
    output:
        vcf="snv_indels/mutect2_postprocess/{sample}_{type}.filtered.norm.vcf.gz",
        tbi="snv_indels/mutect2_postprocess/{sample}_{type}.filtered.norm.vcf.gz.tbi",
    params:
        extra=config.get("bcftools_norm_mutect2", {}).get("extra", "--multiallelics -both --check-ref w"),
    log:
        "logs/bcftools_norm_mutect2/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcftools_norm_mutect2/{sample}_{type}.tsv",
            config.get("bcftools_norm_mutect2", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_norm_mutect2", {}).get("threads", 1)
    resources:
        mem_mb=config.get("bcftools_norm_mutect2", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_norm_mutect2", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_norm_mutect2", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_norm_mutect2", {}).get("threads", 1),
        time=config.get("bcftools_norm_mutect2", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_norm_mutect2", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: normalize filtered Mutect2 calls for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        bcftools norm \
            --fasta-ref {input.ref} \
            {params.extra} \
            --output-type z \
            --output {output.vcf} \
            {input.vcf} &> {log}
        bcftools index --force --tbi {output.vcf} >> {log} 2>&1
        """


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


rule bcftools_norm_clairs_to:
    """Left-align, normalize and split multiallelic ClairS-TO calls."""
    input:
        vcf="snv_indels/clairs_to/{sample}_{type}.snv-indels.vcf.gz",
        tbi="snv_indels/clairs_to/{sample}_{type}.snv-indels.vcf.gz.tbi",
        ref=config["reference"]["fasta"],
        fai=config["reference"]["fai"],
    output:
        vcf="snv_indels/clairs_to_postprocess/{sample}_{type}.filtered.norm.vcf.gz",
        tbi="snv_indels/clairs_to_postprocess/{sample}_{type}.filtered.norm.vcf.gz.tbi",
    params:
        extra=config.get("bcftools_norm_clairs_to", {}).get("extra", "--multiallelics -both --check-ref w"),
    log:
        "logs/bcftools_norm_clairs_to/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcftools_norm_clairs_to/{sample}_{type}.tsv",
            config.get("bcftools_norm_clairs_to", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_norm_clairs_to", {}).get("threads", 1)
    resources:
        mem_mb=config.get("bcftools_norm_clairs_to", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_norm_clairs_to", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_norm_clairs_to", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_norm_clairs_to", {}).get("threads", 1),
        time=config.get("bcftools_norm_clairs_to", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_norm_clairs_to", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: normalize ClairS-TO calls for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        bcftools norm \
            --fasta-ref {input.ref} \
            {params.extra} \
            --output-type z \
            --output {output.vcf} \
            {input.vcf} &> {log}
        bcftools index --force --tbi {output.vcf} >> {log} 2>&1
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
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} {{print $1, 0, $2}}' {input.fai} > {output.bed}
        """


rule bcftools_norm_deepvariant:
    """Left-align, normalize and split multiallelic DeepVariant calls."""
    input:
        vcf="snv_indels/deepvariant/{sample}_{type}.vcf.gz",
        tbi="snv_indels/deepvariant/{sample}_{type}.vcf.gz.tbi",
        ref=config["reference"]["fasta"],
        fai=config["reference"]["fai"],
    output:
        vcf="snv_indels/deepvariant_postprocess/{sample}_{type}.filtered.norm.vcf.gz",
        tbi="snv_indels/deepvariant_postprocess/{sample}_{type}.filtered.norm.vcf.gz.tbi",
    params:
        extra=config.get("bcftools_norm_deepvariant", {}).get("extra", "--multiallelics -both --check-ref w"),
    log:
        "logs/bcftools_norm_deepvariant/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcftools_norm_deepvariant/{sample}_{type}.tsv",
            config.get("bcftools_norm_deepvariant", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_norm_deepvariant", {}).get("threads", 1)
    resources:
        mem_mb=config.get("bcftools_norm_deepvariant", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_norm_deepvariant", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_norm_deepvariant", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_norm_deepvariant", {}).get("threads", 1),
        time=config.get("bcftools_norm_deepvariant", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_norm_deepvariant", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: normalize DeepVariant calls for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        bcftools norm \
            --fasta-ref {input.ref} \
            {params.extra} \
            --output-type z \
            --output {output.vcf} \
            {input.vcf} &> {log}
        bcftools index --force --tbi {output.vcf} >> {log} 2>&1
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

