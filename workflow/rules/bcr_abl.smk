__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule bcr_abl_fusion_counts:
    """Identify and quantify BCR-ABL1 fusion isoforms from spliced alignments."""
    input:
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
        features=lambda wildcards: (
            config.get("bcr_abl", {}).get("features_json", "")
            or os.path.join(workflow.basedir, "resources/bcr_abl/bcr_abl_features.json")
        ),
    output:
        tsv="reports/bcr_abl/{sample}_{type}.fusion_summary.tsv",
        json="reports/bcr_abl/{sample}_{type}.fusion_summary.json",
    params:
        tolerance=lambda wildcards: config.get("bcr_abl", {}).get("tolerance", 150),
        min_mapq=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 1),
        genome_build=lambda wildcards: config.get("reference", {}).get("build", "hg38"),
    log:
        "logs/bcr_abl_fusion_counts/{sample}_{type}.log",
    threads: config.get("bcr_abl", {}).get("threads", 1)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 1),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        python3 {workflow.basedir}/scripts/bcr_abl_fusion_counter.py \
            --bam {input.bam} \
            --sample {wildcards.sample} \
            --genome-build {params.genome_build} \
            --features-json {input.features} \
            --output {output.tsv} \
            --json-output {output.json} \
            --tolerance {params.tolerance} \
            --min-mapq {params.min_mapq} > {log} 2>&1
        """


rule bcr_abl_branford_calling:
    """Force-call BCR-ABL1 kinase domain resistance mutations at high depth."""
    input:
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
        bed=lambda wildcards: (
            config.get("bcr_abl", {}).get("branford_bed", "")
            or os.path.join(
                workflow.basedir,
                f"resources/bcr_abl/branford_positions_{'hg19' if '19' in str(config.get('reference', {}).get('build', '')).lower() else 'hg38'}.bed",
            )
        ),
        features=lambda wildcards: (
            config.get("bcr_abl", {}).get("features_json", "")
            or os.path.join(workflow.basedir, "resources/bcr_abl/bcr_abl_features.json")
        ),
    output:
        all_sites="reports/bcr_abl/{sample}_{type}.branford_all_sites.tsv",
        mutations="reports/bcr_abl/{sample}_{type}.branford_mutations.tsv",
        vcf="reports/bcr_abl/{sample}_{type}.branford_mutations.vcf",
    params:
        min_af=lambda wildcards: get_assay_setting(wildcards.sample, "bcr_abl", "min_af", config.get("bcr_abl", {}).get("min_af", 0.01)),
        min_depth=lambda wildcards: config.get("bcr_abl", {}).get("min_depth", 50),
        min_alt_reads=lambda wildcards: config.get("bcr_abl", {}).get("min_alt_reads", 5),
    log:
        "logs/bcr_abl_branford_calling/{sample}_{type}.log",
    threads: config.get("bcr_abl", {}).get("threads", 1)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 1),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        python3 {workflow.basedir}/scripts/bcr_abl_branford_caller.py \
            --bam {input.bam} \
            --bed {input.bed} \
            --sample {wildcards.sample} \
            --output-all {output.all_sites} \
            --output-mutations {output.mutations} \
            --output-vcf {output.vcf} \
            --features-json {input.features} \
            --min-af {params.min_af} \
            --min-depth {params.min_depth} \
            --min-alt-reads {params.min_alt_reads} > {log} 2>&1
        """


rule bcr_abl_phasing:
    """Single-molecule phasing of compound resistance mutations from long reads."""
    input:
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
        mutations="reports/bcr_abl/{sample}_{type}.branford_mutations.tsv",
    output:
        tsv="reports/bcr_abl/{sample}_{type}.phasing.tsv",
        json="reports/bcr_abl/{sample}_{type}.phasing.json",
    params:
        min_spanning_reads=lambda wildcards: config.get("bcr_abl", {}).get("min_spanning_reads", 10),
        min_mapq=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 1),
    log:
        "logs/bcr_abl_phasing/{sample}_{type}.log",
    threads: config.get("bcr_abl", {}).get("threads", 1)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 1),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        python3 {workflow.basedir}/scripts/bcr_abl_phaser.py \
            --bam {input.bam} \
            --mutations-tsv {input.mutations} \
            --sample {wildcards.sample} \
            --output-tsv {output.tsv} \
            --output-json {output.json} \
            --min-spanning-reads {params.min_spanning_reads} \
            --min-mapq {params.min_mapq} > {log} 2>&1
        """
