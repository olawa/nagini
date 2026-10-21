__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


def _bcr_abl_resource(filename):
    """Return absolute path to a bundled BCR-ABL1 resource file."""
    return os.path.join(workflow.basedir, "resources", "bcr_abl", filename)


def _bcr_abl_build():
    """Only select bundled coordinates for a reference they describe."""
    build = str(config.get("reference_build", "hg38")).strip().lower()
    aliases = {"hg19": "hg19", "grch37": "hg19", "b37": "hg19",
               "hg38": "hg38", "grch38": "hg38"}
    if build not in aliases:
        raise WorkflowError(
            "No bundled BCR-ABL1 resources for reference_build={!r}; "
            "supported builds are hg19/GRCh37 and hg38/GRCh38.".format(build)
        )
    return aliases[build]


def _branford_bed():
    """Pick the branford positions BED that matches the configured reference build."""
    override = config.get("bcr_abl", {}).get("branford_bed", "")
    if override:
        return override
    bed_name = f"branford_positions_{_bcr_abl_build()}.bed"
    return _bcr_abl_resource(bed_name)


def _bcr_abl_exons_bed():
    """Pick the BCR-ABL1 exon coordinates BED that matches the configured reference build."""
    override = config.get("bcr_abl", {}).get("exons_bed", "")
    if override:
        return override
    bed_name = f"bcr_abl1_exons_{_bcr_abl_build()}.bed"
    return _bcr_abl_resource(bed_name)


def _bcr_abl_gtf():
    """Pick the gene annotation GTF for BCR-ABL1 fusion plotting."""
    override = config.get("bcr_abl", {}).get("gtf", "")
    if override and os.path.isfile(override):
        return override
    ref_gtf = config.get("reference", {}).get("gene_annotation", "") or config.get("reference", {}).get("gtf", "")
    if ref_gtf and os.path.isfile(ref_gtf):
        return ref_gtf
    gtf_name = f"bcr_abl_{_bcr_abl_build()}.gtf"
    return _bcr_abl_resource(gtf_name)




rule bcr_abl_fusion_counts:
    """Identify and quantify BCR-ABL1 fusion isoforms from spliced alignments."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        features=lambda wildcards: (
            config.get("bcr_abl", {}).get("features_json", "")
            or _bcr_abl_resource("bcr_abl_features.json")
        ),
    output:
        tsv="reports/bcr_abl/{sample}_{type}.fusion_summary.tsv",
        json="reports/bcr_abl/{sample}_{type}.fusion_summary.json",
        arriba="reports/bcr_abl/{sample}_{type}.arriba_fusions.tsv",
    params:
        tolerance=lambda wildcards: config.get("bcr_abl", {}).get("tolerance", 150),
        min_mapq=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 1),
        genome_build=lambda wildcards: _bcr_abl_build(),
    log:
        "logs/bcr_abl_fusion_counts/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcr_abl_fusion_counts/{sample}_{type}.tsv",
            config.get("bcr_abl", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcr_abl", {}).get("threads", 4)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 4),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcr_abl", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: count BCR-ABL1 fusion-spanning reads for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/bcr_abl_fusion_counter.py"


rule bcr_abl_branford_calling:
    """Force-call BCR-ABL1 kinase domain resistance mutations at high depth."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        bed=lambda wildcards: _branford_bed(),
        features=lambda wildcards: (
            config.get("bcr_abl", {}).get("features_json", "")
            or _bcr_abl_resource("bcr_abl_features.json")
        ),
        ref_fasta=lambda wildcards: config.get("reference", {}).get("fasta", ""),
        known_variants=lambda wildcards: (
            [config["bcr_abl"]["known_variants"]]
            if config.get("bcr_abl", {}).get("known_variants") else []
        ),
    output:
        all_sites="reports/bcr_abl/{sample}_{type}.branford_all_sites.tsv",
        mutations="reports/bcr_abl/{sample}_{type}.branford_mutations.tsv",
        vcf="reports/bcr_abl/{sample}_{type}.branford_mutations.vcf",
    params:
        min_af=lambda wildcards: get_assay_setting(wildcards.sample, "bcr_abl", "min_af", config.get("bcr_abl", {}).get("min_af", 0.01)),
        subclonal_af=lambda wildcards: get_assay_setting(wildcards.sample, "bcr_abl", "subclonal_af", config.get("bcr_abl", {}).get("subclonal_af", 0.002)),
        min_depth=lambda wildcards: config.get("bcr_abl", {}).get("min_depth", 50),
        min_alt_reads=lambda wildcards: config.get("bcr_abl", {}).get("min_alt_reads", 5),
        max_depth=lambda wildcards: config.get("bcr_abl", {}).get("max_depth", 1000000),
        min_base_quality=lambda wildcards: config.get("bcr_abl", {}).get("min_base_quality", 12),
        min_mapping_quality=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 1),
        min_read_qs=lambda wildcards: config.get("bcr_abl", {}).get("min_read_qs", 0),
        min_snr=lambda wildcards: config.get("bcr_abl", {}).get("min_snr", 3.0),
        min_minor_strand_ratio=lambda wildcards: config.get("bcr_abl", {}).get("min_minor_strand_ratio", 0.05),
        min_mean_baseq=lambda wildcards: config.get("bcr_abl", {}).get("min_mean_baseq", 14.0),
        min_kmer_rate=lambda wildcards: config.get("bcr_abl", {}).get("min_kmer_rate", 0.60),
        max_sor=lambda wildcards: config.get("bcr_abl", {}).get("max_sor", 3.0),
        kmer_window=lambda wildcards: config.get("bcr_abl", {}).get("kmer_window", 3),
        extra=lambda wildcards: config.get("bcr_abl", {}).get("branford_caller_extra", ""),
    log:
        "logs/bcr_abl_branford_calling/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcr_abl_branford_calling/{sample}_{type}.tsv",
            config.get("bcr_abl", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcr_abl", {}).get("threads", 4)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 4),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcr_abl", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: force-call BCR-ABL1 kinase domain mutations for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/bcr_abl_branford_caller.py"



rule bcr_abl_phasing:
    """Single-molecule phasing of compound resistance mutations from long reads."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        mutations="reports/bcr_abl/{sample}_{type}.branford_mutations.tsv",
    output:
        tsv="reports/bcr_abl/{sample}_{type}.phasing.tsv",
        json="reports/bcr_abl/{sample}_{type}.phasing.json",
    params:
        min_spanning_reads=lambda wildcards: config.get("bcr_abl", {}).get("min_spanning_reads", 10),
        min_mapq=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 1),
        min_base_quality=lambda wildcards: config.get("bcr_abl", {}).get("min_base_quality", 12),
        min_read_qs=lambda wildcards: config.get("bcr_abl", {}).get("min_read_qs", 0),
    log:
        "logs/bcr_abl_phasing/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcr_abl_phasing/{sample}_{type}.tsv",
            config.get("bcr_abl", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcr_abl", {}).get("threads", 4)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 4),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcr_abl", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: phase BCR-ABL1 resistance mutations for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/bcr_abl_phaser.py"


rule bcr_abl_exon_coverage:
    """Calculate depth and coverage metrics across BCR and ABL1 exons."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        bed=lambda wildcards: _bcr_abl_exons_bed(),
    output:
        tsv="reports/bcr_abl/{sample}_{type}.exon_coverage.tsv",
    params:
        min_mapq=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 1),
        max_depth=lambda wildcards: config.get("bcr_abl", {}).get("max_depth", 1000000),
    log:
        "logs/bcr_abl_exon_coverage/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcr_abl_exon_coverage/{sample}_{type}.tsv",
            config.get("bcr_abl", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcr_abl", {}).get("threads", 4)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 4),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcr_abl", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: measure coverage across BCR and ABL1 exons for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/bcr_abl_exon_coverage.py"



rule bcr_abl_draw_fusions:
    """Generate Arriba fusion diagrams using draw_fusions.R."""
    input:
        fusions="reports/bcr_abl/{sample}_{type}.arriba_fusions.tsv",
        fusion_summary="reports/bcr_abl/{sample}_{type}.fusion_summary.tsv",
        mutations="reports/bcr_abl/{sample}_{type}.branford_mutations.tsv",
        phasing="reports/bcr_abl/{sample}_{type}.phasing.tsv",
        exon_coverage="reports/bcr_abl/{sample}_{type}.exon_coverage.tsv",
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
    output:
        pdf="reports/bcr_abl/{sample}_{type}.fusion_plot.pdf",
        png="reports/bcr_abl/{sample}_{type}.fusion_plot.png",
    params:
        gtf=lambda wildcards: _bcr_abl_gtf(),
        cytobands=lambda wildcards: config.get("reference", {}).get("cytoband", "") or config.get("reference", {}).get("cytobands", ""),
        protein_domains=lambda wildcards: config.get("arriba", {}).get("protein_domains", ""),
        build=lambda wildcards: str(config.get("reference_build", "hg19")).strip().lower(),
        extra=config.get("arriba", {}).get("extra", ""),
    log:
        "logs/bcr_abl_draw_fusions/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcr_abl_draw_fusions/{sample}_{type}.tsv",
            config.get("arriba", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("arriba", {}).get("threads", 1)
    resources:
        mem_mb=config.get("arriba", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("arriba", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("arriba", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("arriba", {}).get("threads", 1),
        time=config.get("arriba", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("arriba", {}).get("container", config.get("ARRIBA_CONTAINER", "docker://uhrigs/arriba:2.4.0"))
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: draw the BCR-ABL1 fusion plots for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        bash {workflow.basedir}/scripts/bcr_abl_draw_fusions.sh \
            --sample {wildcards.sample} \
            --fusions {input.fusions} \
            --fusion-summary {input.fusion_summary} \
            --mutations {input.mutations} \
            --phasing {input.phasing} \
            --exon-coverage {input.exon_coverage} \
            --bam {input.bam} \
            --gtf={params.gtf:q} \
            --cytobands={params.cytobands:q} \
            --protein-domains={params.protein_domains:q} \
            --build {params.build} \
            --extra={params.extra:q} \
            --output-pdf {output.pdf} \
            --output-png {output.png} > {log} 2>&1
        """


rule bcr_abl_bcftools_pileup:
    """Run bcftools mpileup and call across BCR-ABL1 target exons at full depth."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        ref=config["reference"]["fasta"],
        bed=lambda wildcards: _bcr_abl_exons_bed(),
    output:
        vcf="snv_indels/bcftools_mpileup/{sample}_{type}.bcr_abl.vcf.gz",
        tbi="snv_indels/bcftools_mpileup/{sample}_{type}.bcr_abl.vcf.gz.tbi",
    params:
        max_depth=lambda wildcards: config.get("bcr_abl", {}).get("max_depth", 100000),
        min_mapq=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 20),
        min_baseq=lambda wildcards: config.get("bcr_abl", {}).get("min_base_quality", 12),
        prior=lambda wildcards: config.get("bcr_abl", {}).get("pileup_prior", 0.01),
        extra=lambda wildcards: config.get("bcr_abl", {}).get("bcftools_mpileup_extra", ""),
    log:
        "logs/bcr_abl_bcftools_pileup/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcr_abl_bcftools_pileup/{sample}_{type}.tsv",
            config.get("bcr_abl", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcr_abl", {}).get("threads", 4)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 4),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools", {}).get("container", "docker://hydragenetics/common:3.1.1")
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: call BCR-ABL1 target exons with bcftools for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        exec > {log} 2>&1
        bcftools mpileup --threads {threads} -Ou -B -f {input.ref} -R {input.bed} \
            -d {params.max_depth} -q {params.min_mapq} -Q {params.min_baseq} \
            -a FORMAT/AD,FORMAT/DP {params.extra} {input.bam} | \
        bcftools call --threads {threads} -mv -P {params.prior} -Oz -o {output.vcf}

        bcftools index -t --threads {threads} {output.vcf}
        """


rule bcr_abl_rindels:
    """Run Rindels in targeted ONT BCR-ABL1 amplicon mode."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        ref=config["reference"]["fasta"],
    output:
        vcf="snv_indels/rindels/{sample}_{type}.rindels.vcf.gz",
        haplotypes="snv_indels/rindels/{sample}_{type}.haplotypes.tsv",
    params:
        binary=lambda wildcards: _find_rindels_binary(),
        genome=lambda wildcards: "hg19" if "19" in str(config.get("reference_build", "hg19")).lower() else "hg38",
        extra=lambda wildcards: config.get("rindels", {}).get("extra", "--ont --amplicon --bcr-abl"),
    log:
        "logs/bcr_abl_rindels/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcr_abl_rindels/{sample}_{type}.tsv",
            config.get("bcr_abl_rindels", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcr_abl_rindels", {}).get("threads", config.get("rindels", {}).get("threads", 4))
    resources:
        mem_mb=config.get("bcr_abl_rindels", {}).get("mem_mb", config.get("rindels", {}).get("mem_mb", config["default_resources"]["mem_mb"])),
        mem_per_cpu=config.get("bcr_abl_rindels", {}).get("mem_per_cpu", config.get("rindels", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])),
        partition=config.get("bcr_abl_rindels", {}).get("partition", config.get("rindels", {}).get("partition", config["default_resources"]["partition"])),
        threads=config.get("bcr_abl_rindels", {}).get("threads", config.get("rindels", {}).get("threads", 4)),
        time=config.get("bcr_abl_rindels", {}).get("time", config.get("rindels", {}).get("time", config["default_resources"]["time"])),
    message:
        "{rule}: run Rindels in BCR-ABL1 mode for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        python3 {workflow.basedir}/scripts/run_rindels.py bcr-abl \
            --binary={params.binary:q} \
            --bam {input.bam} \
            --ref {input.ref} \
            --vcf {output.vcf} \
            --haplotypes {output.haplotypes} \
            --genome {params.genome} \
            --extra={params.extra:q} \
            --threads {threads} > {log} 2>&1
        """


rule bcr_abl_vardictrust:
    """Run VarDictRust on BCR-ABL1 target exons with ONT-tuned parameters."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        ref=config["reference"]["fasta"],
        bed=lambda wildcards: _bcr_abl_exons_bed(),
    output:
        vcf=temp("snv_indels/vardictrust/{sample}_{type}.bcr_abl.raw.vcf"),
    params:
        min_af=lambda wildcards: get_assay_setting(
            wildcards.sample,
            "bcr_abl",
            "vardict_min_af",
            config.get("vardictrust", {}).get(
                "min_af",
                config.get("bcr_abl", {}).get("vardict_min_af", 0.15),
            ),
        ),
        min_baseq=lambda wildcards: config.get("vardictrust", {}).get("min_base_quality", 12.0),
        max_mismatches=lambda wildcards: config.get("vardictrust", {}).get("max_mismatches", 150),
        col_chr=config.get("vardictrust", {}).get("col_chr", 1),
        col_start=config.get("vardictrust", {}).get("col_start", 2),
        col_end=config.get("vardictrust", {}).get("col_end", 3),
        extra=lambda wildcards: config.get("vardictrust", {}).get("extra", ""),
    log:
        "logs/bcr_abl_vardictrust/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcr_abl_vardictrust/{sample}_{type}.tsv",
            config.get("bcr_abl_vardictrust", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("vardictrust", {}).get("threads", 4)
    resources:
        mem_mb=config.get("vardictrust", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vardictrust", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vardictrust", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vardictrust", {}).get("threads", 4),
        time=config.get("vardictrust", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vardictrust", {}).get("container", "docker://quay.io/biocontainers/vardictrust:1.0.2--hfa8f182_0")
    message:
        "{rule}: call BCR-ABL1 target exons with VarDictRust for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        exec 2> {log}
        vardictrust \
            -G {input.ref} \
            -b {input.bam} \
            -N {wildcards.sample}_{wildcards.type} \
            -c {params.col_chr} -S {params.col_start} -E {params.col_end} \
            -f {params.min_af} \
            -q {params.min_baseq} \
            -m {params.max_mismatches} \
            -U \
            --fisher \
            --threads {threads} \
            {params.extra} \
            {input.bed} > {output.vcf}
        """


rule bcr_abl_vardictrust_sort_index:
    """Sort, compress and index VarDictRust VCF."""
    input:
        vcf="snv_indels/vardictrust/{sample}_{type}.bcr_abl.raw.vcf",
    output:
        vcf="snv_indels/vardictrust/{sample}_{type}.bcr_abl.vcf.gz",
        tbi="snv_indels/vardictrust/{sample}_{type}.bcr_abl.vcf.gz.tbi",
    log:
        "logs/bcr_abl_vardictrust_sort_index/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/bcr_abl_vardictrust_sort_index/{sample}_{type}.tsv",
            config.get("bcr_abl_vardictrust_sort_index", {}).get("benchmark_repeats", 1),
        )
    threads: 2
    resources:
        mem_mb=config["default_resources"]["mem_mb"],
        mem_per_cpu=config["default_resources"]["mem_per_cpu"],
        partition=config["default_resources"]["partition"],
        threads=2,
        time=config["default_resources"]["time"],
    container:
        config.get("bcftools", {}).get("container", "docker://hydragenetics/common:3.1.1")
    message:
        "{rule}: sort and index VarDictRust VCF for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        exec > {log} 2>&1
        bcftools sort -Oz -o {output.vcf} {input.vcf}
        bcftools index -t {output.vcf}
        """







