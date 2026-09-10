__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule cnv_rindels:
    """Run Rindels in CNV-only mode to produce whole-genome binned read counts and ASCAT segments."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        ref=config["reference"]["fasta"],
        gtf=lambda wildcards: get_gene_annotation() or [],
    output:
        # Rindels names its sidecars after the -o base, so -o must be
        # {sample}_{type}.vcf for the segmentation to land on
        # {sample}_{type}.ascat.bed. Passing -o *.vcf.gz would let rindels
        # write BGZF and its own .tbi, but it also changes how it names the
        # sidecars -- with a .gz suffix it keeps the full name for some
        # (X.vcf.gz.report_baf.parquet) and strips it for others (X.baf.bed) --
        # and .ascat.bed is the one file cnv_summary cannot do without.  The
        # plain VCF is compressed by cnv_rindels_bgzip instead, which runs in
        # a container that has htslib; this rule runs on the host, because the
        # rindels binary is a host binary, so it must not depend on bgzip or
        # tabix being installed there. The parquet sidecars stay undeclared:
        # which ones appear depends on the mode and the data, as
        # sv_added_snpengine_del.tsv appearing for only three of five samples
        # shows, and declaring an output that may not appear is what made this
        # rule fail in the first place.
        bins="cnv/rindels/{sample}_{type}.cnv_bins.bed",
        ascat="cnv/rindels/{sample}_{type}.ascat.bed",
        vcf=temp("cnv/rindels/{sample}_{type}.vcf"),
    params:
        binary=lambda wildcards: _find_rindels_binary(),
        bin_size=config.get("rindels", {}).get("cnv_bin_size", 50000),
        # Keep this in step with config.yaml's rindels.extra_cnv: --only-cnv
        # implies --no-snp, which leaves ASCAT without heterozygous markers and
        # produces no .ascat.bed at all.
        extra=config.get("rindels", {}).get("extra_cnv", "--ont --sv-lite"),
        # Gene annotation when the reference build has a GTF installed, so the
        # call set carries genes and consequences rather than coordinates only.
        annotation=lambda wildcards: get_rindels_annotation_args(),
    log:
        "logs/cnv_rindels/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/cnv_rindels/{sample}_{type}.tsv",
            config.get("cnv_rindels", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("rindels", {}).get("threads", 4)
    resources:
        mem_mb=config.get("rindels", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("rindels", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("rindels", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("rindels", {}).get("threads", 4),
        time=config.get("rindels", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: call whole-genome copy number variations with Rindels for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        exec > {log} 2>&1
        RINDELS_BIN="{params.binary}"
        if ! command -v "$RINDELS_BIN" >/dev/null 2>&1 && [ ! -x "$RINDELS_BIN" ]; then
            echo "Error: Rindels binary not found at $RINDELS_BIN"
            exit 1
        fi

        "$RINDELS_BIN" \
            -b {input.bam} \
            -r {input.ref} \
            --cnv-output {output.bins} \
            --cnv-bin-size {params.bin_size} \
            -o {output.vcf} \
            -t {threads} \
            {params.annotation} \
            {params.extra}

        # Rindels skips the ASCAT segmentation when it finds fewer than 100
        # heterozygous sites, which is a legitimate outcome for a shallow or
        # genomically quiet sample. Leave a header-only file so the run
        # reports no segments instead of dying on a missing output;
        # cnv_summary already treats an empty file as zero segments.
        if [ ! -f "{output.ascat}" ]; then
            echo "ASCAT segmentation skipped by Rindels (too few heterozygous sites); see the log."
            printf '%s\n' \
                "# ASCAT-style Heuristic Copy Number Segments" \
                "# NOTE: Rindels skipped segmentation for this sample (too few heterozygous sites)" \
                "# Columns: chrom, start, end, CN_A, CN_B, BAF, logR, n_sites, state, n_het, n_hom, n_total" \
                > "{output.ascat}"
        fi
        """


rule cnv_summary:
    """Summarise Rindels ASCAT CNV segments and annotate with cytobands and target regions."""
    input:
        ascat="cnv/rindels/{sample}_{type}.ascat.bed",
        target_bed=lambda wildcards: get_sample_design_bed(wildcards.sample) or [],
        cytoband=lambda wildcards: get_cytoband_file(wildcards) or [],
    output:
        tsv="reports/cnv_summary/{sample}_{type}.cnv.tsv",
        json="reports/cnv_summary/{sample}_{type}.cnv.json",
    params:
        extra=config.get("cnv_summary", {}).get("extra", ""),
    log:
        "logs/cnv_summary/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/cnv_summary/{sample}_{type}.tsv",
            config.get("cnv_summary", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("cnv_summary", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnv_summary", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnv_summary", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("cnv_summary", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnv_summary", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: summarise copy number alterations for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/cnv_summary.py"
