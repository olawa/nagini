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
        pdf="cnv/rindels/{sample}_{type}_report.pdf",
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
    threads: config.get("cnv_rindels", {}).get("threads", config.get("rindels", {}).get("threads", 4))
    resources:
        mem_mb=config.get("cnv_rindels", {}).get("mem_mb", config.get("rindels", {}).get("mem_mb", config["default_resources"]["mem_mb"])),
        mem_per_cpu=config.get("cnv_rindels", {}).get("mem_per_cpu", config.get("rindels", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"])),
        partition=config.get("cnv_rindels", {}).get("partition", config.get("rindels", {}).get("partition", config["default_resources"]["partition"])),
        threads=config.get("cnv_rindels", {}).get("threads", config.get("rindels", {}).get("threads", 4)),
        time=config.get("cnv_rindels", {}).get("time", config.get("rindels", {}).get("time", config["default_resources"]["time"])),
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
            --report \
            --report-path {output.pdf} \
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

        if [ ! -f "{output.pdf}" ]; then
            echo "Warning: Rindels unified PDF report was not generated; creating empty placeholder."
            touch "{output.pdf}"
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


rule cytoband_centromeres:
    """Centromere intervals for the ideogram, from the build's own cytobands.

    Rindels falls back to a built-in hg38 set when it is given none, which on
    T2T puts the centromeres in the wrong places.
    """
    output:
        bed="references/centromeres_{build}.bed",
    params:
        cytoband=lambda wildcards: get_cytoband_file() or "",
    log:
        "logs/cytoband_centromeres/centromeres_{build}.log",
    threads: 1
    resources:
        mem_mb=config["default_resources"]["mem_mb"],
        mem_per_cpu=config["default_resources"]["mem_per_cpu"],
        partition=config["default_resources"]["partition"],
        threads=1,
        time=config["default_resources"]["time"],
    container:
        config["default_container"]
    conda:
        "../../envs/nagini.yaml"
    wildcard_constraints:
        build=r"[A-Za-z0-9._-]+",
    message:
        "{rule}: derive centromere intervals from the cytoband table"
    script:
        "../scripts/centromeres_from_cytoband.py"


rule cnv_rindels_ideogram:
    """Whole-genome ideogram from the Rindels bins and BAF track.

    The CNV data comes from the same Parquet table the unified PDF report
    renders from, not from cnv_bins.bed, and it is aggregated at the run's own
    bin size rather than the subcommand's 1 Mb default. Drawing the BED at 1 Mb
    while the PDF drew the Parquet at 50 kb is what made the two pictures of
    one sample disagree.

    Both the Parquet and the BAF track are sidecars of cnv_rindels rather than
    declared outputs, so each is used only when it exists and cnv_bins.bed,
    which is declared, remains the dependency.
    """
    input:
        bins="cnv/rindels/{sample}_{type}.cnv_bins.bed",
        centromeres=lambda wildcards: "references/centromeres_{}.bed".format(
            config.get("reference_build", "hg38")
        ),
    output:
        png="reports/cnv_ideogram/{sample}_{type}.ideogram.png",
    params:
        binary=lambda wildcards: _find_rindels_binary(),
        baf=lambda wildcards: f"cnv/rindels/{wildcards.sample}_{wildcards.type}.baf.bed",
        # `-o X.vcf` makes Rindels write its CNV bin table to X.vcf.parquet.
        parquet=lambda wildcards: f"cnv/rindels/{wildcards.sample}_{wildcards.type}.vcf.parquet",
        bin_size=config.get("rindels", {}).get("cnv_bin_size", 50000),
        width=config.get("rindels_ideogram", {}).get("width", 2000),
        height=config.get("rindels_ideogram", {}).get("height", 2400),
        extra=config.get("rindels_ideogram", {}).get("extra", ""),
    log:
        "logs/cnv_rindels_ideogram/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/cnv_rindels_ideogram/{sample}_{type}.tsv",
            config.get("rindels_ideogram", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("rindels_ideogram", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("rindels_ideogram", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("rindels_ideogram", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("rindels_ideogram", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: draw the whole-genome ideogram for {wildcards.sample}_{wildcards.type}"
    shell:
        """
        exec > {log} 2>&1
        RINDELS_BIN="{params.binary}"
        if ! command -v "$RINDELS_BIN" >/dev/null 2>&1 && [ ! -x "$RINDELS_BIN" ]; then
            echo "Error: Rindels binary not found at $RINDELS_BIN"
            exit 1
        fi

        mkdir -p "$(dirname {output.png})"

        # Prefer the table the PDF report renders from; fall back to the BED.
        if [ -s "{params.parquet}" ]; then
            CNV_ARG="--cnv {params.parquet}"
        else
            echo "No CNV Parquet at {params.parquet}; falling back to the binned BED."
            CNV_ARG="--bed {input.bins}"
        fi

        BAF_ARG=""
        if [ -s "{params.baf}" ]; then
            BAF_ARG="--baf-bed {params.baf}"
        else
            echo "No BAF track at {params.baf}; drawing density only."
        fi

        CENTROMERE_ARG=""
        if [ -s "{input.centromeres}" ]; then
            CENTROMERE_ARG="--centromeres {input.centromeres}"
        else
            echo "No centromere BED; Rindels will fall back to its built-in hg38 set."
        fi

        "$RINDELS_BIN" ideogram \
            $CNV_ARG \
            $BAF_ARG \
            $CENTROMERE_ARG \
            --bin-size {params.bin_size} \
            --res $(( {params.bin_size} / 1000 )) \
            --width {params.width} \
            --height {params.height} \
            {params.extra} \
            -o {output.png}
        """
