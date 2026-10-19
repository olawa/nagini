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
        # Rindels has no assembly option: the unified PDF's ideogram gaps are
        # hg38 whatever reference it was given, and `rindels plot --ref-genome`
        # says as much ("hg38 only supported for ideogram gaps currently").
        # The report is therefore switchable, so a T2T run can drop a figure
        # drawn on the wrong coordinates instead of publishing it.
        report=lambda wildcards: "--report" if config.get("rindels", {}).get("report", True) else "",
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
        python3 {workflow.basedir}/scripts/run_rindels.py cnv \
            --binary={params.binary:q} \
            --bam {input.bam} \
            --ref {input.ref} \
            --bins {output.bins} \
            --bin-size {params.bin_size} \
            --vcf {output.vcf} \
            --ascat {output.ascat} \
            --pdf {output.pdf} \
            {params.report} \
            --annotation={params.annotation:q} \
            --extra={params.extra:q} \
            --threads {threads} > {log} 2>&1
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
        gaps=lambda wildcards: get_gaps_bed() or [],
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
        python3 {workflow.basedir}/scripts/run_rindels.py ideogram \
            --binary={params.binary:q} \
            --bins {input.bins} \
            --parquet={params.parquet:q} \
            --baf={params.baf:q} \
            --centromeres={input.centromeres:q} \
            --gaps={input.gaps:q} \
            --bin-size {params.bin_size} \
            --width {params.width} \
            --height {params.height} \
            --extra={params.extra:q} \
            --output {output.png} > {log} 2>&1
        """
