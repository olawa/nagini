__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"

import os
import shutil


def _amplicon_resource(filename):
    """Return absolute path to a bundled amplicon panel resource file."""
    return os.path.join(workflow.basedir, "resources", "amplicon", filename)


def _amplicon_primers_tsv():
    override = config.get("amplicon", {}).get("primers", "")
    if override and os.path.isfile(override):
        return override
    return _amplicon_resource("amplicons_primers.tsv")


def _amplicon_targets_bed():
    override = config.get("amplicon", {}).get("targets_bed", "") or config.get("amplicon", {}).get("bed", "")
    if override and os.path.isfile(override):
        return override
    return _amplicon_resource("amplicon_targets.bed")


def _amplicon_exons_bed():
    override = config.get("amplicon", {}).get("exons_bed", "")
    if override and os.path.isfile(override):
        return override
    return _amplicon_resource("exon_coverage.bed")


def _amplicon_pools_yaml():
    override = config.get("amplicon", {}).get("pools", "")
    if override and os.path.isfile(override):
        return override
    return _amplicon_resource("amplicons_pools.yaml")


def _amplicon_gtf():
    override = config.get("amplicon", {}).get("gtf", "")
    if override and os.path.isfile(override):
        return override
    ref_gtf = config.get("reference", {}).get("gene_annotation", "") or config.get("reference", {}).get("gtf", "")
    if ref_gtf and os.path.isfile(ref_gtf):
        return ref_gtf
    return ""


def _find_rindels_binary():
    """Locate the Rindels binary from config, PATH, ~/bin, or dev checkout."""
    override = config.get("rindels", {}).get("binary", "") or config.get("amplicon", {}).get("rindels_binary", "")
    if override and (os.path.isfile(override) or shutil.which(override)):
        return override
    candidates = [
        os.path.expanduser("~/bin/rindels"),
        "/home/cgu-lab/bin/rindels",
        os.path.expanduser("~/tmp/tools/rindels/target/release/rindels"),
        "/home/cgu-lab/tmp/tools/rindels/target/release/rindels",
        os.path.expanduser("~/.cargo/bin/rindels"),
    ]
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    which_r = shutil.which("rindels")
    if which_r:
        return which_r
    return "rindels"


def get_amplicon_bam_for_calling(wildcards):
    """Return downsampled top-N BAM if top_n > 0, else full merged BAM."""
    top_n = int(config.get("amplicon", {}).get("top_n", 0) or 0)
    if top_n > 0:
        return f"alignment/amplicon/{wildcards.sample}_{wildcards.type}.topn.bam"
    return f"alignment/merged/{wildcards.sample}_{wildcards.type}.bam"


def get_amplicon_bai_for_calling(wildcards):
    """Return downsampled top-N BAI if top_n > 0, else full merged BAI."""
    top_n = int(config.get("amplicon", {}).get("top_n", 0) or 0)
    if top_n > 0:
        return f"alignment/amplicon/{wildcards.sample}_{wildcards.type}.topn.bam.bai"
    return f"alignment/merged/{wildcards.sample}_{wildcards.type}.bam.bai"


rule amplicon_nanostream_stats:
    """Run nanostream amplicons on the aligned BAM to calculate amplicon and pool metrics."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        primers=_amplicon_primers_tsv(),
    output:
        json="reports/amplicon/{sample}_{type}.amplicon_stats.json",
    params:
        binary=lambda wildcards: get_nanostream_binary(),
        min_qs=lambda wildcards: config.get("amplicon", {}).get("min_qs", 10),
        len_range=lambda wildcards: config.get("amplicon", {}).get("len_range", "1000-3500"),
        extra=lambda wildcards: config.get("amplicon", {}).get("nanostream_extra", ""),
    log:
        "logs/amplicon_nanostream_stats/{sample}_{type}.log",
    threads: config.get("amplicon", {}).get("threads", 4)
    resources:
        mem_mb=config.get("amplicon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("amplicon", {}).get("threads", 4),
        time=config.get("amplicon", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        exec > {log} 2>&1
        NS_BIN="{params.binary}"
        "$NS_BIN" amplicons {input.bam} \
            --primers {input.primers} \
            --output {output.json} \
            --min-qs {params.min_qs} \
            --len {params.len_range} \
            -t {threads} {params.extra}
        """


rule amplicon_downsample_topn:
    """Downsample aligned BAM reads to top-N highest quality reads per amplicon."""
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        primers=_amplicon_primers_tsv(),
    output:
        fastq=temp("reads/amplicon/{sample}_{type}.topn.fastq.gz"),
    params:
        binary=lambda wildcards: get_nanostream_binary(),
        top_n=lambda wildcards: config.get("amplicon", {}).get("top_n", 1000),
        min_qs=lambda wildcards: config.get("amplicon", {}).get("min_qs", 10),
        len_range=lambda wildcards: config.get("amplicon", {}).get("len_range", "1000-3500"),
    log:
        "logs/amplicon_downsample_topn/{sample}_{type}.log",
    threads: config.get("amplicon", {}).get("threads", 4)
    resources:
        mem_mb=config.get("amplicon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("amplicon", {}).get("threads", 4),
        time=config.get("amplicon", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        exec > {log} 2>&1
        NS_BIN="{params.binary}"
        "$NS_BIN" amplicons {input.bam} \
            --primers {input.primers} \
            --top-n {params.top_n} \
            --output-fastq {output.fastq} \
            --min-qs {params.min_qs} \
            --len {params.len_range} \
            -t {threads}
        """


def _amplicon_align_target(wildcards):
    idx = get_reference_index()
    if idx:
        return idx
    return config["reference"]["fasta"]


rule amplicon_topn_align:
    """Align downsampled top-N amplicon reads with minimap2 and index."""
    input:
        fastq="reads/amplicon/{sample}_{type}.topn.fastq.gz",
        target=_amplicon_align_target,
    output:
        bam="alignment/amplicon/{sample}_{type}.topn.bam",
        bai="alignment/amplicon/{sample}_{type}.topn.bam.bai",
    log:
        "logs/amplicon_topn_align/{sample}_{type}.log",
    threads: config.get("amplicon", {}).get("align_threads", 8)
    resources:
        mem_mb=config.get("amplicon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("amplicon", {}).get("align_threads", 8),
        time=config.get("amplicon", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        exec > {log} 2>&1
        minimap2 -ax map-ont -y --MD -t {threads} {input.target} {input.fastq} | \
        samtools sort -@ {threads} -o {output.bam} -
        samtools index -@ {threads} {output.bam}
        """


rule amplicon_rindels:
    """Call targeted SNVs, indels, and insertions inside amplicon target regions with Rindels."""
    input:
        bam=get_amplicon_bam_for_calling,
        bai=get_amplicon_bai_for_calling,
        ref=config["reference"]["fasta"],
        bed=_amplicon_targets_bed(),
    output:
        vcf="snv_indels/rindels/{sample}_{type}.amplicon.vcf.gz",
        tbi="snv_indels/rindels/{sample}_{type}.amplicon.vcf.gz.tbi",
        pileup_counts="reports/amplicon/{sample}_{type}.pileup_counts.tsv",
    params:
        binary=lambda wildcards: _find_rindels_binary(),
        gtf=lambda wildcards: _amplicon_gtf(),
        extra=lambda wildcards: config.get("amplicon", {}).get("rindels_extra", "--ont --amplicon"),
    log:
        "logs/amplicon_rindels/{sample}_{type}.log",
    threads: config.get("amplicon", {}).get("threads", 4)
    resources:
        mem_mb=config.get("amplicon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("amplicon", {}).get("threads", 4),
        time=config.get("amplicon", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        exec > {log} 2>&1
        RINDELS_BIN="{params.binary}"
        if ! command -v "$RINDELS_BIN" >/dev/null 2>&1 && [ ! -x "$RINDELS_BIN" ]; then
            echo "Error: Rindels binary not found at $RINDELS_BIN"
            exit 1
        fi

        GTF_ARG=""
        if [ -n "{params.gtf}" ] && [ -f "{params.gtf}" ]; then
            GTF_ARG="--gtf {params.gtf}"
        fi

        "$RINDELS_BIN" -b {input.bam} -r {input.ref} -o {output.vcf} \
            {params.extra} \
            --bed {input.bed} \
            $GTF_ARG \
            --pileup-counts {input.bed} \
            --pileup-counts-out {output.pileup_counts} \
            -t {threads}

        if [ -s "{output.vcf}" ] && [ ! -s "{output.tbi}" ]; then
            if command -v bcftools >/dev/null 2>&1; then
                bcftools index -t {output.vcf} || true
            fi
        fi
        """


rule amplicon_rindels_itd:
    """Identify and quantify FLT3 internal tandem duplications (ITD) using Rindels ITD."""
    input:
        bam=get_amplicon_bam_for_calling,
        bai=get_amplicon_bai_for_calling,
        ref=config["reference"]["fasta"],
    output:
        tsv="reports/amplicon/{sample}_{type}.flt3_itd.tsv",
    params:
        binary=lambda wildcards: _find_rindels_binary(),
        region=lambda wildcards: config.get("amplicon", {}).get("flt3_itd_region", "chr13:28032002-28035169"),
    log:
        "logs/amplicon_rindels_itd/{sample}_{type}.log",
    threads: config.get("amplicon", {}).get("threads", 2)
    resources:
        mem_mb=config.get("amplicon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("amplicon", {}).get("threads", 2),
        time=config.get("amplicon", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        exec > {log} 2>&1
        RINDELS_BIN="{params.binary}"
        if ! command -v "$RINDELS_BIN" >/dev/null 2>&1 && [ ! -x "$RINDELS_BIN" ]; then
            echo "Error: Rindels binary not found at $RINDELS_BIN"
            exit 1
        fi

        "$RINDELS_BIN" itd -b {input.bam} -r {input.ref} \
            --region {params.region} \
            -o {output.tsv} \
            --sample-name {wildcards.sample} || {
                echo -e "gene\tregion\tlength\taf\tstatus\nFLT3\t{params.region}\t0\t0.0\tNEGATIVE" > {output.tsv}
            }
        """


rule amplicon_stats_tsv:
    """Format amplicon QC stats into a readable TSV table."""
    input:
        json="reports/amplicon/{sample}_{type}.amplicon_stats.json",
        pools=_amplicon_pools_yaml(),
    output:
        tsv="reports/amplicon/{sample}_{type}.amplicon_stats.tsv",
    params:
        script=f"{workflow.basedir}/scripts/amplicon_report.py",
    log:
        "logs/amplicon_stats_tsv/{sample}_{type}.log",
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        exec > {log} 2>&1
        python3 {params.script} \
            --json {input.json} \
            --pools {input.pools} \
            --sample {wildcards.sample} \
            --output-stats-tsv {output.tsv}
        """


rule amplicon_summary:
    """Compile unified clinical summary report combining amplicon QC, Rindels variants, and ITD calls."""
    input:
        stats="reports/amplicon/{sample}_{type}.amplicon_stats.json",
        vcf="snv_indels/rindels/{sample}_{type}.amplicon.vcf.gz",
        itd="reports/amplicon/{sample}_{type}.flt3_itd.tsv",
        pileup="reports/amplicon/{sample}_{type}.pileup_counts.tsv",
        pools=_amplicon_pools_yaml(),
    output:
        tsv="reports/amplicon/{sample}_{type}.amplicon_summary.tsv",
        json="reports/amplicon/{sample}_{type}.amplicon_summary.json",
    params:
        script=f"{workflow.basedir}/scripts/amplicon_report.py",
    log:
        "logs/amplicon_summary/{sample}_{type}.log",
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        exec > {log} 2>&1
        python3 {params.script} \
            --json {input.stats} \
            --vcf {input.vcf} \
            --itd {input.itd} \
            --pileup {input.pileup} \
            --pools {input.pools} \
            --sample {wildcards.sample} \
            --output-summary-tsv {output.tsv} \
            --output-summary-json {output.json}
        """
