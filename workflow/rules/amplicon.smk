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


def _amplicon_rindels_binary():
    """Rindels binary for the amplicon rules.

    Deliberately named apart from bcr_abl.smk's ``_find_rindels_binary``: both
    files are included into the same namespace, so a second definition of that
    name would silently replace the first for every rule in the workflow.
    """
    override = config.get("amplicon", {}).get("rindels_binary", "")
    if override and (os.path.isfile(override) or shutil.which(override)):
        return override
    return _find_rindels_binary()


def _amplicon_top_n():
    """Reads kept per amplicon before calling; 0 disables downsampling."""
    return int(config.get("amplicon", {}).get("top_n", 0) or 0)


def get_amplicon_bam_for_calling(wildcards):
    """Return downsampled top-N BAM if top_n > 0, else full merged BAM."""
    if _amplicon_top_n() > 0:
        return f"alignment/amplicon/{wildcards.sample}_{wildcards.type}.topn.bam"
    return f"alignment/merged/{wildcards.sample}_{wildcards.type}.bam"


def get_amplicon_bai_for_calling(wildcards):
    """Return downsampled top-N BAI if top_n > 0, else full merged BAI."""
    if _amplicon_top_n() > 0:
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
    # Deliberately no container by default: this rule runs the nanostream binary
    # resolved by get_nanostream_binary(), and a container would hide a binary
    # found on the host PATH. Set amplicon.nanostream_container together with a
    # binary path inside the image to containerise it.
    container:
        config.get("amplicon", {}).get("nanostream_container", None)
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
        top_n=lambda wildcards: _amplicon_top_n(),
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
    # See amplicon_nanostream_stats: host binary, so no container by default.
    container:
        config.get("amplicon", {}).get("nanostream_container", None)
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
    container:
        config.get("amplicon", {}).get("container", config["default_container"])
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
    """Call targeted SNVs, indels, and insertions inside amplicon target regions with Rindels.

    FLT3-ITD is read out of this VCF by amplicon_summary; the ``rindels itd``
    subcommand is a short-read caller and is not used on ONT data.
    """
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
        binary=lambda wildcards: _amplicon_rindels_binary(),
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
    # Deliberately no container by default: Rindels is a host binary, as in
    # bcr_abl_rindels. Set amplicon.rindels_container plus amplicon.rindels_binary
    # together to run it from an image.
    container:
        config.get("amplicon", {}).get("rindels_container", None)
    shell:
        """
        exec > {log} 2>&1
        RINDELS_BIN="{params.binary}"
        if ! command -v "$RINDELS_BIN" >/dev/null 2>&1 && [ ! -x "$RINDELS_BIN" ]; then
            echo "Error: Rindels binary not found at $RINDELS_BIN"
            exit 1
        fi
        # Record which binary produced this call set: the path is resolved from
        # the host, so the version is not otherwise pinned anywhere.
        echo "Rindels binary: $RINDELS_BIN"
        "$RINDELS_BIN" --version || true

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
            bcftools index -t {output.vcf}
        fi
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
    threads: 1
    resources:
        mem_mb=config.get("amplicon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("amplicon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("amplicon", {}).get("container", config["default_container"])
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
    """Compile the clinical summary: amplicon QC, Rindels variants, derived FLT3-ITD.

    The FLT3-ITD table is derived here from insertions inside the FLT3 exon
    13-15 target of the long-read VCF. When no VCF is available the table is
    written with status NOT_ASSESSED, so a missing call is never delivered as a
    negative result.
    """
    input:
        stats="reports/amplicon/{sample}_{type}.amplicon_stats.json",
        vcf="snv_indels/rindels/{sample}_{type}.amplicon.vcf.gz",
        pileup="reports/amplicon/{sample}_{type}.pileup_counts.tsv",
        pools=_amplicon_pools_yaml(),
    output:
        tsv="reports/amplicon/{sample}_{type}.amplicon_summary.tsv",
        json="reports/amplicon/{sample}_{type}.amplicon_summary.json",
        itd="reports/amplicon/{sample}_{type}.flt3_itd.tsv",
    params:
        script=f"{workflow.basedir}/scripts/amplicon_report.py",
        itd_region=lambda wildcards: config.get("amplicon", {}).get("flt3_itd_region", ""),
        itd_min_length=lambda wildcards: config.get("amplicon", {}).get("flt3_itd_min_length", 6),
        itd_min_vaf=lambda wildcards: config.get("amplicon", {}).get("flt3_itd_min_vaf", 1.0),
    log:
        "logs/amplicon_summary/{sample}_{type}.log",
    threads: 1
    resources:
        mem_mb=config.get("amplicon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("amplicon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("amplicon", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        exec > {log} 2>&1
        REGION_ARG=""
        if [ -n "{params.itd_region}" ]; then
            REGION_ARG="--flt3-itd-region {params.itd_region}"
        fi
        python3 {params.script} \
            --json {input.stats} \
            --vcf {input.vcf} \
            --pileup {input.pileup} \
            --pools {input.pools} \
            --sample {wildcards.sample} \
            $REGION_ARG \
            --flt3-itd-min-length {params.itd_min_length} \
            --flt3-itd-min-vaf {params.itd_min_vaf} \
            --output-itd-tsv {output.itd} \
            --output-summary-tsv {output.tsv} \
            --output-summary-json {output.json}
        """

# ===========================================================================
# neville_mx_amplicon comparison report
# ===========================================================================
# neville_mx_amplicon runs the same 15 amplicons with a VarDict/ClairS-TO/
# DeepSomatic ensemble. These rules render nagini's Rindels calls into
# neville's table shape so the two can be compared sample by sample.
#
# The VEP flags below are neville's, deliberately: a different --pick_order
# selects a different transcript, and the same variant would then show a
# different HGVSc/ENSP in the two reports - a difference in annotation config
# that reads as a difference in calling.


def _amplicon_run_id(sample, sample_type):
    """RunID for the neville RunID column, from units.run_id."""
    configured = str(config.get("amplicon", {}).get("run_id", "") or "").strip()
    if configured:
        return configured
    rows = _unit_rows(sample, sample_type)
    if not rows.empty:
        run_ids = [str(r).strip() for r in rows.get("run_id", pd.Series(dtype=str)).tolist() if str(r).strip()]
        if run_ids:
            return run_ids[0]
    return ""


rule amplicon_norm_rindels:
    """Left-align and split the amplicon calls before annotation.

    neville normalises with vt before VEP. Without the same step an indel can
    sit at a different position in the two reports and look like a discrepancy.
    """
    input:
        vcf="snv_indels/rindels/{sample}_{type}.amplicon.vcf.gz",
        tbi="snv_indels/rindels/{sample}_{type}.amplicon.vcf.gz.tbi",
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("snv_indels/rindels/{sample}_{type}.amplicon.norm.vcf.gz"),
        tbi=temp("snv_indels/rindels/{sample}_{type}.amplicon.norm.vcf.gz.tbi"),
    params:
        extra=config.get("amplicon_norm_rindels", {}).get("extra", "--multiallelics -both --check-ref w"),
    log:
        "logs/amplicon_norm_rindels/{sample}_{type}.log",
    threads: 1
    resources:
        mem_mb=config.get("amplicon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("amplicon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("amplicon", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: normalize amplicon calls for {wildcards.sample}_{wildcards.type}"
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


rule amplicon_vep_annotate:
    """Annotate the amplicon calls with neville's VEP settings."""
    input:
        vcf="snv_indels/rindels/{sample}_{type}.amplicon.norm.vcf.gz",
        tbi="snv_indels/rindels/{sample}_{type}.amplicon.norm.vcf.gz.tbi",
        fasta=config["reference"]["fasta"],
        cache=lambda wildcards: config["vep"]["vep_cache"],
    output:
        vcf=temp("snv_indels/rindels/{sample}_{type}.amplicon.norm.vep.vcf"),
        stats="reports/neville_compat/{sample}_{type}.vep.summary.html",
    params:
        mode=config.get("amplicon_vep", {}).get("mode", config.get("vep", {}).get("mode", "--offline --cache --merged")),
        extra=config.get("amplicon_vep", {}).get(
            "extra",
            "--format vcf --assembly GRCh38 --check_existing --pick --variant_class --everything "
            "--pick_order mane_select,mane_plus_clinical,canonical,biotype,rank,appris,tsl,ccds,length,ensembl,refseq",
        ),
    log:
        "logs/amplicon_vep_annotate/{sample}_{type}.log",
    threads: config.get("amplicon_vep", {}).get("threads", config.get("vep", {}).get("threads", 4))
    resources:
        mem_mb=config.get("amplicon_vep", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon_vep", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon_vep", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("amplicon_vep", {}).get("threads", config.get("vep", {}).get("threads", 4)),
        time=config.get("amplicon_vep", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("amplicon_vep", {}).get("container", config.get("vep", {}).get("container", config["default_container"]))
    conda:
        "../../envs/vep.yaml"
    message:
        "{rule}: VEP annotate amplicon calls for {wildcards.sample}_{wildcards.type}"
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
            {params.extra} &> {log}
        """


rule amplicon_vep_bgzip:
    """Compress and index the VEP-annotated amplicon calls."""
    input:
        vcf="snv_indels/rindels/{sample}_{type}.amplicon.norm.vep.vcf",
    output:
        vcf="snv_indels/rindels/{sample}_{type}.amplicon.norm.vep.vcf.gz",
        tbi="snv_indels/rindels/{sample}_{type}.amplicon.norm.vep.vcf.gz.tbi",
    log:
        "logs/amplicon_vep_bgzip/{sample}_{type}.log",
    threads: 1
    resources:
        mem_mb=config.get("amplicon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("amplicon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("amplicon", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        bcftools view --output-type z --output {output.vcf} {input.vcf} &> {log}
        bcftools index --force --tbi {output.vcf} >> {log} 2>&1
        """


rule amplicon_neville_compat_report:
    """Render the amplicon calls in the neville_mx_amplicon table layout."""
    input:
        vcf="snv_indels/rindels/{sample}_{type}.amplicon.norm.vep.vcf.gz",
        tbi="snv_indels/rindels/{sample}_{type}.amplicon.norm.vep.vcf.gz.tbi",
    output:
        snv_tsv="reports/neville_compat/{sample}_{type}.snv_indels.tsv",
        sv_tsv="reports/neville_compat/{sample}_{type}.svs.tsv",
        xlsx="reports/neville_compat/{sample}_{type}.neville_layout.xlsx",
    params:
        script=f"{workflow.basedir}/scripts/neville_compat_report.py",
        run_id=lambda wildcards: _amplicon_run_id(wildcards.sample, wildcards.type),
        reference=config["reference"]["fasta"],
        callers=config.get("amplicon", {}).get("neville_callers_label", "rindels"),
        min_af=config.get("amplicon", {}).get("neville_min_af", 0.005),
        min_sv_length=config.get("amplicon", {}).get("neville_min_sv_length", 50),
    log:
        "logs/amplicon_neville_compat_report/{sample}_{type}.log",
    threads: 1
    resources:
        mem_mb=config.get("amplicon", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("amplicon", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("amplicon", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("amplicon", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("amplicon", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: render {wildcards.sample}_{wildcards.type} in the neville report layout"
    shell:
        """
        exec > {log} 2>&1
        python3 {params.script} \
            --vcf {input.vcf} \
            --sample {wildcards.sample} \
            --run-id "{params.run_id}" \
            --reference {params.reference} \
            --callers "{params.callers}" \
            --min-af {params.min_af} \
            --min-sv-length {params.min_sv_length} \
            --output-snv-tsv {output.snv_tsv} \
            --output-sv-tsv {output.sv_tsv} \
            --output-xlsx {output.xlsx}
        """
