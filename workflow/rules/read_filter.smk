__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"

import os
import platform
import shutil


def get_nanostream_binary():
    """Resolve the nanostream executable, or fail with what to do about it.

    There is deliberately no Python fallback. A fallback that reimplements the
    filter is worse than no filter at all: it computed mean read quality as the
    arithmetic mean of Phred scores, which overestimates quality on exactly the
    ONT failure mode this filter exists to catch (a good read with a bad tail),
    so the same config silently kept different reads depending on which
    implementation happened to be resolved on that host.

    nanostream computes it as -10*log10(mean(10^(-q/10))) on the FASTQ path and
    reads the qs tag directly on the BAM path.
    """

    configured = config.get("read_filter", {}).get("binary", "")
    if configured:
        expanded = os.path.expanduser(os.path.expandvars(configured))
        if os.path.isabs(expanded):
            if os.path.exists(expanded) and os.access(expanded, os.X_OK):
                return expanded
            raise WorkflowError(
                f"read_filter.binary points at {expanded}, which is not an executable file."
            )
        resolved = shutil.which(configured)
        if resolved:
            return resolved
        raise WorkflowError(f"read_filter.binary='{configured}' is not on PATH.")

    for name in ("nanostream", "nanofilter"):
        resolved = shutil.which(name)
        if resolved:
            return resolved

    # Resolved on the submitting host, so a bundled binary has to be one the
    # compute nodes can also execute; keep it in the repo on shared storage.
    system = platform.system().lower()
    machine = platform.machine().lower()
    arch_alias = {"aarch64": "arm64", "arm64": "arm64", "x86_64": "x86_64", "amd64": "x86_64"}.get(
        machine, machine
    )
    bundled = os.path.join(workflow.basedir, "bin", f"{system}-{arch_alias}", "nanostream")
    if os.path.exists(bundled) and os.access(bundled, os.X_OK):
        return bundled

    raise WorkflowError(
        "read_filter is enabled but nanostream was not found. Put it on PATH, set "
        "read_filter.binary to its path, drop the binary in "
        f"workflow/bin/{system}-{arch_alias}/nanostream, or set read_filter.enabled: "
        "false to skip filtering."
    )


def is_read_filtering_enabled(sample):
    """Check if read pre-filtering is enabled for a given sample/assay."""
    if config.get("read_filter", {}).get("enabled") is False:
        return False
    setting = get_assay_setting(sample, "read_filter", "enabled", None)
    if setting is not None:
        return bool(setting)
    return bool(config.get("read_filter", {}).get("enabled", False))


def get_fastq_for_alignment(wildcards):
    """Return filtered FASTQ if filtering is enabled, else raw unit FASTQ."""
    if is_read_filtering_enabled(wildcards.sample):
        return f"reads/filtered/{wildcards.sample}_{wildcards.type}_{wildcards.processing_unit}_{wildcards.barcode}.fastq.gz".strip()
    return get_unit_input(wildcards, "fastq1")


def get_unmapped_bam_for_alignment(wildcards):
    """Return filtered unmapped BAM if filtering is enabled, else raw BAM."""
    if is_read_filtering_enabled(wildcards.sample):
        return f"reads/filtered/{wildcards.sample}_{wildcards.type}_{wildcards.processing_unit}_{wildcards.barcode}.unmapped.bam".strip()
    return get_unit_input(wildcards, "bam")



rule read_filter_fastq:
    """Filter raw FASTQ reads by QV and length before alignment."""
    input:
        fastq1=lambda wildcards: get_unit_input(wildcards, "fastq1"),
    output:
        fastq1=temp("reads/filtered/{sample}_{type}_{processing_unit}_{barcode}.fastq.gz"),
    params:
        binary=lambda wildcards: get_nanostream_binary(),
        min_qv=lambda wildcards: get_assay_setting(wildcards.sample, "read_filter", "min_qv", config.get("read_filter", {}).get("min_qv", 20)),
        length_range=lambda wildcards: get_assay_setting(wildcards.sample, "read_filter", "length_range", config.get("read_filter", {}).get("length_range", "500-3000")),
        extra=lambda wildcards: config.get("read_filter", {}).get("extra", ""),
    log:
        "logs/read_filter/{sample}_{type}_{processing_unit}_{barcode}.fastq.log",
    threads: config.get("read_filter", {}).get("threads", 4)
    resources:
        mem_mb=config.get("read_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("read_filter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("read_filter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("read_filter", {}).get("threads", 4),
        time=config.get("read_filter", {}).get("time", config["default_resources"]["time"]),
    # Deliberately no container: this rule runs the nanostream binary resolved
    # by get_nanostream_binary(), and a container would hide a binary found on
    # the host PATH. Set read_filter.container plus read_filter.binary together
    # if you want it containerised, with the binary at a path inside the image.
    container:
        config.get("read_filter", {}).get("container", None)
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        {params.binary} filter {input.fastq1} -o {output.fastq1} \
            --qv {params.min_qv} \
            --len {params.length_range} \
            -t {threads} {params.extra} > {log} 2>&1
        """


rule read_filter_unmapped_bam:
    """Filter raw unmapped BAM reads by QV and length before alignment."""
    input:
        bam=lambda wildcards: get_unit_input(wildcards, "bam"),
    output:
        bam=temp("reads/filtered/{sample}_{type}_{processing_unit}_{barcode}.unmapped.bam"),
    params:
        binary=lambda wildcards: get_nanostream_binary(),
        min_qv=lambda wildcards: get_assay_setting(wildcards.sample, "read_filter", "min_qv", config.get("read_filter", {}).get("min_qv", 20)),
        length_range=lambda wildcards: get_assay_setting(wildcards.sample, "read_filter", "length_range", config.get("read_filter", {}).get("length_range", "500-3000")),
        extra=lambda wildcards: config.get("read_filter", {}).get("extra", ""),
    log:
        "logs/read_filter/{sample}_{type}_{processing_unit}_{barcode}.bam.log",
    threads: config.get("read_filter", {}).get("threads", 4)
    resources:
        mem_mb=config.get("read_filter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("read_filter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("read_filter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("read_filter", {}).get("threads", 4),
        time=config.get("read_filter", {}).get("time", config["default_resources"]["time"]),
    # Deliberately no container: this rule runs the nanostream binary resolved
    # by get_nanostream_binary(), and a container would hide a binary found on
    # the host PATH. Set read_filter.container plus read_filter.binary together
    # if you want it containerised, with the binary at a path inside the image.
    container:
        config.get("read_filter", {}).get("container", None)
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        {params.binary} filter {input.bam} -o {output.bam} \
            --qv {params.min_qv} \
            --len {params.length_range} \
            -t {threads} {params.extra} > {log} 2>&1
        """


