__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"

import os
import platform
import shutil


def get_nanofilter_binary():
    """Resolve the nanofilter executable path across platforms."""
    configured = config.get("read_filter", {}).get("binary", "")
    if configured:
        expanded = os.path.expanduser(os.path.expandvars(configured))
        if os.path.isabs(expanded) and os.path.exists(expanded):
            return expanded
        if shutil.which(configured):
            return shutil.which(configured)

    path_bin = shutil.which("nanofilter")
    if path_bin:
        return path_bin

    system = platform.system().lower()
    machine = platform.machine().lower()
    if machine in ("aarch64", "arm64"):
        arch_alias = "arm64"
    elif machine in ("x86_64", "amd64"):
        arch_alias = "x86_64"
    else:
        arch_alias = machine

    bundled_dir = os.path.join(workflow.basedir, "bin", f"{system}-{arch_alias}")
    bundled_exe = os.path.join(bundled_dir, "nanofilter")
    if os.path.exists(bundled_exe) and os.access(bundled_exe, os.X_OK):
        return bundled_exe

    fallback_script = os.path.join(workflow.basedir, "scripts", "nanofilter_fallback.py")
    return f"python3 {fallback_script}"


def is_read_filtering_enabled(sample):
    """Check if read pre-filtering is enabled for a given sample/assay."""
    setting = get_assay_setting(sample, "read_filter", "enabled", None)
    if setting is not None:
        return bool(setting)
    return bool(config.get("read_filter", {}).get("enabled", False))


def get_fastq_for_alignment(wildcards):
    """Return filtered FASTQ if filtering is enabled, else raw unit FASTQ."""
    if is_read_filtering_enabled(wildcards.sample):
        return f"reads/filtered/{wildcards.sample}_{wildcards.type}_{wildcards.processing_unit}_{wildcards.barcode}.fastq.gz"
    return get_unit_input(wildcards, "fastq1")


def get_unmapped_bam_for_alignment(wildcards):
    """Return filtered unmapped BAM if filtering is enabled, else raw BAM."""
    if is_read_filtering_enabled(wildcards.sample):
        return f"reads/filtered/{wildcards.sample}_{wildcards.type}_{wildcards.processing_unit}_{wildcards.barcode}.unmapped.bam"
    return get_unit_input(wildcards, "bam")


rule read_filter_fastq:
    """Filter raw FASTQ reads by QV and length before alignment."""
    input:
        fastq1=lambda wildcards: get_unit_input(wildcards, "fastq1"),
    output:
        fastq1=temp("reads/filtered/{sample}_{type}_{processing_unit}_{barcode}.fastq.gz"),
    params:
        binary=lambda wildcards: get_nanofilter_binary(),
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
        binary=lambda wildcards: get_nanofilter_binary(),
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
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        {params.binary} filter {input.bam} -o {output.bam} \
            --qv {params.min_qv} \
            --len {params.length_range} \
            -t {threads} {params.extra} > {log} 2>&1
        """
