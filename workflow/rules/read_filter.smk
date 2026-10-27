__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


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
    benchmark:
        repeat(
            "benchmarks/read_filter_fastq/{sample}_{type}_{processing_unit}_{barcode}.tsv",
            config.get("read_filter", {}).get("benchmark_repeats", 1),
        )
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
    message:
        "{rule}: filter {input.fastq1} by quality and length"
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
    benchmark:
        repeat(
            "benchmarks/read_filter_unmapped_bam/{sample}_{type}_{processing_unit}_{barcode}.tsv",
            config.get("read_filter", {}).get("benchmark_repeats", 1),
        )
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
    message:
        "{rule}: filter {input.bam} by quality and length"
    shell:
        """
        {params.binary} filter {input.bam} -o {output.bam} \
            --qv {params.min_qv} \
            --len {params.length_range} \
            -t {threads} {params.extra} > {log} 2>&1
        """


