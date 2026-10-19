__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__license__ = "GPL-3"

"""Containerised Dorado basecalling and barcode-level raw-data export."""


_dorado_img = config.get("dorado", {}).get("container")
if _dorado_img in ("", None, "none", "None", False, "local"):
    _dorado_img = None

_dorado_exe = config.get("dorado", {}).get("binary", "")
_dorado_exe = os.path.expanduser(str(_dorado_exe)) if _dorado_exe else "dorado"

_dorado_output_format = str(
    config.get("dorado", {}).get("demux_output_format", "fastq")
).lower()
_dorado_keep_bam = _dorado_output_format in ("bam", "both")
_dorado_keep_fastq = _dorado_output_format in ("fastq", "both")
_dorado_barcode_bam = "basecalling/dorado_demux/{dorado_barcode}.bam"
_dorado_barcode_bam_output = (
    _dorado_barcode_bam if _dorado_keep_bam else temp(_dorado_barcode_bam)
)
_dorado_barcode_fastq = "basecalling/dorado_fastq/{dorado_barcode}.fastq.gz"
_dorado_barcode_fastq_output = (
    _dorado_barcode_fastq
    if _dorado_keep_fastq
    else temp(_dorado_barcode_fastq)
)


rule dorado_basecaller:
    """Basecall a MinKNOW/POD5 directory and retain barcode sequence context."""
    input:
        reads=lambda wildcards: config["dorado"]["input"],
    output:
        bam="basecalling/dorado/basecalled.bam",
    params:
        dorado_bin=_dorado_exe,
        model=config.get("dorado", {}).get("model", "hac"),
        kit=config.get("dorado", {}).get("kit", ""),
        device_arg=(
            f"--device {config['dorado']['device']}"
            if config.get("dorado", {}).get("device")
            else ""
        ),
        models_directory_arg=(
            f"--models-directory {config['dorado']['models_directory']}"
            if config.get("dorado", {}).get("models_directory")
            else ""
        ),
        extra=config.get("dorado", {}).get("basecaller_extra", ""),
    log:
        "logs/dorado/basecalling.log",
    benchmark:
        repeat(
            "benchmarks/dorado/basecalling.tsv",
            config.get("dorado", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("dorado", {}).get("basecaller_threads", config["default_resources"]["threads"])
    resources:
        # A GPU partition alone does not allocate a device.  Sites that submit
        # through a scheduler must also request the device itself; the resource
        # value carries the complete flag so that it can be empty elsewhere.
        gres=config.get("dorado", {}).get(
            "basecaller_gres", config["default_resources"].get("gres", "")
        ),
        mem_mb=config.get("dorado", {}).get("basecaller_mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("dorado", {}).get("basecaller_mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("dorado", {}).get("basecaller_partition", config["default_resources"]["partition"]),
        threads=config.get("dorado", {}).get("basecaller_threads", config["default_resources"]["threads"]),
        time=config.get("dorado", {}).get("basecaller_time", config["default_resources"]["time"]),
    container:
        _dorado_img
    message:
        "{rule}: basecall MinKNOW/POD5 input with Dorado"
    shell:
        r"""
        {params.dorado_bin} basecaller {params.model:q} {input.reads:q} --recursive \
            --no-trim \
            {params.device_arg} {params.models_directory_arg} \
            {params.extra} > {output.bam} 2> {log}
        """


rule dorado_demux:
    """Classify reads and retain Dorado's native demultiplexing directory."""
    input:
        bam="basecalling/dorado/basecalled.bam",
    output:
        raw_dir=temp(directory("basecalling/dorado_demux_raw")),
    params:
        dorado_bin=_dorado_exe,
        kit=config.get("dorado", {}).get("kit", ""),
        extra=lambda wildcards: (
            ("--barcode-both-ends " if (config.get("dorado", {}).get("demux_both_ends", False) or config.get("dorado", {}).get("both_ends", False) or config.get("dorado", {}).get("barcode_both_ends", False)) else "")
            + config.get("dorado", {}).get("demux_extra", "")
        ).strip(),


    log:
        "logs/dorado/demux.log",
    benchmark:
        repeat(
            "benchmarks/dorado/demux.tsv",
            config.get("dorado", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("dorado", {}).get("demux_threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("dorado", {}).get("demux_mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("dorado", {}).get("demux_mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("dorado", {}).get("demux_partition", config["default_resources"]["partition"]),
        threads=config.get("dorado", {}).get("demux_threads", config["default_resources"]["threads"]),
        time=config.get("dorado", {}).get("demux_time", config["default_resources"]["time"]),
    container:
        _dorado_img
    message:
        "{rule}: demultiplex Dorado BAM output by barcode"
    shell:
        r"""
        rm -rf {output.raw_dir:q}
        mkdir -p {output.raw_dir:q}
        {params.dorado_bin} demux --output-dir {output.raw_dir:q} \
            --kit-name {params.kit:q} --threads {threads} \
            --no-trim --emit-summary {params.extra} {input.bam:q} > {log} 2>&1
        """


rule dorado_demux_summary:
    """Copy Dorado's demultiplexing summary out of its raw output directory.

    A separate step because the Dorado image has no Python and the other
    containers have no Dorado."""
    input:
        raw_dir="basecalling/dorado_demux_raw",
    output:
        summary="basecalling/dorado_demux/demux_summary.tsv",
    log:
        "logs/dorado/demux_summary.log",
    benchmark:
        repeat(
            "benchmarks/dorado_demux_summary/dorado_demux_summary.tsv",
            config.get("dorado_demux_summary", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("dorado_demux_summary", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("dorado_demux_summary", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("dorado_demux_summary", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("dorado_demux_summary", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("dorado_demux_summary", {}).get("container", config["default_container"])
    message:
        "{rule}: take the demultiplexing summary from {input.raw_dir}"
    shell:
        """
        python3 {workflow.basedir}/scripts/dorado_demux_collect.py summary \
            --raw-dir {input.raw_dir:q} \
            --output {output.summary:q} > {log} 2>&1
        """


rule dorado_demux_barcode_bam:
    """Normalize Dorado output to one explicit, unaligned BAM per barcode."""
    input:
        raw_dir="basecalling/dorado_demux_raw",
        basecalled="basecalling/dorado/basecalled.bam",
    output:
        bam=_dorado_barcode_bam_output,
    params:
        barcode=lambda wildcards: wildcards.dorado_barcode,
        allow_empty_flag="--allow-empty" if config.get("dorado", {}).get("allow_empty_barcodes", False) else "",
    log:
        "logs/dorado/demux_{dorado_barcode}.log",
    benchmark:
        repeat(
            "benchmarks/dorado/demux_{dorado_barcode}.tsv",
            config.get("dorado", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("dorado_demux_bam", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("dorado_demux_bam", {}).get(
            "mem_mb", config["default_resources"]["mem_mb"]
        ),
        mem_per_cpu=config.get("dorado_demux_bam", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("dorado_demux_bam", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("dorado_demux_bam", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("dorado_demux_bam", {}).get(
            "time", config["default_resources"]["time"]
        ),
    container:
        config.get("dorado_demux_bam", {}).get(
            "container", config["default_container"]
        )
    wildcard_constraints:
        dorado_barcode="barcode[0-9]+",
    message:
        "{rule}: retain unaligned reads assigned to {wildcards.dorado_barcode}"
    shell:
        r"""
        python3 {workflow.basedir}/scripts/dorado_demux_collect.py barcode \
            --raw-dir {input.raw_dir:q} \
            --barcode {params.barcode:q} \
            --basecalled {input.basecalled:q} \
            --output {output.bam:q} \
            {params.allow_empty_flag} > {log} 2>&1
        """


rule dorado_demux_fastq:
    """Optionally export a gzip-compressed FASTQ for a demultiplexed barcode."""
    input:
        bam="basecalling/dorado_demux/{dorado_barcode}.bam",
    output:
        fastq=_dorado_barcode_fastq_output,
    params:
        extra=config.get("samtools_fastq", {}).get("extra", "-n -T MM,ML"),
    log:
        "logs/dorado/fastq_{dorado_barcode}.log",
    benchmark:
        repeat(
            "benchmarks/dorado/fastq_{dorado_barcode}.tsv",
            config.get("dorado_fastq", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("dorado_fastq", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("dorado_fastq", {}).get(
            "mem_mb", config["default_resources"]["mem_mb"]
        ),
        mem_per_cpu=config.get("dorado_fastq", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("dorado_fastq", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("dorado_fastq", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("dorado_fastq", {}).get(
            "time", config["default_resources"]["time"]
        ),
    container:
        config.get("dorado_fastq", {}).get(
            "container", config["default_container"]
        )
    wildcard_constraints:
        dorado_barcode="barcode[0-9]+",
    message:
        "{rule}: export {wildcards.dorado_barcode} as FASTQ"
    shell:
        r"""
        samtools fastq -@ {threads} {params.extra} {input.bam:q} 2> {log} \
            | gzip -c > {output.fastq:q}
        """


rule dorado_demux_qc:
    """Aggregate the per-read demux summary into a barcode yield table.

    Run level rather than per sample, because the question it answers -- which
    barcode got no reads, and how much went unclassified -- is about the
    flowcell.
    """
    input:
        summary="basecalling/dorado_demux/demux_summary.tsv",
        units=config["units"],
    output:
        mqc="qc/dorado_demux/dorado_demux_mqc.tsv",
    params:
        kit=config.get("dorado", {}).get("kit", ""),
    log:
        "logs/dorado/demux_qc.log",
    benchmark:
        repeat(
            "benchmarks/dorado/demux_qc.tsv",
            config.get("dorado", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("dorado_demux_qc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("dorado_demux_qc", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("dorado_demux_qc", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=1,
        time=config.get("dorado_demux_qc", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("dorado_demux_qc", {}).get("container", config["default_container"])
    message:
        "{rule}: summarise reads per barcode for the run report"
    script:
        "../scripts/dorado_demux_summary.py"
