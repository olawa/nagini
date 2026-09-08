__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__license__ = "GPL-3"

"""Align pipeline-produced, demultiplexed Dorado BAMs per Hydra unit."""


_dorado_aligner = str(
    config.get("dorado_alignment", {}).get("aligner", "minimap2")
).strip().lower()
if _dorado_aligner not in ("minimap2", "dorado"):
    raise WorkflowError(
        "Unsupported dorado_alignment.aligner={!r}; use minimap2 or dorado.".format(
            _dorado_aligner
        )
    )

_dorado_alignment_raw = temp(
    "alignment/units_unsorted/"
    "{sample}_{type}_{processing_unit}_{barcode}.unsorted.bam"
)


if _dorado_aligner == "minimap2":

    rule align_dorado_bam_raw:
        """Align the exported per-barcode FASTQ with minimap2."""
        input:
            fastq=lambda wildcards: get_unit_dorado_fastq(wildcards),
            target=lambda wildcards: get_reference_target_for_unit(wildcards),
        output:
            bam=_dorado_alignment_raw,
        params:
            preset=lambda wildcards: get_minimap2_preset_for_unit(wildcards),
            extra=config.get("minimap2_align", {}).get("extra", ""),
        log:
            "logs/align_dorado_bam/{sample}_{type}_{processing_unit}_{barcode}.log",
        benchmark:
            repeat(
                "benchmarks/align_dorado_bam/"
                "{sample}_{type}_{processing_unit}_{barcode}.tsv",
                config.get("dorado_alignment", {}).get("benchmark_repeats", 1),
            )
        threads: config.get("dorado_alignment", {}).get("threads", config["default_resources"]["threads"])
        resources:
            mem_mb=config.get("dorado_alignment", {}).get(
                "mem_mb", config["default_resources"]["mem_mb"]
            ),
            mem_per_cpu=config.get("dorado_alignment", {}).get(
                "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
            ),
            partition=config.get("dorado_alignment", {}).get(
                "partition", config["default_resources"]["partition"]
            ),
            threads=config.get("dorado_alignment", {}).get(
                "threads", config["default_resources"]["threads"]
            ),
            time=config.get("dorado_alignment", {}).get(
                "time", config["default_resources"]["time"]
            ),
        container:
            config.get("dorado_alignment", {}).get(
                "container",
                config.get("minimap2_align", {}).get(
                    "container", config["default_container"]
                ),
            )
        conda:
            "../../envs/nagini.yaml"
        message:
            "{rule}: align Dorado barcode FASTQ {input.fastq} with minimap2"
        shell:
            r"""
            minimap2 -a -t {threads} -x {params.preset:q} \
                    {params.extra} {input.target:q} {input.fastq:q} 2> {log} \
                | samtools view -b -o {output.bam:q} - >> {log} 2>&1
            """


else:
    _dorado_alignment_img = config.get("dorado", {}).get("container")
    if _dorado_alignment_img in ("", None, "none", "None", False, "local"):
        _dorado_alignment_img = None
    _dorado_alignment_exe = config.get("dorado", {}).get("binary", "")
    _dorado_alignment_exe = (
        os.path.expanduser(str(_dorado_alignment_exe))
        if _dorado_alignment_exe
        else "dorado"
    )

    rule align_dorado_bam_raw:
        """Align a demultiplexed BAM directly with Dorado aligner."""
        input:
            bam=lambda wildcards: get_unit_dorado_bam(wildcards),
            target=lambda wildcards: get_reference_target_for_unit(wildcards),
        output:
            bam=_dorado_alignment_raw,
        params:
            dorado_bin=_dorado_alignment_exe,
            extra=config.get("dorado_alignment", {}).get("extra", ""),
        log:
            "logs/align_dorado_bam/{sample}_{type}_{processing_unit}_{barcode}.log",
        benchmark:
            repeat(
                "benchmarks/align_dorado_bam/"
                "{sample}_{type}_{processing_unit}_{barcode}.tsv",
                config.get("dorado_alignment", {}).get("benchmark_repeats", 1),
            )
        threads: config.get("dorado_alignment", {}).get("threads", config["default_resources"]["threads"])
        resources:
            mem_mb=config.get("dorado_alignment", {}).get(
                "mem_mb", config["default_resources"]["mem_mb"]
            ),
            mem_per_cpu=config.get("dorado_alignment", {}).get(
                "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
            ),
            partition=config.get("dorado_alignment", {}).get(
                "partition", config["default_resources"]["partition"]
            ),
            threads=config.get("dorado_alignment", {}).get(
                "threads", config["default_resources"]["threads"]
            ),
            time=config.get("dorado_alignment", {}).get(
                "time", config["default_resources"]["time"]
            ),
        container:
            _dorado_alignment_img
        message:
            "{rule}: align unaligned Dorado BAM {input.bam} with Dorado"
        shell:
            r"""
            {params.dorado_bin} aligner {input.target:q} {input.bam:q} \
                --threads {threads} {params.extra} > {output.bam:q} 2> {log}
            """


rule finalize_dorado_bam_alignment:
    """Add sample metadata, coordinate-sort, and index an aligned Dorado unit."""
    input:
        bam=(
            "alignment/units_unsorted/"
            "{sample}_{type}_{processing_unit}_{barcode}.unsorted.bam"
        ),
    output:
        bam=(
            "alignment/units/"
            "{sample}_{type}_{processing_unit}_{barcode}.from_dorado.bam"
        ),
        bai=(
            "alignment/units/"
            "{sample}_{type}_{processing_unit}_{barcode}.from_dorado.bam.bai"
        ),
    params:
        read_group=lambda wildcards: get_minimap2_read_group(wildcards),
        extra=config.get("finalize_dorado_alignment", {}).get("extra", ""),
    log:
        "logs/finalize_dorado_alignment/"
        "{sample}_{type}_{processing_unit}_{barcode}.log",
    benchmark:
        repeat(
            "benchmarks/finalize_dorado_alignment/"
            "{sample}_{type}_{processing_unit}_{barcode}.tsv",
            config.get("finalize_dorado_alignment", {}).get(
                "benchmark_repeats", 1
            ),
        )
    threads: config.get("finalize_dorado_alignment", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("finalize_dorado_alignment", {}).get(
            "mem_mb", config["default_resources"]["mem_mb"]
        ),
        mem_per_cpu=config.get("finalize_dorado_alignment", {}).get(
            "mem_per_cpu", config["default_resources"]["mem_per_cpu"]
        ),
        partition=config.get("finalize_dorado_alignment", {}).get(
            "partition", config["default_resources"]["partition"]
        ),
        threads=config.get("finalize_dorado_alignment", {}).get(
            "threads", config["default_resources"]["threads"]
        ),
        time=config.get("finalize_dorado_alignment", {}).get(
            "time", config["default_resources"]["time"]
        ),
    container:
        config.get("finalize_dorado_alignment", {}).get(
            "container", config["default_container"]
        )
    message:
        "{rule}: add RG, sort and index {input.bam}"
    shell:
        r"""
        samtools addreplacerg \
            -r '{params.read_group}' \
            -m overwrite_all -u -o - {input.bam:q} 2> {log} \
            | samtools sort -@ {threads} {params.extra} \
                -o {output.bam:q} - >> {log} 2>&1
        samtools index -@ {threads} {output.bam:q} {output.bai:q} \
            >> {log} 2>&1
        """
