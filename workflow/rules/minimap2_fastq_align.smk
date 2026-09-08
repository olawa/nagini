__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule minimap2_fastq_align:
    """Align a FASTQ unit with the Hydra minimap2 index."""
    input:
        fastq1=lambda wildcards: get_fastq_for_alignment(wildcards),
        target=lambda wildcards: get_reference_target_for_unit(wildcards),

    output:
        bam="alignment/units/{sample}_{type}_{processing_unit}_{barcode}.from_fastq.bam",
        bai="alignment/units/{sample}_{type}_{processing_unit}_{barcode}.from_fastq.bam.bai",
    params:
        preset=lambda wildcards: get_minimap2_preset_for_unit(wildcards),
        extra=config.get("minimap2_align", {}).get("extra", ""),
        read_group=lambda wildcards: get_minimap2_read_group(wildcards),
    log:
        "logs/minimap2_fastq_align/{sample}_{type}_{processing_unit}_{barcode}.log",
    benchmark:
        repeat(
            "benchmarks/minimap2_fastq_align/{sample}_{type}_{processing_unit}_{barcode}.tsv",
            config.get("minimap2_fastq_align", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("minimap2_fastq_align", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("minimap2_fastq_align", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("minimap2_fastq_align", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("minimap2_fastq_align", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("minimap2_fastq_align", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("minimap2_fastq_align", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("minimap2_fastq_align", {}).get("container", config.get("minimap2_align", {}).get("container", config["default_container"]))
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: minimap2 {input.fastq1} against {input.target}"
    shell:
        """
        minimap2 -a -t {threads} -x {params.preset} \
            -R '{params.read_group}' {params.extra} \
            {input.target} {input.fastq1} 2> {log} | \
        samtools sort -@ {threads} -o {output.bam} - >> {log} 2>&1
        samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1
        """
