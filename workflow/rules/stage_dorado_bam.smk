__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__license__ = "GPL-3"


rule stage_dorado_bam:
    """Add a GATK-compatible read group and index a pipeline-produced BAM."""
    input:
        demux_dir="basecalling/dorado_demux",
    output:
        bam="alignment/nagini_align/{sample}_{type}_{processing_unit}_{barcode}.dorado_bam.bam",
        bai="alignment/nagini_align/{sample}_{type}_{processing_unit}_{barcode}.dorado_bam.bam.bai",
    params:
        input_bam=lambda wildcards: get_unit_dorado_bam(wildcards),
        read_group=lambda wildcards: get_minimap2_read_group(wildcards),
        extra=config.get("stage_dorado_bam", {}).get("extra", ""),
    log:
        "logs/stage_dorado_bam/{sample}_{type}_{processing_unit}_{barcode}.log",
    benchmark:
        repeat(
            "benchmarks/stage_dorado_bam/{sample}_{type}_{processing_unit}_{barcode}.tsv",
            config.get("stage_dorado_bam", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("stage_dorado_bam", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("stage_dorado_bam", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("stage_dorado_bam", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("stage_dorado_bam", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("stage_dorado_bam", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("stage_dorado_bam", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("stage_dorado_bam", {}).get("container", config["default_container"])
    message:
        "{rule}: add sample metadata and sort/index Dorado output"
    shell:
        """
        samtools addreplacerg \
            -r '{params.read_group}' \
            -m overwrite_all \
            -o {output.bam}.with_rg.bam \
            {params.input_bam} &> {log}
        samtools sort -@ {threads} {params.extra} -o {output.bam} {output.bam}.with_rg.bam >> {log} 2>&1
        rm -f {output.bam}.with_rg.bam
        samtools index -@ {threads} {output.bam} {output.bai} >> {log} 2>&1
        """
