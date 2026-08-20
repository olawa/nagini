__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__license__ = "GPL-3"

"""Containerised Dorado preprocessing for MinKNOW/POD5 input."""


dorado_barcodes = sorted(
    {
        str(value).strip()
        for value in units.get("barcode_name", pd.Series(dtype=str)).tolist()
        if str(value).strip()
    }
)


_dorado_img = config.get("dorado", {}).get("container")
if _dorado_img in ("", None, "none", "None", False, "local"):
    _dorado_img = None

_dorado_exe = config.get("dorado", {}).get("binary", "")
_dorado_exe = os.path.expanduser(str(_dorado_exe)) if _dorado_exe else "dorado"


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
    """Demultiplex basecalled reads and normalize one BAM per barcode."""
    input:
        bam="basecalling/dorado/basecalled.bam",
    output:
        demux_dir=directory("basecalling/dorado_demux"),
    params:
        dorado_bin=_dorado_exe,
        kit=config.get("dorado", {}).get("kit", ""),
        threads=config.get("dorado", {}).get("demux_threads", config["default_resources"]["threads"]),
        barcodes=" ".join(dorado_barcodes),
        extra=config.get("dorado", {}).get("demux_extra", ""),
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
        raw="{output.demux_dir}.raw"
        rm -rf "$raw" "{output.demux_dir}"
        mkdir -p "$raw" "{output.demux_dir}"
        {params.dorado_bin} demux --output-dir "$raw" --kit-name {params.kit} \
            --no-trim --emit-summary {params.extra} {input.bam:q} > {log} 2>&1
        for barcode in {params.barcodes}; do
            # A read loop instead of mapfile keeps this working on the bash 3.x
            # that some container images still ship.
            bam_files=()
            while IFS= read -r bam_file; do
                bam_files+=("$bam_file")
            done < <(find "$raw" -type f \( -path "*/${{barcode}}/*.bam" -o -name "${{barcode}}.bam" -o -name "*${{barcode}}*.bam" \) -print)
            if [[ "${{#bam_files[@]}}" -eq 0 ]]; then
                # A barcode listed in the sample sheet can legitimately receive
                # no reads (failed library, empty well).  Emitting a header-only
                # BAM keeps that sample's own branch failing on its own merits
                # instead of discarding the basecalling for every other sample.
                echo "WARNING: no demultiplexed reads for $barcode; writing an empty BAM" >&2
                echo "WARNING: no demultiplexed reads for $barcode" >> {log}
                samtools view -H {input.bam:q} \
                    | samtools view -b -o "{output.demux_dir}/${{barcode}}.bam" - >> {log} 2>&1
            else
                samtools merge -f -@ {threads} "{output.demux_dir}/${{barcode}}.bam" "${{bam_files[@]}}" >> {log} 2>&1
            fi
            samtools index -@ {threads} "{output.demux_dir}/${{barcode}}.bam" >> {log} 2>&1
        done
        rm -rf "$raw"
        """
