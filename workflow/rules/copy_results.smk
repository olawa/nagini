__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"

"""Copy the deliverables named in the output spec into results/.

An entry in config/output_files.yaml with an `input:` path is produced in
place by its own rule and copied here, so results/ holds the run's
deliverables while the working directories keep the intermediates. Entries
with `input: null` are delivered in place only; read-level data (the merged
and target BAMs, the merged FASTQ, per-base coverage) is deliberately among
them, since copying it would roughly double a run's footprint.

This used to generate one rule per entry by exec'ing decorator source. That
broke on Snakemake 9, which cannot resolve a rule whose source file it has
never seen, and it had carried two latent bugs the whole time: `when:` was
read as a plain config flag, so no conditional entry ever produced a copy
rule, and mem_mb was emitted as a string. One generic rule needs none of it.
"""


rule copy_to_results:
    """Copy one deliverable into results/, leaving the in-place file alone."""
    input:
        # results/<path> comes from <path>; the constraint stops results/
        # from being its own source.
        lambda wildcards: wildcards.path,
    output:
        "results/{path}",
    wildcard_constraints:
        path=r"(?!results/).+",
    params:
        extra=config.get("_copy", {}).get("extra", ""),
    log:
        "logs/copy_to_results/{path}.log",
    benchmark:
        repeat(
            "benchmarks/copy_to_results/{path}.tsv",
            config.get("_copy", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("_copy", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("_copy", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("_copy", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("_copy", {}).get("container", config["default_container"])
    message:
        "{rule}: deliver {input} to {output}"
    shell:
        "(cp --preserve=timestamps {params.extra} {input} {output}) &> {log}"
