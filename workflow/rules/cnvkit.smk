__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"

import os


def _get_ref(key, config):
    """Return a reference file path resolved via reference_build."""
    build = config.get("reference_build", "hg38")
    return config["references"][build][key]


rule cnvkit_targets:
    """Tile the accessible genome into fixed bins, once per reference build.

    `cnvkit.py batch -m wgs` derives its own bins when given no --targets, and
    on T2T that produced an empty target set: batch wrote "target.bed with 0
    regions" and then died in do_reference_flat, where zip(*[]) cannot unpack
    an empty probe list. Running the same two steps by hand on the same FASTA
    and the same bin size yields 24 accessible regions and 62344 bins, so it
    is batch's internal tiling that fails, not access or target.

    The likely reason is structural rather than T2T-specific: hg38's N blocks
    split each chromosome into many accessible regions, while a gapless
    assembly gives one interval per chromosome, up to 248 Mb wide. Generating
    the bins here and passing them in sidesteps that path entirely, and costs
    one run per reference build rather than one per sample.
    """
    input:
        fasta=lambda wildcards: _get_ref("fasta", config),
    output:
        access=temp("references/cnvkit_access_{build}.bed"),
        targets="references/cnvkit_targets_{build}.bed",
    params:
        avg_size=config.get("cnvkit_batch", {}).get("target_avg_size", 50000),
    log:
        "logs/cnvkit_targets/cnvkit_targets_{build}.log",
    threads: 1
    resources:
        mem_mb=config.get("cnvkit_batch", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_batch", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_batch", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("cnvkit_batch", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_batch", {}).get("container", config["default_container"])
    wildcard_constraints:
        build=r"[A-Za-z0-9._-]+",
    message:
        "{rule}: tile the accessible genome into {params.avg_size} bp bins for CNVkit"
    shell:
        """
        exec > {log} 2>&1
        cnvkit.py access {input.fasta} -o {output.access}
        cnvkit.py target {output.access} --split --avg-size {params.avg_size} -o {output.targets}

        if [ ! -s "{output.targets}" ]; then
            echo "Error: cnvkit.py target produced no bins from {output.access}."
            exit 1
        fi
        """


rule cnvkit_batch:
    """Run CNVkit in WGS mode on a tumor-only ONT BAM.

    Bins come from cnvkit_targets rather than from batch's own tiling, which
    produced an empty target set on a gapless assembly; see that rule.

    No Panel of Normals is needed: -n builds an internal flat reference.
    --no-gc skips GC-bias correction (ONT is PCR-free, so correction
    introduces more noise than it removes).
    --count-reads counts read midpoints instead of computing per-base depth,
    which is orders of magnitude faster for long-read data.
    """
    input:
        bam="alignment/merged/{sample}_{type}.bam",
        bai="alignment/merged/{sample}_{type}.bam.bai",
        fasta=lambda wildcards: _get_ref("fasta", config),
        targets=lambda wildcards: "references/cnvkit_targets_{}.bed".format(
            config.get("reference_build", "hg38")
        ),
    output:
        cnr="cnv/cnvkit/{sample}_{type}.cnr",
        cns="cnv/cnvkit/{sample}_{type}.cns",
    params:
        tmpdir="cnv/cnvkit/tmp_{sample}_{type}",
        method=config.get("cnvkit_batch", {}).get("method", "wgs"),
        # No --target-avg-size: the bins arrive through --targets, and passing
        # both would ask batch to re-derive what it was just handed.
        extra=config.get("cnvkit_batch", {}).get("extra", "--count-reads"),
    log:
        "cnv/cnvkit/{sample}_{type}.cnvkit_batch.log",
    benchmark:
        repeat(
            "cnv/cnvkit/{sample}_{type}.cnvkit_batch.benchmark.tsv",
            config.get("cnvkit_batch", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cnvkit_batch", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cnvkit_batch", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cnvkit_batch", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cnvkit_batch", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cnvkit_batch", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cnvkit_batch", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cnvkit_batch", {}).get("container", config["default_container"])
    message:
        "{rule}: CNV calling with CNVkit (WGS, tumor-only) for {wildcards.sample}_{wildcards.type}"
    shell:
        "mkdir -p {params.tmpdir} && "
        "cnvkit.py batch {input.bam} -n -f {input.fasta} "
        "--output-reference {params.tmpdir}/reference.cnn "
        "--targets {input.targets} "
        "-m {params.method} -d {params.tmpdir} -p {threads} {params.extra} "
        "> {log} 2>&1 && "
        "bam_base=$(basename {input.bam} .bam) && "
        "mv {params.tmpdir}/${{bam_base}}.cnr {output.cnr} && "
        "mv {params.tmpdir}/${{bam_base}}.cns {output.cns} && "
        "rm -rf {params.tmpdir}"
