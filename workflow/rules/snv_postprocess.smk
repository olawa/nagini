__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"

import os


def get_vep_custom_annotations(wildcards):
    """Add ClinVar as an optional, indexed VEP custom annotation."""

    clinvar = get_clinvar_vcf(wildcards)
    if not clinvar or not os.path.isfile(f"{clinvar}.tbi"):
        return ""
    return (
        f"--custom {clinvar},short_name=ClinVar,format=vcf,type=exact,"
        "fields=CLNSIG%CLNREVSTAT%CLNDN"
    )


rule create_genome_contigs_bed:
    """Create a whole-genome BED from FASTA index when no design BED is provided."""
    input:
        fai=config["reference"]["fai"],
    output:
        bed="references/genome_contigs.bed",
    threads: 1
    container:
        config["default_container"]
    shell:
        """
        awk 'BEGIN {{OFS="\\t"}} {{print $1, 0, $2}}' {input.fai} > {output.bed}
        """
