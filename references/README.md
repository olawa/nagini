# Reference bundles

The pipeline does not download reference data.  Put prepared, indexed bundles
under `references/hg38/` and/or `references/hs1/`, or replace the paths in the
selected config file.

Each selected build needs at least:

* FASTA and `.fai`
* a sequence dictionary if GATK/Mutect2 is enabled
* a BED file for target-region coverage and Mutect2 interval splitting
* a build-compatible offline VEP cache if
  `annotate_small_variants_with_vep` is enabled

`reference_build` selects either `hg38` or `hs1`.  The workspace contained
known hg38/GRCh38 paths but no hs1/T2T path, so the hs1 entries in the example
config are deliberately placeholders.
