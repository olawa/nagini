# Reference bundles

The pipeline does not download reference data.  Put prepared, indexed bundles
under `references/hg38/` and/or `references/hs1/`, or select a site config that
points to a shared reference directory.

Each selected build needs at least:

* FASTA and `.fai`
* a sequence dictionary if GATK/Mutect2 is enabled
* a BED file for target-region coverage and Mutect2 interval splitting
* a build-compatible offline VEP cache if
  `annotate_small_variants_with_vep` is enabled

The workstation profile expects the following GRCh38 layout under
`/data/references`:

* `GRCh38_no_alt_analysis_set.fasta`, `.fai`, `.dict`, and `.mmi`
* `ReadUntil_Final_Extra_30kb_hg38_v5.bed`
* `human_GRCh38_no_alt_analysis_set.trf.bed` for Sniffles2 tandem repeats
* `clinvar.vcf.gz` and `clinvar.vcf.gz.tbi` for optional VEP ClinVar fields
* `gencode.v46.annotation.sorted.gtf` and `annotations` (reference provenance)
* `vep_cache/` when offline VEP is enabled

The corresponding site config is
`config/site_configs/site_config_workstation.yaml`; its registry is
`config/reference_files/reference_files_workstation.yaml`. The registry is a
manifest, not a downloader. Fill in its checksum fields after the workstation
bundle has been checked.

ClinVar is used only when VEP is enabled and both the VCF and its tabix index
are readable. Cytoband annotation is optional; the workstation config leaves
`cytoband` empty because no `cytoBand.txt`/`cytoBand.txt.gz` was present in the
installed file list. Add that file and set `HG38_CYTOBAND` when it is available.

`reference_build` selects either `hg38` or `hs1`. The checkout-relative local
config keeps hs1 as a placeholder; the workstation config uses the installed
`hs1.fa`, `.fai`, `.dict`, and `.mmi` files.

Select a runtime site with `pipeline_start_scripts/start_quicken.sh --site`
(`local`, `workstation`, or `marvin`).
The launcher combines the shared config with the matching file under
`config/site_configs/` and keeps the corresponding manifest under
`config/reference_files/` in the analysis directory. Fill in the `checksum`
fields in the selected registry when a workstation reference bundle has been
verified.
