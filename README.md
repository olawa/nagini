# Nagini

Nagini is an Oxford Nanopore sequencing analysis pipeline. It turns raw or
pre-aligned reads into normalized alignments, quality-control metrics, small
variant calls, structural variant calls, and a consolidated analysis report.

## Input modes

Nagini accepts four kinds of sample input:

| Input | Processing |
| --- | --- |
| MinKNOW/POD5 folder | Dorado basecalls and demultiplexes reads per barcode; the selected aligner then maps each barcode independently. |
| Existing barcode FASTQ folder | FASTQs are collected recursively across runs and `fastq_pass`/`fastq_fail`, merged per barcode, and aligned directly. |
| FASTQ | Reads are aligned to the selected reference genome. |
| Unmapped BAM | Reads are extracted, aligned, and normalized. |
| Aligned BAM | Existing alignments are retained; the BAM is prepared, sorted, and indexed for downstream variant calling. |

### MinKNOW / POD5 Input & Dorado Demux

For MinKNOW/POD5 runs, the launcher accepts a Dorado model, modified basecalling flags (e.g. `5mCG_5hmCG`), and a barcode kit. It schedules basecalling and demultiplexing as part of the workflow.

If a folder contains both existing FASTQ files and POD5 files, the launcher uses
the FASTQs by default. It only schedules new basecalling when `--input-type
minknow` is explicitly given. Barcode FASTQs can be mapped with `--barcode-map`;
`sample_sheet*.csv/tsv` files found below the input folder are used automatically
when they contain barcode and sample columns. By default both `fastq_pass` and
`fastq_fail` are included; use `--fastq-status pass` for pass reads only.

A barcode map CSV/TSV maps barcode identifiers (e.g. `barcode01`) to clinical sample names and types (`T` for tumor, `N` for normal). Dorado can be run via container (default) or a local binary via `--dorado-binary`.

Example starting from raw MinKNOW POD5 data:

```bash
pipeline_start_scripts/start_nagini.sh \
  --site workstation \
  --input-type minknow \
  --inbox-path /data/minknow_run \
  --barcode-kit SQK-NBD114-24 \
  --basecalling-model "hac,5mCG_5hmCG" \
  --methylation Yes \
  --barcode-map barcode_map.csv \
  --dorado-binary ~/tools/dorado-2.1.1-linux-x64/bin/dorado \
  --analysis-dir /data/analysis/run_01
```

Example starting from existing aligned BAM files:

```bash
pipeline_start_scripts/start_nagini.sh \
  --site workstation \
  --input-type aligned_bam \
  --inbox-path /data/analysis/run_01/alignment/merged \
  --analysis-dir /data/analysis/run_01
```

Example starting from several completed barcode FASTQ runs for one assay:

```bash
pipeline_start_scripts/start_nagini.sh \
  --site workstation \
  --inbox-path /data/bcr_abl_dilution_P2_run_folders \
  --analysis-dir /data/analysis/bcr_abl_dilution_P2 \
  --assay bcr_abl \
  --library-type cdna \
  --reference-build hg19 \
  --barcode-map /data/bcr_abl/barcode_map.csv \
  --basecalling-model sup \
  --cores 32
```

The launcher writes one merged `input_fastq/barcodeXX.fastq.gz` and one
`units.tsv` row per observed barcode. Omit `--barcode-map` when the exported
sample sheets provide the mapping. Use `--input-type minknow` only when the
POD5 data should be basecalled again.

## Analysis & Variant Calling

For each sample, Nagini performs:

- **Alignment & BAM staging:** Minimap2 alignment with `@RG` sample tags and samtools sorting/indexing.
- **Coverage QC:** Mosdepth whole-genome and target-region (BED) coverage analysis.
- **Read & Yield QC:** Cramino read-length, N50, and Q-score statistics.
- **Small Variant Calling:** chosen with `--small-variant-caller auto|clairs_to|deepvariant|deepsomatic|mutect2|both`.
  - **auto** (default): ClairS-TO for tumour (`T`) and cDNA samples, DeepVariant for normal (`N`) samples.
  - **ClairS-TO**: Tumour-only somatic SNV and small-indel caller for ONT R10.4.1, with the platform model derived from the basecalling model. Tags non-somatic variants against its built-in panels of normals (gnomAD, dbSNP, 1000 Genomes, CoLoRSdb).
  - **DeepVariant**: Germline calling with the `ONT_R104` model.
  - **DeepSomatic**: Tumour-only calling with the `ONT_TUMOR_ONLY` model; its built-in panel of normals is applied on hg38. Only runs when selected by name.
  - **GATK Mutect2**: Retained on request; experimental for ONT. `both` runs ClairS-TO and Mutect2.
- **Structural Variant Calling**, genome-wide, with one row per caller in the report:
  - **Sniffles2**, plus a second run in `--mosaic` mode for low-allele-fraction (subclonal) SVs, reported as `Sniffles2-mosaic`.
  - **SAVANA** in tumour-only mode.
  - **Rindels** `--sv-lite` (DEL/DUP/INV/TRA from pair and clip evidence).
  - Translocation/BND extraction, cytogenetic band notation (ISCN-style) and target BED overlap for each. `sv_engine` restricts the set.
- **Copy Number Calling**, genome-wide from binned read counts (50 kb): **Rindels** (ASCAT-style segmentation with BAF, purity and ploidy) and **CNVkit** (`wgs`, read counts, flat reference), with an ideogram and an interactive HTML report.
- **Annotation:** Optional VEP annotation with ClinVar pathogenicity fields and SnpEff.
- **Reporting & Aggregation:** Fluffy-style consolidated multi-tab Excel reports and MultiQC HTML dashboard.

## Assays

What runs for a sample is set by its `assay` (`--assay`, or the `assay` column of the barcode map):

| assay | small variants | SV | CNV | confined to the design BED |
|---|---|---|---|---|
| `wgs` | yes | yes | yes | no |
| `adaptive` | yes | yes | yes | small variants and target QC |
| `adaptive_cnv_sv` | no | yes | yes | target QC |
| `cdna` | no | no | no | – |
| `bcr_abl` | BCR-ABL1 module | no | no | yes |

New assays are declared under `assays:` in the config; see `DEVELOPMENT.md`.

## Target Region & Adaptive Sampling Support

Target BED files (e.g. for ReadUntil or targeted panels) can be provided with `--bed-file /path/to/panel.bed` or placed directly in the input directory where they are automatically discovered.

The default design is `dic_9_20_iamp21` (hg38 and T2T): four blocks on chr9, chr20, chr21 and chr22, about 110 Mb.

For a target-restricted assay with a design BED:
- **Small-variant callers** are confined to the design, and only the chromosomes it touches are run.
- **Mosdepth** computes on-target coverage distributions and threshold fractions (10x, 20x, 30x, 50x, 100x), and the report gives fold enrichment over the genome.
- **Cramino** is also run on the target BAM, since genome-wide read statistics on an adaptive run mostly describe rejected reads.
- A **target BAM** (design ±25 kb) is written for IGV.
- **SV and CNV calling stay genome-wide**, so a partner breakpoint outside the design is not lost. CNV is called from read counts per bin, not from base coverage. Segments and SVs are annotated with the design regions they overlap.

Modified-base tags (`MM`/`ML`) are carried through alignment into the BAM for viewing in IGV; no methylation calling is done.

## Resource Scaling

On workstations and local compute, Nagini dynamically scales compute-heavy sequential rules (`minimap2`, `clairs_to_call`, `sniffles2_call`, `vep`) to the full allocated core budget (e.g. `--cores 48`), avoiding disk I/O thrashing from over-parallelized alignment while maximizing single-sample processing throughput.

## Results

The workflow produces:

- **Alignments:** Sorted, indexed BAM files under `alignment/merged/`;
- **Quality Control:** Cramino read metrics, Mosdepth WGS/target coverage under `qc/`, and an aggregated `qc/multiqc/multiqc_nanopore.html` dashboard;
- **Small Variants:** Filtered, normalized VCFs under `snv_indels/`;
- **Structural Variants:** Sniffles2 (ordinary and mosaic), SAVANA and Rindels VCFs under `cnv_sv/` and `cnv/rindels/`, with SV summaries and translocation tables under `reports/sv_summary*/`;
- **Copy Number:** Rindels segments and PDF under `cnv/rindels/`, CNVkit `.cnr`/`.cns` under `cnv/cnvkit/`, and the summary, ideogram and HTML report under `reports/`;
- **Excel Report:** One consolidated Excel workbook per sample under `reports/{sample}_{type}.nagini_report.xlsx` containing:
  - **Overview:** Sample metadata, caller provenance, sequencing QC KPIs, and variant summary;
  - **SNVs and indels:** Annotated small variants with VAF, depth, Gene, Consequence, and ClinVar significance;
  - **CNVs** and **Ideogram:** Copy-number segments with cytobands and design overlap, and the whole-genome picture;
  - **SVs:** Structural variants from every caller, with cytogenetic band notation and target overlaps;
  - **Translocations:** Breakends and translocations from every caller, with cytoband ISCN notation;
  - **WGS coverage & Target coverage:** Mosdepth coverage tables;
  - **Cramino QC:** Detailed read yield and quality metrics.

Development and reproducibility notes are kept in [`DEVELOPMENT.md`](DEVELOPMENT.md).
