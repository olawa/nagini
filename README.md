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
  --inbox-path /data/analysis/run_01/alignment/nagini_bam \
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
- **Somatic Small Variant Calling:**
  - **ClairS-TO** (default): Deep-learning somatic SNV and small-indel caller optimized for ONT R10.4.1 data with automatic platform model selection (e.g. `ont_r10_dorado_hac_4khz`, `ont_r10_dorado_sup_5khz_ssrs`) and target BED restriction.
  - **GATK Mutect2**: Somatic small-variant calling with merged statistics, filtering, and left-normalization.
  - Caller selection is configurable via `--small-variant-caller clairs_to|mutect2|both`.
- **Structural Variant Calling:** Sniffles2 SV calling with tandem-repeat annotations, translocation/BND extraction, cytogenetic band notation (ISCN-style), and target BED overlap.
- **Annotation:** Optional VEP annotation with ClinVar pathogenicity fields and SnpEff.
- **Reporting & Aggregation:** Fluffy-style consolidated multi-tab Excel reports and MultiQC HTML dashboard.

## Target Region & Adaptive Sampling Support

Target BED files (e.g. for ReadUntil or targeted panels) can be provided with `--bed-file /path/to/panel.bed` or placed directly in the input directory where they are automatically discovered.

When a target BED is configured:
- **ClairS-TO** and **Mutect2** restrict variant calling strictly to the target regions (`--bed_fn` / `-L`).
- **Mosdepth** computes on-target coverage distributions and threshold fractions (10x, 20x, 30x, 50x, 100x).
- **SV Summaries & Excel Reports** annotate and highlight structural variants and translocations overlapping the target genes/regions.

## Resource Scaling

On workstations and local compute, Nagini dynamically scales compute-heavy sequential rules (`minimap2`, `clairs_to_call`, `sniffles2_call`, `vep`) to the full allocated core budget (e.g. `--cores 48`), avoiding disk I/O thrashing from over-parallelized alignment while maximizing single-sample processing throughput.

## Results

The workflow produces:

- **Alignments:** Sorted, indexed BAM files under `alignment/nagini_bam/`;
- **Quality Control:** Cramino read metrics, Mosdepth WGS/target coverage under `qc/`, and an aggregated `qc/multiqc/multiqc_nanopore.html` dashboard;
- **Small Variants:** Filtered, normalized VCFs under `snv_indels/`;
- **Structural Variants:** Sniffles2 VCFs, SV summaries, and translocation tables under `cnv_sv/` and `reports/sv_summary/`;
- **Excel Report:** One consolidated Excel workbook per sample under `reports/{sample}_{type}.nagini_report.xlsx` containing:
  - **Overview:** Sample metadata, caller provenance, sequencing QC KPIs, and variant summary;
  - **SNVs and indels:** Annotated small variants with VAF, depth, Gene, Consequence, and ClinVar significance;
  - **SVs:** Structural variants with cytogenetic band notation and target gene overlaps;
  - **Translocations:** Breakends and translocations with cytoband ISCN notation;
  - **WGS coverage & Target coverage:** Mosdepth coverage tables;
  - **Cramino QC:** Detailed read yield and quality metrics.

Development and reproducibility notes are kept in [`DEVELOPMENT.md`](DEVELOPMENT.md).
