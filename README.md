# Nagini

Nagini is an Oxford Nanopore sequencing analysis pipeline. It turns raw or
pre-aligned reads into normalized alignments, quality-control metrics, small
variant calls, structural variant calls, and a consolidated analysis report.

## Input modes

Nagini accepts four kinds of sample input:

| Input | Processing |
| --- | --- |
| MinKNOW/POD5 folder | Dorado basecalling and barcode demultiplexing are run before the analysis. |
| FASTQ | Reads are aligned to the selected reference. |
| Unmapped BAM | Reads are extracted, aligned, and normalized. |
| Aligned BAM | Existing mapping is retained; the BAM is sorted and indexed for downstream analysis. |

Dorado `bam_pass/barcodeXX` directories are supported with a separate
barcode-to-sample map. The map is used to create the Hydra-style input tables;
the original BAM paths and Dorado barcode identifiers are retained in the
`barcode_name` column.

For a MinKNOW/POD5 folder, the launcher accepts a Dorado model and barcode kit,
then schedules containerized basecalling and demultiplexing as part of the
workflow. The same barcode-to-sample map is used to create Hydra-style input
tables before the DAG is built. Basecalling logs and intermediate BAMs are kept
in the run's `basecalling/` directory.

Dorado is supplied by the pinned Singularity/Apptainer image configured for the
selected site; no local Dorado installation is required.

For example:

```bash
pipeline_start_scripts/start_nagini.sh \
  --input-type minknow \
  --inbox-path /data/minknow_run \
  --barcode-kit SQK-RBK114-24 \
  --basecalling-model sup \
  --barcode-map barcode_map.csv
```

The pipeline supports the hg38 and T2T/hs1 reference builds, configurable
target regions, and sample- and unit-level metadata.

## Analysis

For each sample, Nagini can perform:

- alignment normalization and BAM indexing;
- read and whole-genome coverage QC;
- target-region coverage QC;
- Mutect2 small-variant calling with merged statistics, filtering, and
  left-normalization;
- optional VEP annotation with ClinVar fields, or SnpEff annotation of small variants;
- Sniffles2 structural-variant calling;
- structural-variant annotation and summaries, including translocations, BND
  records, and optional cytoband notation; and
- aggregation of QC results with MultiQC.

## Results

The workflow produces:

- sorted, indexed BAM files under `alignment/nagini_bam/`;
- read and coverage metrics under `qc/`;
- filtered and normalized SNV/indel VCFs under `snv_indels/`;
- structural-variant VCFs and summaries under `cnv_sv/` and `reports/`; and
- one Excel report per sample and input type under `reports/`.

Each Excel report contains separate sheets for an overview, SNVs and indels,
structural variants, translocations, whole-genome coverage, target coverage,
and Cramino QC metrics.

The reports are technical summaries and are not clinical interpretations.

Development and reproducibility notes are kept in
[`DEVELOPMENT.md`](DEVELOPMENT.md).
