# Quicken

Quicken is an Oxford Nanopore sequencing analysis pipeline. It turns raw or
pre-aligned reads into normalized alignments, quality-control metrics, small
variant calls, structural variant calls, and a consolidated analysis report.

## Input modes

Quicken accepts three kinds of sample input:

| Input | Processing |
| --- | --- |
| FASTQ | Reads are aligned to the selected reference. |
| Unmapped BAM | Reads are extracted, aligned, and normalized. |
| Aligned BAM | Existing mapping is retained; the BAM is sorted and indexed for downstream analysis. |

The pipeline supports the hg38 and T2T/hs1 reference builds, configurable
target regions, and sample- and unit-level metadata.

## Analysis

For each sample, Quicken can perform:

- alignment normalization and BAM indexing;
- read and whole-genome coverage QC;
- target-region coverage QC;
- Mutect2 small-variant calling with merged statistics, filtering, and
  left-normalization;
- optional VEP or SnpEff annotation of small variants;
- Sniffles2 structural-variant calling;
- structural-variant annotation and summaries, including translocations and
  BND records; and
- aggregation of QC results with MultiQC.

## Results

The workflow produces:

- sorted, indexed BAM files under `alignment/quicken_bam/`;
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
