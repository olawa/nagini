# Poster text — BCR-ABL1 on nagini

Every number below is taken from the code and resources in this checkout, not
from memory. The commands at the bottom reproduce them.

---

## Methods (ready to paste)

The BCR-ABL1 analysis runs as an *assay* within nagini, the Snakemake pipeline
for Oxford Nanopore data at Clinical Genetics, Uppsala. One value in
`units.tsv` (`assay: bcr_abl`) selects the subgraph: 18 rules and 40 jobs,
against 70 for the WGS route. It is not a separate pipeline and not a fork —
the same rules, containers and references as the rest of the pipeline, on both
marvin and miarka.

Reads from MinKNOW are filtered before alignment on the basecaller's own mean
read quality (the `qs` tag, ≥ 20) and on length (500–3000 bp). The cDNA
amplicon spans introns, so alignment uses `minimap2 -ax splice` against
GRCh37/hg19.

Three independent questions are asked of the same BAM:

**Which fusion?** Reads spanning the breakpoint are counted per isoform. Five
isoforms are distinguished: e13a2 (b2a2) and e14a2 (b3a2), both p210 in CML;
e1a2 (p190, ALL); e19a2 (p230, CNL); and atypical e14a3 (p205). Each read is
one molecule, so the fraction is a direct count rather than a model estimate.

**Which resistance?** 73 Branford variants at 52 genomic positions across 38
codons (244–499) of the ABL1 kinase domain are interrogated. Variants are
reported from 1 % allele fraction, and subclonally down to 0.2 %.

**The same molecule?** A single ONT read spans the entire kinase domain
(~1 kb), so when two resistance mutations are found, cis versus trans is read
off individual reads — without statistical phasing and without a haplotype
assumption. The distinction is clinical: a compound mutation in cis changes the
TKI choice, two clones in trans do not.

Results are delivered as an xlsx workbook built with xlsxwriter in the same
shape as the other CGU pipelines, with VCF, TSV and JSON alongside for
downstream use.

---

## Specification

| | |
|---|---|
| **Pipeline** | nagini, Snakemake 7.32.4 (hydra-genetics conventions) |
| **Analysis selected by** | `assay: bcr_abl` in `units.tsv`, or `--assay bcr_abl` |
| **Scope** | 18 rules / 40 jobs (WGS route: 70) |
| **Reference** | GRCh37 / hg19 |
| **Read filter** | nanostream, `qs` ≥ 20, 500–3000 bp |
| **Alignment** | minimap2 2.26, `-ax splice -y -Y --MD` |
| **Fusion isoforms** | 5 — e13a2, e14a2, e1a2, e19a2, e14a3 |
| **Resistance panel** | 73 variants / 52 positions / 38 codons (244–499) |
| **Limit of detection** | AF ≥ 1 %; subclonal ≥ 0.2 % |
| **Coverage required** | depth ≥ 50, ≥ 5 alternate reads |
| **Phasing** | read-backed cis/trans, ≥ 10 spanning reads |
| **Clusters** | marvin (local) and miarka (UPPMAX), same version |

### Artefact filters in the resistance caller

Every candidate must clear all of them:

| filter | threshold | catches |
|---|---|---|
| Signal-to-noise | alt / next-largest non-reference base ≥ 3.0 | positional noise |
| Strand bias | minor-strand fraction ≥ 0.05 **and** SOR ≤ 3.0 | one-sided artefacts |
| Mean base quality, alt | ≥ 14.0 | low-quality calls |
| Flanking k-mer | ≥ 60 % clean, window ±3 bp | homopolymer / indel noise |
| Mapping quality | ≥ 1 | multi-mapping |

Strand bias requires **both** signals: SOR alone flags true variants on
strand-skewed amplicons.

---

## Limitations worth putting on the poster

Leaving them out makes the poster weaker, not stronger.

- **The signal-to-noise filter is blind to context-specific error.** A
  recurrent ONT error favouring a single substitution at 1.4 % gives roughly
  the same signal-to-noise ratio as a true subclone at 1.4 %. The filter does
  not separate them.
- **The 0.2 % threshold is not yet calibrated against negative samples.** The
  per-position background that should set it has been measured but not
  incorporated.
- **No UMI.** Fractions are raw read counts; PCR duplicates cannot be
  distinguished from independent molecules.

---

## Reproducing the numbers

```bash
# panel size
grep -vc '^#' workflow/resources/bcr_abl/branford_positions_hg19.bed   # 73

# rules and jobs for the assay
snakemake -s workflow/Snakefile --configfile <config> -n --quiet       # 40

# the isoforms
grep -A 6 'ISOFORMS' workflow/scripts/bcr_abl_fusion_counter.py
```
