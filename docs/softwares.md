# Verktyg och Hydra-regler i nagini

`nagini` använder i första hand Hydra-Genetics-regler och wrappers. Lokala
regler finns för inputbranching, SV-sammanfattning och Excelrapportering.

## Lokala regler

| Regel | Verktyg/funktion | Varför lokal? |
| --- | --- | --- |
| `stage_aligned_bam` | samtools sort/index | Kräver uttrycklig no-realign-väg för alignad BAM. |
| `sv_summary` | Python/VCF-parser | Hydra har ingen generell Sniffles2→translocations-tabell; regeln räknar SVTYPE/BND och kan överlappa BED. |
| `nagini_excel_report` | Python XLSX-writer | Samlar QC, coverage, SNV/indel och SV/BND utan externa Python-paket. |

Shell-/script-reglerna har `container`, `conda`, logg, benchmark och de fem
standardresurserna. Conda-fallback finns i `envs/nagini.yaml`. Regler från
Hydra-moduler använder sina containrar eller wrapperns egen miljödefinition.

## Hydra-moduler

| Modul | Tag i nagini | Användning |
| --- | --- | --- |
| alignment | `v0.9.0` | minimap2-index, `minimap2_align` (FASTQ- och uBAM-units), `minimap2_merge` (sammanslagning av units), `samtools_index` (mergad BAM) |
| cnv_sv | `v3.4.0` | Sniffles2, bgzip, tabix |
| qc | `v0.8.0` | cramino, mosdepth (`mosdepth`, `mosdepth_bed`), MultiQC |
| snv_indels | `v2.1.0-nagini` | BED-split, Mutect2, MergeMutectStats (`gatk_mutect2_merge_stats`), FilterMutectCalls (`gatk_mutect2_filter`), bcftools concat, bcftools norm (`bcftools_norm`, för Mutect2, ClairS-TO, DeepVariant och amplikon), ClairS-TO, DeepVariant, bgzip, tabix |
| annotation | `v1.5.0` | SnpEff, simple SV annotation och VEP (`vep`, för Mutect2, ClairS-TO, DeepVariant och amplikon). VEP körs med `--no_stats`, så ingen `*.vep.summary.html` skrivs (D16). |

`v2.1.0-nagini` är en lokal kompatibilitetstagg från Hydra
`snv_indels v2.1.0`; patchen tar bort en oanvänd DeepTrio-include som annars
ger parsefel med Snakemake 7.32.4. Patchfilen finns i
`build/patches/snv_indels_v2.1.0_snakemake7.patch`.

## [mosdepth_genome](url_to_tool)
Introduction to mosdepth_genome

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__mosdepth_genome__mosdepth_genome#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__mosdepth_genome__mosdepth_genome#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__mosdepth_genome#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__mosdepth_genome#
argets#
