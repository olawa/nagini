# Verktyg och Hydra-regler i nagini

`nagini` använder i första hand Hydra-Genetics-regler och wrappers. Lokala
regler finns för inputbranching, en wrapper-kompatibel mosdepth-kedja,
Mutect2-efterbehandling, VEP, SV-sammanfattning och Excelrapportering.

## Lokala regler

| Regel | Verktyg/funktion | Varför lokal? |
| --- | --- | --- |
| `stage_aligned_bam` | samtools sort/index | Kräver uttrycklig no-realign-väg för alignad BAM. |
| `minimap2_fastq_align` | minimap2 + samtools | Input är FASTQ och behöver ett annat query-flöde än Hydra-modulens BAM-kolumn. |
| `minimap2_unmapped_bam_align` | samtools fastq + minimap2 + samtools | Oalignerad BAM konverteras innan mapping. |
| `samtools_merge_bam` | samtools merge/index | Slår ihop flera ONT units till sample/type. |
| `mosdepth_genome` | mosdepth via Hydra/Snakemake-wrapper | Hydra qc v0.8.0 deklarerar inkompatibla outputnamn med den pinnade wrappern; regeln producerar windowed whole-genome coverage. |
| `mosdepth_targets` | mosdepth via Hydra/Snakemake-wrapper | Producerar target-regioner och per-base BED med wrapperns `per_base`-output. |
| `sv_summary` | Python/VCF-parser | Hydra har ingen generell Sniffles2→translocations-tabell; regeln räknar SVTYPE/BND och kan överlappa BED. |
| `gatk_merge_mutect_stats` | GATK MergeMutectStats | Slår ihop statistik från den kromosomindelade Mutect2-körningen. |
| `gatk_filter_mutect_calls` | GATK FilterMutectCalls | Filtrerar den sammanfogade Mutect2-VCF:en med den sammanslagna statistiken. |
| `bcftools_norm_mutect2` | bcftools norm/index | Vänsterjusterar, normaliserar och delar multialleliska varianter. |
| `vep_annotate_mutect2` | Ensembl VEP | Valfri offline-annotering mot en referenskompatibel cache. |
| `nagini_excel_report` | Python XLSX-writer | Samlar QC, coverage, SNV/indel och SV/BND utan externa Python-paket. |

Shell-/script-reglerna har `container`, `conda`, logg, benchmark och de fem
standardresurserna. Mosdepth-reglerna använder i stället den pinnade Hydra-
wrapperns egna miljödefinition. Conda-fallback finns i `envs/nagini.yaml`,
`envs/gatk.yaml` och `envs/vep.yaml`.

## Hydra-moduler

| Modul | Tag i nagini | Användning |
| --- | --- | --- |
| alignment | `v0.9.0` | minimap2-index |
| cnv_sv | `v3.4.0` | Sniffles2, bgzip, tabix |
| qc | `v0.8.0` | cramino, mosdepth, target-BED mosdepth, MultiQC |
| snv_indels | `v2.1.0-nagini` | BED-split, Mutect2, bcftools concat, tabix |
| annotation | `v1.5.0` | SnpEff och simple SV annotation; VEP-kedjan är lokal |

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
