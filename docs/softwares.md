# Verktyg och Hydra-regler i nagini

`nagini` använder i första hand Hydra-Genetics-regler och wrappers. Lokala
regler finns där ingen Hydra-regel gör jobbet: Dorado, Rindels, BCR-ABL1,
CNV-rapporten, inputbranching, SV-sammanfattning och Excelrapportering.

## Hydra-moduler

| Modul | Tag i nagini | Användning |
| --- | --- | --- |
| alignment | `v0.9.0` | minimap2-index, `minimap2_align` (FASTQ- och uBAM-units), `minimap2_merge` (sammanslagning av units), `samtools_index` (mergad BAM) |
| cnv_sv | `v3.4.0` | Sniffles2, bgzip, tabix |
| qc | `v0.8.0` | cramino, mosdepth (`mosdepth`, `mosdepth_bed`), MultiQC |
| snv_indels | `v2.1.0-nagini` | BED-split, Mutect2, MergeMutectStats (`gatk_mutect2_merge_stats`), FilterMutectCalls (`gatk_mutect2_filter`), bcftools concat, bcftools norm (`bcftools_norm`, för Mutect2, ClairS-TO och DeepVariant), ClairS-TO, DeepVariant, bgzip, tabix |
| annotation | `v1.5.0` | SnpEff, simple SV annotation och VEP (`vep`, för Mutect2, ClairS-TO och DeepVariant). VEP körs med `--no_stats`, så ingen `*.vep.summary.html` skrivs (D16). |

`v2.1.0-nagini` är en lokal kompatibilitetstagg från Hydra
`snv_indels v2.1.0`; patchen tar bort en oanvänd DeepTrio-include som annars
ger parsefel med Snakemake 7.32.4. Patchfilen finns i
`build/patches/snv_indels_v2.1.0_snakemake7.patch`.

## Lokala regler

| Regel | Fil | Funktion | Varför lokal? |
| --- | --- | --- | --- |
| `align_dorado_bam_raw` | `align_dorado_bam.smk` | Align the exported per-barcode FASTQ with minimap2. | Dorado-justering per barcode; hydras `minimap2_align` eller Dorados egen aligner väljs per unit. |
| `finalize_dorado_bam_alignment` | `align_dorado_bam.smk` | Add sample metadata, coordinate-sort, and index an aligned Dorado unit. | Dorado-justering per barcode; hydras `minimap2_align` eller Dorados egen aligner väljs per unit. |
| `bcr_abl_fusion_counts` | `bcr_abl.smk` | Identify and quantify BCR-ABL1 fusion isoforms from spliced alignments. | BCR-ABL1-specifik analys utan Hydra-motsvarighet. |
| `bcr_abl_branford_calling` | `bcr_abl.smk` | Force-call BCR-ABL1 kinase domain resistance mutations at high depth. | BCR-ABL1-specifik analys utan Hydra-motsvarighet. |
| `bcr_abl_phasing` | `bcr_abl.smk` | Single-molecule phasing of compound resistance mutations from long reads. | BCR-ABL1-specifik analys utan Hydra-motsvarighet. |
| `bcr_abl_exon_coverage` | `bcr_abl.smk` | Calculate depth and coverage metrics across BCR and ABL1 exons. | BCR-ABL1-specifik analys utan Hydra-motsvarighet. |
| `bcr_abl_draw_fusions` | `bcr_abl.smk` | Generate Arriba fusion diagrams using draw_fusions.R. | BCR-ABL1-specifik analys utan Hydra-motsvarighet. |
| `bcr_abl_bcftools_pileup` | `bcr_abl.smk` | Run bcftools mpileup and call across BCR-ABL1 target exons at full depth. | BCR-ABL1-specifik analys utan Hydra-motsvarighet. |
| `bcr_abl_rindels` | `bcr_abl.smk` | Run Rindels in targeted ONT BCR-ABL1 amplicon mode. | BCR-ABL1-specifik analys utan Hydra-motsvarighet. |
| `cnv_rindels` | `cnv.smk` | Run Rindels in CNV-only mode to produce whole-genome binned read counts and ASCAT segments. | Rindels är en värdbinär utan Hydra-modul. |
| `cnv_summary` | `cnv.smk` | Summarise Rindels ASCAT CNV segments and annotate with cytobands and target regions. | Rindels är en värdbinär utan Hydra-modul. |
| `cytoband_centromeres` | `cnv.smk` | Centromere intervals for the ideogram, from the build's own cytobands. | Rindels är en värdbinär utan Hydra-modul. |
| `cnv_rindels_ideogram` | `cnv.smk` | Whole-genome ideogram from the Rindels bins and BAF track. | Rindels är en värdbinär utan Hydra-modul. |
| `decompress_cytoband` | `cnv_html_report.smk` | Decompress cytobands if compressed, because merge_cnv_json expects plain text. | Hydras reports-modul kräver `general_report` i configen; skripten anropas direkt. |
| `cnvkit_json` | `cnv_html_report.smk` | Convert CNVkit .cnr / .cns output to the intermediate JSON format expected by merge_cnvkit_json. | Hydras reports-modul kräver `general_report` i configen; skripten anropas direkt. |
| `merge_cnvkit_json` | `cnv_html_report.smk` | Merge per-caller CNV JSON(s) with genome context (cytobands, FAI). | Hydras reports-modul kräver `general_report` i configen; skripten anropas direkt. |
| `compile_cnv_html_report` | `cnv_html_report.smk` | Render the final interactive D3.js HTML report from the merged JSON. | Hydras reports-modul kräver `general_report` i configen; skripten anropas direkt. |
| `cnvkit_targets` | `cnvkit.smk` | Tile the accessible genome into fixed bins, once per reference build. | Tumör-only WGS med färdigtilade bins och `--count-reads`; `cnvkit_batch` i cnv_sv kräver en normal_reference. |
| `cnvkit_batch` | `cnvkit.smk` | Run CNVkit in WGS mode on a tumor-only ONT BAM. | Tumör-only WGS med färdigtilade bins och `--count-reads`; `cnvkit_batch` i cnv_sv kräver en normal_reference. |
| `copy_to_results` | `copy_results.smk` | Copy one deliverable into results/, leaving the in-place file alone. | En generisk kopieringsregel i stället för genererade (Snakemake 9, `when:`). |
| `dorado_basecaller` | `dorado.smk` | Basecall a MinKNOW/POD5 directory and retain barcode sequence context. | Inget Hydra-modul för Dorado (`hydra-genetics/long_read` har bara Guppy och pbmm2). |
| `dorado_demux` | `dorado.smk` | Classify reads and retain Dorado's native demultiplexing directory. | Inget Hydra-modul för Dorado (`hydra-genetics/long_read` har bara Guppy och pbmm2). |
| `dorado_demux_summary` | `dorado.smk` | Copy Dorado's demultiplexing summary out of its raw output directory. | Inget Hydra-modul för Dorado (`hydra-genetics/long_read` har bara Guppy och pbmm2). |
| `dorado_demux_barcode_bam` | `dorado.smk` | Normalize Dorado output to one explicit, unaligned BAM per barcode. | Inget Hydra-modul för Dorado (`hydra-genetics/long_read` har bara Guppy och pbmm2). |
| `dorado_demux_fastq` | `dorado.smk` | Optionally export a gzip-compressed FASTQ for a demultiplexed barcode. | Inget Hydra-modul för Dorado (`hydra-genetics/long_read` har bara Guppy och pbmm2). |
| `dorado_demux_qc` | `dorado.smk` | Aggregate the per-read demux summary into a barcode yield table. | Inget Hydra-modul för Dorado (`hydra-genetics/long_read` har bara Guppy och pbmm2). |
| `nagini_excel_report` | `excel_report.smk` | Build a portable multi-sheet workbook from final QC and variant outputs. | Rapporterna är nagini-specifika. |
| `nagini_summary_excel_report` | `excel_report.smk` | Build a multi-sample summary Excel workbook across all samples in the run. | Rapporterna är nagini-specifika. |
| `merged_fastq` | `merged_fastq.smk` | One FASTQ per sample from its unaligned source BAMs, mod tags kept. | Slår ihop uBAM-källor till en FASTQ per prov med modifieringstaggar. |
| `read_filter_fastq` | `read_filter.smk` | Filter raw FASTQ reads by QV and length before alignment. | Kvalitets- och längdfilter med `nanostream`, som saknar Hydra-modul. |
| `read_filter_unmapped_bam` | `read_filter.smk` | Filter raw unmapped BAM reads by QV and length before alignment. | Kvalitets- och längdfilter med `nanostream`, som saknar Hydra-modul. |
| `create_genome_contigs_bed` | `snv_postprocess.smk` | Create a whole-genome BED from FASTA index when no design BED is provided. | Hjälpregel för en BED över hela genomet. |
| `stage_aligned_bam` | `stage_aligned_bam.smk` | Prepare an existing aligned BAM without invoking an aligner. | Kräver uttrycklig no-realign-väg för alignad BAM. |
| `sv_summary` | `sv_summary.smk` | Summarise Sniffles2 SVs and extract translocations/BNDs. | Hydra har ingen generell Sniffles2→translocations-tabell. |
| `sv_summary_rindels` | `sv_summary.smk` | The same summary over the Rindels call set, for a second opinion. | Hydra har ingen generell Sniffles2→translocations-tabell. |
| `sv_summary_savana` | `sv_summary.smk` | Summarise SAVANA somatic SVs and extract translocations/BNDs. | Hydra har ingen generell Sniffles2→translocations-tabell. |
| `pad_target_bed` | `target_bam.smk` | Pad the design BED file by +/- 25 kb (or configurable bp) and merge overlapping intervals. | Läsurval runt designregionerna för adaptiv sampling. |
| `extract_target_bam` | `target_bam.smk` | Extract reads overlapping padded target regions into a target-only BAM. | Läsurval runt designregionerna för adaptiv sampling. |

Varje lokal regel har `input`, `output`, `params`, `log`, `benchmark`, `threads`,
`resources`, `container`, `message` och en av `shell`/`script`/`wrapper`. Undantag:
`params` saknas på `dorado_demux_summary` och `create_genome_contigs_bed`, som inte har några
verktygsflaggor att skicka vidare; `conda` saknas där ingen conda-miljö kan ersätta containern (Dorado, Rindels, Arriba);
övriga regler delar `envs/nagini.yaml`. Rindels-reglerna har ingen `container`, eftersom
Rindels är en värdbinär.

## Regelreferens

Tre regler saknar referensavsnitt eftersom dokumentationspluginet inte kan läsa dem:
`align_dorado_bam_raw` (definieras villkorligt), `copy_to_results` och
`nagini_summary_excel_report` (flerradiga och `unpack`-baserade indata). De finns i tabellen ovan.

### finalize_dorado_bam_alignment
Add sample metadata, coordinate-sort, and index an aligned Dorado unit.

#SNAKEMAKE_RULE_SOURCE__align_dorado_bam__finalize_dorado_bam_alignment#

#SNAKEMAKE_RULE_TABLE__align_dorado_bam__finalize_dorado_bam_alignment#

#CONFIGSCHEMA__finalize_dorado_alignment#

#RESOURCESSCHEMA__finalize_dorado_alignment#

### bcr_abl_fusion_counts
Identify and quantify BCR-ABL1 fusion isoforms from spliced alignments.

#SNAKEMAKE_RULE_SOURCE__bcr_abl__bcr_abl_fusion_counts#

#SNAKEMAKE_RULE_TABLE__bcr_abl__bcr_abl_fusion_counts#

#CONFIGSCHEMA__bcr_abl#

#RESOURCESSCHEMA__bcr_abl#

### bcr_abl_branford_calling
Force-call BCR-ABL1 kinase domain resistance mutations at high depth.

#SNAKEMAKE_RULE_SOURCE__bcr_abl__bcr_abl_branford_calling#

#SNAKEMAKE_RULE_TABLE__bcr_abl__bcr_abl_branford_calling#

#CONFIGSCHEMA__bcr_abl#

#RESOURCESSCHEMA__bcr_abl#

### bcr_abl_phasing
Single-molecule phasing of compound resistance mutations from long reads.

#SNAKEMAKE_RULE_SOURCE__bcr_abl__bcr_abl_phasing#

#SNAKEMAKE_RULE_TABLE__bcr_abl__bcr_abl_phasing#

#CONFIGSCHEMA__bcr_abl#

#RESOURCESSCHEMA__bcr_abl#

### bcr_abl_exon_coverage
Calculate depth and coverage metrics across BCR and ABL1 exons.

#SNAKEMAKE_RULE_SOURCE__bcr_abl__bcr_abl_exon_coverage#

#SNAKEMAKE_RULE_TABLE__bcr_abl__bcr_abl_exon_coverage#

#CONFIGSCHEMA__bcr_abl#

#RESOURCESSCHEMA__bcr_abl#

### bcr_abl_draw_fusions
Generate Arriba fusion diagrams using draw_fusions.R.

#SNAKEMAKE_RULE_SOURCE__bcr_abl__bcr_abl_draw_fusions#

#SNAKEMAKE_RULE_TABLE__bcr_abl__bcr_abl_draw_fusions#

#CONFIGSCHEMA__arriba#

#RESOURCESSCHEMA__arriba#

### bcr_abl_bcftools_pileup
Run bcftools mpileup and call across BCR-ABL1 target exons at full depth.

#SNAKEMAKE_RULE_SOURCE__bcr_abl__bcr_abl_bcftools_pileup#

#SNAKEMAKE_RULE_TABLE__bcr_abl__bcr_abl_bcftools_pileup#

#CONFIGSCHEMA__bcr_abl#

#RESOURCESSCHEMA__bcr_abl#

### bcr_abl_rindels
Run Rindels in targeted ONT BCR-ABL1 amplicon mode.

#SNAKEMAKE_RULE_SOURCE__bcr_abl__bcr_abl_rindels#

#SNAKEMAKE_RULE_TABLE__bcr_abl__bcr_abl_rindels#

#CONFIGSCHEMA__bcr_abl_rindels#

#RESOURCESSCHEMA__bcr_abl_rindels#

### cnv_rindels
Run Rindels in CNV-only mode to produce whole-genome binned read counts and ASCAT segments.

#SNAKEMAKE_RULE_SOURCE__cnv__cnv_rindels#

#SNAKEMAKE_RULE_TABLE__cnv__cnv_rindels#

#CONFIGSCHEMA__cnv_rindels#

#RESOURCESSCHEMA__cnv_rindels#

### cnv_summary
Summarise Rindels ASCAT CNV segments and annotate with cytobands and target regions.

#SNAKEMAKE_RULE_SOURCE__cnv__cnv_summary#

#SNAKEMAKE_RULE_TABLE__cnv__cnv_summary#

#CONFIGSCHEMA__cnv_summary#

#RESOURCESSCHEMA__cnv_summary#

### cytoband_centromeres
Centromere intervals for the ideogram, from the build's own cytobands.

#SNAKEMAKE_RULE_SOURCE__cnv__cytoband_centromeres#

#SNAKEMAKE_RULE_TABLE__cnv__cytoband_centromeres#

#CONFIGSCHEMA__cytoband_centromeres#

#RESOURCESSCHEMA__cytoband_centromeres#

### cnv_rindels_ideogram
Whole-genome ideogram from the Rindels bins and BAF track.

#SNAKEMAKE_RULE_SOURCE__cnv__cnv_rindels_ideogram#

#SNAKEMAKE_RULE_TABLE__cnv__cnv_rindels_ideogram#

#CONFIGSCHEMA__rindels_ideogram#

#RESOURCESSCHEMA__rindels_ideogram#

### decompress_cytoband
Decompress cytobands if compressed, because merge_cnv_json expects plain text.

#SNAKEMAKE_RULE_SOURCE__cnv_html_report__decompress_cytoband#

#SNAKEMAKE_RULE_TABLE__cnv_html_report__decompress_cytoband#

#CONFIGSCHEMA__decompress_cytoband#

#RESOURCESSCHEMA__decompress_cytoband#

### cnvkit_json
Convert CNVkit .cnr / .cns output to the intermediate JSON format expected by merge_cnvkit_json.

#SNAKEMAKE_RULE_SOURCE__cnv_html_report__cnvkit_json#

#SNAKEMAKE_RULE_TABLE__cnv_html_report__cnvkit_json#

#CONFIGSCHEMA__cnv_json#

#RESOURCESSCHEMA__cnv_json#

### merge_cnvkit_json
Merge per-caller CNV JSON(s) with genome context (cytobands, FAI).

#SNAKEMAKE_RULE_SOURCE__cnv_html_report__merge_cnvkit_json#

#SNAKEMAKE_RULE_TABLE__cnv_html_report__merge_cnvkit_json#

#CONFIGSCHEMA__merge_cnv_json#

#RESOURCESSCHEMA__merge_cnv_json#

### compile_cnv_html_report
Render the final interactive D3.js HTML report from the merged JSON.

#SNAKEMAKE_RULE_SOURCE__cnv_html_report__compile_cnv_html_report#

#SNAKEMAKE_RULE_TABLE__cnv_html_report__compile_cnv_html_report#

#CONFIGSCHEMA__cnv_html_report#

#RESOURCESSCHEMA__cnv_html_report#

### cnvkit_targets
Tile the accessible genome into fixed bins, once per reference build.

#SNAKEMAKE_RULE_SOURCE__cnvkit__cnvkit_targets#

#SNAKEMAKE_RULE_TABLE__cnvkit__cnvkit_targets#

#CONFIGSCHEMA__cnvkit_batch#

#RESOURCESSCHEMA__cnvkit_batch#

### cnvkit_batch
Run CNVkit in WGS mode on a tumor-only ONT BAM.

#SNAKEMAKE_RULE_SOURCE__cnvkit__cnvkit_batch#

#SNAKEMAKE_RULE_TABLE__cnvkit__cnvkit_batch#

#CONFIGSCHEMA__cnvkit_batch#

#RESOURCESSCHEMA__cnvkit_batch#

### dorado_basecaller
Basecall a MinKNOW/POD5 directory and retain barcode sequence context.

#SNAKEMAKE_RULE_SOURCE__dorado__dorado_basecaller#

#SNAKEMAKE_RULE_TABLE__dorado__dorado_basecaller#

#CONFIGSCHEMA__dorado#

#RESOURCESSCHEMA__dorado#

### dorado_demux
Classify reads and retain Dorado's native demultiplexing directory.

#SNAKEMAKE_RULE_SOURCE__dorado__dorado_demux#

#SNAKEMAKE_RULE_TABLE__dorado__dorado_demux#

#CONFIGSCHEMA__dorado#

#RESOURCESSCHEMA__dorado#

### dorado_demux_summary
Copy Dorado's demultiplexing summary out of its raw output directory.

#SNAKEMAKE_RULE_SOURCE__dorado__dorado_demux_summary#

#SNAKEMAKE_RULE_TABLE__dorado__dorado_demux_summary#

#CONFIGSCHEMA__dorado_demux_summary#

#RESOURCESSCHEMA__dorado_demux_summary#

### dorado_demux_barcode_bam
Normalize Dorado output to one explicit, unaligned BAM per barcode.

#SNAKEMAKE_RULE_SOURCE__dorado__dorado_demux_barcode_bam#

#SNAKEMAKE_RULE_TABLE__dorado__dorado_demux_barcode_bam#

#CONFIGSCHEMA__dorado_demux_bam#

#RESOURCESSCHEMA__dorado_demux_bam#

### dorado_demux_fastq
Optionally export a gzip-compressed FASTQ for a demultiplexed barcode.

#SNAKEMAKE_RULE_SOURCE__dorado__dorado_demux_fastq#

#SNAKEMAKE_RULE_TABLE__dorado__dorado_demux_fastq#

#CONFIGSCHEMA__dorado_fastq#

#RESOURCESSCHEMA__dorado_fastq#

### dorado_demux_qc
Aggregate the per-read demux summary into a barcode yield table.

#SNAKEMAKE_RULE_SOURCE__dorado__dorado_demux_qc#

#SNAKEMAKE_RULE_TABLE__dorado__dorado_demux_qc#

#CONFIGSCHEMA__dorado_demux_qc#

#RESOURCESSCHEMA__dorado_demux_qc#

### nagini_excel_report
Build a portable multi-sheet workbook from final QC and variant outputs.

#SNAKEMAKE_RULE_SOURCE__excel_report__nagini_excel_report#

#SNAKEMAKE_RULE_TABLE__excel_report__nagini_excel_report#

#CONFIGSCHEMA__nagini_excel_report#

#RESOURCESSCHEMA__nagini_excel_report#

### merged_fastq
One FASTQ per sample from its unaligned source BAMs, mod tags kept.

#SNAKEMAKE_RULE_SOURCE__merged_fastq__merged_fastq#

#SNAKEMAKE_RULE_TABLE__merged_fastq__merged_fastq#

#CONFIGSCHEMA__samtools_fastq#

#RESOURCESSCHEMA__samtools_fastq#

### read_filter_fastq
Filter raw FASTQ reads by QV and length before alignment.

#SNAKEMAKE_RULE_SOURCE__read_filter__read_filter_fastq#

#SNAKEMAKE_RULE_TABLE__read_filter__read_filter_fastq#

#CONFIGSCHEMA__read_filter#

#RESOURCESSCHEMA__read_filter#

### read_filter_unmapped_bam
Filter raw unmapped BAM reads by QV and length before alignment.

#SNAKEMAKE_RULE_SOURCE__read_filter__read_filter_unmapped_bam#

#SNAKEMAKE_RULE_TABLE__read_filter__read_filter_unmapped_bam#

#CONFIGSCHEMA__read_filter#

#RESOURCESSCHEMA__read_filter#

### create_genome_contigs_bed
Create a whole-genome BED from FASTA index when no design BED is provided.

#SNAKEMAKE_RULE_SOURCE__snv_postprocess__create_genome_contigs_bed#

#SNAKEMAKE_RULE_TABLE__snv_postprocess__create_genome_contigs_bed#

#CONFIGSCHEMA__create_genome_contigs_bed#

#RESOURCESSCHEMA__create_genome_contigs_bed#

### stage_aligned_bam
Prepare an existing aligned BAM without invoking an aligner.

#SNAKEMAKE_RULE_SOURCE__stage_aligned_bam__stage_aligned_bam#

#SNAKEMAKE_RULE_TABLE__stage_aligned_bam__stage_aligned_bam#

#CONFIGSCHEMA__stage_aligned_bam#

#RESOURCESSCHEMA__stage_aligned_bam#

### sv_summary
Summarise Sniffles2 SVs and extract translocations/BNDs.

#SNAKEMAKE_RULE_SOURCE__sv_summary__sv_summary#

#SNAKEMAKE_RULE_TABLE__sv_summary__sv_summary#

#CONFIGSCHEMA__sv_summary#

#RESOURCESSCHEMA__sv_summary#

### sv_summary_rindels
The same summary over the Rindels call set, for a second opinion.

#SNAKEMAKE_RULE_SOURCE__sv_summary__sv_summary_rindels#

#SNAKEMAKE_RULE_TABLE__sv_summary__sv_summary_rindels#

#CONFIGSCHEMA__sv_summary#

#RESOURCESSCHEMA__sv_summary#

### sv_summary_savana
Summarise SAVANA somatic SVs and extract translocations/BNDs.

#SNAKEMAKE_RULE_SOURCE__sv_summary__sv_summary_savana#

#SNAKEMAKE_RULE_TABLE__sv_summary__sv_summary_savana#

#CONFIGSCHEMA__sv_summary#

#RESOURCESSCHEMA__sv_summary#

### pad_target_bed
Pad the design BED file by +/- 25 kb (or configurable bp) and merge overlapping intervals.

#SNAKEMAKE_RULE_SOURCE__target_bam__pad_target_bed#

#SNAKEMAKE_RULE_TABLE__target_bam__pad_target_bed#

#CONFIGSCHEMA__pad_target_bed#

#RESOURCESSCHEMA__pad_target_bed#

### extract_target_bam
Extract reads overlapping padded target regions into a target-only BAM.

#SNAKEMAKE_RULE_SOURCE__target_bam__extract_target_bam#

#SNAKEMAKE_RULE_TABLE__target_bam__extract_target_bam#

#CONFIGSCHEMA__extract_target_bam#

#RESOURCESSCHEMA__extract_target_bam#
