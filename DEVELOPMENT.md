# quicken

`quicken` är en Hydra-Genetics/Snakemake-pipeline för Oxford Nanopore-data.
Den stöder FASTQ, oalignerad BAM och alignad BAM, referensbyggena hg38 och
T2T/hs1, minimap2-alignering, coverage-QC, filtrerade och normaliserade
Mutect2-anrop, Sniffles2, valfri VEP-annotering och en samlad Excelrapport.

Författare: Ola Wallerman — `ola.wallerman@scilifelab.uu.se` — GitHub `olawa`.

## Viktigaste designbeslutet: tre inputvägar

`units.tsv` har en explicit `input_type`:

| `input_type` | Pipelineväg | Realignering? |
| --- | --- | --- |
| `fastq` | minimap2 → sorterad BAM → index | Ja |
| `unmapped_bam` | samtools fastq → minimap2 → sorterad BAM → index | Ja |
| `aligned_bam` | samtools sort → index | Nej |

En alignad BAM skickas alltså aldrig till minimap2. Sortering/indexering är
bara normalisering av filen för downstream-regler; den ändrar inte mappingen.

## Hydra-kommandon som användes

Pipeline-scaffolden skapades med:

```bash
hydra-genetics create-pipeline \
  -n quicken \
  -d "Nanopore pipeline for alignment, QC, small variants and structural variants" \
  -a "Ola Wallerman" \
  -e "ola.wallerman@scilifelab.uu.se" \
  -g olawa \
  --version 0.0.1 \
  --min-snakemake-version 7.8.0 \
  -o quicken
```

Egna regler scaffoldades med Hydra `create-rule` och fylldes sedan med
pipeline-specifik logik:

```bash
hydra-genetics create-rule -c stage_aligned_bam -m quicken \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c minimap2_fastq_align -m quicken \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c minimap2_unmapped_bam_align -m quicken \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c merge_bam -t samtools -m quicken \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c mosdepth_genome -m quicken \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c mosdepth_targets -m quicken \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c sv_summary -m quicken \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
```

Kommandona ovan behöver köras från katalogen som innehåller pipelinekatalogen
och med `-o` satt till den katalogen. De är dokumenterade här främst för att
visa hur Hydra scaffoldar en regel; efter scaffolden ska placeholders alltid
ersättas och schema-/dokumentationsfiler uppdateras.

## Indata

Redigera `config/samples.tsv` och `config/units.tsv`. En ONT-rad behöver minst:

```text
sample  type  platform  processing_unit  barcode  input_type  fastq1  fastq2  bam  methylation  run_id  basecalling_model
```

Använd `fastq1` för FASTQ, `bam` för `unmapped_bam`/`aligned_bam`, och sätt
oanvända fält till exempelvis `UNUSED` (inte `NA`, eftersom pandas annars
tolkar det som saknat värde). `barcode` används som teknisk unit-nyckel och
kan vara `NB` om biblioteket inte har en biologisk barcode.

Kolumnen `methylation` dokumenterar om MM/ML-taggar förväntas finnas; den
skapar inte taggar själv. `minimap2_align.extra` skickas vidare verbatim. Den
generella konfigurationen använder `-y -Y --MD`: kommentarer från FASTQ
kopieras, supplementary alignments soft-clippas och MD-taggen skrivs. För
oalignerad BAM försöker `samtools_fastq.extra` föra över `MM,ML` till FASTQ-
kommentarerna. Ändra dessa parametrar när basecaller-/modbaseformatet kräver
det och verifiera alltid MM/ML-taggarna i den färdiga BAM-filen.

`reference_build` i config väljer mellan `references.hg38` och
`references.hs1`. hs1-sökvägarna är avsiktligt placeholders: jag hittade ingen
verifierbar T2T/hs1-sökväg i de befintliga arbetsytans konfigurationer.

## Vad som produceras

* `alignment/quicken_bam/`: en sammanslagen, sorterad BAM och `.bai` per sample/type.
* `qc/cramino/`: läs-/sekvensmängd och read-statistik.
* `qc/mosdepth/`: whole-genome coverage.
* `qc/mosdepth_bed/`: coverage mot `reference.design_bed`, inklusive per-base BED.
* `cnv_sv/sniffles2_call/`: Sniffles2 VCF, `.tbi` och `.snf`.
* `reports/sv_summary/`: TSV/JSON och separat translocations/BND-tabell.
* `snv_indels/bcftools_concat/`: sammanfogad, ofiltrerad per-kromosom Mutect2-VCF.
* `snv_indels/mutect2_postprocess/`: FilterMutectCalls-resultat,
  filtreringsstatistik, normaliserad VCF och, när aktiverat, VEP-annoterad VCF.
* `reports/{sample}_{type}.quicken_report.xlsx`: översikt, SNV/indel, SV/BND,
  WGS-/target-coverage och Cramino-QC i separata blad.
* `qc/multiqc/`: kombinerad QC-rapport.

VEP aktiveras med `annotate_small_variants_with_vep: true` och kräver att
`vep.vep_cache` pekar på en offline-cache som är kompatibel med den valda
referensen. Excelrapporten använder den VEP-annoterade VCF:en när grenen är
aktiv och annars den filtrerade, normaliserade VCF:en. Den äldre SnpEff-grenen
kan fortfarande aktiveras separat med `annotate_small_variants`; SV-annotering
styrs av `annotate_svs`. `config/config_marvin.yaml` aktiverar VEP och
SV-annotering med de referenssökvägar som används i Fluffy-miljön.

## Kör lokalt/offline

Första gången krävs nätverk för att lägga Hydra-moduler och wrappers lokalt:

```bash
cd quicken
./build/vendor_hydra_dependencies.sh
```

Scriptet lägger de pinnade modulerna i `vendor/hydra-genetics/` och
Snakemake-wrappers commit `d3e6a345acb4e03ac0662ade94d9970bb233bdab` i
`vendor/snakemake-wrappers/`. Därefter kan Snakefile och
moduler läsas lokalt. `snv_indels` får en lokal kompatibilitetstagg
`v2.1.0-quicken`, som bygger på Hydra `v2.1.0` men inte inkluderar den
oanvända DeepTrio-regeln som inte kan parsas av den här Snakemake-versionen.

Referenser och containerimages laddas inte automatiskt av pipelinekörningen.
För Singularity/Apptainer kan Hydra användas vid förberedelsen:

```bash
hydra-genetics prepare-environment create-singularity-files \
  -c config/config_marvin.yaml \
  -o vendor/containers
hydra-genetics prepare-environment container-path-update \
  -c config/config_marvin.yaml \
  -n config/config_marvin.local.yaml \
  -p vendor/containers
```

När containers och referensfiler väl finns lokalt är själva DAG-bygget offline.
För en lokal körning med containercache:

```bash
WRAPPER_PREFIX="git+file://$(pwd)/vendor/snakemake-wrappers"
snakemake \
  --profile profiles/local \
  --configfile config/config.yaml \
  --wrapper-prefix "${WRAPPER_PREFIX}"
```

Om Ubuntu-installationen saknar Singularity/Apptainer kan samma conda-fallback
användas med `--use-conda` och utan `profiles/local`; containrar är dock den
rekommenderade reproducerbara vägen när de redan är förberedda.

## Marvin

```bash
WRAPPER_PREFIX="git+file://$(pwd)/vendor/snakemake-wrappers"
snakemake \
  --profile profiles/marvin \
  --configfile config/config_marvin.yaml \
  --wrapper-prefix "${WRAPPER_PREFIX}"
```

`config/config_marvin.yaml` använder de hg38-vägar som redan förekommer i
Fluffy/CGU-konfigurationerna och resursprofilen följer samma Marvin/DRMAA-
mönster. Kontrollera alltid att just ditt Marvinprojekt har rätt bind mounts,
konto och partition innan fullskalig körning.

## Mutect2 på ONT

Mutect2 är med eftersom det uttryckligen efterfrågades och Hydra har en färdig
regel. Denna pipeline använder Hydra-regeln per kromosom, dess BED-split och
Hydras `bcftools_concat`/`tabix`. Per-kromosom-statistik slås därefter ihop med
`MergeMutectStats`, varianterna filtreras med `FilterMutectCalls` och
normaliseras med `bcftools norm`. Mutect2-delen ska ändå betraktas som
experimentell för ONT tills den är validerad mot en ONT-lämplig felmodell,
provtagningsdesign och truth set. Excelrapporten är en teknisk sammanställning,
inte ett kliniskt validerat utlåtande. Sniffles2 är den primära SV-kallaren.

## Verifiering

En liten test-DAG finns under `.tests/integration`. Den kan kontrolleras utan
att starta verktyg:

```bash
snakemake --dry-run \
  -s workflow/Snakefile \
  --configfile .tests/integration/config/config.yaml \
  --cores 1 \
  --wrapper-prefix "git+file://$(pwd)/vendor/snakemake-wrappers"
```

Inputväxeln testas separat med `--config units=.tests/integration/units_aligned.tsv`
respektive `units_unmapped.tsv`. För granskning av hela beroendegrafen finns
`workflow/rulegraph.dot`.
