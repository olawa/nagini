# nagini

`nagini` är en Hydra-Genetics/Snakemake-pipeline för Oxford Nanopore-data.
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
  -n nagini \
  -d "Nanopore pipeline for alignment, QC, small variants and structural variants" \
  -a "Ola Wallerman" \
  -e "ola.wallerman@scilifelab.uu.se" \
  -g olawa \
  --version 0.0.1 \
  --min-snakemake-version 7.8.0 \
  -o nagini
```

Egna regler scaffoldades med Hydra `create-rule` och fylldes sedan med
pipeline-specifik logik:

```bash
hydra-genetics create-rule -c stage_aligned_bam -m nagini \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c minimap2_fastq_align -m nagini \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c minimap2_unmapped_bam_align -m nagini \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c merge_bam -t samtools -m nagini \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c mosdepth_genome -m nagini \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c mosdepth_targets -m nagini \
  -a "Ola Wallerman" -e "ola.wallerman@scilifelab.uu.se"
hydra-genetics create-rule -c sv_summary -m nagini \
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
`references.hs1`. Den checkout-relativa local-configen innehåller placeholders,
medan workstation-configen pekar på den installerade `hs1.fa`-bundlen.

## Site-config och reference registry

Den gemensamma [`config/config.yaml`](config/config.yaml) innehåller pipeline-
logik och Pickett-liknande `{{VARIABLE}}`-placeholders. Miljöberoende värden
ligger i:

* [`config/site_configs/site_config_local.yaml`](config/site_configs/site_config_local.yaml)
  för en checkout-relativ lokal körning;
* [`config/site_configs/site_config_workstation.yaml`](config/site_configs/site_config_workstation.yaml)
  för den delade `/data/references`-mappen på workstation;
* [`config/site_configs/site_config_marvin.yaml`](config/site_configs/site_config_marvin.yaml)
  för Marvin.

De separata registren
[`config/reference_files/reference_files_local.yaml`](config/reference_files/reference_files_local.yaml),
[`config/reference_files/reference_files_workstation.yaml`](config/reference_files/reference_files_workstation.yaml)
och
[`config/reference_files/reference_files_marvin.yaml`](config/reference_files/reference_files_marvin.yaml)
beskriver samma referenser med logiskt namn, path, källa och checksumfält.
Registret är en reproducerbar provenance-snapshot; den aktiva Snakemake-
konfigurationen är den körspecifika `config.nagini.yaml` som startscriptet
skriver.

Site-konfigurationen slås ihop med den gemensamma konfigurationen före start.
Relativa paths löses mot pipeline-roten och placeholders expanderas innan
Snakemake validerar configen. `config/config_marvin.yaml` finns kvar som en
äldre kompatibilitetsfil, men nya körningar bör använda `--site marvin`.

## Vad som produceras

* `alignment/nagini_bam/`: en sammanslagen, sorterad BAM och `.bai` per sample/type.
* `qc/cramino/`: läs-/sekvensmängd och read-statistik.
* `qc/mosdepth/`: whole-genome coverage.
* `qc/mosdepth_bed/`: coverage mot `reference.design_bed`, inklusive per-base BED.
* `cnv_sv/sniffles2_call/`: Sniffles2 VCF, `.tbi` och `.snf`.
* `reports/sv_summary/`: TSV/JSON och separat translocations/BND-tabell.
* `snv_indels/bcftools_concat/`: sammanfogad, ofiltrerad per-kromosom Mutect2-VCF.
* `snv_indels/mutect2_postprocess/`: FilterMutectCalls-resultat,
  filtreringsstatistik, normaliserad VCF och, när aktiverat, VEP-annoterad VCF.
* `reports/{sample}_{type}.nagini_report.xlsx`: översikt, SNV/indel, SV/BND,
  WGS-/target-coverage och Cramino-QC i separata blad.
* `qc/multiqc/`: kombinerad QC-rapport.

VEP aktiveras med `annotate_small_variants_with_vep: true` och kräver att
`vep.vep_cache` pekar på en offline-cache som är kompatibel med den valda
referensen. Sätt `vep.skip: true` under tiden cachen byggs. Om den valda
referensen dessutom har en läsbar, tabix-indexerad `clinvar.vcf.gz` skickas
ClinVar-fälten `CLNSIG`, `CLNREVSTAT` och `CLNDN` som VEP-custom annotation och
visas i Excelrapporten. Excelrapporten använder den VEP-annoterade VCF:en när
grenen är aktiv och annars den filtrerade, normaliserade VCF:en.

En referens kan även ha `tandem_repeats`, som skickas till Sniffles2. En
valfri `cytoband`-fil (UCSC `cytoBand.txt` eller `.gz`) ger cytobandskolumner och
ISCN-liknande notation i SV-sammanfattningen och rapporten. Den äldre
SnpEff-grenen kan fortfarande aktiveras separat med
`annotate_small_variants`; SV-annotering styrs av `annotate_svs`.

## Kör lokalt/offline

Första gången krävs nätverk för att lägga Hydra-moduler och wrappers lokalt:

```bash
cd nagini
./build/vendor_hydra_dependencies.sh
```

Scriptet lägger de pinnade modulerna i `vendor/hydra-genetics/` och
Snakemake-wrappers commit `d3e6a345acb4e03ac0662ade94d9970bb233bdab` i
`vendor/snakemake-wrappers/`. Därefter kan Snakefile och
moduler läsas lokalt. `snv_indels` får en lokal kompatibilitetstagg
`v2.1.0-nagini`, en lokal patch-tagg baserad på Hydra `v2.1.0` som inte inkluderar den
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

För en site-aware körning som samtidigt skapar Hydra-tabeller och sparar alla
konfigurationsfiler i analyskatalogen:

```bash
pipeline_start_scripts/start_nagini.sh \
  --site local \
  --inbox-path /data/run_001 \
  --input-type unmapped_bam \
  --analysis-dir /data/analysis/run_001 \
  --dry-run
```

På Linux-workstationen med den delade referensmappen används i stället
`workstation`-site-configen. Den återanvänder det befintliga GRCh38-minimap2-
indexet och kör första provet utan VEP medan cachen färdigställs:

```bash
pipeline_start_scripts/start_nagini.sh \
  --site workstation \
  --inbox-path /data/run_001 \
  --input-type unmapped_bam \
  --analysis-dir /data/analysis/run_001 \
  --dry-run
```

När `/data/references/vep_cache` är komplett ändras
`annotate_small_variants_with_vep` till `true` och `vep.skip` till `false` i
`config/site_configs/site_config_workstation.yaml` för att aktivera VEP och
ClinVar i samma körning.

Om Ubuntu-installationen saknar Singularity/Apptainer kan samma conda-fallback
användas med `--use-conda` och utan `profiles/local`; containrar är dock den
rekommenderade reproducerbara vägen när de redan är förberedda.

## Marvin

Site-aware launcher (rekommenderad):

```bash
pipeline_start_scripts/start_nagini.sh \
  --site marvin \
  --inbox-path /data/run_001 \
  --input-type unmapped_bam \
  --analysis-dir /projects/wp2/nagini/run_001 \
  --dry-run
```

Detta använder gemensam `config/config.yaml`, Marvin-site-config,
Marvin-reference registry och `profiles/marvin`. Den äldre direkta varianten
med `config/config_marvin.yaml` fungerar fortfarande för kompatibilitet.

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

## Startscript och Hydra-input

`pipeline_start_scripts/start_nagini.sh` följer samma mönster som CGU:s
övriga `pipeline_start_scripts`: den tar en inbox, skapar eller återanvänder
`samples.tsv` och `units.tsv` i analyskatalogen, skriver en körspecifik config
och startar Snakemake med vald site-profil.

`--site local`, `--site workstation` och `--site marvin` väljer automatiskt site-config,
reference registry och Snakemake-profil. De kan ersättas individuellt med
`--site-config`, `--reference-files` och `--profile`. `--library-root` kan
anges för att dokumentera en fysisk `cgu/library`-rot i körmetadata medan
library-migreringen fortfarande är valfri.

För ONT-BAM anropas Hydra-Genetics `create-input-files`, som läser metadata ur
BAM-read groups. Scriptet kompletterar därefter tabellen med Naginis
`input_type` och oanvända FASTQ-kolumner. För FASTQ skapas samma Hydra-stil
direkt; ange en regex som fångar sample-namnet om filnamnet inte följer
standardmönstret.

### Dorado `bam_pass` och barcode-mappning

När Dorado redan har producerat `bam_pass/barcodeXX/*.bam` ska Nagini inte
köra om basecalling eller demultiplexing. Ange i stället en CSV- eller TSV-fil
med minst kolumnerna `barcode` och `sample`; `type` är valfri och använder
`--sample-type` som standard:

```text
barcode,sample,type
barcode01,SAMPLE_A,N
barcode02,SAMPLE_B,N
```

Starta med:

```bash
pipeline_start_scripts/start_nagini.sh \
  --site workstation \
  --inbox-path /data/run_001/bam_pass \
  --input-type unmapped_bam \
  --barcode-map /data/run_001/barcode_map.csv \
  --analysis-dir /data/analysis/run_001 \
  --dry-run
```

Barcode-ID:n kan anges som `01`, `bc01` eller `barcode01` och normaliseras till
`barcode01`. I Hydra-tabellen används `BC` som tekniskt ONT-barcode och det
ursprungliga Dorado-ID:t sparas i `barcode_name`. Alla BAM-filer under ett barcode blir separata units för samma
sample; Nagini slår ihop dem senare. Ommappade barcode-mappar ger ett fel för
att undvika tyst dataförlust. `unclassified` ignoreras eftersom den inte är en
sample-mappning.

Exempel:

```bash
# Oalignerad ONT-BAM; metadata läses från BAM-headern
pipeline_start_scripts/start_nagini.sh \
  --inbox-path /data/run_001 \
  --input-type unmapped_bam \
  --analysis-dir /data/analysis/run_001

# FASTQ; sample fångas från filnamnets första del före underscore
pipeline_start_scripts/start_nagini.sh \
  --inbox-path /data/fastq \
  --input-type fastq \
  --sample-regex '^(?P<sample>[^_]+)_' \
  --analysis-dir /data/analysis/run_001

# Bara kontrollera DAG:en
pipeline_start_scripts/start_nagini.sh \
  --inbox-path /data/run_001 \
  --input-type unmapped_bam \
  --dry-run
```

`--input-type aligned_bam` behåller befintlig mapping; `unmapped_bam` gör
realignering; `fastq` kör minimap2-alignering. Befintliga tabeller återanvänds
om inte `--force-inputs` anges. Scriptet skriver `config.nagini.yaml`,
`site_config.yaml`, `reference_files.yaml` och `run_metadata.yaml` i
analyskatalogen, så att en körning kan granskas eller startas om utan att
förlita sig på den aktuella arbetskatalogen.
