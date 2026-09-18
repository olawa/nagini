# Beslut: PCR-amplikonläget (FLT3/TP53) — 2026-09-14

`harmonized_pipelines/DECISIONS.md` nämner inte nagini, så besluten nedan ligger
här tills nagini harmoniseras. Alla ändrar **levererade filer eller
analysresultat** och är därför märkta **[bekräfta]** enligt arbetsreglerna.

## A1. FLT3-ITD härleds ur långläsnings-VCF, inte med `rindels itd` **[bekräfta]**

`rindels itd` är en kortläsningscaller och används inte på ONT-data. Regeln
`amplicon_rindels_itd` är borttagen. ITD kallas i stället ur insertioner i
`snv_indels/rindels/{sample}_{type}.amplicon.vcf.gz` inom FLT3 exon 13–15
(`chr13:28033864-28034426`, från `amplicon_targets.bed`).

- Tröskel: insertion ≥ `amplicon.flt3_itd_min_length` (6 bp) och VAF ≥
  `amplicon.flt3_itd_min_vaf` (1,0 %), FILTER = PASS.
- AR = alt-läsningar / (djup − alt-läsningar). Läsram = längd % 3.
- Sökregionen är exonintervallet, inte hela 3 kb-amplikonet: en insertion
  utanför duplikationshotspotten är ingen ITD.

Leveransfilen `reports/amplicon/{sample}_{type}.flt3_itd.tsv` finns kvar men har
nya kolumner (`dup_len`, `af_pct`, `ar_category`, `filter`, `depth`, ...).

## A2. Tredje ITD-status: `NOT_ASSESSED` **[bekräfta]**

Tidigare skrev regeln `status=NEGATIVE` när callern kraschade, så ett
verktygsfel levererades som ett negativt svar. Nu gäller:

| status | betyder |
|---|---|
| `DETECTED` | ITD kallad |
| `NEGATIVE` | VCF fanns, ingen kvalificerande insertion |
| `NOT_ASSESSED` | ingen VCF — analysen gjordes aldrig |

`NOT_ASSESSED` visas som **REVIEW** i rapporten, aldrig som PASS.

## A3. Korrigerade hotspot-koordinater (GRCh38) **[bekräfta]**

Bakgrundsbrus och mutationsmatchning läste tidigare fel baser:

| Hotspot | Före | Efter |
|---|---|---|
| TP53 p.R175 | chr17:7675107 (kodon 169) | **chr17:7675088** |
| IDH2 p.R172 | chr15:90088461, "Exon 5" | **chr15:90088606, Exon 4** |
| FLT3 p.D835 | chr13:28018504, ref G | **chr13:28018505, ref C** |
| NPM1 p.W288 | chr5:171410544 | **chr5:171410542** (kodon 288 c.862) |

R248, R273 och IDH1 R132 var redan rätt.

## A4. Hotspot-matchning sker på kodonspann, inte ±2 bp **[bekräfta]**

En variant kopplas till en hotspot om den överlappar kodonets spann; indelar
vidgas med sin egen längd så att vänsterjusterade representationer träffar.
Gen-annotationen måste stämma när VCF:en bär en. NPM1c+-duplikationen har ett
bredare spann (`171410536-171410546`) eftersom den normaliseras några baser
tidigare. Statusraden skriver ut variantens **egen** koordinat, så att en
avvikelse syns i stället för att döljas.

## A5. Filtrerade varianter separeras från den kliniska tabellen **[bekräfta]**

Icke-PASS-anrop hamnade tidigare i `# Section: VARIANTS` och i Excel-avsnitt 2
utan avskiljning. De ligger nu i `# Section: FILTERED_VARIANTS` respektive
avsnitt "2b. FILTERED CALLS – NOT FOR CLINICAL REPORTING".

## A6. AR rapporteras utan ELN-riskklassificering **[bekräfta]**

Kolumnen hette "ELN Risk Classification" med "HIGH AR (>=0.5)". Det är ELN 2017;
ELN 2022 använder inte allelkvoten för riskstratifiering. Kolumnen heter nu
"AR Category" med värdena `AR >= 0.5` / `AR < 0.5`, och en fotnot förklarar det.

## Öppen punkt

NPM1 p.W288: kodonet är c.862_864, men den exakta ankarbasen bör kontrolleras
mot referens-FASTA:n innan brusgolvet för just den raden används kliniskt.
Matchningen av själva NPM1c+-insertionen är inte beroende av den kontrollen.
