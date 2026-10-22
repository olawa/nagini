# Postertext — BCR-ABL1 på nagini

Allt nedan är hämtat ur koden och resurserna i den här checkouten, inte ur minnet.
Siffrorna går att verifiera med kommandona längst ned.

---

## Metod (att klippa in)

BCR-ABL1-analysen körs som en *assay* i nagini, Klinisk Genetik Uppsalas
Snakemake-pipeline för Oxford Nanopore-data. Ett värde i `units.tsv`
(`assay: bcr_abl`) väljer delgrafen: 20 regler och 41 jobb, mot 72 för
WGS-vägen. Det är ingen separat pipeline och ingen avgrening — samma regler,
containrar och referenser som pipelinen i övrigt.

Läsningar från MinKNOW filtreras före alignment på basecallerns egen
medelkvalitet (`qs`-taggen, ≥ 20) och längd (500–3000 bp). cDNA-amplikonet
korsar introner, så alignment sker med `minimap2 -ax splice` mot GRCh37/hg19.

Tre oberoende frågor ställs på samma BAM:

**Vilken fusion?** Läsningar som korsar brytpunkten räknas per isoform. Fem
isoformer skiljs åt: e13a2 (b2a2) och e14a2 (b3a2), båda p210 vid KML; e1a2
(p190, ALL); e19a2 (p230, KNL); samt atypiska e14a3 (p205). Varje läsning är
en molekyl, så fraktionen är en direkt räkning och inte en modellskattning.

**Vilken resistens?** 73 Branford-varianter på 52 genomiska positioner över 38
kodon (244–499) i ABL1:s kinasdomän interrogeras. 38 av dem har en kurerad
TKI-resistensprofil; övriga rapporteras som oklassade eller synonyma, inte som
resistens. Varianter rapporteras från 1 % allelfraktion, och subklonalt ned
till 0,2 %.

**Samma molekyl?** En ONT-läsning spänner hela kinasdomänen (~1 kb), så när två
resistensmutationer hittas avgörs cis/trans genom att läsa av enskilda
läsningar — utan statistisk fasning och utan haplotypantagande. Skillnaden är
klinisk: en sammansatt mutation i cis ändrar TKI-valet, två kloner i trans gör
det inte. Avläsningen är säker när båda mutationerna ligger väl över
positionsbakgrunden; se begränsningarna.

Leverans är en xlsx-arbetsbok byggd med xlsxwriter i samma form som övriga
CGU-pipelines, med VCF, TSV och JSON bredvid för vidare bearbetning.

---

## Specifikationstabell

| | |
|---|---|
| **Pipeline** | nagini, Snakemake 7 (hydra-genetics-konventioner) |
| **Val av analys** | `assay: bcr_abl` i `units.tsv`, eller `--assay bcr_abl` |
| **Omfattning** | 20 regler / 41 jobb (WGS-vägen: 72) |
| **Referens** | GRCh37 / hg19 |
| **Läsfilter** | nanostream, `qs` ≥ 20, 500–3000 bp |
| **Alignment** | minimap2 2.26, `-ax splice -y -Y --MD` |
| **Fusionsisoformer** | 5 — e13a2, e14a2, e1a2, e19a2, e14a3 |
| **Resistenspanel** | 73 varianter / 52 positioner / 38 kodon (244–499) |
| **Kurerad TKI-profil** | 38 av de 73 varianterna |
| **Detektionsgräns** | AF ≥ 1 %; subklonalt ≥ 0,2 % |
| **Täckningskrav** | djup ≥ 50, ≥ 5 alt-läsningar |
| **Fasning** | läsbaserad cis/trans, ≥ 10 spännande läsningar |

### Artefaktfilter i resistenskallningen

Varje kandidat måste passera samtliga:

| filter | tröskel | fångar |
|---|---|---|
| Signal/brus | alt / näst största icke-referensbas ≥ 3,0 | positionsbrus |
| Strängbias | minorsträngandel ≥ 0,05 **och** SOR ≤ 3,0 | ensidiga artefakter |
| Medelbaskvalitet, alt | ≥ 14,0 | lågkvalitativa anrop |
| Flankerande k-mer | ≥ 60 % rena, fönster ±3 bp | homopolymer-/indelbrus |
| Mappningskvalitet | ≥ 1 | multimappning |

Strängbias kräver **båda** signalerna: SOR ensamt flaggar äkta varianter på
strängsneda amplikon.

---

## Begränsningar som bör stå på postern

Att utelämna dem gör postern svagare, inte starkare.

- **Signal/brus-filtret är blint för kontextspecifika fel.** Ett återkommande
  ONT-fel som gynnar en enda substitution vid 1,4 % ger ungefär samma
  signal/brus-kvot som en äkta subklon vid 1,4 %. Filtret separerar dem inte.
- **Tröskeln 0,2 % är ännu inte kalibrerad mot negativa prov.** Den per-positions
  bakgrund som skulle sätta gränsen är mätt men inte inarbetad.
- **Fasningen drar inte av bakgrundsfelet.** Sekvenseringsfel på bara den ena
  positionen ger skenbara enkelmutanta läsningar från vildtypsmolekyler. En ren
  sammansatt mutation på 1 % med 0,3 % positionsfel klassas därför som
  `COEXISTING (Compound + Polyclonal)` — och på 5 % gör den det med 1,4 % fel.
  Nära detektionsgränsen ska en polyklonal tolkning inte tas för given.
- **Ingen UMI.** Fraktionerna är rålästräkningar; PCR-duplikat går inte att
  skilja från oberoende molekyler.

---

## Verifiera siffrorna

```bash
# panelens storlek
grep -vc '^#' workflow/resources/bcr_abl/branford_positions_hg19.bed   # 73

# antal regler och jobb för assayen
snakemake -s workflow/Snakefile --configfile <config> -n --quiet       # 41

# isoformerna
grep -A 6 'ISOFORMS' workflow/scripts/bcr_abl_fusion_counter.py
```
