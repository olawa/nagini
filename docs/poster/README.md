# Postermaterial — BCR-ABL1 på nagini

Figurer och text för postern, plus ett körbart exempel så att påståendet
"det är ingen separat pipeline" går att kontrollera i stället för att tros på.

| fil | vad |
|---|---|
| `fig1_architecture_sv.svg` / `_en.svg` | pipelinefiguren: från FASTQ till levererad arbetsbok |
| `fig2_analysis_sv.svg` / `_en.svg` | analysfiguren: fusionsisoformer, Branford-panelen, cis/trans |
| `metod_sv.md` / `methods_en.md` | metodstycke, specifikationstabell, begränsningar |
| `example/` | minimal konfiguration som kör BCR-ABL-vägen |

Båda språken finns i båda figurerna. De är separata filer och inte lager, så
en ändring i sak måste göras i båda — det är två filer, inte åtta.

Figurerna är handskriven SVG — vektor, så de skalar till posterstorlek utan
förlust, och texten går att redigera direkt i Illustrator eller Inkscape.
De delar palett: blått för BCR/chr22, orange för ABL1/chr9 och de analyser som
kallar varianter, grönt för fasning och leverans. Paletten är validerad för
färgseendevariation på *alla* par (värsta ΔE 9,2 deutan, 24,0 normalseende),
eftersom vilka två element som helst kan jämföras på en poster.

## Kör exemplet

```bash
.venv/bin/snakemake -s workflow/Snakefile \
    --configfile .tests/integration/config/config.yaml docs/poster/example/config.yaml \
    -n
```

Ger 40 jobb och 18 regler, varav fem BCR-ABL-specifika. Samma kommando utan
`docs/poster/example/config.yaml` ger 70 jobb — WGS-vägen. Det är hela skillnaden,
och den ligger i `assay`-kolumnen i `units.tsv`, inte i koden.

## Vad som ska stå på postern

Pinna en tagg på den commit figuren beskriver, så att SHA:t på postern
motsvarar koden som kör i drift:

```bash
git tag -a nagini-bcr_abl-poster -m "Koden som postern beskriver"
```

En tagg, inte en gren: en avgrening skulle börja identisk och sedan divergera,
och versionen du presenterar skulle sluta vara den som analyserar prov.

## Om figurerna behöver ändras

Siffrorna i dem är hämtade ur resurserna, inte inskrivna för hand. Ändras
panelen ändras figuren:

```bash
# 73 varianter / 52 positioner
grep -vc '^#' workflow/resources/bcr_abl/branford_positions_hg19.bed
awk '!/^#/{print $1":"$2}' workflow/resources/bcr_abl/branford_positions_hg19.bed | sort -u | wc -l
```

Kodonstrecken i `fig2_analysis_*.svg` är genererade ur den BED-filen. Byts panelen
ut måste figuren genereras om — den räknar inte själv.

## Språk

Svenska och engelska finns färdiga. Engelska versionen är inte en maskinell
översättning: raderna är omformulerade där engelskan blev bredare än rutan.
Två rader fick kortas i båda språken av samma skäl.

Båda versionerna är renderade och granskade, inte bara validerade som XML.
Det fångade tre kollisioner i den svenska figur 2 som inte syns i koden:
T315I-etiketten ovanpå rubriken, e1a2-linjen genom ett textblock och
cis/trans-raderna ovanpå varandra.
