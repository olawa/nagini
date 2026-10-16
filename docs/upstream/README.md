# Ändringar som hör hemma uppströms

Patcharna här är testade lokalt mot nagini och har skickats som pull requests:

- snakemake-wrappers: https://github.com/snakemake/snakemake-wrappers/pull/5954
- hydra-genetics/fusions: https://github.com/hydra-genetics/fusions/pull/78

| patch | repo | varför |
|---|---|---|
| `snakemake-wrappers_minimap2_samtools_fastq_extra.patch` | `snakemake/snakemake-wrappers` | Valfri parameter `samtools_fastq_extra` i `bio/minimap2/aligner` (standard oförändrat `-T "*"`). **Inte nödvändig för nagini** (se nedan). |
| `hydra-fusions_longread_units.patch` | `hydra-genetics/fusions` | Modulen går inte att ladda med långläsdata (`units.tsv` utan `flowcell`/`lane`), och `arriba_draw_fusion` saknar `params.extra`. |

## minimap2-wrappern och taggar från uBAM

Wrappern `bio/minimap2/aligner` konverterar en oalignerad BAM med
`samtools fastq -T "*"`, så alla taggar följer med till den alignerade BAM:en.
Nagini använde tidigare bara `-T MM,ML`.

Nagini kör nu uBAM-units genom Hydras `minimap2_align` (beslut D17). Skillnaden
mot den gamla lokala regeln är de extra instrumenttaggarna (till exempel `qs`,
`ts`, `ch`, `du`); alignmentfälten och `MM`/`ML` är desamma. Jämfört på en uBAM
med sådana taggar gav båda reglerna identiska alignments, och inga taggar
dubbleras i `map-ont`-läge.

Parametern `samtools_fastq_extra` behövs därför bara om man vill hålla borta
instrumenttaggar (till exempel en stor `mv`-tabell från `--emit-moves`), eller
om en uBAM redan bär taggar som minimap2 själv skriver (`ts` kan krocka i
splice-läge, vilket inte är kontrollerat på riktig Dorado-data). Wrapper-PR:en
är en bonus och ingenting väntar på den.

## fusions

1. Patchen gäller `develop` i `hydra-genetics/fusions` (patchen finns i
   `docs/upstream/hydra-fusions_longread_units.patch`).
2. Ändringen följer samma mönster som alignment-modulen: långläsdata (`ONT`,
   `PACBIO`) indexeras på `processing_unit` och `barcode`, och `units.schema.yaml`
   kräver fälten per plattform.
3. Efter en ny modultagg kan nagini lägga till modulen och använda
   `arriba_draw_fusion` med `params.extra` för `--sampleName`,
   `--transcriptSelection` och `--plotPanels`.

Det som återstår i nagini efter det är referensfilerna. Hydras regel kräver
cytobands, proteindomäner och GTF som indatafiler, medan nagini i dag hämtar
dem ur Arriba-containern. Den lokala regeln `bcr_abl_draw_fusions` behövs
dessutom kvar för PNG-bilden och Python-ritningen som reservväg, eftersom
Hydras regel bara skriver en PDF.
