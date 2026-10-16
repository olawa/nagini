# Ändringar som hör hemma uppströms

Två saker i nagini kan inte lösas med `use rule ... with:` och behöver ändras i
Hydra respektive snakemake-wrappers. Patcharna här är testade lokalt mot nagini och har skickats som pull requests:

- snakemake-wrappers: https://github.com/snakemake/snakemake-wrappers/pull/5954
- hydra-genetics/fusions: https://github.com/hydra-genetics/fusions/pull/78

När de är släppta följer stegen nedan i nagini.

| patch | repo | varför |
|---|---|---|
| `snakemake-wrappers_minimap2_samtools_fastq_extra.patch` | `snakemake/snakemake-wrappers` | Wrappern `bio/minimap2/aligner` konverterar en oalignerad BAM med `samtools fastq -T "*"`. Parametern `samtools_fastq_extra` (standard oförändrat `-T "*"`) låter nagini välja `-T MM,ML`. |
| `hydra-fusions_longread_units.patch` | `hydra-genetics/fusions` | Modulen går inte att ladda med långläsdata (`units.tsv` utan `flowcell`/`lane`), och `arriba_draw_fusion` saknar `params.extra`. |

## Varför `-T MM,ML` spelar roll

Med `-T "*"` följer alla taggar från uBAM:en med till FASTQ-headern, och
`minimap2 -y` skriver tillbaka dem. De som minimap2 själv skriver (`ms`, `nn`,
`tp`, `cm`, `s1`, `s2`, `de`, `rl`, och `ts` i splice-läge) hamnar då två gånger
på samma rad. I wrapperns eget testdata ger standardläget dubbla taggar; med
`samtools_fastq_extra="-T s1,s2"` blir det inga dubbletter. Nagini bär bara
`MM` och `ML`, och det ska den fortsätta göra.

## Så går det till

### snakemake-wrappers

1. (Gjort, se PR ovan.) Patchen gäller `master` i `snakemake/snakemake-wrappers`.
2. Testet `minimap2_ubam_sam_tags` i `bio/minimap2/aligner/test/Snakefile`
   körs med repots egna `test_wrappers.py` (kräver conda).
3. PR-titeln ska följa conventional commits, till exempel
   `feat(minimap2/aligner): add samtools_fastq_extra`. Repots mall säger det.
4. När PR:en är släppt (release-please ger en ny wrapper-tagg, till exempel
   `v9.x.y`) kan Hydra peka på den.

### hydra-genetics/alignment (efter wrapper-släppet)

I `minimap2_align`: byt wrappertaggen `v4.3.0/bio/minimap2/aligner` mot den nya
och lägg till
`samtools_fastq_extra=config.get("minimap2_align", {}).get("samtools_fastq_extra", '-T "*"')`
i `params`. Lägg även till nyckeln i modulens `config.schema.yaml`, skriv en rad
i dokumentationen och släpp en ny modultagg. I nagini ersätts då
`minimap2_unmapped_bam_align` av `use rule minimap2_align from alignment` med
`samtools_fastq_extra: "-T MM,ML"`, och wrappertaggen läggs till i
`build/vendor_hydra_dependencies.sh`.

Under tiden behåller nagini den lokala regeln. Den pipar `samtools fastq` rakt
in i minimap2. En tillfällig FASTQ-fil mellan stegen skulle bli lika stor som
hela readmängden.

### hydra-genetics/fusions

1. (Gjort, se PR ovan.) Patchen gäller `develop` i `hydra-genetics/fusions`:
   (patchen finns i `docs/upstream/hydra-fusions_longread_units.patch`).
2. Ändringen följer samma mönster som alignment-modulen: långläsdata (`ONT`,
   `PACBIO`) indexeras på `processing_unit` och `barcode`, och `units.schema.yaml`
   kräver fälten per plattform.
3. Släpp en ny modultagg. Därefter kan nagini lägga till modulen och använda
   `arriba_draw_fusion` med `params.extra` för `--sampleName`,
   `--transcriptSelection` och `--plotPanels`.

Det som återstår i nagini efter det är referensfilerna. Hydras regel kräver
cytobands, proteindomäner och GTF som indatafiler, medan nagini i dag hämtar
dem ur Arriba-containern. Den lokala regeln `bcr_abl_draw_fusions` behövs
dessutom kvar för PNG-bilden och Python-ritningen som reservväg, eftersom
Hydras regel bara skriver en PDF.
