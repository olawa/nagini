# Pipelineöversikt

```mermaid
flowchart LR
    A[FASTQ] --> M[minimap2]
    B[Oalignerad BAM] --> F[samtools fastq]
    F --> M
    C[Alignad BAM] --> S[sort/index, ingen realignering]
    M --> S2[sort/index]
    S --> X[Sammanslagen BAM]
    S2 --> X
    X --> Q[cramino + mosdepth]
    X --> V[Mutect2 per kromosom]
    V --> VF[Merge stats + FilterMutectCalls]
    VF --> VN[bcftools norm]
    VN --> VA[VEP, valfritt]
    X --> SV[Sniffles2]
    SV --> R[SV-annotering/sammanfattning]
    VN --> E[Excelrapport]
    VA --> E
    R --> E
```

Hydra-reglerna används för minimap2-index, alignment och sammanslagning, coverage, Mutect2 och
dess filtrering, ClairS-TO, DeepVariant, Sniffles2, SAVANA, normalisering, VEP, SnpEff, komprimering/indexering
och MultiQC. Lokala regler hanterar Dorado, Rindels, BCR-ABL1, CNV-rapporten, de tre inputvägarna,
SV-sammanfattningen och Excelrapporten (se [verktyg och regler](softwares.md)). Regelgrafen finns som
[`images/rulegraph.svg`](../images/rulegraph.svg) och i källformat som
[`workflow/rulegraph.dot`](../workflow/rulegraph.dot).
