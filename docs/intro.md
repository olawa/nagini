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

Hydra-reglerna används för minimap2-index, coverage, Mutect2, Sniffles2,
komprimering/indexering, annotering och MultiQC. Lokala regler hanterar de tre
inputvägarna, filtrering och normalisering av Mutect2-resultat, VEP-kopplingen,
SV-sammanfattningen och Excelrapporten. Den genererade regelgrafen finns i
[`workflow/rulegraph.dot`](../workflow/rulegraph.dot).
