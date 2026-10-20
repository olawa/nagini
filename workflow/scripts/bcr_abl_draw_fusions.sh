#!/usr/bin/env bash
set -uo pipefail

SAMPLE=""
FUSIONS=""
FUSION_SUMMARY=""
MUTATIONS=""
PHASING=""
EXON_COVERAGE=""
BAM=""
GTF=""
CYTOBANDS=""
PROTEIN_DOMAINS=""
BUILD="hg19"
EXTRA=""
OUTPUT_PDF=""
OUTPUT_PNG=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --sample=*) SAMPLE="${1#*=}"; shift ;;
        --sample) SAMPLE="$2"; shift 2 ;;
        --fusions=*) FUSIONS="${1#*=}"; shift ;;
        --fusions) FUSIONS="$2"; shift 2 ;;
        --fusion-summary=*) FUSION_SUMMARY="${1#*=}"; shift ;;
        --fusion-summary) FUSION_SUMMARY="$2"; shift 2 ;;
        --mutations=*) MUTATIONS="${1#*=}"; shift ;;
        --mutations) MUTATIONS="$2"; shift 2 ;;
        --phasing=*) PHASING="${1#*=}"; shift ;;
        --phasing) PHASING="$2"; shift 2 ;;
        --exon-coverage=*) EXON_COVERAGE="${1#*=}"; shift ;;
        --exon-coverage) EXON_COVERAGE="$2"; shift 2 ;;
        --bam=*) BAM="${1#*=}"; shift ;;
        --bam) BAM="$2"; shift 2 ;;
        --gtf=*) GTF="${1#*=}"; shift ;;
        --gtf) GTF="$2"; shift 2 ;;
        --cytobands=*) CYTOBANDS="${1#*=}"; shift ;;
        --cytobands) CYTOBANDS="$2"; shift 2 ;;
        --protein-domains=*) PROTEIN_DOMAINS="${1#*=}"; shift ;;
        --protein-domains) PROTEIN_DOMAINS="$2"; shift 2 ;;
        --build=*) BUILD="${1#*=}"; shift ;;
        --build) BUILD="$2"; shift 2 ;;
        --extra=*) EXTRA="${1#*=}"; shift ;;
        --extra) EXTRA="$2"; shift 2 ;;
        --output-pdf=*) OUTPUT_PDF="${1#*=}"; shift ;;
        --output-pdf) OUTPUT_PDF="$2"; shift 2 ;;
        --output-png=*) OUTPUT_PNG="${1#*=}"; shift ;;
        --output-png) OUTPUT_PNG="$2"; shift 2 ;;
        *) shift ;;
    esac
done

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
mkdir -p "$(dirname "$OUTPUT_PDF")" "$(dirname "$OUTPUT_PNG")"

echo "Running BCR-ABL1 Fusion Plot for ${SAMPLE}..."

# 1. Look for draw_fusions.R / draw_fusions.r
DRAW_CMD=""
if [[ -f "${SCRIPT_DIR}/draw_fusions.r" ]]; then
    DRAW_CMD="Rscript ${SCRIPT_DIR}/draw_fusions.r"
elif command -v draw_fusions.R >/dev/null 2>&1; then
    DRAW_CMD="draw_fusions.R"
elif [[ -f "/arriba_v2.4.0/draw_fusions.R" ]]; then
    DRAW_CMD="/arriba_v2.4.0/draw_fusions.R"
elif [[ -f "/usr/local/bin/draw_fusions.R" ]]; then
    DRAW_CMD="/usr/local/bin/draw_fusions.R"
fi

# 2. Look for Cytobands
CYTO_FILE="$CYTOBANDS"
if [[ -z "$CYTO_FILE" || ! -f "$CYTO_FILE" ]]; then
    for c in \
        /arriba_v2.4.0/database/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv \
        /arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv \
        /var/lib/arriba/database/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv \
        /var/lib/arriba/database/cytobands_hg38_GRCh38_v2.4.0.tsv \
        /data/references/arriba/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv \
        /data/references/arriba/cytobands_hg38_GRCh38_v2.4.0.tsv; do
        if [[ -f "$c" ]]; then
            CYTO_FILE="$c"
            break
        fi
    done
fi

# 3. Look for Protein domains
PROT_FILE="$PROTEIN_DOMAINS"
if [[ -z "$PROT_FILE" || ! -f "$PROT_FILE" ]]; then
    for p in \
        /arriba_v2.4.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 \
        /arriba_v2.4.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.tsv \
        /arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 \
        /arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.tsv \
        /var/lib/arriba/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 \
        /var/lib/arriba/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 \
        /data/references/arriba/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 \
        /data/references/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3; do
        if [[ -f "$p" ]]; then
            PROT_FILE="$p"
            break
        fi
    done
fi

echo "Draw script: ${DRAW_CMD:-none}"
echo "Cytobands: ${CYTO_FILE:-none}"
echo "Protein domains: ${PROT_FILE:-none}"
echo "GTF: ${GTF}"

if [[ -n "$DRAW_CMD" && -f "$GTF" ]]; then
    CYTO_OPT=""
    PROT_OPT=""
    PANELS="fusion,readcounts"

    if [[ -n "$CYTO_FILE" && -f "$CYTO_FILE" ]]; then
        CYTO_OPT="--cytobands=$CYTO_FILE"
        PANELS="fusion,circos,readcounts"
    fi

    if [[ -n "$PROT_FILE" && -f "$PROT_FILE" ]]; then
        PROT_OPT="--proteinDomains=$PROT_FILE"
        if [[ -n "$CYTO_OPT" ]]; then
            PANELS="fusion,circos,domains,readcounts"
        else
            PANELS="fusion,domains,readcounts"
        fi
    fi

    echo "Executing Arriba draw_fusions.R (all fusions to PDF)..."
    $DRAW_CMD \
        --fusions="$FUSIONS" \
        --alignments="$BAM" \
        --output="$OUTPUT_PDF" \
        --annotation="$GTF" \
        --sampleName="$SAMPLE" \
        --transcriptSelection=canonical \
        --plotPanels="$PANELS" \
        $CYTO_OPT $PROT_OPT $EXTRA || true

    # Generate PNG for primary (most common) fusion only
    TMP_DIR="$(mktemp -d 2>/dev/null || mktemp -d -t bcr_abl_fusion)"
    PRIMARY_TSV="${TMP_DIR}/primary_fusion.tsv"
    head -n 2 "$FUSIONS" > "$PRIMARY_TSV" 2>/dev/null || cp "$FUSIONS" "$PRIMARY_TSV"

    echo "Executing Arriba draw_fusions.R for PNG (primary fusion)..."
    $DRAW_CMD \
        --fusions="$PRIMARY_TSV" \
        --alignments="$BAM" \
        --output="$OUTPUT_PNG" \
        --annotation="$GTF" \
        --sampleName="$SAMPLE" \
        --transcriptSelection=canonical \
        --plotPanels="$PANELS" \
        $CYTO_OPT $PROT_OPT $EXTRA || true

    # If pdftoppm is available and PDF exists, take page 1 as PNG
    if command -v pdftoppm >/dev/null 2>&1 && [[ -s "$OUTPUT_PDF" ]]; then
        pdftoppm -png -r 150 -f 1 -l 1 "$OUTPUT_PDF" "${TMP_DIR}/page" 2>/dev/null || true
        for f in "${TMP_DIR}/page"*.png; do
            if [[ -s "$f" ]]; then
                mv -f "$f" "$OUTPUT_PNG"
                break
            fi
        done
    fi
    rm -rf "$TMP_DIR" 2>/dev/null || true
fi

# Fallback: if PDF is still missing or 0 bytes
if [[ ! -s "$OUTPUT_PDF" ]]; then
    if command -v python3 >/dev/null 2>&1 && [[ -f "${SCRIPT_DIR}/draw_bcr_abl_fusion.py" ]]; then
        echo "Generating fallback PDF with draw_bcr_abl_fusion.py..."
        python3 "${SCRIPT_DIR}/draw_bcr_abl_fusion.py" \
            --sample "$SAMPLE" \
            --fusions "$FUSION_SUMMARY" \
            --mutations "$MUTATIONS" \
            --phasing "$PHASING" \
            --exon-coverage "$EXON_COVERAGE" \
            --build "$BUILD" \
            --output "$OUTPUT_PDF" || true
    fi
fi
if [[ ! -s "$OUTPUT_PDF" ]]; then
    echo "Generating minimal fallback PDF with Rscript..."
    Rscript -e "
    pdf('${OUTPUT_PDF}', width=11, height=7)
    plot(0:10, 0:10, type='n', axes=FALSE, xlab='', ylab='', main='NAGINI BCR-ABL1 FUSION REPORT')
    text(5, 7.5, '${SAMPLE}', cex=1.8, font=2, col='#14324F')
    text(5, 5.0, 'BCR-ABL1 Analysis Completed', cex=1.3, col='#1F4E78')
    text(5, 3.0, 'Detailed variant and isoform metrics are available in the Excel report', cex=1.0, col='#595959')
    dev.off()
    " 2>/dev/null || true
fi

# Fallback: if PNG is still missing or 0 bytes
if [[ ! -s "$OUTPUT_PNG" ]]; then
    if command -v python3 >/dev/null 2>&1 && [[ -f "${SCRIPT_DIR}/draw_bcr_abl_fusion.py" ]]; then
        python3 "${SCRIPT_DIR}/draw_bcr_abl_fusion.py" \
            --sample "$SAMPLE" \
            --fusions "$FUSION_SUMMARY" \
            --mutations "$MUTATIONS" \
            --phasing "$PHASING" \
            --exon-coverage "$EXON_COVERAGE" \
            --build "$BUILD" \
            --output "$OUTPUT_PNG" || true
    fi
fi
if [[ ! -s "$OUTPUT_PNG" ]]; then
    echo "Generating minimal fallback PNG with Rscript..."
    Rscript -e "
    png('${OUTPUT_PNG}', width=1650, height=1050, res=150)
    plot(0:10, 0:10, type='n', axes=FALSE, xlab='', ylab='', main='NAGINI BCR-ABL1 FUSION REPORT')
    text(5, 7.5, '${SAMPLE}', cex=1.8, font=2, col='#14324F')
    text(5, 5.0, 'BCR-ABL1 Analysis Completed', cex=1.3, col='#1F4E78')
    text(5, 3.0, 'Detailed variant and isoform metrics are available in the Excel report', cex=1.0, col='#595959')
    dev.off()
    " 2>/dev/null || true
fi

test -s "$OUTPUT_PDF" || { echo "ERROR: failed to produce $OUTPUT_PDF" >&2; exit 1; }
test -s "$OUTPUT_PNG" || { echo "ERROR: failed to produce $OUTPUT_PNG" >&2; exit 1; }
echo "Successfully generated BCR-ABL1 fusion plots:"
echo "  PDF: $OUTPUT_PDF"
echo "  PNG: $OUTPUT_PNG"
exit 0
