__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


def _bcr_abl_resource(filename):
    """Return absolute path to a bundled BCR-ABL1 resource file."""
    return os.path.join(workflow.basedir, "resources", "bcr_abl", filename)


def _branford_bed():
    """Pick the branford positions BED that matches the configured reference build."""
    override = config.get("bcr_abl", {}).get("branford_bed", "")
    if override:
        return override
    build = str(config.get("reference_build", "hg38")).strip().lower()
    bed_name = "branford_positions_hg19.bed" if "19" in build else "branford_positions_hg38.bed"
    return _bcr_abl_resource(bed_name)


def _bcr_abl_exons_bed():
    """Pick the BCR-ABL1 exon coordinates BED that matches the configured reference build."""
    override = config.get("bcr_abl", {}).get("exons_bed", "")
    if override:
        return override
    build = str(config.get("reference_build", "hg38")).strip().lower()
    bed_name = "bcr_abl1_exons_hg19.bed" if "19" in build else "bcr_abl1_exons_hg38.bed"
    return _bcr_abl_resource(bed_name)


def _bcr_abl_gtf():
    """Pick the gene annotation GTF for BCR-ABL1 fusion plotting."""
    override = config.get("bcr_abl", {}).get("gtf", "")
    if override and os.path.isfile(override):
        return override
    ref_gtf = config.get("reference", {}).get("gene_annotation", "") or config.get("reference", {}).get("gtf", "")
    if ref_gtf and os.path.isfile(ref_gtf):
        return ref_gtf
    build = str(config.get("reference_build", "hg38")).strip().lower()
    gtf_name = "bcr_abl_hg19.gtf" if "19" in build else "bcr_abl_hg38.gtf"
    return _bcr_abl_resource(gtf_name)


def _find_rindels_binary():
    """Locate the Rindels binary from config, PATH, ~/bin, or dev checkout."""
    override = config.get("rindels", {}).get("binary", "") or config.get("bcr_abl", {}).get("rindels_binary", "")
    if override and (os.path.isfile(override) or shutil.which(override)):
        return override
    candidates = [
        os.path.expanduser("~/bin/rindels"),
        "/home/cgu-lab/bin/rindels",
        os.path.expanduser("~/tmp/tools/rindels/target/release/rindels"),
        "/home/cgu-lab/tmp/tools/rindels/target/release/rindels",
        os.path.expanduser("~/.cargo/bin/rindels"),
    ]

    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    which_r = shutil.which("rindels")
    if which_r:
        return which_r
    return "rindels"




rule bcr_abl_fusion_counts:
    """Identify and quantify BCR-ABL1 fusion isoforms from spliced alignments."""
    input:
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
        features=lambda wildcards: (
            config.get("bcr_abl", {}).get("features_json", "")
            or _bcr_abl_resource("bcr_abl_features.json")
        ),
    output:
        tsv="reports/bcr_abl/{sample}_{type}.fusion_summary.tsv",
        json="reports/bcr_abl/{sample}_{type}.fusion_summary.json",
        arriba="reports/bcr_abl/{sample}_{type}.arriba_fusions.tsv",
    params:
        tolerance=lambda wildcards: config.get("bcr_abl", {}).get("tolerance", 150),
        min_mapq=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 1),
        genome_build=lambda wildcards: config.get("reference_build", "hg38"),
    log:
        "logs/bcr_abl_fusion_counts/{sample}_{type}.log",
    threads: config.get("bcr_abl", {}).get("threads", 4)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 4),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        python3 {workflow.basedir}/scripts/bcr_abl_fusion_counter.py \
            --bam {input.bam} \
            --sample {wildcards.sample} \
            --genome-build {params.genome_build} \
            --features-json {input.features} \
            --output {output.tsv} \
            --json-output {output.json} \
            --arriba-output {output.arriba} \
            --tolerance {params.tolerance} \
            --min-mapq {params.min_mapq} \
            --threads {threads} > {log} 2>&1
        """


rule bcr_abl_branford_calling:
    """Force-call BCR-ABL1 kinase domain resistance mutations at high depth."""
    input:
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
        bed=lambda wildcards: _branford_bed(),
        features=lambda wildcards: (
            config.get("bcr_abl", {}).get("features_json", "")
            or _bcr_abl_resource("bcr_abl_features.json")
        ),
    output:
        all_sites="reports/bcr_abl/{sample}_{type}.branford_all_sites.tsv",
        mutations="reports/bcr_abl/{sample}_{type}.branford_mutations.tsv",
        vcf="reports/bcr_abl/{sample}_{type}.branford_mutations.vcf",
    params:
        min_af=lambda wildcards: get_assay_setting(wildcards.sample, "bcr_abl", "min_af", config.get("bcr_abl", {}).get("min_af", 0.01)),
        min_depth=lambda wildcards: config.get("bcr_abl", {}).get("min_depth", 50),
        min_alt_reads=lambda wildcards: config.get("bcr_abl", {}).get("min_alt_reads", 5),
        max_depth=lambda wildcards: config.get("bcr_abl", {}).get("max_depth", 1000000),
    log:
        "logs/bcr_abl_branford_calling/{sample}_{type}.log",
    threads: config.get("bcr_abl", {}).get("threads", 4)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 4),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        python3 {workflow.basedir}/scripts/bcr_abl_branford_caller.py \
            --bam {input.bam} \
            --bed {input.bed} \
            --sample {wildcards.sample} \
            --output-all {output.all_sites} \
            --output-mutations {output.mutations} \
            --output-vcf {output.vcf} \
            --features-json {input.features} \
            --min-af {params.min_af} \
            --min-depth {params.min_depth} \
            --min-alt-reads {params.min_alt_reads} \
            --max-depth {params.max_depth} \
            --threads {threads} > {log} 2>&1
        """


rule bcr_abl_phasing:
    """Single-molecule phasing of compound resistance mutations from long reads."""
    input:
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
        mutations="reports/bcr_abl/{sample}_{type}.branford_mutations.tsv",
    output:
        tsv="reports/bcr_abl/{sample}_{type}.phasing.tsv",
        json="reports/bcr_abl/{sample}_{type}.phasing.json",
    params:
        min_spanning_reads=lambda wildcards: config.get("bcr_abl", {}).get("min_spanning_reads", 10),
        min_mapq=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 1),
    log:
        "logs/bcr_abl_phasing/{sample}_{type}.log",
    threads: config.get("bcr_abl", {}).get("threads", 4)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 4),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        python3 {workflow.basedir}/scripts/bcr_abl_phaser.py \
            --bam {input.bam} \
            --mutations-tsv {input.mutations} \
            --sample {wildcards.sample} \
            --output-tsv {output.tsv} \
            --output-json {output.json} \
            --min-spanning-reads {params.min_spanning_reads} \
            --min-mapq {params.min_mapq} \
            --threads {threads} > {log} 2>&1
        """


rule bcr_abl_exon_coverage:
    """Calculate depth and coverage metrics across BCR and ABL1 exons."""
    input:
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
        bed=lambda wildcards: _bcr_abl_exons_bed(),
    output:
        tsv="reports/bcr_abl/{sample}_{type}.exon_coverage.tsv",
    params:
        min_mapq=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 1),
        max_depth=lambda wildcards: config.get("bcr_abl", {}).get("max_depth", 1000000),
    log:
        "logs/bcr_abl_exon_coverage/{sample}_{type}.log",
    threads: config.get("bcr_abl", {}).get("threads", 4)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 4),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        python3 {workflow.basedir}/scripts/bcr_abl_exon_coverage.py \
            --bam {input.bam} \
            --bed {input.bed} \
            --sample {wildcards.sample} \
            --output {output.tsv} \
            --min-mapq {params.min_mapq} \
            --max-depth {params.max_depth} \
            --threads {threads} > {log} 2>&1
        """



rule bcr_abl_draw_fusions:
    """Generate Arriba fusion diagrams using draw_fusions.R."""
    input:
        fusions="reports/bcr_abl/{sample}_{type}.arriba_fusions.tsv",
        fusion_summary="reports/bcr_abl/{sample}_{type}.fusion_summary.tsv",
        mutations="reports/bcr_abl/{sample}_{type}.branford_mutations.tsv",
        phasing="reports/bcr_abl/{sample}_{type}.phasing.tsv",
        exon_coverage="reports/bcr_abl/{sample}_{type}.exon_coverage.tsv",
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
    output:
        pdf="reports/bcr_abl/{sample}_{type}.fusion_plot.pdf",
        png="reports/bcr_abl/{sample}_{type}.fusion_plot.png",
    params:
        gtf=lambda wildcards: _bcr_abl_gtf(),
        cytobands=lambda wildcards: config.get("reference", {}).get("cytoband", "") or config.get("reference", {}).get("cytobands", ""),
        protein_domains=lambda wildcards: config.get("arriba", {}).get("protein_domains", ""),
        build=lambda wildcards: str(config.get("reference_build", "hg19")).strip().lower(),
        extra=config.get("arriba", {}).get("extra", ""),
    log:
        "logs/bcr_abl_draw_fusions/{sample}_{type}.log",
    threads: config.get("arriba", {}).get("threads", 1)
    resources:
        mem_mb=config.get("arriba", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("arriba", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("arriba", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("arriba", {}).get("threads", 1),
        time=config.get("arriba", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("arriba", {}).get("container", config.get("ARRIBA_CONTAINER", "docker://uhrigs/arriba:2.4.0"))
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        exec > {log} 2>&1

        echo "Running BCR-ABL1 Fusion Plot for {wildcards.sample}..."

        # 1. Look for draw_fusions.R
        DRAW_CMD=""
        if [ -f "{workflow.basedir}/scripts/draw_fusions.r" ]; then
            DRAW_CMD="Rscript {workflow.basedir}/scripts/draw_fusions.r"
        elif command -v draw_fusions.R >/dev/null 2>&1; then
            DRAW_CMD="draw_fusions.R"
        elif [ -f "/arriba_v2.4.0/draw_fusions.R" ]; then
            DRAW_CMD="/arriba_v2.4.0/draw_fusions.R"
        elif [ -f "/usr/local/bin/draw_fusions.R" ]; then
            DRAW_CMD="/usr/local/bin/draw_fusions.R"
        fi

        # 2. Look for Cytobands
        CYTO_FILE="{params.cytobands}"
        if [ -z "$CYTO_FILE" ] || [ ! -f "$CYTO_FILE" ]; then
            for c in \
                /arriba_v2.4.0/database/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv \
                /arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv \
                /var/lib/arriba/database/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv \
                /var/lib/arriba/database/cytobands_hg38_GRCh38_v2.4.0.tsv \
                /data/references/arriba/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv \
                /data/references/arriba/cytobands_hg38_GRCh38_v2.4.0.tsv
            do
                if [ -f "$c" ]; then
                    CYTO_FILE="$c"
                    break
                fi
            done
        fi

        # 3. Look for Protein domains
        PROT_FILE="{params.protein_domains}"
        if [ -z "$PROT_FILE" ] || [ ! -f "$PROT_FILE" ]; then
            for p in \
                /arriba_v2.4.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 \
                /arriba_v2.4.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.tsv \
                /arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 \
                /arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.tsv \
                /var/lib/arriba/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 \
                /var/lib/arriba/database/protein_domains_hg38_GRCh38_v2.4.0.gff3 \
                /data/references/arriba/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 \
                /data/references/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3
            do
                if [ -f "$p" ]; then
                    PROT_FILE="$p"
                    break
                fi
            done
        fi

        echo "Draw script: $DRAW_CMD"
        echo "Cytobands: $CYTO_FILE"
        echo "Protein domains: $PROT_FILE"
        echo "GTF: {params.gtf}"

        # Try Arriba draw_fusions.R / bundled draw_fusions.r if available
        if [ -n "$DRAW_CMD" ] && [ -f "{params.gtf}" ]; then
            PROT_OPT=""
            CYTO_OPT=""
            PANELS="fusion,readcounts"

            if [ -n "$CYTO_FILE" ] && [ -f "$CYTO_FILE" ]; then
                CYTO_OPT="--cytobands=$CYTO_FILE"
                PANELS="fusion,circos,readcounts"
            fi

            if [ -n "$PROT_FILE" ] && [ -f "$PROT_FILE" ]; then
                PROT_OPT="--proteinDomains=$PROT_FILE"
                if [ -n "$CYTO_OPT" ]; then
                    PANELS="fusion,circos,domains,readcounts"
                else
                    PANELS="fusion,domains,readcounts"
                fi
            fi

            echo "Executing Arriba draw_fusions.R (all fusions to PDF)..."
            $DRAW_CMD \
                --fusions={input.fusions} \
                --alignments={input.bam} \
                --output={output.pdf} \
                --annotation={params.gtf} \
                --sampleName="{wildcards.sample}" \
                --transcriptSelection=canonical \
                --plotPanels=$PANELS \
                $CYTO_OPT $PROT_OPT {params.extra} || true

            # Generate PNG for primary (most common) fusion only
            PRIMARY_TSV="reports/bcr_abl/{wildcards.sample}_{wildcards.type}.primary_fusion.tsv"
            head -n 2 {input.fusions} > "$PRIMARY_TSV" 2>/dev/null || cp {input.fusions} "$PRIMARY_TSV"

            echo "Executing Arriba draw_fusions.R for PNG (primary fusion)..."
            $DRAW_CMD \
                --fusions="$PRIMARY_TSV" \
                --alignments={input.bam} \
                --output={output.png} \
                --annotation={params.gtf} \
                --sampleName="{wildcards.sample}" \
                --transcriptSelection=canonical \
                --plotPanels=$PANELS \
                $CYTO_OPT $PROT_OPT {params.extra} || true

            PDF_FILE="{output.pdf}"
            PNG_FILE="{output.png}"
            P1_PREFIX="reports/bcr_abl/{wildcards.sample}_{wildcards.type}.fusion_plot_p1"

            # If pdftoppm is available, ensure first page of PDF is used as PNG
            if command -v pdftoppm >/dev/null 2>&1 && [ -s "$PDF_FILE" ]; then
                pdftoppm -png -r 150 -f 1 -l 1 "$PDF_FILE" "$P1_PREFIX" 2>/dev/null || true
                for f in "$P1_PREFIX-1.png" "$P1_PREFIX-01.png" "$P1_PREFIX-001.png"; do
                    if [ -s "$f" ]; then
                        mv -f "$f" "$PNG_FILE"
                        break
                    fi
                done
                rm -f "$P1_PREFIX"*.png 2>/dev/null || true
            fi
        fi




        # If PDF was not produced or is 0 bytes, generate publication-ready fusion diagram via python
        if [ ! -s "{output.pdf}" ]; then
            echo "Generating BCR-ABL1 diagram with draw_bcr_abl_fusion.py..."
            python3 {workflow.basedir}/scripts/draw_bcr_abl_fusion.py \
                --sample {wildcards.sample} \
                --fusions {input.fusion_summary} \
                --mutations {input.mutations} \
                --phasing {input.phasing} \
                --exon-coverage {input.exon_coverage} \
                --build {params.build} \
                --output {output.pdf}
        fi

        # If PNG was not produced or is 0 bytes, generate via python
        if [ ! -s "{output.png}" ]; then
            python3 {workflow.basedir}/scripts/draw_bcr_abl_fusion.py \
                --sample {wildcards.sample} \
                --fusions {input.fusion_summary} \
                --mutations {input.mutations} \
                --phasing {input.phasing} \
                --exon-coverage {input.exon_coverage} \
                --build {params.build} \
                --output {output.png} || true
        fi

        test -s "{output.pdf}"
        test -s "{output.png}"
        """


rule bcr_abl_bcftools_pileup:
    """Run bcftools mpileup and call across BCR-ABL1 target exons at full depth."""
    input:
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
        ref=config["reference"]["fasta"],
        bed=lambda wildcards: _bcr_abl_exons_bed(),
    output:
        vcf="snv_indels/bcftools_mpileup/{sample}_{type}.bcr_abl.vcf.gz",
        tbi="snv_indels/bcftools_mpileup/{sample}_{type}.bcr_abl.vcf.gz.tbi",
    params:
        max_depth=lambda wildcards: config.get("bcr_abl", {}).get("max_depth", 100000),
        min_mapq=lambda wildcards: config.get("bcr_abl", {}).get("min_mapq", 20),
        min_baseq=lambda wildcards: config.get("bcr_abl", {}).get("min_base_quality", 12),
        prior=lambda wildcards: config.get("bcr_abl", {}).get("pileup_prior", 0.01),
        extra=lambda wildcards: config.get("bcr_abl", {}).get("bcftools_mpileup_extra", ""),
    log:
        "logs/bcr_abl_bcftools_pileup/{sample}_{type}.log",
    threads: config.get("bcr_abl", {}).get("threads", 2)
    resources:
        mem_mb=config.get("bcr_abl", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcr_abl", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcr_abl", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcr_abl", {}).get("threads", 2),
        time=config.get("bcr_abl", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools", {}).get("container", "docker://hydragenetics/common:3.1.1")
    conda:
        "../../envs/nagini.yaml"
    shell:
        """
        exec > {log} 2>&1
        bcftools mpileup -Ou -f {input.ref} -R {input.bed} \
            -d {params.max_depth} -q {params.min_mapq} -Q {params.min_baseq} \
            -a FORMAT/AD,FORMAT/DP {params.extra} {input.bam} | \
        bcftools call -mv -P {params.prior} -Oz -o {output.vcf}

        bcftools index -t {output.vcf}
        """


rule bcr_abl_rindels:
    """Run Rindels in targeted ONT BCR-ABL1 amplicon mode."""
    input:
        bam="alignment/nagini_bam/{sample}_{type}.bam",
        bai="alignment/nagini_bam/{sample}_{type}.bam.bai",
        ref=config["reference"]["fasta"],
    output:
        vcf="snv_indels/rindels/{sample}_{type}.rindels.vcf.gz",
        tbi="snv_indels/rindels/{sample}_{type}.rindels.vcf.gz.tbi",
        haplotypes="snv_indels/rindels/{sample}_{type}.haplotypes.tsv",
    params:
        binary=lambda wildcards: _find_rindels_binary(),
        genome=lambda wildcards: "hg19" if "19" in str(config.get("reference_build", "hg19")).lower() else "hg38",
        extra=lambda wildcards: config.get("rindels", {}).get("extra", "--ont --amplicon --bcr-abl"),
    log:
        "logs/bcr_abl_rindels/{sample}_{type}.log",
    threads: config.get("rindels", {}).get("threads", 4)
    resources:
        mem_mb=config.get("rindels", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("rindels", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("rindels", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("rindels", {}).get("threads", 4),
        time=config.get("rindels", {}).get("time", config["default_resources"]["time"]),
    shell:
        """
        exec > {log} 2>&1
        RINDELS_BIN="{params.binary}"
        if ! command -v "$RINDELS_BIN" >/dev/null 2>&1 && [ ! -x "$RINDELS_BIN" ]; then
            echo "Error: Rindels binary not found at $RINDELS_BIN"
            exit 1
        fi

        "$RINDELS_BIN" -b {input.bam} -r {input.ref} -o {output.vcf} \
            {params.extra} --bcr-abl-genome {params.genome} \
            --write-haplotypes {output.haplotypes} -t {threads}

        if [ -s "{output.vcf}" ] && [ ! -s "{output.tbi}" ]; then
            if command -v bcftools >/dev/null 2>&1; then
                bcftools index -t {output.vcf} || true
            fi
        fi
        """






