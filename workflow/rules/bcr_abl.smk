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
        if command -v draw_fusions.R >/dev/null 2>&1; then
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
        if [ -n "$DRAW_CMD" ] && [ -f "{params.gtf}" ] && [ -n "$CYTO_FILE" ] && [ -f "$CYTO_FILE" ]; then
            PROT_OPT=""
            PANELS="fusion,circos,readcounts"
            if [ -n "$PROT_FILE" ] && [ -f "$PROT_FILE" ]; then
                PROT_OPT="--proteinDomains=$PROT_FILE"
                PANELS="fusion,circos,domains,readcounts"
            fi

            echo "Executing Arriba draw_fusions.R..."
            $DRAW_CMD \
                --fusions={input.fusions} \
                --alignments={input.bam} \
                --output={output.pdf} \
                --annotation={params.gtf} \
                --cytobands="$CYTO_FILE" \
                --transcriptSelection=canonical \
                --plotPanels=$PANELS \
                $PROT_OPT {params.extra} || true
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

        test -s "{output.pdf}"
        """




