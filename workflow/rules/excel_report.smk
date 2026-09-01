__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"


rule nagini_excel_report:
    """Build a portable multi-sheet workbook from final QC and variant outputs."""
    input:
        snv_vcf=get_report_snv_vcf,
        snv_tbi=get_report_snv_tbi,
        sv=get_report_sv_tsv,
        translocations=get_report_translocations_tsv,
        wgs_coverage="qc/mosdepth/{sample}_{type}.mosdepth.summary.txt",
        target_coverage=get_report_target_coverage,
        target_dist=get_report_target_dist,
        demux_summary=get_report_demux_summary,
        cramino="qc/cramino/{sample}_{type}.txt",
        bcr_abl_fusion=get_report_bcr_abl_fusion,
        bcr_abl_mutations=get_report_bcr_abl_mutations,
        bcr_abl_phasing=get_report_bcr_abl_phasing,
        bcr_abl_exon_coverage=get_report_bcr_abl_exon_coverage,
    output:
        xlsx="reports/{sample}_{type}.nagini_report.xlsx",
    params:
        reference_build=lambda wildcards: config.get("reference", {}).get("build", config.get("reference_build", "hg38")),
        snv_annotation=(
            "VEP" if vep_enabled() else "none"
        ),
        cytoband=get_cytoband_file,
        small_variant_caller=get_small_variant_caller(),
        design_bed=lambda wildcards: config.get("reference", {}).get("design_bed", ""),
        clairs_to_platform=lambda wildcards: config.get("clairs_to_call", {}).get("platform", ""),
        barcode=get_sample_barcode_name,
    log:
        "logs/nagini_excel_report/{sample}_{type}.log",
    benchmark:
        repeat(
            "benchmarks/nagini_excel_report/{sample}_{type}.tsv",
            config.get("nagini_excel_report", {}).get("benchmark_repeats", 1),
        )
    threads: 1
    resources:
        mem_mb=config.get("nagini_excel_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("nagini_excel_report", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("nagini_excel_report", {}).get("partition", config["default_resources"]["partition"]),
        threads=1,
        time=config.get("nagini_excel_report", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("nagini_excel_report", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: compile Excel report for {wildcards.sample}_{wildcards.type}"
    script:
        "../scripts/nagini_report.py"


def get_summary_report_inputs(wildcards):
    res = {
        "reports": [],
        "craminos": [],
        "wgs_coverages": [],
        "snv_vcfs": [],
        "bcr_abl_fusions": [],
        "bcr_abl_mutations": [],
        "bcr_abl_phasings": [],
        "bcr_abl_exons": [],
    }
    for _, row in samples.iterrows():
        s = str(row["sample"]).strip()
        t = str(row["type"]).strip()
        res["reports"].append(f"reports/{s}_{t}.nagini_report.xlsx")
        res["craminos"].append(f"qc/cramino/{s}_{t}.txt")
        res["wgs_coverages"].append(f"qc/mosdepth/{s}_{t}.mosdepth.summary.txt")

        if sample_needs_dna_variant_calling(s):
            vcf_suffix = ".filtered.norm.vep.vcf.gz" if vep_enabled() else ".filtered.norm.vcf.gz"
            res["snv_vcfs"].append(f"snv_indels/clairs_to_postprocess/{s}_{t}{vcf_suffix}")

        if sample_has_assay_module(s, "bcr_abl"):
            res["bcr_abl_fusions"].append(f"reports/bcr_abl/{s}_{t}.fusion_summary.tsv")
            res["bcr_abl_mutations"].append(f"reports/bcr_abl/{s}_{t}.branford_mutations.tsv")
            res["bcr_abl_phasings"].append(f"reports/bcr_abl/{s}_{t}.phasing.tsv")
            res["bcr_abl_exons"].append(f"reports/bcr_abl/{s}_{t}.exon_coverage.tsv")
    return res


rule nagini_summary_excel_report:
    """Build a multi-sample summary Excel workbook across all samples in the run."""
    input:
        unpack(get_summary_report_inputs),
    output:
        xlsx="reports/nagini_summary_report.xlsx",
    params:
        sample_names=lambda wildcards: [str(row["sample"]).strip() for _, row in samples.iterrows()],
        sample_types=lambda wildcards: [str(row["type"]).strip() for _, row in samples.iterrows()],
        barcodes=lambda wildcards: [get_sample_barcode_name_str(str(row["sample"]).strip(), str(row["type"]).strip()) for _, row in samples.iterrows()],
        assays=lambda wildcards: [get_assay(str(row["sample"]).strip()) for _, row in samples.iterrows()],
        reference_build=lambda wildcards: config.get("reference", {}).get("build", config.get("reference_build", "hg38")),
        sequence_id=lambda wildcards: config.get("sequence_id", os.path.basename(os.getcwd())),
    log:
        "logs/nagini_summary_excel_report/summary.log",
    threads: 1
    resources:
        mem_mb=config.get("nagini_excel_report", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        threads=1,
    container:
        config.get("nagini_excel_report", {}).get("container", config["default_container"])
    conda:
        "../../envs/nagini.yaml"
    message:
        "{rule}: compile multi-sample Excel summary report"
    script:
        "../scripts/nagini_summary_report.py"

