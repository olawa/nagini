#!/usr/bin/env bash
set -euo pipefail

# Start Nagini from an ONT inbox while keeping Hydra-style samples.tsv and
# units.tsv in the analysis directory. BAM input is inspected by Hydra's
# create-input-files command; FASTQ input is converted to the same table
# contract by prepare_nagini_inputs.py. Dorado bam_pass/barcodeXX input can
# use --barcode-map to create the tables without calling Hydra's discovery.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

# Activate conda environment 'nagini' if available and not already active
if [[ "${CONDA_DEFAULT_ENV:-}" != "nagini" ]]; then
    if command -v conda &>/dev/null; then
        eval "$(conda shell.bash hook 2>/dev/null || true)"
        conda activate nagini 2>/dev/null || true
    elif [[ -f "${HOME}/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "${HOME}/miniconda3/etc/profile.d/conda.sh"
        conda activate nagini 2>/dev/null || true
    elif [[ -f "${HOME}/anaconda3/etc/profile.d/conda.sh" ]]; then
        source "${HOME}/anaconda3/etc/profile.d/conda.sh"
        conda activate nagini 2>/dev/null || true
    elif [[ -f "/opt/conda/etc/profile.d/conda.sh" ]]; then
        source "/opt/conda/etc/profile.d/conda.sh"
        conda activate nagini 2>/dev/null || true
    elif [[ -f "/home/cgu-lab/miniconda3/etc/profile.d/conda.sh" ]]; then
        source "/home/cgu-lab/miniconda3/etc/profile.d/conda.sh"
        conda activate nagini 2>/dev/null || true
    fi
fi


input_dir=""
analysis_dir="${PWD}"
pipeline_path="${PIPELINE_ROOT}"
base_config=""
profile=""
site="local"
site_config=""
reference_files=""
library_root=""
input_type=""
sample_type="N"
library_type=""
assay=""
reference_build=""
sequence_id=""
sample_regex=""
barcode_map=""
init_sample_sheet=""
barcode="NB"
barcode_kit=""
run_id=""
basecalling_model="dna_r10.4.1_e8.2_400bps_sup@v5.0.0"
fastq_status="all"

dorado_device=""
dorado_models_directory=""
dorado_container=""
dorado_binary=""
dorado_demux_output_format="fastq"
dorado_demux_both_ends=false
dorado_demux_extra=""
dorado_aligner="minimap2"

bed_file=""
small_variant_caller=""
clairs_to_platform=""
clairs_to_container=""
methylation="No"
tumor_content=""
cores=""
filter_reads=""
merged_fastq=""
barcode_map_origin=""
filter_qv=""
filter_len=""
filter_binary=""
force_inputs=false
dry_run=false
unlock=false



usage() {
    cat <<EOF
Usage:
  $(basename "$0") --inbox-path PATH [options]

Required:
  --inbox-path PATH             Directory containing ONT POD5, FASTQ, or BAM files.

Options:
  --analysis-dir PATH           Analysis/output directory (default: current directory).
  --pipeline-path PATH          Nagini checkout (default: this checkout).
  --configfile PATH             Base config (default: config/config.yaml).
  --site NAME                   Runtime site: local, workstation, or marvin (default: local).
  --site-config PATH            Site config (default: config/site_configs/site_config_<site>.yaml).
  --reference-files PATH        Reference registry (default: config/reference_files/reference_files_<site>.yaml).
  --library-root PATH           Optional cgu/library root recorded in the run snapshot.
  --profile PATH                Snakemake profile (default: profiles/<site>).
  --input-type TYPE             minknow, fastq, unmapped_bam, or aligned_bam.
  --sample-type TYPE            N, T, or R (default: N).
  --library-type TYPE           Default library type: wgs, amplicon, cdna, or direct_rna
                                (default: wgs). Decides how reads are aligned.
  --assay NAME                  Default assay: wgs, adaptive, bcr_abl, or one declared in
                                config assays. Decides which analysis runs. A per-sample
                                assay column in the barcode map overrides it.
  --reference-build BUILD       Override config['reference_build'] for this run, e.g. hg19
                                for --assay bcr_abl (default: whatever config.yaml says,
                                normally hg38). One run uses one build for every sample.
  --sequenceid ID               Run identifier; defaults to inbox directory name.
  --sample-regex REGEX          Regex for FASTQ sample extraction (named group: sample).
  --barcode-map PATH            CSV/TSV barcode-to-sample map for barcode-demultiplexed FASTQ/BAM.
                                sample_sheet*.csv/tsv is auto-detected for existing FASTQ when omitted.
  --init-sample-sheet [FILE]    Generate a sample sheet CSV template in the input directory and exit.
  --barcode-kit KIT             Dorado barcode kit used for demultiplexing minknow input.
  --default-barcode BARCODE     Barcode for FASTQ and BAM rows (default: NB).
  --run-id ID                   Run ID for generated FASTQ rows.
  --basecalling-model MODEL     Basecalling model (also used by Dorado for minknow input).
  --dorado-model MODEL          Alias for --basecalling-model.
  --fastq-status all|pass       Existing FASTQ input: include pass and fail reads (default: all),
                                or only fastq_pass reads.
  --dorado-device DEVICE        Dorado device, for example cuda:all or cuda:0.
  --dorado-models-directory DIR Directory containing downloaded Dorado models.
  --dorado-container IMAGE      Override default Dorado Docker/Singularity image (or 'none' for local binary).
  --dorado-binary PATH          Path to local dorado binary (e.g. ~/tmp/tools/dorado-2.1.1-linux-x64/bin/dorado).
  --demux-output-format FORMAT  Per-barcode raw output: bam, fastq, or both (default: fastq).
  --demux-both-ends             Require barcode detection at both ends in Dorado demux.
  --demux-extra FLAGS           Additional flags passed to dorado demux.
  --dorado-aligner ALIGNER      Alignment after demux: minimap2 or dorado (default: minimap2).

  --bed-file PATH               Custom target region BED file (aliases: --bed, --design-bed).
                                (If omitted, any single *.bed in inbox is auto-detected).
  --small-variant-caller CALLER clairs_to, mutect2, or both (default: clairs_to).
  --clairs-to-platform MODEL    ClairS-TO platform model (auto-selected from basecalling model if omitted).
  --clairs-to-container IMAGE   Override default ClairS-TO container image.
  --methylation Yes|No          Whether modified-base tags are expected (default: No).
  --merged-fastq                Also write one FASTQ per sample (merged across flowcells)
                                to reads/nagini_fastq/. Unaligned BAM input only.
  --tumor-content FLOAT         Optional value added to samples.tsv.
  --filter-reads                Enable pre-alignment read filtering (nanofilter).
  --no-filter-reads             Disable pre-alignment read filtering.
  --filter-qv QV                Minimum mean read quality (e.g. 20).
  --filter-len MIN-MAX          Read length range (e.g. 500-3000).
  --filter-binary PATH          Path to nanofilter binary or script.
  --cores N                     Override the profile's core count.
  --force-inputs                Recreate samples.tsv and units.tsv.

  --unlock                      Remove a stale Snakemake directory lock before/instead of running.
  --dry-run                     Build the DAG without executing jobs.
  -h, --help                    Show this help.


Examples:
  $(basename "$0") --inbox-path /data/run_001 --input-type unmapped_bam
  $(basename "$0") --inbox-path /data/fastq --input-type fastq --sample-regex '^(?P<sample>[^_]+)_'
  $(basename "$0") --inbox-path /data/minknow_run --input-type minknow \\
    --barcode-kit SQK-RBK114-24 --basecalling-model sup --barcode-map barcode_map.csv
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_command() {
    command -v "$1" >/dev/null 2>&1 || die "Required command not found: $1"
}

absolute_path() {
    local path="$1"
    if [[ "$path" == /* ]]; then
        printf '%s\n' "$path"
    else
        printf '%s/%s\n' "$PWD" "$path"
    fi
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --inbox-path|--input-dir)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            input_dir="$2"
            shift 2
            ;;
        --analysis-dir|--work-dir)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            analysis_dir="$2"
            shift 2
            ;;
        --pipeline-path|--repo-path)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            pipeline_path="$2"
            shift 2
            ;;
        --configfile)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            base_config="$2"
            shift 2
            ;;
        --site)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            site="$2"
            shift 2
            ;;
        --site-config)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            site_config="$2"
            shift 2
            ;;
        --reference-files)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            reference_files="$2"
            shift 2
            ;;
        --library-root)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            library_root="$2"
            shift 2
            ;;
        --profile)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            profile="$2"
            shift 2
            ;;
        --input-type)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            input_type="$2"
            shift 2
            ;;
        --sample-type)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            sample_type="$2"
            shift 2
            ;;
        --assay)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            assay="$2"
            shift 2
            ;;
        --reference-build)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            reference_build="$2"
            shift 2
            ;;
        --library-type)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            library_type="$2"
            shift 2
            ;;
        --sequenceid|--sequence-id)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            sequence_id="$2"
            shift 2
            ;;
        --sample-regex)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            sample_regex="$2"
            shift 2
            ;;
        --barcode-map)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            barcode_map="$2"
            shift 2
            ;;
        --init-sample-sheet|--create-sample-sheet)
            if [[ $# -ge 2 && ! "$2" =~ ^-- ]]; then
                init_sample_sheet="$2"
                shift 2
            else
                init_sample_sheet="sample_sheet.csv"
                shift 1
            fi
            ;;
        --barcode-kit)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            barcode_kit="$2"
            shift 2
            ;;
        --default-barcode|--barcode)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            barcode="$2"
            shift 2
            ;;
        --run-id)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            run_id="$2"
            shift 2
            ;;
        --basecalling-model)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            basecalling_model="$2"
            shift 2
            ;;
        --dorado-model)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            basecalling_model="$2"
            shift 2
            ;;
        --fastq-status)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            fastq_status="$2"
            shift 2
            ;;
        --dorado-device)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            dorado_device="$2"
            shift 2
            ;;
        --dorado-models-directory)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            dorado_models_directory="$2"
            shift 2
            ;;
        --dorado-container)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            dorado_container="$2"
            shift 2
            ;;
        --dorado-binary|--dorado-path)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            dorado_binary="$2"
            shift 2
            ;;
        --demux-output-format|--dorado-demux-output-format)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            dorado_demux_output_format="$2"
            shift 2
            ;;
        --demux-both-ends|--both-ends)
            dorado_demux_both_ends=true
            shift
            ;;
        --demux-extra)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            dorado_demux_extra="$2"
            shift 2
            ;;
        --dorado-aligner)

            [[ $# -ge 2 ]] || die "$1 requires a value"
            dorado_aligner="$2"
            shift 2
            ;;
        --bed-file|--bed|--design-bed)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            bed_file="$(absolute_path "$2")"
            [[ -f "$bed_file" ]] || die "Target BED file does not exist: $bed_file"
            shift 2
            ;;
        --small-variant-caller)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            small_variant_caller="$2"
            shift 2
            ;;
        --clairs-to-platform)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            clairs_to_platform="$2"
            shift 2
            ;;
        --clairs-to-container)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            clairs_to_container="$2"
            shift 2
            ;;
        --methylation)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            methylation="$2"
            shift 2
            ;;
        --tumor-content)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            tumor_content="$2"
            shift 2
            ;;
        --merged-fastq)
            merged_fastq=true
            shift
            ;;
        --no-merged-fastq)
            merged_fastq=false
            shift
            ;;
        --filter-reads)
            filter_reads=true
            shift
            ;;
        --no-filter-reads)
            filter_reads=false
            shift
            ;;
        --filter-qv)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            filter_qv="$2"
            shift 2
            ;;
        --filter-len)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            filter_len="$2"
            shift 2
            ;;
        --filter-binary)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            filter_binary="$2"
            shift 2
            ;;
        --cores)

            [[ $# -ge 2 ]] || die "$1 requires a value"
            cores="$2"
            shift 2
            ;;
        --force|--force-inputs)
            force_inputs=true
            shift
            ;;
        --unlock)
            unlock=true
            shift
            ;;
        --dry-run)
            dry_run=true
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            usage >&2
            die "Unknown option: $1"
            ;;
    esac
done

if $unlock && [[ -z "$input_dir" ]]; then
    pipeline_path="$(absolute_path "$pipeline_path")"
    analysis_dir="$(absolute_path "$analysis_dir")"
    [[ -n "$profile" ]] || profile="$pipeline_path/profiles/${site}"
    profile="$(absolute_path "$profile")"
    echo "Unlocking analysis directory: $analysis_dir"
    (
        cd "$analysis_dir"
        snakemake \
            --profile "$profile" \
            --snakefile "$pipeline_path/workflow/Snakefile" \
            --unlock
    )
    exit 0
fi

[[ -n "$input_dir" ]] || { usage >&2; die "--inbox-path is required"; }

pipeline_path="$(absolute_path "$pipeline_path")"
input_dir="$(absolute_path "$input_dir")"
analysis_dir="$(absolute_path "$analysis_dir")"
[[ -n "$base_config" ]] || base_config="$pipeline_path/config/config.yaml"
[[ -n "$site_config" ]] || site_config="$pipeline_path/config/site_configs/site_config_${site}.yaml"
[[ -n "$reference_files" ]] || reference_files="$pipeline_path/config/reference_files/reference_files_${site}.yaml"
[[ -n "$profile" ]] || profile="$pipeline_path/profiles/${site}"
base_config="$(absolute_path "$base_config")"
site_config="$(absolute_path "$site_config")"
reference_files="$(absolute_path "$reference_files")"
profile="$(absolute_path "$profile")"
if [[ -n "$library_root" ]]; then
    library_root="$(absolute_path "$library_root")"
fi
if [[ -n "$barcode_map" ]]; then
    barcode_map="$(absolute_path "$barcode_map")"
elif [[ -d "$input_dir" ]]; then
    for candidate in \
        "$input_dir/barcode_map.csv" "$input_dir/barcode_map.tsv" \
        "$input_dir/barcode-map.csv" "$input_dir/barcode-map.tsv" \
        "$input_dir/sample_sheet.csv" "$input_dir/sample_sheet.tsv" \
        "$input_dir/sample-sheet.csv" "$input_dir/sample-sheet.tsv" \
        "$input_dir/samplesheet.csv" "$input_dir/samplesheet.tsv" \
        "$input_dir/sample_map.csv" "$input_dir/sample_map.tsv" \
        "$input_dir/sample-map.csv" "$input_dir/sample-map.tsv" \
        "$input_dir/barcodes.csv" "$input_dir/barcodes.tsv"
    do
        if [[ -f "$candidate" ]]; then
            barcode_map="$(absolute_path "$candidate")"
            barcode_map_origin="auto-detected in inbox"
            break
        fi
    done
fi
if [[ -n "$dorado_models_directory" ]]; then

    dorado_models_directory="$(absolute_path "$dorado_models_directory")"
fi

[[ -d "$input_dir" ]] || die "Input directory does not exist: $input_dir"
[[ -d "$pipeline_path" ]] || die "Pipeline directory does not exist: $pipeline_path"
[[ -f "$base_config" ]] || die "Config file does not exist: $base_config"
[[ -z "$barcode_map" || -f "$barcode_map" ]] || die "Barcode map does not exist: $barcode_map"
case "$site" in
    local|workstation|marvin) ;;
    *) die "Unsupported --site: $site (expected local, workstation, or marvin)" ;;
esac
[[ -f "$site_config" ]] || die "Site config does not exist: $site_config"
[[ -f "$reference_files" ]] || die "Reference registry does not exist: $reference_files"
[[ -f "$pipeline_path/workflow/Snakefile" ]] || die "Snakefile not found below $pipeline_path"
[[ -d "$profile" ]] || die "Snakemake profile does not exist: $profile"
[[ -d "$pipeline_path/vendor/hydra-genetics" ]] || die "Hydra modules are missing. Run $pipeline_path/build/vendor_hydra_dependencies.sh first."
[[ -d "$pipeline_path/vendor/snakemake-wrappers" ]] || die "Snakemake wrappers are missing. Run $pipeline_path/build/vendor_hydra_dependencies.sh first."

require_command python3

if [[ -n "$init_sample_sheet" ]]; then
    init_cmd=(
        python3 "$pipeline_path/pipeline_start_scripts/prepare_nagini_inputs.py"
        --init-sample-sheet "$init_sample_sheet"
        --input-dir "$input_dir"
        --pipeline-root "$pipeline_path"
    )
    [[ -n "$barcode_kit" ]] && init_cmd+=(--barcode-kit "$barcode_kit")
    [[ -n "$sample_type" ]] && init_cmd+=(--sample-type "$sample_type")
    [[ -n "$assay" ]] && init_cmd+=(--assay "$assay")
    [[ -n "$library_type" ]] && init_cmd+=(--library-type "$library_type")
    "${init_cmd[@]}"
    exit 0
fi

require_command snakemake

if [[ -z "$input_type" ]]; then
    bam_count="$(find "$input_dir" -type f -name '*.bam' | wc -l | tr -d ' ')"
    fastq_count="$(find "$input_dir" -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) | wc -l | tr -d ' ')"
    pod5_count="$(find "$input_dir" -type f -iname '*.pod5' | wc -l | tr -d ' ')"
    # A completed Dorado export often contains both POD5 and barcode FASTQ.
    # Prefer the reusable FASTQ unless the caller explicitly selects minknow;
    # this avoids silently spending GPU time on a second basecalling pass.
    if [[ "$fastq_count" -gt 0 ]]; then
        input_type="fastq"
    elif [[ "$pod5_count" -gt 0 ]]; then
        input_type="minknow"
    elif [[ "$bam_count" -gt 0 ]]; then
        input_type="unmapped_bam"
    else
        die "Could not infer input type; pass --input-type minknow, fastq, unmapped_bam, or aligned_bam."
    fi
fi
case "$input_type" in
    minknow|fastq|unmapped_bam|aligned_bam) ;;
    *) die "Unsupported --input-type: $input_type" ;;
esac
[[ "$sample_type" =~ ^(N|T|R)$ ]] || die "--sample-type must be N, T, or R"
[[ "$methylation" =~ ^(Yes|No)$ ]] || die "--methylation must be Yes or No"
case "$dorado_demux_output_format" in
    bam|fastq|both) ;;
    *) die "--demux-output-format must be bam, fastq, or both" ;;
esac
case "$dorado_aligner" in
    minimap2|dorado) ;;
    *) die "--dorado-aligner must be minimap2 or dorado" ;;
esac
case "$fastq_status" in
    all|pass) ;;
    *) die "--fastq-status must be all or pass" ;;
esac
source_input_dir="$input_dir"
source_input_type="$input_type"
if [[ "$input_type" == "minknow" ]]; then
    [[ -n "$barcode_kit" ]] || die "--barcode-kit is required with --input-type minknow"
    if [[ -z "$barcode_map" ]]; then
        echo "NOTE: No --barcode-map provided. Samples will be automatically named per barcode (e.g. bc1, bc10) based on kit $barcode_kit."
    fi
    [[ "$basecalling_model" != "unknown" ]] || die "--basecalling-model (or --dorado-model) is required with --input-type minknow"
    [[ -z "$dorado_models_directory" || -d "$dorado_models_directory" ]] || die "Dorado models directory does not exist: $dorado_models_directory"
fi

mkdir -p "$analysis_dir"
if [[ -z "$sequence_id" ]]; then
    sequence_id="$(basename "$input_dir")"
fi
# A launch banner rather than a flat list of echos: which inbox, which
# analysis directory and which site config were picked is what a rerun is
# reconstructed from, and grouping keeps the three apart at a glance.
section() { printf '\n%s\n' "$1"; }
field() { printf '  %-18s %s\n' "$1" "$2"; }
field_if() { [[ -n "$2" ]] && field "$1" "$2"; return 0; }

section "Input"
field "directory" "$input_dir"
field "type" "$input_type"
field "sequence id" "$sequence_id"
if [[ -n "$barcode_map" ]]; then
    if [[ -n "$barcode_map_origin" ]]; then
        field "barcode map" "$barcode_map ($barcode_map_origin)"
    else
        field "barcode map" "$barcode_map"
    fi
fi
field_if "barcode kit" "$barcode_kit"
field_if "basecalling model" "$basecalling_model"
field_if "methylation" "$methylation"

section "Analysis"
field "directory" "$analysis_dir"
field_if "assay" "${assay:-from barcode map}"
field_if "library type" "$library_type"
field_if "sample type" "$sample_type"
field_if "reference build" "$reference_build"
field "variant caller" "${small_variant_caller:-auto (ClairS-TO + DeepVariant)}"
field_if "design bed" "$bed_file"
field_if "merged fastq" "$merged_fastq"

section "Site"
field "site" "$site"
field "config" "$site_config"
field "reference registry" "$reference_files"
field "profile" "$profile"
field_if "pipeline" "$pipeline_path"
printf '\n'

prepare_cmd=(
    python3 "$pipeline_path/pipeline_start_scripts/prepare_nagini_inputs.py"
    --input-dir "$input_dir"
    --analysis-dir "$analysis_dir"
    --pipeline-root "$pipeline_path"
    --base-config "$base_config"
    --site-config "$site_config"
    --reference-files "$reference_files"
    --input-type "$input_type"
    --sample-type "$sample_type"
    --library-type "$library_type"
    --sequence-id "$sequence_id"
    --barcode "$barcode"
    --run-id "$run_id"
    --basecalling-model "$basecalling_model"
    --fastq-status "$fastq_status"
    --methylation "$methylation"
)
[[ -n "$assay" ]] && prepare_cmd+=(--assay "$assay")
[[ -n "$reference_build" ]] && prepare_cmd+=(--reference-build "$reference_build")
[[ -n "$sample_regex" ]] && prepare_cmd+=(--sample-regex "$sample_regex")
[[ -n "$barcode_map" ]] && prepare_cmd+=(--barcode-map "$barcode_map")
[[ "$input_type" == "minknow" ]] && prepare_cmd+=(--dorado-input "$input_dir" --dorado-model "$basecalling_model" --dorado-kit "$barcode_kit")
[[ "$input_type" == "minknow" ]] && prepare_cmd+=(--dorado-demux-output-format "$dorado_demux_output_format" --dorado-aligner "$dorado_aligner")
[[ "$dorado_demux_both_ends" == "true" ]] && prepare_cmd+=(--dorado-demux-both-ends)
[[ -n "$dorado_demux_extra" ]] && prepare_cmd+=(--dorado-demux-extra "$dorado_demux_extra")
[[ -n "$dorado_device" ]] && prepare_cmd+=(--dorado-device "$dorado_device")

[[ -n "$dorado_models_directory" ]] && prepare_cmd+=(--dorado-models-directory "$dorado_models_directory")
[[ -n "$dorado_container" ]] && prepare_cmd+=(--dorado-container "$dorado_container")
[[ -n "$dorado_binary" ]] && prepare_cmd+=(--dorado-binary "$dorado_binary")
[[ -n "$bed_file" ]] && prepare_cmd+=(--bed-file "$bed_file")
[[ -n "$small_variant_caller" ]] && prepare_cmd+=(--small-variant-caller "$small_variant_caller")
[[ -n "$clairs_to_platform" ]] && prepare_cmd+=(--clairs-to-platform "$clairs_to_platform")
[[ -n "$clairs_to_container" ]] && prepare_cmd+=(--clairs-to-container "$clairs_to_container")
[[ -n "$cores" ]] && prepare_cmd+=(--cores "$cores")
[[ -n "$tumor_content" ]] && prepare_cmd+=(--tumor-content "$tumor_content")
[[ -n "$library_root" ]] && prepare_cmd+=(--library-root "$library_root")
[[ "$merged_fastq" == "true" ]] && prepare_cmd+=(--merged-fastq)
[[ "$merged_fastq" == "false" ]] && prepare_cmd+=(--no-merged-fastq)
[[ "$filter_reads" == "true" ]] && prepare_cmd+=(--filter-reads)
[[ "$filter_reads" == "false" ]] && prepare_cmd+=(--no-filter-reads)
[[ -n "$filter_qv" ]] && prepare_cmd+=(--filter-qv "$filter_qv")
[[ -n "$filter_len" ]] && prepare_cmd+=(--filter-len "$filter_len")
[[ -n "$filter_binary" ]] && prepare_cmd+=(--filter-binary "$filter_binary")
$force_inputs && prepare_cmd+=(--force)
"${prepare_cmd[@]}"


cp "$site_config" "$analysis_dir/site_config.yaml"
cp "$reference_files" "$analysis_dir/reference_files.yaml"

cat > "$analysis_dir/run_metadata.yaml" <<EOF
sequence_id: "$sequence_id"
input_dir: "$input_dir"
input_type: "$input_type"
source_input_dir: "$source_input_dir"
source_input_type: "$source_input_type"
sample_type: "$sample_type"
library_type: "$library_type"
assay: "${assay:-}"
site: "$site"
pipeline_path: "$pipeline_path"
base_config: "$base_config"
site_config: "$site_config"
reference_files: "$reference_files"
barcode_map: "${barcode_map:-}"
barcode_kit: "${barcode_kit:-}"
basecalling_model: "$basecalling_model"
fastq_status: "$fastq_status"
merged_fastq: "${merged_fastq:-}"
dorado_device: "${dorado_device:-}"
dorado_models_directory: "${dorado_models_directory:-}"
dorado_demux_output_format: "$dorado_demux_output_format"
dorado_aligner: "$dorado_aligner"
runtime_config: "$analysis_dir/config.nagini.yaml"
profile: "$profile"
library_root: "${library_root:-}"
EOF

snakemake_cmd=(
    snakemake
    --profile "$profile"
    --snakefile "$pipeline_path/workflow/Snakefile"
    --configfile "$analysis_dir/config.nagini.yaml"
    --wrapper-prefix "git+file://$pipeline_path/vendor/snakemake-wrappers"
)
[[ -n "$cores" ]] && snakemake_cmd+=(--cores "$cores")
$dry_run && snakemake_cmd+=(--dry-run)

if $unlock; then

    echo "Unlocking analysis directory: $analysis_dir"
    (
        cd "$analysis_dir"
        snakemake \
            --profile "$profile" \
            --snakefile "$pipeline_path/workflow/Snakefile" \
            --configfile "$analysis_dir/config.nagini.yaml" \
            --unlock || true
    )
fi

cd "$analysis_dir"
"${snakemake_cmd[@]}"

