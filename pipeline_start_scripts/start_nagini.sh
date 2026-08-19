#!/usr/bin/env bash
set -euo pipefail

# Start Nagini from an ONT inbox while keeping Hydra-style samples.tsv and
# units.tsv in the analysis directory. BAM input is inspected by Hydra's
# create-input-files command; FASTQ input is converted to the same table
# contract by prepare_nagini_inputs.py. Dorado bam_pass/barcodeXX input can
# use --barcode-map to create the tables without calling Hydra's discovery.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

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
sequence_id=""
sample_regex=""
barcode_map=""
barcode="NB"
barcode_kit=""
run_id=""
basecalling_model="unknown"
dorado_device=""
dorado_models_directory=""
dorado_container=""
methylation="No"
tumor_content=""
cores=""
force_inputs=false
dry_run=false

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
  --sequenceid ID               Run identifier; defaults to inbox directory name.
  --sample-regex REGEX          Regex for FASTQ sample extraction (named group: sample).
  --barcode-map PATH            CSV/TSV barcode-to-sample map for barcode-demultiplexed input.
  --barcode-kit KIT             Dorado barcode kit used for demultiplexing minknow input.
  --default-barcode BARCODE     Barcode for FASTQ and BAM rows (default: NB).
  --run-id ID                   Run ID for generated FASTQ rows.
  --basecalling-model MODEL     Basecalling model (also used by Dorado for minknow input).
  --dorado-model MODEL          Alias for --basecalling-model.
  --dorado-device DEVICE        Dorado device, for example cuda:all or cuda:0.
  --dorado-models-directory DIR Directory containing downloaded Dorado models.
  --dorado-container IMAGE      Override default Dorado Docker/Singularity image.
  --methylation Yes|No          Whether modified-base tags are expected (default: No).
  --tumor-content FLOAT         Optional value added to samples.tsv.
  --cores N                     Override the profile's core count.
  --force-inputs                Recreate samples.tsv and units.tsv.
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
        --cores)
            [[ $# -ge 2 ]] || die "$1 requires a value"
            cores="$2"
            shift 2
            ;;
        --force-inputs)
            force_inputs=true
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
require_command snakemake

if [[ -z "$input_type" ]]; then
    bam_count="$(find "$input_dir" -type f -name '*.bam' | wc -l | tr -d ' ')"
    fastq_count="$(find "$input_dir" -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) | wc -l | tr -d ' ')"
    pod5_count="$(find "$input_dir" -type f -iname '*.pod5' | wc -l | tr -d ' ')"
    if [[ "$pod5_count" -gt 0 ]]; then
        input_type="minknow"
    elif [[ "$bam_count" -gt 0 && "$fastq_count" -eq 0 ]]; then
        input_type="unmapped_bam"
    elif [[ "$fastq_count" -gt 0 && "$bam_count" -eq 0 ]]; then
        input_type="fastq"
    else
        die "Could not infer input type; pass --input-type minknow, fastq, unmapped_bam, or aligned_bam."
    fi
fi
case "$input_type" in
    minknow|fastq|unmapped_bam|aligned_bam) ;;
    *) die "Unsupported --input-type: $input_type" ;;
esac
if [[ ("$input_type" == "unmapped_bam" || "$input_type" == "aligned_bam") && -z "$barcode_map" ]]; then
    require_command hydra-genetics
fi
[[ "$sample_type" =~ ^(N|T|R)$ ]] || die "--sample-type must be N, T, or R"
[[ "$methylation" =~ ^(Yes|No)$ ]] || die "--methylation must be Yes or No"

source_input_dir="$input_dir"
source_input_type="$input_type"
if [[ "$input_type" == "minknow" ]]; then
    [[ -n "$barcode_kit" ]] || die "--barcode-kit is required with --input-type minknow"
    [[ -n "$barcode_map" ]] || die "--barcode-map is required with --input-type minknow"
    [[ "$basecalling_model" != "unknown" ]] || die "--basecalling-model (or --dorado-model) is required with --input-type minknow"
    [[ -z "$dorado_models_directory" || -d "$dorado_models_directory" ]] || die "Dorado models directory does not exist: $dorado_models_directory"
fi

mkdir -p "$analysis_dir"
if [[ -z "$sequence_id" ]]; then
    sequence_id="$(basename "$input_dir")"
fi
echo "Input directory: $input_dir"
echo "Analysis directory: $analysis_dir"
echo "Input type: $input_type"
echo "Sequence ID: $sequence_id"
echo "Site: $site"
echo "Site config: $site_config"
echo "Reference registry: $reference_files"
[[ -n "$barcode_map" ]] && echo "Barcode map: $barcode_map"

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
    --sequence-id "$sequence_id"
    --barcode "$barcode"
    --run-id "$run_id"
    --basecalling-model "$basecalling_model"
    --methylation "$methylation"
)
[[ -n "$sample_regex" ]] && prepare_cmd+=(--sample-regex "$sample_regex")
[[ -n "$barcode_map" ]] && prepare_cmd+=(--barcode-map "$barcode_map")
[[ "$input_type" == "minknow" ]] && prepare_cmd+=(--dorado-input "$input_dir" --dorado-model "$basecalling_model" --dorado-kit "$barcode_kit")
[[ -n "$dorado_device" ]] && prepare_cmd+=(--dorado-device "$dorado_device")
[[ -n "$dorado_models_directory" ]] && prepare_cmd+=(--dorado-models-directory "$dorado_models_directory")
[[ -n "$dorado_container" ]] && prepare_cmd+=(--dorado-container "$dorado_container")
[[ -n "$tumor_content" ]] && prepare_cmd+=(--tumor-content "$tumor_content")
[[ -n "$library_root" ]] && prepare_cmd+=(--library-root "$library_root")
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
site: "$site"
pipeline_path: "$pipeline_path"
base_config: "$base_config"
site_config: "$site_config"
reference_files: "$reference_files"
barcode_map: "${barcode_map:-}"
barcode_kit: "${barcode_kit:-}"
basecalling_model: "$basecalling_model"
dorado_device: "${dorado_device:-}"
dorado_models_directory: "${dorado_models_directory:-}"
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

cd "$analysis_dir"
"${snakemake_cmd[@]}"
