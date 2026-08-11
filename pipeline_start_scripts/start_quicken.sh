#!/usr/bin/env bash
set -euo pipefail

# Start Quicken from an ONT inbox while keeping Hydra-style samples.tsv and
# units.tsv in the analysis directory. BAM input is inspected by Hydra's
# create-input-files command; FASTQ input is converted to the same table
# contract by prepare_quicken_inputs.py.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PIPELINE_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"

input_dir=""
analysis_dir="${PWD}"
pipeline_path="${PIPELINE_ROOT}"
base_config="${PIPELINE_ROOT}/config/config.yaml"
profile="${PIPELINE_ROOT}/profiles/local"
input_type=""
sample_type="N"
sequence_id=""
sample_regex=""
barcode="NB"
run_id=""
basecalling_model="unknown"
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
  --inbox-path PATH             Directory containing ONT FASTQ or BAM files.

Options:
  --analysis-dir PATH           Analysis/output directory (default: current directory).
  --pipeline-path PATH          Quicken checkout (default: this checkout).
  --configfile PATH             Base config (default: config/config.yaml).
  --profile PATH                Snakemake profile (default: profiles/local).
  --input-type TYPE             fastq, unmapped_bam, or aligned_bam.
  --sample-type TYPE            N, T, or R (default: N).
  --sequenceid ID               Run identifier; defaults to inbox directory name.
  --sample-regex REGEX          Regex for FASTQ sample extraction (named group: sample).
  --default-barcode BARCODE     Barcode for FASTQ and BAM rows (default: NB).
  --run-id ID                   Run ID for generated FASTQ rows.
  --basecalling-model MODEL     Basecalling model for generated FASTQ rows.
  --methylation Yes|No          Whether modified-base tags are expected (default: No).
  --tumor-content FLOAT         Optional value added to samples.tsv.
  --cores N                     Override the profile's core count.
  --force-inputs                Recreate samples.tsv and units.tsv.
  --dry-run                     Build the DAG without executing jobs.
  -h, --help                    Show this help.

Examples:
  $(basename "$0") --inbox-path /data/run_001 --input-type unmapped_bam
  $(basename "$0") --inbox-path /data/fastq --input-type fastq --sample-regex '^(?P<sample>[^_]+)_'
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
base_config="$(absolute_path "$base_config")"
profile="$(absolute_path "$profile")"

[[ -d "$input_dir" ]] || die "Input directory does not exist: $input_dir"
[[ -d "$pipeline_path" ]] || die "Pipeline directory does not exist: $pipeline_path"
[[ -f "$base_config" ]] || die "Config file does not exist: $base_config"
[[ -f "$pipeline_path/workflow/Snakefile" ]] || die "Snakefile not found below $pipeline_path"
[[ -d "$profile" ]] || die "Snakemake profile does not exist: $profile"
[[ -d "$pipeline_path/vendor/hydra-genetics" ]] || die "Hydra modules are missing. Run $pipeline_path/build/vendor_hydra_dependencies.sh first."
[[ -d "$pipeline_path/vendor/snakemake-wrappers" ]] || die "Snakemake wrappers are missing. Run $pipeline_path/build/vendor_hydra_dependencies.sh first."

require_command python3
require_command snakemake

if [[ -z "$input_type" ]]; then
    bam_count="$(find "$input_dir" -type f -name '*.bam' | wc -l | tr -d ' ')"
    fastq_count="$(find "$input_dir" -type f \( -name '*.fastq' -o -name '*.fastq.gz' -o -name '*.fq' -o -name '*.fq.gz' \) | wc -l | tr -d ' ')"
    if [[ "$bam_count" -gt 0 && "$fastq_count" -eq 0 ]]; then
        input_type="unmapped_bam"
    elif [[ "$fastq_count" -gt 0 && "$bam_count" -eq 0 ]]; then
        input_type="fastq"
    else
        die "Could not infer input type; pass --input-type fastq, unmapped_bam, or aligned_bam."
    fi
fi
case "$input_type" in
    fastq|unmapped_bam|aligned_bam) ;;
    *) die "Unsupported --input-type: $input_type" ;;
esac
if [[ "$input_type" == "unmapped_bam" || "$input_type" == "aligned_bam" ]]; then
    require_command hydra-genetics
fi
[[ "$sample_type" =~ ^(N|T|R)$ ]] || die "--sample-type must be N, T, or R"
[[ "$methylation" =~ ^(Yes|No)$ ]] || die "--methylation must be Yes or No"

mkdir -p "$analysis_dir"
if [[ -z "$sequence_id" ]]; then
    sequence_id="$(basename "$input_dir")"
fi
echo "Input directory: $input_dir"
echo "Analysis directory: $analysis_dir"
echo "Input type: $input_type"
echo "Sequence ID: $sequence_id"

prepare_cmd=(
    python3 "$pipeline_path/pipeline_start_scripts/prepare_quicken_inputs.py"
    --input-dir "$input_dir"
    --analysis-dir "$analysis_dir"
    --pipeline-root "$pipeline_path"
    --base-config "$base_config"
    --input-type "$input_type"
    --sample-type "$sample_type"
    --sequence-id "$sequence_id"
    --barcode "$barcode"
    --run-id "$run_id"
    --basecalling-model "$basecalling_model"
    --methylation "$methylation"
)
[[ -n "$sample_regex" ]] && prepare_cmd+=(--sample-regex "$sample_regex")
[[ -n "$tumor_content" ]] && prepare_cmd+=(--tumor-content "$tumor_content")
$force_inputs && prepare_cmd+=(--force)
"${prepare_cmd[@]}"

cat > "$analysis_dir/run_metadata.yaml" <<EOF
sequence_id: "$sequence_id"
input_dir: "$input_dir"
input_type: "$input_type"
sample_type: "$sample_type"
pipeline_path: "$pipeline_path"
base_config: "$base_config"
EOF

snakemake_cmd=(
    snakemake
    --profile "$profile"
    --snakefile "$pipeline_path/workflow/Snakefile"
    --configfile "$analysis_dir/config.quicken.yaml"
    --wrapper-prefix "git+file://$pipeline_path/vendor/snakemake-wrappers"
)
[[ -n "$cores" ]] && snakemake_cmd+=(--cores "$cores")
$dry_run && snakemake_cmd+=(--dry-run)

cd "$analysis_dir"
"${snakemake_cmd[@]}"
