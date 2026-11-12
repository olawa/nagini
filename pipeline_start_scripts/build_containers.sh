#!/usr/bin/env bash
# Build the container images nagini needs, as listed in containers/manifest.yaml.
#
#   build_containers.sh --dest <CONTAINER_DIR> [--only <sif> ...] [--dry-run]
#
# An installation step: it needs network and apptainer or singularity. A run
# pulls nothing (D10); the launcher stops if an image it needs is missing.
# Images that already exist in --dest are left alone. --only limits the build
# to the named .sif files, for one missing image instead of all of them.
# Images that cannot be pulled (hand-built, such as MELT) are listed with
# their source so they can be copied there by hand.
set -euo pipefail

pipeline="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
dest=""
dry_run=false
only=()

while [[ $# -gt 0 ]]; do
    case "$1" in
        --dest) dest="$2"; shift 2 ;;
        --only) only+=("$2"); shift 2 ;;
        --dry-run) dry_run=true; shift ;;
        -h|--help) sed -n '2,12p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

[[ -n "$dest" ]] || { echo "--dest (CONTAINER_DIR) is required" >&2; exit 1; }

manifest="${pipeline}/containers/manifest.yaml"
[[ -f "$manifest" ]] || { echo "Missing $manifest" >&2; exit 1; }

runtime=""
for candidate in apptainer singularity; do
    command -v "$candidate" >/dev/null 2>&1 && runtime="$candidate" && break
done
if [[ -z "$runtime" && "$dry_run" == false ]]; then
    echo "Neither apptainer nor singularity was found" >&2
    exit 1
fi

mkdir -p "$dest"
built=0; skipped=0; manual=0

# Names in --only that the manifest does not know are a typo, not a no-op.
if [[ ${#only[@]} -gt 0 ]]; then
    unknown="$(python3 - "$manifest" "${only[@]}" <<'PY'
import sys, yaml
known = {entry["sif"] for entry in yaml.safe_load(open(sys.argv[1]))["containers"].values()}
print(" ".join(name for name in sys.argv[2:] if name not in known))
PY
)"
    [[ -z "$unknown" ]] || { echo "Not in the manifest: $unknown" >&2; exit 1; }
fi

while IFS=$'\t' read -r sif source; do
    if [[ ${#only[@]} -gt 0 ]]; then
        wanted=false
        for name in "${only[@]}"; do [[ "$name" == "$sif" ]] && wanted=true; done
        $wanted || continue
    fi
    target="${dest}/${sif}"
    if [[ -f "$target" ]]; then
        skipped=$((skipped + 1))
        continue
    fi
    if [[ "$source" != docker://* ]]; then
        echo "MANUAL: ${sif} is not built here, copy it from ${source}"
        manual=$((manual + 1))
        continue
    fi
    echo "building ${sif}  <-  ${source}"
    if [[ "$dry_run" == false ]]; then
        "$runtime" build --force "$target" "$source"
    fi
    built=$((built + 1))
done < <(python3 - "$manifest" <<'PY'
import sys, yaml
for entry in yaml.safe_load(open(sys.argv[1]))["containers"].values():
    print(f"{entry['sif']}\t{entry['source']}")
PY
)

echo
echo "built: ${built}, already there: ${skipped}, need a manual copy: ${manual}"
[[ $manual -eq 0 ]] || exit 2
