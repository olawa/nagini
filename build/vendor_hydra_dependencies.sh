#!/usr/bin/env bash
set -euo pipefail

# One-time online preparation.  After this script has completed, the module
# Snakefiles and Snakemake wrappers are local and the workflow can be parsed
# without GitHub access.  Reference data and container images are prepared by
# separate Hydra commands documented in README.md.

pipeline_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
module_root="${pipeline_root}/vendor/hydra-genetics"
wrapper_root="${pipeline_root}/vendor/snakemake-wrappers"
wrapper_commit="d3e6a345acb4e03ac0662ade94d9970bb233bdab"

mkdir -p "${module_root}"

for module_tag in \
  "alignment:v0.8.0" \
  "annotation:v1.4.1" \
  "cnv_sv:v3.2.1" \
  "qc:v0.8.0" \
  "snv_indels:v2.1.0"; do
  module="${module_tag%%:*}"
  tag="${module_tag##*:}"
  destination="${module_root}/${module}"
  if [[ -d "${destination}/.git" ]]; then
    echo "Keeping existing ${destination}; checking for ${tag}"
    git -C "${destination}" rev-parse --verify "refs/tags/${tag}" >/dev/null
  else
    echo "Cloning hydra-genetics/${module}@${tag}"
    git clone --branch "${tag}" --depth 1 "https://github.com/hydra-genetics/${module}.git" "${destination}"
  fi
done

# snv_indels v2.1.0 contains an unused DeepTrio rule whose formatted
# ``expand`` expression is not parseable by the Snakemake 7.32 runtime used
# for this pipeline. Keep the Mutect2/BCFtools rules from the pinned module,
# but remove that unused include in a local compatibility commit.
snv_module="${module_root}/snv_indels"
snv_patch="${pipeline_root}/build/patches/snv_indels_v2.1.0_snakemake7.patch"
if git -C "${snv_module}" rev-parse --verify refs/tags/v2.1.0-quicken >/dev/null 2>&1; then
  echo "Keeping existing snv_indels compatibility tag v2.1.0-quicken"
elif git -C "${snv_module}" apply --check "${snv_patch}" >/dev/null 2>&1; then
  git -C "${snv_module}" apply "${snv_patch}"
  git -C "${snv_module}" -c user.name="quicken builder" -c user.email="quicken@localhost" commit -am "Compatibility patch for Snakemake 7" >/dev/null
  git -C "${snv_module}" tag -f "v2.1.0-quicken" >/dev/null
else
  echo "snv_indels is neither pristine v2.1.0 nor tagged v2.1.0-quicken; refusing to overwrite it" >&2
  exit 1
fi

if [[ -d "${wrapper_root}/.git" ]]; then
  echo "Keeping existing ${wrapper_root}"
else
  echo "Fetching pinned snakemake-wrappers commit"
  mkdir -p "${wrapper_root}"
  git -C "${wrapper_root}" init -q
  git -C "${wrapper_root}" remote add origin "https://github.com/snakemake/snakemake-wrappers.git"
  git -C "${wrapper_root}" fetch --depth 1 origin "${wrapper_commit}" >/dev/null
fi
git -C "${wrapper_root}" checkout --detach "${wrapper_commit}" >/dev/null
echo "Using snakemake-wrappers@${wrapper_commit}"

echo "Hydra modules and wrappers are available under vendor/."
