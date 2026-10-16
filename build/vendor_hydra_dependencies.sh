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
  "alignment:v0.9.0" \
  "annotation:v1.5.0" \
  "cnv_sv:v3.4.0" \
  "qc:v0.8.0" \
  "reports:v1.1.1" \
  "snv_indels:v2.1.0"; do
  module="${module_tag%%:*}"
  tag="${module_tag##*:}"
  destination="${module_root}/${module}"
  if [[ -d "${destination}/.git" ]]; then
    echo "Keeping existing ${destination}; checking for ${tag}"
    # Modules are read from their git tag (LocalGitFile), so a pin bump only
    # needs the new tag present, not a checkout of it.
    if ! git -C "${destination}" rev-parse --verify "refs/tags/${tag}" >/dev/null 2>&1; then
      echo "Fetching ${tag} into ${destination}"
      git -C "${destination}" fetch --depth 1 origin "refs/tags/${tag}:refs/tags/${tag}" >/dev/null
    fi
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
if git -C "${snv_module}" rev-parse --verify refs/tags/v2.1.0-nagini >/dev/null 2>&1; then
  echo "Keeping existing snv_indels compatibility tag v2.1.0-nagini"
elif git -C "${snv_module}" rev-parse --verify refs/tags/v2.1.0-quicken >/dev/null 2>&1; then
  # Migrate a vendor checkout created before the pipeline was renamed.
  git -C "${snv_module}" tag -f "v2.1.0-nagini" "refs/tags/v2.1.0-quicken" >/dev/null
  echo "Migrated existing snv_indels tag v2.1.0-quicken to v2.1.0-nagini"
elif git -C "${snv_module}" apply --check "${snv_patch}" >/dev/null 2>&1; then
  git -C "${snv_module}" apply "${snv_patch}"
  git -C "${snv_module}" -c user.name="nagini builder" -c user.email="nagini@localhost" commit -am "Compatibility patch for Snakemake 7" >/dev/null
  git -C "${snv_module}" tag -f "v2.1.0-nagini" >/dev/null
else
  echo "snv_indels is neither pristine v2.1.0 nor tagged v2.1.0-nagini; refusing to overwrite it" >&2
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
# Snakemake wrapper rules reference version tags (for example v1.5.0 for
# GATK Mutect2).  Fetch those tags explicitly; fetching only the pinned
# commit leaves the wrapper files present but makes git show <tag>:<path>
# fail during execution.
for wrapper_tag in 0.80.2 0.85.0 v1.1.0 v1.3.1 v1.5.0 v1.7.0 v1.25.0 v3.9.0 v4.3.0; do
  git -C "${wrapper_root}" fetch --depth 1 origin \
    "refs/tags/${wrapper_tag}:refs/tags/${wrapper_tag}" >/dev/null
done
echo "Using snakemake-wrappers@${wrapper_commit}"

echo "Hydra modules and wrappers are available under vendor/."
