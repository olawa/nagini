__author__ = "Ola Wallerman"
__copyright__ = "Copyright 2026, Ola Wallerman"
__email__ = "ola.wallerman@scilifelab.uu.se"
__license__ = "GPL-3"

"""Pipeline-wide configuration, input handling and output compilation.

The unit table deliberately contains an ``input_type`` column.  This keeps
the three input paths explicit:

* ``fastq`` -> minimap2 alignment
* ``unmapped_bam`` -> samtools fastq, then minimap2 alignment
* ``aligned_bam`` -> copied and indexed, never realigned

Hydra modules are imported after this file, so the selected reference is
available to every module as ``config["reference"]``.
"""

import json
import os
import pathlib
import re
import shlex
import sys
from datetime import datetime

import pandas as pd
import yaml
from snakemake.exceptions import WorkflowError
from snakemake.utils import min_version, validate

from hydra_genetics import min_version as hydra_min_version
from hydra_genetics.utils.misc import export_config_as_file
from hydra_genetics.utils.misc import get_module_snakefile
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import get_samples
from hydra_genetics.utils.software_versions import add_software_version_to_config
from hydra_genetics.utils.software_versions import add_version_files_to_multiqc
from hydra_genetics.utils.software_versions import export_pipeline_version_as_file
from hydra_genetics.utils.software_versions import export_software_version_as_file
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name
from hydra_genetics.utils.software_versions import touch_software_version_file
from hydra_genetics.utils.software_versions import use_container


hydra_min_version("3.0.0")
min_version("7.8.0")


if not workflow.overwrite_configfiles:
    sys.exit("Pass a config file with --configfile or --configfiles.")


_config_template = re.compile(r"\{\{\s*([A-Za-z_][A-Za-z0-9_]*)\s*\}\}")


def _expand_config_templates(value, context):
    if isinstance(value, dict):
        return {key: _expand_config_templates(item, context) for key, item in value.items()}
    if isinstance(value, list):
        return [_expand_config_templates(item, context) for item in value]
    if not isinstance(value, str):
        return value
    expanded = value
    for _ in range(12):
        updated = _config_template.sub(
            lambda match: context.get(match.group(1), match.group(0)), expanded
        )
        if updated == expanded:
            break
        expanded = updated
    return expanded


_config_context = {
    str(key): str(value)
    for key, value in config.items()
    if isinstance(value, (str, int, float, bool))
}
for _ in range(12):
    _changed = False
    for _key, _value in list(_config_context.items()):
        _expanded = _expand_config_templates(_value, _config_context)
        if _expanded != _value:
            _config_context[_key] = _expanded
            _changed = True
    if not _changed:
        break
config = _expand_config_templates(config, _config_context)
_unresolved_config_variables = sorted(
    {
        name
        for match in _config_template.finditer(str(config))
        for name in [match.group(1)]
    }
)
if _unresolved_config_variables:
    raise WorkflowError(
        "Unresolved config variable(s): " + ", ".join(_unresolved_config_variables)
    )

try:
    validate(config, schema="../schemas/config.schema.yaml")
except WorkflowError as error:
    if not error.args[0].lower().startswith("error validating config file"):
        raise
    error_msg = "\n".join(error.args[0].splitlines()[:2])
    parent_schema_path = error.args[0].splitlines()[3].split()[-1]
    if parent_schema_path == "schema:":
        sys.exit(error_msg)
    schema_hierarchy = parent_schema_path.split()[-1]
    schema_section = ".".join(re.findall(r"\['([^']+)'\]", schema_hierarchy)[1::2])
    sys.exit(f"{error_msg} in {schema_section}")


reference_build = config.get("reference_build", "hg38")
references = config.get("references", {})
if reference_build not in references:
    raise WorkflowError(
        f"reference_build={reference_build} is not defined in config['references']. "
        f"Available builds: {', '.join(sorted(references))}"
    )
config["reference"] = references[reference_build]
config["reference"]["build"] = reference_build

def vep_enabled():
    """Whether the report should depend on the VEP-annotated VCF."""

    return bool(
        config.get("annotate_small_variants_with_vep", False)
        and not config.get("vep", {}).get("skip", False)
    )


if vep_enabled() and not config.get("vep", {}).get("vep_cache"):
    raise WorkflowError(
        "annotate_small_variants_with_vep is true, but vep.vep_cache is empty. "
        "Configure a reference-compatible offline VEP cache or disable VEP annotation."
    )


date_string = datetime.now().strftime("%Y%m%d--%H-%M-%S")
pipeline_version = get_pipeline_version(workflow, pipeline_name="quicken")
version_files = touch_pipeline_version_file_name(
    pipeline_version,
    date_string=date_string,
    directory="results/versions/software",
)
update_config = config
if use_container(workflow):
    version_files.append(
        touch_software_version_file(
            config,
            date_string=date_string,
            directory="results/versions/software",
        )
    )
add_version_files_to_multiqc(config, version_files)


onstart:
    export_pipeline_version_as_file(
        pipeline_version,
        date_string=date_string,
        directory="results/versions/software",
    )
    if use_container(workflow):
        update_config, software_info = add_software_version_to_config(config, workflow, False)
        export_software_version_as_file(
            software_info,
            date_string=date_string,
            directory="results/versions/software",
        )
    export_config_as_file(update_config, date_string=date_string, directory="results/versions")


config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")

samples = pd.read_table(config["samples"], dtype=str).fillna("").set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

# Keep this dataframe flat.  Hydra's ONT modules index the same table in
# their own namespace, while the custom input rules need input_type as well.
units = pd.read_table(config["units"], dtype=str).fillna("")
validate(units, schema="../schemas/units.schema.yaml")

with open(config["output"]) as output_file:
    if config["output"].endswith("json"):
        output_spec = json.load(output_file)
    else:
        output_spec = yaml.safe_load(output_file.read())
validate(output_spec, schema="../schemas/output_files.schema.yaml")


wildcard_constraints:
    sample="|".join(re.escape(str(sample)) for sample in samples.index),
    type="N|T|R",
    processing_unit="[A-Za-z0-9_.-]+",
    barcode="[A-Za-z0-9+_.-]+",


def _unit_rows(sample, unit_type, processing_unit=None, barcode=None):
    """Return unit rows for a sample/type and optional ONT unit key."""

    rows = units[(units["sample"] == sample) & (units["type"] == unit_type)]
    if processing_unit is not None:
        rows = rows[rows["processing_unit"] == processing_unit]
    if barcode is not None:
        rows = rows[rows["barcode"] == barcode]
    return rows


def get_unit(wildcards):
    rows = _unit_rows(
        wildcards.sample,
        wildcards.type,
        getattr(wildcards, "processing_unit", None),
        getattr(wildcards, "barcode", None),
    )
    if len(rows) != 1:
        raise WorkflowError(
            "Expected exactly one units.tsv row for "
            f"sample={wildcards.sample}, type={wildcards.type}, "
            f"processing_unit={getattr(wildcards, 'processing_unit', None)}, "
            f"barcode={getattr(wildcards, 'barcode', None)}; found {len(rows)}."
        )
    return rows.iloc[0]


def get_unit_input(wildcards, column):
    row = get_unit(wildcards)
    value = str(row.get(column, "")).strip()
    if value in {"", ".", "NA", "nan", "None"}:
        raise WorkflowError(
            f"units.tsv column {column} is empty for "
            f"{wildcards.sample}/{wildcards.type}/{wildcards.processing_unit}/{wildcards.barcode}."
        )
    return value


def get_unit_input_type(wildcards):
    return str(get_unit(wildcards)["input_type"]).strip().lower()


def get_unit_alignment_path(row):
    suffix = {
        "fastq": "fastq",
        "unmapped_bam": "unmapped_bam",
        "aligned_bam": "aligned_bam",
    }.get(str(row["input_type"]).strip().lower())
    if suffix is None:
        raise WorkflowError(
            f"Unsupported input_type={row['input_type']}; use fastq, unmapped_bam or aligned_bam."
        )
    return "alignment/quicken_align/{}_{}_{}_{}.{}.bam".format(
        row["sample"],
        row["type"],
        row["processing_unit"],
        row["barcode"],
        suffix,
    )


def get_unit_alignment_bams(wildcards):
    rows = _unit_rows(wildcards.sample, wildcards.type)
    if rows.empty:
        raise WorkflowError(f"No units.tsv rows for {wildcards.sample}/{wildcards.type}.")
    return [get_unit_alignment_path(row) for _, row in rows.iterrows()]


def get_reference_index():
    configured_index = config["reference"].get("mmi", "")
    if configured_index and os.path.isfile(configured_index):
        # Reuse a site-managed index instead of rebuilding a multi-gigabyte
        # minimap2 index in every analysis directory.
        return configured_index
    fasta_name = os.path.basename(config["reference"]["fasta"])
    preset = config.get("minimap2_align", {}).get("preset", "map-ont")
    return "alignment/minimap2_index/{}.{}.mmi".format(fasta_name, preset)


def _optional_reference_path(key):
    """Return an optional selected-reference asset only when it is readable."""

    value = str(config.get("reference", {}).get(key, "") or "").strip()
    return value if value and os.path.isfile(value) else ""


def get_clinvar_vcf(wildcards=None):
    return _optional_reference_path("clinvar")


def get_cytoband_file(wildcards=None):
    return _optional_reference_path("cytoband")


def get_tandem_repeats(wildcards=None):
    return _optional_reference_path("tandem_repeats")


# Hydra's cnv_sv module reads this setting when its rule is imported. Keep the
# raw reference path in the config and let the module add the CLI flag itself.
_tandem_repeats = get_tandem_repeats()
if _tandem_repeats:
    config.setdefault("sniffles2_call", {})["tandem_repeats"] = _tandem_repeats


def get_minimap2_read_group(wildcards):
    configured = config.get("minimap2_align", {}).get("read_group", "")
    if configured:
        return configured.format(
            sample=wildcards.sample,
            type=wildcards.type,
            processing_unit=wildcards.processing_unit,
            barcode=wildcards.barcode,
        )
    return (
        r"@RG\tID:{sample}_{type}_{processing_unit}_{barcode}"
        r"\tSM:{sample}_{type}\tPL:ONT"
    ).format(
        sample=wildcards.sample,
        type=wildcards.type,
        processing_unit=wildcards.processing_unit,
        barcode=wildcards.barcode,
    )


def get_reference_chromosomes():
    chromosomes = config["reference"].get("chromosomes", [])
    if chromosomes:
        return list(chromosomes)
    fai = config["reference"].get("fai", f"{config['reference']['fasta']}.fai")
    if not os.path.exists(fai):
        raise WorkflowError(
            "No reference.chromosomes configured and the reference .fai is not available: "
            f"{fai}"
        )
    with open(fai) as fai_file:
        return [line.split("\t", 1)[0] for line in fai_file if line.strip()]


reference_chromosomes = get_reference_chromosomes()


def get_small_variant_report_vcf(wildcards):
    """Return the normalized VCF, preferring VEP annotation when enabled."""
    suffix = ".filtered.norm.vep.vcf.gz" if vep_enabled() else ".filtered.norm.vcf.gz"
    return f"snv_indels/mutect2_postprocess/{wildcards.sample}_{wildcards.type}{suffix}"


def get_small_variant_report_vcf_index(wildcards):
    return f"{get_small_variant_report_vcf(wildcards)}.tbi"


def get_unit_types_for_sample(sample):
    return sorted(set(units.loc[units["sample"] == sample, "type"]))


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec.get("directory", "./"))
    output_files = []
    for file_definition in output_spec["files"]:
        condition = file_definition.get("when")
        enabled = vep_enabled() if condition == "annotate_small_variants_with_vep" else config.get(condition, False)
        if condition and not enabled:
            continue
        allowed_types = file_definition.get("types")
        for sample in get_samples(samples):
            for unit_type in get_unit_types_for_sample(sample):
                if allowed_types and unit_type not in allowed_types:
                    continue
                output_files.append(
                    outdir
                    / file_definition["output"].format(
                        sample=sample,
                        type=unit_type,
                    )
                )
    return output_files
