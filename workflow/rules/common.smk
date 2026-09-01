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
* ``dorado_bam`` -> produced by the pipeline's Dorado basecalling/demux stage

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

# Routing lives in a plain module so that it can be unit tested without
# Snakemake.  workflow.basedir is workflow/, so scripts/ sits next to rules/.
SCRIPTS = os.path.join(workflow.basedir, "scripts")
sys.path.insert(0, SCRIPTS)
import nagini_routing


hydra_min_version("3.3.0")
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
pipeline_version = get_pipeline_version(workflow, pipeline_name="nagini")
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


def _scale_workstation_resources(cfg):
    """Dynamically scale compute-heavy sequential rules to the allocated core budget."""
    site = str(cfg.get("site", "")).strip().lower()
    # Apply dynamic scaling on workstation and local runs (not cluster scheduler environments like marvin)
    if site not in ("workstation", "local", ""):
        return cfg

    cores = None
    if cfg.get("cores"):
        try:
            cores = int(cfg["cores"])
        except (ValueError, TypeError):
            pass
    if not cores and hasattr(workflow, "cores") and workflow.cores:
        try:
            cores = int(workflow.cores)
        except (ValueError, TypeError):
            pass

    if not cores or cores <= 1:
        return cfg

    # Compute-heavy sequential rules that should utilize the full allocated core budget:
    full_core_rules = (
        "minimap2_fastq_align",
        "minimap2_unmapped_bam_align",
        "dorado_alignment",
        "clairs_to_call",
        "sniffles2_call",
        "vep",
    )
    for rule_name in full_core_rules:
        if rule_name in cfg and isinstance(cfg[rule_name], dict):
            cfg[rule_name]["threads"] = cores

    # Moderate parallel / intermediate rules scaled proportionally:
    moderate_core_rules = {
        "samtools_merge_bam": max(4, min(16, cores // 2)),
        "extract_target_bam": max(4, min(16, cores // 2)),
        "mosdepth": max(4, min(16, cores // 2)),
        "mosdepth_bed": max(4, min(16, cores // 2)),
    }
    for rule_name, target_threads in moderate_core_rules.items():
        if rule_name in cfg and isinstance(cfg[rule_name], dict):
            cfg[rule_name]["threads"] = target_threads

    if "dorado" in cfg and isinstance(cfg["dorado"], dict):
        cfg["dorado"]["demux_threads"] = max(4, min(24, cores // 2))

    return cfg


config = _scale_workstation_resources(config)

samples = (
    pd.read_table(config["samples"], dtype=str)
    .fillna("")
    .apply(lambda col: col.str.strip() if hasattr(col, "str") else col)
    .set_index("sample", drop=False)
)
validate(samples, schema="../schemas/samples.schema.yaml")

# Keep this dataframe flat.  Hydra's ONT modules index the same table in
# their own namespace, while the custom input rules need input_type as well.
units = (
    pd.read_table(config["units"], dtype=str)
    .fillna("")
    .apply(lambda col: col.str.strip() if hasattr(col, "str") else col)
)
validate(units, schema="../schemas/units.schema.yaml")

dorado_units = units[units["input_type"].str.lower() == "dorado_bam"]
dorado_barcodes = sorted(
    {
        str(value).strip()
        for value in dorado_units.get("barcode_name", pd.Series(dtype=str)).tolist()
        if str(value).strip()
    }
)
dorado_demux_output_format = str(
    config.get("dorado", {}).get("demux_output_format", "fastq")
).strip().lower()

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
        "dorado_bam": "dorado_bam",
    }.get(str(row["input_type"]).strip().lower())
    if suffix is None:
        raise WorkflowError(
            f"Unsupported input_type={row['input_type']}; use fastq, unmapped_bam, aligned_bam or dorado_bam."
        )
    return "alignment/nagini_align/{}_{}_{}_{}.{}.bam".format(
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


def get_unit_dorado_bam(wildcards):
    """Return the pipeline-owned demultiplexed BAM for a Dorado unit."""

    row = get_unit(wildcards)
    barcode = str(row.get("barcode_name", "")).strip()
    if not barcode:
        raise WorkflowError("Dorado units require barcode_name in units.tsv.")
    return str(row.get("bam", f"basecalling/dorado_demux/{barcode}.bam"))


def get_unit_dorado_fastq(wildcards):
    """Return the exported FASTQ for a pipeline-produced Dorado unit."""

    row = get_unit(wildcards)
    barcode = str(row.get("barcode_name", "")).strip()
    if not barcode:
        raise WorkflowError("Dorado units require barcode_name in units.tsv.")
    return f"basecalling/dorado_fastq/{barcode}.fastq.gz"


def get_reference_index():
    configured_index = config["reference"].get("mmi", "")
    if configured_index and os.path.isfile(configured_index):
        # Reuse a site-managed index instead of rebuilding a multi-gigabyte
        # minimap2 index in every analysis directory.
        return configured_index
    fasta_name = os.path.basename(config["reference"]["fasta"])
    preset = config.get("minimap2_align", {}).get("preset", "map-ont")
    return "alignment/minimap2_index/{}.{}.mmi".format(fasta_name, preset)


def get_index_preset():
    """The preset the available minimap2 index was built with."""

    return str(
        config.get("minimap2_align", {}).get(
            "index_preset", config.get("minimap2_align", {}).get("preset", "map-ont")
        )
    ).strip()


def get_reference_target_for_unit(wildcards):
    """The minimap2 target for a unit: the index when it matches, else the FASTA.

    A prebuilt index carries its own -k/-w, and minimap2 silently keeps them:
    aligning with -x splice against a map-ont index runs the splice model on
    map-ont seeding.  Rather than align on parameters nobody chose, spliced
    units are pointed at the FASTA and minimap2 indexes as it goes.
    """

    preset = get_minimap2_preset_for_unit(wildcards)
    if preset == get_index_preset():
        return get_reference_index()
    return config["reference"]["fasta"]


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


# ---------------------------------------------------------------------------
# Sample routing
# ---------------------------------------------------------------------------
#
# Two independent facts decide what runs for a sample:
#
#   library_type  what the molecules are, and therefore how they align
#   assay         what the sample was sequenced for, and therefore what runs
#
# The logic lives in workflow/scripts/nagini_routing.py so that it is testable
# without Snakemake.  The wrappers below bind it to this workflow's `units`
# table and resolved assay registry.

try:
    assays = nagini_routing.resolve_assays(config)
except nagini_routing.RoutingError as error:
    raise WorkflowError(str(error))

VALID_LIBRARY_TYPES = set(nagini_routing.LIBRARY_TYPES)
VALID_ASSAYS = set(assays)


def _routing_call(function, *args):
    """Turn a routing error into a Snakemake error with the same message."""

    try:
        return function(*args)
    except nagini_routing.RoutingError as error:
        raise WorkflowError(str(error))


def get_assay(sample):
    """Return the assay name for a sample."""

    return _routing_call(nagini_routing.get_assay, units, sample, assays)


def get_assay_definition(sample):
    return _routing_call(nagini_routing.get_assay_definition, units, sample, assays)


def get_library_type(sample):
    """Return the library type for a sample (units.tsv wins over the assay)."""

    return _routing_call(nagini_routing.get_library_type, units, sample, assays)


def sample_needs_dna_variant_calling(sample):
    """Whether SNV/indel calling runs for this sample."""

    return _routing_call(nagini_routing.sample_needs_snv_calling, units, sample, assays)


def sample_needs_sv_calling(sample):
    """Whether structural variant calling runs for this sample."""

    return _routing_call(nagini_routing.sample_needs_sv_calling, units, sample, assays)


def sample_needs_rna_analysis(sample):
    """Whether the library is spliced and needs RNA-aware handling."""

    return nagini_routing.library_type_is_spliced(get_library_type(sample))


def sample_is_target_restricted(sample):
    """Whether analysis for this sample is confined to the design BED."""

    return _routing_call(nagini_routing.sample_is_target_restricted, units, sample, assays)


def sample_has_assay_module(sample, module):
    """Whether an assay-specific analysis module runs for this sample."""

    return _routing_call(nagini_routing.sample_has_module, units, sample, assays, module)


def samples_with_assay_module(module):
    """Every sample whose assay declares this module."""

    return [
        str(sample).strip()
        for sample in get_samples(samples)
        if sample_has_assay_module(str(sample).strip(), module)
    ]


def get_resolved_design_bed():
    """Return the resolved absolute path to design_bed if it exists, else empty string."""
    bed = str(config.get("reference", {}).get("design_bed", "") or "").strip()
    if not bed:
        return ""
    if os.path.isfile(bed):
        return os.path.abspath(bed)
    if not os.path.isabs(bed):
        repo_relative = os.path.normpath(os.path.join(workflow.basedir, "..", bed))
        if os.path.isfile(repo_relative):
            return repo_relative
    return ""


def has_design_bed():
    """Return True if a valid design_bed is configured and exists."""
    return bool(get_resolved_design_bed())


def get_minimap2_preset_for_unit(wildcards):
    """Return the minimap2 preset for a unit, from its library type."""

    return nagini_routing.minimap2_preset(
        get_library_type(wildcards.sample),
        config.get("minimap2_align", {}).get("preset", "map-ont"),
    )


def get_assay_setting(sample, section, key, default=None):
    """A tool setting for a sample, with its assay's override applied."""

    return _routing_call(
        nagini_routing.assay_setting, units, sample, assays, section, key, default
    )


def get_clairs_to_platform(wildcards):
    """Derive the ClairS-TO error model from the sample's basecalling model."""

    configured = str(config.get("clairs_to_call", {}).get("platform", "") or "").strip()
    if configured:
        # An explicit platform in the config is an override, not a suggestion.
        return configured
    model = nagini_routing.get_basecalling_model(units, wildcards.sample)
    platform, _warning = nagini_routing.clairs_to_platform(
        model,
        config.get("clairs_to_call", {}).get("platform_map", {}),
        config.get("clairs_to_call", {}).get("default_platform", ""),
    )
    return platform


def get_clairs_to_snv_min_af(wildcards):
    return get_assay_setting(
        wildcards.sample,
        "clairs_to_call",
        "snv_min_af",
        config.get("clairs_to_call", {}).get("snv_min_af", 0.05),
    )


def get_clairs_to_indel_min_af(wildcards):
    return get_assay_setting(
        wildcards.sample,
        "clairs_to_call",
        "indel_min_af",
        config.get("clairs_to_call", {}).get("indel_min_af", 0.1),
    )


# Fail early on an unroutable sample sheet rather than midway through the DAG,
# and surface a library_type that disagrees with its assay as a warning.
_routing_summary = _routing_call(
    nagini_routing.routing_summary, units, list(get_samples(samples)), assays
)
for _warning in _routing_call(
    nagini_routing.routing_warnings, units, list(get_samples(samples)), assays
):
    sys.stderr.write(f"WARNING: {_warning}\n")

if not str(config.get("clairs_to_call", {}).get("platform", "") or "").strip():
    _seen_models = set()
    for _sample in get_samples(samples):
        _model = nagini_routing.get_basecalling_model(units, str(_sample).strip())
        if _model in _seen_models:
            continue
        _seen_models.add(_model)
        _platform, _model_warning = nagini_routing.clairs_to_platform(
            _model,
            config.get("clairs_to_call", {}).get("platform_map", {}),
            config.get("clairs_to_call", {}).get("default_platform", ""),
        )
        if _model_warning:
            sys.stderr.write(f"WARNING: {_model_warning}\n")


def get_small_variant_caller(sample=None, unit_type=None):
    """Determine the small variant caller for a sample/unit or the run default."""
    if sample is None:
        configured = str(config.get("small_variant_caller", "auto")).strip().lower()
        return configured if configured in ("clairs_to", "deepvariant", "deepsomatic", "mutect2", "both", "auto") else "auto"
    return _routing_call(
        nagini_routing.get_sample_caller,
        units,
        sample,
        assays,
        config.get("small_variant_caller", "auto"),
        unit_type or "",
    )


def deepvariant_enabled(sample=None, unit_type=None):
    caller = get_small_variant_caller(sample, unit_type)
    return caller in ("deepvariant", "both", "auto")


def clairs_to_enabled(sample=None, unit_type=None):
    caller = get_small_variant_caller(sample, unit_type)
    return caller in ("clairs_to", "both", "auto")


def mutect2_enabled(sample=None, unit_type=None):
    caller = get_small_variant_caller(sample, unit_type)
    return caller in ("mutect2", "both")


def get_clairs_to_bed(wildcards):
    """Return the design BED if configured, otherwise whole-genome contigs BED."""
    resolved = get_resolved_design_bed()
    if resolved:
        return resolved
    return "references/genome_contigs.bed"


# For target-restricted panels, alert if design_bed is missing
for _sample in get_samples(samples):
    _s_str = str(_sample).strip()
    if sample_is_target_restricted(_s_str) and sample_needs_dna_variant_calling(_s_str) and not has_design_bed():
        sys.stderr.write(
            f"WARNING: Sample {_s_str} is target_restricted but no reference.design_bed is configured.\n"
        )


def get_sample_type_path(wildcards, directory, suffix):
    """Build report paths from normalized sample/type components.

    Keep path construction in one place.  In particular, never use implicit
    multiline string concatenation for these paths: indentation or a visual
    separator such as ``" _ "`` would become a literal filename character and
    Snakemake would only report the problem later as a missing input.
    """

    sample = str(wildcards.sample).strip()
    unit_type = str(wildcards.type).strip()
    if not sample or not unit_type:
        raise WorkflowError(
            f"Cannot build a report path with empty sample/type: {sample}/{unit_type}."
        )
    if any(character.isspace() for character in directory + suffix):
        raise WorkflowError(
            f"Internal whitespace in report path components: {directory}, {suffix}."
        )
    return f"{directory}/{sample}_{unit_type}{suffix}"


def get_small_variant_report_vcf(wildcards):
    """Return the normalized VCF, preferring VEP annotation when enabled."""
    suffix = ".filtered.norm.vep.vcf.gz" if vep_enabled() else ".filtered.norm.vcf.gz"
    caller = get_small_variant_caller(wildcards.sample, wildcards.type)
    if caller == "deepvariant":
        caller_dir = "snv_indels/deepvariant_postprocess"
    elif caller == "mutect2":
        caller_dir = "snv_indels/mutect2_postprocess"
    else:
        caller_dir = "snv_indels/clairs_to_postprocess"
    return get_sample_type_path(wildcards, caller_dir, suffix)



def get_small_variant_report_vcf_index(wildcards):
    return f"{get_small_variant_report_vcf(wildcards)}.tbi"


def get_report_snv_vcf(wildcards):
    if sample_needs_dna_variant_calling(wildcards.sample):
        return get_small_variant_report_vcf(wildcards)
    return []


def get_report_snv_tbi(wildcards):
    if sample_needs_dna_variant_calling(wildcards.sample):
        return get_small_variant_report_vcf_index(wildcards)
    return []


def get_report_sv_tsv(wildcards):
    if sample_needs_sv_calling(wildcards.sample):
        return get_sample_type_path(wildcards, "reports/sv_summary", ".tsv")
    return []


def get_report_translocations_tsv(wildcards):
    if sample_needs_sv_calling(wildcards.sample):
        return get_sample_type_path(
            wildcards, "reports/sv_summary", ".translocations.tsv"
        )
    return []


def get_report_target_coverage(wildcards):
    if has_design_bed():
        return get_sample_type_path(
            wildcards, "qc/mosdepth_bed", ".mosdepth.summary.txt"
        )
    return []


def get_report_target_dist(wildcards):
    if has_design_bed():
        return get_sample_type_path(
            wildcards, "qc/mosdepth_bed", ".mosdepth.global.dist.txt"
        )
    return []


def get_report_demux_summary(wildcards):
    if dorado_barcodes:
        return "qc/dorado_demux/dorado_demux_mqc.tsv"
    return []


def get_sample_barcode_name(wildcards):
    rows = _unit_rows(wildcards.sample, wildcards.type)
    if not rows.empty:
        barcodes = [
            str(b).strip()
            for b in rows.get("barcode_name", pd.Series(dtype=str)).tolist()
            if str(b).strip()
        ]
        if barcodes:
            return barcodes[0]
        bc = [
            str(b).strip()
            for b in rows.get("barcode", pd.Series(dtype=str)).tolist()
            if str(b).strip()
        ]
        if bc:
            return bc[0]
    return ""


def get_report_bcr_abl_fusion(wildcards):
    if sample_has_assay_module(wildcards.sample, "bcr_abl"):
        return get_sample_type_path(
            wildcards, "reports/bcr_abl", ".fusion_summary.tsv"
        )
    return []


def get_report_bcr_abl_mutations(wildcards):
    if sample_has_assay_module(wildcards.sample, "bcr_abl"):
        return get_sample_type_path(
            wildcards, "reports/bcr_abl", ".branford_mutations.tsv"
        )
    return []


def get_report_bcr_abl_phasing(wildcards):
    if sample_has_assay_module(wildcards.sample, "bcr_abl"):
        return get_sample_type_path(wildcards, "reports/bcr_abl", ".phasing.tsv")
    return []


def get_report_bcr_abl_exon_coverage(wildcards):
    if sample_has_assay_module(wildcards.sample, "bcr_abl"):
        return get_sample_type_path(wildcards, "reports/bcr_abl", ".exon_coverage.tsv")
    return []



def get_unit_types_for_sample(sample):
    return sorted(set(units.loc[units["sample"] == sample, "type"]))


# Conditions understood by `when:` in config/output_files.yaml.  Pipeline-level
# ones are evaluated once; sample-level ones are evaluated per sample.  Any
# other value falls back to a plain config flag lookup, which is how the
# existing annotate_* switches work.
_PIPELINE_CONDITIONS = {
    "annotate_small_variants_with_vep": lambda: vep_enabled(),
    "clairs_to_enabled": lambda: clairs_to_enabled(),
    "clairs_to_vep_enabled": lambda: clairs_to_enabled() and vep_enabled(),
    "deepvariant_enabled": lambda: deepvariant_enabled(),
    "deepvariant_vep_enabled": lambda: deepvariant_enabled() and vep_enabled(),
    "deepsomatic_enabled": lambda: get_small_variant_caller() in ("deepsomatic", "both", "auto"),
    "deepsomatic_vep_enabled": lambda: get_small_variant_caller() in ("deepsomatic", "both", "auto") and vep_enabled(),
    "mutect2_enabled": lambda: mutect2_enabled(),
    "mutect2_vep_enabled": lambda: mutect2_enabled() and vep_enabled(),
    "annotate_small_variants": lambda: bool(config.get("annotate_small_variants", False)),
    "annotate_svs": lambda: bool(config.get("annotate_svs", False)),
    "has_design_bed": lambda: has_design_bed(),
}

_SAMPLE_CONDITIONS = {
    "clairs_to_enabled": lambda s: clairs_to_enabled(s) and sample_needs_dna_variant_calling(s),
    "clairs_to_vep_enabled": lambda s: clairs_to_enabled(s) and sample_needs_dna_variant_calling(s),
    "deepvariant_enabled": lambda s: deepvariant_enabled(s) and sample_needs_dna_variant_calling(s),
    "deepvariant_vep_enabled": lambda s: deepvariant_enabled(s) and sample_needs_dna_variant_calling(s),
    "deepsomatic_enabled": lambda s: sample_needs_dna_variant_calling(s),
    "deepsomatic_vep_enabled": lambda s: sample_needs_dna_variant_calling(s),
    "mutect2_enabled": lambda s: mutect2_enabled(s) and sample_needs_dna_variant_calling(s),
    "mutect2_vep_enabled": lambda s: mutect2_enabled(s) and sample_needs_dna_variant_calling(s),
    "annotate_small_variants": sample_needs_dna_variant_calling,
    "annotate_small_variants_with_vep": sample_needs_dna_variant_calling,
    "sample_needs_dna_variant_calling": sample_needs_dna_variant_calling,
    "sv_calling": sample_needs_sv_calling,
    "annotate_svs": sample_needs_sv_calling,
    "sample_needs_sv_calling": sample_needs_sv_calling,
    "sample_needs_rna_analysis": sample_needs_rna_analysis,
    "target_restricted": sample_is_target_restricted,
}


def _pipeline_condition_met(condition):
    """Whether a `when:` condition holds for the run as a whole."""

    if not condition:
        return True
    if condition in _PIPELINE_CONDITIONS:
        return _PIPELINE_CONDITIONS[condition]()
    if condition.startswith(("assay:", "assay_module:")):
        # Resolved per sample; nothing to decide at pipeline level.
        return True
    if condition in _SAMPLE_CONDITIONS:
        return True
    return bool(config.get(condition, False))


def _sample_condition_met(condition, sample):
    """Whether a `when:` condition holds for one sample."""

    if not condition:
        return True
    if condition.startswith("assay_module:"):
        return sample_has_assay_module(sample, condition.split(":", 1)[1])
    if condition.startswith("assay:"):
        return get_assay(sample) == nagini_routing.normalize(condition.split(":", 1)[1])
    predicate = _SAMPLE_CONDITIONS.get(condition)
    return predicate(sample) if predicate else True


# A typo in `when:` used to fall through to a config lookup and silently drop
# the output.  Fail at parse time instead.
_unknown_conditions = nagini_routing.unknown_conditions(output_spec, assays)
if _unknown_conditions:
    raise WorkflowError(
        "Unknown `when:` condition(s) in {}: {}. Valid values: a config flag, "
        "one of {}, or assay:<name> / assay_module:<name>.".format(
            config["output"],
            ", ".join(_unknown_conditions),
            ", ".join(
                sorted(
                    nagini_routing.PIPELINE_CONDITION_NAMES
                    | nagini_routing.SAMPLE_CONDITION_NAMES
                )
            ),
        )
    )


def compile_output_file_list(wildcards):
    outdir = pathlib.Path(output_spec.get("directory", "./"))
    output_files = []
    for file_definition in output_spec["files"]:
        condition = file_definition.get("when")
        allowed_types = file_definition.get("types")

        if not _pipeline_condition_met(condition):
            continue

        for sample in get_samples(samples):
            clean_sample = str(sample).strip()
            if not _sample_condition_met(condition, clean_sample):
                continue

            for unit_type in get_unit_types_for_sample(sample):
                if allowed_types and unit_type not in allowed_types:
                    continue
                clean_type = str(unit_type).strip()
                template = str(file_definition["output"]).strip()
                output_files.append(
                    outdir
                    / template.format(
                        sample=clean_sample,
                        type=clean_type,
                    )
                )

    if dorado_barcodes:
        output_files.append(outdir / "basecalling/dorado_demux/demux_summary.tsv")
        output_files.append(outdir / "qc/dorado_demux/dorado_demux_mqc.tsv")
        if dorado_demux_output_format in ("bam", "both"):
            output_files.extend(
                outdir / f"basecalling/dorado_demux/{barcode}.bam"
                for barcode in dorado_barcodes
            )
        if dorado_demux_output_format in ("fastq", "both"):
            output_files.extend(
                outdir / f"basecalling/dorado_fastq/{barcode}.fastq.gz"
                for barcode in dorado_barcodes
            )
    return output_files
