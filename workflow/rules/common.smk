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

# A reference entry with no FASTA means the build is declared but not installed
# at this site (the registry says so with an empty path). Catch it here, where
# the message can name the build, rather than letting Snakemake report an empty
# input file for whichever rule happens to be scheduled first.
if not str(config["reference"].get("fasta", "") or "").strip():
    _installed = sorted(
        name
        for name, entry in references.items()
        if isinstance(entry, dict) and str(entry.get("fasta", "") or "").strip()
    )
    raise WorkflowError(
        f"reference_build={reference_build} has no fasta configured, which means it "
        f"is not installed at site {config.get('site', '?')}. Builds available here: "
        f"{', '.join(_installed) or 'none'}. Add the paths to "
        "config/reference_files/reference_files_<site>.yaml, or pick another build."
    )

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
    if not cores and hasattr(workflow, "resource_settings") and getattr(workflow.resource_settings, "cores", None):
        try:
            cores = int(workflow.resource_settings.cores)
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
        # Both rindels rules read this section, and both are a single
        # compute-heavy job per sample.
        "rindels",
        "cnv_rindels",
        "vep",
    )
    for rule_name in full_core_rules:
        if rule_name in cfg and not isinstance(cfg[rule_name], dict):
            continue
        cfg.setdefault(rule_name, {})["threads"] = cores

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


# Every per-unit alignment carries a token naming where its reads came from.
# It is not decoration: one rule per input type writes into alignment/units
# and Snakemake needs their output patterns to differ, while a single run may mix
# input types across units so a ruleorder cannot decide it. Reading as
# provenance ("from_ubam") rather than as a state ("unmapped_bam") keeps the
# name from contradicting the aligned file it is on.
def get_unit_alignment_path(row):
    suffix = {
        "fastq": "from_fastq",
        "unmapped_bam": "from_ubam",
        "aligned_bam": "from_aligned",
        "dorado_bam": "from_dorado",
    }.get(str(row["input_type"]).strip().lower())
    if suffix is None:
        raise WorkflowError(
            f"Unsupported input_type={row['input_type']}; use fastq, unmapped_bam, aligned_bam or dorado_bam."
        )
    return "alignment/units/{}_{}_{}_{}.{}.bam".format(
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


# Input types that carry their reads in a BAM that has not been aligned yet,
# and can therefore be turned back into FASTQ without undoing an alignment.
UNALIGNED_BAM_INPUT_TYPES = ("unmapped_bam", "dorado_bam")


def sample_source_bams_are_unaligned(sample, unit_type=None):
    """Whether every unit of this sample came in as an unaligned BAM."""

    rows = _unit_rows(sample, unit_type) if unit_type else units[
        units["sample"].astype(str).str.strip() == str(sample).strip()
    ]
    if rows.empty:
        return False
    return all(
        str(row["input_type"]).strip().lower() in UNALIGNED_BAM_INPUT_TYPES
        for _, row in rows.iterrows()
    )


def sample_has_merged_fastq(sample):
    """Whether one FASTQ per sample is asked for and can be produced.

    The FASTQ is built from the source BAMs rather than from the merged
    alignment: re-extracting reads from an aligned BAM means undoing the
    reverse-complementing minimap2 did, which the MM/ML modified-base tags
    have to follow exactly.  Going back to the unaligned input sidesteps that
    entirely, and costs one pass over files the run already reads.
    """

    return bool(config.get("merged_fastq", False)) and sample_source_bams_are_unaligned(
        sample
    )


def get_sample_source_bams(wildcards):
    """Every unaligned source BAM for one sample, across all its flowcells."""

    rows = _unit_rows(wildcards.sample, wildcards.type)
    if rows.empty:
        raise WorkflowError(f"No units.tsv rows for {wildcards.sample}/{wildcards.type}.")
    bams = []
    for _, row in rows.iterrows():
        input_type = str(row["input_type"]).strip().lower()
        if input_type not in UNALIGNED_BAM_INPUT_TYPES:
            raise WorkflowError(
                f"{wildcards.sample}/{wildcards.type} has a unit with "
                f"input_type={input_type}; a merged FASTQ needs unaligned BAM input."
            )
        bam = str(row.get("bam", "") or "").strip()
        if not bam:
            raise WorkflowError(
                f"units.tsv column bam is empty for {wildcards.sample}/{wildcards.type}/"
                f"{row['processing_unit']}."
            )
        bams.append(bam)
    return bams


def get_barcode_unit_alignment(wildcards):
    """Return the per-unit alignment for one preserved barcode."""

    barcode = str(wildcards.barcode).strip()
    if "barcode_name" not in units.columns:
        raise WorkflowError("units.tsv has no barcode_name column.")
    rows = _unit_rows(wildcards.sample, wildcards.type)
    rows = rows[rows["barcode_name"].astype(str).str.strip() == barcode]
    if len(rows) != 1:
        raise WorkflowError(
            "Expected exactly one unit for "
            f"{wildcards.sample}/{wildcards.type}/{barcode}; found {len(rows)}."
        )
    return get_unit_alignment_path(rows.iloc[0])


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


def get_gene_annotation(wildcards=None):
    return _optional_reference_path("gene_annotation")


def get_rindels_annotation_args():
    """--gtf and coding-only annotation for Rindels, when a GTF is installed.

    Deliberately not --exonic: that limits the whole analysis to exons, which
    would gut the genome-wide read-count binning the CNV segmentation is built
    on. Annotating instead keeps the call set genome-wide and lets the report
    pick out the coding and splice variants.
    """

    gtf = get_gene_annotation()
    if not gtf:
        return ""
    return f"--gtf {gtf} --annotate --annot-coding-only"


def get_clinvar_vcf(wildcards=None):
    return _optional_reference_path("clinvar")


def get_cytoband_file(wildcards=None):
    return _optional_reference_path("cytoband")


def get_gaps_bed(wildcards=None):
    """Assembly gaps for the ideogram, when the build registers a track.

    chm13 is gapless, so leaving this unset is the correct answer for T2T:
    Rindels then draws none rather than GRCh38's unresolved arms.
    """

    return _optional_reference_path("gaps_bed")


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


_design_chromosomes_cache = []


def get_design_chromosomes():
    """Reference chromosomes that the design BED actually covers.

    A target-restricted sample has nothing to call on a chromosome the design
    never touches, so those shards should not be created at all rather than
    run and produce an empty VCF.  Resolved on first use, since the design BED
    helpers are defined further down this file.
    """

    if _design_chromosomes_cache:
        return _design_chromosomes_cache[0]
    bed = get_resolved_design_bed()
    if not bed:
        _design_chromosomes_cache.append(list(reference_chromosomes))
        return _design_chromosomes_cache[0]
    present = set()
    with open(bed) as handle:
        for line in handle:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            present.add(line.split("\t", 1)[0].strip())
    covered = [name for name in reference_chromosomes if name in present]
    _design_chromosomes_cache.append(covered or list(reference_chromosomes))
    return _design_chromosomes_cache[0]


def sample_calls_within_design(sample):
    """Whether small-variant calling for this sample is confined to the design."""

    return has_design_bed() and sample_is_target_restricted(sample)


def get_snv_shard_chromosomes(sample):
    """The chromosomes a per-chromosome caller should be sharded over."""

    if sample_calls_within_design(sample):
        return get_design_chromosomes()
    return reference_chromosomes


def get_deepvariant_regions_bed(wildcards):
    """The design slice DeepVariant should be confined to, or nothing."""

    if not sample_calls_within_design(wildcards.sample):
        return []
    return f"snv_indels/bed_split/design_bedfile_{wildcards.chr}.bed"


def get_deepvariant_regions(wildcards):
    """--regions for DeepVariant: the design slice, else the whole chromosome.

    The Hydra rule always passes the bare chromosome, so a target-restricted
    sample was called genome-wide while ClairS-TO ran only inside the design.
    DeepVariant's --regions takes a BED, and snv_indels' bed_split already
    produces one slice per chromosome for Mutect2, so reuse that.
    """

    bed = get_deepvariant_regions_bed(wildcards)
    return f" --regions {bed} " if bed else f" --regions {wildcards.chr} "


def get_deepvariant_shards(wildcards, template):
    """The per-chromosome DeepVariant files to concatenate for one sample."""

    return [
        template.format(sample=wildcards.sample, type=wildcards.type, chr=name)
        for name in get_snv_shard_chromosomes(wildcards.sample)
    ]


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


def sample_needs_cnv_calling(sample):
    """Whether whole-genome CNV calling runs for this sample."""
    if not config.get("cnv_calling", True):
        return False
    return _routing_call(nagini_routing.sample_needs_cnv_calling, units, sample, assays)


def cnvkit_enabled(sample):
    """Whether CNVkit runs for this sample.

    CNVkit runs whenever CNV calling is enabled and cnv_engine is not set to
    'rindels_only'. Both Rindels and CNVkit can run in parallel.
    """
    if not sample_needs_cnv_calling(sample):
        return False
    return config.get("cnv_engine", "both") != "rindels_only"


def sniffles2_enabled(sample):
    """Whether Sniffles2 SV calling runs for this sample.

    Sniffles2 runs whenever SV calling is enabled and sv_engine is not set to
    'savana_only'. Both Sniffles2 and SAVANA run in parallel as complements.
    """
    if not sample_needs_sv_calling(sample):
        return False
    return config.get("sv_engine", "both") != "savana_only"


def savana_enabled(sample):
    """Whether SAVANA SV calling runs for this sample.

    SAVANA runs whenever SV calling is enabled and sv_engine is not set to
    'sniffles2_only'. Both Sniffles2 and SAVANA run in parallel as complements.
    """
    if not sample_needs_sv_calling(sample):
        return False
    return config.get("sv_engine", "both") != "sniffles2_only"


def get_savana_genome_version(wildcards=None):
    """Return genome version argument for SAVANA tumor-only mode."""
    configured = config.get("savana_ont_to", {}).get("genome_version", "")
    if configured:
        return configured
    build = str(config.get("reference_build", "hg38")).strip().lower()
    if "19" in build or "37" in build:
        return "1000g_hg19"
    elif "hs1" in build or "t2t" in build or "chm13" in build:
        return "1000g_t2t"
    return "1000g_hg38"


def _find_rindels_binary():
    """Locate the Rindels binary from config, PATH, ~/bin, or dev checkout."""
    override = config.get("rindels", {}).get("binary", "") or config.get("bcr_abl", {}).get("rindels_binary", "")
    if override and (os.path.isfile(override) or shutil.which(override)):
        return override
    candidates = [
        os.path.expanduser("~/bin/rindels"),
        "/home/cgu-lab/bin/rindels",
        os.path.expanduser("~/tmp/tools/rindels/target/release/rindels"),
        "/home/cgu-lab/tmp/tools/rindels/target/release/rindels",
        os.path.expanduser("~/.cargo/bin/rindels"),
    ]
    for c in candidates:
        if os.path.isfile(c) and os.access(c, os.X_OK):
            return c
    which_r = shutil.which("rindels")
    if which_r:
        return which_r
    return "rindels"


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


def sample_has_target_metrics(sample):
    """Whether target-region outputs are meaningful for this sample.

    A design BED existing in the config is not enough: in a run that mixes an
    adaptive-sampling panel with a plain WGS control, the control was never
    enriched for those regions, so a target BAM and target coverage for it
    would be an artefact of the config rather than of the experiment.
    """

    return has_design_bed() and sample_is_target_restricted(sample)


def get_sample_design_bed(sample):
    """The design BED for one sample, empty unless its assay is confined to it.

    Annotations keyed off the design (CNV target overlap, the report header)
    must not claim a WGS sample was sequenced against a panel it shares a
    flowcell with.
    """

    return get_resolved_design_bed() if sample_has_target_metrics(sample) else ""


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
    """Return the design BED for target-restricted samples if configured, otherwise whole-genome contigs BED."""
    if sample_is_target_restricted(wildcards.sample):
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


def validate_design_bed():
    """Reject a malformed design BED here rather than inside a caller.

    A single blank line in the T2T design was enough to make ClairS-TO exit
    with "list index out of range" after the whole alignment stage had already
    run.  The BED is passed verbatim to tools whose parsers assume three
    fields per line, so check it once, at DAG build, where it costs nothing.
    """

    bed = get_resolved_design_bed()
    if not bed:
        return
    problems = []
    with open(bed) as handle:
        for number, line in enumerate(handle, start=1):
            stripped = line.rstrip("\n")
            if stripped.startswith(("#", "track", "browser")):
                continue
            if not stripped.strip():
                problems.append(f"line {number}: blank line")
                continue
            fields = stripped.split("\t")
            if len(fields) < 3:
                problems.append(
                    f"line {number}: {len(fields)} tab-separated field(s), need at least 3"
                )
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                problems.append(f"line {number}: non-integer start/end {fields[1:3]}")
                continue
            if end <= start:
                problems.append(f"line {number}: end {end} is not after start {start}")
    if problems:
        raise WorkflowError(
            "Malformed design BED {}:\n  {}\nTools that take this file verbatim "
            "(ClairS-TO, mosdepth, samtools view -L) fail on it.".format(
                bed, "\n  ".join(problems)
            )
        )


validate_design_bed()


def get_sample_type_path(wildcards, directory, suffix):
    """Build report paths from normalized sample/type components.

    Keep path construction in one place.  In particular, never use implicit
    multiline string concatenation for these paths: indentation or a visual
    separator such as ``" _ "`` would become a literal filename character and
    Snakemake would only report the problem later as a missing input.
    """

    sample = str(wildcards.sample).strip()
    unit_type = str(wildcards.type).strip()
    d = str(directory).strip()
    sfx = str(suffix).strip()
    if not sample or not unit_type:
        raise WorkflowError(
            f"Cannot build a report path with empty sample/type: {sample}/{unit_type}."
        )
    res = f"{d}/{sample}_{unit_type}{sfx}"
    while " /" in res or "/ " in res or " _" in res or "_ " in res or " ." in res:
        res = res.replace(" /", "/").replace("/ ", "/").replace(" _", "_").replace("_ ", "_").replace(" .", ".")
    return res.strip()



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
    if sniffles2_enabled(wildcards.sample):
        return get_sample_type_path(wildcards, "reports/sv_summary", ".tsv")
    return []


def sample_has_rindels_sv(sample):
    """Whether a Rindels structural call set exists for this sample.

    Rindels runs for CNV, and in --sv-lite it also reports DEL/DUP/INV/TRA, so
    the second call set exists exactly when CNV calling does and the mode is
    not the CNV-only one that disables SV.
    """

    if not sample_needs_cnv_calling(sample):
        return False
    extra = str(config.get("rindels", {}).get("extra_cnv", "")).lower()
    return "--only-cnv" not in extra and "--no-sv" not in extra


def get_report_target_cramino(wildcards):
    """Read statistics over the target BAM, for target-restricted samples."""

    if sample_has_target_metrics(wildcards.sample):
        return get_sample_type_path(wildcards, "qc/cramino_target", ".txt")
    return []


def get_report_cnv_ideogram(wildcards):
    """The ideogram PNG for a sample, when Rindels produced the bins."""

    if sample_needs_cnv_calling(wildcards.sample):
        return get_sample_type_path(wildcards, "reports/cnv_ideogram", ".ideogram.png")
    return []


def get_report_rindels_sv_tsv(wildcards):
    if sample_has_rindels_sv(wildcards.sample):
        return get_sample_type_path(wildcards, "reports/sv_summary_rindels", ".tsv")
    return []


def get_report_rindels_translocations_tsv(wildcards):
    if sample_has_rindels_sv(wildcards.sample):
        return get_sample_type_path(
            wildcards, "reports/sv_summary_rindels", ".translocations.tsv"
        )
    return []


def get_report_savana_sv_tsv(wildcards):
    if savana_enabled(wildcards.sample):
        return get_sample_type_path(wildcards, "reports/sv_summary_savana", ".tsv")
    return []


def get_report_savana_translocations_tsv(wildcards):
    if savana_enabled(wildcards.sample):
        return get_sample_type_path(
            wildcards, "reports/sv_summary_savana", ".translocations.tsv"
        )
    return []


def get_report_translocations_tsv(wildcards):
    if sniffles2_enabled(wildcards.sample):
        return get_sample_type_path(
            wildcards, "reports/sv_summary", ".translocations.tsv"
        )
    return []


def get_report_cnv_tsv(wildcards):
    if sample_needs_cnv_calling(wildcards.sample):
        return get_sample_type_path(
            wildcards, "reports/cnv_summary", ".cnv.tsv"
        )
    return []


def get_report_target_coverage(wildcards):
    if sample_has_target_metrics(wildcards.sample):
        return get_sample_type_path(
            wildcards, "qc/mosdepth_bed", ".mosdepth.summary.txt"
        )
    return []


def get_report_target_dist(wildcards):
    if sample_has_target_metrics(wildcards.sample):
        return get_sample_type_path(
            wildcards, "qc/mosdepth_bed", ".mosdepth.global.dist.txt"
        )
    return []


def get_report_demux_summary(wildcards):
    if dorado_barcodes:
        return "qc/dorado_demux/dorado_demux_mqc.tsv"
    return []


def get_sample_dna_number(sample):
    """The lab identifier for a sample, empty when samples.tsv has none.

    Reporting it alongside the pipeline's own sample id lets a run keep the
    barcode-derived names that every file is built from while the workbook
    still says which DNA the result belongs to.
    """

    if "dna_number" not in samples.columns:
        return ""
    key = str(sample).strip()
    if key not in samples.index:
        return ""
    value = samples.loc[key, "dna_number"]
    if hasattr(value, "iloc"):
        value = value.iloc[0]
    return str(value or "").strip()


def get_sample_barcode_name(wildcards):
    return get_sample_barcode_name_str(wildcards.sample, wildcards.type)


def get_sample_barcode_name_str(sample, sample_type):
    rows = _unit_rows(sample, sample_type)
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


def get_report_bcr_abl_plot_png(wildcards):
    if sample_has_assay_module(wildcards.sample, "bcr_abl"):
        return get_sample_type_path(wildcards, "reports/bcr_abl", ".fusion_plot.png")
    return []


def get_report_amplicon_stats(wildcards):
    if sample_has_assay_module(wildcards.sample, "amplicon"):
        return get_sample_type_path(wildcards, "reports/amplicon", ".amplicon_stats.tsv")
    return []


def get_report_amplicon_summary(wildcards):
    if sample_has_assay_module(wildcards.sample, "amplicon"):
        return get_sample_type_path(wildcards, "reports/amplicon", ".amplicon_summary.tsv")
    return []


def get_report_amplicon_itd(wildcards):
    if sample_has_assay_module(wildcards.sample, "amplicon"):
        return get_sample_type_path(wildcards, "reports/amplicon", ".flt3_itd.tsv")
    return []


def get_report_amplicon_vcf(wildcards):
    if sample_has_assay_module(wildcards.sample, "amplicon"):
        return get_sample_type_path(wildcards, "snv_indels/rindels", ".amplicon.vcf.gz")
    return []


def get_report_amplicon_pileup(wildcards):
    if sample_has_assay_module(wildcards.sample, "amplicon"):
        return get_sample_type_path(wildcards, "reports/amplicon", ".pileup_counts.tsv")
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
    "target_metrics": lambda: any(
        sample_has_target_metrics(s) for s in get_samples(samples)
    ),
    "merged_fastq": lambda: any(
        sample_has_merged_fastq(s) for s in get_samples(samples)
    ),
    "rindels_sv": lambda: any(
        sample_has_rindels_sv(s) for s in get_samples(samples)
    ),
    "cnv_calling": lambda: any(sample_needs_cnv_calling(s) for s in get_samples(samples)),
    "cnvkit_enabled": lambda: any(cnvkit_enabled(s) for s in get_samples(samples)),
    "sniffles2_enabled": lambda: any(sniffles2_enabled(s) for s in get_samples(samples)),
    "savana_enabled": lambda: any(savana_enabled(s) for s in get_samples(samples)),
    "enable_bcftools_pileup": lambda: bool(config.get("bcr_abl", {}).get("enable_bcftools_pileup", False) or config.get("enable_bcftools_pileup", False)),
    "enable_rindels": lambda: bool(config.get("bcr_abl", {}).get("enable_rindels", False) or config.get("enable_rindels", False)),
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
    "sniffles2_enabled": sniffles2_enabled,
    "savana_enabled": savana_enabled,
    "cnv_calling": sample_needs_cnv_calling,
    "sample_needs_cnv_calling": sample_needs_cnv_calling,
    "cnvkit_enabled": cnvkit_enabled,
    "sample_needs_rna_analysis": sample_needs_rna_analysis,
    "target_restricted": sample_is_target_restricted,
    "target_metrics": sample_has_target_metrics,
    "merged_fastq": sample_has_merged_fastq,
    "rindels_sv": sample_has_rindels_sv,
    "enable_bcftools_pileup": lambda s: sample_has_assay_module(s, "bcr_abl") and bool(config.get("bcr_abl", {}).get("enable_bcftools_pileup", False) or config.get("enable_bcftools_pileup", False)),
    "enable_rindels": lambda s: sample_has_assay_module(s, "bcr_abl") and bool(config.get("bcr_abl", {}).get("enable_rindels", False) or config.get("enable_rindels", False)),
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

        template = str(file_definition["output"]).strip()
        if "{sample}" in template:
            for sample in get_samples(samples):
                clean_sample = str(sample).strip()
                if not _sample_condition_met(condition, clean_sample):
                    continue

                for unit_type in get_unit_types_for_sample(sample):
                    if allowed_types and unit_type not in allowed_types:
                        continue
                    clean_type = str(unit_type).strip()
                    output_files.append(
                        outdir
                        / template.format(
                            sample=clean_sample,
                            type=clean_type,
                        )
                    )
        else:
            output_files.append(outdir / template)


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
