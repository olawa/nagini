"""Sample routing: which analysis runs for which sample, and how it aligns.

Two facts decide what happens to a sample, and they are genuinely independent:

``library_type``
    What the molecules are and therefore how they must be aligned.  A cDNA
    library needs spliced alignment whatever it was made for.

``assay``
    Which analysis the sample was sequenced for.  A BCR-ABL1 panel needs
    force-called variants at known positions and a fusion breakdown; an
    adaptive-sampling run needs genome-wide calling restricted to a design.

Collapsing the two into one column is what made a cDNA amplicon impossible to
express: it is cDNA (so it must be spliced) *and* it needs variant calling,
and a single enum forces a choice between those.

Adding an assay is a config edit, not a code edit.  ``config["assays"]`` is
merged over the built-in definitions below, so a new panel is declared rather
than implemented::

    assays:
      my_panel:
        library_type: amplicon
        snv_calling: true
        sv_calling: false
        target_restricted: true
        modules: [my_panel]

This module deliberately imports nothing from Snakemake so that it can be
tested directly.
"""

from __future__ import annotations


# Alignment behaviour per library type.  ``spliced`` drives the minimap2
# preset; nothing else in the pipeline should branch on the library type.
LIBRARY_TYPES = {
    "wgs": {"spliced": False, "description": "Genomic DNA, whole genome"},
    "amplicon": {"spliced": False, "description": "Genomic DNA, amplicon or targeted"},
    "cdna": {"spliced": True, "description": "Reverse-transcribed RNA"},
    "direct_rna": {"spliced": True, "description": "Direct RNA sequencing"},
}

SPLICED_PRESET = "splice"

# Every key an assay definition may carry, with the value used when the
# definition leaves it out.
ASSAY_DEFAULTS = {
    "description": "",
    "library_type": "wgs",
    "snv_calling": False,
    "sv_calling": False,
    "cnv_calling": False,
    "target_restricted": False,
    "modules": [],
    # Tool settings that differ for this assay, as {section: {key: value}}.
    # A 20 000x amplicon and a 30x genome cannot share a minimum allele
    # fraction, so the assay carries the difference instead of the site config.
    "config_overrides": {},
}

BUILTIN_ASSAYS = {
    # The three names below match the historical library_type values so that a
    # units.tsv without an assay column keeps routing exactly as it did.
    "wgs": {
        "description": "Whole-genome ONT sequencing",
        "library_type": "wgs",
        "snv_calling": True,
        "sv_calling": True,
        "cnv_calling": True,
        "target_restricted": False,
    },
    "amplicon": {
        "description": "Genomic DNA amplicon or targeted panel",
        "library_type": "amplicon",
        "snv_calling": True,
        "sv_calling": False,
        "target_restricted": True,
        "modules": ["amplicon"],
        "config_overrides": {
            "clairs_to_call": {"snv_min_af": 0.01, "indel_min_af": 0.01},
            "read_filter": {"enabled": True},
        },
    },
    "flt3_tp53": {
        "description": "FLT3/TP53 AML targeted amplicon panel",
        "library_type": "amplicon",
        "snv_calling": False,
        "sv_calling": False,
        "target_restricted": True,
        "modules": ["amplicon"],
        "config_overrides": {
            "read_filter": {"enabled": True},
        },
    },
    "cdna": {
        "description": "cDNA library without an assay-specific analysis",
        "library_type": "cdna",
        "snv_calling": False,
        "sv_calling": False,
        "target_restricted": True,
    },
    "adaptive": {
        "description": "Adaptive sampling (ReadUntil) against a design BED",
        "library_type": "wgs",
        "snv_calling": True,
        "sv_calling": True,
        "cnv_calling": True,
        "target_restricted": True,
    },
    # Adaptive sampling read as a structural workup: what an aneuploidy or
    # rearrangement question needs (coverage, copy number, translocations)
    # without the small-variant callers, which dominate the runtime and answer
    # a question this workup is not asking.
    "adaptive_cnv_sv": {
        "description": "Adaptive sampling, structural workup: coverage, CNV, SV and translocations, no small variants",
        "library_type": "wgs",
        "snv_calling": False,
        "sv_calling": True,
        "cnv_calling": True,
        "target_restricted": True,
    },
    "bcr_abl": {
        "description": "BCR-ABL1 cDNA amplicon: fusion isoforms and kinase-domain resistance",
        "library_type": "cdna",
        "snv_calling": False,
        "sv_calling": False,
        "target_restricted": True,
        "modules": ["bcr_abl"],
        "config_overrides": {
            "read_filter": {
                "enabled": True,
            },
        },

    },

}


DEFAULT_ASSAY = "wgs"

_EMPTY = {"", "nan", "none", "na", "."}


class RoutingError(ValueError):
    """Raised for a sample sheet or config the pipeline cannot route."""


def normalize(value) -> str:
    """Lower-case and strip a cell, mapping the usual empty markers to ''."""

    text = str(value if value is not None else "").strip().lower()
    return "" if text in _EMPTY else text


def _as_bool(value, key: str, assay: str) -> bool:
    if isinstance(value, bool):
        return value
    text = normalize(value)
    if text in ("true", "yes", "1"):
        return True
    if text in ("false", "no", "0"):
        return False
    raise RoutingError(f"assays.{assay}.{key} must be true or false, got {value!r}.")


def resolve_assays(config=None) -> dict:
    """Merge ``config['assays']`` over the built-in definitions and validate."""

    resolved = {}
    for name, definition in BUILTIN_ASSAYS.items():
        merged = dict(ASSAY_DEFAULTS)
        merged.update(definition)
        merged["modules"] = list(merged["modules"])
        merged["config_overrides"] = {
            str(section): dict(values)
            for section, values in (merged.get("config_overrides") or {}).items()
        }
        resolved[name] = merged

    overrides = (config or {}).get("assays") or {}
    if not isinstance(overrides, dict):
        raise RoutingError("config['assays'] must be a mapping of assay name to definition.")

    for raw_name, definition in overrides.items():
        name = normalize(raw_name)
        if not name:
            raise RoutingError("An assay in config['assays'] has an empty name.")
        if definition is None:
            definition = {}
        if not isinstance(definition, dict):
            raise RoutingError(f"assays.{name} must be a mapping, got {type(definition).__name__}.")
        unknown = sorted(set(definition) - set(ASSAY_DEFAULTS))
        if unknown:
            raise RoutingError(
                f"assays.{name} has unsupported key(s): {', '.join(unknown)}. "
                f"Supported keys: {', '.join(sorted(ASSAY_DEFAULTS))}."
            )
        merged = dict(resolved.get(name) or ASSAY_DEFAULTS)
        merged.update(definition)
        merged["modules"] = list(merged.get("modules") or [])
        merged["config_overrides"] = dict(merged.get("config_overrides") or {})
        resolved[name] = merged

    for name, definition in resolved.items():
        library_type = normalize(definition["library_type"])
        if library_type not in LIBRARY_TYPES:
            raise RoutingError(
                f"assays.{name}.library_type={definition['library_type']!r} is not supported. "
                f"Valid library types: {', '.join(sorted(LIBRARY_TYPES))}."
            )
        definition["library_type"] = library_type
        definition["description"] = str(definition.get("description") or "")
        for key in ("snv_calling", "sv_calling", "target_restricted"):
            definition[key] = _as_bool(definition[key], key, name)
        definition["modules"] = [normalize(module) for module in definition["modules"] if normalize(module)]
        overrides_for_assay = definition.get("config_overrides") or {}
        if not isinstance(overrides_for_assay, dict):
            raise RoutingError(f"assays.{name}.config_overrides must be a mapping.")
        for section, values in overrides_for_assay.items():
            if not isinstance(values, dict):
                raise RoutingError(
                    f"assays.{name}.config_overrides.{section} must be a mapping of setting to value."
                )
        definition["config_overrides"] = {
            str(section): dict(values) for section, values in overrides_for_assay.items()
        }
    return resolved


def _column_values(units, sample: str, column: str) -> list:
    """Distinct non-empty values of a units column for one sample."""

    if units is None or column not in getattr(units, "columns", []):
        return []
    rows = units.loc[units["sample"].astype(str).str.strip() == str(sample).strip(), column]
    return sorted({normalize(value) for value in rows.tolist() if normalize(value)})


def get_assay(units, sample: str, assays: dict) -> str:
    """Return the assay for a sample.

    The ``assay`` column wins.  Without it the historical ``library_type``
    value is used when it names an assay, so an existing sample sheet keeps
    behaving the way it did.
    """

    values = _column_values(units, sample, "assay")
    if len(values) > 1:
        raise RoutingError(
            f"Sample {sample} has more than one assay in units.tsv: {', '.join(values)}. "
            "One sample belongs to one assay."
        )
    if values:
        assay = values[0]
        if assay not in assays:
            raise RoutingError(
                f"Sample {sample} has assay={assay!r}, which is not defined. "
                f"Known assays: {', '.join(sorted(assays))}. "
                "Declare it under config['assays'] to add one."
            )
        return assay

    legacy = _column_values(units, sample, "library_type")
    if len(legacy) > 1:
        raise RoutingError(
            f"Sample {sample} has more than one library_type in units.tsv: {', '.join(legacy)}."
        )
    if legacy and legacy[0] in assays:
        return legacy[0]
    if legacy and legacy[0] not in LIBRARY_TYPES:
        raise RoutingError(
            f"Sample {sample} has library_type={legacy[0]!r}, which is not supported. "
            f"Valid library types: {', '.join(sorted(LIBRARY_TYPES))}."
        )
    return DEFAULT_ASSAY


def get_library_type(units, sample: str, assays: dict) -> str:
    """Return the library type, preferring an explicit units.tsv value."""

    values = _column_values(units, sample, "library_type")
    if len(values) > 1:
        raise RoutingError(
            f"Sample {sample} has more than one library_type in units.tsv: {', '.join(values)}."
        )
    if values:
        if values[0] not in LIBRARY_TYPES:
            raise RoutingError(
                f"Sample {sample} has library_type={values[0]!r}, which is not supported. "
                f"Valid library types: {', '.join(sorted(LIBRARY_TYPES))}."
            )
        return values[0]
    return assays[get_assay(units, sample, assays)]["library_type"]


def get_assay_definition(units, sample: str, assays: dict) -> dict:
    return assays[get_assay(units, sample, assays)]


def sample_needs_snv_calling(units, sample: str, assays: dict) -> bool:
    return bool(get_assay_definition(units, sample, assays)["snv_calling"])


def sample_needs_sv_calling(units, sample: str, assays: dict) -> bool:
    return bool(get_assay_definition(units, sample, assays)["sv_calling"])


def sample_needs_cnv_calling(units, sample: str, assays: dict) -> bool:
    return bool(get_assay_definition(units, sample, assays).get("cnv_calling", False))


def sample_is_target_restricted(units, sample: str, assays: dict) -> bool:
    return bool(get_assay_definition(units, sample, assays)["target_restricted"])


def get_sample_caller(units, sample: str, assays: dict, global_caller: str = "auto", unit_type: str = "") -> str:
    """Determine the small variant caller for a given sample and unit type.

    Precedence:
    1. Assay config override: config_overrides.<section>.small_variant_caller or caller
    2. Global config setting: if not 'auto', uses the specified caller
    3. Auto resolution:
       - 'N' (Normal/Germline) -> 'deepvariant'
       - 'T' (Tumor) -> 'clairs_to'
       - amplicon / cdna -> 'clairs_to'
       - wgs -> 'deepvariant'
    """
    overrides = get_assay_definition(units, sample, assays).get("config_overrides", {})
    for sec_data in overrides.values():
        if isinstance(sec_data, dict):
            if "small_variant_caller" in sec_data:
                return normalize(sec_data["small_variant_caller"])
            if "caller" in sec_data:
                return normalize(sec_data["caller"])
    norm_global = normalize(global_caller or "auto")
    if norm_global != "auto":
        return norm_global
    lib = get_library_type(units, sample, assays)
    if lib in ("amplicon", "cdna"):
        return "clairs_to"
    if str(unit_type).upper() == "N":
        return "deepvariant"
    if str(unit_type).upper() == "T":
        return "clairs_to"
    return "deepvariant"





def sample_modules(units, sample: str, assays: dict) -> list:
    return list(get_assay_definition(units, sample, assays)["modules"])


def sample_has_module(units, sample: str, assays: dict, module: str) -> bool:
    return normalize(module) in sample_modules(units, sample, assays)


def library_type_is_spliced(library_type: str) -> bool:
    return bool(LIBRARY_TYPES.get(normalize(library_type), {}).get("spliced", False))


def minimap2_preset(library_type: str, default_preset: str = "map-ont") -> str:
    """Spliced libraries always align spliced, whatever the assay wants."""

    return SPLICED_PRESET if library_type_is_spliced(library_type) else default_preset


def routing_warnings(units, samples, assays: dict) -> list:
    """Non-fatal inconsistencies worth showing at parse time."""

    warnings = []
    for sample in samples:
        assay = get_assay(units, sample, assays)
        declared = _column_values(units, sample, "library_type")
        expected = assays[assay]["library_type"]
        if declared and declared[0] != expected:
            warnings.append(
                f"Sample {sample}: units.tsv says library_type={declared[0]}, but assay "
                f"{assay} expects {expected}. Using {declared[0]}; alignment will follow it."
            )
    return warnings


def routing_summary(units, samples, assays: dict) -> list:
    """One row per sample: what will run and why.  Used for logging."""

    summary = []
    for sample in samples:
        assay = get_assay(units, sample, assays)
        library_type = get_library_type(units, sample, assays)
        definition = assays[assay]
        summary.append(
            {
                "sample": str(sample).strip(),
                "assay": assay,
                "library_type": library_type,
                "preset": minimap2_preset(library_type),
                "snv_calling": bool(definition["snv_calling"]),
                "sv_calling": bool(definition["sv_calling"]),
                "target_restricted": bool(definition["target_restricted"]),
                "modules": list(definition["modules"]),
            }
        )
    return summary


# ---------------------------------------------------------------------------
# Conditions usable in `when:` in config/output_files.yaml
# ---------------------------------------------------------------------------
#
# Declared here so that a typo in the output spec fails at parse time instead
# of silently dropping the output: an unknown `when:` used to fall through to a
# config flag lookup, which is False for anything that does not exist.

PIPELINE_CONDITION_NAMES = {
    "annotate_small_variants_with_vep",
    "clairs_to_enabled",
    "clairs_to_vep_enabled",
    "deepvariant_enabled",
    "deepvariant_vep_enabled",
    "deepsomatic_enabled",
    "deepsomatic_vep_enabled",
    "mutect2_enabled",
    "mutect2_vep_enabled",
    "annotate_small_variants",
    "annotate_svs",
    "has_design_bed",
    "target_metrics",
    "merged_fastq",
    "cnv_calling",
    "cnvkit_enabled",
    "sniffles2_enabled",
    "savana_enabled",
    "rindels_sv",
}

SAMPLE_CONDITION_NAMES = {
    "clairs_to_enabled",
    "clairs_to_vep_enabled",
    "deepvariant_enabled",
    "deepvariant_vep_enabled",
    "deepsomatic_enabled",
    "deepsomatic_vep_enabled",
    "mutect2_enabled",
    "mutect2_vep_enabled",
    "annotate_small_variants",
    "annotate_small_variants_with_vep",
    "sample_needs_dna_variant_calling",
    "sv_calling",
    "annotate_svs",
    "sample_needs_sv_calling",
    "sample_needs_rna_analysis",
    "target_restricted",
    "target_metrics",
    "merged_fastq",
    "cnv_calling",
    "sample_needs_cnv_calling",
    "cnvkit_enabled",
    "sniffles2_enabled",
    "savana_enabled",
    "rindels_sv",
}

# Config flags that are legitimate `when:` values even though they are plain
# booleans rather than named predicates.
CONFIG_FLAG_CONDITIONS = {
    "annotate_small_variants",
    "annotate_svs",
    "sv_calling",
    "cnv_calling",
    "enable_rindels",
    "enable_bcftools_pileup",
}


CONDITION_PREFIXES = ("assay:", "assay_module:")


def is_known_condition(condition, assays: dict | None = None) -> bool:
    """Whether a `when:` value names something the pipeline can evaluate."""

    if not condition:
        return True
    for prefix in CONDITION_PREFIXES:
        if condition.startswith(prefix):
            value = normalize(condition[len(prefix):])
            if not value:
                return False
            if prefix == "assay:" and assays is not None:
                return value in assays
            return True
    return (
        condition in PIPELINE_CONDITION_NAMES
        or condition in SAMPLE_CONDITION_NAMES
        or condition in CONFIG_FLAG_CONDITIONS
    )


def unknown_conditions(output_spec: dict, assays: dict | None = None) -> list:
    """Every `when:` in an output spec that the pipeline cannot evaluate."""

    return sorted(
        {
            str(entry.get("when"))
            for entry in (output_spec or {}).get("files", [])
            if entry.get("when") and not is_known_condition(str(entry["when"]), assays)
        }
    )


def assay_setting(units, sample: str, assays: dict, section: str, key: str, default=None):
    """A tool setting for one sample, with the assay's override applied.

    Precedence is assay override, then the value the caller passes as
    ``default`` (normally the plain config value).
    """

    overrides = get_assay_definition(units, sample, assays).get("config_overrides", {})
    return overrides.get(section, {}).get(key, default)


# ---------------------------------------------------------------------------
# Basecalling model -> ClairS-TO platform
# ---------------------------------------------------------------------------
#
# ClairS-TO ships one error model per basecaller configuration.  Running SUP
# 5 kHz data against a HAC 4 kHz model is a silent mismatch, so the platform is
# derived from the basecalling_model recorded in units.tsv rather than pinned
# in the config.  The mapping is overridable from the config because ONT adds
# model names faster than a pipeline gets released.

CLAIRS_TO_DEFAULT_PLATFORM = "ont_r10_dorado_sup_5khz"

CLAIRS_TO_PLATFORMS = {
    ("r10", "sup", "5khz"): "ont_r10_dorado_sup_5khz",
    ("r10", "hac", "5khz"): "ont_r10_dorado_hac_5khz",
    ("r10", "sup", "4khz"): "ont_r10_dorado_sup_4khz",
    ("r10", "hac", "4khz"): "ont_r10_dorado_hac_4khz",
}

# ONT moved to 5 kHz sampling with the v4.2.0 basecalling models; anything
# older is 4 kHz.  An explicit 4khz/5khz token in the model name always wins.
_FIVE_KHZ_FROM_VERSION = (4, 2)


def parse_basecalling_model(model: str) -> dict:
    """Pull chemistry, accuracy and sampling rate out of a Dorado model name."""

    text = normalize(model)
    result = {"chemistry": "", "accuracy": "", "sampling": "", "version": ""}
    if not text:
        return result

    if "r10" in text:
        result["chemistry"] = "r10"
    elif "r9" in text:
        result["chemistry"] = "r9"

    for accuracy in ("sup", "hac", "fast"):
        if accuracy in text.replace("@", "_").split("_") or f"_{accuracy}" in text:
            result["accuracy"] = accuracy
            break

    if "5khz" in text:
        result["sampling"] = "5khz"
    elif "4khz" in text:
        result["sampling"] = "4khz"

    version = ""
    if "@" in text:
        version = text.rsplit("@", 1)[1].lstrip("v")
    result["version"] = version

    if not result["sampling"] and version:
        parts = []
        for piece in version.split("."):
            digits = "".join(character for character in piece if character.isdigit())
            parts.append(int(digits) if digits else 0)
        while len(parts) < 2:
            parts.append(0)
        result["sampling"] = (
            "5khz" if tuple(parts[:2]) >= _FIVE_KHZ_FROM_VERSION else "4khz"
        )
    return result


def clairs_to_platform(model: str, table: dict | None = None, default: str | None = None):
    """Return ``(platform, warning)`` for a basecalling model.

    An unrecognised model is not an error: ONT names models faster than this
    pipeline is released, and refusing to run would be worse than running with
    a stated assumption.  The caller is expected to surface the warning.
    """

    fallback = default or CLAIRS_TO_DEFAULT_PLATFORM
    lookup = dict(CLAIRS_TO_PLATFORMS)
    for key, value in (table or {}).items():
        if isinstance(key, str):
            parts = tuple(normalize(part) for part in key.split(","))
            if len(parts) == 3:
                lookup[parts] = value
        elif isinstance(key, tuple) and len(key) == 3:
            lookup[tuple(normalize(part) for part in key)] = value

    parsed = parse_basecalling_model(model)
    signature = (parsed["chemistry"], parsed["accuracy"], parsed["sampling"])
    if all(signature) and signature in lookup:
        return lookup[signature], None

    return fallback, (
        f"Basecalling model {model!r} does not map to a ClairS-TO platform "
        f"(read as chemistry={parsed['chemistry'] or '?'}, "
        f"accuracy={parsed['accuracy'] or '?'}, "
        f"sampling={parsed['sampling'] or '?'}). Using {fallback}."
    )


def get_basecalling_model(units, sample: str) -> str:
    """The basecalling model recorded for a sample, empty when absent."""

    values = _column_values(units, sample, "basecalling_model")
    return values[0] if values else ""


# ---------------------------------------------------------------------------
# Barcode kits
# ---------------------------------------------------------------------------
#
# A barcode listed in the sample sheet that the kit cannot produce is a
# transcription slip, and it is cheap to catch before a basecalling job is
# queued rather than after it.

# Kits whose barcode count is not the trailing number in their name.
BARCODE_KIT_SIZES = {
    "EXP-NBD103": 12,
    "EXP-NBD104": 12,
    "EXP-NBD114": 24,
    "EXP-PBC001": 12,
    "EXP-PBC096": 96,
    "SQK-16S024": 24,
    "SQK-LWB001": 12,
    "SQK-PBK004": 12,
    "SQK-RAB204": 12,
    "SQK-RBK004": 12,
    "SQK-RLB001": 12,
}


def barcode_kit_size(kit: str, table: dict | None = None):
    """How many barcodes a kit provides, or None when it cannot be determined.

    Most ONT kit names end in their barcode count (SQK-NBD114-24), so the name
    is read first and the table covers the ones that do not follow it.
    """

    name = str(kit or "").strip().upper()
    if not name:
        return None
    lookup = {str(key).strip().upper(): value for key, value in (table or {}).items()}
    lookup_all = {**BARCODE_KIT_SIZES, **lookup}
    if name in lookup_all:
        return int(lookup_all[name])
    tail = name.rsplit("-", 1)[-1]
    if tail.isdigit():
        size = int(tail)
        # Guard against reading a chemistry version (RBK114) as a count.
        return size if 1 <= size <= 384 else None
    return None


def barcodes_outside_kit(barcodes, kit: str, table: dict | None = None) -> list:
    """Barcodes that the kit cannot produce, as ``[(barcode, number)]``."""

    size = barcode_kit_size(kit, table)
    if size is None:
        return []
    outside = []
    for barcode in barcodes:
        digits = "".join(character for character in str(barcode) if character.isdigit())
        if not digits:
            continue
        number = int(digits)
        if number < 1 or number > size:
            outside.append((str(barcode), number))
    return outside
