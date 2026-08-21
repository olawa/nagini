#!/usr/bin/env python3
"""Create Hydra-style input tables and a run-local Nagini config.

Hydra-Genetics can create ONT tables from BAM read groups. FASTQ and MinKNOW
input do not have enough metadata for Hydra to infer the same table, so this
helper creates the equivalent table using conservative defaults and explicit
CLI options. Dorado's ``bam_pass/barcodeXX`` layout can be ingested directly
when ``--barcode-map`` supplies the biological sample name for each barcode.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

import yaml


FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
BASE_UNIT_COLUMNS = [
    "sample",
    "type",
    "platform",
    "machine",
    "processing_unit",
    "barcode",
    "barcode_name",
    "input_type",
    "fastq1",
    "fastq2",
    "bam",
    "methylation",
    "run_id",
    "basecalling_model",
]
DEFAULT_SAMPLE_REGEX = r"(?P<sample>[^_]+)(?:_.*)?\.(?:fastq|fq)(?:\.gz)?$"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument("--pipeline-root", type=Path, required=True)
    parser.add_argument("--base-config", type=Path, required=True)
    parser.add_argument("--site-config", type=Path)
    parser.add_argument("--reference-files", type=Path)
    parser.add_argument("--library-root", type=Path)
    parser.add_argument(
        "--input-type", choices=("minknow", "fastq", "unmapped_bam", "aligned_bam"), required=True
    )
    parser.add_argument("--sample-type", choices=("N", "T", "R"), default="N")
    parser.add_argument("--sequence-id", default="")
    parser.add_argument("--sample-regex", default=DEFAULT_SAMPLE_REGEX)
    parser.add_argument(
        "--barcode-map",
        type=Path,
        help="CSV/TSV mapping barcode or barcodeXX to sample (optional for Dorado BAM input).",
    )
    parser.add_argument("--barcode", default="NB")
    parser.add_argument("--run-id", default="")
    parser.add_argument("--basecalling-model", default="unknown")
    parser.add_argument("--dorado-input", type=Path)
    parser.add_argument("--dorado-model", default="")
    parser.add_argument("--dorado-kit", default="")
    parser.add_argument("--dorado-device", default="")
    parser.add_argument("--dorado-models-directory", type=Path)
    parser.add_argument("--dorado-container", default="")
    parser.add_argument("--dorado-binary", default="")
    parser.add_argument(
        "--dorado-demux-output-format",
        choices=("bam", "fastq", "both"),
        default="fastq",
    )
    parser.add_argument(
        "--dorado-aligner", choices=("minimap2", "dorado"), default="minimap2"
    )
    parser.add_argument("--methylation", choices=("Yes", "No"), default="No")
    parser.add_argument(
        "--bed-file",
        "--bed",
        "--design-bed",
        dest="bed_file",
        type=Path,
        help="Target region BED file (overrides site default). If omitted, any single *.bed in inbox is auto-detected.",
    )
    parser.add_argument(
        "--small-variant-caller",
        choices=("clairs_to", "mutect2", "both"),
        default="",
        help="Somatic small-variant caller to execute and report (default: clairs_to from config).",
    )
    parser.add_argument(
        "--clairs-to-platform",
        default="",
        help="ClairS-TO platform model (e.g. ont_r10_dorado_hac_4khz). Auto-derived from basecalling model if omitted.",
    )
    parser.add_argument(
        "--clairs-to-container",
        default="",
        help="Override default ClairS-TO container image.",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=None,
        help="Number of CPU cores allocated for the pipeline run.",
    )
    parser.add_argument("--tumor-content", default=None)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


def normalize_barcode(value: object) -> str:
    """Canonicalize barcode05, bc05, barcode5 and 05 to the same key."""

    barcode = str(value or "").strip().lower()
    match = re.fullmatch(r"(?:barcode|bc)?(\d+)", barcode)
    if match:
        return f"barcode{int(match.group(1)):02d}"
    return barcode


def read_barcode_map(path: Path, default_type: str) -> dict[str, dict[str, str]]:
    """Read a barcode-to-sample CSV/TSV mapping with strict unique barcodes."""

    if not path.is_file():
        raise SystemExit(f"Barcode map does not exist: {path}")
    with path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
        except csv.Error:
            dialect = csv.excel_tab if path.suffix.lower() in {".tsv", ".tab"} else csv.excel
        rows = list(csv.DictReader(handle, dialect=dialect))

    if not rows:
        raise SystemExit(f"Barcode map is empty: {path}")
    normalized_rows = []
    for row in rows:
        fields = {
            re.sub(r"[^a-z0-9]", "", str(key or "").lower()): str(value or "").strip()
            for key, value in row.items()
        }
        barcode = normalize_barcode(fields.get("barcode", fields.get("barcodename", "")))
        sample_name = fields.get("sample", fields.get("sampleid", ""))
        sample_type = fields.get("type", fields.get("sampletype", "")) or default_type
        if not barcode or not re.fullmatch(r"barcode\d+", barcode):
            raise SystemExit(f"Invalid or missing barcode in {path}: {row}")
        if not sample_name:
            raise SystemExit(f"Missing sample/sample_id for {barcode} in {path}.")
        if sample_type not in {"N", "T", "R"}:
            raise SystemExit(
                f"Invalid type={sample_type!r} for {barcode}; expected N, T or R."
            )
        normalized_rows.append(
            {"barcode": barcode, "sample": sample_name, "type": sample_type}
        )

    mapping: dict[str, dict[str, str]] = {}
    for row in normalized_rows:
        if row["barcode"] in mapping:
            raise SystemExit(f"Barcode appears more than once in {path}: {row['barcode']}")
        mapping[row["barcode"]] = row
    return mapping


def discover_bams_by_barcode(input_dir: Path) -> dict[str, list[Path]]:
    """Find BAMs below barcodeXX directories or with barcodeXX filenames."""

    bams: defaultdict[str, list[Path]] = defaultdict(list)
    root_barcode = normalize_barcode(input_dir.name)
    if not re.fullmatch(r"barcode\d+", root_barcode):
        root_barcode = ""
    for path in sorted(input_dir.rglob("*.bam")):
        relative_parts = path.relative_to(input_dir).parts[:-1]
        barcode = ""
        for part in reversed(relative_parts):
            candidate = normalize_barcode(part)
            if re.fullmatch(r"barcode\d+", candidate):
                barcode = candidate
                break
        if not barcode:
            # Dorado demux versions may emit ``barcode01.bam`` directly below
            # the output directory instead of creating barcode subdirectories.
            for token in re.findall(r"(?i)(?:barcode|bc)\d+", path.stem):
                candidate = normalize_barcode(token)
                if re.fullmatch(r"barcode\d+", candidate):
                    barcode = candidate
                    break
        barcode = barcode or root_barcode
        if barcode:
            bams[barcode].append(path.resolve())
    return dict(bams)


def make_bam_tables_from_barcode_map(
    args: argparse.Namespace,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build Hydra-style tables from Dorado ``bam_pass/barcodeXX`` output."""

    mapping = read_barcode_map(args.barcode_map.resolve(), args.sample_type)
    bams_by_barcode = discover_bams_by_barcode(args.input_dir)
    if not bams_by_barcode:
        raise SystemExit(f"No BAM files below barcode directories in {args.input_dir}.")

    observed = set(bams_by_barcode)
    expected = set(mapping)
    missing = sorted(expected - observed)
    unmapped = sorted(observed - expected)
    if missing:
        raise SystemExit(
            "Barcode map contains barcodes with no BAM files: " + ", ".join(missing)
        )
    if unmapped:
        raise SystemExit(
            "BAM barcode directories are absent from the barcode map: "
            + ", ".join(unmapped)
            + ". Add them to the map or remove them from the input directory."
        )

    run_id = args.run_id or args.sequence_id or args.input_dir.name
    samples_by_name: dict[str, dict[str, str]] = {}
    units: list[dict[str, str]] = []
    for barcode in sorted(mapping):
        assignment = mapping[barcode]
        sample_name = assignment["sample"]
        samples_by_name.setdefault(sample_name, {"sample": sample_name})
        for index, bam_path in enumerate(sorted(bams_by_barcode[barcode]), start=1):
            units.append(
                {
                    "sample": sample_name,
                    "type": assignment["type"],
                    "platform": "ONT",
                    "machine": "unknown",
                "processing_unit": f"{barcode}_{index:03d}",
                # Hydra's ONT schema accepts flowcell-style values such as
                # BC/NB but rejects Dorado's numeric barcode05 identifiers.
                # Keep the original Dorado identifier in barcode_name and
                # use a valid technical barcode value for Hydra.
                "barcode": "BC",
                "barcode_name": barcode,
                    "input_type": args.input_type,
                    "fastq1": "UNUSED",
                    "fastq2": "UNUSED",
                    "bam": str(bam_path),
                    "methylation": args.methylation,
                    "run_id": run_id,
                    "basecalling_model": args.basecalling_model,
                }
            )
    return sorted(samples_by_name.values(), key=lambda row: row["sample"]), units


def read_yaml(path: Path) -> dict[str, Any]:
    with path.open() as handle:
        value = yaml.safe_load(handle) or {}
    if not isinstance(value, dict):
        raise SystemExit(f"Expected a YAML mapping in {path}.")
    return value


def deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Merge site values over the shared config without mutating either input."""

    merged: dict[str, Any] = dict(base)
    for key, value in override.items():
        if isinstance(merged.get(key), dict) and isinstance(value, dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


TEMPLATE = re.compile(r"\{\{\s*([A-Za-z_][A-Za-z0-9_]*)\s*\}\}")


def expand_templates(value: Any, context: dict[str, str]) -> Any:
    """Expand Pickett-style ``{{VARIABLE}}`` values recursively."""

    if isinstance(value, dict):
        return {key: expand_templates(item, context) for key, item in value.items()}
    if isinstance(value, list):
        return [expand_templates(item, context) for item in value]
    if not isinstance(value, str):
        return value

    expanded = value
    for _ in range(12):
        updated = TEMPLATE.sub(
            lambda match: context.get(match.group(1), match.group(0)), expanded
        )
        if updated == expanded:
            break
        expanded = updated
    unresolved = sorted(set(TEMPLATE.findall(expanded)))
    if unresolved:
        names = ", ".join(unresolved)
        raise SystemExit(f"Unresolved config variable(s) in {value!r}: {names}")
    return expanded


def write_tsv(path: Path, rows: list[dict[str, str]], columns: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows({column: row.get(column, "") for column in columns} for row in rows)


def sample_from_fastq(path: Path, pattern: re.Pattern[str]) -> str:
    match = pattern.search(path.name)
    if not match:
        raise SystemExit(
            f"Cannot extract a sample name from {path.name!r} with --sample-regex."
        )
    if "sample" in match.groupdict():
        sample = match.group("sample")
    elif match.groups():
        sample = match.group(1)
    else:
        sample = match.group(0)
    if not sample:
        raise SystemExit(f"The sample name extracted from {path.name!r} is empty.")
    return sample


def make_fastq_tables(args: argparse.Namespace) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    files = sorted(
        path.resolve()
        for path in args.input_dir.rglob("*")
        if path.is_file() and path.name.lower().endswith(FASTQ_SUFFIXES)
    )
    if not files:
        raise SystemExit(f"No FASTQ files found below {args.input_dir}.")

    pattern = re.compile(args.sample_regex)
    counters: defaultdict[str, int] = defaultdict(int)
    units: list[dict[str, str]] = []
    for path in files:
        sample = sample_from_fastq(path, pattern)
        counters[sample] += 1
        run_id = args.run_id or args.sequence_id or args.input_dir.name
        units.append(
            {
                "sample": sample,
                "type": args.sample_type,
                "platform": "ONT",
                "machine": "unknown",
                "processing_unit": f"PU{counters[sample]}",
                "barcode": args.barcode,
                "input_type": "fastq",
                "fastq1": str(path),
                "fastq2": "UNUSED",
                "bam": "UNUSED",
                "methylation": args.methylation,
                "run_id": run_id,
                "basecalling_model": args.basecalling_model,
            }
        )
    samples = [{"sample": sample} for sample in sorted(counters)]
    return samples, units


def make_minknow_tables(args: argparse.Namespace) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Create deferred units for Dorado outputs produced by Snakemake."""

    mapping = read_barcode_map(args.barcode_map.resolve(), args.sample_type)
    run_id = args.run_id or args.sequence_id or args.input_dir.name
    samples_by_name: dict[str, dict[str, str]] = {}
    units: list[dict[str, str]] = []
    for barcode in sorted(mapping):
        assignment = mapping[barcode]
        sample_name = assignment["sample"]
        samples_by_name.setdefault(sample_name, {"sample": sample_name})
        units.append(
            {
                "sample": sample_name,
                "type": assignment["type"],
                "platform": "ONT",
                "machine": "unknown",
                "processing_unit": f"{barcode}_001",
                "barcode": "BC",
                "barcode_name": barcode,
                "input_type": "dorado_bam",
                "fastq1": "UNUSED",
                "fastq2": "UNUSED",
                "bam": f"basecalling/dorado_demux/{barcode}.bam",
                "methylation": args.methylation,
                "run_id": run_id,
                "basecalling_model": args.basecalling_model,
            }
        )
    return sorted(samples_by_name.values(), key=lambda row: row["sample"]), units


def make_bam_tables(args: argparse.Namespace) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    files = sorted(
        path.resolve()
        for path in args.input_dir.rglob("*")
        if path.is_file() and path.name.lower().endswith(".bam")
    )
    if not files:
        raise SystemExit(f"No BAM files found below {args.input_dir}.")

    pattern = re.compile(args.sample_regex) if args.sample_regex else None
    counters: defaultdict[str, int] = defaultdict(int)
    units: list[dict[str, str]] = []
    samples_by_name: dict[str, dict[str, str]] = {}

    for path in files:
        sample = ""
        sample_type = args.sample_type
        if pattern:
            match = pattern.search(path.name)
            if match:
                if "sample" in match.groupdict():
                    sample = match.group("sample")
                elif match.groups():
                    sample = match.group(1)
                else:
                    sample = match.group(0)
                sample_type = match.groupdict().get("type", args.sample_type)
        if not sample:
            stem = path.stem
            if stem.endswith(("_T", "_N", "_R")):
                sample, sample_type = stem[:-2], stem[-1]
            else:
                sample = stem

        sample = str(sample).strip()
        sample_type = str(sample_type).strip()
        counters[sample] += 1
        run_id = args.run_id or args.sequence_id or args.input_dir.name
        samples_by_name.setdefault(sample, {"sample": sample})
        units.append(
            {
                "sample": sample,
                "type": sample_type,
                "platform": "ONT",
                "machine": "unknown",
                "processing_unit": f"PU{counters[sample]:03d}",
                "barcode": args.barcode,
                "input_type": args.input_type,
                "fastq1": "UNUSED",
                "fastq2": "UNUSED",
                "bam": str(path),
                "methylation": args.methylation,
                "run_id": run_id,
                "basecalling_model": args.basecalling_model,
            }
        )
    return sorted(samples_by_name.values(), key=lambda row: row["sample"]), units


def validate_tables(
    samples: list[dict[str, str]],
    units: list[dict[str, str]],
    allow_deferred_inputs: bool = False,
) -> None:
    sample_ids = {row.get("sample", "") for row in samples}
    if not sample_ids or "" in sample_ids:
        raise SystemExit("samples.tsv must contain non-empty sample identifiers.")
    seen: set[tuple[str, str, str, str]] = set()
    for row in units:
        key = tuple(row.get(column, "") for column in ("sample", "type", "processing_unit", "barcode"))
        if "" in key:
            raise SystemExit(f"Incomplete unit metadata: {row}")
        if key in seen:
            raise SystemExit(f"Duplicate unit metadata: {key}")
        seen.add(key)
        if row["sample"] not in sample_ids:
            raise SystemExit(f"Unit refers to sample absent from samples.tsv: {row['sample']}")
        input_type = row.get("input_type")
        input_path = row.get("fastq1") if input_type == "fastq" else row.get("bam")
        if (
            not input_path
            or input_path == "UNUSED"
            or (not allow_deferred_inputs and not Path(input_path).exists())
        ):
            raise SystemExit(f"Input file does not exist for unit {key}: {input_path}")


def resolve_config_path(value: object, base_dir: Path) -> object:
    if not isinstance(value, str) or not value or value.startswith(("docker://", "http://", "https://")):
        return value
    path = Path(value).expanduser()
    return str(path if path.is_absolute() else (base_dir / path).resolve())


def make_runtime_config(args: argparse.Namespace, analysis_dir: Path) -> Path:
    base_config = args.base_config.resolve()
    config = read_yaml(base_config)
    if args.site_config:
        site_config = args.site_config.resolve()
        if not site_config.is_file():
            raise SystemExit(f"Site config does not exist: {site_config}")
        config = deep_merge(config, read_yaml(site_config))
    reference_registry = None
    if args.reference_files:
        reference_path = args.reference_files.resolve()
        if not reference_path.is_file():
            raise SystemExit(f"Reference registry does not exist: {reference_path}")
        reference_registry = read_yaml(reference_path)

    context: dict[str, str] = {
        str(key): str(value)
        for key, value in config.items()
        if isinstance(value, (str, int, float, bool))
    }
    context.update(
        {
            "PIPELINE_ROOT": str(args.pipeline_root.resolve()),
            "ANALYSIS_DIR": str(analysis_dir.resolve()),
        }
    )
    # Site values can themselves refer to another site variable. Resolve the
    # scalar context first, then expand the complete nested configuration.
    for _ in range(12):
        changed = False
        for key, value in list(context.items()):
            expanded = expand_templates(value, context)
            if expanded != value:
                context[key] = expanded
                changed = True
        if not changed:
            break
    config = expand_templates(config, context)
    if reference_registry:
        registry_site = reference_registry.get("site")
        if registry_site and registry_site != config.get("site"):
            raise SystemExit(
                "Reference registry site does not match selected site: "
                f"{registry_site!r} != {config.get('site')!r}"
            )

    # Nagini's checked-in configs use paths relative to the repository root,
    # even though the YAML files themselves live below config/.
    base_dir = args.pipeline_root.resolve()
    config["samples"] = str((analysis_dir / "samples.tsv").resolve())
    config["units"] = str((analysis_dir / "units.tsv").resolve())
    config["output"] = str(Path(resolve_config_path(config.get("output", "config/output_files.yaml"), base_dir)))
    config["resources"] = str(Path(resolve_config_path(config.get("resources", "config/resources.yaml"), base_dir)))
    config["hydra_local_path"] = str(Path(resolve_config_path(config.get("hydra_local_path", "vendor"), base_dir)))
    if config.get("container_cache"):
        config["container_cache"] = str(resolve_config_path(config["container_cache"], base_dir))

    for key in (
        "REFERENCE_ROOT",
        "CONTAINER_CACHE",
        "HG38_FASTA",
        "HG38_FAI",
        "HG38_DICT",
        "HG38_MMI",
        "HG38_DESIGN_BED",
        "HG38_SV_ANNOTATION_BED",
        "HG38_TANDEM_REPEATS",
        "HG38_CLINVAR",
        "HG38_CYTOBAND",
        "HG38_GTF",
        "HG38_ANNOTATIONS",
        "HS1_FASTA",
        "HS1_FAI",
        "HS1_DICT",
        "HS1_MMI",
        "HS1_DESIGN_BED",
        "HS1_SV_ANNOTATION_BED",
        "HS1_TANDEM_REPEATS",
        "HS1_CLINVAR",
        "HS1_CYTOBAND",
        "HS1_GTF",
        "HS1_ANNOTATIONS",
        "VEP_CACHE",
        "SNPEFF_DB",
        "SIMPLE_SV_PANEL",
        "SIMPLE_SV_FUSION_PAIRS",
    ):
        if config.get(key):
            config[key] = str(resolve_config_path(config[key], base_dir))

    for reference in config.get("references", {}).values():
        if not isinstance(reference, dict):
            continue
        for key in (
            "fasta",
            "fai",
            "dict",
            "mmi",
            "design_bed",
            "sv_annotation_bed",
            "tandem_repeats",
            "clinvar",
            "cytoband",
            "gene_annotation",
            "annotations",
        ):
            if key in reference and reference[key]:
                reference[key] = resolve_config_path(reference[key], base_dir)
    for section, keys in {
        "vep": ("vep_cache",),
        "snpeff": ("db",),
        "simple_sv_annotation": ("panel", "fusion_pairs"),
    }.items():
        for key in keys:
            if section in config and config[section].get(key):
                config[section][key] = resolve_config_path(config[section][key], base_dir)

    config["site_config_path"] = str(args.site_config.resolve()) if args.site_config else ""
    config["reference_files_path"] = (
        str(args.reference_files.resolve()) if args.reference_files else ""
    )
    config["library_root"] = str(args.library_root.resolve()) if args.library_root else ""

    if getattr(args, "input_type", "") == "minknow":
        dorado_settings = {
            "input": str((args.dorado_input or args.input_dir).resolve()),
            "model": args.dorado_model or args.basecalling_model,
            "kit": args.dorado_kit,
            "device": args.dorado_device,
            "models_directory": str(args.dorado_models_directory.resolve())
            if args.dorado_models_directory
            else "",
            "demux_output_format": getattr(
                args, "dorado_demux_output_format", "fastq"
            ),
        }
        if getattr(args, "dorado_container", ""):
            dorado_settings["container"] = args.dorado_container
        if getattr(args, "dorado_binary", ""):
            dorado_settings["binary"] = str(Path(args.dorado_binary).expanduser().resolve())
            if not getattr(args, "dorado_container", ""):
                dorado_settings["container"] = "none"
        config.setdefault("dorado", {}).update(dorado_settings)
        config.setdefault("dorado_alignment", {})["aligner"] = getattr(
            args, "dorado_aligner", "minimap2"
        )

    bed_path = None
    if getattr(args, "bed_file", None):
        resolved_bed = args.bed_file.resolve()
        if not resolved_bed.is_file():
            raise SystemExit(f"Target BED file does not exist: {resolved_bed}")
        bed_path = str(resolved_bed)
    elif getattr(args, "input_dir", None) and args.input_dir.is_dir():
        bed_files = sorted(args.input_dir.glob("*.bed"))
        if len(bed_files) == 1:
            bed_path = str(bed_files[0].resolve())

    if bed_path:
        ref_build = config.get("reference_build", "hg38")
        if "references" in config and isinstance(config["references"], dict):
            if ref_build in config["references"] and isinstance(config["references"][ref_build], dict):
                config["references"][ref_build]["design_bed"] = bed_path
        if "reference" in config and isinstance(config["reference"], dict):
            config["reference"]["design_bed"] = bed_path
        config["HG38_DESIGN_BED"] = bed_path
        config["HS1_DESIGN_BED"] = bed_path

    if getattr(args, "small_variant_caller", ""):
        config["small_variant_caller"] = args.small_variant_caller

    clairs_to_platform = getattr(args, "clairs_to_platform", "")
    if not clairs_to_platform:
        model_str = (
            getattr(args, "dorado_model", "")
            or getattr(args, "basecalling_model", "")
        ).lower()
        if "sup" in model_str:
            clairs_to_platform = "ont_r10_dorado_sup_5khz_ssrs"
        elif "hac" in model_str or "fast" in model_str:
            clairs_to_platform = "ont_r10_dorado_hac_4khz"

    if clairs_to_platform:
        config.setdefault("clairs_to_call", {})["platform"] = clairs_to_platform

    if getattr(args, "clairs_to_container", ""):
        config.setdefault("clairs_to_call", {})["container"] = args.clairs_to_container

    if getattr(args, "cores", None):
        config["cores"] = int(args.cores)

    runtime_config = analysis_dir / "config.nagini.yaml"
    with runtime_config.open("w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return runtime_config


def main() -> None:
    args = parse_args()
    args.input_dir = args.input_dir.resolve()
    args.analysis_dir = args.analysis_dir.resolve()
    args.analysis_dir.mkdir(parents=True, exist_ok=True)
    if not args.input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist: {args.input_dir}")

    samples_path = args.analysis_dir / "samples.tsv"
    units_path = args.analysis_dir / "units.tsv"
    if args.input_type == "fastq":
        if samples_path.exists() and units_path.exists() and not args.force:
            samples, units = read_tsv(samples_path), read_tsv(units_path)
        else:
            samples, units = make_fastq_tables(args)
    elif args.input_type == "minknow":
        if samples_path.exists() and units_path.exists() and not args.force:
            samples, units = read_tsv(samples_path), read_tsv(units_path)
        else:
            samples, units = make_minknow_tables(args)
    elif args.barcode_map:
        if samples_path.exists() and units_path.exists() and not args.force:
            samples, units = read_tsv(samples_path), read_tsv(units_path)
        else:
            samples, units = make_bam_tables_from_barcode_map(args)
    else:
        if samples_path.exists() and units_path.exists() and not args.force:
            samples, units = read_tsv(samples_path), read_tsv(units_path)
        else:
            samples, units = make_bam_tables(args)

    if args.tumor_content is not None:
        for row in samples:
            row["tumor_content"] = args.tumor_content
    validate_tables(samples, units, allow_deferred_inputs=args.input_type == "minknow")

    sample_columns = list(samples[0]) if samples else ["sample"]
    unit_columns = list(BASE_UNIT_COLUMNS)
    for row in units:
        for column in row:
            if column not in unit_columns:
                unit_columns.append(column)
    write_tsv(samples_path, samples, sample_columns)
    write_tsv(units_path, units, unit_columns)
    runtime_config = make_runtime_config(args, args.analysis_dir)
    print(f"Created {samples_path}")
    print(f"Created {units_path}")
    print(f"Created {runtime_config}")


if __name__ == "__main__":
    main()
