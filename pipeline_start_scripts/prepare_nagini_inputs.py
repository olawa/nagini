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
import gzip
import hashlib
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any

import yaml


FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
FASTQ_STATUS_VALUES = ("all", "pass")
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
    "library_type",
    "assay",
]
DEFAULT_SAMPLE_REGEX = r"(?P<sample>[^_]+)(?:_.*)?\.(?:fastq|fq)(?:\.gz)?$"
DEFAULT_BASECALLING_MODEL = "dna_r10.4.1_e8.2_400bps_sup@v5.0.0"


def resolve_library_type_for_assay(
    assay: str, explicit_lib: str = "", default_lib: str = ""
) -> str:
    """Determine library_type from assay defaults unless an explicit type is provided."""
    explicit = str(explicit_lib or "").strip().lower()
    if explicit:
        return explicit
    norm_assay = str(assay or "").strip().lower()
    if norm_assay in ("bcr_abl", "cdna"):
        return "cdna"
    if norm_assay == "amplicon":
        return "amplicon"
    if norm_assay in ("wgs", "adaptive"):
        return "wgs"
    default_clean = str(default_lib or "").strip().lower()
    return default_clean if default_clean else "wgs"



def _pipeline_root(args) -> Path:
    """The checkout to import the routing module from.

    Falls back to this file's own checkout so that the table builders stay
    callable without a fully populated argument namespace.
    """

    root = getattr(args, "pipeline_root", None)
    return Path(root).resolve() if root else Path(__file__).resolve().parents[1]


def _routing_module(pipeline_root: Path):
    """Import the workflow's routing module so the launcher validates the same
    vocabulary the pipeline will use, rather than a second copy of it."""

    scripts = (pipeline_root / "workflow" / "scripts").resolve()
    if str(scripts) not in sys.path:
        sys.path.insert(0, str(scripts))
    import nagini_routing

    return nagini_routing


def known_library_types(pipeline_root: Path) -> set[str]:
    return set(_routing_module(pipeline_root).LIBRARY_TYPES)


def known_assays(pipeline_root: Path, config: dict) -> set[str]:
    routing = _routing_module(pipeline_root)
    try:
        return set(routing.resolve_assays(config))
    except routing.RoutingError as error:
        raise SystemExit(f"Invalid assay registry in the config: {error}")


def parse_args() -> argparse.Namespace:
    repo_root = Path(__file__).resolve().parent.parent
    default_base_config = repo_root / "config" / "config.yaml"
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--analysis-dir", type=Path, default=None)
    parser.add_argument(
        "--init-sample-sheet",
        nargs="?",
        const="sample_sheet.csv",
        default=None,
        help="Generate a sample sheet CSV template for the input directory and exit.",
    )
    parser.add_argument(
        "--pipeline-root",
        type=Path,
        default=repo_root,
        help="Path to Nagini repository root (default: auto-detected).",
    )
    parser.add_argument(
        "--base-config",
        type=Path,
        default=default_base_config,
        help="Path to base config.yaml (default: <pipeline-root>/config/config.yaml).",
    )
    parser.add_argument(
        "--site",
        choices=("local", "workstation", "marvin"),
        default="workstation",
        help="Runtime site for default site-config and reference-files (default: workstation).",
    )
    parser.add_argument("--site-config", type=Path)
    parser.add_argument("--reference-files", type=Path)
    parser.add_argument("--library-root", type=Path)
    parser.add_argument(
        "--input-type", choices=("minknow", "fastq", "unmapped_bam", "aligned_bam", ""), default=""
    )

    parser.add_argument("--sample-type", choices=("N", "T", "R"), default="N")
    parser.add_argument(
        "--reference-build",
        default="",
        help=(
            "Override config['reference_build'] for this run, e.g. hg19 for the "
            "bcr_abl assay. Must be a key already present in config['references']."
        ),
    )
    parser.add_argument(
        "--assay",
        default="",
        help=(
            "Default assay for every sample, for example wgs, adaptive or bcr_abl. "
            "Validated against the assay registry in the config. A per-sample "
            "assay column in the barcode map overrides it."
        ),
    )
    parser.add_argument(
        "--library-type",
        choices=("wgs", "amplicon", "cdna", ""),
        default="",
        help="Default library type when not specified in barcode map (default: derived from assay or wgs).",
    )
    parser.add_argument("--sequence-id", default="")
    parser.add_argument("--sample-regex", default=DEFAULT_SAMPLE_REGEX)
    parser.add_argument(
        "--barcode-map",
        type=Path,
        help=(
            "CSV/TSV mapping barcode or barcodeXX to sample for existing barcode "
            "FASTQ/BAM input. sample_sheet*.csv/tsv is auto-detected when omitted."
        ),
    )
    parser.add_argument("--barcode", default="NB")
    parser.add_argument("--run-id", default="")
    parser.add_argument("--basecalling-model", default=DEFAULT_BASECALLING_MODEL)
    parser.add_argument(
        "--fastq-status",
        choices=FASTQ_STATUS_VALUES,
        default="all",
        help=(
            "For existing barcode FASTQ input, include all fastq_pass and fastq_fail "
            "files (default) or only fastq_pass files."
        ),
    )

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
        "--dorado-demux-both-ends",
        action="store_true",
        default=False,
        help="Require barcode detection at both ends of reads in Dorado demux.",
    )
    parser.add_argument(
        "--dorado-demux-extra",
        default="",
        help="Extra flags passed to dorado demux.",
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
    parser.add_argument(
        "--no-check-references",
        dest="check_references",
        action="store_false",
        help="Skip the existence check on the reference registry.",
    )
    parser.add_argument(
        "--verify-references",
        action="store_true",
        help="Also verify registry checksums. Reads every file, so it is slow.",
    )
    parser.add_argument(
        "--filter-reads",
        dest="filter_reads",
        action="store_true",
        default=None,
        help="Enable pre-alignment read filtering (nanofilter).",
    )
    parser.add_argument(
        "--no-filter-reads",
        dest="filter_reads",
        action="store_false",
        help="Disable pre-alignment read filtering.",
    )
    parser.add_argument(
        "--merged-fastq",
        dest="merged_fastq",
        action="store_true",
        default=None,
        help="Write one FASTQ per sample, merged across flowcells, to reads/nagini_fastq/.",
    )
    parser.add_argument(
        "--no-merged-fastq",
        dest="merged_fastq",
        action="store_false",
        default=None,
        help="Do not write a merged per-sample FASTQ (the default).",
    )
    parser.add_argument(
        "--filter-qv",
        default="",
        help="Minimum mean read quality QV for nanofilter (e.g. 20).",
    )
    parser.add_argument(
        "--filter-len",
        default="",
        help="Read length range for nanofilter (e.g. 500-3000).",
    )
    parser.add_argument(
        "--filter-binary",
        default="",
        help="Path to nanofilter binary or script.",
    )
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    args.input_dir = args.input_dir.expanduser().resolve()
    if args.analysis_dir:
        args.analysis_dir = args.analysis_dir.expanduser().resolve()
    args.pipeline_root = args.pipeline_root.expanduser().resolve()
    args.base_config = args.base_config.expanduser().resolve()

    if args.site_config:
        args.site_config = args.site_config.expanduser().resolve()
    else:
        site_candidate = args.pipeline_root / "config" / "site_configs" / f"site_config_{args.site}.yaml"
        if site_candidate.is_file():
            args.site_config = site_candidate

    if args.reference_files:
        args.reference_files = args.reference_files.expanduser().resolve()
    else:
        ref_candidate = args.pipeline_root / "config" / "reference_files" / f"reference_files_{args.site}.yaml"
        if ref_candidate.is_file():
            args.reference_files = ref_candidate

    if getattr(args, "barcode_map", None):
        args.barcode_map = args.barcode_map.expanduser().resolve()
    if getattr(args, "bed_file", None):
        args.bed_file = Path(args.bed_file).expanduser().resolve()
    if getattr(args, "library_root", None):
        args.library_root = Path(args.library_root).expanduser().resolve()

    return args





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


def format_short_barcode(value: object) -> str:
    """Format barcode01 -> bc1, barcode10 -> bc10, bc01 -> bc1, 01 -> bc1, etc."""
    barcode = str(value or "").strip().lower()
    match = re.search(r"(\d+)", barcode)
    if match:
        return f"bc{int(match.group(1))}"
    return barcode if barcode else "bc"


def make_barcode_sample_name(base_sample: str, barcode: str) -> str:
    """Combine base sample name and barcode into a sample name with barcode info (bc1, bc10, etc.)."""
    short_bc = format_short_barcode(barcode)
    base = str(base_sample or "").strip()
    generic_names = {
        "sample",
        "reads",
        "read",
        "fastq",
        "fastqs",
        "bam",
        "bams",
        "data",
        "inputs",
        "run",
        "runs",
        "inbox",
        "combined",
        "minknow",
        "pass",
        "fail",
    }
    if not base or base.lower() in generic_names:
        return short_bc
    if normalize_barcode(base) == normalize_barcode(barcode) or base.lower() == short_bc:
        return short_bc

    pattern = rf"(?i)[-_]?(?:{re.escape(short_bc)}|{re.escape(normalize_barcode(barcode))})$"
    if re.search(pattern, base):
        prefix = re.sub(pattern, "", base).rstrip("-_")
        return f"{prefix}_{short_bc}" if prefix else short_bc

    return f"{base}_{short_bc}"


def sample_row(sample_name: str, assignment: dict) -> dict[str, str]:
    """A samples.tsv row, carrying the lab identifier when the map gave one."""

    row = {"sample": sample_name}
    dna_number = str(assignment.get("dna_number", "") or "").strip()
    if dna_number:
        row["dna_number"] = dna_number
    return row


def read_barcode_map(
    path: Path,
    default_type: str,
    default_library_type: str = "wgs",
    default_assay: str = "",
    valid_library_types: set[str] | None = None,
) -> dict[str, dict[str, str]]:
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
        barcode = normalize_barcode(
            fields.get("barcode", fields.get("barcodename", fields.get("barcodeid", "")))
        )
        sample_name = next(
            (
                fields.get(key, "")
                for key in ("sample", "sampleid", "samplename", "alias", "name")
                if fields.get(key, "")
            ),
            "",
        )
        sample_type = fields.get("type", fields.get("sampletype", "")) or default_type
        sample_type = sample_type.upper()
        # The lab identifier the report should show. Kept apart from `sample`,
        # which names files and so has to stay filename-safe and stable: a run
        # can be delivered under a real DNA number without renaming anything.
        dna_number = next(
            (
                fields.get(key, "")
                for key in ("dnanumber", "dna", "dnanr", "labid", "label")
                if fields.get(key, "")
            ),
            "",
        )
        raw_lib = fields.get("librarytype", fields.get("library_type", ""))
        assay = str(fields.get("assay", "") or default_assay).strip().lower()
        library_type = resolve_library_type_for_assay(
            assay,
            explicit_lib=raw_lib,
            default_lib=default_library_type,
        )
        if not barcode or not re.fullmatch(r"barcode\d+", barcode):

            raise SystemExit(f"Invalid or missing barcode in {path}: {row}")
        if not sample_name:
            raise SystemExit(f"Missing sample/sample_id for {barcode} in {path}.")
        if sample_type not in {"N", "T", "R"}:
            raise SystemExit(
                f"Invalid type={sample_type!r} for {barcode}; expected N, T or R."
            )
        allowed = valid_library_types or {"wgs", "amplicon", "cdna", "direct_rna"}
        if library_type not in allowed:
            raise SystemExit(
                f"Invalid library_type={library_type!r} for {barcode} in {path}; "
                f"expected one of {', '.join(sorted(allowed))}."
            )
        normalized_rows.append(
            {
                "barcode": barcode,
                "sample": sample_name,
                "type": sample_type,
                "library_type": library_type,
                "assay": assay,
                "dna_number": dna_number,
            }
        )

    mapping: dict[str, dict[str, str]] = {}
    for row in normalized_rows:
        if row["barcode"] in mapping:
            raise SystemExit(f"Barcode appears more than once in {path}: {row['barcode']}")
        mapping[row["barcode"]] = row
    return mapping


def _delimited_header(path: Path) -> set[str]:
    """Return normalized column names without requiring the file to be a map."""

    try:
        with path.open(newline="") as handle:
            sample = handle.read(4096)
            handle.seek(0)
            try:
                dialect = csv.Sniffer().sniff(sample, delimiters=",\t")
            except csv.Error:
                dialect = csv.excel_tab if path.suffix.lower() in {".tsv", ".tab"} else csv.excel
            reader = csv.reader(handle, dialect=dialect)
            header = next(reader, [])
    except (OSError, UnicodeError, StopIteration):
        return set()
    return {
        re.sub(r"[^a-z0-9]", "", str(column or "").lower())
        for column in header
    }


def discover_barcode_map(
    input_dir: Path,
    default_type: str,
    default_library_type: str = "wgs",
    default_assay: str = "",
    valid_library_types: set[str] | None = None,
) -> tuple[dict[str, dict[str, str]] | None, list[Path]]:
    """Find and combine ONT ``sample_sheet*.csv/tsv`` barcode maps.

    MinKNOW/Dorado exports commonly keep the sample sheet beside the FASTQs.
    Treating it as an implicit map makes a directory of existing runs usable
    without first copying the rows into a second hand-written file.  Multiple
    run sheets are allowed when they agree for every barcode.
    """

    candidates = []
    for path in sorted(input_dir.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in {".csv", ".tsv", ".tab"}:
            continue
        stem = path.stem.lower().replace("-", "_")
        if not (
            stem.startswith("sample_sheet")
            or stem.startswith("samplesheet")
            or "barcode_map" in stem
            or "sample_map" in stem
            or "barcodes" in stem
        ):
            continue
        header = _delimited_header(path)
        has_sample = bool(
            header
            & {"sample", "sampleid", "samplename", "alias", "name"}
        )
        if "barcode" in header or "barcodename" in header or "barcodeid" in header:
            if has_sample:
                candidates.append(path)


    if not candidates:
        return None, []

    combined: dict[str, dict[str, str]] = {}
    for path in candidates:
        parsed = read_barcode_map(
            path,
            default_type,
            default_library_type,
            default_assay,
            valid_library_types,
        )
        for barcode, row in parsed.items():
            previous = combined.get(barcode)
            if previous and any(
                previous.get(column, "") != row.get(column, "")
                for column in ("sample", "type", "library_type", "assay")
            ):
                raise SystemExit(
                    f"Barcode {barcode} has conflicting assignments in sample sheets "
                    f"({previous['sample']!r} versus {row['sample']!r}); check {path}."
                )
            combined[barcode] = row
    return combined, candidates


def _barcode_from_path(path: Path, input_dir: Path) -> str:
    """Extract a barcode from a directory component or the filename."""

    root_barcode = normalize_barcode(input_dir.name)
    if not re.fullmatch(r"barcode\d+", root_barcode):
        root_barcode = ""
    relative_parts = path.relative_to(input_dir).parts[:-1]
    for part in reversed(relative_parts):
        candidate = normalize_barcode(part)
        if re.fullmatch(r"barcode\d+", candidate):
            return candidate

    # Dorado demux versions may emit ``barcode01.bam``/FASTQ directly below
    # the output directory instead of creating barcode subdirectories.
    for token in re.findall(r"(?i)(?:barcode|bc)\d+", path.stem):
        candidate = normalize_barcode(token)
        if re.fullmatch(r"barcode\d+", candidate):
            return candidate
    return root_barcode


def discover_bams_by_barcode(input_dir: Path) -> dict[str, list[Path]]:
    """Find BAMs below barcodeXX directories or with barcodeXX filenames."""

    bams: defaultdict[str, list[Path]] = defaultdict(list)
    for path in sorted(input_dir.rglob("*.bam")):
        barcode = _barcode_from_path(path, input_dir)
        if barcode:
            bams[barcode].append(path.resolve())
    return dict(bams)


def _fastq_is_allowed(path: Path, input_dir: Path, status: str) -> bool:
    """Select all FASTQs or only files below a Dorado ``fastq_pass`` tree."""

    if status == "all":
        return True
    parts = {part.lower() for part in path.relative_to(input_dir).parts[:-1]}
    return "fastq_fail" not in parts


def discover_fastqs_by_barcode(
    input_dir: Path,
    status: str = "all",
    exclude_dir: Path | None = None,
) -> dict[str, list[Path]]:
    """Find existing barcode FASTQs across one or more run directories."""

    if status not in FASTQ_STATUS_VALUES:
        raise SystemExit(
            f"Unsupported FASTQ status {status!r}; expected one of {', '.join(FASTQ_STATUS_VALUES)}."
        )
    fastqs: defaultdict[str, list[Path]] = defaultdict(list)
    excluded = exclude_dir.resolve() if exclude_dir else None
    for path in sorted(input_dir.rglob("*")):
        if not path.is_file() or not path.name.lower().endswith(FASTQ_SUFFIXES):
            continue
        resolved = path.resolve()
        if excluded and (resolved == excluded or excluded in resolved.parents):
            continue
        if not _fastq_is_allowed(path, input_dir, status):
            continue
        barcode = _barcode_from_path(path, input_dir)
        if barcode:
            fastqs[barcode].append(resolved)
    return dict(fastqs)


def generate_sample_sheet_template(args: argparse.Namespace) -> Path:
    """Generate a pre-filled sample sheet CSV template for the input directory."""
    input_dir = Path(args.input_dir).resolve()
    if not input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist: {input_dir}")

    target_str = getattr(args, "init_sample_sheet", None) or "sample_sheet.csv"
    target_path = Path(target_str)
    if not target_path.is_absolute():
        target_path = input_dir / target_path

    # Discover barcodes from existing FASTQ or BAM files
    barcodes: list[str] = []
    fastqs = discover_fastqs_by_barcode(input_dir, status="all")
    if fastqs:
        barcodes.extend(fastqs.keys())

    bams = discover_bams_by_barcode(input_dir)
    if bams:
        barcodes.extend(bams.keys())

    kit = str(getattr(args, "dorado_kit", "") or getattr(args, "barcode_kit", "") or "").strip()
    if kit:
        routing = _routing_module(_pipeline_root(args))
        kit_size = routing.barcode_kit_size(kit)
        if kit_size:
            barcodes.extend(f"barcode{i:02d}" for i in range(1, kit_size + 1))

    unique_barcodes = sorted(
        {normalize_barcode(bc) for bc in barcodes if normalize_barcode(bc)},
        key=lambda bc: int(re.search(r"\d+", bc).group(0)) if re.search(r"\d+", bc) else 0,
    )

    if not unique_barcodes:
        unique_barcodes = [f"barcode{i:02d}" for i in range(1, 25)]

    default_type = str(getattr(args, "sample_type", "T") or "T").strip().upper()
    default_assay = str(getattr(args, "assay", "") or "").strip().lower() or "wgs"
    default_lib = resolve_library_type_for_assay(default_assay, getattr(args, "library_type", ""), "wgs")

    target_path.parent.mkdir(parents=True, exist_ok=True)
    with target_path.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["barcode", "sample", "type", "library_type", "assay"])
        for bc in unique_barcodes:
            short_bc = format_short_barcode(bc)
            writer.writerow([bc, f"sample_{short_bc}", default_type, default_lib, default_assay])

    print(f"Created sample sheet template with {len(unique_barcodes)} barcode(s) at: {target_path}")
    return target_path


def discover_fastq_files(
    input_dir: Path,
    status: str = "all",
    exclude_dir: Path | None = None,
) -> list[Path]:
    """Find FASTQ files while applying the same pass/fail selection as barcode mode."""

    if status not in FASTQ_STATUS_VALUES:
        raise SystemExit(
            f"Unsupported FASTQ status {status!r}; expected one of {', '.join(FASTQ_STATUS_VALUES)}."
        )
    excluded = exclude_dir.resolve() if exclude_dir else None
    files = []
    for path in sorted(input_dir.rglob("*")):
        if not path.is_file() or not path.name.lower().endswith(FASTQ_SUFFIXES):
            continue
        resolved = path.resolve()
        if excluded and (resolved == excluded or excluded in resolved.parents):
            continue
        if _fastq_is_allowed(path, input_dir, status):
            files.append(resolved)
    return files


def merge_fastq_files(input_files: list[Path], output: Path) -> None:
    """Concatenate FASTQ inputs into one gzip stream without loading them in memory."""

    if not input_files:
        raise SystemExit(f"Cannot create merged FASTQ {output}: no input files.")
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.name}.{os.getpid()}.tmp")
    try:
        if all(path.name.lower().endswith(".gz") for path in input_files):
            # Concatenated gzip members are valid and are understood by samtools,
            # minimap2 and gzip. This keeps merging large runs effectively I/O-only.
            with temporary.open("wb") as destination:
                for path in input_files:
                    with path.open("rb") as source:
                        shutil.copyfileobj(source, destination, length=1024 * 1024)
        else:
            with gzip.open(temporary, "wb") as destination:
                for path in input_files:
                    opener = gzip.open if path.name.lower().endswith(".gz") else open
                    with opener(path, "rb") as source:
                        shutil.copyfileobj(source, destination, length=1024 * 1024)
        os.replace(temporary, output)
    except OSError:
        temporary.unlink(missing_ok=True)
        raise


def flowcell_token_from_path(path: Path) -> str:
    """The flowcell a MinKNOW BAM came from, as a filename-safe token.

    MinKNOW names its run directory ``<date>_<time>_<position>_<flowcell>_<id>``,
    which is the only place the flowcell is recorded for an already
    demultiplexed BAM.  Naming the processing unit after it makes a merged
    multi-flowcell sample say which flowcell each unit came from, instead of
    a counter whose order depends on how the files happened to sort.
    """

    for parent in path.parents:
        parts = parent.name.split("_")
        if len(parts) >= 5 and re.fullmatch(r"\d{8}", parts[0]) and re.fullmatch(r"\d{4}", parts[1]):
            token = f"{parts[-2]}-{parts[-1]}"
            if re.fullmatch(r"[A-Za-z0-9.-]+", token):
                return token
    return ""


def name_processing_units(bam_paths: list[Path], barcode: str) -> list[tuple[Path, str]]:
    """Pair each BAM with a processing unit that is unique for the sample.

    A flowcell can hold several BAMs for one barcode (MinKNOW splits large
    outputs into ``_0.bam``, ``_1.bam``), so the token alone is not unique;
    the second and later files within one flowcell get a counter.  Inputs with
    no recognisable MinKNOW run directory keep the old positional naming.
    """

    tokens = [(path, flowcell_token_from_path(path)) for path in bam_paths]
    if not all(token for _, token in tokens):
        return [
            (path, f"{barcode}_{index:03d}")
            for index, (path, _) in enumerate(tokens, start=1)
        ]
    seen: defaultdict[str, int] = defaultdict(int)
    named = []
    for path, token in tokens:
        seen[token] += 1
        named.append((path, token if seen[token] == 1 else f"{token}-{seen[token]}"))
    return named


def make_bam_tables_from_barcode_map(
    args: argparse.Namespace,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Build Hydra-style tables from Dorado ``bam_pass/barcodeXX`` output."""

    mapping = None
    if getattr(args, "barcode_map", None):
        mapping = read_barcode_map(
            args.barcode_map.resolve(),
            args.sample_type,
            getattr(args, "library_type", "wgs"),
            getattr(args, "assay", ""),
            known_library_types(_pipeline_root(args)),
        )
    else:
        mapping, map_sources = discover_barcode_map(
            args.input_dir,
            args.sample_type,
            getattr(args, "library_type", ""),
            getattr(args, "assay", ""),
            known_library_types(_pipeline_root(args)),
        )

    bams_by_barcode = discover_bams_by_barcode(args.input_dir)
    if not bams_by_barcode:
        raise SystemExit(f"No BAM files below barcode directories in {args.input_dir}.")

    if not mapping:
        library_type = resolve_library_type_for_assay(
            getattr(args, "assay", ""), getattr(args, "library_type", ""), "wgs"
        )
        base_prefix = args.sequence_id or args.run_id or ""
        mapping = {
            barcode: {
                "barcode": barcode,
                "sample": make_barcode_sample_name(base_prefix, barcode),
                "type": str(args.sample_type).strip().upper(),
                "library_type": library_type,
                "assay": str(getattr(args, "assay", "") or "").strip().lower(),
            }
            for barcode in sorted(bams_by_barcode)
        }
        print(
            "Auto-inferred barcode map for BAM files: "
            + ", ".join(f"{bc}->{m['sample']}" for bc, m in sorted(mapping.items()))
        )

    observed = set(bams_by_barcode)
    expected = set(mapping)
    missing = sorted(expected - observed)
    unmapped = sorted(observed - expected)
    if unmapped:
        print(
            "WARNING: Ignoring BAM barcode directories absent from barcode map: "
            + ", ".join(unmapped),
            file=sys.stderr,
        )
    if missing:
        print(
            "WARNING: Barcode map contains barcodes with no BAM files: "
            + ", ".join(missing),
            file=sys.stderr,
        )
    selected_barcodes = sorted(expected & observed)
    if not selected_barcodes:
        raise SystemExit(
            f"None of the barcodes in the barcode map ({', '.join(sorted(expected))}) "
            f"were found in {args.input_dir}."
        )

    run_id = args.run_id or args.sequence_id or args.input_dir.name
    samples_by_name: dict[str, dict[str, str]] = {}
    units: list[dict[str, str]] = []
    for barcode in selected_barcodes:
        assignment = mapping[barcode]

        sample_name = assignment["sample"]
        samples_by_name.setdefault(sample_name, sample_row(sample_name, assignment))
        for bam_path, unit_name in name_processing_units(
            sorted(bams_by_barcode[barcode]), barcode
        ):
            units.append(
                {
                    "sample": sample_name,
                    "type": assignment["type"],
                    "platform": "ONT",
                    "machine": "unknown",
                    "processing_unit": unit_name,
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
                    "run_id": flowcell_token_from_path(bam_path) or run_id,
                    "basecalling_model": args.basecalling_model,
                    "library_type": assignment.get(
                        "library_type", getattr(args, "library_type", "wgs")
                    ),
                    "assay": assignment.get("assay", "") or getattr(args, "assay", ""),
                }
            )
    return sorted(samples_by_name.values(), key=lambda row: row["sample"]), units


# ---------------------------------------------------------------------------
# Reference registry
# ---------------------------------------------------------------------------
#
# The registry is the inventory of reference assets for a site: where each file
# is, whether it is required, and what it should hash to.  It is flattened into
# the same {{VARIABLE}} names the shared config templates use, so a reference
# path is declared once, in the registry, rather than twice.
#
# Precedence is base config < registry < site config, so a site can still pin
# one path without taking over the whole inventory.

REGISTRY_BUILD_VARIABLES = {
    "fasta": "FASTA",
    "fai": "FAI",
    "dict": "DICT",
    "mmi": "MMI",
    "minimap2_index": "MMI",
    "design_bed": "DESIGN_BED",
    "sv_annotation_bed": "SV_ANNOTATION_BED",
    "tandem_repeats": "TANDEM_REPEATS",
    "clinvar": "CLINVAR",
    "clinvar_index": "CLINVAR_INDEX",
    "cytoband": "CYTOBAND",
    "gaps_bed": "GAPS_BED",
    "gene_annotation": "GTF",
    "annotations": "ANNOTATIONS",
}

REGISTRY_ANNOTATION_VARIABLES = {
    "vep_cache": "VEP_CACHE",
    "snpeff_db": "SNPEFF_DB",
    "simple_sv_panel": "SIMPLE_SV_PANEL",
    "simple_sv_fusion_pairs": "SIMPLE_SV_FUSION_PAIRS",
}


class ReferenceEntry:
    """One asset: the variable it fills, its path, and how to check it."""

    __slots__ = ("variable", "path", "kind", "checksum", "source")

    def __init__(self, variable, path, kind="file", checksum="", source=""):
        self.variable = variable
        self.path = path
        self.kind = (kind or "file").strip().lower()
        self.checksum = (checksum or "").strip()
        self.source = source

    @property
    def required(self) -> bool:
        return self.kind not in ("optional", "") and bool(self.path)


def _registry_entry(variable, value, reference_root, source):
    """Accept either a bare path or a {path, type, url, checksum} mapping."""

    if value is None:
        return None
    if isinstance(value, str):
        path, kind, checksum = value.strip(), "file", ""
    elif isinstance(value, dict):
        path = str(value.get("path", "") or "").strip()
        kind = str(value.get("type", "file") or "file").strip().lower()
        checksum = str(value.get("checksum", "") or "").strip()
    else:
        raise SystemExit(
            f"{source}: expected a path or a mapping, got {type(value).__name__}."
        )
    if path and reference_root and not os.path.isabs(path):
        if not path.startswith(("workflow/", "config/")):
            path = os.path.join(reference_root, path)
    return ReferenceEntry(variable, path, kind, checksum, source)


def flatten_reference_registry(registry: dict, default_root: str = "") -> list[ReferenceEntry]:
    """Turn a registry into the {{VARIABLE}} entries the config templates use.

    Relative ``path`` values are joined to a reference root, so a registry can
    say where an asset sits inside the site's reference tree without repeating
    the tree's location on every line. The root is the registry's own
    ``reference_root`` when it sets one, otherwise the site's
    ``REFERENCE_ROOT`` -- so a site declares that location once.
    """

    if not registry:
        return []
    reference_root = str(
        registry.get("reference_root", "") or default_root or ""
    ).strip()
    entries: list[ReferenceEntry] = []

    for build, assets in (registry.get("reference_builds") or {}).items():
        if not isinstance(assets, dict):
            raise SystemExit(f"reference_builds.{build} must be a mapping.")
        prefix = str(build).strip().upper()
        for key, value in assets.items():
            suffix = REGISTRY_BUILD_VARIABLES.get(str(key).strip().lower())
            if not suffix:
                continue
            entry = _registry_entry(
                f"{prefix}_{suffix}", value, reference_root, f"reference_builds.{build}.{key}"
            )
            if entry:
                entries.append(entry)

    for key, value in (registry.get("annotations") or {}).items():
        variable = REGISTRY_ANNOTATION_VARIABLES.get(str(key).strip().lower())
        if not variable:
            continue
        entry = _registry_entry(variable, value, reference_root, f"annotations.{key}")
        if entry:
            entries.append(entry)
    return entries


def registry_variables(registry: dict, default_root: str = "") -> dict[str, str]:
    """Every variable the registry declares, empty paths included.

    An entry with no path is a statement that the asset is not installed here,
    which has to override the shared config's default rather than fall through
    to it.
    """

    return {
        entry.variable: entry.path
        for entry in flatten_reference_registry(registry, default_root)
    }


def _md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def validate_reference_registry(
    registry: dict,
    verify_checksums: bool = False,
    strict: bool = False,
    default_root: str = "",
) -> None:
    """Check the inventory before Snakemake builds a DAG against it.

    Reports rather than raises by default: the shipped registries describe more
    assets than every site installs, and a launcher that refuses to start over
    an unused reference would be worse than one that says what is missing.
    ``--verify-references`` makes it strict and adds checksums.
    """

    entries = flatten_reference_registry(registry, default_root)
    missing = []
    wrong_kind = []
    bad_checksum = []
    checked = 0

    for entry in entries:
        if not entry.path:
            if entry.required:
                missing.append(f"{entry.source} has no path")
            continue
        path = Path(entry.path)
        if entry.kind == "directory":
            if not path.is_dir():
                (missing if entry.required else wrong_kind).append(
                    f"{entry.source}: not a directory: {path}"
                )
                continue
        elif not path.exists():
            if entry.required:
                missing.append(f"{entry.source}: {path}")
            continue
        elif path.is_dir() and entry.kind == "file":
            wrong_kind.append(f"{entry.source}: expected a file, found a directory: {path}")
            continue

        if verify_checksums and entry.checksum and path.is_file():
            actual = _md5(path)
            checked += 1
            if actual != entry.checksum:
                bad_checksum.append(
                    f"{entry.source}: {path} is {actual}, registry says {entry.checksum}"
                )

    problems = missing + wrong_kind + bad_checksum
    if problems:
        message = "Reference registry problems:\n  " + "\n  ".join(problems)
        if strict or bad_checksum:
            # A wrong checksum is never tolerable: the file is not the one the
            # registry describes, whatever else is true.
            raise SystemExit(message)
        print("WARNING: " + message, file=sys.stderr)
    if verify_checksums:
        print(f"Verified {checked} reference checksum(s).")


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


def infer_barcode_map_from_fastq_files(
    args: argparse.Namespace, files: list[Path]
) -> dict[str, dict[str, str]]:
    """Infer barcode assignments when FASTQ names contain sample and barcode.

    When no sample sheet is provided, each barcode is assigned a sample name
    with barcode info (bc1, bc10, etc.), e.g. ``PBE82562_bc1`` or ``bc1``.
    """

    pattern = re.compile(args.sample_regex)
    input_dir = Path(args.input_dir).resolve()
    samples_by_barcode: defaultdict[str, set[str]] = defaultdict(set)
    for path in files:
        barcode = _barcode_from_path(path, input_dir)
        if not barcode:
            continue
        try:
            sample = sample_from_fastq(path, pattern).strip()
        except SystemExit:
            sample = ""
        samples_by_barcode[barcode].add(sample)

    if not samples_by_barcode:
        return {}

    library_type = resolve_library_type_for_assay(
        getattr(args, "assay", ""), getattr(args, "library_type", ""), "wgs"
    )
    default_prefix = args.sequence_id or args.run_id or ""
    all_observed_samples = sorted(
        {s for samples in samples_by_barcode.values() for s in samples if s}
    )
    common_sample_prefix = (
        all_observed_samples[0] if len(all_observed_samples) == 1 else ""
    )

    mapping = {}
    multi_flowcell = []
    for barcode, samples in sorted(samples_by_barcode.items()):
        non_empty = sorted(set(s for s in samples if s))
        if len(non_empty) > 1:
            multi_flowcell.append(f"{barcode} ({', '.join(non_empty)})")

        if common_sample_prefix:
            base_sample = common_sample_prefix
        elif default_prefix:
            base_sample = default_prefix
        else:
            base_sample = ""
        sample_name = make_barcode_sample_name(base_sample, barcode)
        mapping[barcode] = {
            "barcode": barcode,
            "sample": sample_name,
            "type": str(args.sample_type).strip().upper(),
            "library_type": library_type,
            "assay": str(getattr(args, "assay", "") or "").strip().lower(),
        }
    if multi_flowcell:
        print(
            "Merging same-barcode reads across runs/flowcells for: "
            + "; ".join(multi_flowcell)
        )
    return mapping


def make_barcode_fastq_tables(
    args: argparse.Namespace,
    mapping: dict[str, dict[str, str]],
    fastqs_by_barcode: dict[str, list[Path]] | None = None,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Merge all run/chunk FASTQs for each barcode and create one unit per barcode."""

    merged_dir = Path(args.analysis_dir).resolve() / "input_fastq"
    fastqs_by_barcode = fastqs_by_barcode or discover_fastqs_by_barcode(
        args.input_dir,
        getattr(args, "fastq_status", "all"),
        exclude_dir=merged_dir,
    )
    if not fastqs_by_barcode:
        raise SystemExit(f"No barcode FASTQ files found below {args.input_dir}.")

    observed = set(fastqs_by_barcode)
    expected = set(mapping)
    unmapped = sorted(observed - expected)
    if unmapped:
        print(
            "WARNING: Ignoring barcode FASTQ directories/files absent from barcode map: "
            + ", ".join(unmapped),
            file=sys.stderr,
        )
    missing = sorted(expected - observed)
    if missing:
        print(
            "WARNING: Barcode map contains barcodes with no selected FASTQ files: "
            + ", ".join(missing),
            file=sys.stderr,
        )
    selected_barcodes = sorted(expected & observed)
    if not selected_barcodes:
        raise SystemExit(
            f"None of the barcodes in the barcode map ({', '.join(sorted(expected))}) "
            f"were found in {args.input_dir}."
        )

    run_id = args.run_id or args.sequence_id or args.input_dir.name
    samples_by_name: dict[str, dict[str, str]] = {}
    units: list[dict[str, str]] = []
    for barcode in selected_barcodes:
        assignment = mapping[barcode]

        sample_name = assignment["sample"]
        samples_by_name.setdefault(sample_name, sample_row(sample_name, assignment))
        source_files = sorted(fastqs_by_barcode[barcode])
        short_bc = format_short_barcode(barcode)
        if short_bc in sample_name.lower() or barcode.lower() in sample_name.lower():
            fastq_name = f"{sample_name}.fastq.gz"
        else:
            fastq_name = f"{sample_name}_{short_bc}.fastq.gz"
        merged_fastq = merged_dir / fastq_name
        merge_fastq_files(source_files, merged_fastq)
        units.append(
            {
                "sample": sample_name,
                "type": assignment["type"],
                "platform": "ONT",
                "machine": "unknown",
                "processing_unit": f"{barcode}_merged",
                "barcode": "BC",
                "barcode_name": barcode,
                "input_type": "fastq",
                "fastq1": str(merged_fastq),
                "fastq2": "UNUSED",
                "bam": "UNUSED",
                "methylation": args.methylation,
                "run_id": run_id,
                "basecalling_model": args.basecalling_model,
                "library_type": assignment.get(
                    "library_type", getattr(args, "library_type", "wgs")
                ),
                "assay": assignment.get("assay", "") or getattr(args, "assay", ""),
            }
        )
    return sorted(samples_by_name.values(), key=lambda row: row["sample"]), units


def make_fastq_tables(args: argparse.Namespace) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    merged_dir = Path(getattr(args, "analysis_dir", args.input_dir)).resolve() / "input_fastq"
    files = discover_fastq_files(
        args.input_dir,
        getattr(args, "fastq_status", "all"),
        exclude_dir=merged_dir,
    )
    if not files:
        raise SystemExit(f"No FASTQ files found below {args.input_dir}.")

    mapping = None
    map_sources: list[Path] = []
    if getattr(args, "barcode_map", None):
        mapping = read_barcode_map(
            args.barcode_map.resolve(),
            args.sample_type,
            getattr(args, "library_type", "wgs"),
            getattr(args, "assay", ""),
            known_library_types(_pipeline_root(args)),
        )
        map_sources = [args.barcode_map.resolve()]
    else:
        mapping, map_sources = discover_barcode_map(
            args.input_dir,
            args.sample_type,
            getattr(args, "library_type", ""),
            getattr(args, "assay", ""),
            known_library_types(_pipeline_root(args)),
        )
    inferred_mapping = False
    if not mapping:
        mapping = infer_barcode_map_from_fastq_files(args, files)
        inferred_mapping = bool(mapping)
    if mapping:
        if map_sources and not getattr(args, "barcode_map", None):
            print(
                "Auto-detected barcode map(s): "
                + ", ".join(str(path) for path in map_sources)
            )
        elif inferred_mapping:
            print(
                "Auto-inferred barcode map from FASTQ filenames for: "
                + ", ".join(sorted(mapping))
            )
        return make_barcode_fastq_tables(args, mapping)

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
                "basecalling_model": args.basecalling_model or DEFAULT_BASECALLING_MODEL,
                "library_type": resolve_library_type_for_assay(
                    getattr(args, "assay", ""), getattr(args, "library_type", ""), "wgs"
                ),
                "assay": getattr(args, "assay", ""),
            }
        )

    samples = [{"sample": sample} for sample in sorted(counters)]
    return samples, units


def make_minknow_tables(args: argparse.Namespace) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    """Create deferred units for Dorado outputs produced by Snakemake."""

    mapping = None
    if getattr(args, "barcode_map", None):
        mapping = read_barcode_map(
            args.barcode_map.resolve(),
            args.sample_type,
            getattr(args, "library_type", "wgs"),
            getattr(args, "assay", ""),
            known_library_types(_pipeline_root(args)),
        )
    else:
        mapping, _ = discover_barcode_map(
            args.input_dir,
            args.sample_type,
            getattr(args, "library_type", ""),
            getattr(args, "assay", ""),
            known_library_types(_pipeline_root(args)),
        )

    if not mapping:
        kit = getattr(args, "dorado_kit", "") or getattr(args, "barcode_kit", "")
        routing = _routing_module(_pipeline_root(args))
        kit_size = routing.barcode_kit_size(kit) if kit else None
        barcodes = []
        if kit_size:
            barcodes = [f"barcode{i:02d}" for i in range(1, kit_size + 1)]
        else:
            bams = discover_bams_by_barcode(args.input_dir)
            fqs = discover_fastqs_by_barcode(args.input_dir)
            barcodes = sorted(set(bams) | set(fqs))

        if not barcodes:
            raise SystemExit(
                "Cannot determine barcodes for MinKNOW input without --barcode-map. "
                "Provide --barcode-kit (e.g. SQK-RBK114-24) or --barcode-map."
            )

        library_type = resolve_library_type_for_assay(
            getattr(args, "assay", ""), getattr(args, "library_type", ""), "wgs"
        )
        base_prefix = args.sequence_id or args.run_id or ""
        mapping = {
            barcode: {
                "barcode": barcode,
                "sample": make_barcode_sample_name(base_prefix, barcode),
                "type": str(args.sample_type).strip().upper(),
                "library_type": library_type,
                "assay": str(getattr(args, "assay", "") or "").strip().lower(),
            }
            for barcode in sorted(barcodes)
        }
        print(
            "Auto-inferred barcode map for MinKNOW input: "
            + ", ".join(f"{bc}->{m['sample']}" for bc, m in sorted(mapping.items()))
        )
    run_id = args.run_id or args.sequence_id or args.input_dir.name
    samples_by_name: dict[str, dict[str, str]] = {}
    units: list[dict[str, str]] = []
    for barcode in sorted(mapping):
        assignment = mapping[barcode]
        sample_name = assignment["sample"]
        samples_by_name.setdefault(sample_name, sample_row(sample_name, assignment))
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
                "library_type": assignment.get(
                    "library_type", getattr(args, "library_type", "wgs")
                ),
                "assay": assignment.get("assay", "") or getattr(args, "assay", ""),
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
        barcode_name = ""
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
                barcode_name = match.groupdict().get("barcode", "")
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
        unit_entry = {
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
            "library_type": getattr(args, "library_type", "wgs"),
            "assay": getattr(args, "assay", ""),
        }
        if barcode_name:
            unit_entry["barcode_name"] = str(barcode_name).strip()
        units.append(unit_entry)

    return sorted(samples_by_name.values(), key=lambda row: row["sample"]), units


def can_reuse_existing_tables(
    samples_path: Path,
    units_path: Path,
    expected_input_type: str,
) -> bool:
    """Check if existing samples.tsv and units.tsv can be reused for the requested input type."""
    if not (samples_path.is_file() and units_path.is_file()):
        return False
    try:
        units = read_tsv(units_path)
        if not units:
            return False
        expected_unit_type = "dorado_bam" if expected_input_type == "minknow" else expected_input_type
        if any(row.get("input_type") != expected_unit_type for row in units):
            return False
        return True
    except Exception:
        return False


def detect_basecalling_model(input_dir: Path) -> str | None:
    """Inspect run summary files, BAM headers, or FASTQ read headers to detect the basecalling model."""
    if not input_dir or not input_dir.exists():
        return None

    # 1. Try MinKNOW / Dorado final_summary or report files
    try:
        for pattern in ("final_summary*.txt", "*summary*.txt"):
            for summary_path in input_dir.rglob(pattern):
                if not summary_path.is_file():
                    continue
                try:
                    with summary_path.open("r", errors="ignore") as f:
                        for line in f:
                            m = re.search(r"^(?:basecall_model_version_id|model)=([^\s,;]+)", line.strip())
                            if m:
                                val = m.group(1).strip()
                                if "@" in val or "dna_r" in val:
                                    return val
                except Exception:
                    continue
    except Exception:
        pass

    # 2. Try FASTQ files
    try:
        count = 0
        for path in input_dir.rglob("*"):
            if not path.is_file() or not path.name.lower().endswith(FASTQ_SUFFIXES):
                continue
            count += 1
            if count > 5:
                break
            try:
                opener = gzip.open if str(path).lower().endswith(".gz") else open
                with opener(str(path), "rt", errors="ignore") as handle:
                    for _ in range(50):
                        line = handle.readline()
                        if not line:
                            break
                        if line.startswith("@"):
                            m = re.search(r"(?:model_version_id|basecall_model_version_id)=([^\s,;]+)", line)
                            if m:
                                return m.group(1).strip()
            except Exception:
                continue
    except Exception:
        pass

    # 3. Try BAM files
    try:
        count = 0
        for path in input_dir.rglob("*.bam"):
            if not path.is_file():
                continue
            count += 1
            if count > 5:
                break
            try:
                import pysam
                with pysam.AlignmentFile(str(path), "rb", check_sq=False) as bf:
                    for pg in bf.header.get("PG", []):
                        cl = pg.get("CL", "")
                        m = re.search(r"--model\s+([^\s]+)", cl)
                        if m:
                            return m.group(1).strip()
                    for rg in bf.header.get("RG", []):
                        ds = rg.get("DS", "")
                        m = re.search(r"(?:basecall_model_version_id|model)=([^\s,;]+)", ds)
                        if m:
                            return m.group(1).strip()
            except Exception:
                continue
    except Exception:
        pass

    return None


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

    reference_registry = None
    if args.reference_files:
        reference_path = args.reference_files.resolve()
        if not reference_path.is_file():
            raise SystemExit(f"Reference registry does not exist: {reference_path}")
        reference_registry = read_yaml(reference_path)

    site_values: dict[str, Any] = {}
    if args.site_config:
        site_config = args.site_config.resolve()
        if not site_config.is_file():
            raise SystemExit(f"Site config does not exist: {site_config}")
        site_values = read_yaml(site_config)

    # Where this site keeps its reference tree. Declared once, in the site
    # config, and used to resolve every tree-relative path in the registry, so
    # the registry says what an asset is called and the site says where the
    # tree lives.
    reference_root = str(
        site_values.get("REFERENCE_ROOT", config.get("REFERENCE_ROOT", "")) or ""
    ).strip()

    # The registry may refer to the site's own variables, so that it can write
    # reference_root: "{{REFERENCE_ROOT}}" instead of repeating the location.
    if reference_registry:
        registry_context = {
            str(key): str(value)
            for key, value in {**config, **site_values}.items()
            if isinstance(value, (str, int, float, bool))
        }
        reference_registry = expand_templates(reference_registry, registry_context)
        reference_root = str(
            reference_registry.get("reference_root", "") or reference_root
        ).strip()

    # base config defaults < registry inventory < site config. The registry is
    # where a reference path is declared; the site config can still pin one
    # without taking over the whole inventory.
    if reference_registry:
        config = deep_merge(
            config, registry_variables(reference_registry, reference_root)
        )

    if site_values:
        config = deep_merge(config, site_values)

    # An explicit --reference-build wins over the config default: a single
    # nagini run uses one build for every sample, so switching assay (bcr_abl
    # needs hg19; wgs/adaptive need hg38) means switching this per invocation.
    reference_build = getattr(args, "reference_build", "") or ""
    if reference_build:
        available = sorted((config.get("references") or {}))
        if reference_build not in available:
            raise SystemExit(
                f"--reference-build {reference_build!r} is not defined in config['references']. "
                f"Available: {', '.join(available)}."
            )
        config["reference_build"] = reference_build

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
        if getattr(args, "check_references", True):
            _verify = bool(getattr(args, "verify_references", False))
            validate_reference_registry(
                reference_registry,
                verify_checksums=_verify,
                strict=_verify,
                default_root=reference_root,
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
        if getattr(args, "dorado_demux_both_ends", False):
            dorado_settings["demux_both_ends"] = True
        if getattr(args, "dorado_demux_extra", ""):
            dorado_settings["demux_extra"] = args.dorado_demux_extra
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

    if getattr(args, "merged_fastq", None) is not None:
        config["merged_fastq"] = bool(args.merged_fastq)

    read_filter_overrides = {}
    if getattr(args, "filter_reads", None) is not None:
        read_filter_overrides["enabled"] = bool(args.filter_reads)
    if getattr(args, "filter_qv", ""):
        try:
            read_filter_overrides["min_qv"] = float(args.filter_qv)
        except ValueError:
            read_filter_overrides["min_qv"] = args.filter_qv
    if getattr(args, "filter_len", ""):
        read_filter_overrides["length_range"] = str(args.filter_len)
    if getattr(args, "filter_binary", ""):
        read_filter_overrides["binary"] = str(args.filter_binary)
    if read_filter_overrides:
        config.setdefault("read_filter", {}).update(read_filter_overrides)

    runtime_config = analysis_dir / "config.nagini.yaml"

    with runtime_config.open("w") as handle:
        yaml.safe_dump(config, handle, sort_keys=False)
    return runtime_config


def validate_assays(args: argparse.Namespace, units: list[dict[str, str]]) -> None:
    """Reject an unknown assay here rather than midway through the DAG."""

    requested = sorted(
        {
            str(row.get("assay", "") or "").strip().lower()
            for row in units
            if str(row.get("assay", "") or "").strip()
        }
    )
    if not requested:
        return
    config = read_yaml(args.base_config.resolve())
    if args.site_config and args.site_config.is_file():
        config = deep_merge(config, read_yaml(args.site_config.resolve()))
    valid = known_assays(_pipeline_root(args), config)
    unknown = [assay for assay in requested if assay not in valid]
    if unknown:
        raise SystemExit(
            "Unknown assay(s): {}. Known assays: {}. Declare a new one under "
            "config['assays'] in the pipeline or site config.".format(
                ", ".join(unknown), ", ".join(sorted(valid))
            )
        )
    print(f"Assays in this run: {', '.join(requested)}")

    # A per-sample assay column wins over --assay by design, which is easy to
    # forget when the map is auto-detected in the inbox: the flag looks like it
    # took effect and the run quietly does the other workup.
    flag = str(getattr(args, "assay", "") or "").strip().lower()
    if flag and flag not in requested:
        print(
            f"NOTE: --assay {flag} was overridden by the assay column in the "
            f"barcode map; this run uses {', '.join(requested)}. Remove that "
            "column, or set it to the assay you want.",
            file=sys.stderr,
        )


def validate_barcode_kit(args: argparse.Namespace, units: list[dict[str, str]]) -> None:
    """Catch a barcode the kit cannot produce before anything is queued."""

    kit = str(getattr(args, "dorado_kit", "") or "").strip()
    if not kit:
        return
    routing = _routing_module(_pipeline_root(args))
    size = routing.barcode_kit_size(kit)
    if size is None:
        print(
            f"NOTE: barcode count for kit {kit} is unknown; barcodes are not checked.",
            file=sys.stderr,
        )
        return
    barcodes = sorted(
        {
            str(row.get("barcode_name", "") or "").strip()
            for row in units
            if str(row.get("barcode_name", "") or "").strip()
        }
    )
    outside = routing.barcodes_outside_kit(barcodes, kit)
    if outside:
        raise SystemExit(
            "Barcode(s) outside kit {} ({} barcodes): {}. Check the barcode map "
            "or the --barcode-kit value.".format(
                kit, size, ", ".join(name for name, _ in outside)
            )
        )
    if barcodes:
        print(f"{len(barcodes)} barcode(s) checked against kit {kit} ({size} barcodes).")


def warn_about_model_download(args: argparse.Namespace) -> None:
    """A model that has to be fetched fails after the GPU job has started.

    Say so at launch instead, where it costs nothing to fix.
    """

    if getattr(args, "input_type", "") != "minknow":
        return
    model = str(getattr(args, "dorado_model", "") or getattr(args, "basecalling_model", "") or "").strip()
    directory = getattr(args, "dorado_models_directory", None)
    if not directory:
        print(
            f"WARNING: no --dorado-models-directory given, so Dorado will try to "
            f"download {model or 'the basecalling model'} when the job starts. On a "
            "host without outbound network that fails after basecalling is queued.",
            file=sys.stderr,
        )
        return
    directory = Path(directory)
    if model and "/" not in model and not (directory / model).exists():
        available = sorted(child.name for child in directory.iterdir() if child.is_dir())
        print(
            "WARNING: model {!r} is not present in {}. Dorado will try to download it. "
            "Models found there: {}".format(
                model, directory, ", ".join(available[:8]) or "none"
            ),
            file=sys.stderr,
        )


def main() -> None:
    args = parse_args()
    args.input_dir = args.input_dir.resolve()
    if not args.input_dir.is_dir():
        raise SystemExit(f"Input directory does not exist: {args.input_dir}")

    if getattr(args, "init_sample_sheet", None):
        generate_sample_sheet_template(args)
        return

    if not args.analysis_dir:
        raise SystemExit("--analysis-dir is required.")
    args.analysis_dir = args.analysis_dir.resolve()
    args.analysis_dir.mkdir(parents=True, exist_ok=True)
    if not args.input_type:
        raise SystemExit("--input-type is required.")
    if (
        not getattr(args, "basecalling_model", None)
        or str(args.basecalling_model).strip().lower() in ("unknown", "", "none")
        or args.basecalling_model == DEFAULT_BASECALLING_MODEL
    ):
        detected_model = detect_basecalling_model(args.input_dir)
        if detected_model:
            if getattr(args, "basecalling_model", None) != detected_model:
                print(f"Auto-detected basecalling model from input files: {detected_model}")
            args.basecalling_model = detected_model
        elif not getattr(args, "basecalling_model", None):
            args.basecalling_model = DEFAULT_BASECALLING_MODEL

    samples_path = args.analysis_dir / "samples.tsv"
    units_path = args.analysis_dir / "units.tsv"
    reuse_existing = (
        not args.force
        and can_reuse_existing_tables(samples_path, units_path, args.input_type)
    )
    if args.input_type == "fastq":
        if reuse_existing:
            samples, units = read_tsv(samples_path), read_tsv(units_path)
        else:
            samples, units = make_fastq_tables(args)
    elif args.input_type == "minknow":
        if reuse_existing:
            samples, units = read_tsv(samples_path), read_tsv(units_path)
        else:
            samples, units = make_minknow_tables(args)
    elif args.barcode_map or discover_bams_by_barcode(args.input_dir):
        if reuse_existing:
            samples, units = read_tsv(samples_path), read_tsv(units_path)
        else:
            samples, units = make_bam_tables_from_barcode_map(args)
    else:
        if reuse_existing:
            samples, units = read_tsv(samples_path), read_tsv(units_path)
        else:
            samples, units = make_bam_tables(args)

    if args.tumor_content is not None:
        for row in samples:
            row["tumor_content"] = args.tumor_content
    validate_tables(samples, units, allow_deferred_inputs=args.input_type == "minknow")
    validate_assays(args, units)
    validate_barcode_kit(args, units)
    warn_about_model_download(args)

    # The union, not the first row's keys: an optional column such as
    # dna_number may be absent from the first sample and present on a later
    # one, and taking only the first row's keys would drop it for everyone.
    sample_columns = ["sample"]
    for row in samples:
        for column in row:
            if column not in sample_columns:
                sample_columns.append(column)
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
