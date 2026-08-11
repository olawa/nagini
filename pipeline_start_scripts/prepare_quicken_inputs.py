#!/usr/bin/env python3
"""Create Hydra-style input tables and a run-local Quicken config.

Hydra-Genetics can create ONT tables from BAM read groups. FASTQ input does
not have enough metadata for Hydra to infer the same table, so this helper
creates the equivalent table using conservative defaults and explicit CLI
options.
"""

from __future__ import annotations

import argparse
import csv
import re
import subprocess
from collections import defaultdict
from pathlib import Path

import yaml


FASTQ_SUFFIXES = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
BASE_UNIT_COLUMNS = [
    "sample",
    "type",
    "platform",
    "machine",
    "processing_unit",
    "barcode",
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
    parser.add_argument(
        "--input-type", choices=("fastq", "unmapped_bam", "aligned_bam"), required=True
    )
    parser.add_argument("--sample-type", choices=("N", "T", "R"), default="N")
    parser.add_argument("--sequence-id", default="")
    parser.add_argument("--sample-regex", default=DEFAULT_SAMPLE_REGEX)
    parser.add_argument("--barcode", default="NB")
    parser.add_argument("--run-id", default="")
    parser.add_argument("--basecalling-model", default="unknown")
    parser.add_argument("--methylation", choices=("Yes", "No"), default="No")
    parser.add_argument("--tumor-content", default=None)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle, delimiter="\t")]


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


def run_hydra_create_input_files(args: argparse.Namespace, analysis_dir: Path) -> None:
    command = [
        "hydra-genetics",
        "create-input-files",
        "--platform",
        "ONT",
        "--sample-type",
        args.sample_type,
        "--directory",
        str(args.input_dir),
        "--default-barcode",
        args.barcode,
    ]
    if args.force:
        command.append("--force")
    subprocess.run(command, cwd=analysis_dir, check=True)


def normalize_bam_tables(
    args: argparse.Namespace, analysis_dir: Path
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    samples_path = analysis_dir / "samples.tsv"
    units_path = analysis_dir / "units.tsv"
    if not samples_path.exists() or not units_path.exists():
        run_hydra_create_input_files(args, analysis_dir)
    if not samples_path.exists() or not units_path.exists():
        raise SystemExit("Hydra did not create samples.tsv and units.tsv.")

    samples = read_tsv(samples_path)
    units = read_tsv(units_path)
    if not units:
        raise SystemExit("units.tsv contains no rows.")
    for row in units:
        row["input_type"] = args.input_type
        row["fastq1"] = row.get("fastq1") or "UNUSED"
        row["fastq2"] = row.get("fastq2") or "UNUSED"
        bam_path = Path(row["bam"]).expanduser()
        if not bam_path.is_absolute():
            bam_path = analysis_dir / bam_path
        row["bam"] = str(bam_path.resolve())
        row.setdefault("machine", "unknown")
        row.setdefault("methylation", args.methylation)
        row.setdefault("run_id", args.run_id or args.sequence_id or args.input_dir.name)
        row.setdefault("basecalling_model", args.basecalling_model)
    return samples, units


def validate_tables(samples: list[dict[str, str]], units: list[dict[str, str]]) -> None:
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
        if not input_path or input_path == "UNUSED" or not Path(input_path).exists():
            raise SystemExit(f"Input file does not exist for unit {key}: {input_path}")


def resolve_config_path(value: object, base_dir: Path) -> object:
    if not isinstance(value, str) or not value or value.startswith(("docker://", "http://", "https://")):
        return value
    path = Path(value).expanduser()
    return str(path if path.is_absolute() else (base_dir / path).resolve())


def make_runtime_config(args: argparse.Namespace, analysis_dir: Path) -> Path:
    base_config = args.base_config.resolve()
    with base_config.open() as handle:
        config = yaml.safe_load(handle) or {}
    # Quicken's checked-in configs use paths relative to the repository root,
    # even though the YAML files themselves live below config/.
    base_dir = args.pipeline_root.resolve()
    config["samples"] = str((analysis_dir / "samples.tsv").resolve())
    config["units"] = str((analysis_dir / "units.tsv").resolve())
    config["output"] = str(Path(resolve_config_path(config.get("output", "config/output_files.yaml"), base_dir)))
    config["resources"] = str(Path(resolve_config_path(config.get("resources", "config/resources.yaml"), base_dir)))
    config["hydra_local_path"] = str(Path(resolve_config_path(config.get("hydra_local_path", "vendor"), base_dir)))

    for reference in config.get("references", {}).values():
        if not isinstance(reference, dict):
            continue
        for key in ("fasta", "fai", "dict", "design_bed", "sv_annotation_bed"):
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

    runtime_config = analysis_dir / "config.quicken.yaml"
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
    else:
        samples, units = normalize_bam_tables(args, args.analysis_dir)

    if args.tumor_content is not None:
        for row in samples:
            row["tumor_content"] = args.tumor_content
    validate_tables(samples, units)

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
