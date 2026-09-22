#!/usr/bin/env python3
"""Rename samples in a finished Nagini analysis directory.

A run started from placeholder barcode names carries them everywhere: in file
names, in the tables, in BAM read groups and in the sample column of every VCF.
Renaming only the files would produce a ``DNA12345_T.bam`` whose read group
still says ``gustaf1_T``, which is worse for traceability than leaving it
alone, so this renames the contents too and says what it could not reach.

What it touches, in order of how safe the operation is:

* **paths** -- files and directories whose name contains an old id.
* **text** -- TSV, CSV, JSON, YAML, BED, HTML, plain VCF and the like.
* **bgzipped text** -- ``.vcf.gz``/``.bed.gz``, rewritten and re-indexed,
  because rewriting a BGZF file invalidates its ``.tbi``.
* **BAM headers** -- ``@RG`` ``ID``/``SM`` via ``samtools reheader``.

Excel workbooks are deliberately not patched: they are the clinical
deliverable, and rebuilding them from the renamed tables is both safer and
cheap. ``--delete-reports`` removes them so a rerun regenerates them.

Nothing is written without ``--apply``; the default is a dry run.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import re
import shutil
import subprocess
import sys
from pathlib import Path

import yaml


# Extensions rewritten as plain text. Everything else is left alone unless a
# more specific handler below claims it.
TEXT_SUFFIXES = (
    ".bed",
    ".cnr",
    ".cns",
    ".csv",
    ".html",
    ".json",
    ".log",
    ".md",
    ".tsv",
    ".txt",
    ".vcf",
    ".yaml",
    ".yml",
)

# Rewritten through gzip, then re-indexed when an index exists beside them.
GZIPPED_TEXT_SUFFIXES = (".bed.gz", ".csv.gz", ".tsv.gz", ".txt.gz", ".vcf.gz")

ALIGNMENT_SUFFIXES = (".bam", ".cram")

# Never rewritten as text: a rename would corrupt them, and they are either
# regenerated or carried along unchanged.
SKIP_SUFFIXES = (".bai", ".crai", ".csi", ".tbi", ".arrow", ".parquet", ".png", ".pdf", ".snf")

REPORT_SUFFIXES = (".xlsx",)


def read_id_map(path: Path, genome: str = "") -> dict[str, str]:
    """Read ``old-id,new-id`` pairs, with or without a header row."""

    if not path.is_file():
        raise SystemExit(f"Id map does not exist: {path}")

    with path.open(newline="") as handle:
        sample = handle.read(4096)
        handle.seek(0)
        try:
            dialect = csv.Sniffer().sniff(sample, delimiters=",\t;")
        except csv.Error:
            dialect = csv.excel_tab if path.suffix.lower() in {".tsv", ".tab"} else csv.excel
        rows = [row for row in csv.reader(handle, dialect) if row]

    if not rows:
        raise SystemExit(f"Id map is empty: {path}")

    # A header is optional: skip the first row only when it names the columns
    # rather than holding an id pair.
    first = [str(field).strip().lower() for field in rows[0][:2]]
    if first and first[0] in {"old", "old_id", "old-id", "from", "sample", "old_sample"}:
        rows = rows[1:]

    mapping: dict[str, str] = {}
    for number, row in enumerate(rows, start=1):
        if len(row) < 2:
            raise SystemExit(f"{path}: row {number} has {len(row)} field(s), need old-id,new-id")
        old, new = str(row[0]).strip(), str(row[1]).strip()
        if not old or not new:
            raise SystemExit(f"{path}: row {number} has an empty id")
        if genome:
            new = f"{new}_{genome}"
        if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]*", new):
            raise SystemExit(
                f"{path}: new id {new!r} is not usable in a file name; "
                "use letters, digits, dot, dash or underscore."
            )
        if old in mapping:
            raise SystemExit(f"{path}: {old} appears more than once")
        mapping[old] = new

    collisions = [new for new in set(mapping.values()) if list(mapping.values()).count(new) > 1]
    if collisions:
        raise SystemExit(
            "Several samples map to the same new id: " + ", ".join(sorted(collisions))
        )
    overlapping = sorted(set(mapping) & set(mapping.values()))
    if overlapping:
        raise SystemExit(
            "An id is both an old and a new name, which makes the result depend "
            "on the order of the rewrite: " + ", ".join(overlapping)
        )
    return mapping


def build_pattern(mapping: dict[str, str]) -> re.Pattern:
    """One regex matching any old id as a whole token.

    The boundaries matter: without them ``gustaf1`` also matches inside
    ``gustaf10``. Longest first so a shared prefix cannot win.
    """

    alternatives = "|".join(
        re.escape(old) for old in sorted(mapping, key=len, reverse=True)
    )
    return re.compile(rf"(?<![A-Za-z0-9]){alternatives}(?![A-Za-z0-9])")


def substitute(text: str, pattern: re.Pattern, mapping: dict[str, str]) -> str:
    return pattern.sub(lambda match: mapping[match.group(0)], text)


def genome_from_config(analysis_dir: Path) -> str:
    """The reference build this run used, from its own config snapshot."""

    config_path = analysis_dir / "config.nagini.yaml"
    if not config_path.is_file():
        return ""
    try:
        config = yaml.safe_load(config_path.read_text()) or {}
    except yaml.YAMLError:
        return ""
    reference = config.get("reference", {})
    build = reference.get("build") if isinstance(reference, dict) else ""
    return str(build or config.get("reference_build", "") or "").strip()


def classify(path: Path) -> str:
    """Which rewrite, if any, applies to this file."""

    name = path.name.lower()
    if name.endswith(SKIP_SUFFIXES):
        return "skip"
    if name.endswith(REPORT_SUFFIXES):
        return "report"
    if name.endswith(ALIGNMENT_SUFFIXES):
        return "alignment"
    if name.endswith(GZIPPED_TEXT_SUFFIXES):
        return "gzip"
    if name.endswith(TEXT_SUFFIXES):
        return "text"
    return "other"


def rewrite_text(path: Path, pattern, mapping, apply_changes: bool) -> bool:
    try:
        original = path.read_text()
    except (OSError, UnicodeDecodeError):
        return False
    updated = substitute(original, pattern, mapping)
    if updated == original:
        return False
    if apply_changes:
        path.write_text(updated)
    return True


def rewrite_gzip(path: Path, pattern, mapping, apply_changes: bool) -> bool:
    """Rewrite a gzipped text file, keeping BGZF where bgzip is available.

    A plain ``gzip`` rewrite of a ``.vcf.gz`` would still read fine but would
    no longer be block-compressed, so tabix could not index it.
    """

    try:
        with gzip.open(path, "rt") as handle:
            original = handle.read()
    except (OSError, UnicodeDecodeError, EOFError):
        return False
    updated = substitute(original, pattern, mapping)
    if updated == original:
        return False
    if not apply_changes:
        return True

    temporary = path.with_name(path.name + ".rename.tmp")
    if shutil.which("bgzip"):
        with temporary.open("wb") as out:
            process = subprocess.run(
                ["bgzip", "-c"], input=updated.encode(), stdout=out, check=False
            )
        if process.returncode != 0:
            temporary.unlink(missing_ok=True)
            raise SystemExit(f"bgzip failed while rewriting {path}")
    else:
        with gzip.open(temporary, "wt") as handle:
            handle.write(updated)
    temporary.replace(path)

    # The old index describes byte offsets that no longer exist.
    for index_suffix in (".tbi", ".csi"):
        index = path.with_name(path.name + index_suffix)
        if index.exists():
            index.unlink()
            if shutil.which("tabix"):
                subprocess.run(["tabix", "-f", "-p", "vcf", str(path)], check=False)
    return True


def reheader_alignment(path: Path, pattern, mapping, apply_changes: bool) -> bool:
    """Rewrite @RG ID/SM in a BAM or CRAM header with samtools."""

    if not shutil.which("samtools"):
        return False
    header = subprocess.run(
        ["samtools", "view", "-H", str(path)], capture_output=True, text=True, check=False
    )
    if header.returncode != 0:
        return False
    updated = substitute(header.stdout, pattern, mapping)
    if updated == header.stdout:
        return False
    if not apply_changes:
        return True

    header_file = path.with_name(path.name + ".header.tmp")
    header_file.write_text(updated)
    temporary = path.with_name(path.name + ".rename.tmp")
    with temporary.open("wb") as out:
        process = subprocess.run(
            ["samtools", "reheader", str(header_file), str(path)], stdout=out, check=False
        )
    header_file.unlink(missing_ok=True)
    if process.returncode != 0:
        temporary.unlink(missing_ok=True)
        raise SystemExit(f"samtools reheader failed for {path}")
    temporary.replace(path)
    # Reheadering shifts every record, so the old index is wrong.
    for index_suffix in (".bai", ".csi", ".crai"):
        index = path.with_name(path.name + index_suffix)
        if index.exists():
            index.unlink()
            subprocess.run(["samtools", "index", str(path)], check=False)
    return True


def rename_paths(root: Path, pattern, mapping, apply_changes: bool) -> list[tuple[Path, Path]]:
    """Rename every file and directory whose name carries an old id.

    Bottom-up, so renaming a directory cannot invalidate the paths of the
    entries inside it that were queued first.
    """

    renames = []
    for path in sorted(root.rglob("*"), key=lambda p: len(p.parts), reverse=True):
        new_name = substitute(path.name, pattern, mapping)
        if new_name == path.name:
            continue
        target = path.with_name(new_name)
        if target.exists():
            raise SystemExit(f"Refusing to overwrite an existing path: {target}")
        renames.append((path, target))
        if apply_changes:
            path.rename(target)
    return renames


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Rename samples in a finished Nagini analysis directory.",
    )
    parser.add_argument("--analysis-dir", type=Path, required=True)
    parser.add_argument(
        "--map",
        dest="id_map",
        type=Path,
        required=True,
        help="CSV/TSV of old-id,new-id, with or without a header row.",
    )
    parser.add_argument(
        "--genome",
        default="",
        help="Append a reference build to every new id, for example hg38. "
        "Use 'auto' to take it from the run's config.nagini.yaml.",
    )
    parser.add_argument(
        "--apply",
        action="store_true",
        help="Perform the changes. Without it nothing is written.",
    )
    parser.add_argument(
        "--no-contents",
        dest="contents",
        action="store_false",
        help="Rename paths only, leaving the old id inside the files.",
    )
    parser.add_argument(
        "--bams",
        action="store_true",
        help="Also rewrite @RG in BAM/CRAM headers. Needs samtools, rewrites "
        "the whole file and re-indexes it.",
    )
    parser.add_argument(
        "--delete-reports",
        action="store_true",
        help="Delete the Excel workbooks so a rerun regenerates them from the "
        "renamed tables instead of carrying a patched binary.",
    )
    args = parser.parse_args()

    analysis_dir = args.analysis_dir.expanduser().resolve()
    if not analysis_dir.is_dir():
        raise SystemExit(f"Analysis directory does not exist: {analysis_dir}")

    genome = args.genome.strip()
    if genome.lower() == "auto":
        genome = genome_from_config(analysis_dir)
        if not genome:
            raise SystemExit(
                "--genome auto could not read a reference build from "
                f"{analysis_dir / 'config.nagini.yaml'}; pass the build explicitly."
            )
        print(f"Reference build from the run config: {genome}")

    mapping = read_id_map(args.id_map.expanduser().resolve(), genome)
    pattern = build_pattern(mapping)

    print(f"Analysis directory: {analysis_dir}")
    for old, new in sorted(mapping.items()):
        print(f"  {old} -> {new}")
    if not args.apply:
        print("\nDry run. Nothing is written; add --apply to perform the changes.\n")

    changed_files: list[Path] = []
    pending_reports: list[Path] = []
    unreachable: list[Path] = []

    if args.contents:
        for path in sorted(analysis_dir.rglob("*")):
            if not path.is_file():
                continue
            kind = classify(path)
            if kind == "skip":
                continue
            if kind == "report":
                try:
                    stale = any(old.encode() in path.read_bytes() for old in mapping)
                except OSError:
                    stale = False
                if stale:
                    pending_reports.append(path)
                    # Delete here, before the paths move: collecting them now
                    # and unlinking after rename_paths looked for names that no
                    # longer exist.
                    if args.delete_reports and args.apply:
                        path.unlink()
                continue
            if kind == "alignment":
                if args.bams:
                    if reheader_alignment(path, pattern, mapping, args.apply):
                        changed_files.append(path)
                elif shutil.which("samtools"):
                    header = subprocess.run(
                        ["samtools", "view", "-H", str(path)],
                        capture_output=True,
                        text=True,
                        check=False,
                    )
                    if header.returncode == 0 and pattern.search(header.stdout):
                        unreachable.append(path)
                continue
            if kind == "gzip":
                if rewrite_gzip(path, pattern, mapping, args.apply):
                    changed_files.append(path)
                continue
            if kind == "text":
                if rewrite_text(path, pattern, mapping, args.apply):
                    changed_files.append(path)

    renames = rename_paths(analysis_dir, pattern, mapping, args.apply)

    print(f"Paths renamed: {len(renames)}")
    for source, target in renames[:20]:
        print(f"  {source.relative_to(analysis_dir)} -> {target.name}")
    if len(renames) > 20:
        print(f"  ... {len(renames) - 20} more")
    print(f"Files with rewritten contents: {len(changed_files)}")

    if pending_reports:
        # Keep the warnings after the summary they qualify.
        sys.stdout.flush()
        print(
            f"\n{len(pending_reports)} Excel workbook(s) still contain an old id. "
            "They are not patched:",
            file=sys.stderr,
        )
        for path in pending_reports:
            print(f"  {path.relative_to(analysis_dir)}", file=sys.stderr)
        if args.delete_reports and args.apply:
            print("Deleted; rerun the pipeline to rebuild them.", file=sys.stderr)
        else:
            print(
                "Delete them and rerun the pipeline to rebuild them from the "
                "renamed tables, or pass --delete-reports.",
                file=sys.stderr,
            )

    if unreachable:
        sys.stdout.flush()
        print(
            f"\n{len(unreachable)} alignment file(s) still carry an old id in "
            "their @RG header. Pass --bams to rewrite them:",
            file=sys.stderr,
        )
        for path in unreachable[:10]:
            print(f"  {path.relative_to(analysis_dir)}", file=sys.stderr)

    if args.apply:
        print("\nDone.")


if __name__ == "__main__":
    main()
