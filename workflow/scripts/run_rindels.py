#!/usr/bin/env python3
"""Run the Rindels binary for nagini's three Rindels steps.

``cnv``       whole-genome CNV calling; leaves a placeholder when Rindels skips
              the ASCAT segmentation or the PDF report
``ideogram``  whole-genome ideogram; uses the optional sidecars that exist
``bcr-abl``   targeted BCR-ABL1 calling; indexes the VCF when it can

The Rindels binary is a host binary, so none of this may depend on tools that
only the containers have.  Standard library only.  These steps were shell
blocks in cnv.smk and bcr_abl.smk, each with its own copy of the binary check.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Callable, List, Optional, Sequence

ASCAT_PLACEHOLDER = (
    "# ASCAT-style Heuristic Copy Number Segments\n"
    "# NOTE: Rindels skipped segmentation for this sample (too few heterozygous sites)\n"
    "# Columns: chrom, start, end, CN_A, CN_B, BAF, logR, n_sites, state, n_het, n_hom, n_total\n"
)


class RindelsError(RuntimeError):
    pass


def check_binary(binary: str, which: Callable[[str], Optional[str]] = shutil.which) -> None:
    """Fail clearly if the binary is neither on PATH nor an executable file."""

    if which(binary) or (os.path.isfile(binary) and os.access(binary, os.X_OK)):
        return
    raise RindelsError(f"Rindels binary not found at {binary}")


def nonempty(path: str) -> bool:
    return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0


def run_checked(run: Callable, command: Sequence[str]) -> None:
    print("+ " + " ".join(shlex.quote(str(part)) for part in command), flush=True)
    result = run(list(command))
    if result.returncode != 0:
        raise RindelsError(f"Rindels exited with status {result.returncode}")


# ---------------------------------------------------------------------------
# cnv
# ---------------------------------------------------------------------------
def cnv_command(args: argparse.Namespace) -> List[str]:
    command = [
        args.binary,
        "-b", args.bam,
        "-r", args.ref,
        "--cnv-output", args.bins,
        "--cnv-bin-size", str(args.bin_size),
        "-o", args.vcf,
    ]
    if args.report:
        command += ["--report", "--report-path", args.pdf]
    command += ["-t", str(args.threads)]
    command += shlex.split(args.annotation or "")
    command += shlex.split(args.extra or "")
    return command


def run_cnv(args: argparse.Namespace, run: Callable = subprocess.run, which=shutil.which) -> None:
    check_binary(args.binary, which)
    run_checked(run, cnv_command(args))
    ensure_compressed_vcf(args.vcf, run, which)

    # Rindels skips the ASCAT segmentation when it finds fewer than 100
    # heterozygous sites, a legitimate outcome for a shallow or genomically
    # quiet sample.  Leave a header-only file so the run reports no segments
    # instead of dying on a missing output; cnv_summary treats an empty file
    # as zero segments.
    if not os.path.isfile(args.ascat):
        print("ASCAT segmentation skipped by Rindels (too few heterozygous sites); see the log.", flush=True)
        Path(args.ascat).write_text(ASCAT_PLACEHOLDER)
    if not os.path.isfile(args.pdf):
        print("Rindels unified PDF report not produced or switched off; leaving a placeholder.", flush=True)
        Path(args.pdf).touch()


# ---------------------------------------------------------------------------
# ideogram
# ---------------------------------------------------------------------------
def ideogram_command(args: argparse.Namespace, say: Callable[[str], None] = print) -> List[str]:
    command = [args.binary, "ideogram"]

    # Prefer the table the PDF report renders from; fall back to the BED.
    if nonempty(args.parquet):
        command += ["--cnv", args.parquet]
    else:
        say(f"No CNV Parquet at {args.parquet}; falling back to the binned BED.")
        command += ["--bed", args.bins]

    if nonempty(args.baf):
        command += ["--baf-bed", args.baf]
    else:
        say(f"No BAF track at {args.baf}; drawing density only.")

    if nonempty(args.centromeres):
        command += ["--centromeres", args.centromeres]
    else:
        say("No centromere BED; Rindels decides from the data's own lengths.")

    # Gaps are registered per build and usually absent: a gapless assembly such
    # as chm13 has none to draw, and Rindels no longer substitutes GRCh38's when
    # the lengths say the assembly is something else.
    if nonempty(args.gaps):
        command += ["--gaps", args.gaps]

    command += [
        "--bin-size", str(args.bin_size),
        "--res", str(args.bin_size // 1000),
        "--width", str(args.width),
        "--height", str(args.height),
    ]
    command += shlex.split(args.extra or "")
    command += ["-o", args.output]
    return command


def run_ideogram(args: argparse.Namespace, run: Callable = subprocess.run, which=shutil.which) -> None:
    check_binary(args.binary, which)
    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    run_checked(run, ideogram_command(args, say=lambda m: print(m, flush=True)))


# ---------------------------------------------------------------------------
# bcr-abl
# ---------------------------------------------------------------------------
def bcr_abl_command(args: argparse.Namespace) -> List[str]:
    command = [args.binary, "-b", args.bam, "-r", args.ref, "-o", args.vcf]
    command += shlex.split(args.extra or "")
    command += [
        "--bcr-abl-genome", args.genome,
        "--write-haplotypes", args.haplotypes,
        "-t", str(args.threads),
    ]
    return command


def ensure_compressed_vcf(vcf: str, run: Callable, which: Callable[[str], Optional[str]]) -> None:
    """Ensure that if a .vcf.gz was requested, it is non-empty and properly compressed.

    If Rindels falls back to writing an uncompressed .vcf because bgzip is missing or failed,
    it can leave behind a 0-byte .vcf.gz while .vcf has the variant records. If so, compress
    the uncompressed .vcf into .vcf.gz.
    """
    if not vcf or not vcf.endswith(".gz"):
        return

    uncompressed = vcf[:-3]
    vcf_empty = (not os.path.exists(vcf)) or os.path.getsize(vcf) == 0
    raw_exists = os.path.exists(uncompressed) and os.path.getsize(uncompressed) > 0

    if vcf_empty and raw_exists:
        print(f"Rindels left uncompressed VCF at {uncompressed}; compressing to {vcf}...", flush=True)
        if which("bcftools"):
            res = run(["bcftools", "view", "-Oz", "-o", vcf, uncompressed])
            if getattr(res, "returncode", 0) == 0 and os.path.exists(vcf) and os.path.getsize(vcf) > 0:
                return
        if which("bgzip"):
            try:
                res = run(["bgzip", "-f", uncompressed])
                bgz_out = f"{uncompressed}.gz"
                if os.path.exists(bgz_out) and bgz_out != vcf:
                    os.replace(bgz_out, vcf)
                if os.path.exists(vcf) and os.path.getsize(vcf) > 0:
                    return
            except Exception:
                pass
        import gzip
        with open(uncompressed, "rb") as f_in, gzip.open(vcf, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)


def index_vcf(vcf: str, run: Callable, which: Callable[[str], Optional[str]]) -> None:
    """Index the VCF if bcftools or tabix is available.  The index is optional
    and nothing downstream needs it, so a failure is a warning."""

    if not nonempty(vcf):
        return
    if which("bcftools"):
        command = ["bcftools", "index", "-t", vcf]
    elif which("tabix"):
        command = ["tabix", "-p", "vcf", vcf]
    else:
        print("WARNING: neither bcftools nor tabix is available; the VCF is not indexed.", flush=True)
        return
    result = run(command)
    if result.returncode != 0:
        print(f"WARNING: indexing {vcf} failed (status {result.returncode}); continuing without an index.", flush=True)


def run_bcr_abl(args: argparse.Namespace, run: Callable = subprocess.run, which=shutil.which) -> None:
    check_binary(args.binary, which)
    run_checked(run, bcr_abl_command(args))
    ensure_compressed_vcf(args.vcf, run, which)
    index_vcf(args.vcf, run, which)


# ---------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="mode", required=True)

    def common(p):
        p.add_argument("--binary", required=True)
        p.add_argument("--extra", default="")

    cnv = sub.add_parser("cnv")
    common(cnv)
    cnv.add_argument("--bam", required=True)
    cnv.add_argument("--ref", required=True)
    cnv.add_argument("--bins", required=True, help="CNV bin BED to write")
    cnv.add_argument("--bin-size", type=int, required=True)
    cnv.add_argument("--vcf", required=True)
    cnv.add_argument("--ascat", required=True, help="ASCAT segment BED Rindels writes next to --vcf")
    cnv.add_argument("--pdf", required=True)
    cnv.add_argument("--report", action="store_true", help="Ask Rindels for the unified PDF report")
    cnv.add_argument("--annotation", default="", help="Annotation arguments, e.g. '--gtf x --annotate'")
    cnv.add_argument("--threads", type=int, default=1)

    ideo = sub.add_parser("ideogram")
    common(ideo)
    ideo.add_argument("--bins", required=True)
    ideo.add_argument("--parquet", default="")
    ideo.add_argument("--baf", default="")
    ideo.add_argument("--centromeres", default="")
    ideo.add_argument("--gaps", default="")
    ideo.add_argument("--bin-size", type=int, required=True)
    ideo.add_argument("--width", type=int, required=True)
    ideo.add_argument("--height", type=int, required=True)
    ideo.add_argument("--output", required=True)

    bcr = sub.add_parser("bcr-abl")
    common(bcr)
    bcr.add_argument("--bam", required=True)
    bcr.add_argument("--ref", required=True)
    bcr.add_argument("--vcf", required=True)
    bcr.add_argument("--haplotypes", required=True)
    bcr.add_argument("--genome", required=True)
    bcr.add_argument("--threads", type=int, default=1)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    handlers = {"cnv": run_cnv, "ideogram": run_ideogram, "bcr-abl": run_bcr_abl}
    try:
        handlers[args.mode](args)
    except RindelsError as error:
        print(f"Error: {error}", flush=True)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
