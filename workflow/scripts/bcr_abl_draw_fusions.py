#!/usr/bin/env python3
"""Draw the BCR-ABL1 fusion plots (PDF and PNG) for one sample.

Arriba's draw_fusions.R is tried first.  It needs the gene annotation, and it
draws more panels when cytobands and protein domains are available.  When it
cannot run or leaves an empty file, draw_bcr_abl_fusion.py produces the
diagram from the fusion summary instead, so that the two deliverables always
exist.  The step fails only if a deliverable is still empty at the end.

This used to be a 150-line shell block in bcr_abl.smk.  Standard library only,
because it runs inside the Arriba container.
"""

from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Callable, Dict, List, Optional, Sequence

HERE = Path(__file__).resolve().parent

# Where draw_fusions.R lives in the images and installs this has run in, after
# the copy bundled next to this script and anything on PATH.
DRAW_SCRIPT_CANDIDATES = (
    "/arriba_v2.4.0/draw_fusions.R",
    "/usr/local/bin/draw_fusions.R",
)
CYTOBAND_CANDIDATES = (
    "/arriba_v2.4.0/database/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv",
    "/arriba_v2.4.0/database/cytobands_hg38_GRCh38_v2.4.0.tsv",
    "/var/lib/arriba/database/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv",
    "/var/lib/arriba/database/cytobands_hg38_GRCh38_v2.4.0.tsv",
    "/data/references/arriba/cytobands_hg19_hs37d5_GRCh37_v2.4.0.tsv",
    "/data/references/arriba/cytobands_hg38_GRCh38_v2.4.0.tsv",
)
PROTEIN_DOMAIN_CANDIDATES = (
    "/arriba_v2.4.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3",
    "/arriba_v2.4.0/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.tsv",
    "/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.gff3",
    "/arriba_v2.4.0/database/protein_domains_hg38_GRCh38_v2.4.0.tsv",
    "/var/lib/arriba/database/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3",
    "/var/lib/arriba/database/protein_domains_hg38_GRCh38_v2.4.0.gff3",
    "/data/references/arriba/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3",
    "/data/references/arriba/protein_domains_hg38_GRCh38_v2.4.0.gff3",
)


# Arriba names its database files after the assembly. The image ships both
# builds side by side, so a candidate list ordered hg19-first silently gave an
# hg38 run hg19 cytobands and protein domains.
BUILD_TOKENS = {
    "hg19": ("hg19", "grch37", "hs37d5"),
    "hg38": ("hg38", "grch38"),
}


def candidates_for_build(candidates: Sequence[str], build: str) -> List[str]:
    """The candidates named for this assembly, in their original order.

    An unknown build gets none: no reference file beats one from the wrong
    assembly, because the plot then simply omits that panel.
    """

    tokens = BUILD_TOKENS.get(str(build).strip().lower())
    if not tokens:
        return []
    return [c for c in candidates if any(t in os.path.basename(c).lower() for t in tokens)]


def is_nonempty(path: str) -> bool:
    return os.path.isfile(path) and os.path.getsize(path) > 0


def first_existing(
    configured: str,
    candidates: Sequence[str],
    isfile: Callable[[str], bool] = os.path.isfile,
) -> str:
    """The configured file if it exists, else the first candidate that does."""

    if configured and isfile(configured):
        return configured
    for candidate in candidates:
        if isfile(candidate):
            return candidate
    return ""


def find_draw_command(
    bundled: Path = HERE / "draw_fusions.r",
    isfile: Callable[[str], bool] = os.path.isfile,
    which: Callable[[str], Optional[str]] = shutil.which,
) -> List[str]:
    """The command that runs draw_fusions.R, or an empty list if none exists."""

    if isfile(str(bundled)):
        return ["Rscript", str(bundled)]
    on_path = which("draw_fusions.R")
    if on_path:
        return ["draw_fusions.R"]
    for candidate in DRAW_SCRIPT_CANDIDATES:
        if isfile(candidate):
            return [candidate]
    return []


def choose_panels(cytobands: str, protein_domains: str) -> str:
    """The plot panels draw_fusions.R can fill with what is available."""

    if cytobands and protein_domains:
        return "fusion,circos,domains,readcounts"
    if cytobands:
        return "fusion,circos,readcounts"
    if protein_domains:
        return "fusion,domains,readcounts"
    return "fusion,readcounts"


def build_draw_args(
    draw_cmd: Sequence[str],
    fusions: str,
    bam: str,
    output: str,
    gtf: str,
    sample: str,
    cytobands: str,
    protein_domains: str,
    extra: str = "",
) -> List[str]:
    args = list(draw_cmd) + [
        f"--fusions={fusions}",
        f"--alignments={bam}",
        f"--output={output}",
        f"--annotation={gtf}",
        f"--sampleName={sample}",
        "--transcriptSelection=canonical",
        f"--plotPanels={choose_panels(cytobands, protein_domains)}",
    ]
    if cytobands:
        args.append(f"--cytobands={cytobands}")
    if protein_domains:
        args.append(f"--proteinDomains={protein_domains}")
    args.extend(shlex.split(extra or ""))
    return args


def primary_fusion_table(source: str, destination: str) -> None:
    """Header plus the first (most supported) fusion, for the single-fusion PNG."""

    with open(source) as handle:
        lines = [line for _, line in zip(range(2), handle)]
    with open(destination, "w") as handle:
        handle.writelines(lines)


def draw_fusion_plots(
    args: argparse.Namespace,
    run: Callable[..., "subprocess.CompletedProcess"] = subprocess.run,
    which: Callable[[str], Optional[str]] = shutil.which,
    isfile: Callable[[str], bool] = os.path.isfile,
    bundled: Path = HERE / "draw_fusions.r",
    fallback_script: Path = HERE / "draw_bcr_abl_fusion.py",
    python: str = sys.executable,
    log=print,
) -> int:
    """Produce args.output_pdf and args.output_png.  Returns a process exit code."""

    log(f"Running BCR-ABL1 fusion plot for {args.sample}...")

    draw_cmd = find_draw_command(bundled, isfile, which)
    cytobands = first_existing(args.cytobands, candidates_for_build(CYTOBAND_CANDIDATES, args.build), isfile)
    protein_domains = first_existing(
        args.protein_domains, candidates_for_build(PROTEIN_DOMAIN_CANDIDATES, args.build), isfile
    )
    log(f"Draw script: {' '.join(draw_cmd) if draw_cmd else '(none found)'}")
    log(f"Cytobands: {cytobands or '(none)'}")
    log(f"Protein domains: {protein_domains or '(none)'}")
    log(f"GTF: {args.gtf}")

    def attempt(label: str, command: List[str], tolerate: bool = True) -> bool:
        log(f"Executing {label}...")
        result = run(command)
        if result.returncode != 0:
            log(f"{label} exited with status {result.returncode}"
                f"{'; continuing' if tolerate else ''}")
        return result.returncode == 0

    if draw_cmd and args.gtf and isfile(args.gtf):
        common = dict(
            bam=args.bam, gtf=args.gtf, sample=args.sample,
            cytobands=cytobands, protein_domains=protein_domains, extra=args.extra,
        )
        attempt(
            "Arriba draw_fusions.R (all fusions to PDF)",
            build_draw_args(draw_cmd, fusions=args.fusions, output=args.output_pdf, **common),
        )
        with tempfile.TemporaryDirectory(prefix="bcr_abl_primary_") as tmp:
            primary = os.path.join(tmp, "primary_fusion.tsv")
            primary_fusion_table(args.fusions, primary)
            attempt(
                "Arriba draw_fusions.R for PNG (primary fusion)",
                build_draw_args(draw_cmd, fusions=primary, output=args.output_png, **common),
            )

        # The R device writes the whole PDF; take its first page as the PNG so
        # that the picture in the Excel report is the primary fusion.
        if which("pdftoppm") and is_nonempty(args.output_pdf):
            with tempfile.TemporaryDirectory(prefix="bcr_abl_page_") as tmp:
                prefix = os.path.join(tmp, "page")
                run(["pdftoppm", "-png", "-r", "150", "-f", "1", "-l", "1", args.output_pdf, prefix])
                for produced in sorted(Path(tmp).glob("page*.png")):
                    if produced.stat().st_size > 0:
                        shutil.move(str(produced), args.output_png)
                        break
    else:
        log("draw_fusions.R or the GTF is missing; using the Python diagram.")

    def fallback(output: str, tolerate: bool) -> None:
        command = [
            python, str(fallback_script),
            "--sample", args.sample,
            "--fusions", args.fusion_summary,
            "--mutations", args.mutations,
            "--phasing", args.phasing,
            "--exon-coverage", args.exon_coverage,
            "--build", args.build,
            "--output", output,
        ]
        log(f"Generating BCR-ABL1 diagram with draw_bcr_abl_fusion.py -> {output}")
        result = run(command)
        if result.returncode != 0 and not tolerate:
            raise RuntimeError(f"draw_bcr_abl_fusion.py failed for {output} (status {result.returncode})")

    if not is_nonempty(args.output_pdf):
        fallback(args.output_pdf, tolerate=False)
    if not is_nonempty(args.output_png):
        fallback(args.output_png, tolerate=True)

    missing = [p for p in (args.output_pdf, args.output_png) if not is_nonempty(p)]
    if missing:
        log("No plot was produced: " + ", ".join(missing))
        return 1
    return 0


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--sample", required=True)
    parser.add_argument("--fusions", required=True, help="Arriba-format fusion table")
    parser.add_argument("--fusion-summary", required=True, help="fusion_summary.tsv for the Python diagram")
    parser.add_argument("--mutations", default="")
    parser.add_argument("--phasing", default="")
    parser.add_argument("--exon-coverage", default="")
    parser.add_argument("--bam", required=True)
    parser.add_argument("--gtf", default="")
    parser.add_argument("--cytobands", default="")
    parser.add_argument("--protein-domains", default="")
    parser.add_argument("--build", default="hg19")
    parser.add_argument("--extra", default="", help="Extra draw_fusions.R arguments")
    parser.add_argument("--output-pdf", required=True)
    parser.add_argument("--output-png", required=True)
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    for output in (args.output_pdf, args.output_png):
        Path(output).parent.mkdir(parents=True, exist_ok=True)
    try:
        return draw_fusion_plots(args, log=lambda message: print(message, flush=True))
    except RuntimeError as error:
        print(f"Error: {error}", flush=True)
        return 1


if __name__ == "__main__":
    sys.exit(main())
