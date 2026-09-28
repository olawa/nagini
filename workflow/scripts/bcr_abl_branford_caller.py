#!/usr/bin/env python3
"""Force-call BCR-ABL1 kinase domain resistance mutations (Branford panel).

High-depth pileup calling at the Branford ABL1 kinase-domain positions,
excluding reference skips (introns in cDNA spliced alignments), computing VAFs
and annotating TKI resistance profiles.

Every interrogated position is written to the all-sites TSV with its status and
filter_reason, so the output answers "was T315I excluded?" and not only "was
T315I found?". The mutations TSV and VCF carry DETECTED and SUBCLONAL only.

Artifact filters, and what each one can and cannot see
-----------------------------------------------------
At ~46 000x, Nanopore's substitution noise puts hundreds of alternate reads on
every position, so a bare VAF threshold calls noise. The filters below narrow
that, but they are not a substitute for a per-position error baseline built
from negative samples.

SNR (--min-snr)
    ad_alt divided by the largest of the two *other* non-reference bases.
    Catches noise that spreads evenly across all three substitutions.
    It is blind to the opposite and more common case: a context-specific
    Nanopore error produces the same substitution every time, leaving the other
    two bases clean, which yields a *high* SNR. At 46 000x a recurrent 1.4 %
    error scores SNR 25.8x against a true 1.4 % subclone's 29.3x -- the two are
    not separable by this statistic. Only a background panel separates them.

STRAND_BIAS (--min-minor-strand-ratio, --max-sor)
    Requires both a one-sided alt distribution and a departure from the
    reference's own strand balance. Requiring both matters on amplicons: when
    the whole amplicon is strand-skewed the reference is one-sided too, and a
    real variant inherits that skew without being an artifact.

LOW_BASEQ (--min-base-quality, --min-mean-baseq)
    --min-base-quality drops individual base calls during the pileup, so it
    also changes the depth that VAF is computed against. depth_unfiltered and
    ad_alt_unfiltered are reported alongside, so a downstream comparison can
    use counts that do not depend on this setting.

DIRTY_CONTEXT (--min-kmer-rate, --kmer-window)
    Fraction of alt reads whose flanking bases match the reference. Targets
    homopolymer and indel-adjacent errors. Deviations that coincide with
    another panel position are ignored, because a compound mutation in cis is
    not a dirty context -- 25 of the 52 Branford positions sit within 3 bp of
    another one, so without that exemption this filter would reject exactly
    the compound mutations the phaser exists to resolve.

Thresholds are conventional starting points, not values derived from these
data. Calibrate them against the per-position background once it exists.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pysam

BASES = ("A", "C", "G", "T")

# Built-in TKI resistance lookup table
DEFAULT_TKI_PROFILES = {
    "T315I": {"Imatinib": "Resistant", "Dasatinib": "Resistant", "Nilotinib": "Resistant", "Bosutinib": "Resistant", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "T315A": {"Imatinib": "Resistant", "Dasatinib": "Resistant", "Nilotinib": "Sensitive", "Bosutinib": "Resistant", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "Y253H": {"Imatinib": "Resistant", "Dasatinib": "Intermediate", "Nilotinib": "Resistant", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "Y253F": {"Imatinib": "Resistant", "Dasatinib": "Intermediate", "Nilotinib": "Resistant", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "E255K": {"Imatinib": "Resistant", "Dasatinib": "Intermediate", "Nilotinib": "Resistant", "Bosutinib": "Intermediate", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "E255V": {"Imatinib": "Resistant", "Dasatinib": "Intermediate", "Nilotinib": "Resistant", "Bosutinib": "Intermediate", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "F359V": {"Imatinib": "Resistant", "Dasatinib": "Sensitive", "Nilotinib": "Resistant", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "F359I": {"Imatinib": "Resistant", "Dasatinib": "Sensitive", "Nilotinib": "Resistant", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "F359C": {"Imatinib": "Resistant", "Dasatinib": "Sensitive", "Nilotinib": "Resistant", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "G250E": {"Imatinib": "Resistant", "Dasatinib": "Intermediate", "Nilotinib": "Intermediate", "Bosutinib": "Intermediate", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "F317L": {"Imatinib": "Intermediate", "Dasatinib": "Resistant", "Nilotinib": "Sensitive", "Bosutinib": "Intermediate", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "F317V": {"Imatinib": "Intermediate", "Dasatinib": "Resistant", "Nilotinib": "Sensitive", "Bosutinib": "Intermediate", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "F317I": {"Imatinib": "Intermediate", "Dasatinib": "Resistant", "Nilotinib": "Sensitive", "Bosutinib": "Intermediate", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "F317C": {"Imatinib": "Intermediate", "Dasatinib": "Resistant", "Nilotinib": "Sensitive", "Bosutinib": "Intermediate", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "V299L": {"Imatinib": "Sensitive", "Dasatinib": "Resistant", "Nilotinib": "Sensitive", "Bosutinib": "Resistant", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "M244V": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "L248V": {"Imatinib": "Resistant", "Dasatinib": "Intermediate", "Nilotinib": "Intermediate", "Bosutinib": "Intermediate", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "Q252H": {"Imatinib": "Intermediate", "Dasatinib": "Intermediate", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "Q252R": {"Imatinib": "Intermediate", "Dasatinib": "Intermediate", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "D276G": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "E279K": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "F311I": {"Imatinib": "Intermediate", "Dasatinib": "Intermediate", "Nilotinib": "Intermediate", "Bosutinib": "Intermediate", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "F311L": {"Imatinib": "Intermediate", "Dasatinib": "Intermediate", "Nilotinib": "Intermediate", "Bosutinib": "Intermediate", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "M351T": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "E355G": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "V379I": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "L387M": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Intermediate", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "L387F": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Intermediate", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "H396R": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "H396P": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "E450K": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "E450G": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
    "E459K": {"Imatinib": "Intermediate", "Dasatinib": "Sensitive", "Nilotinib": "Sensitive", "Bosutinib": "Sensitive", "Ponatinib": "Sensitive", "Asciminib": "Sensitive"},
}

AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}


def short_aa_notation(aa_change: str) -> str:
    """Convert 'p.Thr315Ile' to 'T315I'."""
    text = aa_change.strip()
    if text.startswith("p."):
        text = text[2:]
    for three, one in AA3_TO_AA1.items():
        text = text.replace(three, one)
    return text


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path, help="Input aligned BAM file")
    parser.add_argument("--bed", required=True, type=Path, help="Branford positions BED file")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--output-all", required=True, type=Path, help="Output TSV with all interrogated sites")
    parser.add_argument("--output-mutations", required=True, type=Path, help="Output TSV with detected mutations")
    parser.add_argument("--output-vcf", type=Path, default=None, help="Optional output VCF file")
    parser.add_argument("--features-json", type=Path, default=None, help="Optional BCR-ABL features JSON database")
    parser.add_argument("--reference-fasta", type=Path, default=None, help="Optional reference FASTA for k-mer validation")
    parser.add_argument("--min-af", type=float, default=0.01, help="Minimum AF threshold for mutation calling (default 0.01 / 1 percent)")
    parser.add_argument("--subclonal-af", type=float, default=0.002, help="Subclonal AF threshold (default 0.002 / 0.2 percent)")
    parser.add_argument("--min-depth", type=int, default=50, help="Minimum total coverage depth (default 50)")
    parser.add_argument("--min-alt-reads", type=int, default=5, help="Minimum alternate supporting reads (default 5)")
    parser.add_argument("--max-depth", type=int, default=1_000_000, help="Maximum pileup depth")
    parser.add_argument("--min-base-quality", type=int, default=12, help="Minimum base quality (default 12)")
    parser.add_argument("--min-mapping-quality", type=int, default=1, help="Minimum mapping quality (default 1)")
    parser.add_argument("--min-snr", type=float, default=3.0, help="Minimum signal-to-noise ratio, ad_alt / largest other non-reference base. Blind to context-specific errors that favour one substitution; see the module docstring (default 3.0)")
    parser.add_argument("--min-minor-strand-ratio", type=float, default=0.05, help="Minimum fraction on minor strand when ad_alt >= 20 (default 0.05)")
    parser.add_argument("--max-sor", type=float, default=3.0, help="Maximum Strand Odds Ratio; alt strand balance is only called biased when it also departs from the reference's own balance (default 3.0)")
    parser.add_argument("--min-mean-baseq", type=float, default=14.0, help="Minimum mean base quality for Alt reads (default 14.0)")
    parser.add_argument("--min-kmer-rate", type=float, default=0.60, help="Minimum kmer reference clean rate (default 0.60 / 60 percent)")
    parser.add_argument("--kmer-window", type=int, default=3, help="Flanking kmer window size in bp (default 3)")
    parser.add_argument("--known-variants", type=Path, default=None, help="Optional BED/VCF of known polymorphisms whose flanking mismatches should not count as dirty context. Panel positions are always exempt.")
    parser.add_argument("--include-artifacts", action="store_true", help="Include ARTIFACT calls in mutations TSV")
    parser.add_argument("--threads", type=int, default=1, help="Number of worker threads for BAM decompression")
    return parser.parse_args()



def read_branford_bed(bed_path: Path) -> List[Dict[str, Any]]:
    rows = []
    with bed_path.open() as handle:
        for line_idx, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t")
            if len(fields) < 7:
                fields = text.split()
            if len(fields) < 7:
                continue

            chrom, start, end, name, ref, alt, aa_change = fields[:7]
            start_pos = int(start)
            pos_1based = start_pos + 1
            short_name = short_aa_notation(aa_change)

            rows.append({
                "line": line_idx,
                "chrom": chrom,
                "pos": pos_1based,
                "name": name,
                "ref": ref.upper(),
                "alt": alt.upper(),
                "aa_change": aa_change,
                "short_mutation": short_name,
            })
    return rows


def load_exempt_positions(path: Optional[Path]) -> Dict[str, set]:
    """Read extra positions whose flanking mismatches should not count as dirty.

    Accepts a BED (chrom, start0, end) or a VCF (chrom, pos). Used for known
    polymorphisms: a germline SNP beside a panel position is a real base, not a
    context error, and without this it would reject that position in every
    sample carrying the SNP.
    """

    exempt: Dict[str, set] = {}
    if not path or not Path(path).is_file():
        return exempt
    for line in Path(path).read_text().splitlines():
        text = line.strip()
        if not text or text.startswith(("#", "track", "browser")):
            continue
        fields = text.split("\t")
        if len(fields) < 2:
            continue
        chrom = fields[0]
        try:
            if str(path).lower().endswith((".vcf", ".vcf.gz")):
                zero_based = int(fields[1]) - 1
            else:
                zero_based = int(fields[1])
        except ValueError:
            continue
        exempt.setdefault(chrom, set()).add(zero_based)
    return exempt


def flank_masks(
    chrom: str, pos0: int, window: int, exempt: Dict[str, set]
) -> Tuple[Tuple[bool, ...], Tuple[bool, ...]]:
    """Which flanking offsets may differ from the reference without penalty."""

    sites = exempt.get(chrom, set())
    left = tuple((pos0 - window + i) in sites for i in range(window))
    right = tuple((pos0 + 1 + i) in sites for i in range(window))
    return left, right


def context_is_clean(
    read_seq: str, ref_seq: str, mask: Tuple[bool, ...]
) -> bool:
    """Compare a flank to the reference, forgiving exempt offsets."""

    if len(read_seq) != len(ref_seq):
        return False
    for index, (read_base, ref_base) in enumerate(zip(read_seq, ref_seq)):
        if read_base == ref_base:
            continue
        if index < len(mask) and mask[index]:
            continue
        return False
    return True


def count_positions_batch(
    bam: pysam.AlignmentFile,
    chrom: str,
    positions_0based: set[int],
    max_depth: int,
    min_baseq: int,
    min_mapq: int,
    ref_seqs_by_pos: Optional[Dict[Tuple[str, int], Tuple[str, str]]] = None,
    orig_chrom: Optional[str] = None,
) -> Dict[int, Dict[str, Any]]:
    """Extract forward/reverse base counts, base qualities, and k-mer context across coordinates."""
    results: Dict[int, Dict[str, Any]] = {}
    if not positions_0based:
        return results

    for pos_0based in positions_0based:
        counts: Dict[str, Any] = {f"{base}_{strand}": 0 for base in BASES for strand in ("fwd", "rev")}
        counts.update({
            "refskip": 0,
            "deletion": 0,
            "low_baseq": 0,
            "base_qual_sum": {b: 0 for b in BASES},
            "base_qual_count": {b: 0 for b in BASES},
            "kmer_clean_count": {b: 0 for b in BASES},
            "kmer_total_count": {b: 0 for b in BASES},
            # Counted before --min-base-quality is applied. The quality filter
            # changes the denominator that VAF is computed against, so anything
            # comparing these calls against another platform needs counts that
            # do not move when the threshold is tuned.
            "unfiltered": {b: 0 for b in BASES},
        })
        results[pos_0based] = counts

    min_pos = min(positions_0based)
    max_pos = max(positions_0based) + 1

    for column in bam.pileup(
        chrom,
        min_pos,
        max_pos,
        truncate=True,
        stepper="nofilter",
        max_depth=max_depth,
        min_base_quality=0,
        ignore_overlaps=False,
    ):
        ref_pos = column.reference_pos
        if ref_pos not in positions_0based:
            continue

        counts = results[ref_pos]
        lookup_key = (orig_chrom or chrom, ref_pos)
        has_kmer = bool(ref_seqs_by_pos and lookup_key in ref_seqs_by_pos)
        left_mask: Tuple[bool, ...] = ()
        right_mask: Tuple[bool, ...] = ()
        if has_kmer and ref_seqs_by_pos is not None:
            ref_left, ref_right, left_mask, right_mask = ref_seqs_by_pos[lookup_key]
            k_len = len(ref_left)
        else:
            ref_left, ref_right, k_len = "", "", 0

        for read in column.pileups:
            aln = read.alignment
            if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail:
                continue
            if (aln.mapping_quality or 0) < min_mapq:
                continue

            # In cDNA spliced alignments, skip introns
            if read.is_refskip:
                counts["refskip"] += 1
                continue

            if read.is_del or read.query_position is None:
                counts["deletion"] += 1
                continue

            bq = aln.query_qualities[read.query_position] if aln.query_qualities is not None else 30

            base = aln.query_sequence[read.query_position].upper()
            if base not in BASES:
                continue
            counts["unfiltered"][base] += 1

            if min_baseq > 0 and bq < min_baseq:
                counts["low_baseq"] += 1
                continue

            strand = "rev" if aln.is_reverse else "fwd"
            counts[f"{base}_{strand}"] += 1
            counts["base_qual_sum"][base] += bq
            counts["base_qual_count"][base] += 1

            # Flanking k-mer reference context check
            if k_len > 0:
                qpos = read.query_position
                qseq = aln.query_sequence
                if qpos >= k_len and (qpos + 1 + k_len) <= len(qseq):
                    counts["kmer_total_count"][base] += 1
                    read_left = qseq[qpos - k_len : qpos].upper()
                    read_right = qseq[qpos + 1 : qpos + 1 + k_len].upper()
                    if context_is_clean(read_left, ref_left, left_mask) and context_is_clean(
                        read_right, ref_right, right_mask
                    ):
                        counts["kmer_clean_count"][base] += 1

    return results


def count_position(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos_1based: int,
    max_depth: int,
    min_baseq: int,
    min_mapq: int,
) -> Dict[str, Any]:
    """Extract forward/reverse base counts at a single coordinate (wrapper for count_positions_batch)."""
    target_0based = pos_1based - 1
    batch_res = count_positions_batch(bam, chrom, {target_0based}, max_depth, min_baseq, min_mapq)
    return batch_res.get(target_0based, {f"{b}_{s}": 0 for b in BASES for s in ("fwd", "rev")})



def main() -> None:
    args = parse_args()
    if not args.bam.is_file():
        sys.exit(f"BAM file not found: {args.bam}")

    sites = read_branford_bed(args.bed)
    if not sites:
        sys.exit(f"No sites parsed from {args.bed}")

    # Load custom TKI profile if provided
    tki_profiles = dict(DEFAULT_TKI_PROFILES)
    if args.features_json and args.features_json.is_file():
        try:
            with args.features_json.open() as h:
                features = json.load(h)
                if "tki_profiles" in features:
                    tki_profiles.update(features["tki_profiles"])
        except Exception as e:
            sys.stderr.write(f"Warning: could not load TKI profiles from {args.features_json}: {e}\n")

    # Pre-extract reference k-mer context if reference fasta is provided
    # A mismatch in the flank is only evidence of a context error when it is
    # not itself something we expect to vary. Panel positions are the obvious
    # case: 25 of the 52 Branford sites lie within 3 bp of another one, so a
    # compound mutation in cis would otherwise mark both of its own positions
    # dirty and be rejected before the phaser ever sees it.
    exempt_positions: Dict[str, set] = load_exempt_positions(
        getattr(args, "known_variants", None)
    )
    for site in sites:
        exempt_positions.setdefault(site["chrom"], set()).add(site["pos"] - 1)

    ref_seqs_by_pos: Dict[Tuple[str, int], Tuple[str, str, Tuple[bool, ...], Tuple[bool, ...]]] = {}
    if args.reference_fasta and args.reference_fasta.is_file():
        try:
            ref_fasta = pysam.FastaFile(str(args.reference_fasta))
            fasta_refs = set(ref_fasta.references)
            for site in sites:
                c = site["chrom"]
                if c not in fasta_refs:
                    alt_c = c[3:] if c.startswith("chr") else f"chr{c}"
                    if alt_c in fasta_refs:
                        c = alt_c
                if c in fasta_refs:
                    p0 = site["pos"] - 1
                    kw = args.kmer_window
                    c_len = ref_fasta.get_reference_length(c)
                    if p0 >= kw and (p0 + 1 + kw) <= c_len:
                        left_seq = ref_fasta.fetch(c, p0 - kw, p0).upper()
                        right_seq = ref_fasta.fetch(c, p0 + 1, p0 + 1 + kw).upper()
                        left_mask, right_mask = flank_masks(
                            site["chrom"], p0, kw, exempt_positions
                        )
                        ref_seqs_by_pos[(site["chrom"], p0)] = (
                            left_seq, right_seq, left_mask, right_mask,
                        )
            ref_fasta.close()
        except Exception as e:
            sys.stderr.write(f"Warning: could not load reference fasta {args.reference_fasta} for k-mer check: {e}\n")

    bam = pysam.AlignmentFile(str(args.bam), "rb", threads=getattr(args, "threads", 1))

    # Contig renaming check (e.g. chr9 vs 9)
    bam_contigs = set(bam.references)

    # Group all target sites by resolved BAM contig for single-pass scanning
    chrom_positions: Dict[str, set[int]] = {}
    chrom_alias_map: Dict[str, str] = {}
    for site in sites:
        chrom = site["chrom"]
        orig_chrom = chrom
        if chrom not in bam_contigs:
            alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            if alt_chrom in bam_contigs:
                chrom = alt_chrom
        chrom_alias_map[chrom] = orig_chrom
        if chrom not in chrom_positions:
            chrom_positions[chrom] = set()
        chrom_positions[chrom].add(site["pos"] - 1)

    cache: Dict[Tuple[str, int], Dict[str, Any]] = {}
    try:
        for chrom, pos_set in chrom_positions.items():
            if chrom in bam_contigs:
                orig_chrom = chrom_alias_map.get(chrom, chrom)
                batch_counts = count_positions_batch(
                    bam, chrom, pos_set, args.max_depth, args.min_base_quality, args.min_mapping_quality,
                    ref_seqs_by_pos=ref_seqs_by_pos, orig_chrom=orig_chrom
                )
                for pos_0based, counts in batch_counts.items():
                    cache[(chrom, pos_0based + 1)] = counts
            else:
                for pos_0based in pos_set:
                    c: Dict[str, Any] = {f"{b}_{s}": 0 for b in BASES for s in ("fwd", "rev")}
                    c.update({
                        "refskip": 0, "deletion": 0, "low_baseq": 0,
                        "base_qual_sum": {b: 0 for b in BASES},
                        "base_qual_count": {b: 0 for b in BASES},
                        "kmer_clean_count": {b: 0 for b in BASES},
                        "kmer_total_count": {b: 0 for b in BASES},
                    })
                    cache[(chrom, pos_0based + 1)] = c

        all_results: List[Dict[str, Any]] = []
        detected_mutations: List[Dict[str, Any]] = []

        for site in sites:
            chrom = site["chrom"]
            if chrom not in bam_contigs:
                alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
                if alt_chrom in bam_contigs:
                    chrom = alt_chrom

            pos = site["pos"]
            ref = site["ref"]
            alt = site["alt"]
            short_mut = site["short_mutation"]

            counts = cache.get((chrom, pos), {f"{b}_{s}": 0 for b in BASES for s in ("fwd", "rev")})
            per_base = {b: counts.get(f"{b}_fwd", 0) + counts.get(f"{b}_rev", 0) for b in BASES}
            depth = sum(per_base.values())

            ad_ref = per_base.get(ref, 0)
            ad_alt = per_base.get(alt, 0)
            vaf = (ad_alt / depth) if depth > 0 else 0.0

            unfiltered = counts.get("unfiltered", {})
            depth_unfiltered = sum(unfiltered.values())
            ad_alt_unfiltered = unfiltered.get(alt, 0)
            vaf_unfiltered = (
                (ad_alt_unfiltered / depth_unfiltered) if depth_unfiltered > 0 else 0.0
            )

            ad_ref_fwd = counts.get(f"{ref}_fwd", 0)
            ad_ref_rev = counts.get(f"{ref}_rev", 0)
            ad_alt_fwd = counts.get(f"{alt}_fwd", 0)
            ad_alt_rev = counts.get(f"{alt}_rev", 0)

            # Strand balance & Strand Odds Ratio (SOR)
            alt_minor_strand_ratio = round(min(ad_alt_fwd, ad_alt_rev) / max(ad_alt, 1), 3)
            ref_has_both_strands = (ad_ref_fwd > 0 and ad_ref_rev > 0)
            r1 = (ad_alt_fwd + 1) * (ad_ref_rev + 1)
            r2 = (ad_alt_rev + 1) * (ad_ref_fwd + 1)
            sor = round(math.log(max(r1 / r2, r2 / r1)), 2)

            # Mean base quality for Ref and Alt
            bq_sum_alt = counts.get("base_qual_sum", {}).get(alt, 0)
            bq_cnt_alt = counts.get("base_qual_count", {}).get(alt, 0)
            mean_baseq_alt = round(bq_sum_alt / max(bq_cnt_alt, 1), 1) if bq_cnt_alt > 0 else 0.0

            bq_sum_ref = counts.get("base_qual_sum", {}).get(ref, 0)
            bq_cnt_ref = counts.get("base_qual_count", {}).get(ref, 0)
            mean_baseq_ref = round(bq_sum_ref / max(bq_cnt_ref, 1), 1) if bq_cnt_ref > 0 else 0.0

            # Flanking k-mer reference clean percentage
            k_clean = counts.get("kmer_clean_count", {}).get(alt, 0)
            k_tot = counts.get("kmer_total_count", {}).get(alt, 0)
            kmer_clean_percent = round((k_clean / k_tot) * 100.0, 1) if k_tot > 0 else 100.0

            # Error baseline from non-candidate alternative bases
            other_bases = {b: per_base[b] for b in BASES if b != ref and b != alt}
            ad_other = sum(other_bases.values())
            ad_other_max = max(other_bases.values()) if other_bases else 0
            vaf_other_percent = round((ad_other / depth) * 100.0, 3) if depth > 0 else 0.0
            other_bases_str = ",".join(f"{b}:{c}" for b, c in sorted(other_bases.items()))

            # Signal-to-Noise Ratio (SNR)
            if ad_other_max == 0:
                snr = float(ad_alt) if ad_alt > 0 else 0.0
            else:
                snr = round(ad_alt / max(ad_other_max, 1), 2)

            # Quality filtering and status classification
            if depth < args.min_depth:
                status = "LOW_DEPTH"
                filter_reason = f"DEPTH<{args.min_depth}"
            elif ad_alt < args.min_alt_reads:
                status = "WILDTYPE"
                filter_reason = "NO_ALT"
            else:
                fails = []
                if snr < args.min_snr:
                    fails.append(f"LOW_SNR({snr:.1f}x<{args.min_snr:.1f}x)")
                # One-sidedness alone is not bias: on a strand-skewed amplicon the
                # reference is one-sided too, and a real variant inherits that skew.
                # SOR compares the alt balance against the reference's own, so both
                # have to agree before the call is rejected.
                if (
                    ref_has_both_strands
                    and ad_alt >= 20
                    and alt_minor_strand_ratio < args.min_minor_strand_ratio
                    and sor > args.max_sor
                ):
                    fails.append(
                        f"STRAND_BIAS(ratio={alt_minor_strand_ratio:.3f}<{args.min_minor_strand_ratio:.2f},"
                        f"SOR={sor:.2f}>{args.max_sor:.1f})"
                    )
                if mean_baseq_alt < args.min_mean_baseq:
                    fails.append(f"LOW_BASEQ({mean_baseq_alt:.1f}<{args.min_mean_baseq:.1f})")
                if k_tot > 0 and kmer_clean_percent < (args.min_kmer_rate * 100.0):
                    fails.append(f"DIRTY_CONTEXT({kmer_clean_percent:.1f}%<{(args.min_kmer_rate*100):.0f}%)")

                if fails:
                    status = "ARTIFACT"
                    filter_reason = ";".join(fails)
                elif vaf >= args.min_af:
                    status = "DETECTED"
                    filter_reason = "PASS"
                elif vaf >= args.subclonal_af:
                    status = "SUBCLONAL"
                    filter_reason = "PASS"
                else:
                    status = "WILDTYPE"
                    filter_reason = "LOW_AF"

            # TKI responses
            tki_info = tki_profiles.get(short_mut, {})
            imatinib = tki_info.get("Imatinib", "Unknown")
            dasatinib = tki_info.get("Dasatinib", "Unknown")
            nilotinib = tki_info.get("Nilotinib", "Unknown")
            bosutinib = tki_info.get("Bosutinib", "Unknown")
            ponatinib = tki_info.get("Ponatinib", "Unknown")
            asciminib = tki_info.get("Asciminib", "Unknown")

            record = {
                "sample": args.sample,
                "chrom": site["chrom"],
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "mutation": short_mut,
                "aa_change": site["aa_change"],
                "name": site["name"],
                "depth": depth,
                "ad_ref": ad_ref,
                "ad_alt": ad_alt,
                "vaf_percent": round(vaf * 100.0, 3),
                "depth_unfiltered": depth_unfiltered,
                "ad_alt_unfiltered": ad_alt_unfiltered,
                "vaf_unfiltered_percent": round(vaf_unfiltered * 100.0, 3),
                "status": status,
                "filter_reason": filter_reason,
                "snr": snr,
                "ad_other": ad_other,
                "vaf_other_percent": vaf_other_percent,
                "other_bases": other_bases_str,
                "alt_minor_strand_ratio": alt_minor_strand_ratio,
                "sor": sor,
                "mean_baseq_ref": mean_baseq_ref,
                "mean_baseq_alt": mean_baseq_alt,
                "kmer_clean_percent": kmer_clean_percent,
                "imatinib": imatinib,
                "dasatinib": dasatinib,
                "nilotinib": nilotinib,
                "bosutinib": bosutinib,
                "ponatinib": ponatinib,
                "asciminib": asciminib,
                "ad_ref_fwd": ad_ref_fwd,
                "ad_ref_rev": ad_ref_rev,
                "ad_alt_fwd": ad_alt_fwd,
                "ad_alt_rev": ad_alt_rev,
                "n_refskip": counts.get("refskip", 0),
                "n_deletion": counts.get("deletion", 0),
            }
            all_results.append(record)

            if status in ("DETECTED", "SUBCLONAL"):
                detected_mutations.append(record)
            elif args.include_artifacts and status == "ARTIFACT":
                detected_mutations.append(record)

    finally:
        bam.close()

    # Sort results
    detected_mutations.sort(key=lambda r: (r["status"] != "DETECTED", -r["vaf_percent"]))

    # Output headers
    headers = [
        "sample", "chrom", "pos", "ref", "alt", "mutation", "aa_change", "name",
        "depth", "ad_ref", "ad_alt", "vaf_percent", "status", "filter_reason",
        "depth_unfiltered", "ad_alt_unfiltered", "vaf_unfiltered_percent",
        "snr", "ad_other", "vaf_other_percent", "other_bases",
        "alt_minor_strand_ratio", "sor",
        "mean_baseq_ref", "mean_baseq_alt", "kmer_clean_percent",
        "imatinib", "dasatinib", "nilotinib", "bosutinib", "ponatinib", "asciminib",
        "ad_ref_fwd", "ad_ref_rev", "ad_alt_fwd", "ad_alt_rev", "n_refskip", "n_deletion",
    ]

    # Write All Sites TSV
    args.output_all.parent.mkdir(parents=True, exist_ok=True)
    with args.output_all.open("w") as handle:
        handle.write(f"# sample={args.sample}\n")
        handle.write(f"# total_sites={len(all_results)}\n")
        handle.write(f"# detected_mutations={len([m for m in detected_mutations if m['status'] == 'DETECTED'])}\n")
        handle.write(f"# subclonal_mutations={len([m for m in detected_mutations if m['status'] == 'SUBCLONAL'])}\n")
        handle.write("\t".join(headers) + "\n")
        for rec in all_results:
            handle.write("\t".join(str(rec.get(h, "")) for h in headers) + "\n")

    # Write Detected Mutations TSV
    args.output_mutations.parent.mkdir(parents=True, exist_ok=True)
    with args.output_mutations.open("w") as handle:
        handle.write(f"# sample={args.sample}\n")
        handle.write(f"# min_af={args.min_af}\n")
        handle.write(f"# subclonal_af={args.subclonal_af}\n")
        handle.write(f"# min_snr={args.min_snr}\n")
        handle.write(f"# max_sor={args.max_sor}\n")
        handle.write(f"# min_depth={args.min_depth}\n")
        handle.write(f"# min_base_quality={args.min_base_quality}\n")
        handle.write(f"# min_mean_baseq={args.min_mean_baseq}\n")
        handle.write(f"# min_kmer_rate={args.min_kmer_rate} kmer_window={args.kmer_window}\n")
        handle.write("\t".join(headers) + "\n")
        for rec in detected_mutations:
            handle.write("\t".join(str(rec.get(h, "")) for h in headers) + "\n")

    # Write VCF if requested
    if args.output_vcf:
        args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
        with args.output_vcf.open("w") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write("##source=nagini_bcr_abl_branford_caller\n")
            handle.write('##INFO=<ID=AA,Number=1,Type=String,Description="Amino acid change">\n')
            handle.write('##INFO=<ID=MUT,Number=1,Type=String,Description="Short mutation notation">\n')
            handle.write('##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS transcript annotation">\n')
            handle.write('##INFO=<ID=SNR,Number=1,Type=Float,Description="Signal-to-noise ratio vs non-candidate error bases">\n')
            handle.write('##INFO=<ID=SOR,Number=1,Type=Float,Description="Strand Odds Ratio">\n')
            handle.write('##INFO=<ID=BQ_ALT,Number=1,Type=Float,Description="Mean base quality of ALT reads">\n')
            handle.write('##FILTER=<ID=PASS,Description="All quality filters passed">\n')
            handle.write('##FILTER=<ID=ARTIFACT,Description="Failed error baseline, strand bias, or base quality filters">\n')
            handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">\n')
            handle.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (REF, ALT)">\n')
            handle.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">\n')
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + args.sample + "\n")
            for rec in detected_mutations:
                af_dec = rec["vaf_percent"] / 100.0
                info_field = f"AA={rec['aa_change']};MUT={rec['mutation']};HGVS={rec['name']};SNR={rec['snr']};SOR={rec['sor']};BQ_ALT={rec['mean_baseq_alt']}"
                gt = "0/1" if rec["status"] in ("DETECTED", "SUBCLONAL") else "./."
                vcf_filter = "PASS" if rec["status"] in ("DETECTED", "SUBCLONAL") else "ARTIFACT"
                format_val = f"{gt}:{rec['depth']}:{rec['ad_ref']},{rec['ad_alt']}:{af_dec:.4f}"
                handle.write(f"{rec['chrom']}\t{rec['pos']}\t{rec['mutation']}\t{rec['ref']}\t{rec['alt']}\t.\t{vcf_filter}\t{info_field}\tGT:DP:AD:AF\t{format_val}\n")

    print(
        f"BCR-ABL1 Branford Caller [{args.sample}]: evaluated {len(sites)} resistance sites; "
        f"found {len(detected_mutations)} candidate mutations.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

