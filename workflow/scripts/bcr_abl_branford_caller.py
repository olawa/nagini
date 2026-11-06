#!/usr/bin/env python3
"""Force-call BCR-ABL1 kinase domain resistance mutations (Branford panel).

High-depth base counting at the Branford ABL1 kinase-domain positions, from
one pass over the reads (a CIGAR walk per read, not a pileup), excluding
reference skips (introns in cDNA spliced alignments), computing VAFs and
annotating TKI resistance profiles.

Every interrogated position is written to the all-sites TSV with its status and
filter_reason, so the output answers "was T315I excluded?" and not only "was
T315I found?". The mutations TSV and VCF retain detected/subclonal evidence,
including explicitly unresolved codons; unresolved VCF records have a named
filter and a missing genotype. --include-artifacts also retains failed calls.

Artifact filters, and what each one can and cannot see
-----------------------------------------------------
At ~46 000x, Nanopore's substitution noise puts hundreds of alternate reads on
every position, so a bare VAF threshold calls noise. The filters below narrow
that, but they are not a substitute for a per-position error baseline built
from negative samples.

SNR (--min-snr)
    ad_alt divided by the largest observed non-candidate base count. Other
    panel substitutions are exempted only when they have independent support
    and exceed the weakest observed non-reference base by --min-snr, and pass
    the base-quality, strand-bias, and context checks. This
    keeps a balanced spectrum across all three alternatives conservative.
    Catches noise that spreads evenly across the remaining substitutions.
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

READ QUALITY (--min-read-qs)
    Drops whole reads whose mean quality (the qs value) is below the cut-off,
    before any counting. Off by default. It applies to the aligned reads, so a
    stricter value than the pre-alignment filter can be tried without aligning
    again; reads the pre-alignment filter removed cannot be brought back.

LOW_BASEQ (--min-base-quality, --min-mean-baseq)
    --min-base-quality drops individual base calls during the pileup, so it
    also changes the depth that VAF is computed against. depth_unfiltered and
    ad_alt_unfiltered are reported alongside, so a downstream comparison can
    use counts that do not depend on this setting.

FLANKING CONTEXT (--kmer-window)
    With --reference-fasta, a read whose bases within --kmer-window of the site
    differ from the reference is not counted at all -- not for the reference,
    not for any alternate. Such reads are what a context error looks like, and
    excluding them from both sides keeps the VAF unbiased: true alternate reads
    carry random flank errors as often as reference reads do. They are reported
    as depth_context_excluded / ad_alt_context_excluded, and kmer_clean_percent
    still shows how clean the alternate reads were before exclusion.

    This replaced rejecting a whole call when under 60 % of its alternate reads
    were clean (DIRTY_CONTEXT, --min-kmer-rate). K357E (c.1069A>G) shows why:
    on the reverse strand the G of AAG|AAAAA is read one base late, giving a
    G>A + A>G pair at 0.2-0.3 % in every sample. A true 0.3 % subclone there
    would have had about half its alternate reads dirty and been rejected;
    counting only clean reads keeps the subclone and drops the artefact.

    A flank base that is a real allele is not a context error, and three kinds
    are forgiven -- each only as the specific base, never "any base at that
    position" (D25):

    * another panel position's ALT, when that ALT is itself called in the
      sample. 38 of the 52 Branford positions sit within 3 bp of another one,
      so a compound mutation in cis would otherwise discard its own reads.
      Which ALTs are called is settled in two rounds: a provisional one that
      forgives every panel ALT, then the final one that forgives only those
      the provisional round called. Noise at a neighbouring panel position
      therefore still makes a read dirty.
    * a substitution the sample carries beside the site
      (--flank-variant-min-af, default 1 %, on both strands, with the
      ordinary base-quality check). A germline SNV next to a panel position
      used to make every read dirty (LOW_DEPTH), and a non-panel change in cis
      with a panel mutation hid the mutation. These are listed in the
      flank_alleles_forgiven column and the file header. Only reads whose
      flank is one gapless aligned block are used to find them, so a
      homopolymer indel cannot pose as a substitution.
    * any base at a --known-variants position.

FUSION READS (--fusion-region, --fusion-reads-only)
    A read is a BCR::ABL1 read when it also aligns in the fusion-partner region
    (its SA tag, i.e. the BCR part). With --fusion-reads-only only those reads
    are counted; resistance mutations matter only in the fusion. Reads without
    a BCR part are 3-9 % of kinase-domain reads in the dilution series and carry
    the true variants at the same VAF, but E499E stood at 2.8 % in them against
    0.14 % in fusion reads: they include ABL1 molecules that are not BCR::ABL1.
    Both read sets are reported separately (*_fusion, *_nonfusion columns)
    whatever the setting.

CODONS
    Two variant positions in one codon (two called panel positions, or a called
    panel position and a sample variant from the flank check) make the single-
    site annotations unreliable, and the codon is read as whole haplotypes.
    A haplotype is credible, and can be reported, only when every base that
    differs from the reference is an allele the site calls accepted (a called
    panel ALT or a sample variant), and it is reported only when at least one
    of them is a called panel ALT. Errors at the codon's other bases are noise,
    as at single sites; before D25 they formed reportable haplotypes of their
    own that were also exempt from the noise estimate.

Thresholds are conventional starting points, not values derived from these
data. Calibrate them against the per-position background once it exists.
"""


import argparse
import bisect
import shlex
import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Set

import pysam

from bcr_abl_annotations import load_features, load_tki_profiles
from bcr_abl_read_quality import ReadQualityFilter

BASES = ("A", "C", "G", "T")

AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}

CODON_TABLE = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L", "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*", "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W",
    "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L", "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q", "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M", "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K", "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R",
    "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V", "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E", "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}

CDNA_RE = re.compile(r"c\.\((\d+)([ACGT])>([ACGT])\)")


def short_aa_notation(aa_change: str) -> str:
    """Convert 'p.Thr315Ile' to 'T315I'."""
    text = aa_change.strip()
    if text.startswith("p."):
        text = text[2:]
    for three, one in AA3_TO_AA1.items():
        text = text.replace(three, one)
    return text


def parse_args(argv=None) -> argparse.Namespace:
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
    parser.add_argument("--max-depth", type=int, default=1_000_000, help="No longer used: reads are scanned once rather than piled up, so there is no depth cap. Accepted so existing invocations keep working.")
    parser.add_argument("--min-base-quality", type=int, default=12, help="Minimum base quality (default 12)")
    parser.add_argument("--min-mapping-quality", type=int, default=1, help="Minimum mapping quality (default 1)")
    parser.add_argument("--min-read-qs", type=float, default=0.0, help="Minimum mean read quality (the qs value), applied to the aligned reads; 0 is off (default 0). Reads below the pre-alignment filter are already gone, so this can only be stricter.")
    parser.add_argument("--min-snr", type=float, default=3.0, help="Minimum signal-to-noise ratio; supported panel competitors are exempt only if they exceed the weakest non-reference base by this ratio (default 3.0)")
    parser.add_argument("--min-minor-strand-ratio", type=float, default=0.05, help="Minimum fraction on minor strand when ad_alt >= 20 (default 0.05)")
    parser.add_argument("--max-sor", type=float, default=3.0, help="Maximum Strand Odds Ratio; alt strand balance is only called biased when it also departs from the reference's own balance (default 3.0)")
    parser.add_argument("--min-mean-baseq", type=float, default=14.0, help="Minimum mean base quality for Alt reads (default 14.0)")
    parser.add_argument("--min-kmer-rate", type=float, default=0.60, help="No longer used: reads with a dirty flanking context are now excluded from the counts instead of a call being rejected below this clean rate. Accepted so existing invocations keep working.")
    parser.add_argument("--kmer-window", type=int, default=3, help="Flanking kmer window size in bp (default 3)")
    parser.add_argument("--known-variants", type=Path, default=None, help="Optional BED/VCF of known polymorphisms; any base at these positions is forgiven in the flanking-context check.")
    parser.add_argument("--flank-variant-min-af", type=float, default=0.01, help="A substitution beside a panel position at this fraction or more, on both strands and with ordinary base quality, is treated as a real allele in the flanking-context check rather than an error (default 0.01)")
    parser.add_argument("--fusion-region", default=None, help="Fusion-partner region, e.g. chr22:23522552-23660224 (BCR, hg19). A read that also aligns here (SA tag) is a fusion read. Without it every read is counted as non-fusion.")
    parser.add_argument("--fusion-reads-only", action="store_true", help="Count only fusion reads (requires --fusion-region)")
    parser.add_argument("--include-artifacts", action="store_true", help="Include ARTIFACT calls in mutations TSV")
    parser.add_argument("--threads", type=int, default=1, help="Number of worker threads for BAM decompression")
    return parser.parse_args(argv)


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
                raise ValueError(f"Malformed BED row {bed_path}:{line_idx}: expected at least 7 columns")

            chrom, start, end, name, ref, alt, aa_change = fields[:7]
            try:
                start_pos = int(start)
                end_pos = int(end)
            except ValueError as exc:
                raise ValueError(f"Invalid BED coordinate at {bed_path}:{line_idx}") from exc
            if start_pos < 0 or end_pos != start_pos + 1:
                raise ValueError(f"Expected one-base BED interval at {bed_path}:{line_idx}")
            if ref.upper() not in BASES or alt.upper() not in BASES or ref.upper() == alt.upper():
                raise ValueError(f"Invalid SNV REF/ALT at {bed_path}:{line_idx}")
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
                "cdna_pos": int(CDNA_RE.search(name).group(1)) if CDNA_RE.search(name) else None,
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


# A flank mask entry says which bases a read may show at one flanking offset
# without its context counting as dirty: True forgives any base (a
# --known-variants position), a set forgives exactly those bases, and an empty
# set or False forgives nothing.
MaskEntry = Any
NO_BASE: frozenset = frozenset()


def forgives(entry: MaskEntry, base: str) -> bool:
    """Whether a flank mask entry forgives ``base`` at its offset."""

    if entry is True:
        return True
    if not entry:
        return False
    return base in entry


def flank_masks(
    chrom: str, pos0: int, window: int, allowed: Dict[str, Any]
) -> Tuple[Tuple[MaskEntry, ...], Tuple[MaskEntry, ...]]:
    """The mask entries for the ``window`` bases either side of ``pos0``.

    ``allowed`` maps a contig to {position: entry}. A plain set of positions
    (the form before D25) forgives any base at them.
    """

    sites = allowed.get(chrom, {})
    if isinstance(sites, (set, frozenset)):
        sites = {position: True for position in sites}
    left = tuple(sites.get(pos0 - window + i, NO_BASE) for i in range(window))
    right = tuple(sites.get(pos0 + 1 + i, NO_BASE) for i in range(window))
    return left, right


def context_is_clean(
    read_seq: str, ref_seq: str, mask: Tuple[MaskEntry, ...]
) -> bool:
    """Compare a flank to the reference, forgiving what the mask allows."""

    if len(read_seq) != len(ref_seq):
        return False
    for index, (read_base, ref_base) in enumerate(zip(read_seq, ref_seq)):
        if read_base == ref_base:
            continue
        if index < len(mask) and forgives(mask[index], read_base):
            continue
        return False
    return True


def flank_offsets(window: int) -> List[int]:
    """Genomic offsets of the flank, left then right: -window..-1, 1..window."""

    return [*range(-window, 0), *range(1, window + 1)]


def masks_from_offsets(allowed: Dict[int, MaskEntry], window: int) -> Tuple[Tuple[MaskEntry, ...], Tuple[MaskEntry, ...]]:
    left = tuple(allowed.get(offset, NO_BASE) for offset in range(-window, 0))
    right = tuple(allowed.get(offset, NO_BASE) for offset in range(1, window + 1))
    return left, right


def offsets_from_masks(left: Tuple[MaskEntry, ...], right: Tuple[MaskEntry, ...]) -> Dict[int, MaskEntry]:
    allowed = {index - len(left): entry for index, entry in enumerate(left)}
    allowed.update({index + 1: entry for index, entry in enumerate(right)})
    return allowed


def variant_classification(aa_change: str, profile: Dict[str, str]) -> Tuple[str, str]:
    """Separate sequence consequence from curated resistance significance."""
    aa = short_aa_notation(aa_change)
    match = re.fullmatch(r"([A-Z*])(\d+)([A-Z*])", aa)
    consequence = "SYNONYMOUS" if match and match.group(1) == match.group(3) else (
        "MISSENSE" if match else "OTHER"
    )
    if consequence == "SYNONYMOUS":
        return consequence, "SYNONYMOUS"
    if any(value in ("Resistant", "Intermediate") for value in profile.values()):
        return consequence, "CURATED_RESISTANCE"
    return consequence, "UNCLASSIFIED"


def codon_groups(sites: List[Dict[str, Any]]) -> Dict[Tuple[str, int], List[Dict[str, Any]]]:
    """Group panel substitutions by transcript codon (cDNA coordinates)."""
    grouped: Dict[Tuple[str, int], List[Dict[str, Any]]] = {}
    for site in sites:
        cpos = site.get("cdna_pos")
        if cpos is not None:
            grouped.setdefault((site["chrom"], (cpos - 1) // 3), []).append(site)
    return grouped


def infer_codon_genomic_positions(group: List[Dict[str, Any]]) -> Optional[Tuple[str, List[int]]]:
    """Infer a codon's genomic positions from its BED anchors (0-based).

    Assumes a plus-strand transcript, as ABL1 is. Anchors must agree on one
    cDNA-to-genome offset. One anchor is enough: no Branford position lies
    within 3 bp of an exon boundary, so every panel codon is contiguous, and if
    one were not, spliced reads would carry the skip inside it,
    collect_codon_haplotypes would drop them all, and the codon would come out
    unresolved rather than wrong.
    """
    c_to_g = {site["cdna_pos"]: site["pos"] - 1 for site in group if site.get("cdna_pos") is not None}
    if not c_to_g:
        return None
    offsets = {g - c for c, g in c_to_g.items()}
    if len(offsets) != 1:
        return None
    codon_start = (group[0]["cdna_pos"] - 1) // 3 * 3 + 1
    return group[0]["chrom"], [c_to_g.get(c, c + next(iter(offsets))) for c in range(codon_start, codon_start + 3)]


def collect_codon_haplotypes(
    bam: pysam.AlignmentFile,
    chrom: str,
    positions: List[int],
    min_baseq: int,
    min_mapq: int,
    read_filter: Optional[ReadQualityFilter] = None,
    reference_context: Optional[Tuple[str, str, Tuple[MaskEntry, ...], Tuple[MaskEntry, ...]]] = None,
    fusion_region: Optional[Tuple[frozenset, int, int]] = None,
    fusion_only: bool = False,
) -> Dict[Tuple[str, bool], Dict[str, Any]]:
    """Count full codons once per read, retaining supplementary-only ABL signal."""
    target = sorted(set(positions))
    by_read = {}
    if not target:
        return {}
    for aln in bam.fetch(chrom, target[0], target[-1] + 1):
        if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail or aln.is_secondary:
            continue
        if (aln.mapping_quality or 0) < min_mapq:
            continue
        if read_filter is not None and not read_filter.passes(aln):
            continue
        if fusion_only and not is_fusion_read(aln, fusion_region):
            continue
        # A codon is only interpretable when all three bases come from one
        # uninterrupted aligned segment: no insertion, deletion or reference
        # skip may fall between its bases.
        if _indel_or_skip_spans_interval(aln, target[0], target[-1]):
            continue
        qname_key = (aln.query_name, aln.is_read1, aln.is_read2)
        rank = (0 if not aln.is_supplementary else 1, -(aln.mapping_quality or 0), aln.reference_start)
        bases: Dict[int, Tuple[str, int, bool]] = {}
        query_positions = {}
        qualities = aln.query_qualities
        sequence = (aln.query_sequence or "").upper()
        for rpos, kind, qpos, _, _ in cigar_events(aln, target, 0):
            if kind != "base":
                continue
            bq = qualities[qpos] if qualities is not None else 30
            base = sequence[qpos]
            if base in BASES and bq >= min_baseq:
                bases[rpos] = (base, bq, aln.is_reverse)
                query_positions[rpos] = qpos
        # Rank only alignments that can provide a complete usable codon.
        if not set(target).issubset(bases):
            continue
        clean_context = None
        if reference_context is not None:
            left, right, left_mask, right_mask = reference_context
            q_start, q_end = min(query_positions.values()), max(query_positions.values()) + 1
            if q_start >= len(left) and q_end + len(right) <= len(sequence):
                clean_context = context_is_clean(sequence[q_start - len(left):q_start], left, left_mask) and context_is_clean(
                    sequence[q_end:q_end + len(right)], right, right_mask
                )
        old = by_read.get(qname_key)
        if old is None or rank < old[0]:
            by_read[qname_key] = (rank, bases, clean_context)
    out: Dict[Tuple[str, bool], Dict[str, Any]] = {}
    for _, bases, clean_context in by_read.values():
        # Codon positions are passed in transcript order and can be genomic +1.
        if not all(pos in bases for pos in positions):
            continue
        codon = "".join(bases[pos][0] for pos in positions)
        reverse = bases[positions[0]][2]
        item = out.setdefault((codon, reverse), {"count": 0, "bq_sum": 0, "kmer_total": 0, "kmer_clean": 0, "excluded": 0})
        # The clean rate is an audit of all reads, as for a single site...
        if clean_context is not None:
            item["kmer_total"] += 1
            item["kmer_clean"] += int(clean_context)
        # ...but a dirty flank makes the read evidence for no haplotype.
        if clean_context is False:
            item["excluded"] += 1
            continue
        item["count"] += 1
        item["bq_sum"] += sum(bases[pos][1] for pos in positions) / 3
    return out


def allele_quality(args, ref_fwd, ref_rev, alt_fwd, alt_rev, bq_sum, bq_count, k_clean=0, k_total=0):
    """The same non-SNR quality requirements apply to sites and full haplotypes."""
    count = alt_fwd + alt_rev
    minor = round(min(alt_fwd, alt_rev) / max(count, 1), 3)
    r1 = (alt_fwd + 1) * (ref_rev + 1)
    r2 = (alt_rev + 1) * (ref_fwd + 1)
    sor = round(math.log(max(r1 / r2, r2 / r1)), 2)
    baseq = round(bq_sum / bq_count, 1) if bq_count else 0.0
    clean_percent = round(100 * k_clean / k_total, 1) if k_total else None
    failures = []
    if ref_fwd > 0 and ref_rev > 0 and count >= 20 and minor < args.min_minor_strand_ratio and sor > args.max_sor:
        failures.append(f"STRAND_BIAS(ratio={minor:.3f}<{args.min_minor_strand_ratio:.2f},SOR={sor:.2f}>{args.max_sor:.1f})")
    if baseq < args.min_mean_baseq:
        failures.append(f"LOW_BASEQ({baseq:.1f}<{args.min_mean_baseq:.1f})")
    return {
        "alt_minor_strand_ratio": minor, "sor": sor, "mean_baseq_alt": baseq,
        "kmer_clean_percent": clean_percent if clean_percent is not None else "NA",
    }, failures


def unresolved_codon_records(all_results, group_pos, reason):
    """Retain marginal evidence when complete codon support cannot resolve it."""
    records = []
    for raw in all_results:
        if (raw["chrom"], raw["pos"]) not in group_pos:
            continue
        raw["codon_context_status"] = "CODON_CONTEXT_UNRESOLVED"
        if raw["status"] in ("DETECTED", "SUBCLONAL"):
            records.append({
                **raw, "record_type": "CODON_CONTEXT_UNRESOLVED",
                "mutation": "CODON_CONTEXT_UNRESOLVED", "aa_change": "",
                "consequence": "OTHER", "resistance_class": "UNCLASSIFIED",
                "filter_reason": reason,
            })
    return records


def _indel_or_skip_spans_interval(aln: pysam.AlignedSegment, first: int, last: int) -> bool:
    """Return whether an insertion, deletion, or reference skip crosses codon bases."""
    ref_pos = aln.reference_start
    for op, length in aln.cigartuples or ():
        if op == 1:  # insertion between reference bases
            if first < ref_pos <= last:
                return True
        elif op in (2, 3):  # deletion or reference skip
            if ref_pos <= last and ref_pos + length > first:
                return True
            ref_pos += length
        elif op in (0, 7, 8):  # aligned bases
            ref_pos += length
    return False


FUSION, NONFUSION = "fusion", "nonfusion"


def parse_region(text: Optional[str]) -> Optional[Tuple[frozenset, int, int]]:
    """'chr22:23522552-23660224' -> ({'chr22', '22'}, 23522552, 23660224); 1-based, inclusive."""

    if not text:
        return None
    match = re.fullmatch(r"\s*([^:\s]+):(\d[\d,]*)-(\d[\d,]*)\s*", str(text))
    if not match:
        raise ValueError(f"Invalid region {text!r}: expected chrom:start-end")
    chrom = match.group(1)
    start, end = (int(match.group(index).replace(",", "")) for index in (2, 3))
    if end < start:
        raise ValueError(f"Invalid region {text!r}: end before start")
    alias = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
    return frozenset({chrom, alias}), start, end


def is_fusion_read(aln: pysam.AlignedSegment, region: Optional[Tuple[frozenset, int, int]]) -> bool:
    """Whether the read also aligns in the fusion-partner region, by its SA tag."""

    if region is None or not aln.has_tag("SA"):
        return False
    names, start, end = region
    for entry in str(aln.get_tag("SA")).split(";"):
        fields = entry.split(",")
        if len(fields) < 2 or fields[0] not in names:
            continue
        try:
            position = int(fields[1])
        except ValueError:
            continue
        if start <= position <= end:
            return True
    return False


class SiteEvidence:
    """What one read set says about one panel position, before any flank policy.

    Reads are kept by base, strand and flank signature -- the flank bases that
    differ from the reference, as (genomic offset, base) -- so any policy for
    which differences to forgive can be applied afterwards without reading the
    BAM again. A signature of None means the flank could not be checked (no
    reference, or the read ends inside the window); () means it is clean.
    """

    __slots__ = ("alleles", "refskip", "deletion", "insertion", "low_baseq", "unfiltered",
                 "flank_reads", "flank_alt", "flank_alt_bq")

    def __init__(self) -> None:
        self.alleles: Dict[Tuple[str, str, Optional[Tuple]], List[int]] = {}
        self.refskip = 0
        self.deletion = 0
        self.insertion = 0
        self.low_baseq = 0
        self.unfiltered = {base: 0 for base in BASES}
        # Substitutions beside the site, counted only in reads whose window is
        # one gapless aligned block (read offset == genomic offset) and only at
        # flank bases that pass --min-base-quality: what the sample carries
        # there, as opposed to what a homopolymer indel makes it look like.
        self.flank_reads = {"fwd": 0, "rev": 0}
        self.flank_alt: Dict[Tuple[int, str, str], int] = {}
        self.flank_alt_bq: Dict[Tuple[int, str], int] = {}


def empty_counts() -> Dict[str, Any]:
    counts: Dict[str, Any] = {f"{base}_{strand}": 0 for base in BASES for strand in ("fwd", "rev")}
    counts.update({
        "refskip": 0,
        "deletion": 0,
        "insertion": 0,
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
        # Reads whose flanking bases disagree with the reference: not
        # counted for any allele (see FLANKING CONTEXT in the module docs).
        "context_excluded": {b: 0 for b in BASES},
    })
    return counts


def cigar_events(aln: pysam.AlignedSegment, targets: List[int], window: int) -> List[Tuple]:
    """Where each sorted 0-based target falls in one alignment.

    Returns (target, kind, query_pos, insertion_follows, gapless_window) with
    kind "base", "del" or "skip" -- what a pileup column would say about this
    read, without walking every column in between. gapless_window is whether
    the ``window`` bases either side sit in the same gapless aligned block.
    """

    events: List[Tuple] = []
    ref = aln.reference_start
    if ref is None or ref < 0:
        return events
    count = len(targets)
    index = bisect.bisect_left(targets, ref)
    if index >= count:
        return events
    ops = aln.cigartuples or ()
    cigar = aln.cigarstring or ""
    if "=" in cigar or "X" in cigar:
        merged: List[Tuple[int, int]] = []
        for op, length in ops:
            op = 0 if op in (7, 8) else op
            if merged and op == 0 and merged[-1][0] == 0:
                merged[-1] = (0, merged[-1][1] + length)
            else:
                merged.append((op, length))
        ops = merged
    query = 0
    last_op = len(ops) - 1
    for op_index, (op, length) in enumerate(ops):
        if op == 0:
            end = ref + length
            while index < count and targets[index] < end:
                target = targets[index]
                insertion = target == end - 1 and op_index < last_op and ops[op_index + 1][0] == 1
                gapless = target - window >= ref and target + window < end
                events.append((target, "base", query + target - ref, insertion, gapless))
                index += 1
            ref = end
            query += length
        elif op in (2, 3):
            end = ref + length
            while index < count and targets[index] < end:
                events.append((targets[index], "del" if op == 2 else "skip", None, False, False))
                index += 1
            ref = end
        elif op in (1, 4):
            query += length
        if index >= count:
            break
    return events


def scan_positions(
    bam: pysam.AlignmentFile,
    chrom: str,
    positions_0based: Set[int],
    min_baseq: int,
    min_mapq: int,
    contexts: Optional[Dict[int, Tuple[str, str]]] = None,
    read_filter: Optional[ReadQualityFilter] = None,
    fusion_region: Optional[Tuple[frozenset, int, int]] = None,
) -> Dict[int, Dict[str, SiteEvidence]]:
    """Read every alignment over the panel once and record what it says at each position.

    This replaces a pileup over the whole span. A cDNA read spans the ABL1
    introns as reference skips, so a pileup visits ~17 000 columns -- nearly all
    intronic -- at full depth to reach 52; one fetch with a CIGAR walk per read
    gives the same counts in about a tenth of the time (D25).

    Each read is counted once per position. When a read has several usable
    alignments over a position, the primary wins over a supplementary, then
    the higher mapping quality, then the first in BAM order -- as the pileup
    implementation did.
    """

    targets = sorted(positions_0based)
    evidence = {position: {FUSION: SiteEvidence(), NONFUSION: SiteEvidence()} for position in targets}
    if not targets:
        return evidence
    contexts = contexts or {}
    window = max((len(left) for left, _ in contexts.values()), default=0)
    first, last = targets[0], targets[-1] + 1

    def usable(aln) -> bool:
        if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail or aln.is_secondary:
            return False
        if (aln.mapping_quality or 0) < min_mapq:
            return False
        return read_filter is None or read_filter.passes(aln)

    # Which reads have more than one usable alignment here; only those need
    # ranking. (This pass also settles the read-quality filter for every
    # alignment of a read before any is counted.)
    seen: Dict[Tuple, int] = {}
    for aln in bam.fetch(chrom, first, last):
        if usable(aln):
            key = (aln.query_name, aln.is_read1, aln.is_read2)
            seen[key] = seen.get(key, 0) + 1
    several = {key for key, number in seen.items() if number > 1}

    def tally(aln, sequence, qualities, fusion, event) -> None:
        target, kind, qpos, insertion, gapless = event
        site = evidence[target][FUSION if fusion else NONFUSION]
        if kind == "skip":
            site.refskip += 1
            return
        if kind == "del":
            site.deletion += 1
            return
        base = sequence[qpos]
        if base not in BASES:
            return
        site.unfiltered[base] += 1
        if insertion:
            site.insertion += 1
        strand = "rev" if aln.is_reverse else "fwd"
        signature = None
        context = contexts.get(target)
        if context is not None:
            ref_left, ref_right = context
            k = len(ref_left)
            if qpos >= k and qpos + 1 + k <= len(sequence):
                read_left = sequence[qpos - k:qpos]
                read_right = sequence[qpos + 1:qpos + 1 + k]
                if read_left == ref_left and read_right == ref_right:
                    signature = ()
                else:
                    signature = tuple(
                        (i - k, b) for i, (b, r) in enumerate(zip(read_left, ref_left)) if b != r
                    ) + tuple(
                        (i + 1, b) for i, (b, r) in enumerate(zip(read_right, ref_right)) if b != r
                    )
                if gapless:
                    site.flank_reads[strand] += 1
                    for offset, flank_base in signature:
                        flank_bq = qualities[qpos + offset] if qualities is not None else 30
                        if flank_base in BASES and flank_bq >= min_baseq:
                            key = (offset, flank_base, strand)
                            site.flank_alt[key] = site.flank_alt.get(key, 0) + 1
                            site.flank_alt_bq[(offset, flank_base)] = site.flank_alt_bq.get((offset, flank_base), 0) + flank_bq
        bq = qualities[qpos] if qualities is not None else 30
        if min_baseq > 0 and bq < min_baseq:
            site.low_baseq += 1
            return
        slot = site.alleles.get((base, strand, signature))
        if slot is None:
            site.alleles[(base, strand, signature)] = [1, bq]
        else:
            slot[0] += 1
            slot[1] += bq

    best: Dict[Tuple, Tuple] = {}
    for aln in bam.fetch(chrom, first, last):
        if not usable(aln):
            continue
        events = cigar_events(aln, targets, window)
        if not events:
            continue
        key = (aln.query_name, aln.is_read1, aln.is_read2)
        if key in several:
            rank = (0 if not aln.is_supplementary else 1, -(aln.mapping_quality or 0))
            for event in events:
                previous = best.get((event[0], key))
                if previous is None or rank < previous[0]:
                    best[(event[0], key)] = (rank, aln, event)
            continue
        sequence = (aln.query_sequence or "").upper()
        qualities = aln.query_qualities
        fusion = is_fusion_read(aln, fusion_region)
        for event in events:
            tally(aln, sequence, qualities, fusion, event)
    for _, aln, event in best.values():
        tally(aln, (aln.query_sequence or "").upper(), aln.query_qualities,
              is_fusion_read(aln, fusion_region), event)
    return evidence


def materialise_counts(evidence: List[SiteEvidence], allowed: Dict[int, MaskEntry]) -> Dict[str, Any]:
    """Counts for the given read sets with flank differences forgiven as ``allowed`` says."""

    counts = empty_counts()
    for site in evidence:
        counts["refskip"] += site.refskip
        counts["deletion"] += site.deletion
        counts["insertion"] += site.insertion
        counts["low_baseq"] += site.low_baseq
        for base in BASES:
            counts["unfiltered"][base] += site.unfiltered[base]
        for (base, strand, signature), (number, bq_sum) in site.alleles.items():
            if signature is not None:
                counts["kmer_total_count"][base] += number
                if not all(forgives(allowed.get(offset, NO_BASE), flank_base) for offset, flank_base in signature):
                    counts["context_excluded"][base] += number
                    continue
                counts["kmer_clean_count"][base] += number
            counts[f"{base}_{strand}"] += number
            counts["base_qual_sum"][base] += bq_sum
            counts["base_qual_count"][base] += number
    return counts


def count_positions_batch(
    bam: pysam.AlignmentFile,
    chrom: str,
    positions_0based: Set[int],
    max_depth: int,
    min_baseq: int,
    min_mapq: int,
    ref_seqs_by_pos: Optional[Dict[Tuple[str, int], Tuple[str, str, Tuple, Tuple]]] = None,
    orig_chrom: Optional[str] = None,
    read_filter: Optional[ReadQualityFilter] = None,
) -> Dict[int, Dict[str, Any]]:
    """Base counts at each position with fixed flank masks (all reads).

    ``max_depth`` is accepted for compatibility; the read scan has no cap.
    """

    contexts: Dict[int, Tuple[str, str]] = {}
    allowed: Dict[int, Dict[int, MaskEntry]] = {}
    for position in positions_0based:
        key = (orig_chrom or chrom, position)
        if ref_seqs_by_pos and key in ref_seqs_by_pos:
            left, right, left_mask, right_mask = ref_seqs_by_pos[key]
            contexts[position] = (left, right)
            allowed[position] = offsets_from_masks(left_mask, right_mask)
    evidence = scan_positions(bam, chrom, positions_0based, min_baseq, min_mapq, contexts, read_filter)
    return {
        position: materialise_counts(list(groups.values()), allowed.get(position, {}))
        for position, groups in evidence.items()
    }


def count_position(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos_1based: int,
    max_depth: int,
    min_baseq: int,
    min_mapq: int,
    read_filter: Optional[ReadQualityFilter] = None,
) -> Dict[str, Any]:
    """Extract forward/reverse base counts at a single coordinate (wrapper for count_positions_batch)."""
    target_0based = pos_1based - 1
    batch_res = count_positions_batch(
        bam, chrom, {target_0based}, max_depth, min_baseq, min_mapq, read_filter=read_filter
    )
    return batch_res.get(target_0based, {f"{b}_{s}": 0 for b in BASES for s in ("fwd", "rev")})


def screen_flank_variants(
    args,
    sites: List[Dict[str, Any]],
    evidence: Dict[Tuple[str, int], Dict[str, SiteEvidence]],
    read_sets: Tuple[str, ...],
    reference_flanks: Dict[Tuple[str, int], Tuple[str, str]],
    panel_alleles: Dict[Tuple[str, int], set],
) -> Dict[Tuple[str, int], Dict[str, Dict[str, Any]]]:
    """Substitutions the sample carries beside panel positions (D25, fix 1).

    A base at least --flank-variant-min-af of the gapless reads, with
    --min-alt-reads support and passing the strand and base-quality checks a
    called allele must pass, is a real allele -- a germline SNV on the fusion
    allele, or a change in cis with a panel mutation -- not a context error.
    Panel ALTs are left out: whether they are forgiven depends on whether they
    are called (fix 2).
    """

    found: Dict[Tuple[str, int], Dict[str, Dict[str, Any]]] = {}
    for key in sorted({(site["chrom"], site["pos"] - 1) for site in sites}):
        groups = evidence.get(key)
        flanks = reference_flanks.get(key)
        if not groups or not flanks:
            continue
        chrom, position = key
        ref_left, ref_right = flanks
        k = len(ref_left)
        reads = {strand: sum(groups[name].flank_reads[strand] for name in read_sets) for strand in ("fwd", "rev")}
        depth = reads["fwd"] + reads["rev"]
        if depth == 0:
            continue
        alt: Dict[Tuple[int, str, str], int] = {}
        alt_bq: Dict[Tuple[int, str], int] = {}
        for name in read_sets:
            for item, number in groups[name].flank_alt.items():
                alt[item] = alt.get(item, 0) + number
            for item, total in groups[name].flank_alt_bq.items():
                alt_bq[item] = alt_bq.get(item, 0) + total
        for offset in flank_offsets(k):
            ref_base = ref_left[offset + k] if offset < 0 else ref_right[offset - 1]
            if ref_base not in BASES:
                continue
            genomic = position + offset
            mismatched = {
                strand: sum(alt.get((offset, base, strand), 0) for base in BASES) for strand in ("fwd", "rev")
            }
            for base in BASES:
                if base == ref_base or base in panel_alleles.get((chrom, genomic), ()):
                    continue
                alt_fwd, alt_rev = alt.get((offset, base, "fwd"), 0), alt.get((offset, base, "rev"), 0)
                number = alt_fwd + alt_rev
                fraction = number / depth
                if number < args.min_alt_reads or fraction < args.flank_variant_min_af:
                    continue
                _, failures = allele_quality(
                    args, reads["fwd"] - mismatched["fwd"], reads["rev"] - mismatched["rev"],
                    alt_fwd, alt_rev, alt_bq.get((offset, base), 0), number,
                )
                if failures:
                    continue
                previous = found.get((chrom, genomic), {}).get(base)
                if previous is None or fraction > previous["af"]:
                    found.setdefault((chrom, genomic), {})[base] = {
                        "ref": ref_base, "af": fraction, "alt_reads": number, "depth": depth,
                    }
    return found


def flank_allowed(
    chrom: str,
    position: int,
    window: int,
    known: Dict[str, set],
    panel_alleles: Dict[Tuple[str, int], set],
    called: Set[Tuple[str, int, str]],
    sample_variants: Dict[Tuple[str, int], Dict[str, Dict[str, Any]]],
) -> Dict[int, MaskEntry]:
    """Which flank bases are real alleles rather than context errors (D25)."""

    allowed: Dict[int, MaskEntry] = {}
    known_here = known.get(chrom, set())
    for offset in flank_offsets(window):
        genomic = position + offset
        if genomic in known_here:
            allowed[offset] = True
            continue
        bases = {alt for alt in panel_alleles.get((chrom, genomic), ()) if (chrom, genomic, alt) in called}
        bases |= set(sample_variants.get((chrom, genomic), {}))
        allowed[offset] = frozenset(bases)
    return allowed


def describe_forgiven(
    chrom: str,
    position: int,
    allowed: Dict[int, MaskEntry],
    reference_flanks: Dict[Tuple[str, int], Tuple[str, str]],
    sample_variants: Dict[Tuple[str, int], Dict[str, Dict[str, Any]]],
) -> str:
    """The forgiven flank alleles as 1-based 'pos:REF>ALT(source)', for the audit table."""

    flanks = reference_flanks.get((chrom, position))
    if not flanks:
        return ""
    ref_left, ref_right = flanks
    k = len(ref_left)
    described = []
    for offset in sorted(allowed):
        entry = allowed[offset]
        if not entry:
            continue
        genomic = position + offset
        ref_base = ref_left[offset + k] if offset < 0 else ref_right[offset - 1]
        if entry is True:
            described.append(f"{genomic + 1}:{ref_base}>N(known)")
            continue
        for base in sorted(entry):
            variant = sample_variants.get((chrom, genomic), {}).get(base)
            source = f"sample,{variant['af'] * 100:.2f}%" if variant else "panel"
            described.append(f"{genomic + 1}:{ref_base}>{base}({source})")
    return ",".join(described)


def evaluate_site(
    args,
    site: Dict[str, Any],
    counts: Dict[str, Any],
    panel_alts_by_pos: Dict[Tuple[str, int], set],
    tki_profiles: Dict[str, Dict[str, str]],
    panel_sites_expected: int,
) -> Dict[str, Any]:
    """Status, filters and annotation of one panel allele from its counts."""

    pos = site["pos"]
    ref = site["ref"]
    alt = site["alt"]
    short_mut = site["short_mutation"]

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

    quality_by_base = {
        b: allele_quality(
            args, ad_ref_fwd, ad_ref_rev,
            counts.get(f"{b}_fwd", 0), counts.get(f"{b}_rev", 0),
            counts.get("base_qual_sum", {}).get(b, 0),
            counts.get("base_qual_count", {}).get(b, 0),
            counts.get("kmer_clean_count", {}).get(b, 0),
            counts.get("kmer_total_count", {}).get(b, 0),
        ) for b in BASES
    }
    quality, quality_failures = quality_by_base[alt]
    bq_sum_ref = counts.get("base_qual_sum", {}).get(ref, 0)
    bq_cnt_ref = counts.get("base_qual_count", {}).get(ref, 0)
    mean_baseq_ref = round(bq_sum_ref / bq_cnt_ref, 1) if bq_cnt_ref else 0.0

    # A panel ALT is a credible independent allele only if it has the
    # ordinary minimum support/frequency, passes quality checks, and dominates the weakest
    # observed non-reference base by the configured SNR. Without this
    # gate, a balanced three-base error spectrum makes every panel ALT
    # exempt itself from the local noise estimate.
    candidate_panel_alts = panel_alts_by_pos.get((site["chrom"], pos), set())
    nonref_floor = min((per_base[b] for b in BASES if b != ref), default=0)
    credible_panel_alts = {
        b for b in candidate_panel_alts
        if b != ref and per_base[b] >= args.min_alt_reads
        and depth > 0 and (per_base[b] / depth) >= args.subclonal_af
        and per_base[b] >= args.min_snr * max(nonref_floor, 1)
        and not quality_by_base[b][1]
    }
    other_bases = {
        b: per_base[b] for b in BASES
        if b != ref and b != alt and b not in credible_panel_alts
    }
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
        fails.extend(quality_failures)

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

    tki_info = tki_profiles.get(short_mut, {})
    consequence, resistance_class = variant_classification(short_mut, tki_info)

    return {
        "record_type": "SITE_AUDIT",
        "panel_sites_expected": panel_sites_expected,
        "sample": args.sample,
        "chrom": site["chrom"],
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "mutation": short_mut,
        "aa_change": site["aa_change"],
        "consequence": consequence,
        "resistance_class": resistance_class,
        "codon_context_status": "SINGLE_SITE_ANNOTATION",
        "name": site["name"],
        "depth": depth,
        "ad_ref": ad_ref,
        "ad_alt": ad_alt,
        "vaf_percent": round(vaf * 100.0, 3),
        "count_A": per_base.get("A", 0),
        "count_C": per_base.get("C", 0),
        "count_G": per_base.get("G", 0),
        "count_T": per_base.get("T", 0),
        "pct_A": round((per_base.get("A", 0) / depth * 100.0), 3) if depth > 0 else 0.0,
        "pct_C": round((per_base.get("C", 0) / depth * 100.0), 3) if depth > 0 else 0.0,
        "pct_G": round((per_base.get("G", 0) / depth * 100.0), 3) if depth > 0 else 0.0,
        "pct_T": round((per_base.get("T", 0) / depth * 100.0), 3) if depth > 0 else 0.0,
        "depth_unfiltered": depth_unfiltered,
        "ad_alt_unfiltered": ad_alt_unfiltered,
        "vaf_unfiltered_percent": round(vaf_unfiltered * 100.0, 3),
        "status": status,
        "filter_reason": filter_reason,
        "snr": snr,
        "ad_other": ad_other,
        "vaf_other_percent": vaf_other_percent,
        "other_bases": other_bases_str,
        "credible_panel_alt_competitors": ",".join(sorted(credible_panel_alts - {alt})),
        "alt_minor_strand_ratio": quality["alt_minor_strand_ratio"],
        "sor": quality["sor"],
        "mean_baseq_ref": mean_baseq_ref,
        "mean_baseq_alt": quality["mean_baseq_alt"],
        "kmer_clean_percent": quality["kmer_clean_percent"],
        "imatinib": tki_info.get("Imatinib", "Unknown"),
        "dasatinib": tki_info.get("Dasatinib", "Unknown"),
        "nilotinib": tki_info.get("Nilotinib", "Unknown"),
        "bosutinib": tki_info.get("Bosutinib", "Unknown"),
        "ponatinib": tki_info.get("Ponatinib", "Unknown"),
        "asciminib": tki_info.get("Asciminib", "Unknown"),
        "ad_ref_fwd": ad_ref_fwd,
        "ad_ref_rev": ad_ref_rev,
        "ad_alt_fwd": ad_alt_fwd,
        "ad_alt_rev": ad_alt_rev,
        "n_refskip": counts.get("refskip", 0),
        "n_deletion": counts.get("deletion", 0),
        "n_insertion": counts.get("insertion", 0),
        "depth_context_excluded": sum(counts.get("context_excluded", {}).values()),
        "ad_alt_context_excluded": counts.get("context_excluded", {}).get(alt, 0),
    }


def main(argv=None) -> None:
    args = parse_args(argv)
    if not args.bam.is_file():
        sys.exit(f"BAM file not found: {args.bam}")
    fusion_region = parse_region(args.fusion_region)
    if args.fusion_reads_only and fusion_region is None:
        sys.exit("--fusion-reads-only needs --fusion-region")

    sites = read_branford_bed(args.bed)
    if not sites:
        sys.exit(f"No sites parsed from {args.bed}")

    # An explicit annotation resource is part of the requested analysis and
    # must not silently fall back to a different clinical profile set.
    features = load_features(args.features_json)
    tki_profiles = load_tki_profiles(args.features_json)

    kw = args.kmer_window
    known_positions: Dict[str, set] = load_exempt_positions(getattr(args, "known_variants", None))
    panel_alleles: Dict[Tuple[str, int], set] = {}
    for site in sites:
        panel_alleles.setdefault((site["chrom"], site["pos"] - 1), set()).add(site["alt"])
    every_panel_allele = {(chrom, position, alt) for (chrom, position), alts in panel_alleles.items() for alt in alts}

    # The reference bases beside each panel position, for the flanking-context
    # check. Keyed by the BED contig name and 0-based position.
    reference_flanks: Dict[Tuple[str, int], Tuple[str, str]] = {}
    fasta_refs: set = set()
    if args.reference_fasta:
        try:
            ref_fasta = pysam.FastaFile(str(args.reference_fasta))
            fasta_refs = set(ref_fasta.references)
            for site in sites:
                c = site["chrom"]
                if c not in fasta_refs:
                    alt_c = c[3:] if c.startswith("chr") else f"chr{c}"
                    if alt_c in fasta_refs:
                        c = alt_c
                if c not in fasta_refs:
                    raise ValueError(f"BED contig {site['chrom']} is absent from FASTA")
                p0 = site["pos"] - 1
                fasta_ref = ref_fasta.fetch(c, p0, p0 + 1).upper()
                if fasta_ref != site["ref"]:
                    raise ValueError(
                        f"BED REF mismatch at {site['chrom']}:{site['pos']}: "
                        f"BED={site['ref']} FASTA={fasta_ref}"
                    )
                c_len = ref_fasta.get_reference_length(c)
                if p0 >= kw and (p0 + 1 + kw) <= c_len:
                    reference_flanks[(site["chrom"], p0)] = (
                        ref_fasta.fetch(c, p0 - kw, p0).upper(),
                        ref_fasta.fetch(c, p0 + 1, p0 + 1 + kw).upper(),
                    )
            ref_fasta.close()
        except Exception as e:
            if isinstance(e, ValueError):
                raise
            raise RuntimeError(f"Could not validate reference FASTA {args.reference_fasta}: {e}") from e

    bam = pysam.AlignmentFile(str(args.bam), "rb", threads=getattr(args, "threads", 1))
    bam_contigs = set(bam.references)

    def bam_contig(chrom: str) -> str:
        if chrom not in bam_contigs:
            alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            if alt_chrom in bam_contigs:
                return alt_chrom
        return chrom

    # Group all target sites by resolved BAM contig for single-pass scanning
    chrom_positions: Dict[str, Set[int]] = {}
    chrom_alias_map: Dict[str, str] = {}
    for site in sites:
        chrom = bam_contig(site["chrom"])
        chrom_alias_map[chrom] = site["chrom"]
        chrom_positions.setdefault(chrom, set()).add(site["pos"] - 1)

    read_filter = ReadQualityFilter(args.min_read_qs)
    read_sets = (FUSION,) if args.fusion_reads_only else (FUSION, NONFUSION)
    evidence: Dict[Tuple[str, int], Dict[str, SiteEvidence]] = {}
    try:
        for chrom, pos_set in chrom_positions.items():
            if chrom not in bam_contigs:
                continue
            orig_chrom = chrom_alias_map.get(chrom, chrom)
            contexts = {
                position: reference_flanks[(orig_chrom, position)]
                for position in pos_set if (orig_chrom, position) in reference_flanks
            }
            scanned = scan_positions(
                bam, chrom, pos_set, args.min_base_quality, args.min_mapping_quality,
                contexts=contexts, read_filter=read_filter, fusion_region=fusion_region,
            )
            for position, groups in scanned.items():
                evidence[(orig_chrom, position)] = groups

        # Fix 1: real substitutions beside panel positions, found from the reads.
        sample_variants = screen_flank_variants(
            args, sites, evidence, read_sets, reference_flanks, panel_alleles
        )

        panel_sites_expected = len({(s["chrom"], s["pos"]) for s in sites})
        panel_alts_by_pos: Dict[Tuple[str, int], set] = {}
        for panel_site in sites:
            panel_alts_by_pos.setdefault((panel_site["chrom"], panel_site["pos"]), set()).add(panel_site["alt"])

        def site_counts(site, called, groups):
            key = (site["chrom"], site["pos"] - 1)
            allowed = flank_allowed(site["chrom"], site["pos"] - 1, kw, known_positions,
                                    panel_alleles, called, sample_variants)
            site_evidence = evidence.get(key)
            if site_evidence is None:
                return empty_counts(), allowed
            return materialise_counts([site_evidence[name] for name in groups], allowed), allowed

        # Fix 2, round one: forgive every panel ALT at a neighbouring position
        # and see which ALTs are called...
        provisional_called = set()
        for site in sites:
            counts, _ = site_counts(site, every_panel_allele, read_sets)
            record = evaluate_site(args, site, counts, panel_alts_by_pos, tki_profiles, panel_sites_expected)
            if record["status"] in ("DETECTED", "SUBCLONAL"):
                provisional_called.add((site["chrom"], site["pos"] - 1, site["alt"]))

        # ...round two: forgive only those. A neighbour's ALT that is noise in
        # this sample keeps the reads it appears in out of the count.
        all_results: List[Dict[str, Any]] = []
        detected_mutations: List[Dict[str, Any]] = []
        for site in sites:
            counts, allowed = site_counts(site, provisional_called, read_sets)
            record = evaluate_site(args, site, counts, panel_alts_by_pos, tki_profiles, panel_sites_expected)
            record["flank_alleles_forgiven"] = describe_forgiven(
                site["chrom"], site["pos"] - 1, allowed, reference_flanks, sample_variants
            )
            for name in (FUSION, NONFUSION):
                if fusion_region is None:
                    record.update({f"depth_{name}": "NA", f"ad_alt_{name}": "NA", f"vaf_{name}_percent": "NA"})
                    continue
                part, _ = site_counts(site, provisional_called, (name,))
                part_depth = sum(part[f"{b}_{s}"] for b in BASES for s in ("fwd", "rev"))
                part_alt = part[f"{site['alt']}_fwd"] + part[f"{site['alt']}_rev"]
                record.update({
                    f"depth_{name}": part_depth,
                    f"ad_alt_{name}": part_alt,
                    f"vaf_{name}_percent": round(100.0 * part_alt / part_depth, 3) if part_depth else 0.0,
                })
            all_results.append(record)

            if record["status"] in ("DETECTED", "SUBCLONAL"):
                # Keep SITE_AUDIT in all-sites.tsv for nucleotide-level
                # traceability, but give downstream variant consumers a
                # distinct record type.
                detected_mutations.append({**record, "record_type": "SITE_VARIANT"})
            elif args.include_artifacts and record["status"] == "ARTIFACT":
                detected_mutations.append({**record, "record_type": "SITE_VARIANT"})

        # Two variant positions in one codon -- two called panel positions, or
        # a called panel position and a sample variant beside it -- make the
        # individual BED amino-acid annotations marginal statements that may
        # not describe the observed haplotype. Replace those mutation rows with
        # read-supported codon haplotypes when a trustworthy reference codon is
        # available. The complete nucleotide-level table above remains intact.
        called_alleles = {
            (r["chrom"], r["pos"] - 1, r["alt"]) for r in all_results
            if r["status"] in ("DETECTED", "SUBCLONAL")
        }
        detected_positions = {(chrom, position + 1) for chrom, position, _ in called_alleles}
        for (_, codon_number), group in codon_groups(sites).items():
            detected_group_positions = sorted({
                s["pos"] for s in group if (s["chrom"], s["pos"]) in detected_positions
            })
            if not detected_group_positions:
                continue
            genomic = infer_codon_genomic_positions(group)
            codon_chrom = group[0]["chrom"]
            variant_positions = {p - 1 for p in detected_group_positions}
            if genomic:
                variant_positions |= {p for p in genomic[1] if sample_variants.get((codon_chrom, p))}
            if len(variant_positions) < 2:
                continue
            group_pos = {(s["chrom"], s["pos"]) for s in group if s["pos"] in detected_group_positions}
            for raw in all_results:
                if (raw["chrom"], raw["pos"]) in group_pos:
                    raw["codon_context_status"] = "CODON_CONTEXT_REQUIRED"
                    raw["resistance_class"] = "UNCLASSIFIED"

            reference_codon = None
            fasta_contig = codon_chrom
            if genomic and args.reference_fasta:
                codon_chrom, codon_positions = genomic
                if fasta_contig not in fasta_refs:
                    fasta_contig = codon_chrom[3:] if codon_chrom.startswith("chr") else f"chr{codon_chrom}"
                if fasta_contig in fasta_refs:
                    with pysam.FastaFile(str(args.reference_fasta)) as fasta:
                        reference_codon = "".join(fasta.fetch(fasta_contig, p, p + 1).upper() for p in codon_positions)
                    if len(reference_codon) != 3 or any(b not in BASES for b in reference_codon):
                        reference_codon = None
                    for s in group:
                        c_index = (s["cdna_pos"] - 1) % 3
                        if reference_codon and reference_codon[c_index] != s["ref"]:
                            raise ValueError(
                                f"BED REF disagrees with reconstructed codon at {s['chrom']}:{s['pos']}"
                            )
            if not genomic or not reference_codon:
                detected_mutations.extend(unresolved_codon_records(
                    all_results, group_pos, "CODON_CONTEXT_UNRESOLVED"
                ))
                detected_mutations = [r for r in detected_mutations if not (r.get("record_type") in {"SITE_AUDIT", "SITE_VARIANT"} and (r["chrom"], r["pos"]) in group_pos)]
                continue

            reference_context = None
            with pysam.FastaFile(str(args.reference_fasta)) as fasta:
                first, last = min(codon_positions), max(codon_positions)
                if first >= kw and last + 1 + kw <= fasta.get_reference_length(fasta_contig):
                    left_mask, _ = masks_from_offsets(flank_allowed(
                        codon_chrom, first, kw, known_positions, panel_alleles, called_alleles, sample_variants), kw)
                    _, right_mask = masks_from_offsets(flank_allowed(
                        codon_chrom, last, kw, known_positions, panel_alleles, called_alleles, sample_variants), kw)
                    reference_context = (
                        fasta.fetch(fasta_contig, first - kw, first).upper(),
                        fasta.fetch(fasta_contig, last + 1, last + 1 + kw).upper(),
                        left_mask, right_mask,
                    )
            codon_start = min(codon_positions) + 1
            bam_chrom = bam_contig(codon_chrom)
            if bam_chrom not in bam_contigs:
                hap_counts = {}
            else:
                hap_counts = collect_codon_haplotypes(
                    bam, bam_chrom, codon_positions, args.min_base_quality, args.min_mapping_quality,
                    read_filter=read_filter, reference_context=reference_context,
                    fusion_region=fusion_region, fusion_only=args.fusion_reads_only,
                )
            state_counts: Dict[str, Dict[str, int]] = {}
            for (codon, reverse), values in hap_counts.items():
                state_counts.setdefault(codon, {"fwd": 0, "rev": 0, "bq_sum": 0, "kmer_total": 0, "kmer_clean": 0, "excluded": 0})
                strand = "rev" if reverse else "fwd"
                state_counts[codon][strand] += values["count"]
                state_counts[codon]["bq_sum"] += values["bq_sum"]
                state_counts[codon]["kmer_total"] += values["kmer_total"]
                state_counts[codon]["kmer_clean"] += values["kmer_clean"]
                state_counts[codon]["excluded"] += values.get("excluded", 0)
            codon_depth = sum(v["fwd"] + v["rev"] for v in state_counts.values())
            codon_excluded = sum(v["excluded"] for v in state_counts.values())
            unresolved_reason = None
            if codon_depth < args.min_depth:
                unresolved_reason = f"LOW_CODON_DEPTH({codon_depth}<{args.min_depth})"
            else:
                for site in group:
                    if not any(r["pos"] == site["pos"] and r["alt"] == site["alt"] and
                               r["status"] in ("DETECTED", "SUBCLONAL") for r in all_results
                               if (r["chrom"], r["pos"]) in group_pos):
                        continue
                    index = (site["cdna_pos"] - 1) % 3
                    support = sum(v["fwd"] + v["rev"] for c, v in state_counts.items() if c[index] == site["alt"])
                    if support < args.min_alt_reads:
                        unresolved_reason = "INSUFFICIENT_CODON_ALT_SUPPORT"
                        break
            if unresolved_reason:
                detected_mutations = [r for r in detected_mutations if not (
                    r.get("record_type") == "SITE_VARIANT" and (r["chrom"], r["pos"]) in group_pos
                )]
                detected_mutations.extend(unresolved_codon_records(all_results, group_pos, unresolved_reason))
                continue

            # Fix 3: a haplotype is made of accepted alleles or it is noise.
            # Accepted at each codon base: the panel ALTs the site calls
            # accepted there, and sample variants from the flank check.
            panel_called = [
                {alt for chrom, position, alt in called_alleles if chrom == codon_chrom and position == p}
                for p in codon_positions
            ]
            accepted = [
                panel_called[i] | set(sample_variants.get((codon_chrom, p), {}))
                for i, p in enumerate(codon_positions)
            ]

            def haplotype_kind(codon: str) -> str:
                differences = [i for i in range(3) if codon[i] != reference_codon[i]]
                if not all(codon[i] in accepted[i] for i in differences):
                    return "noise"
                if any(codon[i] in panel_called[i] for i in differences):
                    return "reportable"
                return "background"

            ref_values = state_counts.get(reference_codon, {"fwd": 0, "rev": 0, "bq_sum": 0})
            ref_count = ref_values["fwd"] + ref_values["rev"]
            codon_quality = {
                c: allele_quality(args, ref_values["fwd"], ref_values["rev"], v["fwd"], v["rev"],
                                  v["bq_sum"], v["fwd"] + v["rev"], v["kmer_clean"], v["kmer_total"])
                for c, v in state_counts.items()
            }
            nonref_codon_counts = {c: v["fwd"] + v["rev"] for c, v in state_counts.items() if c != reference_codon}
            credible_codon_haplotypes = {
                c for c, count in nonref_codon_counts.items()
                if haplotype_kind(c) != "noise"
                and count >= args.min_alt_reads and count / codon_depth >= args.subclonal_af
                and not codon_quality[c][1]
            }
            passing_haplotypes = 0
            ref_aa = CODON_TABLE.get(reference_codon)
            for alternate_codon, values in state_counts.items():
                if alternate_codon == reference_codon or haplotype_kind(alternate_codon) != "reportable":
                    continue
                alt_count = values["fwd"] + values["rev"]
                if not alt_count:
                    continue
                alt_aa = CODON_TABLE.get(alternate_codon)
                aa_mut = f"{ref_aa}{codon_number + 1}{alt_aa}" if ref_aa and alt_aa else ""
                aa_change = f"p.{next((k for k,v in AA3_TO_AA1.items() if v == ref_aa), ref_aa)}{codon_number + 1}{next((k for k,v in AA3_TO_AA1.items() if v == alt_aa), alt_aa)}" if aa_mut else ""
                profile = tki_profiles.get(aa_mut, {})
                consequence, resistance_class = variant_classification(aa_mut, profile)
                quality, quality_failures = codon_quality[alternate_codon]
                other_codon_alt = max(
                    (count for c, count in nonref_codon_counts.items()
                     if c != alternate_codon and c not in credible_codon_haplotypes),
                    default=0,
                )
                snr = round(alt_count / max(other_codon_alt, 1), 2)
                af = alt_count / codon_depth if codon_depth else 0.0
                if codon_depth < args.min_depth:
                    hap_status, hap_reason = "LOW_DEPTH", f"DEPTH<{args.min_depth}"
                elif alt_count < args.min_alt_reads:
                    hap_status, hap_reason = "WILDTYPE", "NO_ALT"
                elif snr < args.min_snr or quality_failures:
                    failures = ([f"LOW_SNR({snr:.1f}x<{args.min_snr:.1f}x)"] if snr < args.min_snr else [])
                    hap_status, hap_reason = "ARTIFACT", ";".join(failures + quality_failures)
                elif af >= args.min_af:
                    hap_status, hap_reason = "DETECTED", "PASS"
                elif af >= args.subclonal_af:
                    hap_status, hap_reason = "SUBCLONAL", "PASS"
                else:
                    hap_status, hap_reason = "WILDTYPE", "LOW_AF"
                record = {
                    "record_type": "CODON_HAPLOTYPE", "panel_sites_expected": panel_sites_expected,
                    "sample": args.sample, "chrom": codon_chrom, "pos": codon_start,
                    "ref": reference_codon, "alt": alternate_codon, "mutation": aa_mut,
                    "aa_change": aa_change, "name": f"NM_005157.6:c.({codon_number * 3 + 1}_{codon_number * 3 + 3}{reference_codon}>{alternate_codon})",
                    "depth": codon_depth, "ad_ref": state_counts.get(reference_codon, {}).get("fwd", 0) + state_counts.get(reference_codon, {}).get("rev", 0),
                    "ad_alt": alt_count, "vaf_percent": round(af * 100, 3), "depth_unfiltered": codon_depth,
                    "ad_alt_unfiltered": alt_count, "vaf_unfiltered_percent": round(af * 100, 3),
                    "status": hap_status, "filter_reason": hap_reason, "snr": snr, "ad_other": other_codon_alt,
                    "vaf_other_percent": round(other_codon_alt / codon_depth * 100, 3) if codon_depth else 0,
                    "other_bases": "", "credible_panel_alt_competitors": "", **quality,
                    "mean_baseq_ref": round(ref_values["bq_sum"] / ref_count, 1) if ref_count else "NA",
                    "consequence": consequence, "resistance_class": resistance_class,
                    "codon_context_status": "RESOLVED", "imatinib": profile.get("Imatinib", "Unknown"),
                    "dasatinib": profile.get("Dasatinib", "Unknown"), "nilotinib": profile.get("Nilotinib", "Unknown"),
                    "bosutinib": profile.get("Bosutinib", "Unknown"), "ponatinib": profile.get("Ponatinib", "Unknown"),
                    "asciminib": profile.get("Asciminib", "Unknown"), "ad_ref_fwd": state_counts.get(reference_codon, {}).get("fwd", 0),
                    "ad_ref_rev": state_counts.get(reference_codon, {}).get("rev", 0), "ad_alt_fwd": values["fwd"],
                    "ad_alt_rev": values["rev"], "n_refskip": "NA", "n_deletion": "NA", "n_insertion": "NA",
                    "depth_context_excluded": codon_excluded, "ad_alt_context_excluded": values["excluded"],
                    "flank_alleles_forgiven": "",
                    "depth_fusion": "NA", "ad_alt_fusion": "NA", "vaf_fusion_percent": "NA",
                    "depth_nonfusion": "NA", "ad_alt_nonfusion": "NA", "vaf_nonfusion_percent": "NA",
                }
                if hap_status in ("DETECTED", "SUBCLONAL"):
                    passing_haplotypes += 1
                if hap_status in ("DETECTED", "SUBCLONAL") or (args.include_artifacts and hap_status == "ARTIFACT"):
                    detected_mutations.append(record)
            detected_mutations = [r for r in detected_mutations if not (r.get("record_type") in {"SITE_AUDIT", "SITE_VARIANT"} and (r["chrom"], r["pos"]) in group_pos)]
            for raw in all_results:
                if (raw["chrom"], raw["pos"]) in group_pos:
                    raw["codon_context_status"] = "RESOLVED" if passing_haplotypes else "CODON_QC_FAILED"

    finally:
        bam.close()

    # Sort results
    detected_mutations.sort(key=lambda r: (r["status"] != "DETECTED", -r["vaf_percent"]))

    print(read_filter.summary(), flush=True)

    read_set_text = "fusion" if args.fusion_reads_only else "all"
    fusion_text = args.fusion_region if fusion_region is not None else "none"
    fusion_depth = sum(r["depth_fusion"] for r in all_results if isinstance(r.get("depth_fusion"), int))
    nonfusion_depth = sum(r["depth_nonfusion"] for r in all_results if isinstance(r.get("depth_nonfusion"), int))
    fusion_fraction = (
        fusion_depth / (fusion_depth + nonfusion_depth) if fusion_region is not None and fusion_depth + nonfusion_depth else None
    )
    fraction_text = f"{fusion_fraction:.3f}" if fusion_fraction is not None else "NA"
    if args.fusion_reads_only and fusion_fraction is not None and fusion_fraction < 0.5:
        print(
            f"WARNING: only {fusion_fraction:.1%} of kinase-domain reads are fusion reads; "
            "check that the BAM keeps supplementary alignments and that --fusion-region matches the reference.",
            file=sys.stderr,
        )
    sample_variant_text = ";".join(
        f"{chrom}:{position + 1}:{info['ref']}>{base}:{info['af'] * 100:.2f}%"
        for (chrom, position), bases in sorted(sample_variants.items()) for base, info in sorted(bases.items())
    ) or "none"

    # Output headers
    headers = [
        "record_type", "panel_sites_expected", "sample", "chrom", "pos", "ref", "alt", "mutation", "aa_change", "name",
        "consequence", "resistance_class", "codon_context_status",
        "depth", "ad_ref", "ad_alt", "vaf_percent", "status", "filter_reason",
        "count_A", "count_C", "count_G", "count_T",
        "pct_A", "pct_C", "pct_G", "pct_T",
        "depth_unfiltered", "ad_alt_unfiltered", "vaf_unfiltered_percent",
        "snr", "ad_other", "vaf_other_percent", "other_bases", "credible_panel_alt_competitors",
        "alt_minor_strand_ratio", "sor",
        "mean_baseq_ref", "mean_baseq_alt", "kmer_clean_percent",
        "imatinib", "dasatinib", "nilotinib", "bosutinib", "ponatinib", "asciminib",
        "ad_ref_fwd", "ad_ref_rev", "ad_alt_fwd", "ad_alt_rev", "n_refskip", "n_deletion", "n_insertion",
        "depth_context_excluded", "ad_alt_context_excluded",
        "flank_alleles_forgiven",
        "depth_fusion", "ad_alt_fusion", "vaf_fusion_percent",
        "depth_nonfusion", "ad_alt_nonfusion", "vaf_nonfusion_percent",
    ]

    # Write All Sites TSV
    args.output_all.parent.mkdir(parents=True, exist_ok=True)
    with args.output_all.open("w") as handle:
        handle.write(f"# sample={args.sample}\n")
        handle.write(f"# min_read_qs={args.min_read_qs:g} reads_dropped={len(read_filter.dropped)}\n")
        handle.write(f"# reads={read_set_text} fusion_region={fusion_text} fusion_read_fraction={fraction_text}\n")
        handle.write(f"# flank_variant_min_af={args.flank_variant_min_af} flank_sample_variants={sample_variant_text}\n")
        handle.write(f"# total_sites={len(all_results)}\n")
        handle.write(f"# detected_mutations={len([m for m in detected_mutations if m['status'] == 'DETECTED'])}\n")
        handle.write(f"# subclonal_mutations={len([m for m in detected_mutations if m['status'] == 'SUBCLONAL'])}\n")
        handle.write(f"# panel_sites_expected={panel_sites_expected}\n")
        handle.write(f"# tki_profile_version={features.get('version', 'unspecified')}\n")
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
        handle.write(f"# min_read_qs={args.min_read_qs:g} reads_dropped={len(read_filter.dropped)}\n")
        handle.write(f"# min_kmer_rate={args.min_kmer_rate} kmer_window={args.kmer_window}\n")
        handle.write(f"# reads={read_set_text} fusion_region={fusion_text} fusion_read_fraction={fraction_text}\n")
        handle.write(f"# flank_variant_min_af={args.flank_variant_min_af} flank_sample_variants={sample_variant_text}\n")
        handle.write(f"# tki_profile_source={features.get('version', 'unspecified')}\n")
        handle.write(f"# panel_sites_expected={panel_sites_expected}\n")
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
            handle.write('##FILTER=<ID=CODON_CONTEXT_UNRESOLVED,Description="Insufficient evidence to resolve the complete codon">\n')
            handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">\n')
            handle.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (REF, ALT)">\n')
            handle.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">\n')
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + args.sample + "\n")
            for rec in sorted(detected_mutations, key=lambda r: (r["chrom"], r["pos"], r["ref"], r["alt"])):
                af_dec = rec["vaf_percent"] / 100.0
                info_field = f"AA={rec['aa_change']};MUT={rec['mutation']};HGVS={rec['name']};SNR={rec['snr']};SOR={rec['sor']};BQ_ALT={rec['mean_baseq_alt']}"
                gt = "0/1" if rec["status"] in ("DETECTED", "SUBCLONAL") else "./."
                vcf_filter = "PASS" if rec["status"] in ("DETECTED", "SUBCLONAL") else "ARTIFACT"
                if rec.get("record_type") == "CODON_CONTEXT_UNRESOLVED":
                    gt, vcf_filter = "./.", "CODON_CONTEXT_UNRESOLVED"
                format_val = f"{gt}:{rec['depth']}:{rec['ad_ref']},{rec['ad_alt']}:{af_dec:.4f}"
                handle.write(f"{rec['chrom']}\t{rec['pos']}\t{rec['mutation']}\t{rec['ref']}\t{rec['alt']}\t.\t{vcf_filter}\t{info_field}\tGT:DP:AD:AF\t{format_val}\n")

    print(
        f"BCR-ABL1 Branford Caller [{args.sample}]: evaluated {len(sites)} resistance sites; "
        f"found {len(detected_mutations)} candidate mutations.",
        file=sys.stderr,
    )

def snakemake_argv(smk) -> list:
    """The command-line arguments the bcr_abl_branford_calling rule stands for."""
    p = smk.params
    argv = [
        "--bam", smk.input.bam,
        "--bed", smk.input.bed,
        "--sample", smk.wildcards.sample,
        "--output-all", smk.output.all_sites,
        "--output-mutations", smk.output.mutations,
        "--output-vcf", smk.output.vcf,
        "--features-json", smk.input.features,
    ]
    if smk.input.ref_fasta:
        argv += ["--reference-fasta", smk.input.ref_fasta]
    argv += [
        "--min-af", p.min_af,
        "--subclonal-af", p.subclonal_af,
        "--min-depth", p.min_depth,
        "--min-alt-reads", p.min_alt_reads,
        "--max-depth", p.max_depth,
        "--min-base-quality", p.min_base_quality,
        "--min-mapping-quality", p.min_mapping_quality,
        "--min-read-qs", p.min_read_qs,
        "--min-snr", p.min_snr,
        "--min-minor-strand-ratio", p.min_minor_strand_ratio,
        "--min-mean-baseq", p.min_mean_baseq,
        "--min-kmer-rate", p.min_kmer_rate,
        "--max-sor", p.max_sor,
    ]
    known_variants = list(smk.input.known_variants)
    if known_variants:
        argv += ["--known-variants", known_variants[0]]
    argv += ["--flank-variant-min-af", p.flank_variant_min_af]
    if p.fusion_region:
        argv += ["--fusion-region", p.fusion_region]
    if p.fusion_reads_only:
        argv += ["--fusion-reads-only"]
    argv += ["--kmer-window", p.kmer_window]
    argv += shlex.split(str(p.extra or ""))
    argv += ["--threads", smk.threads]
    return [str(arg) for arg in argv]


if __name__ == "__main__":
    if "snakemake" in globals():  # run by a Snakemake `script:` rule
        try:
            from snakemake_cli import run
        except ImportError:
            sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
            from snakemake_cli import run
        sys.exit(run(main, snakemake_argv(snakemake), snakemake.log))
    main()
