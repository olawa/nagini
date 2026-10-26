#!/usr/bin/env python3
"""Sweep read quality scores (QS) across BCR-ABL1 kinase domain Branford panel positions.

Reads the BAM file over the ABL1 kinase domain in a single pass, extracts the read
`qs` tag (integer Q-score), stratifies per-base counts across distinct QS values,
and outputs:
1. A multi-page PDF with:
   - Page 1: Titration / Sensitivity overview curve (Median Depth, Mean Noise, Mutation VAF vs QS)
   - Pages 2..N: Full 3-panel stacked bar plot for each distinct QS score / threshold
2. An optional tidy TSV with all counts, depth, and VAFs across all sites and QS bins.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import numpy as np
    import pysam
except ImportError as err:
    # Raise, do not sys.exit: a SystemExit at import time aborts a whole pytest
    # session with INTERNALERROR and zero tests run, where an ImportError only
    # fails the modules that import this one.
    raise ImportError(
        f"missing required package ({err}). Ensure pysam, matplotlib, and numpy are installed."
    ) from err

from plot_branford_background import build_plot_figure, BASE_COLORS, STATUS_COLORS
from bcr_abl_read_quality import ReadQualityFilter

AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q", "Glu": "E",
    "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F",
    "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V", "Ter": "*",
}


def short_aa(aa: str) -> str:
    text = aa.strip()
    if text.startswith("p."):
        text = text[2:]
    for three, one in AA3_TO_AA1.items():
        text = text.replace(three, one)
    return text


def load_branford_bed(bed_path: Path) -> List[Dict[str, Any]]:
    sites: List[Dict[str, Any]] = []
    with bed_path.open() as handle:
        for line_idx, line in enumerate(handle, start=1):
            text = line.strip()
            if not text or text.startswith(("#", "track", "browser")):
                continue
            fields = text.split("\t")
            if len(fields) < 7:
                fields = text.split()
            if len(fields) < 6:
                continue
            chrom = fields[0]
            start_pos = int(fields[1])
            name = fields[3]
            ref = fields[4].upper()
            alt = fields[5].upper()
            aa_change = fields[6] if len(fields) >= 7 else name
            short_name = short_aa(aa_change)
            sites.append({
                "chrom": chrom,
                "pos0": start_pos,
                "pos1": start_pos + 1,
                "name": name,
                "ref": ref,
                "alt": alt,
                "aa_change": aa_change,
                "mutation": short_name,
            })
    return sites


def detect_genome_build_from_bam(bam_path: Path, threads: int = 1) -> str:
    """Auto-detect whether BAM was aligned to hg19 or hg38 from reference contig lengths."""
    try:
        with pysam.AlignmentFile(str(bam_path), "rb", threads=threads) as bam:
            for c in ("chr9", "9"):
                if c in bam.references:
                    length = bam.get_reference_length(c)
                    if length == 141213431 or length > 140000000:
                        return "hg19"
                    elif length == 138394717 or (135000000 < length < 140000000):
                        return "hg38"
    except Exception as exc:
        sys.stderr.write(f"Warning: could not detect build from BAM: {exc}\n")
    return "hg19"


def find_bed_file(build: str) -> Optional[Path]:
    """Find branford BED file for the given build in standard nagini locations."""
    candidates = [
        Path(f"workflow/resources/bcr_abl/branford_positions_{build}.bed"),
        Path(__file__).resolve().parent.parent / "resources" / "bcr_abl" / f"branford_positions_{build}.bed",
    ]
    for c in candidates:
        if c.is_file():
            return c
    return None


def empty_qs_counts() -> Dict[str, int]:
    return {
        "A": 0, "C": 0, "G": 0, "T": 0,
        "A_fwd": 0, "A_rev": 0,
        "C_fwd": 0, "C_rev": 0,
        "G_fwd": 0, "G_rev": 0,
        "T_fwd": 0, "T_rev": 0,
        "del": 0, "del_fwd": 0, "del_rev": 0,
        "ins": 0, "ins_fwd": 0, "ins_rev": 0,
        "total": 0, "total_fwd": 0, "total_rev": 0,
    }


CUMULATIVE_COUNTS = "_cumulative_counts"


def collect_qs_pileup(
    bam_path: Path,
    sites: List[Dict[str, Any]],
    min_base_quality: int = 0,
    min_mapping_quality: int = 0,
    cumulative_thresholds: Optional[List[int]] = None,
    threads: int = 1,
) -> Tuple[str, Dict[Any, Any]]:
    """Collect integer QS bins, evaluating cumulative filters on exact read QS."""
    target_pos = {s["pos0"] for s in sites}
    by_pos: Dict[Any, Any] = {
        s["pos0"]: defaultdict(empty_qs_counts)
        for s in sites
    }
    # Keep the public helper useful for callers that collect and then request
    # the standard QS range directly. main() passes only thresholds it needs.
    thresholds = sorted(set(range(20, 31) if cumulative_thresholds is None else cumulative_thresholds))
    by_pos[CUMULATIVE_COUNTS] = {
        pos: {threshold: empty_qs_counts() for threshold in thresholds}
        for pos in target_pos
    }
    if not target_pos:
        return "", by_pos

    with pysam.AlignmentFile(str(bam_path), "rb", threads=threads) as bam:
        ref_chrom = sites[0]["chrom"]
        chrom = ref_chrom
        if chrom not in bam.references:
            alt_c = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            if alt_c in bam.references:
                chrom = alt_c
            else:
                sys.stderr.write(f"Warning: neither {chrom} nor {alt_c} in BAM references!\n")
                return ref_chrom, by_pos

        min_pos = min(target_pos)
        max_pos = max(target_pos) + 1

        print(f"[Pileup] Scanning BAM {bam_path.name} across {chrom}:{min_pos:,}-{max_pos:,}...", file=sys.stderr)
        total_reads_seen = 0
        qs_value_counts: Counter = Counter()
        quality_values = ReadQualityFilter()
        dropped_by_threshold: Dict[int, Set[str]] = {q: set() for q in thresholds}

        for col in bam.pileup(
            chrom, min_pos, max_pos, truncate=True,
            stepper="nofilter", max_depth=2000000,
            min_base_quality=0, ignore_overlaps=False
        ):
            pos0 = col.reference_pos
            if pos0 not in target_pos:
                continue

            record_type = Tuple[Tuple[int, int], Any, Optional[float]]
            best_by_read: Dict[Tuple[str, bool, bool], record_type] = {}
            best_by_read_threshold: Dict[int, Dict[Tuple[str, bool, bool], record_type]] = {
                threshold: {} for threshold in thresholds
            }
            for read in col.pileups:
                aln = read.alignment
                if aln.is_unmapped or aln.is_duplicate or aln.is_qcfail or aln.is_secondary:
                    continue
                if (aln.mapping_quality or 0) < min_mapping_quality:
                    continue
                # Evaluate every alignment before refskip/base-quality
                # filtering, matching the caller's read-level filter order.
                qs_value = quality_values.value(aln)
                eligible_thresholds: List[int] = []
                for threshold in thresholds:
                    dropped = dropped_by_threshold[threshold]
                    if aln.query_name in dropped:
                        continue
                    if qs_value is not None and qs_value < threshold:
                        dropped.add(aln.query_name)
                        continue
                    eligible_thresholds.append(threshold)
                key = (aln.query_name, aln.is_read1, aln.is_read2)
                rank = (0 if not aln.is_supplementary else 1, -(aln.mapping_quality or 0))
                record = (rank, read, qs_value)
                if key not in best_by_read or rank < best_by_read[key][0]:
                    best_by_read[key] = record
                for threshold in eligible_thresholds:
                    previous = best_by_read_threshold[threshold].get(key)
                    if previous is None or rank < previous[0]:
                        best_by_read_threshold[threshold][key] = record

            def selected_observation(record):
                _, read, qs_value = record
                aln = read.alignment
                if read.is_refskip:
                    return None
                if read.is_del or read.query_position is None:
                    return "DEL", "rev" if aln.is_reverse else "fwd", False, qs_value
                qpos = read.query_position
                if aln.query_qualities is not None and aln.query_qualities[qpos] < min_base_quality:
                    return None
                base = aln.query_sequence[qpos].upper()
                if base not in ("A", "C", "G", "T"):
                    return None
                return base, "rev" if aln.is_reverse else "fwd", bool(getattr(read, "indel", 0) > 0), qs_value

            def add_observation(counts: Dict[str, int], base: str, strand: str, has_ins: bool) -> None:
                if base == "DEL":
                    counts["del"] += 1
                    counts[f"del_{strand}"] += 1
                else:
                    counts[base] += 1
                    counts[f"{base}_{strand}"] += 1
                    counts["total"] += 1
                    counts[f"total_{strand}"] += 1
                    if has_ins:
                        counts["ins"] += 1
                        counts[f"ins_{strand}"] += 1

            for record in best_by_read.values():
                observation = selected_observation(record)
                if observation is None:
                    continue
                base, strand, has_ins, qs_value = observation
                # Integer bins have bounded cardinality even when every read
                # has a distinct fractional tag. Cumulative counts above use
                # the original value before this display-only binning.
                qs_bin = math.floor(qs_value) if qs_value is not None else None
                counts = by_pos[pos0][qs_bin]
                add_observation(counts, base, strand, has_ins)
                if base != "DEL":
                    total_reads_seen += 1

                if qs_value is not None:
                    qs_value_counts[qs_value] += 1

            for threshold, observations in best_by_read_threshold.items():
                counts = by_pos[CUMULATIVE_COUNTS][pos0][threshold]
                for record in observations.values():
                    observation = selected_observation(record)
                    if observation is None:
                        continue
                    base, strand, has_ins, _ = observation
                    add_observation(counts, base, strand, has_ins)

        if total_reads_seen > 0:
            min_q = min(qs_value_counts) if qs_value_counts else "N/A"
            max_q = max(qs_value_counts) if qs_value_counts else "N/A"
            med_q = weighted_median(qs_value_counts) if qs_value_counts else "N/A"
            print(
                f"[Pileup] Collected {total_reads_seen:,} base calls across target sites. "
                f"Read QS range: {min_q} to {max_q} (median: {med_q}).",
                file=sys.stderr
            )
        else:
            sys.stderr.write(
                f"[Warning] 0 reads observed across {chrom}:{min_pos:,}-{max_pos:,}! "
                f"Check BAM contig names and reference genome build (hg19 vs hg38).\n"
            )

    return ref_chrom, by_pos


def weighted_median(counts: Counter) -> float:
    """Compute the median of repeated QS values without storing each read value."""
    total = sum(counts.values())
    if total == 0:
        return float("nan")
    ordered = sorted(counts.items())
    middle_low = (total - 1) // 2
    middle_high = total // 2
    cumulative = 0
    values = []
    for value, count in ordered:
        next_cumulative = cumulative + count
        if cumulative <= middle_low < next_cumulative:
            values.append(value)
        if cumulative <= middle_high < next_cumulative:
            values.append(value)
        if len(values) == 2 or (total % 2 and values):
            break
        cumulative = next_cumulative
    return float(sum(values) / len(values))


def make_bins(qs_min: int, qs_max: int, mode: str) -> List[Dict[str, Any]]:
    bins = []
    if mode in ("binned", "both"):
        # Distinct integer bins
        bins.append({
            "type": "binned",
            "label": f"QS < {qs_min}",
            "short_label": f"<{qs_min}",
            "filter_fn": lambda q: q is not None and q < qs_min,
        })
        for q in range(qs_min, qs_max + 1):
            bins.append({
                "type": "binned",
                "label": f"QS = {q}",
                "short_label": f"{q}",
                # Stored bins correspond to floor(QS); cumulative thresholds
                # are evaluated on exact values during collection.
                "filter_fn": (lambda target_q: lambda value: value is not None and math.floor(value) == target_q)(q),
            })
        bins.append({
            "type": "binned",
            "label": f"QS > {qs_max}",
            "short_label": f">{qs_max}",
            "filter_fn": lambda q: q is not None and q >= qs_max + 1,
        })
        bins.append({
            "type": "binned",
            "label": "QS unknown",
            "short_label": "unknown",
            "filter_fn": lambda q: q is None,
        })

    if mode in ("cumulative", "both"):
        # Cumulative thresholds
        for q in range(qs_min, qs_max + 1):
            bins.append({
                "type": "cumulative",
                "label": f"QS >= {q}",
                "short_label": f">={q}",
                # Unknown QS values are retained, as in ReadQualityFilter.
                "filter_fn": (lambda target_q: lambda value: value is None or value >= target_q)(q),
                "min_qs": float(q),
            })
    return bins


def rows_for_bin(
    sites: List[Dict[str, Any]],
    by_pos: Dict[Any, Any],
    bin_def: Dict[str, Any],
    sample: str,
    strand: str = "all",
) -> List[Dict[str, str]]:
    rows = []
    filter_fn = bin_def["filter_fn"]
    min_qs = bin_def.get("min_qs")

    for s in sites:
        pos0 = s["pos0"]
        counts_for_pos = by_pos.get(pos0, {})
        pos_filter_fn = filter_fn
        cumulative = by_pos.get(CUMULATIVE_COUNTS, {}).get(pos0, {})
        if min_qs is not None and int(min_qs) in cumulative:
            counts_for_pos = {None: cumulative[int(min_qs)]}
            pos_filter_fn = lambda _q: True

        def count(field: str) -> int:
            return sum(c.get(field, 0) for q, c in counts_for_pos.items() if pos_filter_fn(q))

        c_a_fwd = count("A_fwd")
        c_a_rev = count("A_rev")
        c_a = c_a_fwd + c_a_rev

        c_c_fwd = count("C_fwd")
        c_c_rev = count("C_rev")
        c_c = c_c_fwd + c_c_rev

        c_g_fwd = count("G_fwd")
        c_g_rev = count("G_rev")
        c_g = c_g_fwd + c_g_rev

        c_t_fwd = count("T_fwd")
        c_t_rev = count("T_rev")
        c_t = c_t_fwd + c_t_rev

        c_del_fwd = count("del_fwd")
        c_del_rev = count("del_rev")
        c_del = c_del_fwd + c_del_rev

        c_ins_fwd = count("ins_fwd")
        c_ins_rev = count("ins_rev")
        c_ins = c_ins_fwd + c_ins_rev

        depth_fwd = c_a_fwd + c_c_fwd + c_g_fwd + c_t_fwd
        depth_rev = c_a_rev + c_c_rev + c_g_rev + c_t_rev
        depth_all = depth_fwd + depth_rev

        ref = s["ref"]
        alt = s["alt"]

        base_counts_all = {"A": c_a, "C": c_c, "G": c_g, "T": c_t}
        base_counts_fwd = {"A": c_a_fwd, "C": c_c_fwd, "G": c_g_fwd, "T": c_t_fwd}
        base_counts_rev = {"A": c_a_rev, "C": c_c_rev, "G": c_g_rev, "T": c_t_rev}

        ad_ref_all = base_counts_all.get(ref, 0)
        ad_alt_all = base_counts_all.get(alt, 0)
        ad_ref_fwd = base_counts_fwd.get(ref, 0)
        ad_alt_fwd = base_counts_fwd.get(alt, 0)
        ad_ref_rev = base_counts_rev.get(ref, 0)
        ad_alt_rev = base_counts_rev.get(alt, 0)

        vaf_all_pct = (ad_alt_all / depth_all * 100.0) if depth_all > 0 else 0.0
        vaf_fwd_pct = (ad_alt_fwd / depth_fwd * 100.0) if depth_fwd > 0 else 0.0
        vaf_rev_pct = (ad_alt_rev / depth_rev * 100.0) if depth_rev > 0 else 0.0

        alt_minor_strand_ratio = round(min(ad_alt_fwd, ad_alt_rev) / max(ad_alt_all, 1), 3)
        r1 = (ad_alt_fwd + 1) * (ad_ref_rev + 1)
        r2 = (ad_alt_rev + 1) * (ad_ref_fwd + 1)
        sor = round(math.log(max(r1 / r2, r2 / r1)), 2)

        if strand == "fwd":
            plot_depth = depth_fwd
            p_a, p_c, p_g, p_t = c_a_fwd, c_c_fwd, c_g_fwd, c_t_fwd
            plot_del, plot_ins = c_del_fwd, c_ins_fwd
            ad_ref = ad_ref_fwd
            ad_alt = ad_alt_fwd
            vaf_pct = vaf_fwd_pct
        elif strand == "rev":
            plot_depth = depth_rev
            p_a, p_c, p_g, p_t = c_a_rev, c_c_rev, c_g_rev, c_t_rev
            plot_del, plot_ins = c_del_rev, c_ins_rev
            ad_ref = ad_ref_rev
            ad_alt = ad_alt_rev
            vaf_pct = vaf_rev_pct
        else:
            plot_depth = depth_all
            p_a, p_c, p_g, p_t = c_a, c_c, c_g, c_t
            plot_del, plot_ins = c_del, c_ins
            ad_ref = ad_ref_all
            ad_alt = ad_alt_all
            vaf_pct = vaf_all_pct

        r = {
            "record_type": "SITE_AUDIT",
            "sample": sample,
            "strand": strand.upper(),
            "mutation": s["mutation"],
            "aa_change": s["aa_change"],
            "name": s["name"],
            "chrom": s["chrom"],
            "pos": str(s["pos1"]),
            "ref": ref,
            "alt": alt,
            "depth": str(plot_depth),
            "depth_fwd": str(depth_fwd),
            "depth_rev": str(depth_rev),
            "count_A": str(p_a),
            "count_C": str(p_c),
            "count_G": str(p_g),
            "count_T": str(p_t),
            "count_A_fwd": str(c_a_fwd),
            "count_A_rev": str(c_a_rev),
            "count_C_fwd": str(c_c_fwd),
            "count_C_rev": str(c_c_rev),
            "count_G_fwd": str(c_g_fwd),
            "count_G_rev": str(c_g_rev),
            "count_T_fwd": str(c_t_fwd),
            "count_T_rev": str(c_t_rev),
            "count_del": str(plot_del),
            "count_del_fwd": str(c_del_fwd),
            "count_del_rev": str(c_del_rev),
            "count_ins": str(plot_ins),
            "count_ins_fwd": str(c_ins_fwd),
            "count_ins_rev": str(c_ins_rev),
            "pct_A": f"{(p_a / plot_depth * 100.0) if plot_depth > 0 else 0.0:.2f}",
            "pct_C": f"{(p_c / plot_depth * 100.0) if plot_depth > 0 else 0.0:.2f}",
            "pct_G": f"{(p_g / plot_depth * 100.0) if plot_depth > 0 else 0.0:.2f}",
            "pct_T": f"{(p_t / plot_depth * 100.0) if plot_depth > 0 else 0.0:.2f}",
            "ad_ref": str(ad_ref),
            "ad_ref_fwd": str(ad_ref_fwd),
            "ad_ref_rev": str(ad_ref_rev),
            "ad_alt": str(ad_alt),
            "ad_alt_fwd": str(ad_alt_fwd),
            "ad_alt_rev": str(ad_alt_rev),
            "vaf_percent": f"{vaf_pct:.2f}",
            "vaf_fwd_percent": f"{vaf_fwd_pct:.2f}",
            "vaf_rev_percent": f"{vaf_rev_pct:.2f}",
            "alt_minor_strand_ratio": f"{alt_minor_strand_ratio:.3f}",
            "sor": f"{sor:.2f}",
            "status": (
                "LOW_DEPTH" if plot_depth == 0 else
                ("DETECTED" if vaf_pct >= 1.0 else ("SUBCLONAL" if vaf_pct >= 0.2 else "NEGATIVE"))
            ),
            "filter_reason": "LOW_DEPTH" if plot_depth == 0 else ("PASS" if vaf_pct >= 0.2 else "NOISE"),
            "qs_bin": bin_def["label"],
            "qs_type": bin_def["type"],
        }
        rows.append(r)
    return rows


def row_allele_key(row: Dict[str, str]) -> Tuple[str, int, str, str]:
    return row["chrom"], int(row["pos"]) - 1, row["ref"].upper(), row["alt"].upper()


def resolve_named_mutations(names: str, sites: List[Dict[str, Any]]) -> Set[Tuple[str, int, str, str]]:
    """Panel alleles for the mutations a user names, e.g. "Y253H,T315I".

    Accepts the short form (Y253H), the HGVS protein form (p.Tyr253His) and
    any case. A name stands for every panel allele that produces it, so F317L
    covers each nucleotide change listed for it. An unknown name is an error,
    not a silent no-op: a typo would otherwise fall back to treating the real
    mutation as background.
    """

    by_name: Dict[str, Set[Tuple[str, int, str, str]]] = {}
    for site in sites:
        key = (site["chrom"], site["pos0"], site["ref"], site["alt"])
        for alias in (site["mutation"], short_aa(site["aa_change"])):
            by_name.setdefault(alias.strip().upper(), set()).add(key)

    requested = [n.strip() for n in re.split(r"[,\s]+", names or "") if n.strip()]
    if not requested:
        raise ValueError("was given but names no mutation")
    alleles: Set[Tuple[str, int, str, str]] = set()
    unknown = []
    for name in requested:
        found = by_name.get(short_aa(name).upper())
        if found:
            alleles |= found
        else:
            unknown.append(name)
    if unknown:
        known = ", ".join(sorted(by_name, key=lambda n: (int(re.sub(r"\D", "", n) or 0), n)))
        raise ValueError(f"names something not in the panel: {', '.join(unknown)}. Panel mutations: {known}")
    return alleles


def row_allele_label(row: Dict[str, str]) -> str:
    chrom, pos0, ref, alt = row_allele_key(row)
    return f'{row["mutation"]} ({chrom}:{pos0 + 1} {ref}>{alt})'


def unique_position_rows(rows: List[Dict[str, str]]) -> List[Dict[str, str]]:
    seen: Set[Tuple[str, int]] = set()
    result = []
    for row in rows:
        position = (row["chrom"], int(row["pos"]) - 1)
        if position not in seen:
            seen.add(position)
            result.append(row)
    return result


def background_noise_values(
    rows: List[Dict[str, str]],
    candidate_alleles: Set[Tuple[str, int, str, str]],
    strand: str,
) -> List[float]:
    """Return total non-reference rates once per covered, non-candidate position."""
    candidates_by_position = {(chrom, pos0) for chrom, pos0, _, _ in candidate_alleles}
    seen_positions: Set[Tuple[str, int]] = set()
    values: List[float] = []
    for row in rows:
        position = (row["chrom"], int(row["pos"]) - 1)
        if position in seen_positions:
            continue
        seen_positions.add(position)
        if position in candidates_by_position:
            continue

        depth = int(row["depth"] if strand == "all" else row[f"depth_{strand}"])
        if depth == 0:
            continue
        ref = row["ref"]
        ref_count = int(row[f"count_{ref}"] if strand == "all" else row[f"count_{ref}_{strand}"])
        values.append((depth - ref_count) / depth * 100.0)
    return values


def build_summary_figure(
    bin_summaries: List[Dict[str, Any]],
    sample: str,
    min_depth: float = 0.0,
) -> plt.Figure:
    """Build multi-panel summary curve across QS bins with strand and indel metrics.

    A noise or VAF value is drawn only when the median depth behind it reaches
    ``min_depth``; strand curves are judged on their own strand's depth, which
    on these amplicons is far from half the total. At 6x a single read moves a
    VAF by 17 points, and those bins had been drawn as points on the same curve
    as the 9,000x bins. The depth panel is always drawn: it is what shows why a
    bin is missing, and masked bins are shaded on every panel.
    """
    labels = [b["short_label"] for b in bin_summaries]
    x = np.arange(len(labels))
    depths = [b["median_depth"] for b in bin_summaries]
    depths_fwd = [b.get("median_depth_fwd", 0.0) for b in bin_summaries]
    depths_rev = [b.get("median_depth_rev", 0.0) for b in bin_summaries]

    def plotted(values):
        return [np.nan if value is None else value for value in values]

    min_depth = float(min_depth or 0.0)

    def masked(values, depth_list):
        """NaN where the depth behind a value is below min_depth."""
        return [
            np.nan if value is None or (min_depth > 0 and float(depth or 0) < min_depth) else value
            for value, depth in zip(values, depth_list)
        ]

    masked_bins = [i for i, depth in enumerate(depths) if min_depth > 0 and float(depth or 0) < min_depth]

    noises = masked([b["mean_background_noise"] for b in bin_summaries], depths)
    noises_fwd = masked([b.get("mean_background_noise_fwd") for b in bin_summaries], depths_fwd)
    noises_rev = masked([b.get("mean_background_noise_rev") for b in bin_summaries], depths_rev)
    del_rates = masked([b.get("mean_del_rate") for b in bin_summaries], depths)
    ins_rates = masked([b.get("mean_ins_rate") for b in bin_summaries], depths)

    mut_vafs: Dict[str, List[Optional[float]]] = {}
    mut_vafs_fwd: Dict[str, List[Optional[float]]] = {}
    mut_vafs_rev: Dict[str, List[Optional[float]]] = {}
    for b in bin_summaries:
        for mut_name, vaf in b["mutation_vafs"].items():
            mut_vafs.setdefault(mut_name, []).append(vaf)
        for mut_name, vaf in b.get("mutation_vafs_fwd", {}).items():
            mut_vafs_fwd.setdefault(mut_name, []).append(vaf)
        for mut_name, vaf in b.get("mutation_vafs_rev", {}).items():
            mut_vafs_rev.setdefault(mut_name, []).append(vaf)

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(max(14, len(labels) * 0.45), 10.5), sharex=True)
    fig.suptitle(f"BCR-ABL1 Kinase Domain | Quality Score (QS) Sweep Summary ({sample})", fontsize=14, fontweight="bold", y=0.98)

    # 1. Median Depth vs QS
    ax1.plot(x, depths, marker="o", color="#3470a3", linewidth=2.5, label="Median Depth (Combined)")
    ax1.fill_between(x, depths, color="#3470a3", alpha=0.12)
    ax1.plot(x, depths_fwd, marker=".", color="#0288d1", linewidth=1.5, linestyle="--", label="Forward Strand (+)")
    ax1.plot(x, depths_rev, marker=".", color="#1a237e", linewidth=1.5, linestyle=":", label="Reverse Strand (-)")
    ax1.set_ylabel("Median Depth (x)", fontsize=11, fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.5)
    import matplotlib.ticker as ticker
    ax1.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y, _: f"{int(y):,}"))
    ax1.set_title("Sequencing Coverage Retention across QS Bins", fontsize=11, loc="left", fontweight="semibold")
    ax1.legend(loc="upper right", framealpha=0.9, fontsize=9)

    # 2. Mean Background Noise (%) vs QS (Split per strand + indels)
    ax2.plot(x, noises, marker="s", color="#d95f02", linewidth=2.5, label="Mean Non-Ref Noise (%) Combined")
    ax2.plot(x, noises_fwd, marker="^", color="#0288d1", linewidth=1.6, linestyle="--", label="Forward Strand (+) Noise (%)")
    ax2.plot(x, noises_rev, marker="v", color="#2e7d32", linewidth=1.6, linestyle="--", label="Reverse Strand (-) Noise (%)")
    if any(math.isfinite(r) and r > 0 for r in del_rates):
        ax2.plot(x, del_rates, marker="x", color="#757575", linewidth=1.2, linestyle=":", label="Mean Deletion Rate (%)")
    if any(math.isfinite(r) and r > 0 for r in ins_rates):
        ax2.plot(x, ins_rates, marker="+", color="#8e24aa", linewidth=1.2, linestyle=":", label="Mean Insertion Rate (%)")
    ax2.axhline(1.0, color="gray", linestyle=":", linewidth=1, label="1.0% Threshold")
    ax2.axhline(0.2, color="lightgray", linestyle=":", linewidth=0.8, label="0.2% Subclonal")
    ax2.set_ylabel("Background Rate (%)", fontsize=11, fontweight="bold")
    all_noise_pts = noises + noises_fwd + noises_rev + del_rates + ins_rates
    finite_noise_pts = [value for value in all_noise_pts if math.isfinite(value)]
    max_noise = max(finite_noise_pts) if finite_noise_pts else 0.5
    ax2.set_ylim(0, max(1.5, max_noise * 1.25))
    ax2.grid(True, linestyle="--", alpha=0.5)
    ax2.set_title("Background Noise & Strand Bias across Negative Sites", fontsize=11, loc="left", fontweight="semibold")
    ax2.legend(loc="upper right", framealpha=0.9, fontsize=9)

    # 3. Detected Mutations VAF (%) vs QS
    colors = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628"]
    if mut_vafs:
        for idx, (m_name, v_list) in enumerate(mut_vafs.items()):
            col = colors[idx % len(colors)]
            ax3.plot(x, masked(v_list, depths), marker="^", color=col, linewidth=2.2, label=f"{m_name} (Combined)")
            if m_name in mut_vafs_fwd:
                ax3.plot(x, masked(mut_vafs_fwd[m_name], depths_fwd), color=col, linestyle="--", linewidth=1.2, alpha=0.75, label=f"{m_name} (Fwd +)")
            if m_name in mut_vafs_rev:
                ax3.plot(x, masked(mut_vafs_rev[m_name], depths_rev), color=col, linestyle=":", linewidth=1.2, alpha=0.75, label=f"{m_name} (Rev -)")
    else:
        ax3.text(0.5, 0.5, "No variants with VAF >= 1.0% detected in sample", ha="center", va="center", transform=ax3.transAxes)

    ax3.set_ylabel("Mutation VAF (%)", fontsize=11, fontweight="bold")
    # Scale to what is drawn: a masked 7 % at 150x must not set the axis.
    all_vafs = [v for v_list in mut_vafs.values() for v in masked(v_list, depths) if math.isfinite(v)]
    max_vaf = max(all_vafs) if all_vafs else 5.0
    ax3.set_ylim(0, max(5.0, max_vaf * 1.25))
    ax3.grid(True, linestyle="--", alpha=0.5)
    ax3.set_title("Detected Resistance Mutation Stability across QS Bins", fontsize=11, loc="left", fontweight="semibold")
    if mut_vafs:
        ax3.legend(loc="upper right", framealpha=0.9, fontsize=9)

    for axis in (ax1, ax2, ax3):
        for i in masked_bins:
            axis.axvspan(i - 0.5, i + 0.5, color="#bdbdbd", alpha=0.25, linewidth=0)
    if masked_bins:
        ax2.axvspan(np.nan, np.nan, color="#bdbdbd", alpha=0.25, label=f"Median depth < {min_depth:,.0f}x: not plotted")
        ax2.legend(loc="upper right", framealpha=0.9, fontsize=9)

    ax3.set_xticks(x)
    ax3.set_xticklabels(labels, rotation=45, ha="right", fontsize=9)
    ax3.set_xlabel("Read Quality Score (QS) Bin / Cutoff", fontsize=11, fontweight="bold")

    plt.tight_layout()
    return fig


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", "-b", required=True, type=Path, help="Input aligned BAM file")
    parser.add_argument(
        "--bed", type=Path, default=None,
        help="Branford positions BED file (default: auto-detected hg19 or hg38 matching BAM)"
    )
    parser.add_argument("--output-pdf", "-o", required=True, type=Path, help="Path to output multi-page PDF")
    parser.add_argument("--output-tsv", type=Path, default=None, help="Optional path to output tidy summary TSV")
    parser.add_argument("--sample", default=None, help="Sample name (defaults to BAM filename stem)")
    parser.add_argument(
        "--mode", choices=["binned", "cumulative", "both"], default="binned",
        help="Sweep mode: 'binned' (distinct integer bins: <20, 20, 21... 30, >30), "
             "'cumulative' (>=20, >=21... >=30), or 'both' (default: 'binned')"
    )
    parser.add_argument("--qs-min", type=int, default=20, help="Minimum QS value to evaluate (default: 20)")
    parser.add_argument("--qs-max", type=int, default=30, help="Maximum QS value to evaluate (default: 30)")
    parser.add_argument("--min-base-quality", type=int, default=0, help="Minimum base quality to count base (default: 0)")
    parser.add_argument("--min-mapping-quality", type=int, default=0, help="Minimum mapping quality (default: 0)")
    parser.add_argument("--zoom-max", type=float, default=5.0, help="Y-axis maximum for background zoom panel (default: 5.0)")
    parser.add_argument("--color-mode", choices=["bases", "relative"], default="bases", help="Coloring mode for bar plots")
    parser.add_argument(
        "--strand", choices=["all", "both", "fwd", "rev"], default="all",
        help="Strand breakdown for stacked bar pages: 'all' (combined only), 'both' (combined, fwd, and rev), "
             "'fwd' (forward strand only), or 'rev' (reverse strand only) (default: 'all')"
    )
    parser.add_argument(
        "--split-strands", action="store_true", default=False,
        help="Shortcut to output separate forward (+) and reverse (-) stacked bar pages for all bins (equivalent to --strand both)"
    )
    parser.add_argument(
        "--threads", "-t", type=int, default=1,
        help="Number of worker threads for BAM decompression (default: 1)"
    )
    parser.add_argument(
        "--min-plot-depth", type=float, default=500.0,
        help="Draw a bin's noise and VAF on the summary page only when the median depth behind it reaches "
             "this; strand curves use their own strand's depth. Masked bins are shaded. 0 draws everything "
             "(default: 500)"
    )
    parser.add_argument(
        "--mutations", default=None,
        help="The sample's true mutations, comma-separated (e.g. Y253H,T315I or p.Tyr253His). "
             "They are tracked on the summary page and kept out of the background. Without it, any "
             "panel allele at >= 1%% in the unfiltered reads is assumed real, which can turn a noisy "
             "position into a 'mutation' and drop it from the background it belongs to."
    )
    parsed = parser.parse_args(argv)
    if parsed.split_strands:
        parsed.strand = "both"
    return parsed


def main():
    args = parse_args()
    if args.split_strands:
        args.strand = "both"

    if not args.bam.is_file():
        sys.exit(f"Error: BAM file not found: {args.bam}")

    threads = getattr(args, "threads", 1)
    detected_build = detect_genome_build_from_bam(args.bam, threads=threads)
    print(f"[BAM] Detected reference genome build: {detected_build}", file=sys.stderr)

    bed_file = args.bed
    if bed_file is None:
        bed_file = find_bed_file(detected_build)
        if bed_file is None or not bed_file.is_file():
            sys.exit(f"Error: could not find branford_positions_{detected_build}.bed. Please specify --bed.")
        print(f"[BED] Using auto-detected BED file: {bed_file}", file=sys.stderr)
    else:
        if not bed_file.is_file():
            sys.exit(f"Error: BED file not found: {bed_file}")
        sites_test = load_branford_bed(bed_file)
        if sites_test:
            first_pos = sites_test[0]["pos0"]
            is_hg38_coords = (130000000 < first_pos < 131500000)
            is_hg19_coords = (133000000 < first_pos < 134500000)
            if detected_build == "hg19" and is_hg38_coords:
                auto_bed = find_bed_file("hg19")
                if auto_bed and auto_bed.is_file():
                    print(
                        f"[Warning] BAM is hg19 (chr9 ~141.2Mb), but provided BED '{bed_file}' has hg38 coordinates (~130.8Mb). "
                        f"Automatically switching to {auto_bed}!",
                        file=sys.stderr
                    )
                    bed_file = auto_bed
            elif detected_build == "hg38" and is_hg19_coords:
                auto_bed = find_bed_file("hg38")
                if auto_bed and auto_bed.is_file():
                    print(
                        f"[Warning] BAM is hg38 (chr9 ~138.4Mb), but provided BED '{bed_file}' has hg19 coordinates (~133.7Mb). "
                        f"Automatically switching to {auto_bed}!",
                        file=sys.stderr
                    )
                    bed_file = auto_bed

    sample = args.sample or args.bam.stem.replace("_T", "").replace(".target", "").replace(".bam", "")
    sites = load_branford_bed(bed_file)
    if not sites:
        sys.exit(f"Error: no valid sites loaded from BED {bed_file}")

    print(f"[Init] Loaded {len(sites)} Branford resistance positions from {bed_file.name}.", file=sys.stderr)
    bins = make_bins(args.qs_min, args.qs_max, args.mode)
    chrom, by_pos = collect_qs_pileup(
        bam_path=args.bam,
        sites=sites,
        min_base_quality=args.min_base_quality,
        min_mapping_quality=args.min_mapping_quality,
        cumulative_thresholds=[int(b["min_qs"]) for b in bins if "min_qs" in b],
        threads=threads,
    )

    print(f"[Sweep] Evaluating {len(bins)} QS bins ({args.mode} mode)...", file=sys.stderr)

    bin_summaries: List[Dict[str, Any]] = []
    all_tsv_rows: List[Dict[str, str]] = []
    pages_to_render: List[Tuple[Dict[str, Any], str, List[Dict[str, str]]]] = []

    # Identify primary mutations across the whole sample (to track on summary curve)
    # First check an unconstrained pool
    all_reads_fn = lambda q: True
    baseline_rows = rows_for_bin(sites, by_pos, {"filter_fn": all_reads_fn, "label": "ALL", "type": "all"}, sample, strand="all")
    named = getattr(args, "mutations", None)
    if named:
        try:
            candidate_mutations: Set[Tuple[str, int, str, str]] = resolve_named_mutations(named, sites)
        except ValueError as err:
            sys.exit(f"Error: --mutations {err}")
        source = "given with --mutations"
    else:
        candidate_mutations = {
            row_allele_key(r) for r in baseline_rows
            if float(r["vaf_percent"]) >= 1.0
        }
        if not candidate_mutations:
            candidate_mutations = {
                row_allele_key(r) for r in baseline_rows if float(r["vaf_percent"]) >= 0.5
            }
        source = "assumed from >= 1% in the unfiltered reads; pass --mutations to state them"
    chosen = sorted(row_allele_label(r) for r in baseline_rows if row_allele_key(r) in candidate_mutations)
    print(f"[Sweep] Mutations ({source}): {', '.join(chosen) or 'none'}", file=sys.stderr)

    for b in bins:
        rows_all = rows_for_bin(sites, by_pos, b, sample, strand="all")
        rows_fwd = rows_for_bin(sites, by_pos, b, sample, strand="fwd")
        rows_rev = rows_for_bin(sites, by_pos, b, sample, strand="rev")

        unique_rows_all = unique_position_rows(rows_all)
        depths_all = [float(r["depth"]) for r in unique_rows_all]
        depths_fwd = [float(r["depth_fwd"]) for r in unique_rows_all]
        depths_rev = [float(r["depth_rev"]) for r in unique_rows_all]

        med_depth_all = float(np.median(depths_all)) if depths_all else 0.0
        med_depth_fwd = float(np.median(depths_fwd)) if depths_fwd else 0.0
        med_depth_rev = float(np.median(depths_rev)) if depths_rev else 0.0

        # Background noise is the sum of all non-reference bases, once per
        # unique genomic position. Candidate alleles exclude their position
        # from this negative-site summary.
        noise_all = background_noise_values(rows_all, candidate_mutations, "all")
        noise_fwd = background_noise_values(rows_fwd, candidate_mutations, "fwd")
        noise_rev = background_noise_values(rows_rev, candidate_mutations, "rev")
        mean_noise_all = float(np.mean(noise_all)) if noise_all else None
        mean_noise_fwd = float(np.mean(noise_fwd)) if noise_fwd else None
        mean_noise_rev = float(np.mean(noise_rev)) if noise_rev else None

        del_rates: List[float] = []
        ins_rates: List[float] = []
        for row in unique_rows_all:
            depth = int(row["depth"])
            deletions = int(row["count_del"])
            if depth + deletions > 0:
                del_rates.append(deletions / (depth + deletions) * 100.0)
            if depth > 0:
                ins_rates.append(int(row["count_ins"]) / depth * 100.0)
        mean_del_rate = float(np.mean(del_rates)) if del_rates else None
        mean_ins_rate = float(np.mean(ins_rates)) if ins_rates else None

        # VAF of tracked candidate mutations
        mut_vafs_all = {
            row_allele_label(r): (float(r["vaf_percent"]) if float(r["depth"]) > 0 else None)
            for r in rows_all
            if row_allele_key(r) in candidate_mutations
        }
        mut_vafs_fwd = {
            row_allele_label(r): (float(r["vaf_fwd_percent"]) if float(r["depth_fwd"]) > 0 else None)
            for r in rows_all
            if row_allele_key(r) in candidate_mutations
        }
        mut_vafs_rev = {
            row_allele_label(r): (float(r["vaf_rev_percent"]) if float(r["depth_rev"]) > 0 else None)
            for r in rows_all
            if row_allele_key(r) in candidate_mutations
        }

        bin_summary = {
            "label": b["label"],
            "short_label": b["short_label"],
            "type": b["type"],
            "median_depth": med_depth_all,
            "median_depth_fwd": med_depth_fwd,
            "median_depth_rev": med_depth_rev,
            "mean_background_noise": mean_noise_all,
            "mean_background_noise_fwd": mean_noise_fwd,
            "mean_background_noise_rev": mean_noise_rev,
            "mean_del_rate": mean_del_rate,
            "mean_ins_rate": mean_ins_rate,
            "mutation_vafs": mut_vafs_all,
            "mutation_vafs_fwd": mut_vafs_fwd,
            "mutation_vafs_rev": mut_vafs_rev,
        }
        bin_summaries.append(bin_summary)

        if args.strand == "all":
            pages_to_render.append((b, "Combined", rows_all))
            all_tsv_rows.extend(rows_all)
        elif args.strand == "both":
            pages_to_render.append((b, "Combined", rows_all))
            pages_to_render.append((b, "Forward Strand (+)", rows_fwd))
            pages_to_render.append((b, "Reverse Strand (-)", rows_rev))
            all_tsv_rows.extend(rows_all)
            all_tsv_rows.extend(rows_fwd)
            all_tsv_rows.extend(rows_rev)
        elif args.strand == "fwd":
            pages_to_render.append((b, "Forward Strand (+)", rows_fwd))
            all_tsv_rows.extend(rows_fwd)
        elif args.strand == "rev":
            pages_to_render.append((b, "Reverse Strand (-)", rows_rev))
            all_tsv_rows.extend(rows_rev)

    # Generate multi-page PDF
    args.output_pdf.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(args.output_pdf) as pdf:
        print(f"[PDF] Rendering summary curve page...", file=sys.stderr)
        min_plot_depth = float(getattr(args, "min_plot_depth", 0) or 0)
        summary_fig = build_summary_figure(bin_summaries, sample, min_depth=min_plot_depth)
        pdf.savefig(summary_fig, dpi=300)
        plt.close(summary_fig)

        for b, page_strand_label, rows in pages_to_render:
            med_depth = int(np.median([float(r["depth"]) for r in rows]))
            title = f"BCR-ABL1 Kinase Domain ({sample}) | {b['label']} ({page_strand_label}) | Median Depth: {med_depth:,}x"
            if min_plot_depth > 0 and med_depth < min_plot_depth:
                title += f" | below {min_plot_depth:,.0f}x: not interpreted"
            print(f"[PDF] Rendering page: {b['label']} ({page_strand_label}, Depth {med_depth:,}x)...", file=sys.stderr)

            fig = build_plot_figure(
                rows=rows,
                sample=sample,
                title=title,
                zoom_max=args.zoom_max,
                color_mode=args.color_mode,
                include_depth=True,
                depth_position="top",
            )
            if fig is not None:
                pdf.savefig(fig, dpi=300)
                plt.close(fig)

    print(f"[PDF] Successfully created multi-page QS sweep report: {args.output_pdf}", file=sys.stderr)

    # Optional TSV output
    if args.output_tsv:
        headers = [
            "sample", "qs_bin", "qs_type", "strand", "mutation", "aa_change", "chrom", "pos", "ref", "alt",
            "depth", "depth_fwd", "depth_rev",
            "count_A", "count_A_fwd", "count_A_rev",
            "count_C", "count_C_fwd", "count_C_rev",
            "count_G", "count_G_fwd", "count_G_rev",
            "count_T", "count_T_fwd", "count_T_rev",
            "count_del", "count_del_fwd", "count_del_rev",
            "count_ins", "count_ins_fwd", "count_ins_rev",
            "pct_A", "pct_C", "pct_G", "pct_T",
            "ad_ref", "ad_ref_fwd", "ad_ref_rev",
            "ad_alt", "ad_alt_fwd", "ad_alt_rev",
            "vaf_percent", "vaf_fwd_percent", "vaf_rev_percent",
            "alt_minor_strand_ratio", "sor",
            "status", "filter_reason"
        ]
        args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
        with args.output_tsv.open("w") as handle:
            handle.write("\t".join(headers) + "\n")
            for r in all_tsv_rows:
                line = [r.get(h, "") for h in headers]
                handle.write("\t".join(str(val) for val in line) + "\n")
        print(f"[TSV] Saved tidy QS sweep counts summary: {args.output_tsv}", file=sys.stderr)


if __name__ == "__main__":
    main()
