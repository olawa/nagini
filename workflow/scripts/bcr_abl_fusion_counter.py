#!/usr/bin/env python3
"""Count BCR-ABL1 fusion isoforms from spliced/split Oxford Nanopore alignments.

Analyzes long reads spanning BCR (chr22) and ABL1 (chr9) to identify and quantify
junction isoforms (e14a2/b3a2, e13a2/b2a2, e1a2, e19a2, and atypical variants).
"""


import argparse
import json
import math
import os
import statistics
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Set, Tuple

import pysam

try:
    from bcr_abl_read_quality import mean_read_quality, QS_TAG
except ImportError:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        from bcr_abl_read_quality import mean_read_quality, QS_TAG
    except ImportError:
        QS_TAG = "qs"
        def mean_read_quality(qualities):
            if not qualities:
                return None
            mean_error = sum(10.0 ** (-q / 10.0) for q in qualities) / len(qualities)
            return -10.0 * math.log10(mean_error) if mean_error > 0 else float(max(qualities))


# Exon-junction coordinates (1-based genomic position of the splice site).
#
# A "bcr_exons" value is the last base of that BCR exon (the donor site: on a
# spliced cDNA read, this base is immediately followed by a gap to the ABL1
# acceptor). An "abl_exons" value is the first base of that ABL1 exon (the
# acceptor site).
#
# hg19 values, 1-based, checked against the reference sequence (UCSC and
# Ensembl GRCh37 agree): every acceptor is preceded by AG, every donor
# followed by GT.
#   BCR exon 1  (NM_004327.4)  chr22:23523148-23524426  -> donor  23524426
#   BCR exon 13 (NM_004327.4)  chr22:23631704-23631808  -> donor  23631808
#   BCR exon 14 (NM_004327.4)  chr22:23632526-23632600  -> donor  23632600
#   ABL1 exon 2 (NM_005157.6)  chr9:133729451-133729624 -> acceptor 133729451
#   ABL1 exon 3 (NM_005157.6)  chr9:133730188-133730483 -> acceptor 133730188
#
# The acceptors were one base late until 2026-10-22. They came from pickett's
# exon BED, whose start column holds 1-based coordinates; adding one to
# "convert" them made the acceptor the exon's second base. Donors come from the
# end column, which is the same number in both conventions, so they were right.
#
# The previous table had e13/e14 swapped (e13 carried exon 14's donor and vice
# versa) and abl a2/a3 swapped, with a2 additionally off by ~140 kb (it
# pointed near the ABL1 gene start rather than exon 2). Every fusion call from
# that table was misclassified: e14a2 read as e13a2 and vice versa.
#
# hg38 values are these hg19 values shifted by the two fixed, independently
# verified offsets below -- not separately curated, so the two builds cannot
# drift apart the way they did before.
#
# e6 and e19 are not in pickett's exon BED. They are the 1-based ends of BCR
# exons 6 and 19 in Ensembl GRCh37 (ENST00000305877), each followed by GT in
# hg19; e1, e13 and e14 agree with the same model exactly. Until 2026-10-26 e6
# stood at 23595687 -- in intron 1, 299 bp before exon 2 and 18 kb from exon 6,
# so a true e6a2 could never be named -- and e19 was 36 bp short of its donor.

_CHR9_HG19_TO_HG38_OFFSET = 2875387   # ABL1 locus; see regenerate_branford_bed.py
_CHR22_HG19_TO_HG38_OFFSET = 342187   # BCR locus; derived from rindels'
                                       # independently curated breakpoint
                                       # clusters (m-BCR/M-BCR/mu-BCR), all
                                       # three agreeing on this exact offset.

_BCR_EXONS_HG19 = {
    "e1": 23524426,
    "e6": 23613779,
    "e13": 23631808,
    "e14": 23632600,
    "e19": 23654023,
}
_ABL_EXONS_HG19 = {
    "a2": 133729451,
    "a3": 133730188,
}

# Genomic bounds for BCR and ABL1 exon bodies / amplicon targets:
_BCR_EXON_REGIONS_HG19 = {
    "e1": (23524100, 23524550),   # spans minor primer (23524226) and donor (23524426)
    "e6": (23613650, 23613900),   # spans exon 6 body (23613719-23613779) and donor
    "e13": (23631650, 23631900),  # spans exon 13 body and donor (23631808)
    "e14": (23632480, 23632700),  # spans major primer (23632550) and donor (23632600)
    "e19": (23653850, 23654150),
}
_ABL_EXON_REGIONS_HG19 = {
    "a2": (133729451, 133729624),  # exon 2 body (174 bp)
    "a3": (133730188, 133730483),  # exon 3 body (296 bp)
}

KNOWN_JUNCTIONS = {
    "hg19": {
        "bcr_exons": dict(_BCR_EXONS_HG19),
        "abl_exons": dict(_ABL_EXONS_HG19),
        "bcr_regions": dict(_BCR_EXON_REGIONS_HG19),
        "abl_regions": dict(_ABL_EXON_REGIONS_HG19),
        "bcr_chroms": {"chr22", "22"},
        "abl_chroms": {"chr9", "9"},
    },
    "hg38": {
        "bcr_exons": {
            name: pos - _CHR22_HG19_TO_HG38_OFFSET for name, pos in _BCR_EXONS_HG19.items()
        },
        "abl_exons": {
            name: pos - _CHR9_HG19_TO_HG38_OFFSET for name, pos in _ABL_EXONS_HG19.items()
        },
        "bcr_regions": {
            name: (s - _CHR22_HG19_TO_HG38_OFFSET, e - _CHR22_HG19_TO_HG38_OFFSET)
            for name, (s, e) in _BCR_EXON_REGIONS_HG19.items()
        },
        "abl_regions": {
            name: (s - _CHR9_HG19_TO_HG38_OFFSET, e - _CHR9_HG19_TO_HG38_OFFSET)
            for name, (s, e) in _ABL_EXON_REGIONS_HG19.items()
        },
        "bcr_chroms": {"chr22", "22"},
        "abl_chroms": {"chr9", "9"},
    },
}

ISOFORM_METADATA = {
    "e14a2": {"alt_name": "b3a2", "protein": "p210", "cluster": "M-BCR", "desc": "Major CML fusion (e14a2)"},
    "e13a2": {"alt_name": "b2a2", "protein": "p210", "cluster": "M-BCR", "desc": "Major CML fusion (e13a2)"},
    "e1a2": {"alt_name": "", "protein": "p190", "cluster": "m-BCR", "desc": "Minor ALL fusion (e1a2)"},
    "e19a2": {"alt_name": "", "protein": "p230", "cluster": "mu-BCR", "desc": "Micro CNL fusion (e19a2)"},
    "e14a3": {"alt_name": "b3a3", "protein": "p205", "cluster": "M-BCR", "desc": "Atypical CML fusion (e14a3)"},
    "e13a3": {"alt_name": "b2a3", "protein": "p205", "cluster": "M-BCR", "desc": "Atypical CML fusion (e13a3)"},
    "e1a3": {"alt_name": "", "protein": "p185", "cluster": "m-BCR", "desc": "Atypical ALL fusion (e1a3)"},
    "e6a2": {"alt_name": "", "protein": "p195", "cluster": "atypical", "desc": "Atypical fusion (e6a2)"},
}


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path, help="Input aligned BAM file")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--genome-build", default="hg38", choices=["hg38", "hg19", "grch38", "grch37"], help="Reference genome build")
    parser.add_argument("--features-json", type=Path, default=None, help="Optional BCR-ABL features JSON database")
    parser.add_argument("--output", required=True, type=Path, help="Output summary TSV path")
    parser.add_argument("--json-output", type=Path, default=None, help="Optional output JSON summary path")
    parser.add_argument("--arriba-output", type=Path, default=None, help="Optional Arriba-compatible fusions TSV path for draw_fusions.R")
    parser.add_argument("--tolerance", type=int, default=150, help="Coordinate tolerance in bp to match exon junction")
    parser.add_argument("--min-mapq", type=int, default=1, help="Minimum mapping quality for read segments")
    parser.add_argument("--threads", type=int, default=1, help="Number of worker threads for BAM decompression")
    return parser.parse_args(argv)


def detect_genome_build(bam: pysam.AlignmentFile, requested_build: str) -> str:
    norm = requested_build.lower().strip()
    if norm in ("hg19", "grch37"):
        return "hg19"
    if norm in ("hg38", "grch38"):
        return "hg38"

    for contig in bam.references:
        if contig in ("chr9", "9"):
            length = bam.get_reference_length(contig)
            if length > 140_000_000:
                return "hg19" if length == 141213431 else "hg38"
    return "hg38"


def assign_exon(pos: int, exon_map: Dict[str, int], tolerance: int) -> Optional[str]:
    best_exon = None
    min_dist = tolerance + 1
    for exon_name, exon_pos in exon_map.items():
        dist = abs(pos - exon_pos)
        if dist <= tolerance and dist < min_dist:
            min_dist = dist
            best_exon = exon_name
    return best_exon


def assign_bcr_exon(
    bcr_junction_pos: int,
    bcr_segs: Sequence[Dict[str, Any]],
    bcr_exons: Dict[str, int],
    bcr_regions: Dict[str, Tuple[int, int]],
    tolerance: int,
) -> Tuple[Optional[str], int]:
    """Assign BCR exon by junction coordinate or amplicon target region overlap."""
    matched = assign_exon(bcr_junction_pos, bcr_exons, tolerance)
    if matched:
        return matched, bcr_exons.get(matched, bcr_junction_pos)

    bcr_start_pos = min(s["start"] for s in bcr_segs)
    matched_start = assign_exon(bcr_start_pos, bcr_exons, tolerance)
    if matched_start:
        return matched_start, bcr_exons.get(matched_start, bcr_start_pos)

    # Amplicon target region overlap:
    # For e1 (minor breakpoint), the entire amplicon fragment on BCR is ~200 bp.
    # Reads with soft-clipping near the chimeric junction or alignment endpoints
    # within the e1 target window unambiguously derive from e1.
    for exon_name, (r_start, r_end) in bcr_regions.items():
        if (r_start - tolerance) <= bcr_junction_pos <= (r_end + tolerance):
            return exon_name, bcr_exons.get(exon_name, bcr_junction_pos)
        for s in bcr_segs:
            s_start, s_end = s["start"], s["end"]
            overlap = max(0, min(s_end, r_end) - max(s_start, r_start))
            if overlap >= 30:
                return exon_name, bcr_exons.get(exon_name, bcr_junction_pos)

    return None, bcr_junction_pos


def assign_abl_exon(
    abl_junction_pos: int,
    abl_segs: Sequence[Dict[str, Any]],
    abl_exons: Dict[str, int],
    abl_regions: Dict[str, Tuple[int, int]],
    tolerance: int,
) -> Tuple[Optional[str], int]:
    """Assign ABL1 exon by junction coordinate and presence/absence of Exon 2 vs Exon 3."""
    matched = assign_exon(abl_junction_pos, abl_exons, tolerance)
    if matched:
        return matched, abl_exons.get(matched, abl_junction_pos)

    abl_end_pos = max(s["end"] for s in abl_segs)
    matched_end = assign_exon(abl_end_pos, abl_exons, tolerance)
    if matched_end:
        return matched_end, abl_exons.get(matched_end, abl_end_pos)

    # Check exon presence across aligned blocks / segments:
    # ABL1 exon 2 is 174 bp; ABL1 exon 3 is 296 bp.
    # In e1a3 / e13a3 / e14a3 fusions, exon 2 is deleted/skipped.
    # In e1a2 / e13a2 / e14a2 fusions, exon 2 is present.
    a2_start, a2_end = abl_regions.get("a2", (0, 0))
    a3_start, a3_end = abl_regions.get("a3", (0, 0))

    has_a2_overlap = False
    has_a3_overlap = False

    for s in abl_segs:
        blocks = s.get("blocks") or [(s["start"], s["end"])]
        for b_start, b_end in blocks:
            if a2_start and max(0, min(b_end, a2_end) - max(b_start, a2_start)) >= 20:
                has_a2_overlap = True
            if a3_start and max(0, min(b_end, a3_end) - max(b_start, a3_start)) >= 20:
                has_a3_overlap = True

    if has_a2_overlap:
        return "a2", abl_exons.get("a2", abl_junction_pos)
    if has_a3_overlap:
        return "a3", abl_exons.get("a3", abl_junction_pos)

    if a3_start and (a3_start - tolerance) <= abl_junction_pos <= (a3_end + tolerance):
        return "a3", abl_exons.get("a3", abl_junction_pos)

    return None, abl_junction_pos


def compute_n50(lengths: Sequence[int]) -> int:
    """Calculate N50 read length across a sequence of lengths."""
    if not lengths:
        return 0
    sorted_lens = sorted(lengths, reverse=True)
    half_total = sum(sorted_lens) / 2.0
    cum = 0
    for l in sorted_lens:
        cum += l
        if cum >= half_total:
            return l
    return sorted_lens[-1]


def extract_fusion_events(
    bam_path: Path, build: str, tolerance: int, min_mapq: int, threads: int = 1
) -> Tuple[Dict[str, int], List[Dict[str, Any]], int, Dict[str, Any]]:
    """Inspect reads for chimeric BCR-ABL1 alignments, quantify isoforms, and compute read stats."""
    coords = KNOWN_JUNCTIONS[build]
    bcr_chroms = coords["bcr_chroms"]
    abl_chroms = coords["abl_chroms"]
    bcr_exons = coords["bcr_exons"]
    abl_exons = coords["abl_exons"]
    bcr_regions = coords.get("bcr_regions", {})
    abl_regions = coords.get("abl_regions", {})

    bcr_start_bound = 20_000_000
    bcr_end_bound = 26_000_000
    abl_start_bound = 128_000_000 if build == "hg38" else 130_000_000
    abl_end_bound = 136_000_000

    bam = pysam.AlignmentFile(str(bam_path), "rb", threads=threads)
    read_segments: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    primary_reads: Dict[str, Tuple[int, Optional[float]]] = {}
    on_target_read_names: Set[str] = set()
    total_reads_scanned = 0

    try:
        for aln in bam.fetch(until_eof=True):
            total_reads_scanned += 1
            qname = aln.query_name

            # Primary alignment represents the unique sequenced read
            if not aln.is_secondary and not aln.is_supplementary:
                qlen = aln.infer_query_length()
                if not qlen:
                    qlen = len(aln.query_sequence) if aln.query_sequence else aln.query_length
                if qlen:
                    qs = None
                    if aln.has_tag(QS_TAG):
                        try:
                            qs = float(aln.get_tag(QS_TAG))
                        except (ValueError, TypeError):
                            qs = None
                    if qs is None and aln.query_qualities is not None and len(aln.query_qualities) > 0:
                        qs = mean_read_quality(aln.query_qualities)
                    primary_reads[qname] = (qlen, qs)

            if aln.is_unmapped or aln.is_duplicate or (aln.mapping_quality or 0) < min_mapq:
                continue

            ref_name = aln.reference_name
            ref_start = aln.reference_start
            ref_end = aln.reference_end or ref_start

            is_bcr = (ref_name in bcr_chroms) and (bcr_start_bound <= ref_end and ref_start <= bcr_end_bound)
            is_abl = (ref_name in abl_chroms) and (abl_start_bound <= ref_end and ref_start <= abl_end_bound)

            if not (is_bcr or is_abl):
                continue

            qstart = aln.query_alignment_start or 0
            qend = aln.query_alignment_end or 0
            try:
                blocks = [(s + 1, e) for s, e in aln.get_blocks()]
            except Exception:
                blocks = [(ref_start + 1, ref_end)]

            on_target_read_names.add(qname)
            read_segments[qname].append({
                "chrom": ref_name,
                "start": ref_start + 1,  # 1-based
                "end": ref_end,
                "is_reverse": aln.is_reverse,
                "mapq": aln.mapping_quality,
                "is_bcr": is_bcr,
                "is_abl": is_abl,
                "cigar": aln.cigarstring,
                "qstart": qstart,
                "qend": qend,
                "blocks": blocks,
            })
    finally:
        bam.close()

    all_lengths: List[int] = []
    all_accuracies: List[float] = []
    on_target_lengths: List[int] = []
    on_target_accuracies: List[float] = []

    for qname, (qlen, qs) in primary_reads.items():
        all_lengths.append(qlen)
        acc = ((1.0 - 10.0 ** (-qs / 10.0)) * 100.0) if qs is not None else None
        if acc is not None:
            all_accuracies.append(acc)
        if qname in on_target_read_names:
            on_target_lengths.append(qlen)
            if acc is not None:
                on_target_accuracies.append(acc)

    for qname in on_target_read_names:
        if qname not in primary_reads:
            seg_len = max((s["end"] - s["start"] + 1) for s in read_segments[qname]) if read_segments[qname] else 0
            if seg_len:
                on_target_lengths.append(seg_len)
                all_lengths.append(seg_len)

    read_stats = {
        "total_reads": len(all_lengths),
        "yield_gb": sum(all_lengths) / 1e9,
        "mean_q": round(sum(all_accuracies) / len(all_accuracies), 2) if all_accuracies else 0.0,
        "n50_bp": compute_n50(all_lengths),
        "median_bp": round(float(statistics.median(all_lengths)), 2) if all_lengths else 0.0,
        "on_target_reads": len(on_target_lengths),
        "on_target_accuracy": round(sum(on_target_accuracies) / len(on_target_accuracies), 2) if on_target_accuracies else 0.0,
        "on_target_n50_bp": compute_n50(on_target_lengths),
        "on_target_median_bp": round(float(statistics.median(on_target_lengths)), 2) if on_target_lengths else 0.0,
    }

    isoform_counts: Dict[str, int] = defaultdict(int)
    raw_events: List[Dict[str, Any]] = []

    for qname, segments in read_segments.items():
        has_bcr = any(s["is_bcr"] for s in segments)
        has_abl = any(s["is_abl"] for s in segments)

        if not (has_bcr and has_abl):
            continue

        bcr_segs = [s for s in segments if s["is_bcr"]]
        abl_segs = [s for s in segments if s["is_abl"]]

        # Select best chimeric junction pair (closest in query alignment space)
        best_bcr_seg = bcr_segs[0]
        best_abl_seg = abl_segs[0]
        if len(bcr_segs) > 1 or len(abl_segs) > 1:
            min_qdist = float("inf")
            for b_seg in bcr_segs:
                for a_seg in abl_segs:
                    qdist = max(0, b_seg["qstart"] - a_seg["qend"], a_seg["qstart"] - b_seg["qend"])
                    if qdist < min_qdist:
                        min_qdist = qdist
                        best_bcr_seg = b_seg
                        best_abl_seg = a_seg

        bcr_junction_pos = best_bcr_seg["end"] if not best_bcr_seg["is_reverse"] else best_bcr_seg["start"]
        abl_junction_pos = best_abl_seg["start"] if not best_abl_seg["is_reverse"] else best_abl_seg["end"]

        if not bcr_junction_pos:
            bcr_junction_pos = max(s["end"] for s in bcr_segs)
        if not abl_junction_pos:
            abl_junction_pos = min(s["start"] for s in abl_segs)

        assigned_bcr_exon, bcr_junction_pos = assign_bcr_exon(
            bcr_junction_pos, bcr_segs, bcr_exons, bcr_regions, tolerance
        )
        assigned_abl_exon, abl_junction_pos = assign_abl_exon(
            abl_junction_pos, abl_segs, abl_exons, abl_regions, tolerance
        )

        if assigned_bcr_exon and assigned_abl_exon:
            isoform = f"{assigned_bcr_exon}{assigned_abl_exon}"
        else:
            isoform = "atypical_bcr_abl"

        isoform_counts[isoform] += 1
        raw_events.append({
            "read_name": qname,
            "isoform": isoform,
            "bcr_pos": bcr_junction_pos,
            "abl_pos": abl_junction_pos,
            "bcr_exon": assigned_bcr_exon or "unknown",
            "abl_exon": assigned_abl_exon or "unknown",
        })

    return isoform_counts, raw_events, total_reads_scanned, read_stats



def write_outputs(
    sample: str,
    isoform_counts: Dict[str, int],
    raw_events: List[Dict[str, Any]],
    output_path: Path,
    json_path: Optional[Path],
    build: str,
    read_stats: Optional[Dict[str, Any]] = None,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    total_fusion_reads = sum(isoform_counts.values())

    rows: List[Dict[str, Any]] = []

    # Include all standard isoforms plus any detected rare/atypical ones
    all_isoforms = list(ISOFORM_METADATA.keys())
    for iso in sorted(isoform_counts.keys()):
        if iso not in all_isoforms:
            all_isoforms.append(iso)

    # Determine primary isoform
    primary_isoform = None
    if total_fusion_reads > 0:
        sorted_counts = sorted(isoform_counts.items(), key=lambda x: x[1], reverse=True)
        if sorted_counts and sorted_counts[0][1] >= 5:
            primary_isoform = sorted_counts[0][0]

    for iso in all_isoforms:
        count = isoform_counts.get(iso, 0)
        pct = (count / total_fusion_reads * 100.0) if total_fusion_reads > 0 else 0.0
        meta = ISOFORM_METADATA.get(iso, {
            "alt_name": "",
            "protein": "unknown",
            "cluster": "atypical",
            "desc": f"Atypical fusion isoform ({iso})",
        })

        if count == 0:
            status = "NOT_DETECTED"
        elif iso == primary_isoform:
            status = "PRIMARY"
        elif pct >= 1.0 or count >= 5:
            status = "SECONDARY"
        else:
            status = "LOW_CONFIDENCE"

        if iso.startswith("e") and "a" in iso[1:]:
            parts = iso.split("a", 1)
            bcr_exon = parts[0]
            abl_exon = "a" + parts[1]
        else:
            bcr_exon = ""
            abl_exon = ""

        rows.append({
            "sample": sample,
            "fusion_isoform": iso,
            "alt_name": meta["alt_name"],
            "protein_form": meta["protein"],
            "cluster": meta["cluster"],
            "bcr_exon": bcr_exon,
            "abl_exon": abl_exon,
            "read_count": count,
            "fraction_percent": round(pct, 2),
            "status": status,
            "description": meta["desc"],
        })

    rows.sort(key=lambda r: (r["status"] != "PRIMARY", r["status"] != "SECONDARY", -r["read_count"], r["fusion_isoform"]))

    # Write TSV
    headers = [
        "sample", "fusion_isoform", "alt_name", "protein_form", "cluster",
        "bcr_exon", "abl_exon", "read_count", "fraction_percent", "status", "description"
    ]
    with output_path.open("w") as handle:
        handle.write(f"# sample={sample}\n")
        handle.write(f"# genome_build={build}\n")
        handle.write(f"# total_fusion_reads={total_fusion_reads}\n")
        if read_stats:
            for k in (
                "total_reads", "yield_gb", "mean_q", "n50_bp", "median_bp",
                "on_target_reads", "on_target_accuracy", "on_target_n50_bp", "on_target_median_bp"
            ):
                if k in read_stats:
                    handle.write(f"# {k}={read_stats[k]}\n")
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(str(row[h]) for h in headers) + "\n")

    # Write JSON if requested
    if json_path:
        json_path.parent.mkdir(parents=True, exist_ok=True)
        summary_payload = {
            "sample": sample,
            "genome_build": build,
            "total_fusion_reads": total_fusion_reads,
            "primary_isoform": primary_isoform,
            "read_stats": read_stats or {},
            "isoforms": rows,
            "raw_event_count": len(raw_events),
        }
        with json_path.open("w") as handle:
            json.dump(summary_payload, handle, indent=2)



def write_arriba_tsv(
    sample: str,
    isoform_counts: Dict[str, int],
    raw_events: List[Dict[str, Any]],
    arriba_path: Path,
    build: str,
) -> None:
    """Format detected fusions into Arriba standard TSV for draw_fusions.R."""
    arriba_path.parent.mkdir(parents=True, exist_ok=True)
    headers = [
        "#gene1", "gene2", "strand1(gene/fusion)", "strand2(gene/fusion)",
        "breakpoint1", "breakpoint2", "site1", "site2", "type",
        "direction1", "direction2", "split_reads1", "split_reads2",
        "discordant_mates", "coverage1", "coverage2", "confidence",
        "closest_genomic_breakpoint1", "closest_genomic_breakpoint2",
        "filters", "fusion_transcript", "reading_frame", "peptide_sequence",
        "read_identifiers",
    ]
    rows = []
    events_by_iso: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for ev in raw_events:
        events_by_iso[ev["isoform"]].append(ev)

    for iso, count in sorted(isoform_counts.items(), key=lambda x: -x[1]):
        if count <= 0:
            continue
        evs = events_by_iso.get(iso, [])
        bcr_pos_list = [e["bcr_pos"] for e in evs]
        abl_pos_list = [e["abl_pos"] for e in evs]
        iso_bcr = iso.split("a", 1)[0] if (iso.startswith("e") and "a" in iso[1:]) else "e14"
        iso_abl = "a" + iso.split("a", 1)[1] if (iso.startswith("e") and "a" in iso[1:]) else "a2"
        bcr_pos = sorted(bcr_pos_list)[len(bcr_pos_list) // 2] if bcr_pos_list else KNOWN_JUNCTIONS[build]["bcr_exons"].get(iso_bcr, 23632600)
        abl_pos = sorted(abl_pos_list)[len(abl_pos_list) // 2] if abl_pos_list else KNOWN_JUNCTIONS[build]["abl_exons"].get(iso_abl, 133729451)
        bcr_exon = evs[0]["bcr_exon"] if evs else iso_bcr
        abl_exon = evs[0]["abl_exon"] if evs else iso_abl

        bcr_chrom = "chr22"
        abl_chrom = "chr9"

        read_ids = ",".join([e["read_name"] for e in evs[:50]])
        confidence = "high" if count >= 5 else "medium"

        row = [
            "BCR",
            "ABL1",
            "+/+",
            "+/+",
            f"{bcr_chrom}:{bcr_pos}",
            f"{abl_chrom}:{abl_pos}",
            f"exon/{bcr_exon}",
            f"exon/{abl_exon}",
            "translocation",
            "downstream",
            "upstream",
            str(count),
            "0",
            "0",
            str(count),
            str(count),
            confidence,
            ".",
            ".",
            "none",
            f"BCR({bcr_exon}):ABL1({abl_exon})",
            "in-frame",
            ".",
            read_ids,
        ]
        rows.append(row)

    with arriba_path.open("w") as handle:
        handle.write("\t".join(headers) + "\n")
        for row in rows:
            handle.write("\t".join(row) + "\n")


def main(argv=None) -> None:
    args = parse_args(argv)
    if not args.bam.is_file():
        sys.exit(f"BAM file not found: {args.bam}")

    bam = pysam.AlignmentFile(str(args.bam), "rb", threads=getattr(args, "threads", 1))
    build = detect_genome_build(bam, args.genome_build)
    bam.close()

    isoform_counts, raw_events, total_scanned, read_stats = extract_fusion_events(
        args.bam, build, args.tolerance, args.min_mapq, threads=getattr(args, "threads", 1)
    )

    write_outputs(
        args.sample,
        isoform_counts,
        raw_events,
        args.output,
        args.json_output,
        build,
        read_stats=read_stats,
    )

    if args.arriba_output:
        write_arriba_tsv(
            args.sample,
            isoform_counts,
            raw_events,
            args.arriba_output,
            build,
        )
    print(
        f"BCR-ABL1 Fusion Counter [{args.sample}]: identified {sum(isoform_counts.values())} "
        f"fusion-spanning reads across {len(isoform_counts)} isoforms.",
        file=sys.stderr,
    )

def snakemake_argv(smk) -> list:
    """The command-line arguments the bcr_abl_fusion_counts rule stands for."""
    p = smk.params
    argv = [
        "--bam", smk.input.bam,
        "--sample", smk.wildcards.sample,
        "--genome-build", p.genome_build,
        "--features-json", smk.input.features,
        "--output", smk.output.tsv,
        "--json-output", smk.output.json,
        "--arriba-output", smk.output.arriba,
        "--tolerance", p.tolerance,
        "--min-mapq", p.min_mapq,
        "--threads", smk.threads,
    ]
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
