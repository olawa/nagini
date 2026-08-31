#!/usr/bin/env python3
"""Count BCR-ABL1 fusion isoforms from spliced/split Oxford Nanopore alignments.

Analyzes long reads spanning BCR (chr22) and ABL1 (chr9) to identify and quantify
junction isoforms (e14a2/b3a2, e13a2/b2a2, e1a2, e19a2, and atypical variants).
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from collections import defaultdict
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pysam

# Exon-junction coordinates (1-based genomic position of the splice site).
#
# A "bcr_exons" value is the last base of that BCR exon (the donor site: on a
# spliced cDNA read, this base is immediately followed by a gap to the ABL1
# acceptor). An "abl_exons" value is the first base of that ABL1 exon (the
# acceptor site).
#
# hg19 values are derived directly from pickett_bcr_abl_pipeline's
# config/annotations_hg19/bcr_abl1.bed (the exon BED already used for
# Illumina BCR-ABL1 reporting at CGU), which gives:
#   BCR exon 1  (NM_004327.4)  chr22:23523148-23524426  -> donor  23524426
#   BCR exon 13 (NM_004327.4)  chr22:23631704-23631808  -> donor  23631808
#   BCR exon 14 (NM_004327.4)  chr22:23632526-23632600  -> donor  23632600
#   ABL1 exon 2 (NM_005157.6)  chr9:133729451-133729624 -> acceptor 133729452
#   ABL1 exon 3 (NM_005157.6)  chr9:133730188-133730483 -> acceptor 133730189
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
# e6/e19 are not present in pickett's exon BED (only exons 1, 13, 14 are
# annotated there for the isoforms it reports on) and have not been
# cross-checked against RefSeq. Treat e6a2/e19a2/e19a3 calls as provisional
# until these two are verified against an authoritative BCR exon model.

_CHR9_HG19_TO_HG38_OFFSET = 2875387   # ABL1 locus; see regenerate_branford_bed.py
_CHR22_HG19_TO_HG38_OFFSET = 342187   # BCR locus; derived from rindels'
                                       # independently curated breakpoint
                                       # clusters (m-BCR/M-BCR/mu-BCR), all
                                       # three agreeing on this exact offset.

_BCR_EXONS_HG19 = {
    "e1": 23524426,
    "e6": 23595687,   # UNVERIFIED: not present in pickett's exon BED.
    "e13": 23631808,
    "e14": 23632600,
    "e19": 23653987,  # UNVERIFIED: not present in pickett's exon BED.
}
_ABL_EXONS_HG19 = {
    "a2": 133729452,
    "a3": 133730189,
}

KNOWN_JUNCTIONS = {
    "hg19": {
        "bcr_exons": dict(_BCR_EXONS_HG19),
        "abl_exons": dict(_ABL_EXONS_HG19),
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


def parse_args() -> argparse.Namespace:
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
    return parser.parse_args()


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


def extract_fusion_events(
    bam_path: Path, build: str, tolerance: int, min_mapq: int, threads: int = 1
) -> Tuple[Dict[str, int], List[Dict[str, Any]], int]:
    """Inspect reads for chimeric BCR-ABL1 alignments and identify junction isoforms."""
    coords = KNOWN_JUNCTIONS[build]
    bcr_chroms = coords["bcr_chroms"]
    abl_chroms = coords["abl_chroms"]
    bcr_exons = coords["bcr_exons"]
    abl_exons = coords["abl_exons"]

    bam = pysam.AlignmentFile(str(bam_path), "rb", threads=threads)
    read_segments: Dict[str, List[Dict[str, Any]]] = defaultdict(list)


    total_reads_scanned = 0
    try:
        # Fetch alignments from BCR and ABL1 regions
        active_contigs = set(bam.references)
        query_regions = []
        for chrom in bcr_chroms:
            if chrom in active_contigs:
                query_regions.append((chrom, 20_000_000, 26_000_000))
        for chrom in abl_chroms:
            if chrom in active_contigs:
                query_regions.append((chrom, 128_000_000 if build == "hg38" else 130_000_000, 136_000_000))

        seen_alignments = set()
        for chrom, start, end in query_regions:
            for aln in bam.fetch(chrom, start, end):
                if aln.is_unmapped or aln.is_duplicate or (aln.mapping_quality or 0) < min_mapq:
                    continue
                
                aln_key = (aln.query_name, aln.reference_name, aln.reference_start, aln.is_supplementary)
                if aln_key in seen_alignments:
                    continue
                seen_alignments.add(aln_key)
                total_reads_scanned += 1

                is_bcr = aln.reference_name in bcr_chroms
                is_abl = aln.reference_name in abl_chroms
                if not (is_bcr or is_abl):
                    continue

                read_segments[aln.query_name].append({
                    "chrom": aln.reference_name,
                    "start": aln.reference_start + 1,  # 1-based
                    "end": aln.reference_end,
                    "is_reverse": aln.is_reverse,
                    "mapq": aln.mapping_quality,
                    "is_bcr": is_bcr,
                    "is_abl": is_abl,
                    "cigar": aln.cigarstring,
                })
    finally:
        bam.close()

    isoform_counts: Dict[str, int] = defaultdict(int)
    raw_events: List[Dict[str, Any]] = []

    for qname, segments in read_segments.items():
        has_bcr = any(s["is_bcr"] for s in segments)
        has_abl = any(s["is_abl"] for s in segments)

        if not (has_bcr and has_abl):
            continue

        bcr_segs = [s for s in segments if s["is_bcr"]]
        abl_segs = [s for s in segments if s["is_abl"]]

        # Find closest BCR junction coordinate and ABL junction coordinate
        bcr_junction_pos = max(s["end"] for s in bcr_segs)
        abl_junction_pos = min(s["start"] for s in abl_segs)

        assigned_bcr_exon = assign_exon(bcr_junction_pos, bcr_exons, tolerance)
        assigned_abl_exon = assign_exon(abl_junction_pos, abl_exons, tolerance)

        if not assigned_bcr_exon:
            # Check if start coordinate matches an exon
            bcr_start_pos = min(s["start"] for s in bcr_segs)
            assigned_bcr_exon = assign_exon(bcr_start_pos, bcr_exons, tolerance)
            if assigned_bcr_exon:
                bcr_junction_pos = bcr_start_pos

        if not assigned_abl_exon:
            abl_end_pos = max(s["end"] for s in abl_segs)
            assigned_abl_exon = assign_exon(abl_end_pos, abl_exons, tolerance)
            if assigned_abl_exon:
                abl_junction_pos = abl_end_pos

        if assigned_bcr_exon and assigned_abl_exon:
            isoform = f"{assigned_bcr_exon}{assigned_abl_exon}"
        elif assigned_bcr_exon:
            isoform = f"{assigned_bcr_exon}a2_candidate"
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

    return isoform_counts, raw_events, total_reads_scanned


def write_outputs(
    sample: str,
    isoform_counts: Dict[str, int],
    raw_events: List[Dict[str, Any]],
    output_path: Path,
    json_path: Optional[Path],
    build: str,
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

        bcr_exon = iso.split("a")[0] if "a" in iso else ""
        abl_exon = "a" + iso.split("a")[1] if "a" in iso else ""

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
        bcr_pos = sorted(bcr_pos_list)[len(bcr_pos_list) // 2] if bcr_pos_list else KNOWN_JUNCTIONS[build]["bcr_exons"].get(iso.split("a")[0], 23632600)
        abl_pos = sorted(abl_pos_list)[len(abl_pos_list) // 2] if abl_pos_list else KNOWN_JUNCTIONS[build]["abl_exons"].get("a" + iso.split("a")[1] if "a" in iso else "a2", 133729452)
        bcr_exon = evs[0]["bcr_exon"] if evs else (iso.split("a")[0] if "a" in iso else "e14")
        abl_exon = evs[0]["abl_exon"] if evs else ("a" + iso.split("a")[1] if "a" in iso else "a2")

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


def main() -> None:
    args = parse_args()
    if not args.bam.is_file():
        sys.exit(f"BAM file not found: {args.bam}")

    bam = pysam.AlignmentFile(str(args.bam), "rb", threads=getattr(args, "threads", 1))
    build = detect_genome_build(bam, args.genome_build)
    bam.close()

    isoform_counts, raw_events, total_scanned = extract_fusion_events(
        args.bam, build, args.tolerance, args.min_mapq, threads=getattr(args, "threads", 1)
    )


    write_outputs(
        args.sample,
        isoform_counts,
        raw_events,
        args.output,
        args.json_output,
        build,
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


if __name__ == "__main__":
    main()

