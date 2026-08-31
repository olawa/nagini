#!/usr/bin/env python3
"""Force-call BCR-ABL1 kinase domain resistance mutations (Branford panel).

Performs sensitive, high-depth pileup calling at ~50 key ABL1 kinase domain codons,
excluding reference skips (intron skips in cDNA spliced alignments), calculating
exact VAFs, and annotating TKI resistance profiles (Imatinib, Dasatinib, Nilotinib,
Bosutinib, Ponatinib, Asciminib).
"""

from __future__ import annotations

import argparse
import json
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
    parser.add_argument("--min-af", type=float, default=0.01, help="Minimum AF threshold for mutation calling (default 0.01 / 1 percent)")
    parser.add_argument("--min-depth", type=int, default=50, help="Minimum total coverage depth (default 50)")
    parser.add_argument("--min-alt-reads", type=int, default=5, help="Minimum alternate supporting reads (default 5)")
    parser.add_argument("--max-depth", type=int, default=1_000_000, help="Maximum pileup depth")
    parser.add_argument("--min-base-quality", type=int, default=0, help="Minimum base quality")
    parser.add_argument("--min-mapping-quality", type=int, default=1, help="Minimum mapping quality")
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


def count_positions_batch(
    bam: pysam.AlignmentFile,
    chrom: str,
    positions_0based: set[int],
    max_depth: int,
    min_baseq: int,
    min_mapq: int,
) -> Dict[int, Dict[str, int]]:
    """Extract forward/reverse base counts across a set of coordinates in a single pileup traversal."""
    results: Dict[int, Dict[str, int]] = {}
    if not positions_0based:
        return results

    for pos_0based in positions_0based:
        results[pos_0based] = {f"{base}_{strand}": 0 for base in BASES for strand in ("fwd", "rev")}
        results[pos_0based].update({"refskip": 0, "deletion": 0, "low_baseq": 0})

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

            if min_baseq > 0 and aln.query_qualities is not None:
                if aln.query_qualities[read.query_position] < min_baseq:
                    counts["low_baseq"] += 1
                    continue

            base = aln.query_sequence[read.query_position].upper()
            if base not in BASES:
                continue

            strand = "rev" if aln.is_reverse else "fwd"
            counts[f"{base}_{strand}"] += 1

    return results


def count_position(
    bam: pysam.AlignmentFile,
    chrom: str,
    pos_1based: int,
    max_depth: int,
    min_baseq: int,
    min_mapq: int,
) -> Dict[str, int]:
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

    bam = pysam.AlignmentFile(str(args.bam), "rb", threads=getattr(args, "threads", 1))

    # Contig renaming check (e.g. chr9 vs 9)
    bam_contigs = set(bam.references)

    # Group all target sites by resolved BAM contig for single-pass scanning
    chrom_positions: Dict[str, set[int]] = {}
    for site in sites:
        chrom = site["chrom"]
        if chrom not in bam_contigs:
            alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
            if alt_chrom in bam_contigs:
                chrom = alt_chrom
        if chrom not in chrom_positions:
            chrom_positions[chrom] = set()
        chrom_positions[chrom].add(site["pos"] - 1)

    cache: Dict[Tuple[str, int], Dict[str, int]] = {}
    try:
        for chrom, pos_set in chrom_positions.items():
            if chrom in bam_contigs:
                batch_counts = count_positions_batch(
                    bam, chrom, pos_set, args.max_depth, args.min_base_quality, args.min_mapping_quality
                )
                for pos_0based, counts in batch_counts.items():
                    cache[(chrom, pos_0based + 1)] = counts
            else:
                for pos_0based in pos_set:
                    c = {f"{b}_{s}": 0 for b in BASES for s in ("fwd", "rev")}
                    c.update({"refskip": 0, "deletion": 0, "low_baseq": 0})
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

            # Classification: only candidate resistance mutations meeting thresholds are DETECTED
            if depth < args.min_depth:
                status = "LOW_DEPTH"
            elif vaf >= args.min_af and ad_alt >= args.min_alt_reads:
                status = "DETECTED"
            else:
                status = "WILDTYPE"

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
                "ad_ref_fwd": counts[f"{ref}_fwd"],
                "ad_ref_rev": counts[f"{ref}_rev"],
                "ad_alt_fwd": counts[f"{alt}_fwd"],
                "ad_alt_rev": counts[f"{alt}_rev"],
                "n_refskip": counts.get("refskip", 0),
                "n_deletion": counts.get("deletion", 0),
                "status": status,
                "imatinib": imatinib,
                "dasatinib": dasatinib,
                "nilotinib": nilotinib,
                "bosutinib": bosutinib,
                "ponatinib": ponatinib,
                "asciminib": asciminib,
            }
            all_results.append(record)

            if status == "DETECTED":
                detected_mutations.append(record)


    finally:
        bam.close()

    # Sort results
    detected_mutations.sort(key=lambda r: (r["status"] != "DETECTED", -r["vaf_percent"]))

    # Write All Sites TSV
    headers = [
        "sample", "chrom", "pos", "ref", "alt", "mutation", "aa_change", "name",
        "depth", "ad_ref", "ad_alt", "vaf_percent", "status",
        "imatinib", "dasatinib", "nilotinib", "bosutinib", "ponatinib", "asciminib",
        "ad_ref_fwd", "ad_ref_rev", "ad_alt_fwd", "ad_alt_rev", "n_refskip", "n_deletion",
    ]

    args.output_all.parent.mkdir(parents=True, exist_ok=True)
    with args.output_all.open("w") as handle:
        handle.write(f"# sample={args.sample}\n")
        handle.write(f"# total_sites={len(all_results)}\n")
        handle.write(f"# detected_mutations={len([m for m in detected_mutations if m['status'] == 'DETECTED'])}\n")
        handle.write("\t".join(headers) + "\n")
        for rec in all_results:
            handle.write("\t".join(str(rec.get(h, "")) for h in headers) + "\n")

    # Write Detected Mutations TSV
    args.output_mutations.parent.mkdir(parents=True, exist_ok=True)
    with args.output_mutations.open("w") as handle:
        handle.write(f"# sample={args.sample}\n")
        handle.write(f"# min_af={args.min_af}\n")
        handle.write(f"# min_depth={args.min_depth}\n")
        handle.write("\t".join(headers) + "\n")
        for rec in detected_mutations:
            handle.write("\t".join(str(rec.get(h, "")) for h in headers) + "\n")

    # Write VCF if requested
    if args.output_vcf:
        args.output_vcf.parent.mkdir(parents=True, exist_ok=True)
        with args.output_vcf.open("w") as handle:
            handle.write("##fileformat=VCFv4.2\n")
            handle.write(f"##source=nagini_bcr_abl_branford_caller\n")
            handle.write('##INFO=<ID=AA,Number=1,Type=String,Description="Amino acid change">\n')
            handle.write('##INFO=<ID=MUT,Number=1,Type=String,Description="Short mutation notation">\n')
            handle.write('##INFO=<ID=HGVS,Number=1,Type=String,Description="HGVS transcript annotation">\n')
            handle.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            handle.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Total depth">\n')
            handle.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths (REF, ALT)">\n')
            handle.write('##FORMAT=<ID=AF,Number=A,Type=Float,Description="Allele fraction">\n')
            handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + args.sample + "\n")
            for rec in detected_mutations:
                af_dec = rec["vaf_percent"] / 100.0
                info_field = f"AA={rec['aa_change']};MUT={rec['mutation']};HGVS={rec['name']}"
                gt = "0/1" if rec["status"] == "DETECTED" else "./."
                format_val = f"{gt}:{rec['depth']}:{rec['ad_ref']},{rec['ad_alt']}:{af_dec:.4f}"
                handle.write(f"{rec['chrom']}\t{rec['pos']}\t{rec['mutation']}\t{rec['ref']}\t{rec['alt']}\t.\tPASS\t{info_field}\tGT:DP:AD:AF\t{format_val}\n")

    print(
        f"BCR-ABL1 Branford Caller [{args.sample}]: evaluated {len(sites)} resistance sites; "
        f"found {len(detected_mutations)} candidate mutations.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
