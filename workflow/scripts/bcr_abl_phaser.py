#!/usr/bin/env python3
"""Single-molecule read-backed phasing of BCR-ABL1 compound resistance mutations.

Oxford Nanopore reads span the entire ABL1 kinase domain (~1 kb) or entire fusion
amplicon (~2 kb). This script inspects single reads spanning pairs of detected
resistance mutations to determine whether they occur on the SAME transcript molecule
(CIS / Compound mutation) or on DIFFERENT transcript molecules (TRANS / Polyclonal clones).
"""

from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pysam

BASES = ("A", "C", "G", "T")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path, help="Input aligned BAM file")
    parser.add_argument("--mutations-tsv", required=True, type=Path, help="Detected mutations TSV from branford caller")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--output-tsv", required=True, type=Path, help="Output phasing TSV path")
    parser.add_argument("--output-json", type=Path, default=None, help="Optional output phasing JSON path")
    parser.add_argument("--min-spanning-reads", type=int, default=10, help="Minimum spanning reads required for robust phasing")
    parser.add_argument("--min-mapq", type=int, default=1, help="Minimum mapping quality")
    parser.add_argument("--threads", type=int, default=1, help="Number of worker threads for BAM decompression")
    return parser.parse_args()


def load_detected_mutations(tsv_path: Path) -> List[Dict[str, Any]]:
    mutations = []
    if not tsv_path.is_file():
        return mutations

    with tsv_path.open() as handle:
        header = None
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split("\t")
            if header is None:
                header = fields
                continue
            row = dict(zip(header, fields))
            if row.get("status") in ("DETECTED", "SUBCLONAL", "PASS"):
                try:
                    mutations.append({
                        "chrom": row["chrom"],
                        "pos": int(row["pos"]),
                        "ref": row["ref"].upper(),
                        "alt": row["alt"].upper(),
                        "mutation": row.get("mutation", f"{row['ref']}{row['pos']}{row['alt']}"),
                        "aa_change": row.get("aa_change", ""),
                        "vaf_percent": float(row.get("vaf_percent", 0.0)),
                        "depth": int(row.get("depth", 0)),
                    })
                except (ValueError, KeyError):
                    continue
    return mutations


def get_base_at_pos(read: pysam.AlignedSegment, target_1based: int) -> Optional[str]:
    """Return the base in query sequence at a 1-based reference position, or None if skipped/del/out of bounds.
    
    Uses direct CIGAR stepping for maximum performance without allocating pair tuples.
    """
    target_0based = target_1based - 1

    # Quick bounding box check
    if target_0based < read.reference_start or target_0based >= read.reference_end:
        return None

    cigartuples = read.cigartuples
    if not cigartuples:
        return None

    seq = read.query_sequence
    if not seq:
        return None

    ref_cur = read.reference_start
    query_cur = 0

    for op, length in cigartuples:
        if op in (0, 7, 8):  # Match (M), Equal (=), Diff (X)
            if ref_cur <= target_0based < ref_cur + length:
                offset = target_0based - ref_cur
                return seq[query_cur + offset].upper()
            ref_cur += length
            query_cur += length
        elif op == 1:  # Insertion (I)
            query_cur += length
        elif op == 2:  # Deletion (D)
            if ref_cur <= target_0based < ref_cur + length:
                return "DEL"
            ref_cur += length
        elif op == 3:  # RefSkip (N - intron skip)
            if ref_cur <= target_0based < ref_cur + length:
                return "REFSKIP"
            ref_cur += length
        elif op == 4:  # Soft clip (S)
            query_cur += length
        elif op in (5, 6):  # Hard clip (H), Pad (P)
            pass

        if ref_cur > target_0based:
            break

    return None



def phase_mutation_pair(
    bam: pysam.AlignmentFile,
    mut1: Dict[str, Any],
    mut2: Dict[str, Any],
    sample: str,
    min_spanning: int,
    min_mapq: int,
) -> Dict[str, Any]:
    """Interrogate all reads covering both mut1 and mut2 and compute haplotype co-occurrence."""
    chrom = mut1["chrom"]
    bam_contigs = set(bam.references)
    if chrom not in bam_contigs:
        alt_chrom = chrom[3:] if chrom.startswith("chr") else f"chr{chrom}"
        if alt_chrom in bam_contigs:
            chrom = alt_chrom

    pos1 = mut1["pos"]
    pos2 = mut2["pos"]
    if pos1 > pos2:
        pos1, pos2 = pos2, pos1
        mut1, mut2 = mut2, mut1

    ref1, alt1 = mut1["ref"], mut1["alt"]
    ref2, alt2 = mut2["ref"], mut2["alt"]

    n_11 = 0  # CIS (Alt1 + Alt2)
    n_10 = 0  # TRANS (Alt1 + Ref2)
    n_01 = 0  # TRANS (Ref1 + Alt2)
    n_00 = 0  # WT (Ref1 + Ref2)
    n_other = 0

    spanning_read_names = []

    try:
        # Fetch reads overlapping the entire interval
        for read in bam.fetch(chrom, pos1 - 1, pos2):
            if read.is_unmapped or read.is_duplicate or read.is_qcfail or (read.mapping_quality or 0) < min_mapq:
                continue

            # Verify the read actually spans both positions
            if read.reference_start > pos1 - 1 or read.reference_end < pos2:
                continue

            base1 = get_base_at_pos(read, pos1)
            base2 = get_base_at_pos(read, pos2)

            if base1 is None or base2 is None or base1 == "DEL" or base2 == "DEL":
                continue

            spanning_read_names.append(read.query_name)

            is_alt1 = (base1 == alt1)
            is_ref1 = (base1 == ref1)
            is_alt2 = (base2 == alt2)
            is_ref2 = (base2 == ref2)

            if is_alt1 and is_alt2:
                n_11 += 1
            elif is_alt1 and is_ref2:
                n_10 += 1
            elif is_ref1 and is_alt2:
                n_01 += 1
            elif is_ref1 and is_ref2:
                n_00 += 1
            else:
                n_other += 1

    except Exception as e:
        sys.stderr.write(f"Warning: error phasing {mut1['mutation']} and {mut2['mutation']}: {e}\n")

    total_spanning = n_11 + n_10 + n_01 + n_00 + n_other
    total_mutant_reads = n_11 + n_10 + n_01
    cis_fraction = (n_11 / total_mutant_reads) if total_mutant_reads > 0 else 0.0

    # Phasing classification
    if total_spanning < min_spanning:
        phase_state = "LOW_SPANNING_DEPTH"
        clinical_interp = f"Täckningen över båda positionerna ({total_spanning} läsningar) är för låg för säker fasning."
    elif total_mutant_reads < 3:
        phase_state = "INSUFFICIENT_MUTANT_READS"
        clinical_interp = f"För få mutanta läsningar ({total_mutant_reads}) för att fastställa fasningskonfiguration."
    elif cis_fraction >= 0.70:
        phase_state = "CIS (Compound)"
        clinical_interp = (
            f"Sammansatt mutation i CIS på samma BCR-ABL1-transkriptmolekyl "
            f"({n_11}/{total_mutant_reads} mutanta läsningar, {cis_fraction*100:.1f}%). "
            f"Hög risk för multiresistens mot TKI (inkl. Ponatinib/Asciminib om T315I ingår)."
        )
    elif cis_fraction <= 0.15 and (n_10 >= 2 or n_01 >= 2):
        phase_state = "TRANS (Polyclonal)"
        clinical_interp = (
            f"Polyklonala mutationer i TRANS på separata BCR-ABL1-kloner/molekyler "
            f"({n_10} läsningar för {mut1['mutation']}, {n_01} läsningar för {mut2['mutation']}). "
            f"Möjlighet till behandlingssvar vid sekventiell eller kombinerad terapi."
        )
    elif 0.15 < cis_fraction < 0.70:
        phase_state = "COEXISTING (Compound + Polyclonal)"
        clinical_interp = (
            f"Blandad klonarkitektur: både sammansatta molekyler i CIS ({n_11} läsningar) "
            f"och separata kloner i TRANS ({n_10} + {n_01} läsningar) påvisade."
        )
    else:
        phase_state = "AMBIGUOUS"
        clinical_interp = "Tvetydig fasningssignal; manuell inspektion rekommenderas."

    return {
        "sample": sample,
        "mutation_1": mut1["mutation"],
        "pos_1": pos1,
        "ref_1": ref1,
        "alt_1": alt1,
        "mutation_2": mut2["mutation"],
        "pos_2": pos2,
        "ref_2": ref2,
        "alt_2": alt2,
        "phase_state": phase_state,
        "cis_fraction_percent": round(cis_fraction * 100.0, 2),
        "cis_reads_11": n_11,
        "trans_reads_10": n_10,
        "trans_reads_01": n_01,
        "wt_reads_00": n_00,
        "other_reads": n_other,
        "total_spanning_reads": total_spanning,
        "clinical_interpretation": clinical_interp,
    }


def main() -> None:
    args = parse_args()
    if not args.bam.is_file():
        sys.exit(f"BAM file not found: {args.bam}")

    detected_muts = load_detected_mutations(args.mutations_tsv)
    
    phased_pairs: List[Dict[str, Any]] = []

    if len(detected_muts) >= 2:
        bam = pysam.AlignmentFile(str(args.bam), "rb", threads=getattr(args, "threads", 1))
        try:

            for m1, m2 in itertools.combinations(detected_muts, 2):
                pair_result = phase_mutation_pair(
                    bam, m1, m2, args.sample, args.min_spanning_reads, args.min_mapq
                )
                phased_pairs.append(pair_result)
        finally:
            bam.close()

    # Write Output TSV
    headers = [
        "sample", "mutation_1", "mutation_2", "phase_state", "cis_fraction_percent",
        "cis_reads_11", "trans_reads_10", "trans_reads_01", "wt_reads_00", "other_reads",
        "total_spanning_reads", "pos_1", "pos_2", "ref_1", "alt_1", "ref_2", "alt_2",
        "clinical_interpretation",
    ]

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w") as handle:
        handle.write(f"# sample={args.sample}\n")
        handle.write(f"# detected_mutations={len(detected_muts)}\n")
        handle.write(f"# phased_pairs={len(phased_pairs)}\n")
        handle.write("\t".join(headers) + "\n")
        for pair in phased_pairs:
            handle.write("\t".join(str(pair.get(h, "")) for h in headers) + "\n")

    # Write JSON if requested
    if args.output_json:
        args.output_json.parent.mkdir(parents=True, exist_ok=True)
        payload = {
            "sample": args.sample,
            "detected_mutations_count": len(detected_muts),
            "phased_pairs_count": len(phased_pairs),
            "phased_pairs": phased_pairs,
        }
        with args.output_json.open("w") as handle:
            json.dump(payload, handle, indent=2)

    print(
        f"BCR-ABL1 Phaser [{args.sample}]: evaluated {len(phased_pairs)} mutation pairs from "
        f"{len(detected_muts)} detected resistance mutations.",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
