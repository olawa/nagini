#!/usr/bin/env python3
"""Single-molecule read-backed phasing of BCR-ABL1 compound resistance mutations.

Oxford Nanopore reads span the entire ABL1 kinase domain (~1 kb) or entire fusion
amplicon (~2 kb). This script inspects single reads spanning pairs of detected
resistance mutations to determine whether they occur on the SAME transcript molecule
(CIS / Compound mutation) or on DIFFERENT transcript molecules (TRANS / Polyclonal clones).
"""


import argparse
import itertools
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import pysam

BASES = ("A", "C", "G", "T")
AA3_TO_AA1 = {
    "Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C", "Gln": "Q",
    "Glu": "E", "Gly": "G", "His": "H", "Ile": "I", "Leu": "L", "Lys": "K",
    "Met": "M", "Phe": "F", "Pro": "P", "Ser": "S", "Thr": "T", "Trp": "W",
    "Tyr": "Y", "Val": "V", "Ter": "*",
}


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", required=True, type=Path, help="Input aligned BAM file")
    parser.add_argument("--mutations-tsv", required=True, type=Path, help="Detected mutations TSV from branford caller")
    parser.add_argument("--sample", required=True, help="Sample identifier")
    parser.add_argument("--output-tsv", required=True, type=Path, help="Output phasing TSV path")
    parser.add_argument("--output-json", type=Path, default=None, help="Optional output phasing JSON path")
    parser.add_argument("--min-spanning-reads", type=int, default=10, help="Minimum spanning reads required for robust phasing")
    parser.add_argument("--min-mapq", type=int, default=1, help="Minimum mapping quality")
    parser.add_argument("--min-base-quality", type=int, default=12, help="Minimum base quality at both mutation positions")
    parser.add_argument("--threads", type=int, default=1, help="Number of worker threads for BAM decompression")
    return parser.parse_args(argv)


def load_detected_mutations(tsv_path: Path) -> List[Dict[str, Any]]:
    mutations = []
    if not tsv_path.is_file():
        raise FileNotFoundError(f"Mutation TSV not found: {tsv_path}")

    with tsv_path.open() as handle:
        header = None
        for line in handle:
            text = line.strip()
            if not text or text.startswith("#"):
                continue
            fields = text.split("\t")
            if header is None:
                header = fields
                required = {"chrom", "pos", "ref", "alt", "status"}
                missing = required - set(header)
                if missing:
                    raise ValueError(f"Mutation TSV is missing required columns: {', '.join(sorted(missing))}")
                continue
            if len(fields) != len(header):
                raise ValueError(f"Malformed mutation TSV row: expected {len(header)} fields, found {len(fields)}")
            row = dict(zip(header, fields))
            if row.get("status") in ("DETECTED", "SUBCLONAL", "PASS"):
                try:
                    mutation = {
                        "chrom": row["chrom"],
                        "pos": int(row["pos"]),
                        "ref": row["ref"].upper(),
                        "alt": row["alt"].upper(),
                        "mutation": row.get("mutation", f"{row['ref']}{row['pos']}{row['alt']}"),
                        "aa_change": row.get("aa_change", ""),
                        "record_type": row.get("record_type", ""),
                        "consequence": row.get("consequence", ""),
                        "resistance_class": row.get("resistance_class", ""),
                        "vaf_percent": float(row.get("vaf_percent", 0.0)),
                        "depth": int(row.get("depth", 0)),
                    }
                    consequence = mutation["consequence"].upper()
                    resistance_class = mutation["resistance_class"].upper()
                    record_type = mutation["record_type"].upper()
                    legacy_synonymous = _is_legacy_synonymous(
                        mutation["mutation"] or mutation["aa_change"]
                    ) or _is_legacy_synonymous(mutation["aa_change"])
                    if (record_type in {"SITE_AUDIT", "CODON_CONTEXT_UNRESOLVED"} or
                            consequence in {"SYNONYMOUS", "CODON_CONTEXT_UNRESOLVED"} or
                            resistance_class == "SYNONYMOUS" or
                            (not consequence and not resistance_class and legacy_synonymous)):
                        continue
                    mutations.append(mutation)
                except (ValueError, KeyError) as exc:
                    raise ValueError(f"Malformed detected mutation row: {text}") from exc
    if header is None:
        raise ValueError(f"Mutation TSV has no header: {tsv_path}")
    return mutations


def _is_legacy_synonymous(value: str) -> bool:
    match = re.search(r"(?:p\.)?([A-Z][a-z]{2}|[A-Z])\d+([A-Z][a-z]{2}|[A-Z])", value)
    if not match:
        return False
    ref = AA3_TO_AA1.get(match.group(1), match.group(1))
    alt = AA3_TO_AA1.get(match.group(2), match.group(2))
    return ref == alt


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
    min_base_quality: int = 12,
) -> Dict[str, Any]:
    """Interrogate all reads covering both mut1 and mut2 and compute haplotype co-occurrence."""
    chrom = mut1["chrom"]
    base_fields: Dict[str, Any] = {
        "sample": sample,
        "mutation_1": mut1["mutation"], "pos_1": mut1["pos"], "ref_1": mut1["ref"], "alt_1": mut1["alt"],
        "mutation_2": mut2["mutation"], "pos_2": mut2["pos"], "ref_2": mut2["ref"], "alt_2": mut2["alt"],
    }
    # Distinct loci on different contigs, and competing alleles at the same
    # locus, cannot be interpreted as an ordinary two-site haplotype.
    if mut1["chrom"] != mut2["chrom"]:
        return _unphaseable_result(base_fields, "UNPHASEABLE_DIFFERENT_CONTIG")
    if len(mut1["ref"]) != len(mut1["alt"]) or len(mut2["ref"]) != len(mut2["alt"]):
        return _unphaseable_result(base_fields, "UNPHASEABLE_INDEL")
    end1 = mut1["pos"] + len(mut1["ref"]) - 1
    end2 = mut2["pos"] + len(mut2["ref"]) - 1
    if not mut1["ref"] or not mut2["ref"] or set(mut1["ref"] + mut1["alt"] + mut2["ref"] + mut2["alt"]) - set(BASES):
        return _unphaseable_result(base_fields, "UNPHASEABLE_NON_SNV_ALLELE")
    if not (end1 < mut2["pos"] or end2 < mut1["pos"]):
        return _unphaseable_result(base_fields, "UNPHASEABLE_SAME_POSITION")
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
    overlapping_alignments = 0
    quality_filtered_alignments = 0
    unusable_base_alignments = 0
    duplicate_alignments = 0
    conflicting_molecules = 0
    observations: Dict[str, Tuple[Tuple[int, int, int], Optional[str], Optional[str]]] = {}

    try:
        # Fetch reads overlapping the entire interval
        interval_end = max(pos1 + len(ref1) - 1, pos2 + len(ref2) - 1)
        for read in bam.fetch(chrom, pos1 - 1, interval_end):
            if read.is_unmapped or read.is_secondary or read.is_duplicate or read.is_qcfail or (read.mapping_quality or 0) < min_mapq:
                continue

            # Verify the read actually spans both positions
            if read.reference_start > pos1 - 1 or read.reference_end < pos2 + len(ref2) - 1:
                continue

            base1 = _allele_at(read, pos1, len(ref1))
            base2 = _allele_at(read, pos2, len(ref2))

            overlapping_alignments += 1
            if (base1 is None or base2 is None or set(base1 + base2) - set(BASES)):
                unusable_base_alignments += 1
                continue
            qualities1 = _qualities_at(read, pos1, len(ref1))
            qualities2 = _qualities_at(read, pos2, len(ref2))
            if (not qualities1 or not qualities2 or min(qualities1 + qualities2) < min_base_quality):
                quality_filtered_alignments += 1
                continue

            # One molecule can have multiple primary/supplementary records.
            # Prefer primary, then higher MAPQ, then higher minimum base Q;
            # retain supplementary records when no spanning primary exists.
            rank = (int(read.is_supplementary), -int(read.mapping_quality or 0), -min(qualities1 + qualities2))
            prior = observations.get(read.query_name)
            if prior is not None:
                duplicate_alignments += 1
                if rank > prior[0]:
                    continue
                if rank == prior[0]:
                    if prior[1] is None or (base1, base2) != (prior[1], prior[2]):
                        if prior[1] is not None:
                            conflicting_molecules += 1
                        observations[read.query_name] = (rank, None, None)
                    continue
            observations[read.query_name] = (rank, base1, base2)

        for _, base1, base2 in observations.values():
            if base1 is None or base2 is None:
                continue
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
        raise RuntimeError(f"Error phasing {mut1['mutation']} and {mut2['mutation']}: {e}") from e

    # Only informative, high-quality REF/ALT observations contribute to the
    # depth denominator. Other canonical bases remain visible as a QC count.
    total_spanning = n_11 + n_10 + n_01 + n_00
    total_mutant_reads = n_11 + n_10 + n_01
    cis_fraction = (n_11 / total_mutant_reads) if total_mutant_reads > 0 else 0.0
    cis_vaf = (n_11 / total_spanning * 100.0) if total_spanning > 0 else 0.0
    cis_given_alt1 = n_11 / (n_11 + n_10) if n_11 + n_10 else 0.0
    cis_given_alt2 = n_11 / (n_11 + n_01) if n_11 + n_01 else 0.0

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
            f"Mutationerna observerades tillsammans på samma BCR-ABL1-läsningar "
            f"({n_11}/{total_mutant_reads} mutanta läsningar, {cis_fraction*100:.1f}%)."
        )
    elif (cis_given_alt1 <= 0.15 and cis_given_alt2 <= 0.15 and n_10 >= 2 and n_01 >= 2):
        phase_state = "TRANS (Polyclonal)"
        clinical_interp = (
            f"Mutant-only reads were observed for both mutations "
            f"({n_10} reads with only {mut1['mutation']}, {n_01} with only {mut2['mutation']}); "
            f"co-occurrence was uncommon among reads carrying either mutation. "
            f"Read-level phasing does not establish separate cell clones."
        )
    elif 0.15 < cis_fraction < 0.70:
        phase_state = "COEXISTING (Compound + Polyclonal)"
        clinical_interp = (
            f"Both co-occurring reads ({n_11}) and reads carrying one mutation only "
            f"({n_10} + {n_01}) were observed."
        )
    else:
        phase_state = "AMBIGUOUS"
        clinical_interp = "Ambiguous read-level phasing evidence; interpretation is uncertain."

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
        "cis_vaf_percent": round(cis_vaf, 2),
        "cis_reads_11": n_11,
        "trans_reads_10": n_10,
        "trans_reads_01": n_01,
        "wt_reads_00": n_00,
        "other_reads": n_other,
        "total_spanning_reads": total_spanning,
        "usable_spanning_reads": total_spanning,
        "overlapping_alignments": overlapping_alignments,
        "quality_filtered_alignments": quality_filtered_alignments,
        "unusable_base_alignments": unusable_base_alignments,
        "duplicate_alignments": duplicate_alignments,
        "conflicting_molecules": conflicting_molecules,
        "cis_fraction_given_alt1_percent": round(cis_given_alt1 * 100.0, 2),
        "cis_fraction_given_alt2_percent": round(cis_given_alt2 * 100.0, 2),
        "clinical_interpretation": clinical_interp,
    }


def _quality_at_pos(read: pysam.AlignedSegment, target_1based: int) -> Optional[int]:
    """Return query base quality at a reference coordinate, respecting CIGAR."""
    query_pos = _query_position_at(read, target_1based)
    if query_pos is None:
        return None
    qualities = read.query_qualities
    # Match the caller's established policy for BAMs with no quality array.
    return int(qualities[query_pos]) if qualities is not None else 30


def _query_position_at(read: pysam.AlignedSegment, target_1based: int) -> Optional[int]:
    target_0based = target_1based - 1
    if target_0based < read.reference_start or target_0based >= read.reference_end:
        return None
    query_pos = 0
    ref_pos = read.reference_start
    for op, length in read.cigartuples or ():
        if op in (0, 7, 8):
            if ref_pos <= target_0based < ref_pos + length:
                return query_pos + target_0based - ref_pos
            ref_pos += length
            query_pos += length
        elif op == 1 or op == 4:
            query_pos += length
        elif op in (2, 3):
            if ref_pos <= target_0based < ref_pos + length:
                return None
            ref_pos += length
    return None


def _allele_at(read: pysam.AlignedSegment, target_1based: int, length: int) -> Optional[str]:
    positions = [_query_position_at(read, pos) for pos in range(target_1based, target_1based + length)]
    if any(pos is None for pos in positions):
        return None
    # Consecutive reference bases must also be consecutive in the query. This
    # rejects an insertion inside a codon MNV instead of joining across it.
    if any(positions[index] != positions[0] + index for index in range(length)):
        return None
    seq = read.query_sequence
    if not seq:
        return None
    bases = [seq[pos].upper() for pos in positions]
    if any(base is None or base in ("DEL", "REFSKIP") for base in bases):
        return None
    return "".join(bases)


def _qualities_at(read: pysam.AlignedSegment, target_1based: int, length: int) -> List[int]:
    qualities = [_quality_at_pos(read, pos) for pos in range(target_1based, target_1based + length)]
    if any(q is None for q in qualities):
        return []
    return [int(q) for q in qualities]


def _unphaseable_result(fields: Dict[str, Any], state: str) -> Dict[str, Any]:
    return {
        **fields, "phase_state": state, "cis_fraction_percent": 0.0, "cis_vaf_percent": 0.0,
        "cis_reads_11": 0, "trans_reads_10": 0, "trans_reads_01": 0, "wt_reads_00": 0,
        "other_reads": 0, "total_spanning_reads": 0, "usable_spanning_reads": 0,
        "overlapping_alignments": 0, "quality_filtered_alignments": 0,
        "unusable_base_alignments": 0, "duplicate_alignments": 0,
        "conflicting_molecules": 0,
        "cis_fraction_given_alt1_percent": 0.0, "cis_fraction_given_alt2_percent": 0.0,
        "clinical_interpretation": "These mutation records do not represent two distinct, phaseable genomic positions.",
    }


def main(argv=None) -> None:
    args = parse_args(argv)
    if not args.bam.is_file():
        sys.exit(f"BAM file not found: {args.bam}")

    detected_muts = load_detected_mutations(args.mutations_tsv)
    
    phased_pairs: List[Dict[str, Any]] = []

    # Open the BAM even when fewer than two mutations were supplied, so a
    # corrupt or unreadable input cannot become an apparently clean result.
    bam = pysam.AlignmentFile(str(args.bam), "rb", threads=getattr(args, "threads", 1))
    try:
        if len(detected_muts) >= 2:
            for m1, m2 in itertools.combinations(detected_muts, 2):
                pair_result = phase_mutation_pair(
                    bam, m1, m2, args.sample, args.min_spanning_reads, args.min_mapq, args.min_base_quality
                )
                phased_pairs.append(pair_result)
    finally:
        bam.close()

    # Write Output TSV
    headers = [
        "sample", "mutation_1", "mutation_2", "phase_state", "cis_fraction_percent", "cis_vaf_percent",
        "cis_reads_11", "trans_reads_10", "trans_reads_01", "wt_reads_00", "other_reads",
        "total_spanning_reads", "usable_spanning_reads", "overlapping_alignments", "quality_filtered_alignments",
        "unusable_base_alignments", "duplicate_alignments", "conflicting_molecules",
        "cis_fraction_given_alt1_percent",
        "cis_fraction_given_alt2_percent", "pos_1", "pos_2", "ref_1", "alt_1", "ref_2", "alt_2",
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

def snakemake_argv(smk) -> list:
    """The command-line arguments the bcr_abl_phasing rule stands for."""
    p = smk.params
    argv = [
        "--bam", smk.input.bam,
        "--mutations-tsv", smk.input.mutations,
        "--sample", smk.wildcards.sample,
        "--output-tsv", smk.output.tsv,
        "--output-json", smk.output.json,
        "--min-spanning-reads", p.min_spanning_reads,
        "--min-mapq", p.min_mapq,
        "--min-base-quality", p.min_base_quality,
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
