#!/usr/bin/env python3
"""Parse nanostream amplicon stats, Rindels variant calls, and FLT3-ITD detection
into unified summary TSV/JSON reports for the AML targeted amplicon panel.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import yaml

DEFAULT_POOLS = {
    1: ["TP53_3kb_A1", "TP53_D2", "TP53_F2", "TP53_I4", "FLT3_TKD_3kb_1"],
    2: ["IDH1_1", "TP53_3kb_B1", "TP53_G4", "TP53_J3", "NPM1_5"],
    3: ["IDH2_1", "TP53_C10", "TP53_E1", "TP53_H4", "FLT3_ITD_3kb_1"],
}

AMPLICON_TO_GENE = {
    "TP53_3kb_A1": "TP53",
    "TP53_3kb_B1": "TP53",
    "TP53_C10": "TP53",
    "TP53_D2": "TP53",
    "TP53_J3": "TP53",
    "TP53_E1": "TP53",
    "TP53_F2": "TP53",
    "TP53_G4": "TP53",
    "TP53_H4": "TP53",
    "TP53_I4": "TP53",
    "FLT3_TKD_3kb_1": "FLT3",
    "FLT3_ITD_3kb_1": "FLT3",
    "IDH1_1": "IDH1",
    "IDH2_1": "IDH2",
    "NPM1_5": "NPM1",
}


def load_pool_definitions(pools_path: Optional[str] = None) -> tuple[Dict[int, List[str]], Dict[str, str]]:
    """Load pool and gene mappings from YAML or use defaults."""
    if pools_path and os.path.isfile(pools_path):
        try:
            with open(pools_path) as fh:
                data = yaml.safe_load(fh)
                pools = {int(k): v for k, v in data.get("pools", {}).items()}
                genes = data.get("amplicon_genes", AMPLICON_TO_GENE)
                return pools, genes
        except Exception:
            pass
    return DEFAULT_POOLS, AMPLICON_TO_GENE


def parse_amplicon_stats_json(
    json_path: str, pools_path: Optional[str] = None
) -> List[Dict[str, Any]]:
    """Parse nanostream amplicons JSON into a list of structured row dictionaries."""
    if not json_path or not os.path.isfile(json_path):
        return []

    try:
        with open(json_path) as fh:
            data = json.load(fh)
    except Exception as e:
        sys.stderr.write(f"Warning: failed to parse {json_path}: {e}\n")
        return []

    pools, genes = load_pool_definitions(pools_path)
    amp_to_pool = {}
    for pool_id, amps in pools.items():
        for amp in amps:
            amp_to_pool[amp] = pool_id

    amplicons_dict = data.get("amplicons", {})
    total_reads = sum(v.get("count", 0) for v in amplicons_dict.values())

    # Calculate pool totals
    pool_totals: Dict[int, int] = {}
    for name, stats in amplicons_dict.items():
        p = amp_to_pool.get(name, 0)
        pool_totals[p] = pool_totals.get(p, 0) + stats.get("count", 0)

    rows = []
    for name, stats in sorted(amplicons_dict.items()):
        count = stats.get("count", 0)
        pool_id = amp_to_pool.get(name, 0)
        p_total = pool_totals.get(pool_id, 0)

        pct_pool = (count / p_total * 100.0) if p_total > 0 else 0.0
        pct_total = (count / total_reads * 100.0) if total_reads > 0 else 0.0

        status = "PASS" if count >= 500 else ("LOW_COVERAGE" if count >= 100 else "FAIL")

        rows.append({
            "amplicon": name,
            "gene": genes.get(name, "Unknown"),
            "pool": str(pool_id) if pool_id > 0 else "N/A",
            "read_count": count,
            "pct_pool": round(pct_pool, 1),
            "pct_total": round(pct_total, 1),
            "median_length": stats.get("median_length", 0),
            "avg_qs": round(stats.get("avg_qs", 0.0), 1),
            "chrom": stats.get("chrom", ""),
            "start": stats.get("start", ""),
            "end": stats.get("end", ""),
            "status": status,
        })
    return rows


def parse_amplicon_vcf(vcf_path: str) -> List[Dict[str, Any]]:
    """Parse variants from Rindels VCF file."""
    if not vcf_path or not os.path.isfile(vcf_path):
        return []

    opener = gzip.open if vcf_path.endswith(".gz") else open
    variants = []

    try:
        with opener(vcf_path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 8:
                    continue

                chrom, pos, var_id, ref, alt, qual, filt, info = fields[:8]
                fmt = fields[8] if len(fields) > 8 else ""
                sample_data = fields[9] if len(fields) > 9 else ""

                dp = 0
                ad_alt = 0
                vaf = 0.0

                if fmt and sample_data:
                    fmt_keys = fmt.split(":")
                    sample_vals = sample_data.split(":")
                    fmt_map = dict(zip(fmt_keys, sample_vals))
                    try:
                        dp = int(fmt_map.get("DP", 0))
                    except ValueError:
                        pass
                    try:
                        ad_str = fmt_map.get("AD", "")
                        if "," in ad_str:
                            ad_alt = int(ad_str.split(",")[1])
                        elif ad_str:
                            ad_alt = int(ad_str)
                    except ValueError:
                        pass
                    try:
                        vaf = float(fmt_map.get("AF", fmt_map.get("VAF", 0.0)))
                        if vaf > 1.0:
                            vaf /= 100.0
                    except ValueError:
                        if dp > 0:
                            vaf = ad_alt / dp

                # Parse ANN or gene info
                gene = ""
                exon = ""
                consequence = ""
                hgvs_c = ""
                hgvs_p = ""

                info_map = {}
                for item in info.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info_map[k] = v

                gene = info_map.get("GENE", info_map.get("SYMBOL", ""))
                exon = info_map.get("EXON", "")
                consequence = info_map.get("CONSEQUENCE", "")
                hgvs_c = info_map.get("HGVS_C", info_map.get("HGVSC", ""))
                hgvs_p = info_map.get("HGVS_P", info_map.get("HGVSP", ""))
                ann = info_map.get("ANN", "")
                csq = info_map.get("CSQ", "")
                if ann:
                    first_ann = ann.split(",")[0].split("|")
                    if len(first_ann) >= 2 and not consequence:
                        consequence = first_ann[1]
                    if len(first_ann) >= 4 and not gene:
                        gene = first_ann[3]
                    if len(first_ann) >= 10 and not hgvs_c:
                        hgvs_c = first_ann[9]
                    if len(first_ann) >= 11 and not hgvs_p:
                        hgvs_p = first_ann[10]
                    if len(first_ann) >= 9 and not exon:
                        exon = first_ann[8]
                elif csq:
                    first_csq = csq.split(",")[0].split("|")
                    if len(first_csq) >= 2 and not consequence:
                        consequence = first_csq[1]
                    if len(first_csq) >= 4 and not gene:
                        gene = first_csq[3]

                variants.append({
                    "chrom": chrom,
                    "pos": int(pos),
                    "ref": ref,
                    "alt": alt,
                    "filter": filt,
                    "gene": gene,
                    "exon": exon,
                    "consequence": consequence,
                    "hgvs_c": hgvs_c,
                    "hgvs_p": hgvs_p,
                    "depth": dp,
                    "alt_depth": ad_alt,
                    "vaf": round(vaf * 100.0, 2),
                })
    except Exception as e:
        sys.stderr.write(f"Warning: failed parsing VCF {vcf_path}: {e}\n")

    return variants


CLINICAL_HOTSPOTS = [
    {"gene": "TP53", "codon": "p.R175", "exon": "Exon 5", "chrom": "chr17", "pos": 7675107, "ref": "C", "clinical": "Hotspot mutation (dominant negative)"},
    {"gene": "TP53", "codon": "p.R248", "exon": "Exon 7", "chrom": "chr17", "pos": 7674220, "ref": "C", "clinical": "DNA contact hotspot (loss of transactivation)"},
    {"gene": "TP53", "codon": "p.R273", "exon": "Exon 8", "chrom": "chr17", "pos": 7673802, "ref": "C", "clinical": "DNA contact hotspot (loss of transactivation)"},
    {"gene": "FLT3", "codon": "p.D835", "exon": "Exon 20", "chrom": "chr13", "pos": 28018504, "ref": "G", "clinical": "TKD activation loop mutation (TKI resistance)"},
    {"gene": "NPM1", "codon": "p.W288", "exon": "Exon 11", "chrom": "chr5", "pos": 171410544, "ref": "T", "clinical": "Exon 11 4bp frameshift insertion (NPM1c+)"},
    {"gene": "IDH1", "codon": "p.R132", "exon": "Exon 4", "chrom": "chr2", "pos": 208248388, "ref": "C", "clinical": "Active site neomorphic 2-HG hotspot"},
    {"gene": "IDH2", "codon": "p.R140", "exon": "Exon 4", "chrom": "chr15", "pos": 90088702, "ref": "C", "clinical": "Active site neomorphic 2-HG hotspot"},
    {"gene": "IDH2", "codon": "p.R172", "exon": "Exon 5", "chrom": "chr15", "pos": 90088461, "ref": "C", "clinical": "Active site neomorphic 2-HG hotspot"},
]


def parse_amplicon_pileup_counts(path: str) -> Dict[tuple[str, int], Dict[str, Any]]:
    """Parse Rindels pileup counts TSV for hotspot noise auditing."""
    if not path or not os.path.isfile(path):
        return {}
    pileup_map = {}
    opener = gzip.open if path.endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
            reader = None
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if reader is None:
                    reader = [p.lower() for p in parts]
                    continue
                row = dict(zip(reader, parts))
                chrom = row.get("chrom", "")
                try:
                    pos = int(row.get("pos", 0))
                except (ValueError, TypeError):
                    continue
                pileup_map[(chrom, pos)] = row
    except Exception as e:
        sys.stderr.write(f"Warning: failed parsing pileup {path}: {e}\n")
    return pileup_map


def audit_hotspots(
    pileup_map: Dict[tuple[str, int], Dict[str, Any]],
    variant_rows: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Audit background noise floor across clinical AML hotspots."""
    results = []
    for hs in CLINICAL_HOTSPOTS:
        chrom = hs["chrom"]
        pos = hs["pos"]
        ref_b = hs["ref"]
        called_var = next(
            (v for v in variant_rows if v.get("chrom") == chrom and abs(int(v.get("pos", 0)) - pos) <= 2),
            None
        )
        row_data = pileup_map.get((chrom, pos))
        if row_data:
            try:
                dp = int(row_data.get("filtered_depth", row_data.get("raw_depth", 0)))
            except (ValueError, TypeError):
                dp = 0
            base_ref = row_data.get("ref", ref_b).upper()
            try:
                other_sum = sum(int(row_data.get(b, 0)) for b in ("a", "c", "g", "t") if b.upper() != base_ref)
                other_sum += int(row_data.get("ins", 0)) + int(row_data.get("del", 0))
            except (ValueError, TypeError):
                other_sum = 0
            noise_pct = (other_sum / dp * 100.0) if dp > 0 else 0.0

            if called_var:
                mut_str = called_var.get("hgvs_p") or f"{called_var.get('ref')}>{called_var.get('alt')}"
                vaf_call = called_var.get("vaf", 0.0)
                status = f"MUTATION DETECTED ({mut_str}, {vaf_call}%)"
            elif dp < 500:
                status = f"LOW DEPTH ({dp}x)"
            elif noise_pct > 1.5:
                status = f"ELEVATED NOISE ({noise_pct:.2f}%)"
            else:
                status = f"CLEAN (Noise {noise_pct:.2f}%)"

            results.append({
                "gene": hs["gene"],
                "codon": hs["codon"],
                "exon": hs["exon"],
                "chrom": chrom,
                "pos": pos,
                "ref": base_ref,
                "depth": dp,
                "non_ref_reads": other_sum,
                "noise_pct": round(noise_pct, 2),
                "clinical": hs["clinical"],
                "status": status,
            })
        else:
            if called_var:
                mut_str = called_var.get("hgvs_p") or f"{called_var.get('ref')}>{called_var.get('alt')}"
                vaf_call = called_var.get("vaf", 0.0)
                status = f"MUTATION DETECTED ({mut_str}, {vaf_call}%)"
            else:
                status = "NOT ASSESSED"
            results.append({
                "gene": hs["gene"],
                "codon": hs["codon"],
                "exon": hs["exon"],
                "chrom": chrom,
                "pos": pos,
                "ref": ref_b,
                "depth": 0,
                "non_ref_reads": 0,
                "noise_pct": 0.0,
                "clinical": hs["clinical"],
                "status": status,
            })
    return results


def parse_flt3_itd_tsv(itd_path: str) -> List[Dict[str, Any]]:
    """Parse FLT3-ITD TSV file from rindels itd."""
    if not itd_path or not os.path.isfile(itd_path):
        return []

    rows = []
    try:
        with open(itd_path) as fh:
            header = None
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if header is None:
                    header = [p.lower() for p in parts]
                    continue
                row_dict = dict(zip(header, parts))
                chrom = row_dict.get("chrom", "")
                start = row_dict.get("dup_start", "")
                end = row_dict.get("dup_end", "")
                coords = f"{chrom}:{start}-{end}" if chrom and start and end else row_dict.get("region", row_dict.get("interval", "chr13:28032002-28035169"))
                itd_len_str = row_dict.get("dup_len") or row_dict.get("itd_length") or row_dict.get("length") or "0"
                try:
                    itd_len = int(float(itd_len_str))
                except (ValueError, TypeError):
                    itd_len = 0

                af_str = row_dict.get("af", row_dict.get("junction_af", "0.0"))
                try:
                    af_val = float(af_str)
                    if af_val > 1.0:
                        af_val_pct = af_val
                    else:
                        af_val_pct = af_val * 100.0
                except (ValueError, TypeError):
                    af_val_pct = 0.0

                is_det = (str(row_dict.get("status", "")).upper() == "DETECTED") or (af_val_pct > 1.0)
                status_str = "DETECTED" if is_det else "NEGATIVE"

                alt_u = int(row_dict.get("alt_unique_reads", 0))
                ref_u = int(row_dict.get("ref_unique_reads", 0))
                if ref_u > 0:
                    ar_val = alt_u / ref_u
                elif is_det and af_val_pct > 0 and af_val_pct < 100.0:
                    ar_val = af_val_pct / (100.0 - af_val_pct)
                else:
                    ar_val = 0.0

                in_frame_raw = str(row_dict.get("in_frame", "")).strip().lower()
                if in_frame_raw in ("true", "1", "yes"):
                    frame_str = "In-frame"
                elif in_frame_raw in ("false", "0", "no"):
                    frame_str = "Frameshift"
                else:
                    frame_str = "In-frame" if (itd_len % 3 == 0 and itd_len > 0) else ("Frameshift" if itd_len > 0 else "-")

                if is_det:
                    eln_risk = "HIGH AR (>=0.5)" if ar_val >= 0.5 else "LOW AR (<0.5)"
                else:
                    eln_risk = "Wildtype (Normal)"

                rows.append({
                    "gene": row_dict.get("gene", "FLT3"),
                    "region": coords,
                    "itd_length": str(itd_len),
                    "dup_len": str(itd_len),
                    "af": af_str,
                    "af_pct": round(af_val_pct, 2),
                    "ar": round(ar_val, 2) if is_det else "-",
                    "reading_frame": frame_str,
                    "eln_risk": eln_risk,
                    "status": status_str,
                    "inserted_seq": row_dict.get("inserted_seq", row_dict.get("seq", "")),
                    "in_frame": row_dict.get("in_frame", ""),
                    "confidence": row_dict.get("confidence", ""),
                    "alt_unique_reads": row_dict.get("alt_unique_reads", "0"),
                    "ref_unique_reads": row_dict.get("ref_unique_reads", "0"),
                })
    except Exception as e:
        sys.stderr.write(f"Warning: failed parsing ITD TSV {itd_path}: {e}\n")

    return rows


def write_amplicon_stats_tsv(rows: List[Dict[str, Any]], out_path: str) -> None:
    """Write amplicon QC table to TSV."""
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w") as fh:
        headers = [
            "amplicon", "gene", "pool", "read_count", "pct_pool", "pct_total",
            "median_length", "avg_qs", "chrom", "start", "end", "status"
        ]
        fh.write("\t".join(headers) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(h, "")) for h in headers) + "\n")


def write_amplicon_summary_tsv(
    stats_rows: List[Dict[str, Any]],
    variant_rows: List[Dict[str, Any]],
    itd_rows: List[Dict[str, Any]],
    out_path: str,
    sample: str = "",
    hotspot_rows: Optional[List[Dict[str, Any]]] = None,
) -> None:
    """Write overall clinical summary TSV."""
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write(f"# Sample: {sample}\n")
        fh.write("# Section: FLT3_ITD\n")
        fh.write("gene\tstatus\tlength\taf_pct\tar\treading_frame\teln_risk\tregion\n")
        if itd_rows:
            for itd in itd_rows:
                fh.write(
                    f"{itd.get('gene', 'FLT3')}\t{itd.get('status', 'NEGATIVE')}\t"
                    f"{itd.get('itd_length', 0)}\t{itd.get('af_pct', itd.get('af', 0))}\t"
                    f"{itd.get('ar', '-')}\t{itd.get('reading_frame', '-')}\t"
                    f"{itd.get('eln_risk', '-')}\t{itd.get('region', '')}\n"
                )
        else:
            fh.write("FLT3\tNEGATIVE\t0\t0.0\t-\t-\tWildtype (Normal)\tchr13:28032002-28035169\n")

        fh.write("\n# Section: VARIANTS\n")
        var_headers = ["chrom", "pos", "gene", "exon", "consequence", "ref", "alt", "vaf", "depth", "alt_depth", "filter"]
        fh.write("\t".join(var_headers) + "\n")
        for v in variant_rows:
            fh.write("\t".join(str(v.get(h, "")) for h in var_headers) + "\n")

        if hotspot_rows:
            fh.write("\n# Section: CLINICAL_HOTSPOT_BACKGROUND_NOISE\n")
            hs_headers = ["gene", "codon", "exon", "chrom", "pos", "ref", "depth", "non_ref_reads", "noise_pct", "status"]
            fh.write("\t".join(hs_headers) + "\n")
            for h in hotspot_rows:
                fh.write("\t".join(str(h.get(hdr, "")) for hdr in hs_headers) + "\n")

        fh.write("\n# Section: AMPLICON_COVERAGE\n")
        cov_headers = ["amplicon", "gene", "pool", "read_count", "pct_pool", "status"]
        fh.write("\t".join(cov_headers) + "\n")
        for s in stats_rows:
            fh.write("\t".join(str(s.get(h, "")) for h in cov_headers) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--json", type=str, help="Input nanostream amplicons stats JSON")
    parser.add_argument("--vcf", type=str, help="Input Rindels VCF file")
    parser.add_argument("--itd", type=str, help="Input Rindels ITD TSV file")
    parser.add_argument("--pileup", type=str, default="", help="Input Rindels pileup counts TSV")
    parser.add_argument("--pools", type=str, help="Pools definition YAML")
    parser.add_argument("--sample", type=str, default="", help="Sample name")
    parser.add_argument("--output-stats-tsv", type=str, help="Output amplicon stats TSV")
    parser.add_argument("--output-summary-tsv", type=str, help="Output clinical summary TSV")
    parser.add_argument("--output-summary-json", type=str, help="Output clinical summary JSON")

    args = parser.parse_args()

    stats_rows = parse_amplicon_stats_json(args.json, args.pools) if args.json else []
    variant_rows = parse_amplicon_vcf(args.vcf) if args.vcf else []
    itd_rows = parse_flt3_itd_tsv(args.itd) if args.itd else []
    pileup_map = parse_amplicon_pileup_counts(args.pileup) if args.pileup else {}
    hotspot_rows = audit_hotspots(pileup_map, variant_rows) if pileup_map else []

    if args.output_stats_tsv:
        write_amplicon_stats_tsv(stats_rows, args.output_stats_tsv)

    if args.output_summary_tsv:
        write_amplicon_summary_tsv(
            stats_rows, variant_rows, itd_rows, args.output_summary_tsv, args.sample, hotspot_rows=hotspot_rows
        )

    if args.output_summary_json:
        os.makedirs(os.path.dirname(os.path.abspath(args.output_summary_json)), exist_ok=True)
        with open(args.output_summary_json, "w") as fh:
            json.dump({
                "sample": args.sample,
                "amplicon_stats": stats_rows,
                "variants": variant_rows,
                "flt3_itd": itd_rows,
                "hotspots": hotspot_rows,
            }, fh, indent=2)


if __name__ == "__main__":
    main()
