#!/usr/bin/env python3
"""Parse nanostream amplicon stats and Rindels variant calls into unified
summary TSV/JSON reports for the AML targeted amplicon panel.

This module is the single source of truth for amplicon parsing: nagini_report.py
and nagini_summary_report.py import from here rather than keeping their own
copies, so the Excel report and the delivered TSV/JSON cannot drift apart.

FLT3-ITD is derived from the long-read Rindels VCF (insertions inside the FLT3
exon 13-15 target). The ``rindels itd`` subcommand is a short-read caller and is
deliberately not used on ONT data.
"""

from __future__ import annotations

import argparse
import gzip
import json
import os
import sys
from typing import Any, Dict, List, Optional, Sequence, Tuple

try:  # pyyaml is optional: the bundled defaults cover the production panel.
    import yaml
except ImportError:  # pragma: no cover - exercised only where pyyaml is absent
    yaml = None

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

# Read-count thresholds for the per-amplicon QC verdict.
AMPLICON_PASS_READS = 500
AMPLICON_LOW_READS = 100

# FLT3-ITD is called from insertions inside the exon 13-15 target interval of
# workflow/resources/amplicon/amplicon_targets.bed, not across the whole 3 kb
# amplicon: an insertion outside the duplication hotspot is not an ITD.
FLT3_ITD_REGION = ("chr13", 28033864, 28034426)
FLT3_ITD_MIN_LENGTH = 6
FLT3_ITD_MIN_VAF_PCT = 1.0

STATUS_DETECTED = "DETECTED"
STATUS_NEGATIVE = "NEGATIVE"
STATUS_NOT_ASSESSED = "NOT_ASSESSED"

SYMBOLIC_ALT_PREFIX = "<"


# ---------------------------------------------------------------------------
# Clinical hotspots
# ---------------------------------------------------------------------------
# GRCh38 coordinates. ``pos`` is the base the background-noise floor is read at;
# ``span`` is the inclusive interval a called variant must touch to be linked to
# this hotspot, and covers the whole codon so that any of its three bases match.
# A called indel is expanded by its own length before the overlap test, because
# VCF left-alignment can place it several bases before the codon.
#
# R172K, D835Y and the NPM1c+ duplication agree with the canonical coordinates
# used by neville_mx_amplicon (results_report_create_tables.py), which calls the
# same panel.
CLINICAL_HOTSPOTS = [
    {
        "gene": "TP53", "codon": "p.R175", "exon": "Exon 5",
        "chrom": "chr17", "pos": 7675088, "ref": "C", "span": (7675087, 7675089),
        "clinical": "Hotspot mutation (dominant negative)",
    },
    {
        "gene": "TP53", "codon": "p.R248", "exon": "Exon 7",
        "chrom": "chr17", "pos": 7674220, "ref": "C", "span": (7674219, 7674221),
        "clinical": "DNA contact hotspot (loss of transactivation)",
    },
    {
        "gene": "TP53", "codon": "p.R273", "exon": "Exon 8",
        "chrom": "chr17", "pos": 7673802, "ref": "C", "span": (7673801, 7673803),
        "clinical": "DNA contact hotspot (loss of transactivation)",
    },
    {
        "gene": "FLT3", "codon": "p.D835", "exon": "Exon 20",
        "chrom": "chr13", "pos": 28018505, "ref": "C", "span": (28018503, 28018505),
        "clinical": "TKD activation loop mutation (TKI resistance)",
    },
    {
        # Codon 288 is c.862_864; the NPM1c+ type A duplication is c.860_863dup
        # and is left-aligned to chr5:171410539 C>CTCTG, hence the wider span.
        "gene": "NPM1", "codon": "p.W288", "exon": "Exon 11",
        "chrom": "chr5", "pos": 171410542, "ref": "T", "span": (171410536, 171410546),
        "clinical": "Exon 11 4bp frameshift insertion (NPM1c+)",
    },
    {
        "gene": "IDH1", "codon": "p.R132", "exon": "Exon 4",
        "chrom": "chr2", "pos": 208248388, "ref": "C", "span": (208248387, 208248389),
        "clinical": "Active site neomorphic 2-HG hotspot",
    },
    {
        "gene": "IDH2", "codon": "p.R140", "exon": "Exon 4",
        "chrom": "chr15", "pos": 90088702, "ref": "C", "span": (90088701, 90088703),
        "clinical": "Active site neomorphic 2-HG hotspot",
    },
    {
        "gene": "IDH2", "codon": "p.R172", "exon": "Exon 4",
        "chrom": "chr15", "pos": 90088606, "ref": "C", "span": (90088605, 90088607),
        "clinical": "Active site neomorphic 2-HG hotspot",
    },
]


def _open_text(path: str):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") \
        else open(path, "r", encoding="utf-8", errors="replace")


def _to_int(value: Any, default: int = 0) -> int:
    """Best-effort int, tolerating '.', '', None and float-formatted strings."""
    try:
        return int(float(str(value).strip()))
    except (ValueError, TypeError, AttributeError):
        return default


def _to_float(value: Any, default: float = 0.0) -> float:
    """Best-effort float, tolerating '.', '', 'NA', '>100' and None."""
    try:
        s = str(value).strip()
        if s.startswith(">"):
            s = s[1:].strip()
        return float(s)
    except (ValueError, TypeError, AttributeError):
        return default


def _as_percent(value: Any) -> float:
    """Normalise an allele fraction to a percentage.

    Values at or below 1.0 are read as fractions, anything larger as a
    percentage already. Both conventions occur in the wild.
    """
    val = _to_float(value)
    if val <= 0.0:
        return 0.0
    return val if val > 1.0 else val * 100.0


def load_pool_definitions(pools_path: Optional[str] = None) -> Tuple[Dict[int, List[str]], Dict[str, str]]:
    """Load pool and gene mappings from YAML or use defaults."""
    if pools_path and os.path.isfile(pools_path) and yaml is not None:
        try:
            with open(pools_path) as fh:
                data = yaml.safe_load(fh) or {}
            pools = {int(k): v for k, v in data.get("pools", {}).items()}
            genes = data.get("amplicon_genes", AMPLICON_TO_GENE)
            if pools:
                return pools, genes
        except Exception as e:
            sys.stderr.write(f"Warning: failed to parse pools YAML {pools_path}: {e}\n")
    elif pools_path and yaml is None:
        sys.stderr.write("Warning: pyyaml unavailable, using bundled pool defaults\n")
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
    panel_reads = sum(_to_int(v.get("count", 0)) for v in amplicons_dict.values())

    pool_totals: Dict[int, int] = {}
    for name, stats in amplicons_dict.items():
        p = amp_to_pool.get(name, 0)
        pool_totals[p] = pool_totals.get(p, 0) + _to_int(stats.get("count", 0))

    rows = []
    for name, stats in sorted(amplicons_dict.items()):
        count = _to_int(stats.get("count", 0))
        pool_id = amp_to_pool.get(name, 0)
        p_total = pool_totals.get(pool_id, 0)

        pct_pool = (count / p_total * 100.0) if p_total > 0 else 0.0
        # Fraction of reads assigned to a panel amplicon, not of all run reads.
        pct_panel = (count / panel_reads * 100.0) if panel_reads > 0 else 0.0

        if count >= AMPLICON_PASS_READS:
            status = "PASS"
        elif count >= AMPLICON_LOW_READS:
            status = "LOW_COVERAGE"
        else:
            status = "FAIL"

        rows.append({
            "amplicon": name,
            "gene": genes.get(name, "Unknown"),
            "pool": str(pool_id) if pool_id > 0 else "N/A",
            "read_count": count,
            "pct_pool": round(pct_pool, 1),
            "pct_total": round(pct_panel, 1),
            "median_length": stats.get("median_length", 0),
            "avg_qs": round(_to_float(stats.get("avg_qs", 0.0)), 1),
            "chrom": stats.get("chrom", ""),
            "start": stats.get("start", ""),
            "end": stats.get("end", ""),
            "status": status,
        })
    return rows


def _annotation_from_info(info_map: Dict[str, str]) -> Dict[str, str]:
    """Pull gene/consequence/HGVS out of plain INFO keys or an ANN/CSQ field."""
    out = {
        "gene": info_map.get("GENE", info_map.get("SYMBOL", "")),
        "exon": info_map.get("EXON", ""),
        "consequence": info_map.get("CONSEQUENCE", ""),
        "hgvs_c": info_map.get("HGVS_C", info_map.get("HGVSC", "")),
        "hgvs_p": info_map.get("HGVS_P", info_map.get("HGVSP", "")),
    }
    ann = info_map.get("ANN", "")
    csq = info_map.get("CSQ", "")
    raw = ann or csq
    if not raw:
        return out
    fields = raw.split(",")[0].split("|")

    def take(index, key):
        if len(fields) > index and not out[key] and fields[index]:
            out[key] = fields[index]

    take(1, "consequence")
    take(3, "gene")
    if ann:
        take(8, "exon")
        take(9, "hgvs_c")
        take(10, "hgvs_p")
    return out


def _split_field(value: str, index: int) -> str:
    """Return element ``index`` of a comma-separated FORMAT value, or the value."""
    if "," not in value:
        return value
    parts = value.split(",")
    return parts[index] if 0 <= index < len(parts) else ""


def parse_amplicon_vcf(vcf_path: str) -> List[Dict[str, Any]]:
    """Parse variants from a Rindels VCF into one row per ALT allele."""
    if not vcf_path or not os.path.isfile(vcf_path):
        return []

    variants: List[Dict[str, Any]] = []
    try:
        with _open_text(vcf_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue

                chrom, pos, _var_id, ref, alt, _qual, filt, info = fields[:8]
                fmt = fields[8] if len(fields) > 8 else ""
                sample_data = fields[9] if len(fields) > 9 else ""

                fmt_map: Dict[str, str] = {}
                if fmt and sample_data:
                    fmt_map = dict(zip(fmt.split(":"), sample_data.split(":")))

                info_map: Dict[str, str] = {}
                for item in info.split(";"):
                    if "=" in item:
                        k, v = item.split("=", 1)
                        info_map[k] = v
                    elif item:
                        info_map[item] = "True"

                ann = _annotation_from_info(info_map)

                dp = _to_int(fmt_map.get("DP"), 0) or _to_int(info_map.get("DP"), 0)
                ad_values = [p for p in fmt_map.get("AD", "").split(",") if p != ""]

                # One row per ALT allele so a multi-allelic record never shows
                # "A,T" against the depth of a single allele.
                for alt_index, alt_allele in enumerate(alt.split(",")):
                    if not alt_allele or alt_allele == ".":
                        continue

                    ad_alt = _to_int(ad_values[alt_index + 1], 0) if len(ad_values) > alt_index + 1 \
                        else (_to_int(ad_values[0], 0) if len(ad_values) == 1 else 0)

                    af_raw = _split_field(fmt_map.get("AF", fmt_map.get("VAF", "")), alt_index)
                    if af_raw in ("", ".", None):
                        af_raw = _split_field(info_map.get("AF", ""), alt_index)

                    if af_raw not in ("", ".", None):
                        vaf = _to_float(af_raw)
                        if vaf > 1.0:
                            vaf /= 100.0
                    elif dp > 0:
                        # No AF/VAF anywhere: fall back to AD/DP rather than
                        # silently reporting a real variant at 0 %.
                        vaf = ad_alt / dp
                    else:
                        vaf = 0.0

                    svlen = _to_int(_split_field(info_map.get("SVLEN", ""), alt_index), 0)
                    if alt_allele.startswith(SYMBOLIC_ALT_PREFIX):
                        length_change = svlen
                    else:
                        length_change = len(alt_allele) - len(ref)

                    variants.append({
                        "chrom": chrom,
                        "pos": _to_int(pos),
                        "ref": ref,
                        "alt": alt_allele,
                        "filter": filt,
                        "gene": ann["gene"],
                        "exon": ann["exon"],
                        "consequence": ann["consequence"],
                        "hgvs_c": ann["hgvs_c"],
                        "hgvs_p": ann["hgvs_p"],
                        "depth": dp,
                        "alt_depth": ad_alt,
                        "vaf": round(vaf * 100.0, 2),
                        "length_change": length_change,
                    })
    except Exception as e:
        sys.stderr.write(f"Warning: failed parsing VCF {vcf_path}: {e}\n")

    return variants


def is_pass_filter(value: Any) -> bool:
    """Whether a VCF FILTER value counts as a passing call."""
    return str(value or "PASS").strip().upper() in ("PASS", ".", "")


def parse_amplicon_pileup_counts(path: str) -> Dict[Tuple[str, int], Dict[str, Any]]:
    """Parse the Rindels pileup counts TSV for hotspot noise auditing.

    Header keys are lower-cased so callers are insensitive to the caller's
    column casing.
    """
    if not path or not os.path.isfile(path):
        return {}
    pileup_map: Dict[Tuple[str, int], Dict[str, Any]] = {}
    try:
        with _open_text(path) as fh:
            header = None
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if header is None:
                    header = [p.strip().lower() for p in parts]
                    continue
                row = dict(zip(header, parts))
                chrom = row.get("chrom", "")
                if not chrom:
                    continue
                try:
                    pos = int(row["pos"])
                except (ValueError, TypeError, KeyError):
                    continue
                pileup_map[(chrom, pos)] = row
    except Exception as e:
        sys.stderr.write(f"Warning: failed parsing pileup {path}: {e}\n")
    return pileup_map


def _hotspot_span(hotspot: Dict[str, Any]) -> Tuple[int, int]:
    span = hotspot.get("span")
    if span and len(span) == 2:
        return int(span[0]), int(span[1])
    pos = int(hotspot["pos"])
    return pos - 1, pos + 1


def variant_matches_hotspot(variant: Dict[str, Any], hotspot: Dict[str, Any]) -> bool:
    """Whether a called variant belongs to this hotspot codon.

    Matching is by overlap with the hotspot codon span rather than a fixed
    distance, and an indel is widened to the left by its own length (capped at
    12 bp) so that a left-aligned representation still overlaps the codon.
    When the VCF carries a gene annotation it must agree, so a neighbouring
    variant in another gene is never relabelled.
    """
    if str(variant.get("chrom", "")) != hotspot["chrom"]:
        return False

    variant_gene = str(variant.get("gene", "") or "").strip().upper()
    if variant_gene and variant_gene != hotspot["gene"].upper():
        return False

    start, end = _hotspot_span(hotspot)
    pos = _to_int(variant.get("pos"), 0)
    if pos <= 0:
        return False

    ref_len = max(len(str(variant.get("ref", "") or "")), 1)
    indel_len = min(abs(_to_int(variant.get("length_change"), 0)), 12)
    var_start = pos - indel_len
    var_end = pos + ref_len - 1
    return var_start <= end and var_end >= start


def _find_called_variant(variant_rows, hotspot):
    matches = [v for v in variant_rows if variant_matches_hotspot(v, hotspot)]
    if not matches:
        return None
    # Prefer a passing call, then the highest allele fraction.
    matches.sort(key=lambda v: (is_pass_filter(v.get("filter")), _to_float(v.get("vaf"))), reverse=True)
    return matches[0]


def _pileup_row_for_hotspot(pileup_map, hotspot):
    """Pileup row at the hotspot base, else anywhere in its codon span."""
    chrom = hotspot["chrom"]
    row = pileup_map.get((chrom, int(hotspot["pos"])))
    if row is not None:
        return row, int(hotspot["pos"])
    start, end = _hotspot_span(hotspot)
    for pos in range(start, end + 1):
        row = pileup_map.get((chrom, pos))
        if row is not None:
            return row, pos
    return None, int(hotspot["pos"])


def _describe_call(called_var: Dict[str, Any]) -> str:
    """Human label for a detected hotspot call, including its own coordinate."""
    label = called_var.get("hgvs_p") or called_var.get("hgvs_c") \
        or f"{called_var.get('ref')}>{called_var.get('alt')}"
    vaf_call = called_var.get("vaf", 0.0)
    pos = called_var.get("pos", "")
    suffix = "" if is_pass_filter(called_var.get("filter")) else f", {called_var.get('filter')}"
    return f"MUTATION DETECTED ({label} @ {pos}, {vaf_call}%{suffix})"


def audit_hotspots(
    pileup_map: Dict[Tuple[str, int], Dict[str, Any]],
    variant_rows: List[Dict[str, Any]],
) -> List[Dict[str, Any]]:
    """Audit background noise floor across clinical AML hotspots."""
    results = []
    for hs in CLINICAL_HOTSPOTS:
        called_var = _find_called_variant(variant_rows, hs)
        row_data, read_pos = _pileup_row_for_hotspot(pileup_map, hs)

        if row_data is None:
            results.append({
                "gene": hs["gene"], "codon": hs["codon"], "exon": hs["exon"],
                "chrom": hs["chrom"], "pos": hs["pos"], "ref": hs["ref"],
                "depth": 0, "non_ref_reads": 0, "noise_pct": 0.0,
                "clinical": hs["clinical"],
                "status": _describe_call(called_var) if called_var else "NOT ASSESSED",
            })
            continue

        dp = _to_int(row_data.get("filtered_depth"), 0) or _to_int(row_data.get("raw_depth"), 0)
        base_ref = str(row_data.get("ref", hs["ref"]) or hs["ref"]).upper()
        other_sum = sum(
            _to_int(row_data.get(b.lower()), 0) for b in ("A", "C", "G", "T") if b != base_ref
        )
        other_sum += _to_int(row_data.get("ins"), 0) + _to_int(row_data.get("del"), 0)
        noise_pct = (other_sum / dp * 100.0) if dp > 0 else 0.0

        if called_var:
            status = _describe_call(called_var)
        elif dp < AMPLICON_PASS_READS:
            status = f"LOW DEPTH ({dp}x)"
        elif noise_pct > 1.5:
            status = f"ELEVATED NOISE ({noise_pct:.2f}%)"
        else:
            status = f"CLEAN (Noise {noise_pct:.2f}%)"

        results.append({
            "gene": hs["gene"], "codon": hs["codon"], "exon": hs["exon"],
            "chrom": hs["chrom"], "pos": read_pos, "ref": base_ref,
            "depth": dp, "non_ref_reads": other_sum, "noise_pct": round(noise_pct, 2),
            "clinical": hs["clinical"], "status": status,
        })
    return results


# ---------------------------------------------------------------------------
# FLT3-ITD
# ---------------------------------------------------------------------------

def parse_region(region: Any, default=FLT3_ITD_REGION) -> Tuple[str, int, int]:
    """Parse 'chr13:28033864-28034426' into (chrom, start, end)."""
    text = str(region or "").strip()
    if not text or ":" not in text:
        return default
    chrom, _, coords = text.partition(":")
    start_s, _, end_s = coords.replace(",", "").partition("-")
    start, end = _to_int(start_s, -1), _to_int(end_s, -1)
    if not chrom or start < 0 or end < start:
        return default
    return chrom, start, end


def _itd_row(
    gene="FLT3", chrom="", start="", end="", itd_len=0, af_pct=0.0, ar: Any = 0.0,
    status=STATUS_NEGATIVE, filt="", depth=0, alt_depth=0, inserted_seq="", region_str="",
) -> Dict[str, Any]:
    detected = status == STATUS_DETECTED
    ar_num = _to_float(ar)
    if detected:
        frame = "In-frame" if itd_len % 3 == 0 else "Frameshift"
        ar_cat = "AR >= 0.5" if (ar_num >= 0.5 or str(ar).startswith(">")) else "AR < 0.5"
    else:
        frame = "-"
        ar_cat = "Wildtype (Normal)" if status == STATUS_NEGATIVE else "-"
    if detected:
        if isinstance(ar, str) and ar.startswith(">"):
            ar_display = ar
        elif ar_num >= 999.0:
            ar_display = ">100"
        else:
            ar_display = round(ar_num, 2)
    else:
        ar_display = "-"
    return {
        "gene": gene,
        "chrom": chrom,
        "dup_start": str(start),
        "dup_end": str(end),
        "region": region_str or (f"{chrom}:{start}-{end}" if chrom and start != "" else ""),
        "itd_length": str(itd_len),
        "dup_len": str(itd_len),
        "af": f"{af_pct / 100.0:.6f}" if af_pct else "0.0",
        "af_pct": round(af_pct, 2),
        "ar": ar_display,
        "in_frame": "true" if (detected and itd_len % 3 == 0) else ("false" if detected else ""),
        "reading_frame": frame,
        # ELN 2017 used an allelic-ratio cut-off; ELN 2022 dropped it from risk
        # stratification. The category is reported, but not as an ELN risk group.
        "eln_risk": ar_cat,
        "ar_category": ar_cat,
        "status": status,
        "filter": filt,
        "depth": str(depth),
        "alt_depth": str(alt_depth),
        "inserted_seq": inserted_seq,
        "confidence": "",
        "alt_unique_reads": str(alt_depth),
        "ref_unique_reads": str(max(depth - alt_depth, 0)),
    }


def derive_flt3_itd_from_variants(
    variant_rows: Sequence[Dict[str, Any]],
    region: Any = FLT3_ITD_REGION,
    min_length: int = FLT3_ITD_MIN_LENGTH,
    min_vaf_pct: float = FLT3_ITD_MIN_VAF_PCT,
    vcf_available: bool = True,
) -> List[Dict[str, Any]]:
    """Derive FLT3-ITD calls from insertions in the Rindels long-read VCF.

    ``rindels itd`` is a short-read caller and is not used on ONT data; on long
    reads the duplication is carried by the alignment itself and surfaces as an
    insertion inside the FLT3 exon 13-15 target.

    Returns a NOT_ASSESSED row when no VCF was available, so that a failed or
    missing call is never delivered as a negative result.
    """
    chrom, start, end = parse_region(region)
    region_str = f"{chrom}:{start}-{end}"

    if not vcf_available:
        return [_itd_row(status=STATUS_NOT_ASSESSED, chrom=chrom, region_str=region_str)]

    candidates = []
    for v in variant_rows:
        if str(v.get("chrom", "")) != chrom:
            continue
        pos = _to_int(v.get("pos"), 0)
        if not (start <= pos <= end):
            continue
        ins_len = _to_int(v.get("length_change"), 0)
        if ins_len < min_length:
            continue
        if not is_pass_filter(v.get("filter")):
            continue
        vaf = _to_float(v.get("vaf"), 0.0)
        if vaf < min_vaf_pct:
            continue
        candidates.append((vaf, ins_len, pos, v))

    if not candidates:
        return [_itd_row(status=STATUS_NEGATIVE, chrom=chrom, region_str=region_str)]

    # Report the dominant clone first.
    candidates.sort(key=lambda c: (c[0], c[1]), reverse=True)

    rows = []
    for vaf, ins_len, pos, v in candidates:
        depth = _to_int(v.get("depth"), 0)
        alt_depth = _to_int(v.get("alt_depth"), 0)
        wt = depth - alt_depth
        if wt > 0 and alt_depth > 0:
            ar = alt_depth / wt
        elif alt_depth > 0 or vaf >= 100.0:
            # Homozygous ITD / LOH: no wild-type allele detected
            ar = 999.0
        elif 0.0 < vaf < 100.0:
            ar = vaf / (100.0 - vaf)
        else:
            ar = 0.0
        alt_seq = str(v.get("alt", "") or "")
        inserted = alt_seq[len(str(v.get("ref", "") or "")):] if not alt_seq.startswith(SYMBOLIC_ALT_PREFIX) else ""
        rows.append(_itd_row(
            gene="FLT3", chrom=chrom, start=pos, end=pos + ins_len,
            itd_len=ins_len, af_pct=vaf, ar=ar, status=STATUS_DETECTED,
            filt=str(v.get("filter", "")), depth=depth, alt_depth=alt_depth,
            inserted_seq=inserted, region_str=region_str,
        ))
    return rows


ITD_TSV_COLUMNS = [
    "gene", "chrom", "dup_start", "dup_end", "dup_len", "af", "af_pct", "ar",
    "in_frame", "reading_frame", "ar_category", "status", "filter", "depth",
    "alt_depth", "inserted_seq", "region",
]


def write_flt3_itd_tsv(rows: List[Dict[str, Any]], out_path: str) -> None:
    """Write the derived FLT3-ITD table."""
    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(out_path, "w") as fh:
        fh.write("\t".join(ITD_TSV_COLUMNS) + "\n")
        for r in rows:
            fh.write("\t".join(str(r.get(c, "")) for c in ITD_TSV_COLUMNS) + "\n")


def parse_flt3_itd_tsv(itd_path: str) -> List[Dict[str, Any]]:
    """Read back a FLT3-ITD TSV written by this module.

    Tolerates columns being absent or non-numeric: a malformed field degrades
    that one value rather than dropping the row or raising.
    """
    if not itd_path or not os.path.isfile(itd_path):
        return []

    rows: List[Dict[str, Any]] = []
    try:
        with _open_text(itd_path) as fh:
            header = None
            for line in fh:
                line = line.rstrip("\n").strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if header is None:
                    header = [p.strip().lower() for p in parts]
                    continue
                row_dict = dict(zip(header, parts))

                chrom = row_dict.get("chrom", "")
                start = row_dict.get("dup_start", "")
                end = row_dict.get("dup_end", "")
                coords = f"{chrom}:{start}-{end}" if chrom and start and end else \
                    row_dict.get("region", row_dict.get("interval", ""))

                itd_len = _to_int(
                    row_dict.get("dup_len") or row_dict.get("itd_length") or row_dict.get("length"), 0
                )
                af_pct = _as_percent(
                    row_dict.get("af_pct") if row_dict.get("af_pct") not in (None, "", ".")
                    else row_dict.get("af", row_dict.get("junction_af", 0.0))
                )

                raw_status = str(row_dict.get("status", "")).strip().upper()
                if raw_status in (STATUS_DETECTED, STATUS_NEGATIVE, STATUS_NOT_ASSESSED):
                    status_str = raw_status
                elif raw_status:
                    # An explicit status from the caller is respected, never
                    # overridden by an allele fraction.
                    status_str = STATUS_NEGATIVE
                else:
                    status_str = STATUS_DETECTED if af_pct > FLT3_ITD_MIN_VAF_PCT else STATUS_NEGATIVE
                is_det = status_str == STATUS_DETECTED

                alt_u = _to_int(row_dict.get("alt_unique_reads", row_dict.get("alt_depth")), 0)
                ref_u = _to_int(row_dict.get("ref_unique_reads"), 0)
                if ref_u <= 0:
                    ref_u = max(_to_int(row_dict.get("depth"), 0) - alt_u, 0)
                if row_dict.get("ar") not in (None, "", ".", "-"):
                    ar_str = str(row_dict.get("ar")).strip()
                    if ar_str.startswith(">"):
                        ar_val = 999.0
                    else:
                        ar_val = _to_float(ar_str, 0.0)
                elif ref_u > 0 and alt_u > 0:
                    ar_val = alt_u / ref_u
                elif is_det and (alt_u > 0 or af_pct >= 100.0):
                    ar_val = 999.0
                elif is_det and 0.0 < af_pct < 100.0:
                    ar_val = af_pct / (100.0 - af_pct)
                else:
                    ar_val = 0.0

                in_frame_raw = str(row_dict.get("in_frame", "")).strip().lower()
                if in_frame_raw in ("true", "1", "yes"):
                    frame_str = "In-frame"
                elif in_frame_raw in ("false", "0", "no"):
                    frame_str = "Frameshift"
                elif row_dict.get("reading_frame"):
                    frame_str = row_dict["reading_frame"]
                else:
                    frame_str = "In-frame" if (itd_len % 3 == 0 and itd_len > 0) \
                        else ("Frameshift" if itd_len > 0 else "-")

                if is_det:
                    ar_cat = "AR >= 0.5" if (ar_val >= 0.5 or str(row_dict.get("ar", "")).startswith(">")) else "AR < 0.5"
                elif status_str == STATUS_NEGATIVE:
                    ar_cat = "Wildtype (Normal)"
                else:
                    ar_cat = "-"

                if is_det:
                    if str(row_dict.get("ar", "")).startswith(">") or ar_val >= 999.0:
                        ar_out = ">100"
                    else:
                        ar_out = round(ar_val, 2)
                else:
                    ar_out = "-"

                rows.append({
                    "gene": row_dict.get("gene", "FLT3"),
                    "region": coords,
                    "itd_length": str(itd_len),
                    "dup_len": str(itd_len),
                    "af": row_dict.get("af", "0.0"),
                    "af_pct": round(af_pct, 2),
                    "ar": ar_out,
                    "reading_frame": frame_str,
                    "eln_risk": ar_cat,
                    "ar_category": ar_cat,
                    "status": status_str,
                    "filter": row_dict.get("filter", ""),
                    "depth": row_dict.get("depth", ""),
                    "alt_depth": row_dict.get("alt_depth", ""),
                    "inserted_seq": row_dict.get("inserted_seq", row_dict.get("seq", "")),
                    "in_frame": row_dict.get("in_frame", ""),
                    "confidence": row_dict.get("confidence", ""),
                    "alt_unique_reads": str(alt_u),
                    "ref_unique_reads": str(ref_u),
                })
    except Exception as e:
        sys.stderr.write(f"Warning: failed parsing ITD TSV {itd_path}: {e}\n")

    return rows


def itd_is_detected(row: Dict[str, Any]) -> bool:
    """Whether an ITD row is a positive call. Never infers from AF alone."""
    return str(row.get("status", "")).strip().upper() == STATUS_DETECTED


def itd_is_assessed(row: Dict[str, Any]) -> bool:
    """Whether the ITD call was actually made (vs. caller failure/missing input)."""
    return str(row.get("status", "")).strip().upper() != STATUS_NOT_ASSESSED


# ---------------------------------------------------------------------------
# Writers
# ---------------------------------------------------------------------------

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
        fh.write("gene\tstatus\tlength\taf_pct\tar\treading_frame\tar_category\tregion\n")
        if itd_rows:
            for itd in itd_rows:
                fh.write(
                    f"{itd.get('gene', 'FLT3')}\t{itd.get('status', STATUS_NEGATIVE)}\t"
                    f"{itd.get('itd_length', 0)}\t{itd.get('af_pct', 0)}\t"
                    f"{itd.get('ar', '-')}\t{itd.get('reading_frame', '-')}\t"
                    f"{itd.get('ar_category', itd.get('eln_risk', '-'))}\t{itd.get('region', '')}\n"
                )
        else:
            fh.write(
                f"FLT3\t{STATUS_NOT_ASSESSED}\t0\t0.0\t-\t-\t-\t"
                f"{FLT3_ITD_REGION[0]}:{FLT3_ITD_REGION[1]}-{FLT3_ITD_REGION[2]}\n"
            )

        # PASS calls first; filtered calls are listed but clearly separated.
        passing = [v for v in variant_rows if is_pass_filter(v.get("filter"))]
        filtered = [v for v in variant_rows if not is_pass_filter(v.get("filter"))]

        var_headers = ["chrom", "pos", "gene", "exon", "consequence", "ref", "alt", "vaf", "depth", "alt_depth", "filter"]
        fh.write("\n# Section: VARIANTS\n")
        fh.write("\t".join(var_headers) + "\n")
        for v in passing:
            fh.write("\t".join(str(v.get(h, "")) for h in var_headers) + "\n")

        fh.write("\n# Section: FILTERED_VARIANTS\n")
        fh.write("# Not for clinical reporting: calls the variant caller filtered out.\n")
        fh.write("\t".join(var_headers) + "\n")
        for v in filtered:
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
    parser.add_argument("--pileup", type=str, default="", help="Input Rindels pileup counts TSV")
    parser.add_argument("--pools", type=str, help="Pools definition YAML")
    parser.add_argument("--sample", type=str, default="", help="Sample name")
    parser.add_argument(
        "--flt3-itd-region", type=str, default="",
        help="Region searched for FLT3-ITD insertions (default: FLT3 exon 13-15 target)",
    )
    parser.add_argument(
        "--flt3-itd-min-length", type=int, default=FLT3_ITD_MIN_LENGTH,
        help="Minimum inserted length counted as an ITD",
    )
    parser.add_argument(
        "--flt3-itd-min-vaf", type=float, default=FLT3_ITD_MIN_VAF_PCT,
        help="Minimum allele fraction (percent) counted as an ITD",
    )
    parser.add_argument("--output-stats-tsv", type=str, help="Output amplicon stats TSV")
    parser.add_argument("--output-itd-tsv", type=str, help="Output derived FLT3-ITD TSV")
    parser.add_argument("--output-summary-tsv", type=str, help="Output clinical summary TSV")
    parser.add_argument("--output-summary-json", type=str, help="Output clinical summary JSON")

    args = parser.parse_args()

    stats_rows = parse_amplicon_stats_json(args.json, args.pools) if args.json else []
    vcf_available = bool(args.vcf and os.path.isfile(args.vcf))
    variant_rows = parse_amplicon_vcf(args.vcf) if args.vcf else []
    pileup_map = parse_amplicon_pileup_counts(args.pileup) if args.pileup else {}
    hotspot_rows = audit_hotspots(pileup_map, variant_rows) if pileup_map else []

    itd_rows = derive_flt3_itd_from_variants(
        variant_rows,
        region=args.flt3_itd_region or FLT3_ITD_REGION,
        min_length=args.flt3_itd_min_length,
        min_vaf_pct=args.flt3_itd_min_vaf,
        vcf_available=vcf_available,
    )

    if args.output_stats_tsv:
        write_amplicon_stats_tsv(stats_rows, args.output_stats_tsv)

    if args.output_itd_tsv:
        write_flt3_itd_tsv(itd_rows, args.output_itd_tsv)

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
