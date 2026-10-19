"""Shared, conservative interpretation helpers for BCR-ABL1 report views."""

import re


PHASE_LABELS = {
    "CIS (COMPOUND)": "CIS (Compound)",
    "TRANS (POLYCLONAL)": "TRANS (Polyclonal)",
    "COEXISTING (COMPOUND + POLYCLONAL)": "COEXISTING (Compound + Polyclonal)",
    "LOW_SPANNING_DEPTH": "LOW_SPANNING_DEPTH",
    "INSUFFICIENT_MUTANT_READS": "INSUFFICIENT_MUTANT_READS",
    "AMBIGUOUS": "AMBIGUOUS",
    "UNPHASEABLE_SAME_POSITION": "UNPHASEABLE_SAME_POSITION",
    "UNPHASEABLE_DIFFERENT_CONTIG": "UNPHASEABLE_DIFFERENT_CONTIG",
    "UNPHASEABLE_INDEL": "UNPHASEABLE_INDEL",
    "UNPHASEABLE_NON_SNV_ALLELE": "UNPHASEABLE_NON_SNV_ALLELE",
}


def phase_class(value):
    """Return a phase class only for known exact caller states."""
    state = re.sub(r"\s+", " ", str(value or "").strip().upper())
    return {
        "CIS (COMPOUND)": "cis",
        "TRANS (POLYCLONAL)": "trans",
        "COEXISTING (COMPOUND + POLYCLONAL)": "coexisting",
        "LOW_SPANNING_DEPTH": "low_depth",
        "INSUFFICIENT_MUTANT_READS": "insufficient",
        "AMBIGUOUS": "ambiguous",
        "UNPHASEABLE_SAME_POSITION": "unphaseable",
        "UNPHASEABLE_DIFFERENT_CONTIG": "unphaseable",
        "UNPHASEABLE_INDEL": "unphaseable",
        "UNPHASEABLE_NON_SNV_ALLELE": "unphaseable",
    }.get(state, "unknown")


def summarize_phasing(rows):
    """Summarize every pair by its exact status; never infer TRANS by default."""
    groups = {key: [] for key in ("cis", "trans", "coexisting", "low_depth", "insufficient", "ambiguous", "unphaseable", "unknown")}
    for row in rows or []:
        cls = phase_class(row.get("phase_state"))
        pair = f"{row.get('mutation_1', '?')}+{row.get('mutation_2', '?')}"
        groups[cls].append(pair)
    labels = []
    for cls, title in (("cis", "CIS"), ("trans", "TRANS"), ("coexisting", "COEXISTING"),
                       ("low_depth", "LOW_SPANNING_DEPTH"), ("insufficient", "INSUFFICIENT_MUTANT_READS"),
                       ("ambiguous", "AMBIGUOUS"), ("unphaseable", "UNPHASEABLE"), ("unknown", "UNKNOWN_PHASE_STATE")):
        if groups[cls]:
            labels.append(f"{title} ({', '.join(groups[cls])})")
    return "; ".join(labels) if labels else "NOT_ASSESSED"


def is_resistance_associated(row):
    """Conservatively exclude synonymous and explicitly non-resistance rows."""
    resistance_class = str(row.get("resistance_class", "")).strip().upper().replace("-", "_")
    consequence = str(row.get("consequence", row.get("variant_consequence", ""))).strip().lower()
    if resistance_class:
        return resistance_class in {"RESISTANCE", "RESISTANCE_ASSOCIATED", "RESISTANCE_MUTATION", "CURATED_RESISTANCE"}
    if resistance_class in {"SYNONYMOUS", "UNCLASSIFIED", "OTHER"}:
        return False
    if consequence:
        if consequence in {"synonymous", "synonymous_variant", "silent", "non_resistance", "other"}:
            return False
        if consequence in {"missense", "missense_variant", "protein_altering", "resistance_associated"}:
            return True
    # Backward-compatible fallback for existing panel files: equal amino-acid
    # reference/alternate letters identify synonymous rows such as L387L.
    aa = str(row.get("mutation", "") or row.get("aa_change", "")).strip()
    match = re.match(r"(?:p\.)?([A-Z*])\d+([A-Z*=])$", aa, re.IGNORECASE)
    return bool(match and match.group(1).upper() != match.group(2).upper())


def active_resistance_mutations(rows):
    return [r for r in rows or [] if str(r.get("status", "")).upper() in {"DETECTED", "SUBCLONAL"}
            and is_resistance_associated(r)]


def panel_assessment(all_sites_rows):
    """Assess only when nonempty all-sites output has every row at callable depth."""
    rows = all_sites_rows or []
    if not rows:
        return {"state": "NOT_ASSESSED", "positions": 0, "callable_positions": 0}
    by_pos = {}
    for row in rows:
        pos = str(row.get("pos", row.get("position", row.get("site", "")))).strip()
        chrom = str(row.get("chrom", "")).strip()
        key = f"{chrom}:{pos}" if pos else str(row.get("mutation", row.get("aa_change", "")))
        by_pos.setdefault(key, []).append(row)
    callable_count = 0
    for items in by_pos.values():
        passed = True
        for row in items:
            status = str(row.get("status", row.get("filter_status", ""))).strip().upper()
            reason = str(row.get("filter_reason", "")).strip().upper()
            depth = row.get("depth", row.get("dp", ""))
            try:
                depth_ok = float(depth) > 0
            except (TypeError, ValueError):
                depth_ok = False
            if status not in {"PASS", "ASSESSED", "DETECTED", "SUBCLONAL", "REFERENCE", "WILDTYPE"} or "LOW_DEPTH" in reason or status == "LOW_DEPTH" or not depth_ok:
                passed = False
                break
        if passed:
            callable_count += 1
    try:
        expected = int(float(rows[0].get("panel_sites_expected", "")))
    except (TypeError, ValueError):
        expected = 0
    complete = expected > 0 and expected == len(by_pos) and callable_count == len(by_pos)
    state = "ASSESSED" if complete else ("INCOMPLETE" if callable_count < len(by_pos) else "COMPLETENESS_UNVERIFIED")
    return {"state": state, "positions": len(by_pos), "callable_positions": callable_count, "expected_positions": expected}

