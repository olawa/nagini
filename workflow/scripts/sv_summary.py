"""Create compact TSV/JSON summaries of the structural variants in a VCF.

This is intentionally a small, dependency-light report step.  It recognises
symbolic SVs and breakend notation, and optionally adds overlapping BED column
4 values (usually gene names).  Full clinical SV interpretation remains the
responsibility of a reference- and panel-specific annotator.
"""

import collections
import bisect
import gzip
import json
import os
import re


def open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path, "r")


def parse_info(value):
    result = {}
    for item in value.split(";"):
        if not item:
            continue
        if "=" in item:
            key, val = item.split("=", 1)
            result[key] = val
        else:
            result[item] = True
    return result


def remote_chrom_from_alt(alt):
    match = re.search(r"[\[\]]([^:\[\]]+):\d+", alt or "")
    return match.group(1) if match else ""


def remote_position_from_alt(alt):
    match = re.search(r"[\[\]][^:\[\]]+:(\d+)", alt or "")
    return int(match.group(1)) if match else 0


class CytobandLookup:
    """Lookup UCSC cytobands from a plain or gzip-compressed BED-like file."""

    def __init__(self, path=""):
        self.bands = collections.defaultdict(list)
        if path and os.path.isfile(path):
            self._load(path)

    def _load(self, path):
        with open_text(path) as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 4:
                    continue
                try:
                    start, end = int(fields[1]), int(fields[2])
                except ValueError:
                    continue
                self.bands[fields[0]].append((end, fields[3], start))
        for chrom in self.bands:
            self.bands[chrom].sort()

    def get_band(self, chrom, position):
        if not self.bands or not chrom or not position:
            return ""
        keys = (chrom, chrom.removeprefix("chr"), f"chr{chrom.removeprefix('chr')}")
        entries = next((self.bands[key] for key in keys if key in self.bands), None)
        if not entries:
            return ""
        position = int(position)
        index = bisect.bisect_left(entries, (position, "", 0))
        if index >= len(entries):
            return ""
        end, band, start = entries[index]
        return band if start <= position <= end else ""


def cytogenetic_notation(svtype, chrom, position, end, chrom2, partner_position, lookup):
    """Create a compact ISCN-style notation when cytobands are available."""

    if not lookup.bands:
        return ""
    chrom1 = str(chrom).removeprefix("chr")
    band1 = lookup.get_band(chrom, position)
    if chrom2 and chrom2 != chrom:
        band2 = lookup.get_band(chrom2, partner_position or end)
        prefix = f"t({chrom1};{str(chrom2).removeprefix('chr')})"
        return f"{prefix}({band1};{band2})" if band1 and band2 else prefix
    band2 = lookup.get_band(chrom, end)
    if not band1 or not band2:
        return ""
    if str(svtype).upper() == "DEL":
        return f"del({chrom1})({band1}{band2})"
    return f"{svtype}({chrom1})({band1};{band2})"


def read_bed(path):
    if not path:
        return {}
    intervals = collections.defaultdict(list)
    with open_text(path) as bed_file:
        for line in bed_file:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            try:
                start, end = int(fields[1]), int(fields[2])
            except ValueError:
                continue
            label = fields[3] if len(fields) > 3 else "."
            intervals[fields[0]].append((start, end, label))
    return intervals


def overlapping_labels(intervals, chrom, start, end):
    labels = set()
    for bed_start, bed_end, label in intervals.get(chrom, []):
        if bed_start < end and start < bed_end:
            labels.add(label)
    return sorted(labels)


def main():
    snakemake_input = snakemake.input
    vcf_path = str(snakemake_input.vcf)
    annotation_bed_value = getattr(snakemake_input, "annotation_bed", "")
    annotation_bed = str(annotation_bed_value) if annotation_bed_value else ""
    bed_intervals = read_bed(annotation_bed)
    cytoband_value = getattr(snakemake_input, "cytoband", "")
    cytoband_path = str(cytoband_value) if cytoband_value else ""
    cytobands = CytobandLookup(cytoband_path)
    # Which caller produced this VCF, so two callers' rows can share a table.
    caller = str(getattr(snakemake.params, "caller", "") or "Sniffles2")

    records = []
    counts = collections.Counter()
    with open_text(vcf_path) as vcf_file:
        for line in vcf_file:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, _record_id, ref, alt, qual, filt, info_text = fields[:8]
            info = parse_info(info_text)
            svtype = str(info.get("SVTYPE", ""))
            if not svtype:
                # Sniffles2 emits only structural records, but a caller run in
                # a combined mode (rindels --sv-lite) writes SNVs and indels
                # into the same VCF. Without this, every one of those millions
                # of records became an "OTHER" row in the report.
                if "[" not in alt and "]" not in alt:
                    continue
                svtype = "BND"
            try:
                position = int(pos)
            except ValueError:
                position = 0
            try:
                end = int(info.get("END", position))
            except (TypeError, ValueError):
                end = position
            try:
                svlen = int(str(info.get("SVLEN", "0")).split(",")[0])
            except ValueError:
                svlen = 0
            chrom2 = str(info.get("CHR2", "")) or remote_chrom_from_alt(alt)
            partner_position = remote_position_from_alt(alt)
            genes = overlapping_labels(bed_intervals, chrom, max(position - 1, 0), max(end, position))
            is_translocation = svtype in {"BND", "TRA"} or (chrom2 and chrom2 != chrom)
            band = cytobands.get_band(chrom, position)
            band2 = cytobands.get_band(
                chrom2 or chrom, partner_position or end
            )
            record = {
                "sample": str(snakemake.wildcards.sample),
                "type": str(snakemake.wildcards.type),
                "caller": caller,
                "svtype": svtype,
                "chrom": chrom,
                "pos": position,
                "end": end,
                "chrom2": chrom2,
                # The partner coordinate was computed for the cytoband lookup
                # and the ISCN notation but never written out, so a
                # translocation row named the partner chromosome without
                # saying where on it. `end` equals `pos` for a breakend, so it
                # could not stand in.
                "pos2": partner_position or (end if chrom2 and chrom2 != chrom else 0),
                "svlen": svlen,
                "qual": qual,
                "filter": filt,
                "genes": ";".join(genes) if genes else ".",
                "translocation": bool(is_translocation),
                "Start_Band": band,
                "End_Band": band2,
                "Cytogenetic_Notation": cytogenetic_notation(
                    svtype,
                    chrom,
                    position,
                    end,
                    chrom2,
                    partner_position,
                    cytobands,
                ),
            }
            records.append(record)
            counts[svtype] += 1

    columns = [
        "sample",
        "type",
        "caller",
        "svtype",
        "chrom",
        "pos",
        "end",
        "chrom2",
        "pos2",
        "svlen",
        "qual",
        "filter",
        "genes",
        "translocation",
        "Start_Band",
        "End_Band",
        "Cytogenetic_Notation",
    ]
    os.makedirs(os.path.dirname(str(snakemake.output.tsv)), exist_ok=True)
    with open(snakemake.output.tsv, "w") as tsv_file:
        tsv_file.write("\t".join(columns) + "\n")
        for record in records:
            tsv_file.write("\t".join(str(record[column]) for column in columns) + "\n")

    translocations = [record for record in records if record["translocation"]]
    with open(snakemake.output.translocations, "w") as trans_file:
        trans_file.write("\t".join(columns) + "\n")
        for record in translocations:
            trans_file.write("\t".join(str(record[column]) for column in columns) + "\n")

    summary = {
        "sample": str(snakemake.wildcards.sample),
        "type": str(snakemake.wildcards.type),
        "caller": caller,
        "input_vcf": vcf_path,
        "total_records": len(records),
        "svtype_counts": dict(sorted(counts.items())),
        "translocation_records": len(translocations),
        "annotation_bed": annotation_bed or None,
        "cytoband": cytoband_path or None,
    }
    with open(snakemake.output.json, "w") as json_file:
        json.dump(summary, json_file, indent=2, sort_keys=True)


if "snakemake" in globals():
    main()
