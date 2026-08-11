"""Create compact TSV/JSON summaries from a Sniffles2-style VCF.

This is intentionally a small, dependency-light report step.  It recognises
symbolic SVs and breakend notation, and optionally adds overlapping BED column
4 values (usually gene names).  Full clinical SV interpretation remains the
responsibility of a reference- and panel-specific annotator.
"""

import collections
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
                svtype = "BND" if "[" in alt or "]" in alt else "OTHER"
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
            genes = overlapping_labels(bed_intervals, chrom, max(position - 1, 0), max(end, position))
            is_translocation = svtype in {"BND", "TRA"} or (chrom2 and chrom2 != chrom)
            record = {
                "sample": str(snakemake.wildcards.sample),
                "type": str(snakemake.wildcards.type),
                "svtype": svtype,
                "chrom": chrom,
                "pos": position,
                "end": end,
                "chrom2": chrom2,
                "svlen": svlen,
                "qual": qual,
                "filter": filt,
                "genes": ";".join(genes) if genes else ".",
                "translocation": bool(is_translocation),
            }
            records.append(record)
            counts[svtype] += 1

    columns = [
        "sample",
        "type",
        "svtype",
        "chrom",
        "pos",
        "end",
        "chrom2",
        "svlen",
        "qual",
        "filter",
        "genes",
        "translocation",
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
        "input_vcf": vcf_path,
        "total_records": len(records),
        "svtype_counts": dict(sorted(counts.items())),
        "translocation_records": len(translocations),
        "annotation_bed": annotation_bed or None,
    }
    with open(snakemake.output.json, "w") as json_file:
        json.dump(summary, json_file, indent=2, sort_keys=True)


main()
