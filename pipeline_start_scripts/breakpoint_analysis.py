#!/usr/bin/env python3
"""
Analysera SA-tag split reads mellan chr9 och chr20 over flera samples.

Tar BAM/CRAM-filer direkt (t.ex. nagini:s target-BAM:ar) eller redan
extraherade SAM-filer, och:

  1. Parsar varje alignment med en SA-tagg som lankar chr9 <-> chr20.
  2. Normaliserar varje brytpunkt till (pos_chr9, pos_chr20), oavsett om
     den primara alignmenten lag pa chr9 eller chr20 - dvs for varje prov
     summeras stod fran bada riktningarna automatiskt for samma brytpunkt.
  3. Klustrar brytpunkter som ligger inom ett fonster (--window, default 5 bp)
     av varandra.
  4. Skriver ut TVA tabeller:

     a) PER-SAMPLE-tabell (lang format): en rad per (prov, brytpunkt) med
        stod for det provet, rankad fallande pa stod.
        -> default: <output>_per_sample.tsv

     b) CROSS-SAMPLE-tabell (brett format): en rad per brytpunktskluster
        med en kolumn per prov + totalt stod + antal prov som stodjer
        klustret, sa att aterkommande positioner syns direkt.
        -> default: <output>_recurring.tsv

     Aven kluster/rader med lagt stod (t.ex. 1 read) tas med och flaggas,
     sa att du sjalv kan bedoma svaga kandidater istallet for att de
     filtreras bort i tysthet.

BAM-input las genom `samtools view -h` och begransas till --chrom-a/--chrom-b
med samtools egen regionslagning, vilket kraver ett index. Det racker: en
split read har exakt en icke-supplementary alignment, och den bar SA-taggen
med partnerns koordinat aven om partnern sjalv ligger utanfor regionen. Darfor
fungerar en target-BAM (designen +/-25 kb) lika bra som hela genomet, och gar
avsevart fortare.

Exempel:
    python3 breakpoint_analysis.py \\
        alignment/target/gustaf*_T.target.bam \\
        --window 5 --min-support 2 --min-samples-recurring 2 \\
        -o breakpoints

Om --sample-regex inte matchar en fil anvands filnamnet utan andelser
(.bam/.cram/.sam, .target och _chr9_to_chr20-suffix) som sample-namn.
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
from collections import defaultdict
from contextlib import contextmanager

SA_RE = re.compile(r"SA:Z:(\S+)")
CHROM_A_DEFAULT = "chr9"
CHROM_B_DEFAULT = "chr20"
ALIGNMENT_SUFFIXES = (".bam", ".cram")


@contextmanager
def alignment_lines(path, regions, reference=None):
    """Yield SAM text lines from a BAM/CRAM or a plain SAM file.

    BAM is read through `samtools view -h` rather than a Python BAM parser so
    that this stays a stdlib-only script, and so that the region filter runs in
    samtools against the index instead of over every record here.
    """

    if path.lower().endswith(ALIGNMENT_SUFFIXES):
        if not shutil.which("samtools"):
            raise SystemExit(
                f"{path} is a BAM/CRAM but samtools is not on PATH. "
                "Install samtools, or pass a SAM file instead."
            )
        index_suffixes = (".bai", ".csi") if path.lower().endswith(".bam") else (".crai",)
        if regions and not any(os.path.exists(path + suffix) for suffix in index_suffixes):
            raise SystemExit(
                f"{path} has no index, which the region filter needs. "
                f"Run `samtools index {path}`, or pass --regions ''."
            )
        command = ["samtools", "view", "-h"]
        if reference:
            command += ["-T", reference]
        command += [path, *regions]
        process = subprocess.Popen(
            command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        try:
            yield process.stdout
        finally:
            process.stdout.close()
            stderr = process.stderr.read()
            process.stderr.close()
            if process.wait() != 0:
                raise SystemExit(f"samtools failed on {path}:\n{stderr.strip()}")
    else:
        with open(path) as handle:
            yield handle


def parse_alignment_file(path, chrom_a, chrom_b, regions, reference=None):
    """Yield (pos_a, pos_b, qname) for varje alignment vars SA-tagg lankar
    chrom_a <-> chrom_b. Supplementary alignments (flag 0x800) hoppas over
    for att undvika dubbelrakning - varje fragment ger da max en rad per
    SA-lank i denna fil."""
    seen = set()
    with alignment_lines(path, regions, reference) as lines:
        for line in lines:
            if line.startswith("@") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                continue
            qname, flag_s, rname, pos_s = fields[0], fields[1], fields[2], fields[3]
            try:
                flag = int(flag_s)
                pos = int(pos_s)
            except ValueError:
                continue

            if flag & 0x800:  # supplementary -> hoppa over, se docstring
                continue
            if rname not in (chrom_a, chrom_b):
                continue

            sa_match = SA_RE.search(line)
            if not sa_match:
                continue

            for sa_entry in sa_match.group(1).split(";"):
                if not sa_entry:
                    continue
                sa_fields = sa_entry.split(",")
                if len(sa_fields) < 2:
                    continue
                sa_rname, sa_pos_s = sa_fields[0], sa_fields[1]
                try:
                    sa_pos = int(sa_pos_s)
                except ValueError:
                    continue

                pair = {}
                pair[rname] = pos
                pair[sa_rname] = sa_pos
                if chrom_a in pair and chrom_b in pair and len(pair) == 2:
                    key = (qname, pair[chrom_a], pair[chrom_b])
                    if key in seen:
                        continue
                    seen.add(key)
                    yield pair[chrom_a], pair[chrom_b], qname


def sample_name_for(path, sample_regex):
    if sample_regex:
        m = re.search(sample_regex, path)
        if m:
            return m.group(1) if m.groups() else m.group(0)
    base = path.split("/")[-1]
    base = re.sub(r"\.(bam|cram|sam)$", "", base, flags=re.IGNORECASE)
    # nagini writes alignment/target/<sample>_<type>.target.bam; the type stays
    # in the name, since a run may carry both N and T of the same sample.
    base = re.sub(r"\.target$", "", base)
    base = re.sub(r"(_chr\d+_to_chr\d+)$", "", base)
    return base


class UnionFind:
    def __init__(self, n):
        self.parent = list(range(n))
        self.rank = [0] * n

    def find(self, x):
        while self.parent[x] != x:
            self.parent[x] = self.parent[self.parent[x]]
            x = self.parent[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            ra, rb = rb, ra
        self.parent[rb] = ra
        if self.rank[ra] == self.rank[rb]:
            self.rank[ra] += 1


def cluster_records(records, window):
    """Single-linkage-klustring: tva brytpunkter slas ihop om bade pos_a
    och pos_b ligger inom `window` bp av varandra. Sorterar pa pos_a sa att
    varje record bara behover jamforas mot foljande records inom fonstret
    (undviker O(n^2) och undviker kantfel fran fast-grid-binning dar t.ex.
    999 och 1002 hade hamnat i olika bins trots 3 bp avstand).

    OBS: single-linkage kan i teorin kedja ihop en lang rad av brytpunkter
    som ligger tatt efter varandra (A nara B nara C aven om A och C ligger
    langre isar an window). For ett litet fonster (t.ex. 5 bp) pa verkliga
    breakpoint-data ar detta ovanligt, men var uppmarksam om pos_a_range/
    pos_b_range i outputen blir ovantat brett."""
    n = len(records)
    uf = UnionFind(n)
    order = sorted(range(n), key=lambda i: records[i][0])  # sortera pa pos_a
    for oi, i in enumerate(order):
        pos_a_i, pos_b_i = records[i][0], records[i][1]
        for j in order[oi + 1:]:
            pos_a_j = records[j][0]
            if pos_a_j - pos_a_i > window:
                break  # sorterat pa pos_a -> inga fler kan vara inom fonstret
            if abs(pos_b_i - records[j][1]) <= window:
                uf.union(i, j)

    groups = defaultdict(list)
    for i in range(n):
        groups[uf.find(i)].append(i)
    return list(groups.values())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("alignments", nargs="+", help="BAM/CRAM- eller SAM-filer, alla samples")
    ap.add_argument("--chrom-a", default=CHROM_A_DEFAULT)
    ap.add_argument("--chrom-b", default=CHROM_B_DEFAULT)
    ap.add_argument("--regions", default=None,
                    help="Kommaseparerade regioner att lasa ur BAM (default: --chrom-a och "
                         "--chrom-b). Tom strang laser hela filen utan index.")
    ap.add_argument("--reference", default=None, help="Referens-FASTA, kravs for CRAM")
    ap.add_argument("--window", type=int, default=5, help="Klustringsfonster i bp (default 5)")
    ap.add_argument("--min-support", type=int, default=2,
                    help="Troskel for is_confident-flaggan pa read-stod (default 2)")
    ap.add_argument("--min-samples-recurring", type=int, default=2,
                    help="Antal prov som maste stodja ett kluster for att flaggas recurring (default 2)")
    ap.add_argument("--sample-regex", default=None,
                    help="Regex med en capture-grupp for att extrahera sample-namn ur filnamnet")
    ap.add_argument("-o", "--output", default="breakpoints",
                    help="Output-prefix. Skriver <prefix>_per_sample.tsv och <prefix>_recurring.tsv")
    args = ap.parse_args(argv)

    if args.regions is None:
        regions = [args.chrom_a, args.chrom_b]
    else:
        regions = [r for r in args.regions.split(",") if r.strip()]

    # Samla forst in ALLA enskilda brytpunkts-observationer (en per read),
    # klustra dem sedan tillsammans over samtliga prov/filer.
    records = []  # list of (pos_a, pos_b, sample, qname)
    samples_seen = set()

    for path in args.alignments:
        sample = sample_name_for(path, args.sample_regex)
        samples_seen.add(sample)
        n = 0
        for pos_a, pos_b, qname in parse_alignment_file(
            path, args.chrom_a, args.chrom_b, regions, args.reference
        ):
            records.append((pos_a, pos_b, sample, qname))
            n += 1
        print(f"[info] {path} (sample={sample}): {n} split-read-lankar", file=sys.stderr)

    samples_sorted = sorted(samples_seen)

    groups = cluster_records(records, args.window)

    # cluster index -> sample -> set(qnames)  (dedupar inom prov+kluster)
    clusters = defaultdict(lambda: defaultdict(set))
    cluster_repr = {}
    for gi, idxs in enumerate(groups):
        pos_a_vals = sorted(records[i][0] for i in idxs)
        pos_b_vals = sorted(records[i][1] for i in idxs)
        cluster_repr[gi] = {
            "pos_a_repr": pos_a_vals[len(pos_a_vals) // 2],
            "pos_a_range": f"{pos_a_vals[0]}-{pos_a_vals[-1]}",
            "pos_b_repr": pos_b_vals[len(pos_b_vals) // 2],
            "pos_b_range": f"{pos_b_vals[0]}-{pos_b_vals[-1]}",
        }
        for i in idxs:
            _, _, sample, qname = records[i]
            clusters[gi][sample].add(qname)

    # --- Tabell A: per-sample, langformat, rankad pa stod -----------------
    per_sample_rows = []
    for ck, per_sample in clusters.items():
        repr_ = cluster_repr[ck]
        for sample, qnames in per_sample.items():
            support = len(qnames)
            if support == 0:
                continue
            per_sample_rows.append({
                "sample": sample,
                "chrom_a": args.chrom_a,
                "pos_a_repr": repr_["pos_a_repr"],
                "pos_a_range": repr_["pos_a_range"],
                "chrom_b": args.chrom_b,
                "pos_b_repr": repr_["pos_b_repr"],
                "pos_b_range": repr_["pos_b_range"],
                "support": support,
                "is_confident": support >= args.min_support,
            })
    per_sample_rows.sort(key=lambda r: r["support"], reverse=True)

    per_sample_path = f"{args.output}_per_sample.tsv"
    with open(per_sample_path, "w") as out:
        header = ["rank", "sample", "chrom_a", "pos_a_repr", "pos_a_range",
                  "chrom_b", "pos_b_repr", "pos_b_range", "support", "is_confident"]
        out.write("\t".join(header) + "\n")
        for i, r in enumerate(per_sample_rows, start=1):
            out.write("\t".join([
                str(i), r["sample"], r["chrom_a"], str(r["pos_a_repr"]), r["pos_a_range"],
                r["chrom_b"], str(r["pos_b_repr"]), r["pos_b_range"],
                str(r["support"]), str(r["is_confident"]),
            ]) + "\n")

    # --- Tabell B: over prov (brett format), rankad pa totalt stod --------
    recurring_rows = []
    for ck, per_sample in clusters.items():
        repr_ = cluster_repr[ck]
        counts = {s: len(per_sample.get(s, ())) for s in samples_sorted}
        total = sum(counts.values())
        supporting_samples = [s for s in samples_sorted if counts[s] > 0]
        n_samples_with_support = len(supporting_samples)

        recurring_rows.append({
            "chrom_a": args.chrom_a,
            "pos_a_repr": repr_["pos_a_repr"],
            "pos_a_range": repr_["pos_a_range"],
            "chrom_b": args.chrom_b,
            "pos_b_repr": repr_["pos_b_repr"],
            "pos_b_range": repr_["pos_b_range"],
            "total_support": total,
            "n_samples_with_support": n_samples_with_support,
            "recurring": n_samples_with_support >= args.min_samples_recurring,
            "is_confident": total >= args.min_support,
            "samples": ";".join(f"{s}:{counts[s]}" for s in supporting_samples),
            "per_sample": counts,
        })
    recurring_rows.sort(key=lambda r: r["total_support"], reverse=True)

    recurring_path = f"{args.output}_recurring.tsv"
    with open(recurring_path, "w") as out:
        header = ["rank", "chrom_a", "pos_a_repr", "pos_a_range", "chrom_b", "pos_b_repr", "pos_b_range",
                  "total_support", "n_samples_with_support", "recurring", "is_confident", "samples"] + samples_sorted
        out.write("\t".join(header) + "\n")
        for i, r in enumerate(recurring_rows, start=1):
            line = [
                str(i), r["chrom_a"], str(r["pos_a_repr"]), r["pos_a_range"],
                r["chrom_b"], str(r["pos_b_repr"]), r["pos_b_range"],
                str(r["total_support"]), str(r["n_samples_with_support"]),
                str(r["recurring"]), str(r["is_confident"]), r["samples"],
            ] + [str(r["per_sample"][s]) for s in samples_sorted]
            out.write("\t".join(line) + "\n")

    n_confident = sum(1 for r in recurring_rows if r["is_confident"])
    n_recurring = sum(1 for r in recurring_rows if r["recurring"])
    print(f"[info] {len(recurring_rows)} brytpunktskluster totalt "
          f"({n_confident} med stod >= {args.min_support}, "
          f"{n_recurring} aterkommande i >= {args.min_samples_recurring} prov)", file=sys.stderr)
    print(f"[info] Skrev {per_sample_path} ({len(per_sample_rows)} rader) och "
          f"{recurring_path} ({len(recurring_rows)} rader)", file=sys.stderr)


if __name__ == "__main__":
    main()
