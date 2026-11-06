import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import pysam

HERE = Path(__file__).parent
spec = importlib.util.spec_from_file_location("bcr_abl_branford_caller", HERE / "bcr_abl_branford_caller.py")
caller = importlib.util.module_from_spec(spec)
spec.loader.exec_module(caller)
phaser_spec = importlib.util.spec_from_file_location("bcr_abl_phaser", HERE / "bcr_abl_phaser.py")
phaser = importlib.util.module_from_spec(phaser_spec)
phaser_spec.loader.exec_module(phaser)


def write_bam(path, rows, length=200):
    header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"SN": "chr9", "LN": length}]}
    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for i, row in enumerate(rows):
            aln = pysam.AlignedSegment()
            aln.query_name = row.get("name", f"read{i}")
            aln.query_sequence = row["seq"]
            aln.query_qualities = row.get("quals", [30] * len(row["seq"]))
            aln.reference_id = 0
            aln.reference_start = row.get("start", 0)
            aln.cigar = row.get("cigar", [(0, len(row["seq"]))])
            aln.mapping_quality = row.get("mapq", 60)
            aln.flag = row.get("flag", 0)
            for tag, value in row.get("tags", {}).items():
                aln.set_tag(tag, value)
            out.write(aln)
    pysam.index(str(path))


def run_caller(tmp, bam, bed_text, extra=()):
    bed = tmp / "panel.bed"
    all_sites = tmp / "all.tsv"
    mutations = tmp / "mutations.tsv"
    vcf = tmp / "mutations.vcf"
    bed.write_text(bed_text)
    argv = ["caller", "--bam", str(bam), "--bed", str(bed), "--sample", "S",
            "--output-all", str(all_sites), "--output-mutations", str(mutations),
            "--output-vcf", str(vcf), "--min-depth", "10", "--min-alt-reads", "5",
            "--min-snr", "3", *extra]
    with patch.object(caller.sys, "argv", argv):
        caller.main()
    def read(path):
        lines = [line for line in path.read_text().splitlines() if line and not line.startswith("#")]
        return list(csv.DictReader(lines, delimiter="\t"))
    return read(all_sites), read(mutations), vcf


class FlankingContextCounting(unittest.TestCase):
    """Reads with a dirty flank are evidence for no allele.

    Modelled on K357E (c.1069A>G): on the reverse strand the G of AAG|AAAAA is
    read one base late, so the read shows G>A at -1 and A>G at the site. In
    every sample this sits at 0.2-0.3 % of depth.
    """

    SITE = 100
    LEFT, RIGHT = "GAGTACCTGGAGAAG", "AAAACTTCATCCACA"
    BED = "chr9\t100\t101\tNM_005157.6:c.(1069A>G)\tA\tG\tp.Lys357Glu\n"

    def setUp(self):
        import random
        rng = random.Random(3)
        filler = lambda n: "".join(rng.choice("ACGT") for _ in range(n))  # noqa: E731
        self.ref = filler(self.SITE - len(self.LEFT)) + self.LEFT + "A" + self.RIGHT
        self.ref += filler(200 - len(self.ref))
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.tmp = Path(self.directory.name)
        self.fasta = self.tmp / "ref.fa"
        self.fasta.write_text(">chr9\n" + self.ref + "\n")
        pysam.faidx(str(self.fasta))

    def read(self, i, alt=False, artefact=False, flank_error=False, artefact_qual=18):
        seq = list(self.ref[25:175])
        quals = [30] * 150
        site = self.SITE - 25
        if alt or artefact:
            seq[site] = "G"
        if artefact:
            seq[site - 1] = "A"                     # the G read one base late
            quals[site - 1] = quals[site] = artefact_qual   # and called with less confidence
        if flank_error:
            seq[site + 2] = "C" if seq[site + 2] != "C" else "G"
        return {"name": f"r{i}", "seq": "".join(seq), "start": 25, "quals": quals,
                "flag": 16 if (artefact or i % 3 == 0) else 0}

    def call(self, rows, reference=True):
        bam = self.tmp / "reads.bam"
        write_bam(bam, rows)
        extra = ("--reference-fasta", str(self.fasta)) if reference else ()
        audit, variants, _ = run_caller(self.tmp, bam, self.BED, extra)
        return audit[0]

    def test_a_subclone_is_called_although_the_artefact_sits_on_top_of_it(self):
        """12 clean G (0.30 %) + 10 shifted G. Before, 22 alt reads with 55 % clean
        context failed DIRTY_CONTEXT and the subclone was lost."""
        rows = [self.read(i, alt=True) for i in range(12)]
        rows += [self.read(100 + i, artefact=True) for i in range(10)]
        rows += [self.read(1000 + i) for i in range(3978)]
        row = self.call(rows)
        self.assertEqual(row["status"], "SUBCLONAL")
        self.assertEqual(row["ad_alt"], "12")
        self.assertEqual(row["depth"], "3990")
        self.assertEqual(row["ad_alt_context_excluded"], "10")
        self.assertAlmostEqual(float(row["vaf_percent"]), 12 / 3990 * 100, places=2)
        self.assertAlmostEqual(float(row["kmer_clean_percent"]), 54.5, places=1)   # audit of all 22

    def test_the_artefact_alone_is_not_called(self):
        """The G>A rides on most ALT reads, so it is judged as a linked
        neighbour -- and fails both haplotype checks, as in the dilution
        series: one strand, and base quality 0.6 of the reference's."""
        rows = [self.read(i, alt=True) for i in range(4)]
        rows += [self.read(100 + i, artefact=True) for i in range(10)]
        rows += [self.read(1000 + i) for i in range(3986)]
        row = self.call(rows)
        self.assertNotIn(row["status"], ("DETECTED", "SUBCLONAL", "REVIEW"))
        self.assertEqual(row["ad_alt"], "4")
        self.assertIn("rejected(STRAND(", row["linked_neighbour"])
        self.assertIn(";BASEQ(", row["linked_neighbour"])

    def test_the_artefact_at_full_base_quality_is_for_review(self):
        """One-stranded but called confidently: one check of two fails, and a
        person looks at it rather than the caller deciding."""
        rows = [self.read(i, alt=True) for i in range(4)]
        rows += [self.read(100 + i, artefact=True, artefact_qual=30) for i in range(10)]
        rows += [self.read(1000 + i) for i in range(3986)]
        row = self.call(rows)
        self.assertEqual(row["status"], "REVIEW")
        self.assertTrue(row["filter_reason"].startswith("LINKED_NEIGHBOUR_REVIEW(STRAND("), row["filter_reason"])

    def test_excluding_on_both_sides_keeps_the_vaf_unbiased(self):
        """Every 12th read has a flank error, whatever its allele. Excluding
        them only from the alternate side would read 5 % as about 4.6 %."""
        rows = [self.read(i, alt=i < 200, flank_error=i % 12 == 0) for i in range(4000)]
        row = self.call(rows)
        self.assertAlmostEqual(float(row["vaf_percent"]), 5.0, delta=0.1)
        self.assertGreater(int(row["depth_context_excluded"]), 300)

    def test_without_a_reference_nothing_is_excluded(self):
        rows = [self.read(i, alt=True) for i in range(12)]
        rows += [self.read(100 + i, artefact=True) for i in range(10)]
        rows += [self.read(1000 + i) for i in range(3978)]
        row = self.call(rows, reference=False)
        self.assertEqual(row["ad_alt"], "22")
        self.assertEqual(row["depth"], "4000")
        self.assertEqual(row["depth_context_excluded"], "0")


class CallerRegressions(unittest.TestCase):
    def test_single_site_variant_copy_is_loadable_by_phaser_and_audit_is_preserved(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            rows = []
            for i, base in enumerate("A" * 20 + "G" * 80):
                seq = list("N" * 150); seq[49] = base
                rows.append({"seq": "".join(seq), "name": f"r{i}"})
            bam = tmp / "reads.bam"; write_bam(bam, rows)
            all_rows, variants, _ = run_caller(
                tmp, bam, "chr9\t49\t50\tNM:c.(763G>A)\tG\tA\tp.Glu255Lys\n"
            )
            self.assertEqual(all_rows[0]["record_type"], "SITE_AUDIT")
            self.assertEqual(variants[0]["record_type"], "SITE_VARIANT")
            loaded = phaser.load_detected_mutations(tmp / "mutations.tsv")
            self.assertEqual([(r["mutation"], r["record_type"]) for r in loaded], [("E255K", "SITE_VARIANT")])

    def test_two_supported_panel_alt_alleles_do_not_mask_each_other_as_snr_noise(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            rows = []
            for i, base in enumerate("T" * 60 + "A" * 20 + "C" * 20):
                seq = list("N" * 150); seq[49] = base
                rows.append({"seq": "".join(seq), "name": f"r{i}"})
            bam = tmp / "reads.bam"; write_bam(bam, rows)
            all_rows, variants, _ = run_caller(
                tmp, bam,
                "chr9\t49\t50\tNM:c.(949T>A)\tT\tA\tp.Phe317Ile\n"
                "chr9\t49\t50\tNM:c.(949T>C)\tT\tC\tp.Phe317Leu\n",
            )
            self.assertEqual([r["status"] for r in all_rows], ["DETECTED", "DETECTED"])
            self.assertEqual({r["mutation"] for r in variants}, {"F317I", "F317L"})
            self.assertEqual(all_rows[0]["credible_panel_alt_competitors"], "C")

    def test_same_codon_cis_reads_emit_aggregate_protein_haplotype(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            ref = list("N" * 200); ref[99:102] = list("ACT")
            fasta = tmp / "ref.fa"; fasta.write_text(">chr9\n" + "".join(ref) + "\n"); pysam.faidx(str(fasta))
            rows = []
            for i in range(20):
                seq = list("N" * 150); seq[99:102] = list("GTT"); rows.append({"seq": "".join(seq)})
            for i in range(80):
                seq = list("N" * 150); seq[99:102] = list("ACT"); rows.append({"seq": "".join(seq)})
            bam = tmp / "reads.bam"; write_bam(bam, rows)
            bed = (
                "chr9\t99\t100\tNM_005157.6:c.(943A>G)\tA\tG\tp.Thr315Ala\n"
                "chr9\t100\t101\tNM_005157.6:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
            )
            all_rows, variants, vcf = run_caller(tmp, bam, bed, ("--reference-fasta", str(fasta)))
            self.assertEqual(len(all_rows), 2)  # nucleotide audit is retained
            self.assertTrue(all(r["codon_context_status"] == "RESOLVED" for r in all_rows))
            self.assertEqual(len(variants), 1)
            self.assertEqual(variants[0]["record_type"], "CODON_HAPLOTYPE")
            self.assertEqual((variants[0]["mutation"], variants[0]["ref"], variants[0]["alt"]), ("T315V", "ACT", "GTT"))
            self.assertEqual(variants[0]["ad_alt"], "20")
            vcf_rows = [line.split("\t") for line in vcf.read_text().splitlines() if line and not line.startswith("#")]
            self.assertEqual((vcf_rows[0][3], vcf_rows[0][4]), ("ACT", "GTT"))

    def test_supported_independent_codon_haplotypes_do_not_mask_each_other(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            ref = list("N" * 200); ref[99:102] = list("ACT")
            fasta = tmp / "ref.fa"; fasta.write_text(">chr9\n" + "".join(ref) + "\n"); pysam.faidx(str(fasta))
            rows = []
            for codon, count in (("ACT", 60), ("GCT", 20), ("ATT", 20)):
                for _ in range(count):
                    seq = list("N" * 150); seq[99:102] = list(codon); rows.append({"seq": "".join(seq)})
            bam = tmp / "reads.bam"; write_bam(bam, rows)
            bed = (
                "chr9\t99\t100\tNM_005157.6:c.(943A>G)\tA\tG\tp.Thr315Ala\n"
                "chr9\t100\t101\tNM_005157.6:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
            )
            all_rows, variants, _ = run_caller(tmp, bam, bed, ("--reference-fasta", str(fasta)))
            self.assertEqual([row["status"] for row in all_rows], ["DETECTED", "DETECTED"])
            self.assertEqual({(row["alt"], row["ad_alt"], row["status"]) for row in variants}, {
                ("GCT", "20", "DETECTED"), ("ATT", "20", "DETECTED"),
            })

    def test_three_balanced_panel_alternates_remain_conservative(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            rows = []
            for i, base in enumerate("T" * 970 + "A" * 10 + "C" * 10 + "G" * 10):
                seq = list("N" * 150); seq[49] = base
                rows.append({"seq": "".join(seq), "name": f"r{i}"})
            bam = tmp / "reads.bam"; write_bam(bam, rows)
            bed = "".join(
                f"chr9\t49\t50\tNM:c.(949T>{alt})\tT\t{alt}\tp.Phe317{aa}\n"
                for alt, aa in (("A", "I"), ("C", "L"), ("G", "V"))
            )
            all_rows, variants, _ = run_caller(tmp, bam, bed)
            self.assertEqual([row["status"] for row in all_rows], ["ARTIFACT"] * 3)
            self.assertTrue(all("LOW_SNR" in row["filter_reason"] for row in all_rows))
            self.assertEqual(variants, [])

    def test_codon_haplotype_collection_excludes_insertions_and_reference_skips(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            # An insertion between target bases makes the aligned bases G-T-T
            # look like GTT if the inserted query base is silently ignored.
            insertion = list("N" * 151); insertion[99:103] = list("GATT")
            # N at the middle target is likewise ambiguous even though both
            # flanking bases are present in the query.
            skipped = list("N" * 149); skipped[99] = "G"; skipped[100:102] = list("TT")
            rows = [
                {"seq": "".join(insertion), "cigar": [(0, 100), (1, 1), (0, 50)]},
                {"seq": "".join(skipped), "cigar": [(0, 100), (3, 1), (0, 49)]},
            ]
            bam_path = tmp / "reads.bam"; write_bam(bam_path, rows, length=200)
            with pysam.AlignmentFile(str(bam_path), "rb") as bam:
                self.assertEqual(caller.collect_codon_haplotypes(bam, "chr9", [99, 100, 101], 12, 1), {})

    def test_synonymous_detection_is_auditable_and_not_classed_as_resistance(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            rows = []
            for i, base in enumerate("A" * 20 + "G" * 80):
                seq = list("N" * 150); seq[49] = base; rows.append({"seq": "".join(seq), "name": f"r{i}"})
            bam = tmp / "reads.bam"; write_bam(bam, rows)
            all_rows, variants, _ = run_caller(
                tmp, bam, "chr9\t49\t50\tNM:c.(1161G>A)\tG\tA\tp.Leu387Leu\n"
            )
            self.assertEqual(all_rows[0]["status"], "DETECTED")
            self.assertEqual(all_rows[0]["consequence"], "SYNONYMOUS")
            self.assertEqual(all_rows[0]["resistance_class"], "SYNONYMOUS")
            self.assertEqual(variants[0]["mutation"], "L387L")

    def test_alignment_dedup_excludes_secondary_but_keeps_supplementary_only_signal(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d); rows = []
            rows.append({"name": "one", "seq": "N" * 49 + "A" + "N" * 100, "flag": 0})
            rows.extend({"name": "one", "seq": "N" * 49 + "A" + "N" * 100, "flag": 2048} for _ in range(5))
            rows.extend({"name": "secondary", "seq": "N" * 49 + "A" + "N" * 100, "flag": 256} for _ in range(4))
            rows.append({"name": "supp_only", "seq": "N" * 49 + "A" + "N" * 100, "flag": 2048})
            bam = tmp / "reads.bam"; write_bam(bam, rows)
            with pysam.AlignmentFile(str(bam), "rb") as handle:
                counts = caller.count_position(handle, "chr9", 50, 1000, 12, 1)
            self.assertEqual(counts["A_fwd"], 2)

    def test_explicit_reference_ref_mismatch_fails(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d); rows = []
            for i in range(20):
                seq = list("N" * 150); seq[49] = "A"; rows.append({"seq": "".join(seq)})
            bam = tmp / "reads.bam"; write_bam(bam, rows)
            fasta = tmp / "ref.fa"; fasta.write_text(">chr9\n" + "N" * 49 + "C" + "N" * 150 + "\n"); pysam.faidx(str(fasta))
            bed = tmp / "panel.bed"; bed.write_text("chr9\t49\t50\tNM:c.(763G>A)\tG\tA\tp.Glu255Lys\n")
            with patch.object(caller.sys, "argv", ["caller", "--bam", str(bam), "--bed", str(bed), "--sample", "S", "--output-all", str(tmp / "a.tsv"), "--output-mutations", str(tmp / "m.tsv"), "--reference-fasta", str(fasta)]):
                with self.assertRaisesRegex(ValueError, "BED REF mismatch"):
                    caller.main()


    def test_indel_handling_and_reporting(self):
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            rows = []
            # 20 reference reads (G at pos 50, 0-based 49)
            for _ in range(20):
                seq = list("N" * 150)
                seq[49] = "G"
                rows.append({"seq": "".join(seq), "cigar": [(0, 150)]})
            # 5 reads with a 2bp deletion spanning pos 50 (cigar: 48M 2D 100M)
            for _ in range(5):
                seq = list("N" * 148)
                rows.append({"seq": "".join(seq), "cigar": [(0, 48), (2, 2), (0, 100)]})
            # 6 reads with a 3bp insertion after pos 50 (cigar: 50M 3I 97M -> 150 bp query)
            for _ in range(6):
                seq = list("N" * 150)
                seq[49] = "G"
                rows.append({"seq": "".join(seq), "cigar": [(0, 50), (1, 3), (0, 97)]})
            bam = tmp / "reads.bam"
            write_bam(bam, rows)
            all_rows, _, _ = run_caller(
                tmp, bam, "chr9\t49\t50\tNM:c.(763G>A)\tG\tA\tp.Glu255Lys\n"
            )
            self.assertEqual(len(all_rows), 1)
            # 20 ref + 6 with ins = 26 base calls at pos 50
            self.assertEqual(int(all_rows[0]["depth"]), 26)
            self.assertEqual(int(all_rows[0]["n_deletion"]), 5)
            self.assertEqual(int(all_rows[0]["n_insertion"]), 6)


def pileup_reference(bam, chrom, positions_0based, min_baseq, min_mapq, ref_seqs_by_pos=None, read_filter=None):
    """The pileup implementation the read scan replaced (D25), kept verbatim
    apart from its signature, so the two can be compared on awkward reads."""
    results = {}
    for pos_0based in positions_0based:
        counts = caller.empty_counts()
        results[pos_0based] = counts
    min_pos = min(positions_0based)
    max_pos = max(positions_0based) + 1
    for column in bam.pileup(chrom, min_pos, max_pos, truncate=True, stepper="nofilter",
                             max_depth=1_000_000, min_base_quality=0, ignore_overlaps=False):
        ref_pos = column.reference_pos
        if ref_pos not in positions_0based:
            continue
        counts = results[ref_pos]
        lookup_key = (chrom, ref_pos)
        if ref_seqs_by_pos and lookup_key in ref_seqs_by_pos:
            ref_left, ref_right, left_mask, right_mask = ref_seqs_by_pos[lookup_key]
            k_len = len(ref_left)
        else:
            ref_left, ref_right, left_mask, right_mask, k_len = "", "", (), (), 0
        best_by_read = {}
        for read in column.pileups:
            aln = read.alignment
            if (aln.is_unmapped or aln.is_duplicate or aln.is_qcfail or aln.is_secondary
                    or (aln.mapping_quality or 0) < min_mapq):
                continue
            if read_filter is not None and not read_filter.passes(aln):
                continue
            key = (aln.query_name, aln.is_read1, aln.is_read2)
            rank = (0 if not aln.is_supplementary else 1, -(aln.mapping_quality or 0))
            previous = best_by_read.get(key)
            if previous is None or rank < previous[0]:
                best_by_read[key] = (rank, read)
        for _, read in best_by_read.values():
            aln = read.alignment
            if read.is_refskip:
                counts["refskip"] += 1
                continue
            if read.is_del or read.query_position is None:
                counts["deletion"] += 1
                continue
            bq = aln.query_qualities[read.query_position] if aln.query_qualities is not None else 30
            base = aln.query_sequence[read.query_position].upper()
            if base not in caller.BASES:
                continue
            counts["unfiltered"][base] += 1
            if getattr(read, "indel", 0) > 0:
                counts["insertion"] += 1
            if min_baseq > 0 and bq < min_baseq:
                counts["low_baseq"] += 1
                continue
            clean = None
            if k_len > 0:
                qpos = read.query_position
                qseq = aln.query_sequence
                if qpos >= k_len and (qpos + 1 + k_len) <= len(qseq):
                    counts["kmer_total_count"][base] += 1
                    clean = caller.context_is_clean(qseq[qpos - k_len:qpos].upper(), ref_left, left_mask) and \
                        caller.context_is_clean(qseq[qpos + 1:qpos + 1 + k_len].upper(), ref_right, right_mask)
                    if clean:
                        counts["kmer_clean_count"][base] += 1
            if clean is False:
                counts["context_excluded"][base] += 1
                continue
            strand = "rev" if aln.is_reverse else "fwd"
            counts[f"{base}_{strand}"] += 1
            counts["base_qual_sum"][base] += bq
            counts["base_qual_count"][base] += 1
    return results


class ReadScanMatchesPileup(unittest.TestCase):
    """D25: the read scan must count exactly what the pileup counted."""

    def test_identical_counts_on_awkward_reads(self):
        import random
        rng = random.Random(11)
        ref = "".join(rng.choice("ACGT") for _ in range(400))
        targets = {120, 122, 123, 130, 200, 260}
        rows = []
        for i in range(600):
            start = rng.randrange(0, 150)
            length = rng.randrange(60, 230)
            seq, cigar, pos = [], [], start
            while len(seq) < length and pos < 399:
                roll = rng.random()
                if roll < 0.03:
                    gap = rng.randrange(1, 4)
                    cigar.append((2, gap)); pos += gap
                elif roll < 0.05:
                    ins = "".join(rng.choice("ACGT") for _ in range(rng.randrange(1, 3)))
                    seq.extend(ins); cigar.append((1, len(ins)))
                elif roll < 0.06 and pos < 330:
                    cigar.append((3, 40)); pos += 40
                else:
                    run = min(rng.randrange(3, 25), 399 - pos)
                    for offset in range(run):
                        base = ref[pos + offset]
                        if rng.random() < 0.06:
                            base = rng.choice("ACGT")
                        seq.append(base)
                    cigar.append((0, run)); pos += run
            if not seq or cigar[0][0] not in (0,) or cigar[-1][0] not in (0,):
                continue
            if rng.random() < 0.1:   # = and X instead of M
                cigar = [(7, n) if op == 0 else (op, n) for op, n in cigar]
            quals = [rng.choice((5, 11, 12, 20, 35)) for _ in seq]
            flag = 16 if rng.random() < 0.5 else 0
            name = f"r{i}"
            if rng.random() < 0.05:
                flag |= 256           # secondary: never counted
            rows.append({"name": name, "seq": "".join(seq), "quals": quals, "start": start,
                         "cigar": cigar, "flag": flag, "mapq": rng.choice((0, 1, 60))})
            if rng.random() < 0.08:   # a second alignment of the same read over the same span
                rows.append({"name": name, "seq": "".join(seq), "quals": quals, "start": start,
                             "cigar": cigar, "flag": (flag & 16) | 2048, "mapq": rng.choice((1, 30, 60))})
        rows.sort(key=lambda row: row["start"])
        with tempfile.TemporaryDirectory() as d:
            tmp = Path(d)
            bam = tmp / "awkward.bam"
            write_bam(bam, rows, length=400)
            contexts = {}
            for t in targets:
                left_mask = tuple((t - 3 + i) in targets for i in range(3))
                right_mask = tuple((t + 1 + i) in targets for i in range(3))
                contexts[("chr9", t)] = (ref[t - 3:t], ref[t + 1:t + 4], left_mask, right_mask)
            with pysam.AlignmentFile(str(bam)) as handle:
                expected = pileup_reference(handle, "chr9", targets, 12, 1, contexts)
            with pysam.AlignmentFile(str(bam)) as handle:
                observed = caller.count_positions_batch(handle, "chr9", targets, 0, 12, 1, contexts, "chr9")
        for t in sorted(targets):
            self.assertEqual(observed[t], expected[t], f"position {t}")
        self.assertGreater(sum(expected[t]["deletion"] + expected[t]["refskip"] for t in targets), 0)
        self.assertGreater(sum(expected[t]["insertion"] + sum(expected[t]["context_excluded"].values()) for t in targets), 0)


class NearbyVariants(unittest.TestCase):
    """D25 fixes 1-3, on the cases the review found."""

    def setUp(self):
        import random
        rng = random.Random(7)
        ref = "".join(rng.choice("ACGT") for _ in range(200))
        self.ref = ref[:99] + "GAG" + "C" + ref[103:]      # 99 G, 100 A, 101 G, 102 C
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.tmp = Path(self.directory.name)
        self.fasta = self.tmp / "ref.fa"
        self.fasta.write_text(">chr9\n" + self.ref + "\n")
        pysam.faidx(str(self.fasta))

    def read(self, i, changes=(), quals=None, reverse=None, tags=None):
        seq = list(self.ref[25:175])
        for position, base in changes:
            seq[position - 25] = base
        row = {"name": f"r{i}", "seq": "".join(seq), "start": 25,
               "flag": 16 if (i % 2 if reverse is None else reverse) else 0}
        if quals:
            q = [30] * 150
            for position, value in quals:
                q[position - 25] = value
            row["quals"] = q
        if tags:
            row["tags"] = tags
        return row

    def call(self, rows, bed, extra=()):
        bam = self.tmp / "reads.bam"
        write_bam(bam, rows)
        return run_caller(self.tmp, bam, bed, ("--reference-fasta", str(self.fasta), *extra))

    SITE = "chr9\t100\t101\tNM:c.(1A>G)\tA\tG\tp.Lys1Glu\n"

    def test_a_mutation_in_cis_with_a_non_panel_change_is_called(self):
        rows = [self.read(i, [(100, "G"), (102, "T")]) for i in range(200)]
        rows += [self.read(1000 + i) for i in range(1800)]
        audit, variants, _ = self.call(rows, self.SITE)
        self.assertEqual(audit[0]["status"], "DETECTED")
        self.assertEqual(audit[0]["ad_alt"], "200")
        self.assertEqual(audit[0]["flank_alleles_forgiven"], "103:C>T(sample,10.00%)")

    def test_a_germline_snv_beside_the_site_no_longer_empties_it(self):
        rows = [self.read(i, [(100, "G"), (102, "T")]) for i in range(200)]
        rows += [self.read(1000 + i, [(102, "T")]) for i in range(1800)]
        audit, _, _ = self.call(rows, self.SITE)
        self.assertEqual(audit[0]["depth"], "2000")      # was 0 -> LOW_DEPTH
        self.assertEqual(audit[0]["status"], "DETECTED")

    def test_only_the_neighbours_panel_alt_is_forgiven(self):
        """A non-panel base at a neighbouring panel position is still an error."""
        bed = self.SITE + "chr9\t102\t103\tNM:c.(3C>G)\tC\tG\tp.Lys1Glu\n"
        rows = [self.read(i, [(100, "G"), (102, "T")], reverse=True) for i in range(30)]
        rows += [self.read(1000 + i) for i in range(1970)]
        audit, variants, _ = self.call(rows, bed)
        self.assertEqual(audit[0]["ad_alt"], "0")
        self.assertEqual(audit[0]["ad_alt_context_excluded"], "30")
        self.assertEqual(variants, [])

    def test_the_neighbours_panel_alt_is_forgiven_only_when_it_is_called(self):
        """Compound in cis: both called. Same reads with the neighbour's ALT at
        low base quality: the neighbour is not called, so the two-round check
        does not forgive it. It is then judged as a linked neighbour, and its
        base quality fails one of the two haplotype checks: REVIEW, where it
        used to disappear without a word."""
        bed = self.SITE + "chr9\t102\t103\tNM:c.(3C>G)\tC\tG\tp.Lys1Glu\n"
        rows = [self.read(i, [(100, "G"), (102, "G")]) for i in range(40)]
        rows += [self.read(1000 + i) for i in range(1960)]
        audit, _, _ = self.call(rows, bed)
        self.assertEqual([row["status"] for row in audit], ["DETECTED", "DETECTED"])
        self.assertIn("103:C>G(panel)", audit[0]["flank_alleles_forgiven"])

        rows = [self.read(i, [(100, "G"), (102, "G")], quals=[(102, 13)]) for i in range(40)]
        rows += [self.read(1000 + i) for i in range(1960)]
        audit, _, _ = self.call(rows, bed)
        self.assertEqual(audit[1]["status"], "ARTIFACT")
        self.assertIn("LOW_BASEQ", audit[1]["filter_reason"])
        self.assertNotIn("(panel)", audit[0]["flank_alleles_forgiven"])
        self.assertEqual(audit[0]["status"], "REVIEW")
        self.assertTrue(audit[0]["filter_reason"].startswith("LINKED_NEIGHBOUR_REVIEW(BASEQ("), audit[0]["filter_reason"])
        self.assertTrue(audit[0]["linked_neighbour"].startswith("103:C>G reads=40"), audit[0]["linked_neighbour"])

    def test_codon_noise_at_a_non_panel_base_is_not_a_haplotype(self):
        bed = ("chr9\t99\t100\tNM_005157.6:c.(763G>A)\tG\tA\tp.Glu255Lys\n"
               "chr9\t100\t101\tNM_005157.6:c.(764A>T)\tA\tT\tp.Glu255Val\n")
        rows = [self.read(i, [(99, "A")]) for i in range(100)]
        rows += [self.read(200 + i, [(100, "T")]) for i in range(100)]
        for j, base in enumerate("ACT"):          # 0.4 % each at the third base
            rows += [self.read(400 + 10 * j + i, [(101, base)]) for i in range(8)]
        rows += [self.read(1000 + i) for i in range(2000 - len(rows))]
        _, variants, _ = self.call(rows, bed, ("--include-artifacts",))
        self.assertEqual(sorted(row["mutation"] for row in variants), ["E255K", "E255V"])
        self.assertTrue(all(row["record_type"] == "CODON_HAPLOTYPE" for row in variants))

    def test_a_sample_variant_in_the_codon_is_read_as_one_haplotype(self):
        """T315I (ACT>ATT) with a non-panel T>G at the codon's third base in cis
        is ATG, Met: reported as T315M, not as T315I."""
        self.ref = self.ref[:99] + "ACT" + self.ref[102:]
        self.fasta.write_text(">chr9\n" + self.ref + "\n")
        pysam.faidx(str(self.fasta))
        bed = "chr9\t100\t101\tNM_005157.6:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
        rows = [self.read(i, [(100, "T"), (101, "G")]) for i in range(200)]
        rows += [self.read(1000 + i) for i in range(1800)]
        audit, variants, _ = self.call(rows, bed)
        self.assertEqual(audit[0]["status"], "DETECTED")
        self.assertEqual([(row["record_type"], row["mutation"], row["alt"]) for row in variants],
                         [("CODON_HAPLOTYPE", "T315M", "ATG")])

    def neighbour_base(self, position):
        return "A" if self.ref[position] != "A" else "C"

    def test_a_low_level_mutation_in_cis_with_a_rare_neighbour_is_called(self):
        """D25, linked neighbours: 20 of 4000 reads (0.5 %) carry the panel ALT
        and a non-panel change two bases upstream, outside the codon. Below the
        1 % sample-variant threshold, every ALT read used to be excluded."""
        other = self.neighbour_base(98)
        rows = [self.read(i, [(100, "G"), (98, other)]) for i in range(20)]
        rows += [self.read(1000 + i) for i in range(3980)]
        audit, variants, _ = self.call(rows, self.SITE)
        self.assertEqual(audit[0]["status"], "SUBCLONAL", audit[0]["filter_reason"])
        self.assertEqual(audit[0]["ad_alt"], "20")
        self.assertEqual(audit[0]["flank_alleles_forgiven"], f"99:{self.ref[98]}>{other}(linked,0.50%)")
        self.assertTrue(audit[0]["linked_neighbour"].startswith(f"99:{self.ref[98]}>{other} reads=20 share=1.00"))
        self.assertEqual([(row["record_type"], row["status"]) for row in variants], [("SITE_VARIANT", "SUBCLONAL")])

    def test_a_low_level_codon_change_in_two_bases_is_one_haplotype(self):
        """T315M (ACT>ATG) in one subclone at 0.5 %: c.945T>G is not a panel
        ALT and is below 1 %, so both bases were lost; now it is T315M."""
        self.ref = self.ref[:99] + "ACT" + self.ref[102:]
        self.fasta.write_text(">chr9\n" + self.ref + "\n")
        pysam.faidx(str(self.fasta))
        bed = "chr9\t100\t101\tNM_005157.6:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
        rows = [self.read(i, [(100, "T"), (101, "G")]) for i in range(20)]
        rows += [self.read(1000 + i) for i in range(3980)]
        audit, variants, _ = self.call(rows, bed)
        self.assertEqual(audit[0]["status"], "SUBCLONAL")
        self.assertEqual([(row["record_type"], row["mutation"], row["status"]) for row in variants],
                         [("CODON_HAPLOTYPE", "T315M", "SUBCLONAL")])

    def test_review_rows_reach_the_mutation_table_and_vcf(self):
        """Low base quality at both bases, both strands: REVIEW, not a call."""
        other = self.neighbour_base(98)
        rows = [self.read(i, [(100, "G"), (98, other)], quals=[(100, 15), (98, 15)]) for i in range(20)]
        rows += [self.read(1000 + i) for i in range(3980)]
        audit, variants, vcf = self.call(rows, self.SITE)
        self.assertEqual(audit[0]["status"], "REVIEW")
        self.assertEqual([(row["record_type"], row["status"]) for row in variants], [("SITE_VARIANT", "REVIEW")])
        record = [line for line in vcf.read_text().splitlines() if not line.startswith("#")][0].split("\t")
        self.assertEqual(record[6], "REVIEW")
        self.assertTrue(record[9].startswith("./."))

    def test_noise_beside_a_clonal_mutation_is_not_linked(self):
        """A random error on 0.5 % of reads also sits on ALT reads when the
        ALT is half the reads -- 10 of them here, past --min-alt-reads -- but
        on only 0.5 % of the ALT reads: no linked neighbour, and no false
        two-base haplotype in the codon."""
        rows = []
        for i in range(4000):
            changes = [(100, "G")] if i % 2 == 0 else []
            if i % 400 < 2:
                changes.append((102, "T"))
            rows.append(self.read(i, changes))
        audit, variants, _ = self.call(rows, self.SITE)
        self.assertEqual(audit[0]["status"], "DETECTED")
        self.assertEqual(audit[0]["linked_neighbour"], "")
        self.assertEqual([row["record_type"] for row in variants], ["SITE_VARIANT"])

    def test_soft_clipped_bases_are_not_a_linked_neighbour(self):
        """Review 2026-11-06, F1: read offsets reach past a 5' soft-clip into
        the flank. The window is not gapless, so a clipped base that differs
        is not a substitution and cannot carry the ALT reads."""
        rows = [self.read(1000 + i) for i in range(3980)]
        clipped = "A" if self.ref[98] != "A" else "C"
        for i in range(20):
            seq = list(self.ref[99:190])
            seq[1] = "G"                                    # the ALT at 100
            clip = self.ref[89:98] + clipped               # last clipped base reads as 'position 98'
            rows.append({"name": f"c{i}", "seq": clip + "".join(seq), "start": 99,
                         "cigar": [(4, 10), (0, 91)], "flag": 16 if i % 2 else 0})
        audit, variants, _ = self.call(rows, self.SITE)
        self.assertEqual(audit[0]["status"], "WILDTYPE")
        self.assertEqual(audit[0]["ad_alt_context_excluded"], "20")
        self.assertEqual(audit[0]["linked_neighbour"], "")
        self.assertEqual(variants, [])

    def test_a_review_pair_in_one_codon_is_reviewed_as_that_codon(self):
        """F2: T315M (ACT>ATG) failing the base-quality check is REVIEW -- as
        T315M, ACT>ATG at the codon start, not as T315I."""
        self.ref = self.ref[:99] + "ACT" + self.ref[102:]
        self.fasta.write_text(">chr9\n" + self.ref + "\n")
        pysam.faidx(str(self.fasta))
        bed = "chr9\t100\t101\tNM_005157.6:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
        rows = [self.read(i, [(100, "T"), (101, "G")], quals=[(100, 15), (101, 15)]) for i in range(20)]
        rows += [self.read(1000 + i) for i in range(3980)]
        audit, variants, vcf = self.call(rows, bed)
        self.assertEqual((audit[0]["mutation"], audit[0]["status"]), ("T315I", "REVIEW"))   # the panel allele's audit
        self.assertEqual(audit[0]["codon_context_status"], "LINKED_NEIGHBOUR_IN_CODON")
        self.assertEqual([(v["record_type"], v["mutation"], v["aa_change"], v["ref"], v["alt"], v["status"]) for v in variants],
                         [("CODON_HAPLOTYPE", "T315M", "p.Thr315Met", "ACT", "ATG", "REVIEW")])
        record = [line.split("\t") for line in vcf.read_text().splitlines() if not line.startswith("#")]
        self.assertEqual([r[1:7] for r in record], [["100", "T315M", "ACT", "ATG", ".", "REVIEW"]])

    def test_an_uncalled_panel_alt_can_be_a_linked_neighbour(self):
        """F3: A (c.4) and B (c.2, two bases upstream) are both panel ALTs, on
        the same 20 reads. B fails on its own (noise on a third base at B), so
        the two-round check does not forgive it -- and A used to lose all 20
        reads. Judged as a linked neighbour, the pair passes."""
        b_ref = self.ref[98]
        b_alt = "A" if b_ref != "A" else "C"
        noise = next(base for base in "ACGT" if base not in (b_ref, b_alt))
        bed = (f"chr9\t98\t99\tNM_1:c.(2{b_ref}>{b_alt})\t{b_ref}\t{b_alt}\tp.Xaa1Yaa\n"
               "chr9\t100\t101\tNM_1:c.(4A>G)\tA\tG\tp.Lys2Glu\n")
        rows = [self.read(i, [(100, "G"), (98, b_alt)]) for i in range(20)]
        rows += [self.read(500 + i, [(98, noise)]) for i in range(10)]
        rows += [self.read(1000 + i) for i in range(3970)]
        audit, variants, _ = self.call(rows, bed)
        by_name = {row["mutation"]: row for row in audit}
        self.assertEqual(by_name["Xaa1Yaa"]["status"], "ARTIFACT")
        self.assertEqual(by_name["K2E"]["status"], "SUBCLONAL", by_name["K2E"]["filter_reason"])
        self.assertEqual(by_name["K2E"]["ad_alt"], "20")
        self.assertIn(f"99:{b_ref}>{b_alt}(linked,", by_name["K2E"]["flank_alleles_forgiven"])

    def test_a_low_quality_flank_base_is_neither_allele_in_the_screen(self):
        """F4: 36 reads carry T beside the site at full base quality and 500
        at base quality 5. The 500 are neither T nor the reference there, so
        the sample variant is 36 of 3500 (1.03 %), not 36 of 4000 (0.9 %)."""
        rows = [self.read(i, [(102, "T")], quals=[(102, 5)] if i >= 36 else None) for i in range(536)]
        rows += [self.read(1000 + i) for i in range(3464)]
        audit, _, _ = self.call(rows, self.SITE)
        self.assertIn("103:C>T(sample,1.03%)", audit[0]["flank_alleles_forgiven"])

    def test_fusion_reads_are_counted_separately_and_can_be_the_only_ones(self):
        sa = {"SA": "chr22,23632400,+,200M100S,60,0;"}
        rows = [self.read(i, [(100, "G")] if i < 50 else (), tags=sa) for i in range(1000)]
        rows += [self.read(2000 + i, [(100, "G")] if i < 50 else ()) for i in range(100)]
        region = ("--fusion-region", "chr22:23522552-23660224")
        audit, _, _ = self.call(rows, self.SITE, region)
        row = audit[0]
        self.assertEqual((row["depth"], row["ad_alt"]), ("1100", "100"))
        self.assertEqual((row["depth_fusion"], row["ad_alt_fusion"], row["vaf_fusion_percent"]), ("1000", "50", "5.0"))
        self.assertEqual((row["depth_nonfusion"], row["ad_alt_nonfusion"], row["vaf_nonfusion_percent"]), ("100", "50", "50.0"))
        audit, _, _ = self.call(rows, self.SITE, (*region, "--fusion-reads-only"))
        self.assertEqual((audit[0]["depth"], audit[0]["vaf_percent"]), ("1000", "5.0"))
        header = (self.tmp / "all.tsv").read_text()
        self.assertIn("# reads=fusion fusion_region=chr22:23522552-23660224 fusion_read_fraction=0.909", header)

    def test_fusion_reads_only_needs_a_region(self):
        rows = [self.read(i) for i in range(100)]
        with self.assertRaises(SystemExit):
            self.call(rows, self.SITE, ("--fusion-reads-only",))


if __name__ == "__main__":
    unittest.main()
