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


if __name__ == "__main__":
    unittest.main()
