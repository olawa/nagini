import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import pysam

HERE = Path(__file__).parent


def load_module(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FUSION_COUNTER = load_module("bcr_abl_fusion_counter")
BRANFORD_CALLER = load_module("bcr_abl_branford_caller")
PHASER = load_module("bcr_abl_phaser")
EXON_COVERAGE = load_module("bcr_abl_exon_coverage")
try:
    REPORT = load_module("nagini_report")
except ImportError:
    REPORT = None




def create_synthetic_bam(bam_path: Path, reads: list[dict], references: list[tuple[str, int]]):
    """Create a synthetic BAM file for testing."""
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": name, "LN": length} for name, length in references]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out_bam:
        for r_idx, r in enumerate(reads):
            aln = pysam.AlignedSegment()
            aln.query_name = r.get("name", f"read_{r_idx}")
            aln.query_sequence = r["seq"]
            aln.query_qualities = r.get("quals") or [r.get("qual", 30)] * len(r["seq"])
            aln.reference_id = [ref[0] for ref in references].index(r["chrom"])
            aln.reference_start = r["start"]  # 0-based
            aln.cigar = r.get("cigar", [(0, len(r["seq"]))])  # default match
            aln.mapping_quality = r.get("mapq", 60)
            aln.is_supplementary = r.get("supplementary", False)
            aln.is_reverse = r.get("is_reverse", False)
            out_bam.write(aln)
    pysam.index(str(bam_path))


class TestBcrAblBranfordCaller(unittest.TestCase):
    def test_short_aa_notation(self):
        self.assertEqual(BRANFORD_CALLER.short_aa_notation("p.Thr315Ile"), "T315I")
        self.assertEqual(BRANFORD_CALLER.short_aa_notation("p.Glu255Lys"), "E255K")
        self.assertEqual(BRANFORD_CALLER.short_aa_notation("p.Tyr253His"), "Y253H")
        self.assertEqual(BRANFORD_CALLER.short_aa_notation("p.Val299Leu"), "V299L")

    def test_branford_calling_on_synthetic_bam(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "test.bam"
            bed_path = tmp / "branford.bed"
            all_sites_tsv = tmp / "all_sites.tsv"
            muts_tsv = tmp / "mutations.tsv"
            vcf_path = tmp / "mutations.vcf"

            # Create test BED: T315 (pos 100) and E255 (pos 50) on test chr9
            bed_path.write_text(
                "chr9\t49\t50\tNM_1:c.(763G>A)\tG\tA\tp.Glu255Lys\n"
                "chr9\t99\t100\tNM_1:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
            )

            # Reads spanning pos 50 (index 49) and pos 100 (index 99)
            # Reference: 150 bp chr9
            # 80 reads with WT (pos 50 = G, pos 100 = C)
            # 20 reads with T315I mutation (pos 100 = T)
            reads = []
            for i in range(80):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "C"
                reads.append({"name": f"wt_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})
            for i in range(20):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "T"
                reads.append({"name": f"mut_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})

            create_synthetic_bam(bam_path, reads, [("chr9", 200)])

            # Run caller main
            BRANFORD_CALLER.sys.argv = [
                "bcr_abl_branford_caller.py",
                "--bam", str(bam_path),
                "--bed", str(bed_path),
                "--sample", "TEST_SAMPLE",
                "--output-all", str(all_sites_tsv),
                "--output-mutations", str(muts_tsv),
                "--output-vcf", str(vcf_path),
                "--min-af", "0.05",
                "--min-depth", "50",
            ]
            BRANFORD_CALLER.main()

            self.assertTrue(all_sites_tsv.is_file())
            self.assertTrue(muts_tsv.is_file())
            self.assertTrue(vcf_path.is_file())

            mut_content = muts_tsv.read_text().splitlines()
            records = [line for line in mut_content if not line.startswith("#") and line.strip()]
            self.assertEqual(len(records), 2)  # Header + 1 detected mutation
            self.assertIn("T315I", records[1])
            self.assertIn("DETECTED", records[1])
            self.assertIn("20.0", records[1])  # 20% VAF

            all_lines = [line for line in all_sites_tsv.read_text().splitlines() if not line.startswith("#") and line.strip()]
            hdr = all_lines[0].split("\t")
            self.assertIn("count_A", hdr)
            self.assertIn("count_T", hdr)
            rows_by_mut = {row["mutation"]: row for row in [dict(zip(hdr, l.split("\t"))) for l in all_lines[1:]]}
            self.assertIn("T315I", rows_by_mut)
            t315i_row = rows_by_mut["T315I"]
            self.assertEqual(t315i_row["count_T"], "20")
            self.assertEqual(t315i_row["count_C"], "80")
            self.assertEqual(t315i_row["pct_T"], "20.0")
            self.assertEqual(t315i_row["pct_C"], "80.0")

    def test_noise_filtering_snr_and_strand_bias(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "test_noise.bam"
            bed_path = tmp / "branford_noise.bed"
            all_sites_tsv = tmp / "all_sites_noise.tsv"
            muts_tsv = tmp / "muts_noise.tsv"

            # 2 sites on chr9:
            # Site 1: E255 (pos 50, G>A) -> will have background noise (A=12, C=10, T=11) out of 1000 reads -> SNR ~ 1.09x -> ARTIFACT
            # Site 2: T315 (pos 100, C>T) -> will have 25 alt reads ALL forward (0 reverse) with Ref balanced (400 fwd, 400 rev) -> STRAND_BIAS -> ARTIFACT
            bed_path.write_text(
                "chr9\t49\t50\tNM_1:c.(763G>A)\tG\tA\tp.Glu255Lys\n"
                "chr9\t99\t100\tNM_1:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
            )

            reads = []
            # Site 1: 967 WT reads (G at 49, C at 99) balanced fwd/rev
            for i in range(484):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "C"
                reads.append({"name": f"wt_fwd_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq), "is_reverse": False})
            for i in range(483):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "C"
                reads.append({"name": f"wt_rev_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq), "is_reverse": True})

            # Site 1 noise: 12 A (candidate alt), 10 C (other 1), 11 T (other 2)
            for i in range(12):
                seq = list("N" * 150)
                seq[49] = "A"
                seq[99] = "C"
                reads.append({"name": f"noise_a_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq), "is_reverse": i % 2 == 0})
            for i in range(10):
                seq = list("N" * 150)
                seq[49] = "C"
                seq[99] = "C"
                reads.append({"name": f"noise_c_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq), "is_reverse": i % 2 == 0})
            for i in range(11):
                seq = list("N" * 150)
                seq[49] = "T"
                seq[99] = "C"
                reads.append({"name": f"noise_t_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq), "is_reverse": i % 2 == 0})

            # Site 2 strand bias: 25 T (candidate alt) ALL forward, 0 reverse
            for i in range(25):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "T"
                reads.append({"name": f"sb_mut_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq), "is_reverse": False})

            create_synthetic_bam(bam_path, reads, [("chr9", 200)])

            BRANFORD_CALLER.sys.argv = [
                "bcr_abl_branford_caller.py",
                "--bam", str(bam_path),
                "--bed", str(bed_path),
                "--sample", "NOISE_SAMPLE",
                "--output-all", str(all_sites_tsv),
                "--output-mutations", str(muts_tsv),
                "--min-af", "0.01",
                "--min-depth", "50",
                "--min-alt-reads", "5",
                "--min-snr", "3.0",
                "--min-minor-strand-ratio", "0.05",
            ]
            BRANFORD_CALLER.main()

            # muts_tsv should have NO detected mutations (both filtered as ARTIFACT)
            mut_content = [l for l in muts_tsv.read_text().splitlines() if not l.startswith("#") and l.strip()]
            self.assertEqual(len(mut_content), 1)  # Only header line, no detected mutations

            # all_sites_tsv should list both with ARTIFACT and specific filter_reasons
            all_content = [l for l in all_sites_tsv.read_text().splitlines() if not l.startswith("#") and l.strip()]
            header = all_content[0].split("\t")
            status_idx = header.index("status")
            reason_idx = header.index("filter_reason")

            rec_e255 = all_content[1].split("\t")
            rec_t315 = all_content[2].split("\t")

            self.assertEqual(rec_e255[status_idx], "ARTIFACT")
            self.assertIn("LOW_SNR", rec_e255[reason_idx])

            self.assertEqual(rec_t315[status_idx], "ARTIFACT")
            self.assertIn("STRAND_BIAS", rec_t315[reason_idx])


class TestBranfordContextAndStrandFilters(unittest.TestCase):
    """The artifact filters must reject artifacts without rejecting biology."""

    REF = ("ACGT" * 50)  # 200 bp, known at every position

    def _write_reference(self, tmp: Path) -> Path:
        fasta = tmp / "ref.fa"
        fasta.write_text(
            ">chr9\n"
            + "\n".join(self.REF[i:i + 60] for i in range(0, len(self.REF), 60))
            + "\n"
        )
        pysam.faidx(str(fasta))
        return fasta

    def _run_caller(self, tmp: Path, bam: Path, bed: Path, fasta: Path, extra=()):
        all_sites = tmp / "all_sites.tsv"
        muts = tmp / "muts.tsv"
        BRANFORD_CALLER.sys.argv = [
            "bcr_abl_branford_caller.py",
            "--bam", str(bam), "--bed", str(bed), "--sample", "S",
            "--output-all", str(all_sites), "--output-mutations", str(muts),
            "--reference-fasta", str(fasta),
            "--min-depth", "50", "--min-alt-reads", "5",
            *extra,
        ]
        BRANFORD_CALLER.main()
        rows = []
        header = None
        for line in all_sites.read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split("\t")
            if header is None:
                header = fields
                continue
            rows.append(dict(zip(header, fields)))
        return rows

    def test_kmer_filter_exempts_a_neighbouring_panel_position(self):
        """A compound mutation in cis is not a dirty context.

        25 of the 52 real Branford positions sit within 3 bp of another one, so
        without this exemption two mutations on the same molecule would mark
        each other's flank dirty and both would be rejected before the phaser
        could resolve whether they are in cis.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = self._write_reference(tmp)
            # 0-based 99 and 100 are adjacent; both are panel positions.
            ref_a, ref_b = self.REF[99], self.REF[100]
            alt_a = "A" if ref_a != "A" else "C"
            alt_b = "A" if ref_b != "A" else "C"
            bed = tmp / "panel.bed"
            bed.write_text(
                f"chr9\t99\t100\tNM_1:c.(1{ref_a}>{alt_a})\t{ref_a}\t{alt_a}\tp.Leu387Phe\n"
                f"chr9\t100\t101\tNM_1:c.(2{ref_b}>{alt_b})\t{ref_b}\t{alt_b}\tp.Met388Leu\n"
            )

            reads = []
            for i in range(300):                      # wild type
                reads.append({"name": f"wt{i}", "chrom": "chr9", "start": 0,
                              "seq": self.REF[:150], "is_reverse": i % 2 == 0})
            for i in range(200):                      # both variants, same molecule
                seq = list(self.REF[:150])
                seq[99], seq[100] = alt_a, alt_b
                reads.append({"name": f"cis{i}", "chrom": "chr9", "start": 0,
                              "seq": "".join(seq), "is_reverse": i % 2 == 0})

            bam = tmp / "cis.bam"
            create_synthetic_bam(bam, reads, [("chr9", 200)])
            rows = self._run_caller(tmp, bam, bed, fasta)

        for row in rows:
            self.assertEqual(
                row["kmer_clean_percent"], "100.0",
                f"{row['pos']}: the neighbouring panel position should be exempt",
            )
            self.assertNotIn("DIRTY_CONTEXT", row["filter_reason"])
            self.assertEqual(row["status"], "DETECTED", row["filter_reason"])

    def test_kmer_filter_still_rejects_a_genuinely_dirty_context(self):
        """A mismatch that is not a known position still counts as dirty."""

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = self._write_reference(tmp)
            ref_base = self.REF[99]
            alt_base = "A" if ref_base != "A" else "C"
            bed = tmp / "panel.bed"
            bed.write_text(
                f"chr9\t99\t100\tNM_1:c.(1{ref_base}>{alt_base})\t{ref_base}\t{alt_base}\tp.Thr315Ile\n"
            )

            reads = []
            for i in range(300):
                reads.append({"name": f"wt{i}", "chrom": "chr9", "start": 0,
                              "seq": self.REF[:150], "is_reverse": i % 2 == 0})
            for i in range(200):
                seq = list(self.REF[:150])
                seq[99] = alt_base
                # 2 bp downstream, not a panel position: genuine context error.
                seq[101] = "A" if self.REF[101] != "A" else "G"
                reads.append({"name": f"dirty{i}", "chrom": "chr9", "start": 0,
                              "seq": "".join(seq), "is_reverse": i % 2 == 0})

            bam = tmp / "dirty.bam"
            create_synthetic_bam(bam, reads, [("chr9", 200)])
            rows = self._run_caller(tmp, bam, bed, fasta)

        self.assertEqual(rows[0]["kmer_clean_percent"], "0.0")
        self.assertEqual(rows[0]["status"], "ARTIFACT")
        self.assertIn("DIRTY_CONTEXT", rows[0]["filter_reason"])

    def test_strand_bias_spares_a_variant_on_a_strand_skewed_amplicon(self):
        """When the whole amplicon is one-sided the reference is one-sided too.

        The minor-strand ratio alone would reject a real variant here; SOR sees
        that the alt balance matches the reference's own and lets it through.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = self._write_reference(tmp)
            ref_base = self.REF[99]
            alt_base = "A" if ref_base != "A" else "C"
            bed = tmp / "panel.bed"
            bed.write_text(
                f"chr9\t99\t100\tNM_1:c.(1{ref_base}>{alt_base})\t{ref_base}\t{alt_base}\tp.Thr315Ile\n"
            )

            reads = []
            for i in range(380):                       # reference, 97% forward
                reads.append({"name": f"wt_f{i}", "chrom": "chr9", "start": 0,
                              "seq": self.REF[:150], "is_reverse": False})
            for i in range(12):
                reads.append({"name": f"wt_r{i}", "chrom": "chr9", "start": 0,
                              "seq": self.REF[:150], "is_reverse": True})
            for i in range(100):                       # alt, same 97% skew
                seq = list(self.REF[:150]); seq[99] = alt_base
                reads.append({"name": f"alt_f{i}", "chrom": "chr9", "start": 0,
                              "seq": "".join(seq), "is_reverse": False})
            for i in range(3):
                seq = list(self.REF[:150]); seq[99] = alt_base
                reads.append({"name": f"alt_r{i}", "chrom": "chr9", "start": 0,
                              "seq": "".join(seq), "is_reverse": True})

            bam = tmp / "skew.bam"
            create_synthetic_bam(bam, reads, [("chr9", 200)])
            rows = self._run_caller(tmp, bam, bed, fasta)

        row = rows[0]
        self.assertLess(float(row["alt_minor_strand_ratio"]), 0.05,
                        "the ratio alone should look one-sided here")
        self.assertLessEqual(float(row["sor"]), 3.0,
                             "SOR should see the alt matches the reference's skew")
        self.assertNotIn("STRAND_BIAS", row["filter_reason"])
        self.assertEqual(row["status"], "DETECTED", row["filter_reason"])

    def test_one_sided_alt_against_balanced_reference_is_still_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = self._write_reference(tmp)
            ref_base = self.REF[99]
            alt_base = "A" if ref_base != "A" else "C"
            bed = tmp / "panel.bed"
            bed.write_text(
                f"chr9\t99\t100\tNM_1:c.(1{ref_base}>{alt_base})\t{ref_base}\t{alt_base}\tp.Thr315Ile\n"
            )

            reads = []
            for i in range(400):
                reads.append({"name": f"wt{i}", "chrom": "chr9", "start": 0,
                              "seq": self.REF[:150], "is_reverse": i % 2 == 0})
            for i in range(60):                        # every alt read forward
                seq = list(self.REF[:150]); seq[99] = alt_base
                reads.append({"name": f"alt{i}", "chrom": "chr9", "start": 0,
                              "seq": "".join(seq), "is_reverse": False})

            bam = tmp / "onesided.bam"
            create_synthetic_bam(bam, reads, [("chr9", 200)])
            rows = self._run_caller(tmp, bam, bed, fasta)

        self.assertEqual(rows[0]["status"], "ARTIFACT")
        self.assertIn("STRAND_BIAS", rows[0]["filter_reason"])

    def test_unfiltered_counts_do_not_move_with_the_base_quality_filter(self):
        """VAF depends on --min-base-quality; the comparison columns must not."""

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            fasta = self._write_reference(tmp)
            ref_base = self.REF[99]
            alt_base = "A" if ref_base != "A" else "C"
            bed = tmp / "panel.bed"
            bed.write_text(
                f"chr9\t99\t100\tNM_1:c.(1{ref_base}>{alt_base})\t{ref_base}\t{alt_base}\tp.Thr315Ile\n"
            )

            reads = []
            for i in range(300):
                reads.append({"name": f"wt{i}", "chrom": "chr9", "start": 0,
                              "seq": self.REF[:150], "is_reverse": i % 2 == 0, "qual": 30})
            for i in range(60):                        # good-quality alt
                seq = list(self.REF[:150]); seq[99] = alt_base
                reads.append({"name": f"good{i}", "chrom": "chr9", "start": 0,
                              "seq": "".join(seq), "is_reverse": i % 2 == 0, "qual": 30})
            for i in range(40):                        # alt below the Q12 floor
                seq = list(self.REF[:150]); seq[99] = alt_base
                reads.append({"name": f"lowq{i}", "chrom": "chr9", "start": 0,
                              "seq": "".join(seq), "is_reverse": i % 2 == 0, "qual": 5})

            bam = tmp / "baseq.bam"
            create_synthetic_bam(bam, reads, [("chr9", 200)])
            strict = self._run_caller(tmp, bam, bed, fasta, ("--min-base-quality", "12"))[0]
            loose = self._run_caller(tmp, bam, bed, fasta, ("--min-base-quality", "0"))[0]

        # The filtered counts move with the threshold ...
        self.assertEqual(strict["ad_alt"], "60")
        self.assertEqual(loose["ad_alt"], "100")
        # ... the unfiltered ones do not, whichever threshold was used.
        self.assertEqual(strict["ad_alt_unfiltered"], "100")
        self.assertEqual(loose["ad_alt_unfiltered"], "100")
        self.assertEqual(strict["depth_unfiltered"], loose["depth_unfiltered"])


class TestBcrAblPhaser(unittest.TestCase):

    def test_single_molecule_cis_phasing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "test_phasing.bam"
            muts_tsv = tmp / "mutations.tsv"
            phasing_tsv = tmp / "phasing.tsv"
            phasing_json = tmp / "phasing.json"

            # 2 mutations: E255K at pos 50 (G>A), T315I at pos 100 (C>T)
            muts_tsv.write_text(
                "sample\tchrom\tpos\tref\talt\tmutation\taa_change\tname\tdepth\tad_ref\tad_alt\tvaf_percent\tstatus\n"
                "SAMPLE_1\tchr9\t50\tG\tA\tE255K\tp.Glu255Lys\tNM_1\t100\t80\t20\t20.0\tDETECTED\n"
                "SAMPLE_1\tchr9\t100\tC\tT\tT315I\tp.Thr315Ile\tNM_1\t100\t80\t20\t20.0\tDETECTED\n"
            )

            # Create 20 reads with BOTH E255K (A at 49) and T315I (T at 99) -> CIS
            # and 80 WT reads (G at 49, C at 99)
            reads = []
            for i in range(80):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "C"
                reads.append({"name": f"wt_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})
            for i in range(20):
                seq = list("N" * 150)
                seq[49] = "A"
                seq[99] = "T"
                reads.append({"name": f"cis_mut_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})

            create_synthetic_bam(bam_path, reads, [("chr9", 200)])

            PHASER.sys.argv = [
                "bcr_abl_phaser.py",
                "--bam", str(bam_path),
                "--mutations-tsv", str(muts_tsv),
                "--sample", "SAMPLE_1",
                "--output-tsv", str(phasing_tsv),
                "--output-json", str(phasing_json),
                "--min-spanning-reads", "10",
            ]
            PHASER.main()

            self.assertTrue(phasing_tsv.is_file())
            content = phasing_tsv.read_text().splitlines()
            records = [line for line in content if not line.startswith("#") and line.strip()]
            self.assertEqual(len(records), 2)  # Header + 1 pair
            self.assertIn("CIS (Compound)", records[1])

            with phasing_json.open() as h:
                data = json.load(h)
            self.assertEqual(data["phased_pairs_count"], 1)
            pair = data["phased_pairs"][0]
            self.assertEqual(pair["cis_reads_11"], 20)
            self.assertEqual(pair["trans_reads_10"], 0)
            self.assertEqual(pair["trans_reads_01"], 0)
            self.assertEqual(pair["cis_fraction_percent"], 100.0)
            self.assertEqual(pair["cis_vaf_percent"], 20.0)
            self.assertEqual(pair["phase_state"], "CIS (Compound)")

    def test_single_molecule_trans_phasing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "test_trans.bam"
            muts_tsv = tmp / "mutations.tsv"
            phasing_tsv = tmp / "phasing.tsv"
            phasing_json = tmp / "phasing.json"

            muts_tsv.write_text(
                "sample\tchrom\tpos\tref\talt\tmutation\taa_change\tname\tdepth\tad_ref\tad_alt\tvaf_percent\tstatus\n"
                "SAMPLE_2\tchr9\t50\tG\tA\tE255K\tp.Glu255Lys\tNM_1\t100\t90\t10\t10.0\tDETECTED\n"
                "SAMPLE_2\tchr9\t100\tC\tT\tT315I\tp.Thr315Ile\tNM_1\t100\t90\t10\t10.0\tDETECTED\n"
            )

            # 10 reads with E255K alone (A at 49, C at 99)
            # 10 reads with T315I alone (G at 49, T at 99)
            # 80 reads with WT (G at 49, C at 99) -> TRANS
            reads = []
            for i in range(80):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "C"
                reads.append({"name": f"wt_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})
            for i in range(10):
                seq = list("N" * 150)
                seq[49] = "A"
                seq[99] = "C"
                reads.append({"name": f"mut1_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})
            for i in range(10):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "T"
                reads.append({"name": f"mut2_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})

            create_synthetic_bam(bam_path, reads, [("chr9", 200)])

            PHASER.sys.argv = [
                "bcr_abl_phaser.py",
                "--bam", str(bam_path),
                "--mutations-tsv", str(muts_tsv),
                "--sample", "SAMPLE_2",
                "--output-tsv", str(phasing_tsv),
                "--output-json", str(phasing_json),
            ]
            PHASER.main()

            with phasing_json.open() as h:
                data = json.load(h)
            pair = data["phased_pairs"][0]
            self.assertEqual(pair["cis_reads_11"], 0)
            self.assertEqual(pair["trans_reads_10"], 10)
            self.assertEqual(pair["trans_reads_01"], 10)
            self.assertEqual(pair["cis_fraction_percent"], 0.0)
            self.assertEqual(pair["cis_vaf_percent"], 0.0)
            self.assertEqual(pair["phase_state"], "TRANS (Polyclonal)")


class TestBcrAblFusionCounter(unittest.TestCase):
    def test_assign_exon(self):
        exon_map = {"e1": 100, "e13": 500, "e14": 1000}
        self.assertEqual(FUSION_COUNTER.assign_exon(105, exon_map, tolerance=50), "e1")
        self.assertEqual(FUSION_COUNTER.assign_exon(510, exon_map, tolerance=50), "e13")
        self.assertIsNone(FUSION_COUNTER.assign_exon(750, exon_map, tolerance=50))

    def test_assign_bcr_exon_minor_breakpoint(self):
        bcr_exons = FUSION_COUNTER.KNOWN_JUNCTIONS["hg19"]["bcr_exons"]
        bcr_regions = FUSION_COUNTER.KNOWN_JUNCTIONS["hg19"]["bcr_regions"]

        # Case 1: Exact donor match (23524426)
        exon, pos = FUSION_COUNTER.assign_bcr_exon(
            23524426, [{"start": 23524226, "end": 23524426}], bcr_exons, bcr_regions, tolerance=150
        )
        self.assertEqual(exon, "e1")

        # Case 2: Soft-clipped BCR segment ending early at 23524270 (156 bp from donor)
        exon, pos = FUSION_COUNTER.assign_bcr_exon(
            23524270, [{"start": 23524226, "end": 23524270}], bcr_exons, bcr_regions, tolerance=150
        )
        self.assertEqual(exon, "e1")
        self.assertEqual(pos, 23524426)

    def test_assign_abl_exon_a2_vs_a3(self):
        abl_exons = FUSION_COUNTER.KNOWN_JUNCTIONS["hg19"]["abl_exons"]
        abl_regions = FUSION_COUNTER.KNOWN_JUNCTIONS["hg19"]["abl_regions"]

        # Case 1: Exon 3 without Exon 2 (e1a3 / p185)
        # Alignment starts inside Exon 3 at 133730250 (>150 bp from acceptor 133730188)
        exon, pos = FUSION_COUNTER.assign_abl_exon(
            133730250, [{"start": 133730250, "end": 133758000, "blocks": [(133730250, 133730483), (133738150, 133738422)]}],
            abl_exons, abl_regions, tolerance=150
        )
        self.assertEqual(exon, "a3")
        self.assertEqual(pos, 133730188)

        # Case 2: Exon 2 present (e1a2 / e14a2 / p210)
        exon, pos = FUSION_COUNTER.assign_abl_exon(
            133729451, [{"start": 133729451, "end": 133758000, "blocks": [(133729451, 133729624), (133730188, 133730483)]}],
            abl_exons, abl_regions, tolerance=150
        )
        self.assertEqual(exon, "a2")
        self.assertEqual(pos, 133729451)

    def test_e1a3_fusion_chimeric_bam_calling(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "e1a3_test.bam"
            out_tsv = tmp / "fusion_summary.tsv"
            out_json = tmp / "fusion_summary.json"

            # Create synthetic e1a3 chimeric reads:
            # Segment 1 on chr22: BCR exon 1 (start 23524225, end 23524300, soft-clipped)
            # Segment 2 on chr9: ABL1 exon 3 (start 133730187, end 133732000)
            reads = []
            for i in range(20):
                reads.append({
                    "name": f"mol_{i}", "chrom": "chr9", "start": 133730187,
                    "seq": "A" * 500, "qual": 30, "supplementary": False
                })
            for i in range(20):
                reads.append({
                    "name": f"mol_{i}", "chrom": "chr22", "start": 23524225,
                    "seq": "C" * 75, "qual": 30, "supplementary": True
                })

            refs = [("chr9", 141213431), ("chr22", 51304566)]
            create_synthetic_bam(bam_path, reads, refs)

            counts, events, scanned, stats = FUSION_COUNTER.extract_fusion_events(
                bam_path=bam_path, build="hg19", tolerance=150, min_mapq=0
            )

            self.assertEqual(counts.get("e1a3"), 20)
            self.assertEqual(counts.get("atypical_bcr_abl", 0), 0)
            self.assertNotIn("e1a2_candidate", counts)

            FUSION_COUNTER.write_outputs(
                sample="SAMPLE_BC08",
                isoform_counts=counts,
                raw_events=events,
                output_path=out_tsv,
                json_path=out_json,
                build="hg19",
                read_stats=stats,
            )

            lines = [l for l in out_tsv.read_text().splitlines() if not l.startswith("#")]
            header = lines[0].split("\t")
            row1 = dict(zip(header, lines[1].split("\t")))
            self.assertEqual(row1["fusion_isoform"], "e1a3")
            self.assertEqual(row1["bcr_exon"], "e1")
            self.assertEqual(row1["abl_exon"], "a3")
            self.assertEqual(row1["protein_form"], "p185")
            self.assertEqual(row1["cluster"], "m-BCR")
            self.assertEqual(row1["status"], "PRIMARY")
            self.assertEqual(row1["read_count"], "20")

    def test_read_stats_computation_on_bam(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "stats_test.bam"
            reads = [
                {"name": "r1", "chrom": "chr9", "start": 10, "seq": "A" * 1000, "qual": 25},
                {"name": "r2", "chrom": "chr9", "start": 20, "seq": "C" * 1500, "qual": 30},
                {"name": "r3", "chrom": "chr9", "start": 30, "seq": "G" * 2000, "qual": 35},
                # supplementary alignment of r3 (should NOT count towards total_reads or yield)
                {"name": "r3", "chrom": "chr9", "start": 40, "seq": "G" * 500, "qual": 35, "supplementary": True},
            ]
            create_synthetic_bam(bam_path, reads, [("chr9", 3000)])
            _, _, total_scanned, stats = FUSION_COUNTER.extract_fusion_events(
                bam_path=bam_path,
                build="hg19",
                tolerance=50,
                min_mapq=0,
            )
            self.assertEqual(stats["total_reads"], 3)
            self.assertEqual(stats["median_bp"], 1500.0)
            self.assertEqual(stats["n50_bp"], 1500)
            self.assertAlmostEqual(stats["yield_gb"], 4500 / 1e9, places=6)



@unittest.skipIf(REPORT is None or getattr(REPORT, "xlsxwriter", None) is None, "xlsxwriter not available")
class TestBcrAblReportIntegration(unittest.TestCase):
    def test_generates_workbook_with_bcr_abl_tab(self):

        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            out_xlsx = tmp / "report.xlsx"
            fusion_tsv = tmp / "fusion.tsv"
            muts_tsv = tmp / "muts.tsv"
            phasing_tsv = tmp / "phasing.tsv"

            fusion_tsv.write_text(
                "sample\tfusion_isoform\talt_name\tprotein_form\tcluster\tbcr_exon\tabl_exon\tread_count\tfraction_percent\tstatus\tdescription\n"
                "SAMPLE_1\te14a2\tb3a2\tp210\tM-BCR\te14\ta2\t1500\t95.5\tPRIMARY\tMajor CML fusion\n"
            )
            muts_tsv.write_text(
                "sample\tchrom\tpos\tref\talt\tmutation\taa_change\tname\tdepth\tad_ref\tad_alt\tvaf_percent\tstatus\timatinib\tdasatinib\tnilotinib\tbosutinib\tponatinib\tasciminib\n"
                "SAMPLE_1\tchr9\t130872896\tC\tT\tT315I\tp.Thr315Ile\tNM_1\t5000\t4000\t1000\t20.0\tDETECTED\tResistant\tResistant\tResistant\tResistant\tSensitive\tSensitive\n"
            )
            phasing_tsv.write_text(
                "sample\tmutation_1\tmutation_2\tphase_state\tcis_fraction_percent\tcis_vaf_percent\tcis_reads_11\ttrans_reads_10\ttrans_reads_01\twt_reads_00\tother_reads\ttotal_spanning_reads\tpos_1\tpos_2\tref_1\talt_1\tref_2\talt_2\tclinical_interpretation\n"
            )

            exon_cov_tsv = tmp / "exon_cov.tsv"
            exon_cov_tsv.write_text(
                "sample\texon\tchrom\tstart\tend\tlength_bp\tmean_depth\tmedian_depth\tmin_depth\tmax_depth\tpct_gte_100x\tpct_gte_500x\tpct_gte_1000x\tspanning_reads\n"
                "SAMPLE_1\tABL1_Exon2\tchr9\t130854065\t130854237\t173\t4500.0\t4500.0\t3200\t5800\t100.0\t100.0\t100.0\t4500\n"
            )

            REPORT.generate_report(
                output=str(out_xlsx),
                sample="SAMPLE_1",
                sample_type="T",
                reference_build="hg38",
                annotation="none",
                snv_vcf="",
                sv_tsv="",
                translocations_tsv="",
                wgs_coverage="",
                target_coverage="",
                cramino="",
                bcr_abl_fusion=str(fusion_tsv),
                bcr_abl_mutations=str(muts_tsv),
                bcr_abl_phasing=str(phasing_tsv),
                bcr_abl_exon_coverage=str(exon_cov_tsv),
            )

            self.assertTrue(out_xlsx.is_file())
            self.assertGreater(out_xlsx.stat().st_size, 1000)


class TestBcrAblExonCoverage(unittest.TestCase):
    def test_calculate_exon_coverage_on_synthetic_bam(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "test_exon.bam"
            bed_path = tmp / "exons.bed"
            out_tsv = tmp / "exon_cov.tsv"

            bed_path.write_text(
                "chr9\t10\t50\tABL1_Exon2\n"
                "chr9\t60\t100\tABL1_Exon3\n"
            )

            # 50 reads covering Exon 2 (10..50)
            reads = []
            for i in range(50):
                reads.append({
                    "name": f"read_e2_{i}",
                    "chrom": "chr9",
                    "start": 10,
                    "seq": "A" * 40,
                })
            # 20 reads covering Exon 3 (60..100)
            for i in range(20):
                reads.append({
                    "name": f"read_e3_{i}",
                    "chrom": "chr9",
                    "start": 60,
                    "seq": "C" * 40,
                })

            create_synthetic_bam(bam_path, reads, [("chr9", 200)])

            exons = EXON_COVERAGE.read_exon_bed(bed_path)
            self.assertEqual(len(exons), 2)

            results = EXON_COVERAGE.calculate_exon_coverage(
                bam_path=bam_path,
                exons=exons,
                sample="EXON_TEST",
            )
            self.assertEqual(len(results), 2)
            self.assertEqual(results[0]["exon"], "ABL1_Exon2")
            self.assertEqual(results[0]["mean_depth"], 50.0)
            self.assertEqual(results[0]["spanning_reads"], 50)
            self.assertEqual(results[1]["exon"], "ABL1_Exon3")
            self.assertEqual(results[1]["mean_depth"], 20.0)
            self.assertEqual(results[1]["spanning_reads"], 20)

            EXON_COVERAGE.write_exon_coverage_tsv(results, out_tsv)
            self.assertTrue(out_tsv.is_file())


class TestBcrAblArribaOutput(unittest.TestCase):
    def test_write_arriba_tsv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            arriba_tsv = tmp / "arriba.tsv"

            isoform_counts = {"e14a2": 500, "e13a2": 50}
            raw_events = [
                {"isoform": "e14a2", "bcr_pos": 23290413, "abl_pos": 130854065, "bcr_exon": "e14", "abl_exon": "a2", "read_name": "r1"},
                {"isoform": "e13a2", "bcr_pos": 23289621, "abl_pos": 130854065, "bcr_exon": "e13", "abl_exon": "a2", "read_name": "r2"},
            ]

            FUSION_COUNTER.write_arriba_tsv(
                sample="TEST_S",
                isoform_counts=isoform_counts,
                raw_events=raw_events,
                arriba_path=arriba_tsv,
                build="hg38",
            )

            self.assertTrue(arriba_tsv.is_file())
            content = arriba_tsv.read_text().splitlines()
            self.assertGreaterEqual(len(content), 3)
            self.assertTrue(content[0].startswith("#gene1\tgene2"))
            self.assertIn("BCR\tABL1", content[1])
            self.assertIn("chr22:23290413", content[1])
            self.assertIn("chr9:130854065", content[1])


class TestBcrAblReport(unittest.TestCase):
    def test_build_bcr_abl_sheet_with_metrics(self):
        if REPORT is None:
            self.skipTest("nagini_report module could not be imported")
        fusion_rows = [{"fusion_isoform": "e14a2", "status": "PRIMARY"}]
        mutation_rows = [
            {
                "mutation": "T315I",
                "aa_change": "p.Thr315Ile",
                "pos": 133748283,
                "ref": "C",
                "alt": "T",
                "depth": 46000,
                "ad_alt": 650,
                "vaf_percent": 1.41,
                "status": "DETECTED",
                "snr": 12.5,
                "alt_minor_strand_ratio": 0.42,
                "mean_baseq_alt": 28.4,
                "imatinib": "Resistant",
            },
            {
                "mutation": "Y253H",
                "aa_change": "p.Tyr253His",
                "pos": 133738357,
                "ref": "T",
                "alt": "C",
                "depth": 45000,
                "ad_alt": 200,
                "vaf_percent": 0.44,
                "status": "SUBCLONAL",
                "snr": 5.2,
                "alt_minor_strand_ratio": 0.38,
                "mean_baseq_alt": 25.1,
                "imatinib": "Resistant",
            },
        ]
        phasing_rows = []
        exon_rows = []

        layout = REPORT.build_bcr_abl_sheet(
            sample="TEST_BCR_SAMPLE",
            fusion_rows=fusion_rows,
            mutation_rows=mutation_rows,
            phasing_rows=phasing_rows,
            exon_rows=exon_rows,
        )
        self.assertIsInstance(layout, list)

        # Check header
        sec2_header_row = [row for row in layout if len(row) > 0 and row[0][0] == "Mutation"][0]
        header_names = [col[0] for col in sec2_header_row]
        self.assertIn("SNR (x)", header_names)
        self.assertIn("Strand Ratio", header_names)
        self.assertIn("Mean BaseQ", header_names)

        # Check rows
        t315i_row = [row for row in layout if len(row) > 0 and row[0][0] == "T315I"][0]
        y253h_row = [row for row in layout if len(row) > 0 and row[0][0] == "Y253H"][0]

        snr_idx = header_names.index("SNR (x)")
        strand_idx = header_names.index("Strand Ratio")
        baseq_idx = header_names.index("Mean BaseQ")
        status_idx = header_names.index("Status")

        self.assertEqual(t315i_row[snr_idx][0], "12.5x")
        self.assertEqual(t315i_row[strand_idx][0], "0.42")
        self.assertEqual(t315i_row[baseq_idx][0], "28.4")
        self.assertEqual(t315i_row[status_idx][0], "DETECTED")
        self.assertEqual(t315i_row[status_idx][1], 4)

        self.assertEqual(y253h_row[status_idx][0], "SUBCLONAL")
        self.assertEqual(y253h_row[status_idx][1], 6)

    def test_orthogonal_caller_supports_t315i_mnv(self):
        if REPORT is None:
            self.skipTest("nagini_report module could not be imported")
        fusion_rows = [{"fusion_isoform": "e14a2", "status": "PRIMARY"}]
        # Branford reports MNV at 133748282 (ACT>ATT)
        mutation_rows = [
            {
                "mutation": "T315I",
                "aa_change": "p.Thr315Ile",
                "pos": 133748282,
                "ref": "ACT",
                "alt": "ATT",
                "depth": 46000,
                "ad_alt": 650,
                "vaf_percent": 1.41,
                "status": "DETECTED",
                "snr": 12.5,
                "alt_minor_strand_ratio": 0.42,
                "mean_baseq_alt": 28.4,
                "imatinib": "Resistant",
            },
        ]
        # Orthogonal callers report SNV at 133748283 (C>T)
        other_caller_snvs = [
            {
                "Caller": "VarDict",
                "Chrom": "chr9",
                "Position": 133748283,
                "Reference": "C",
                "Alternate": "T",
                "VAF": 0.014,
                "Filter": "PASS",
                "HGVSp": "p.Thr315Ile",
            },
            {
                "Caller": "bcftools",
                "Chrom": "chr9",
                "Position": 133748283,
                "Reference": "C",
                "Alternate": "T",
                "VAF": 0.014,
                "Filter": "PASS",
            },
        ]

        layout = REPORT.build_bcr_abl_sheet(
            sample="TEST_BCR_SAMPLE",
            fusion_rows=fusion_rows,
            mutation_rows=mutation_rows,
            phasing_rows=[],
            exon_rows=[],
            other_caller_snvs=other_caller_snvs,
        )

        sec2_header_row = [row for row in layout if len(row) > 0 and row[0][0] == "Mutation"][0]
        header_names = [col[0] for col in sec2_header_row]
        support_idx = header_names.index("Supporting Callers")

        t315i_row = [row for row in layout if len(row) > 0 and row[0][0] == "T315I"][0]
        support_text = t315i_row[support_idx][0]
        self.assertIn("VarDict (1.40%)", support_text)
        self.assertIn("bcftools (1.40%)", support_text)
        self.assertNotIn("Branford only", support_text)

        # Verify that 133748283 is NOT listed in "ADDITIONAL ABL1 EXON VARIANTS"
        additional_section = [row for row in layout if len(row) > 0 and "ADDITIONAL ABL1 EXON VARIANTS" in str(row[0][0])]
        self.assertEqual(len(additional_section), 0)


if __name__ == "__main__":
    unittest.main()


