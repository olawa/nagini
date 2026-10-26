#!/usr/bin/env python3
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

import pysam

try:
    import matplotlib  # noqa: F401  (the sweep and the plotter need it)
except ImportError:
    matplotlib = None

if matplotlib is not None:
    import sweep_bcr_abl_qs as QS
else:
    QS = None


@unittest.skipIf(matplotlib is None, "matplotlib is required for the QS sweep and background plots")
class TestSweepBcrAblQs(unittest.TestCase):
    def test_load_branford_bed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bed_path = tmp / "branford_test.bed"
            bed_path.write_text(
                "chr9\t130862969\t130862970\tNM_005157.6:c.(757T>C)\tT\tC\tp.Tyr253His\n"
                "chr9\t130872895\t130872896\tNM_005157.6:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
            )
            sites = QS.load_branford_bed(bed_path)
            self.assertEqual(len(sites), 2)
            self.assertEqual(sites[0]["mutation"], "Y253H")
            self.assertEqual(sites[0]["ref"], "T")
            self.assertEqual(sites[0]["alt"], "C")
            self.assertEqual(sites[0]["pos0"], 130862969)
            self.assertEqual(sites[1]["mutation"], "T315I")

    def test_collect_qs_pileup_and_generate_report(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bed_path = tmp / "branford_test.bed"
            bed_path.write_text(
                "chr9\t100\t101\tNM_1\tA\tG\tp.Lys100Arg\n"
                "chr9\t200\t201\tNM_2\tC\tT\tp.Thr200Ile\n"
            )
            sites = QS.load_branford_bed(bed_path)

            bam_path = tmp / "test.bam"
            header = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 1000, "SN": "chr9"}]}

            # Build synthetic BAM
            with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out_bam:
                # Read 1: covers pos 100 with A (ref), qs = 19.5 (int 19 -> <20)
                a1 = pysam.AlignedSegment()
                a1.query_name = "r1"
                a1.query_sequence = "A" * 50
                a1.query_qualities = pysam.qualitystring_to_array("I" * 50)
                a1.reference_id = 0
                a1.reference_start = 80
                a1.cigar = ((0, 50),)
                a1.mapping_quality = 60
                a1.set_tag("qs", 19.5, value_type="f")
                out_bam.write(a1)

                # Read 2: covers pos 100 with G (alt), qs = 22.1 (int 22)
                a2 = pysam.AlignedSegment()
                a2.query_name = "r2"
                # at pos 100 (offset 100-90 = 10): put G
                seq = list("A" * 50)
                seq[10] = "G"
                a2.query_sequence = "".join(seq)
                a2.query_qualities = pysam.qualitystring_to_array("I" * 50)
                a2.reference_id = 0
                a2.reference_start = 90
                a2.cigar = ((0, 50),)
                a2.mapping_quality = 60
                a2.set_tag("qs", 22.1, value_type="f")
                out_bam.write(a2)

                # Read 3: covers pos 100 with G (alt) on REVERSE strand, qs = 32.8 (int 32 -> >30)
                a3 = pysam.AlignedSegment()
                a3.query_name = "r3"
                seq3 = list("A" * 50)
                seq3[10] = "G"
                a3.query_sequence = "".join(seq3)
                a3.query_qualities = pysam.qualitystring_to_array("I" * 50)
                a3.reference_id = 0
                a3.reference_start = 90
                a3.cigar = ((0, 50),)
                a3.mapping_quality = 60
                a3.flag = 16  # reverse strand
                a3.set_tag("qs", 32.8, value_type="f")
                out_bam.write(a3)

                # Read 4: covers pos 100 with a deletion, qs = 22.0
                a4 = pysam.AlignedSegment()
                a4.query_name = "r4"
                a4.query_sequence = "A" * 49
                a4.query_qualities = pysam.qualitystring_to_array("I" * 49)
                a4.reference_id = 0
                a4.reference_start = 90
                # deletion of 1bp at reference pos 100 (offset 10)
                a4.cigar = ((0, 10), (2, 1), (0, 39))
                a4.mapping_quality = 60
                a4.set_tag("qs", 22.0, value_type="f")
                out_bam.write(a4)

            pysam.index(str(bam_path))

            chrom, by_pos = QS.collect_qs_pileup(bam_path, sites)
            self.assertEqual(chrom, "chr9")
            pos100 = by_pos[100]
            # Verify stratification and strand counts
            self.assertEqual(pos100[19]["A"], 1)
            self.assertEqual(pos100[19]["A_fwd"], 1)
            self.assertEqual(pos100[22]["G"], 1)
            self.assertEqual(pos100[22]["G_fwd"], 1)
            self.assertEqual(pos100[22]["del"], 1)
            self.assertEqual(pos100[32]["G"], 1)
            self.assertEqual(pos100[32]["G_rev"], 1)

            # Test bins
            bins = QS.make_bins(qs_min=20, qs_max=30, mode="binned")
            labels = [b["label"] for b in bins]
            self.assertIn("QS < 20", labels)
            self.assertIn("QS = 22", labels)
            self.assertIn("QS > 30", labels)

            # Test rows_for_bin
            bin_22 = [b for b in bins if b["label"] == "QS = 22"][0]
            rows_22 = QS.rows_for_bin(sites, by_pos, bin_22, "TEST_SAMPLE", strand="all")
            site100_r22 = [r for r in rows_22 if r["pos"] == "101"][0]
            self.assertEqual(site100_r22["depth"], "1")
            self.assertEqual(site100_r22["count_G"], "1")
            self.assertEqual(site100_r22["count_G_fwd"], "1")
            self.assertEqual(site100_r22["count_G_rev"], "0")
            self.assertEqual(site100_r22["count_del"], "1")
            self.assertEqual(site100_r22["vaf_percent"], "100.00")
            self.assertEqual(site100_r22["vaf_fwd_percent"], "100.00")
            self.assertEqual(site100_r22["vaf_rev_percent"], "0.00")

            # Test rows_for_bin with strand="fwd" and strand="rev"
            rows_22_fwd = QS.rows_for_bin(sites, by_pos, bin_22, "TEST_SAMPLE", strand="fwd")
            self.assertEqual(rows_22_fwd[0]["depth"], "1")
            rows_22_rev = QS.rows_for_bin(sites, by_pos, bin_22, "TEST_SAMPLE", strand="rev")
            self.assertEqual(rows_22_rev[0]["depth"], "0")

            # Test PDF and TSV generation
            out_pdf = tmp / "output_sweep.pdf"
            out_tsv = tmp / "output_sweep.tsv"

            argv = [
                "--bam", str(bam_path),
                "--bed", str(bed_path),
                "--output-pdf", str(out_pdf),
                "--output-tsv", str(out_tsv),
                "--qs-min", "20",
                "--qs-max", "25",
                "--mode", "both",
                "--split-strands",
                "--sample", "TEST_SAMPLE",
            ]
            parsed_args = QS.parse_args(argv)
            self.assertEqual(parsed_args.sample, "TEST_SAMPLE")
            self.assertEqual(parsed_args.strand, "both")

            # Run main workflow
            with patch.object(QS, "parse_args", return_value=parsed_args):
                QS.main()

            self.assertTrue(out_pdf.is_file())
            self.assertGreater(out_pdf.stat().st_size, 5000)
            # Verify it starts with PDF magic bytes
            with out_pdf.open("rb") as f:
                magic = f.read(4)
                self.assertEqual(magic, b"%PDF")

            self.assertTrue(out_tsv.is_file())
            tsv_lines = out_tsv.read_text().splitlines()
            self.assertGreater(len(tsv_lines), 5)
            self.assertIn("qs_bin", tsv_lines[0])
            self.assertIn("strand", tsv_lines[0])
            self.assertIn("vaf_fwd_percent", tsv_lines[0])
            self.assertIn("alt_minor_strand_ratio", tsv_lines[0])
            self.assertIn("count_del", tsv_lines[0])

    def test_main_summary_tracks_alleles_total_noise_and_valid_deletion_rates(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bed_path = tmp / "panel.bed"
            bed_path.write_text(
                "chr9\t100\t101\tNM:c.(949T>C)\tT\tC\tp.Phe317Leu\n"
                "chr9\t102\t103\tNM:c.(951C>G)\tC\tG\tp.Phe317Leu\n"
                "chr9\t104\t105\tNM:c.(763T>C)\tT\tC\tp.Glu255Lys\n"
                "chr9\t104\t105\tNM:c.(763T>C)\tT\tC\tp.Glu255LysDuplicate\n"
                "chr9\t106\t107\tNM:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
                "chr9\t200\t201\tNM:c.(757A>G)\tA\tG\tp.Tyr253Cys\n"
            )
            bam_path = tmp / "summary.bam"
            header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 1000, "SN": "chr9"}]}
            with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out_bam:
                for i in range(100):
                    seq = list("A" * 50)
                    seq[20] = "C" if i < 20 else "T"  # F317L, 20%
                    seq[22] = "G" if i < 10 else "C"  # same protein label, 10%
                    seq[24] = "A" if i < 10 else ("G" if i < 20 else "T")  # 20% total nonref noise
                    seq[26] = "C"
                    aln = pysam.AlignedSegment()
                    aln.query_name = f"r{i}"
                    if i < 90:
                        # Delete reference offset 26 (chr9:107).
                        seq.pop(26)
                        aln.cigar = ((0, 26), (2, 1), (0, 23))
                    else:
                        aln.cigar = ((0, 50),)
                    aln.query_sequence = "".join(seq)
                    aln.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
                    aln.reference_id = 0
                    aln.reference_start = 80
                    aln.mapping_quality = 60
                    aln.set_tag("qs", 25.5, value_type="f")
                    out_bam.write(aln)
            pysam.index(str(bam_path))

            args = SimpleNamespace(
                split_strands=False, bam=bam_path, bed=bed_path, sample="S",
                min_base_quality=0, min_mapping_quality=0, qs_min=20, qs_max=20,
                mode="cumulative", strand="all", output_pdf=tmp / "summary.pdf",
                output_tsv=None, zoom_max=5, color_mode="bases",
            )
            captured = []
            captured_rows = []

            def capture_summary(summaries, sample, **_):
                captured.extend(summaries)
                return QS.plt.figure()

            def capture_rows(**kwargs):
                captured_rows.extend(kwargs["rows"])
                return None

            with (
                patch.object(QS, "parse_args", return_value=args),
                patch.object(QS, "build_summary_figure", side_effect=capture_summary),
                patch.object(QS, "build_plot_figure", side_effect=capture_rows),
                patch.object(QS, "PdfPages"),
            ):
                QS.main()

            summary = captured[0]
            self.assertEqual(len(summary["mutation_vafs"]), 2)
            self.assertEqual(sorted(summary["mutation_vafs"].values()), [10.0, 20.0])
            self.assertIn("chr9:101 T>C", " ".join(summary["mutation_vafs"]))
            self.assertIn("chr9:103 C>G", " ".join(summary["mutation_vafs"]))
            self.assertAlmostEqual(summary["mean_background_noise"], 10.0)
            self.assertIsNone(summary["mean_background_noise_rev"])
            self.assertAlmostEqual(summary["mean_del_rate"], 22.5)
            uncovered = next(row for row in captured_rows if row["pos"] == "201")
            self.assertEqual(uncovered["status"], "LOW_DEPTH")
            self.assertEqual(uncovered["filter_reason"], "LOW_DEPTH")

            fig = QS.build_summary_figure([summary], "S")
            QS.plt.close(fig)

    def test_fractional_qs_uses_bounded_bins_and_exact_cumulative_thresholds(self):
        with tempfile.TemporaryDirectory() as directory:
            tmp = Path(directory)
            bam_path = tmp / "fractional.bam"
            header = {"HD": {"SO": "coordinate"}, "SQ": [{"SN": "chr9", "LN": 300}]}
            with pysam.AlignmentFile(str(bam_path), "wb", header=header) as output:
                for index in range(1000):
                    read = pysam.AlignedSegment()
                    read.query_name = f"read{index}"
                    read.query_sequence = "A" * 150
                    read.query_qualities = [30] * 150
                    read.reference_id = 0
                    read.reference_start = 0
                    read.mapping_quality = 60
                    read.cigartuples = [(0, 150)]
                    read.set_tag("qs", 20 + index / 500, "f")
                    output.write(read)
            pysam.index(str(bam_path))
            sites = [{"chrom": "chr9", "pos0": pos, "pos1": pos + 1,
                      "ref": "G", "alt": "A", "name": "NM:c.(763G>A)",
                      "mutation": "E255K", "aa_change": "p.Glu255Lys"} for pos in (49, 59)]
            _, counts = QS.collect_qs_pileup(bam_path, sites, cumulative_thresholds=[20, 21, 22])
            for pos in (49, 59):
                self.assertEqual(set(counts[pos]), {20, 21})
                self.assertEqual(counts[pos][20]["A"], 500)
                self.assertEqual(counts[pos][21]["A"], 500)
            bins = QS.make_bins(20, 22, "cumulative")
            depths = [int(QS.rows_for_bin(sites, counts, bin_def, "S")[0]["depth"]) for bin_def in bins]
            self.assertEqual(depths, [1000, 500, 0])

    def test_detect_genome_build_and_find_bed(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            # Create BAM with hg19 chr9 length
            bam19 = tmp / "hg19.bam"
            header19 = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 141213431, "SN": "chr9"}]}
            with pysam.AlignmentFile(str(bam19), "wb", header=header19) as f:
                pass
            self.assertEqual(QS.detect_genome_build_from_bam(bam19), "hg19")

            # Create BAM with hg38 chr9 length
            bam38 = tmp / "hg38.bam"
            header38 = {"HD": {"VN": "1.0"}, "SQ": [{"LN": 138394717, "SN": "chr9"}]}
            with pysam.AlignmentFile(str(bam38), "wb", header=header38) as f:
                pass
            self.assertEqual(QS.detect_genome_build_from_bam(bam38), "hg38")

            # Check find_bed_file resolves standard beds
            bed19 = QS.find_bed_file("hg19")
            self.assertIsNotNone(bed19)
            self.assertTrue(bed19.is_file())
            self.assertIn("hg19", bed19.name)



@unittest.skipIf(matplotlib is None, "matplotlib is required for the QS sweep and background plots")
class NamedMutationsTests(unittest.TestCase):
    """--mutations states the truth instead of guessing it from >= 1 %."""

    # Five panel positions. F317L appears twice: one amino-acid change, two
    # nucleotide changes, as in the real panel.
    BED = (
        "chr9\t100\t101\tNM:c.(758T>C)\tT\tC\tp.Tyr253His\n"
        "chr9\t102\t103\tNM:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
        "chr9\t104\t105\tNM:c.(949T>C)\tT\tC\tp.Phe317Leu\n"
        "chr9\t106\t107\tNM:c.(951C>G)\tC\tG\tp.Phe317Leu\n"
        "chr9\t108\t109\tNM:c.(1052T>C)\tT\tC\tp.Met351Thr\n"
    )
    REF = {100: "T", 102: "C", 104: "T", 106: "C", 108: "T"}

    def _sites(self, tmp):
        bed = Path(tmp) / "panel.bed"
        bed.write_text(self.BED)
        return bed, QS.load_branford_bed(bed)

    def test_names_resolve_in_every_accepted_form(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, sites = self._sites(tmp)
            y253h = {("chr9", 100, "T", "C")}
            for spelling in ("Y253H", "y253h", "p.Tyr253His", " Y253H "):
                self.assertEqual(QS.resolve_named_mutations(spelling, sites), y253h, spelling)
            both = QS.resolve_named_mutations("Y253H, T315I", sites)
            self.assertEqual(both, y253h | {("chr9", 102, "C", "T")})

    def test_a_name_covers_every_allele_that_produces_it(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, sites = self._sites(tmp)
            self.assertEqual(
                QS.resolve_named_mutations("F317L", sites),
                {("chr9", 104, "T", "C"), ("chr9", 106, "C", "G")},
            )

    def test_none_declares_a_negative_control(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, sites = self._sites(tmp)
            for spelling in ("none", "None", " NONE "):
                self.assertEqual(QS.resolve_named_mutations(spelling, sites), set(), spelling)

    def test_a_negative_control_keeps_every_position_in_the_background(self):
        """Without 'none', the sweep would guess -- here M351T, which is noise."""
        with tempfile.TemporaryDirectory() as tmp:
            control = self._run(tmp, "none")
        self.assertEqual(control["mutation_vafs"], {})

    def test_an_unknown_name_is_an_error_not_a_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            _, sites = self._sites(tmp)
            with self.assertRaises(ValueError) as raised:
                QS.resolve_named_mutations("Y253H,T315X", sites)
            self.assertIn("T315X", str(raised.exception))
            self.assertIn("T315I", str(raised.exception))  # and it lists what exists

    def _run(self, tmp, mutations):
        """Y253H is real at 3 %; M351T is a negative position that is noisy
        only in the poor reads, enough to pass 1 % in the unfiltered pool."""

        import random
        random.seed(7)
        bed, _ = self._sites(tmp)
        bam = Path(tmp) / "s.bam"
        header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 1000, "SN": "chr9"}]}
        with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
            for i in range(4000):
                poor = i < 400
                seq = ["A"] * 30
                for pos, base in self.REF.items():
                    seq[pos - 90] = base
                if i % 100 < 3:
                    seq[100 - 90] = "C"                       # Y253H, 3 %
                if random.random() < (0.15 if poor else 0.003):
                    seq[108 - 90] = "C"                       # M351T allele, noise only
                aln = pysam.AlignedSegment()
                aln.query_name = f"r{i}"
                aln.query_sequence = "".join(seq)
                aln.query_qualities = pysam.qualitystring_to_array("I" * 30)
                aln.reference_id, aln.reference_start, aln.mapping_quality = 0, 90, 60
                aln.cigar = ((0, 30),)
                aln.flag = 16 if i % 2 else 0
                aln.set_tag("qs", 8.0 if poor else 15.0, value_type="f")
                out.write(aln)
        pysam.index(str(bam))
        args = SimpleNamespace(
            split_strands=False, bam=bam, bed=bed, sample="S", min_base_quality=0,
            min_mapping_quality=0, qs_min=11, qs_max=11, mode="cumulative", strand="all",
            output_pdf=Path(tmp) / "o.pdf", output_tsv=None, zoom_max=5, color_mode="bases",
            mutations=mutations,
        )
        summaries = []
        with (
            patch.object(QS, "parse_args", return_value=args),
            patch.object(QS, "build_summary_figure",
                         side_effect=lambda s, n, **_: (summaries.extend(s), QS.plt.figure())[1]),
            patch.object(QS, "build_plot_figure", return_value=None),
            patch.object(QS, "PdfPages"),
        ):
            QS.main()
        return summaries[0]

    def test_named_mutations_keep_a_noisy_negative_in_the_background(self):
        with tempfile.TemporaryDirectory() as tmp:
            named = self._run(tmp, "Y253H")
        tracked = " ".join(named["mutation_vafs"])
        self.assertIn("Y253H", tracked)
        self.assertNotIn("M351T", tracked)

    def test_without_names_the_noisy_negative_is_taken_for_a_mutation(self):
        """The behaviour --mutations exists to avoid, kept as the default."""

        with tempfile.TemporaryDirectory() as tmp:
            guessed = self._run(tmp, None)
        with tempfile.TemporaryDirectory() as tmp:
            named = self._run(tmp, "Y253H")
        self.assertIn("M351T", " ".join(guessed["mutation_vafs"]))
        # Leaving the hotspot out of the background makes it look cleaner.
        self.assertLess(guessed["mean_background_noise"], named["mean_background_noise"])

    def test_a_typo_stops_the_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(SystemExit) as stopped:
                self._run(tmp, "Y253X")
        self.assertIn("Y253X", str(stopped.exception.code))



@unittest.skipIf(matplotlib is None, "matplotlib is required for the QS sweep and background plots")
class DepthMaskTests(unittest.TestCase):
    """Bins too thin to mean anything are shaded, not drawn as data."""

    # Three bins: 6x (QS 26 on a real run), 2,000x with a thin reverse strand,
    # and 9,000x. The 6x bin carries the kind of spike it produced for real.
    BINS = [
        {"short_label": "26", "median_depth": 6, "median_depth_fwd": 5, "median_depth_rev": 1,
         "mean_background_noise": 0.35, "mean_background_noise_fwd": 0.4, "mean_background_noise_rev": 0.0,
         "mean_del_rate": 0.1, "mean_ins_rate": 0.0,
         "mutation_vafs": {"T315I": 16.7}, "mutation_vafs_fwd": {"T315I": 20.0}, "mutation_vafs_rev": {"T315I": 0.0}},
        {"short_label": "11", "median_depth": 2000, "median_depth_fwd": 1600, "median_depth_rev": 400,
         "mean_background_noise": 1.1, "mean_background_noise_fwd": 1.2, "mean_background_noise_rev": 0.9,
         "mean_del_rate": 0.8, "mean_ins_rate": 0.3,
         "mutation_vafs": {"T315I": 3.6}, "mutation_vafs_fwd": {"T315I": 3.5}, "mutation_vafs_rev": {"T315I": 3.9}},
        {"short_label": "18", "median_depth": 9000, "median_depth_fwd": 7000, "median_depth_rev": 2000,
         "mean_background_noise": 0.33, "mean_background_noise_fwd": 0.35, "mean_background_noise_rev": 0.25,
         "mean_del_rate": 0.4, "mean_ins_rate": 0.1,
         "mutation_vafs": {"T315I": 3.9}, "mutation_vafs_fwd": {"T315I": 3.9}, "mutation_vafs_rev": {"T315I": 3.9}},
    ]

    @staticmethod
    def _line(ax, label_start):
        import math
        line = next(l for l in ax.get_lines() if l.get_label().startswith(label_start))
        return ["nan" if math.isnan(v) else round(float(v), 2) for v in line.get_ydata()]

    def test_thin_bins_are_not_drawn(self):
        fig = QS.build_summary_figure(self.BINS, "S", min_depth=500)
        try:
            _, ax2, ax3 = fig.axes[:3]
            self.assertEqual(self._line(ax2, "Mean Non-Ref Noise"), ["nan", 1.1, 0.33])
            self.assertEqual(self._line(ax3, "T315I (Combined)"), ["nan", 3.6, 3.9])
        finally:
            QS.plt.close(fig)

    def test_strand_curves_are_judged_on_their_own_depth(self):
        """The reverse strand has 400x at QS 11 although the bin has 2,000x."""
        fig = QS.build_summary_figure(self.BINS, "S", min_depth=500)
        try:
            ax2, ax3 = fig.axes[1], fig.axes[2]
            self.assertEqual(self._line(ax2, "Reverse Strand"), ["nan", "nan", 0.25])
            self.assertEqual(self._line(ax2, "Forward Strand"), ["nan", 1.2, 0.35])
            self.assertEqual(self._line(ax3, "T315I (Rev -)"), ["nan", "nan", 3.9])
        finally:
            QS.plt.close(fig)

    def test_the_depth_panel_is_always_drawn(self):
        fig = QS.build_summary_figure(self.BINS, "S", min_depth=500)
        try:
            self.assertEqual(self._line(fig.axes[0], "Median Depth (Combined)"), [6, 2000, 9000])
        finally:
            QS.plt.close(fig)

    def test_a_masked_spike_does_not_set_the_vaf_axis(self):
        fig = QS.build_summary_figure(self.BINS, "S", min_depth=500)
        try:
            self.assertLess(fig.axes[2].get_ylim()[1], 10)    # 16.7 % at 6x is not on the plot
        finally:
            QS.plt.close(fig)

    def test_masked_bins_are_shaded_and_labelled(self):
        fig = QS.build_summary_figure(self.BINS, "S", min_depth=500)
        try:
            labels = [h.get_label() for h in fig.axes[1].get_legend().legend_handles]
            self.assertTrue(any("500x" in label for label in labels), labels)
        finally:
            QS.plt.close(fig)

    def test_zero_draws_everything(self):
        fig = QS.build_summary_figure(self.BINS, "S")
        try:
            self.assertEqual(self._line(fig.axes[1], "Mean Non-Ref Noise"), [0.35, 1.1, 0.33])
        finally:
            QS.plt.close(fig)

    def test_every_bin_is_drawn_by_default(self):
        """Masking is opt-in: --min-plot-depth 500 to shade thin bins."""
        args = QS.parse_args(["--bam", "x.bam", "--output-pdf", "x.pdf"])
        self.assertEqual(args.min_plot_depth, 0.0)



@unittest.skipIf(matplotlib is None, "matplotlib is required for the QS sweep and background plots")
class SummaryPagePerModeTests(unittest.TestCase):
    """--mode both draws binned and cumulative on separate summary pages."""

    def _calls(self, mode):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            bed = tmp / "p.bed"
            bed.write_text("chr9\t100\t101\tNM:c.(944C>T)\tC\tT\tp.Thr315Ile\n")
            bam = tmp / "s.bam"
            header = {"HD": {"VN": "1.0", "SO": "coordinate"}, "SQ": [{"LN": 1000, "SN": "chr9"}]}
            with pysam.AlignmentFile(str(bam), "wb", header=header) as out:
                for i in range(200):
                    seq = list("A" * 30)
                    seq[10] = "C"
                    aln = pysam.AlignedSegment()
                    aln.query_name = f"r{i}"
                    aln.query_sequence = "".join(seq)
                    aln.query_qualities = pysam.qualitystring_to_array("I" * 30)
                    aln.reference_id, aln.reference_start, aln.mapping_quality = 0, 90, 60
                    aln.cigar = ((0, 30),)
                    aln.set_tag("qs", 14.0 + (i % 6), value_type="f")
                    out.write(aln)
            pysam.index(str(bam))
            args = SimpleNamespace(
                split_strands=False, bam=bam, bed=bed, sample="S", min_base_quality=0,
                min_mapping_quality=0, qs_min=14, qs_max=18, mode=mode, strand="all",
                output_pdf=tmp / "o.pdf", output_tsv=None, zoom_max=5, color_mode="bases",
            )
            calls = []

            def capture(summaries, sample, **kwargs):
                calls.append(([b["type"] for b in summaries], kwargs.get("mode")))
                return QS.plt.figure()

            with (
                patch.object(QS, "parse_args", return_value=args),
                patch.object(QS, "build_summary_figure", side_effect=capture),
                patch.object(QS, "build_plot_figure", return_value=None),
                patch.object(QS, "PdfPages"),
            ):
                QS.main()
            return calls

    def test_both_modes_get_a_page_each(self):
        calls = self._calls("both")
        self.assertEqual([mode for _, mode in calls], ["binned", "cumulative"])
        for types, mode in calls:
            self.assertEqual(set(types), {mode}, "a page must hold one kind of bin only")

    def test_a_single_mode_still_gets_one_page(self):
        for mode in ("binned", "cumulative"):
            calls = self._calls(mode)
            self.assertEqual(len(calls), 1, mode)
            self.assertEqual(calls[0][1], mode)

    def test_the_page_says_which_kind_it_is(self):
        bins = [{"short_label": ">=14", "type": "cumulative", "median_depth": 1000,
                 "mean_background_noise": 0.2, "mutation_vafs": {}}]
        fig = QS.build_summary_figure(bins, "S", mode="cumulative")
        try:
            self.assertIn("cumulative", fig._suptitle.get_text())
            self.assertIn("QS >= N", fig.axes[2].get_xlabel())
            for axis in (fig.axes[0], fig.axes[2]):
                self.assertIn("Thresholds", axis.get_title(loc="left"))
            # It used to give the >= 1 % rule's wording even for a negative control.
            self.assertIn("No mutations tracked", " ".join(t.get_text() for t in fig.axes[2].texts))
        finally:
            QS.plt.close(fig)


if __name__ == "__main__":
    unittest.main()

