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

            def capture_summary(summaries, sample):
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


if __name__ == "__main__":
    unittest.main()

