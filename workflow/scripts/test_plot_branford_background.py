import tempfile
import unittest
from pathlib import Path

import plot_branford_background as PBB


class TestPlotBranfordBackground(unittest.TestCase):
    def test_load_all_sites_and_write_tsv_and_plot(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            all_sites_path = tmp / "test_sample_T.branford_all_sites.tsv"
            out_png = tmp / "test_plot.png"
            out_tsv = tmp / "test_counts.tsv"

            # Create mock all_sites TSV
            content = (
                "# sample=TEST_SAMPLE\n"
                "# min_read_qs=0\n"
                "record_type\tpanel_sites_expected\tsample\tchrom\tpos\tref\talt\tmutation\taa_change\tname\t"
                "consequence\tresistance_class\tcodon_context_status\tdepth\tad_ref\tad_alt\tvaf_percent\tstatus\tfilter_reason\t"
                "count_A\tcount_C\tcount_G\tcount_T\tpct_A\tpct_C\tpct_G\tpct_T\n"
                "SITE_AUDIT\t2\tTEST_SAMPLE\tchr9\t130862975\tG\tA\tE255K\tp.Glu255Lys\tNM_1\t"
                "MISSENSE\tCURATED_RESISTANCE\tSINGLE_SITE_ANNOTATION\t1000\t980\t15\t1.5\tDETECTED\tPASS\t"
                "15\t3\t980\t2\t1.5\t0.3\t98.0\t0.2\n"
                "SITE_AUDIT\t2\tTEST_SAMPLE\tchr9\t130872895\tC\tT\tT315I\tp.Thr315Ile\tNM_1\t"
                "MISSENSE\tCURATED_RESISTANCE\tSINGLE_SITE_ANNOTATION\t2000\t1900\t80\t4.0\tDETECTED\tPASS\t"
                "8\t1900\t12\t80\t0.4\t95.0\t0.6\t4.0\n"
            )
            all_sites_path.write_text(content)

            sample_name, rows = PBB.load_all_sites(all_sites_path)
            self.assertEqual(sample_name, "TEST_SAMPLE")
            self.assertEqual(len(rows), 2)

            PBB.write_summary_tsv(rows, out_tsv, sample_name)
            self.assertTrue(out_tsv.is_file())
            tsv_lines = out_tsv.read_text().splitlines()
            self.assertEqual(len(tsv_lines), 3)  # Header + 2 data rows
            self.assertIn("count_A\tcount_C\tcount_G\tcount_T", tsv_lines[0])
            self.assertIn("E255K", tsv_lines[1])
            self.assertIn("T315I", tsv_lines[2])

            # Generate plot
            PBB.generate_plot(
                rows=rows,
                output_png=out_png,
                sample=sample_name,
                zoom_max=5.0,
                color_mode="bases",
            )
            self.assertTrue(out_png.is_file())
            self.assertGreater(out_png.stat().st_size, 1000)

            # Also test relative mode
            out_png_rel = tmp / "test_plot_rel.png"
            PBB.generate_plot(
                rows=rows,
                output_png=out_png_rel,
                sample=sample_name,
                zoom_max=5.0,
                color_mode="relative",
            )
            self.assertTrue(out_png_rel.is_file())
            self.assertGreater(out_png_rel.stat().st_size, 1000)

    def test_load_rindels_pileup_counts_and_plot(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            pileup_path = tmp / "sample_test_T.pileup_counts.tsv"
            out_png = tmp / "test_rindels_plot.png"
            out_tsv = tmp / "test_rindels_counts.tsv"

            header = (
                "name\tchrom\tpos\tref\talt\traw_depth\tfiltered_depth\t"
                "A\tC\tG\tT\tA_fwd\tA_rev\tC_fwd\tC_rev\tG_fwd\tG_rev\tT_fwd\tT_rev\t"
                "INS\tINS_fwd\tINS_rev\tDEL\tDEL_fwd\tDEL_rev\tDEL_span\tS\tinsertions\tdeletions\tis_capped\tcap_value\n"
            )
            row1 = "p.Glu255Lys\tchr9\t130862976\tG\tA\t1000\t1000\t15\t3\t980\t2\t8\t7\t2\t1\t490\t490\t1\t1\t0\t0\t0\t0\t0\t0\t0\t0\t\t\tfalse\t0\n"
            row2 = "p.Thr315Ile\tchr9\t130872896\tC\tT\t2000\t2000\t8\t1900\t12\t80\t4\t4\t950\t950\t6\t6\t40\t40\t0\t0\t0\t0\t0\t0\t0\t0\t\t\tfalse\t0\n"
            pileup_path.write_text(header + row1 + row2)

            sample_name, rows = PBB.load_all_sites(pileup_path)
            self.assertEqual(sample_name, "sample_test_T")
            self.assertEqual(len(rows), 2)

            PBB.write_summary_tsv(rows, out_tsv, sample_name)
            self.assertTrue(out_tsv.is_file())
            lines = out_tsv.read_text().splitlines()
            self.assertEqual(len(lines), 3)
            self.assertIn("p.Glu255Lys", lines[1])
            self.assertIn("1000", lines[1])

            PBB.generate_plot(rows, out_png, sample=sample_name, include_depth=True, depth_position="top")
            self.assertTrue(out_png.is_file())
            self.assertGreater(out_png.stat().st_size, 1000)

            # Test bottom position
            out_png_bot = tmp / "test_rindels_bot.png"
            PBB.generate_plot(rows, out_png_bot, sample=sample_name, include_depth=True, depth_position="bottom")
            self.assertTrue(out_png_bot.is_file())
            self.assertGreater(out_png_bot.stat().st_size, 1000)

            # Test no depth panel
            out_png_no_depth = tmp / "test_rindels_no_depth.png"
            PBB.generate_plot(rows, out_png_no_depth, sample=sample_name, include_depth=False)
            self.assertTrue(out_png_no_depth.is_file())
            self.assertGreater(out_png_no_depth.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
