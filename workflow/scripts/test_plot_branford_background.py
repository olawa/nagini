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


if __name__ == "__main__":
    unittest.main()
