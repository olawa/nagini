#!/usr/bin/env python3
"""Tests for pipeline_start_scripts/plot_bcr_abl_dilution_series.py."""

import tempfile
import unittest
from pathlib import Path

import pipeline_start_scripts.plot_bcr_abl_dilution_series as plot_tool


class TestPlotBcrAblDilutionSeries(unittest.TestCase):

    def test_parse_vaf(self):
        self.assertEqual(plot_tool.parse_vaf("11,10%"), 11.1)
        self.assertEqual(plot_tool.parse_vaf("9.859"), 9.859)
        self.assertIsNone(plot_tool.parse_vaf("n d"))
        self.assertIsNone(plot_tool.parse_vaf("-"))

    def test_parse_table7(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_file = Path(tmp) / "table7.txt"
            tmp_file.write_text(
                "Tabell 7\tT315I\tY253H\tT315I\tY253H\tT315I\tY253H\t\t\n"
                "ratio\t10%\t10%\t5%\t5%\t1%\t1%\t\t\n"
                "10%\t11,10%\t9,40%\t5,20%\t4,90%\t1,30%\tn d\tingen bakgrund\n"
                "0,50%\t10,30%\t9,50%\t4%\t3,30%\t1,20%\t1,50%\tingen bakgrund\n",
                encoding="utf-8"
            )
            res = plot_tool.parse_table7(tmp_file)
            self.assertEqual(res.get(("10%", 10.0, "T315I")), 11.1)
            self.assertEqual(res.get(("10%", 10.0, "Y253H")), 9.4)
            self.assertEqual(res.get(("10%", 5.0, "T315I")), 5.2)
            self.assertEqual(res.get(("0.5%", 1.0, "T315I")), 1.2)
            self.assertNotIn(("10%", 1.0, "Y253H"), res)  # was n d

    def test_build_dataset_and_generate_figure(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            out_png = tmp_dir / "dilution_test.png"

            miseq_mock = {
                ("10%", 10.0, "T315I"): 11.1,
                ("10%", 10.0, "Y253H"): 9.4,
                ("5%", 5.0, "T315I"): 5.0,
                ("0.5%", 1.0, "T315I"): 1.2,
                ("0.1%", 5.0, "T315I"): 18.4,
            }
            ont_mock = {
                ("10%", 10.0, "T315I"): 10.2,
                ("10%", 10.0, "Y253H"): 9.1,
                ("5%", 5.0, "T315I"): 4.8,
                ("0.5%", 1.0, "T315I"): 1.1,
                ("0.1%", 5.0, "T315I"): 16.5,
            }

            dataset = plot_tool.build_combined_dataset(miseq_mock, ont_mock)
            self.assertEqual(len(dataset), 5)

            plot_tool.create_two_panel_figure(dataset, out_png, dpi=100)
            self.assertTrue(out_png.is_file())
            self.assertGreater(out_png.stat().st_size, 1000)

            single_png = tmp_dir / "dilution_single.png"
            plot_tool.create_single_panel_figure(dataset, single_png, dpi=100)
            self.assertTrue(single_png.is_file())
            self.assertGreater(single_png.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
