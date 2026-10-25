#!/usr/bin/env python3
"""Tests for pipeline_start_scripts/plot_bcr_abl_technology_comparison.py."""

import tempfile
import unittest
from pathlib import Path

import pipeline_start_scripts.plot_bcr_abl_technology_comparison as plot_tool


class TestPlotBcrAblTechnologyComparison(unittest.TestCase):

    def test_normalize_aa(self):
        self.assertEqual(plot_tool.normalize_aa("p.Met351Thr"), "M351T")
        self.assertEqual(plot_tool.normalize_aa("T315I"), "T315I")
        self.assertEqual(plot_tool.normalize_aa(""), "")

    def test_parse_vaf(self):
        self.assertEqual(plot_tool.parse_vaf("96.97%"), 96.97)
        self.assertEqual(plot_tool.parse_vaf("95,3"), 95.3)
        self.assertIsNone(plot_tool.parse_vaf("-"))
        self.assertIsNone(plot_tool.parse_vaf("ses ej"))

    def test_end_to_end_figure_generation(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_dir = Path(tmp)
            truth_file = tmp_dir / "table6.txt"
            truth_file.write_text(
                "Tabell 6\t\tNGS Sequel\tNGS Sequel\tNGS Sequel\tNGS MiSeq\tNGS MiSeq\tNGS MiSeq\n"
                "Provnr\tSample ID\tratio\taa-variant\tallelfrekv %\ttranskript\taa-variant\tallelfrekv %\ttranskript\n"
                "SYN_01\tBCRABL-bc01\t5,3\tM351T\t97\te14-a2\tM351T\t95,3\te14-a2\n"
                "SYN_02\tBCRABL-bc06\t109\tT315I\t40,7\te13-a2\tT315I\t42,4\te13-a2\n"
                "SYN_03\tBCRABL-bc07\t2,3\tE459K\t1,5\te1-a2\tE459K\t1,7\te1-a2\n",
                encoding="utf-8"
            )

            results_file = tmp_dir / "results.tsv"
            results_file.write_text(
                "DNA Number\tSample ID\tType\tCaller\tGene\tChr\tPosition\tRef\tAlt\tVAF\tDepth\tAlt Reads\tConsequence / AA Change\n"
                "SYN_01\tBCRABL-bc01\tT\tBranford\tABL1\tchr9\t133748391\tT\tC\t97.02%\t10580\t10265\tM351T\n"
                "SYN_01\tBCRABL-bc01\tT\tVarDictRust\tABL1\tchr9\t133748391\tT\tC\t95.10%\t10616\t301,10099\t\n"
                "SYN_02\tBCRABL-bc06\tT\tBranford\tABL1\tchr9\t133748283\tC\tT\t43.73%\t29144\t12746\tT315I\n"
                "SYN_03\tBCRABL-bc07\tT\tBranford\tABL1\tchr9\t133753906\tG\tA\t2.01%\t44869\t901\tE459K\n",
                encoding="utf-8"
            )

            bed_file = tmp_dir / "branford.bed"
            bed_file.write_text(
                "chr9\t133748390\t133748391\t.\tT\tC\tp.Met351Thr\n"
                "chr9\t133748282\t133748283\t.\tC\tT\tp.Thr315Ile\n"
                "chr9\t133753905\t133753906\t.\tG\tA\tp.Glu459Lys\n",
                encoding="utf-8"
            )

            out_png = tmp_dir / "test_fig.png"
            rc = plot_tool.main([
                "-t", str(truth_file),
                "-v", str(results_file),
                "-b", str(bed_file),
                "-o", str(out_png),
            ])
            self.assertEqual(rc, 0)
            self.assertTrue(out_png.is_file())
            self.assertGreater(out_png.stat().st_size, 10000)


if __name__ == "__main__":
    unittest.main()
