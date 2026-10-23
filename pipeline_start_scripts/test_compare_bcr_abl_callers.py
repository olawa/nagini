#!/usr/bin/env python3
"""Tests for pipeline_start_scripts/compare_bcr_abl_callers.py."""

import csv
import gzip
import os
import tempfile
import unittest
from pathlib import Path

import pipeline_start_scripts.compare_bcr_abl_callers as comp


class TestCompareBcrAblCallers(unittest.TestCase):

    def test_normalize_aa(self):
        self.assertEqual(comp.normalize_aa("p.Met351Thr"), "M351T")
        self.assertEqual(comp.normalize_aa("p.Thr315Ile"), "T315I")
        self.assertEqual(comp.normalize_aa("p.Tyr253His"), "Y253H")
        self.assertEqual(comp.normalize_aa("M351T"), "M351T")
        self.assertEqual(comp.normalize_aa("ingen variant"), "ingen variant")
        self.assertEqual(comp.normalize_aa(""), "")

    def test_parse_vaf_percent(self):
        self.assertEqual(comp.parse_vaf_percent("96.97%"), 96.97)
        self.assertEqual(comp.parse_vaf_percent("95,3"), 95.3)
        self.assertEqual(comp.parse_vaf_percent("0.946"), 94.6)
        self.assertEqual(comp.parse_vaf_percent("0.015"), 1.5)
        self.assertEqual(comp.parse_vaf_percent(0.963), 96.3)
        self.assertEqual(comp.parse_vaf_percent(42.4), 42.4)
        self.assertIsNone(comp.parse_vaf_percent("-"))
        self.assertIsNone(comp.parse_vaf_percent("ses ej"))
        self.assertIsNone(comp.parse_vaf_percent(""))
        self.assertIsNone(comp.parse_vaf_percent(None))

        # Test sub-1% VAFs with scale="percent" vs scale="fraction"
        self.assertEqual(comp.parse_vaf_percent("0,652", scale="percent"), 0.65)
        self.assertEqual(comp.parse_vaf_percent("0,6", scale="percent"), 0.6)
        self.assertEqual(comp.parse_vaf_percent("0.65", scale="percent"), 0.65)
        self.assertEqual(comp.parse_vaf_percent("0.65%", scale="percent"), 0.65)
        self.assertEqual(comp.parse_vaf_percent("0.65%", scale="fraction"), 0.65)
        self.assertEqual(comp.parse_vaf_percent("0.0065", scale="fraction"), 0.65)
        self.assertEqual(comp.parse_vaf_percent("0.65", scale="fraction"), 65.0)


    def test_parse_table6_real_content(self):
        sample_txt = (
            "Tabell 6\t\tNGS Sequel\tNGS Sequel\tNGS Sequel\tNGS MiSeq\tNGS MiSeq\tNGS MiSeq\n"
            "Provnr\tratio\taa-variant\tallelfrekv %\ttranskript\taa-variant\tallelfrekv %\ttranskript\n"
            "R22-01099\t5,3\tM351T\t97\te14-a2\tM351T\t95,3\te14-a2\n"
            "R21-01511\t102\tY253H/T315I\t99,9/100\te1-a2\tY253H/T315I\t99,9/97,9\te1-a2\n"
            "R21-02135\t64\tingen variant\t \te1-a2\tingen variant\t \te1-a2\n"
        )
        with tempfile.TemporaryDirectory() as tmp:
            p = Path(tmp) / "table6.txt"
            p.write_text(sample_txt, encoding="utf-8")
            records = comp.parse_table6(p)
            self.assertEqual(len(records), 4)

            # R22-01099
            r1 = records[0]
            self.assertEqual(r1["provnr"], "R22-01099")
            self.assertEqual(r1["variant"], "M351T")
            self.assertEqual(r1["sequel_af"], 97.0)
            self.assertEqual(r1["miseq_af"], 95.3)
            self.assertEqual(r1["miseq_iso"], "e14a2")

            # R21-01511 (split into Y253H and T315I)
            r2 = records[1]
            self.assertEqual(r2["provnr"], "R21-01511")
            self.assertEqual(r2["variant"], "Y253H")
            self.assertEqual(r2["sequel_af"], 99.9)
            self.assertEqual(r2["miseq_af"], 99.9)

            r3 = records[2]
            self.assertEqual(r3["provnr"], "R21-01511")
            self.assertEqual(r3["variant"], "T315I")
            self.assertEqual(r3["sequel_af"], 100.0)
            self.assertEqual(r3["miseq_af"], 97.9)

            # R21-02135
            r4 = records[3]
            self.assertEqual(r4["provnr"], "R21-02135")
            self.assertEqual(r4["variant"], "ingen variant")
            self.assertIsNone(r4["sequel_af"])

    def test_build_comparison_with_mock_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            run_dir = Path(tmp)
            (run_dir / "reports" / "bcr_abl").mkdir(parents=True)
            (run_dir / "snv_indels" / "vardictrust").mkdir(parents=True)

            # Mock samples.tsv
            samples_tsv = run_dir / "samples.tsv"
            samples_tsv.write_text("sample\tdna_number\nBCRABL-bc01\tR22-01099\n", encoding="utf-8")

            # Mock units.tsv
            units_tsv = run_dir / "units.tsv"
            units_tsv.write_text("sample\tbarcode\nBCRABL-bc01\tbarcode01\n", encoding="utf-8")

            # Mock Branford mutations TSV
            b_tsv = run_dir / "reports" / "bcr_abl" / "BCRABL-bc01_T.branford_mutations.tsv"
            b_tsv.write_text(
                "chrom\tpos\tref\talt\taa_change\tvaf_percent\tstatus\n"
                "chr9\t133748391\tT\tC\tM351T\t96.97\tDETECTED\n",
                encoding="utf-8"
            )

            # Mock VarDictRust VCF
            vd_vcf = run_dir / "snv_indels" / "vardictrust" / "BCRABL-bc01_T.bcr_abl.vcf"
            vd_vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr9\t133748391\t.\tT\tC\t100\tPASS\t.\tGT:DP:VD:AF\t0/1:12244:11587:0.946\n",
                encoding="utf-8"
            )

            # Mock Rindels VCF
            r_vcf = run_dir / "reports" / "bcr_abl" / "BCRABL-bc01_T.rindels.vcf"
            r_vcf.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n"
                "chr9\t133748391\t.\tT\tC\t100\tPASS\t.\tGT:DP:AD:AF\t0/1:12761:381,12293:0.963\n",
                encoding="utf-8"
            )

            # Mock BED
            bed = run_dir / "branford.bed"
            bed.write_text("chr9\t133748390\t133748391\t.\tT\tC\tp.Met351Thr\n", encoding="utf-8")
            aa_to_coords, coords_to_aa = comp.load_branford_bed(bed)

            truth_records = [{
                "provnr": "R22-01099",
                "ratio": "5,3",
                "variant": "M351T",
                "sequel_af": 97.0,
                "sequel_iso": "e14a2",
                "miseq_af": 95.3,
                "miseq_iso": "e14a2",
            }]

            rows = comp.build_comparison(
                truth_records=truth_records,
                run_dir=run_dir,
                aa_to_coords=aa_to_coords,
                coords_to_aa=coords_to_aa,
            )

            self.assertEqual(len(rows), 1)
            r = rows[0]
            self.assertEqual(r["dna_number"], "R22-01099")
            self.assertEqual(r["sample_id"], "BCRABL-bc01")
            self.assertEqual(r["barcode"], "barcode01")
            self.assertEqual(r["variant"], "M351T")
            self.assertEqual(r["miseq_af"], 95.3)
            self.assertEqual(r["sequel_af"], 97.0)
            self.assertEqual(r["branford_af"], 96.97)
            self.assertEqual(r["vardict_af"], 94.6)
            self.assertEqual(r["rindels_af"], 96.3)

            # Test write TSV and XLSX
            tsv_out = run_dir / "out.tsv"
            comp.write_tsv(rows, tsv_out)
            self.assertTrue(tsv_out.is_file())
            content = tsv_out.read_text(encoding="utf-8")
            self.assertIn("96.97%", content)
            self.assertIn("94.60%", content)
            self.assertIn("96.30%", content)

            try:
                import openpyxl
                has_openpyxl = True
            except ImportError:
                has_openpyxl = False

            if has_openpyxl:
                xlsx_out = run_dir / "out.xlsx"
                comp.write_excel(rows, xlsx_out)
                self.assertTrue(xlsx_out.is_file())


if __name__ == "__main__":
    unittest.main()
