"""Exercise resource selection and report dependencies from the actual rules."""

import os
import sys
import unittest
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))


def common_section(title):
    """The helpers common.smk holds for one rule file, as source text."""

    source = (ROOT / "workflow/rules/common.smk").read_text()
    start = source.index(f"# --- {title} (")
    return source[start:source.index("\n# --- ", start + 1)]


class WorkflowContractTest(unittest.TestCase):
    def test_reference_aliases_never_silently_select_wrong_bundled_build(self):
        source = common_section("BCR-ABL1 resources")
        for configured, expected in (("hg19", "hg19"), ("GRCh37", "hg19"),
                                     ("b37", "hg19"), ("hg38", "hg38"), ("GRCh38", "hg38")):
            scope = {"os": os, "config": {"reference_build": configured},
                     "WorkflowError": ValueError, "workflow": SimpleNamespace(basedir="workflow")}
            exec(compile(source, "common.smk", "exec"), scope)
            self.assertTrue(scope["_branford_bed"]().endswith(f"branford_positions_{expected}.bed"))
            self.assertTrue(scope["_bcr_abl_exons_bed"]().endswith(f"bcr_abl1_exons_{expected}.bed"))
            self.assertTrue(scope["_bcr_abl_gtf"]().endswith(f"bcr_abl_{expected}.gtf"))
        scope["config"]["reference_build"] = "hs1"
        with self.assertRaisesRegex(ValueError, "No bundled BCR-ABL1 resources"):
            scope["_branford_bed"]()

    def test_all_sites_inputs_stay_aligned_in_a_mixed_assay_run(self):
        source = common_section("summary Excel report")
        scope = {
            "samples": ["wgs_sample", "bcr_sample"],
            "get_samples": lambda samples: samples,
            "get_unit_types_for_sample": lambda sample: ["T"],
            "clean_report_path": lambda value: str(value).strip(),
            "sample_needs_dna_variant_calling": lambda sample: False,
            "sample_has_assay_module": lambda sample, module: sample == "bcr_sample" and module == "bcr_abl",
        }
        for name in ("sample_has_target_metrics", "sniffles2_enabled", "sample_has_rindels_sv",
                     "savana_enabled", "sample_needs_cnv_calling",
                     "bcr_abl_vardictrust_enabled", "bcr_abl_rindels_enabled", "bcr_abl_bcftools_enabled"):
            scope[name] = lambda sample: False
        exec(compile(source, "common.smk", "exec"), scope)
        expected = "reports/bcr_abl/bcr_sample_T.branford_all_sites.tsv"
        self.assertEqual(scope["get_summary_report_paths"]("bcr_abl_all_sites"), ["", expected])
        self.assertEqual(scope["get_summary_report_inputs"](None)["bcr_abl_all_sites"], [expected])
        self.assertIn("bcr_abl_all_sites", scope["_SUMMARY_OPTIONAL_KEYS"])
        self.assertIn("bcr_abl_vardictrust_vcfs", scope["_SUMMARY_OPTIONAL_KEYS"])
        self.assertIn("bcr_abl_rindels_vcfs", scope["_SUMMARY_OPTIONAL_KEYS"])
        self.assertIn("bcr_abl_bcftools_vcfs", scope["_SUMMARY_OPTIONAL_KEYS"])

    def test_orthogonal_callers_rules_and_contract(self):
        import nagini_routing
        self.assertNotIn("enable_vardictrust", nagini_routing.CONFIG_FLAG_CONDITIONS)
        self.assertIn("enable_all_callers", nagini_routing.CONFIG_FLAG_CONDITIONS)
        self.assertIn("enable_rindels", nagini_routing.CONFIG_FLAG_CONDITIONS)
        self.assertIn("enable_bcftools_pileup", nagini_routing.CONFIG_FLAG_CONDITIONS)

        bcr_smk = (ROOT / "workflow/rules/bcr_abl.smk").read_text()
        self.assertNotIn("rule bcr_abl_vardictrust:", bcr_smk)
        self.assertIn("rule bcr_abl_bcftools_pileup:", bcr_smk)
        self.assertIn("rule bcr_abl_rindels:", bcr_smk)

        output_files = (ROOT / "config/output_files.yaml").read_text()
        self.assertNotIn("snv_indels/vardictrust/", output_files)
        self.assertIn("snv_indels/rindels/{sample}_{type}.rindels.vcf.gz", output_files)



class ReadFilterContractTest(unittest.TestCase):
    """The nanostream arguments each unit gets, from common.smk itself."""

    def scope(self, read_filter, rows):
        import shlex
        import shutil
        import subprocess
        import platform

        scope = {
            "os": os, "shlex": shlex, "shutil": shutil, "subprocess": subprocess,
            "platform": platform, "WorkflowError": ValueError,
            "workflow": SimpleNamespace(basedir="workflow"),
            "config": {"read_filter": read_filter},
            "get_assay_setting": lambda sample, section, key, default: default,
        }
        exec(compile(common_section("read filtering"), "common.smk", "exec"), scope)
        scope["get_unit"] = lambda wc: rows[wc.processing_unit]
        scope["get_unit_input_type"] = lambda wc: rows[wc.processing_unit]["input_type"]
        return scope

    @staticmethod
    def wc(unit):
        return SimpleNamespace(sample="S", type="N", processing_unit=unit, barcode="BC")

    ROWS = {
        "fq": {"sample": "S", "input_type": "fastq", "barcode_name": "barcode13",
               "sequencing_summary": "/run/input_summary/barcode13.sequencing_summary.tsv"},
        "bare": {"sample": "S", "input_type": "fastq", "barcode_name": "", "sequencing_summary": ""},
        "dorado": {"sample": "S", "input_type": "dorado_bam", "barcode_name": "barcode13",
                   "sequencing_summary": ""},
        "aligned": {"sample": "S", "input_type": "aligned_bam", "barcode_name": "barcode13",
                    "sequencing_summary": ""},
    }

    def test_both_ends_alone_selects_without_thresholds(self):
        s = self.scope({"enabled": False, "both_ends": True, "min_barcode_score": 70}, self.ROWS)
        self.assertEqual(s["get_read_filter_thresholds"](self.wc("fq")), "--qv 0")
        self.assertEqual(
            s["get_read_selection_args"](self.wc("fq")),
            "--sequencing-summary /run/input_summary/barcode13.sequencing_summary.tsv "
            "--barcode barcode13 --both-ends --min-barcode-score 70",
        )
        self.assertTrue(s["read_filter_step_enabled"]("S"))
        # Dorado demultiplexed its own reads, both ends included.
        self.assertEqual(s["get_read_selection_args"](self.wc("dorado")), "")
        with self.assertRaisesRegex(ValueError, "sequencing_summary"):
            s["get_read_selection_args"](self.wc("bare"))

    def test_filtering_uses_the_summary_for_quality_when_there_is_one(self):
        s = self.scope({"enabled": True, "min_qv": 14, "length_range": "800-2500"}, self.ROWS)
        self.assertEqual(s["get_read_filter_thresholds"](self.wc("fq")), "--qv 14 --len 800-2500")
        self.assertEqual(
            s["get_read_selection_args"](self.wc("fq")),
            "--sequencing-summary /run/input_summary/barcode13.sequencing_summary.tsv --barcode barcode13",
        )
        self.assertEqual(s["get_read_selection_args"](self.wc("bare")), "")

    def test_an_aligned_bam_is_realigned_only_when_reads_are_filtered_or_selected(self):
        aligned = self.ROWS["aligned"]
        self.assertFalse(self.scope({"enabled": False}, self.ROWS)["unit_is_realigned"](aligned))
        self.assertTrue(self.scope({"enabled": True}, self.ROWS)["unit_is_realigned"](aligned))
        self.assertTrue(self.scope({"enabled": False, "both_ends": True}, self.ROWS)["unit_is_realigned"](aligned))
        self.assertFalse(self.scope({"enabled": True}, self.ROWS)["unit_is_realigned"](self.ROWS["fq"]))

    def test_an_aligned_bam_handed_in_as_ubam_is_always_unaligned_first(self):
        rows = dict(self.ROWS, ubam_aligned={"sample": "S", "input_type": "unmapped_bam",
                                             "bam_aligned": "yes", "barcode_name": "", "sequencing_summary": ""})
        s = self.scope({"enabled": False}, rows)
        s["get_unit_input"] = lambda wc, column: "raw.bam"
        self.assertTrue(s["unit_bam_is_aligned"](rows["ubam_aligned"]))
        self.assertIn("reads/filtered/", s["get_unmapped_bam_for_alignment"](self.wc("ubam_aligned")))
        self.assertEqual(s["get_read_filter_thresholds"](self.wc("ubam_aligned")), "--qv 0")

    def test_the_bam_filter_always_writes_unaligned_reads(self):
        rules = (ROOT / "workflow/rules/read_filter.smk").read_text()
        bam_rule = rules[rules.index("rule read_filter_unmapped_bam:"):]
        self.assertIn("--unaligned", bam_rule)
        self.assertIn("{params.selection}", bam_rule)


if __name__ == "__main__":
    unittest.main()
