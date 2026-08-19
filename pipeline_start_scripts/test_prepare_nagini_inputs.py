import importlib.util
import tempfile
import unittest
from pathlib import Path

import yaml


SCRIPT = Path(__file__).with_name("prepare_nagini_inputs.py")
SPEC = importlib.util.spec_from_file_location("prepare_nagini_inputs", SCRIPT)
PREPARE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PREPARE)


class PrepareNaginiInputsTest(unittest.TestCase):
    def test_site_merge_and_placeholder_expansion(self):
        base = {"A": "base", "nested": {"keep": True, "replace": "base"}}
        site = {"B": "site", "nested": {"replace": "site"}}
        merged = PREPARE.deep_merge(base, site)
        self.assertEqual(merged["nested"], {"keep": True, "replace": "site"})

        expanded = PREPARE.expand_templates(
            {"path": "{{ROOT}}/{{NAME}}", "items": ["{{NAME}}"]},
            {"ROOT": "/tmp", "NAME": "asset"},
        )
        self.assertEqual(expanded, {"path": "/tmp/asset", "items": ["asset"]})

    def test_runtime_config_uses_selected_site(self):
        pipeline_root = Path(__file__).parents[1]
        args = type("Args", (), {})()
        args.base_config = pipeline_root / "config/config.yaml"
        args.site_config = pipeline_root / "config/site_configs/site_config_marvin.yaml"
        args.reference_files = pipeline_root / "config/reference_files/reference_files_marvin.yaml"
        args.library_root = None
        args.pipeline_root = pipeline_root

        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = PREPARE.make_runtime_config(args, Path(tmpdir))
            config = yaml.safe_load(runtime.read_text())

        self.assertEqual(config["site"], "marvin")
        self.assertEqual(config["resources"], str(pipeline_root / "config/resources_marvin.yaml"))
        self.assertTrue(config["references"]["hg38"]["fasta"].startswith("/beegfs-storage/"))
        self.assertTrue(config["annotate_small_variants_with_vep"])
        self.assertTrue(config["annotate_svs"])
        self.assertIn("--everything", config["vep"]["extra"])
        self.assertEqual(config["site_config_path"], str(args.site_config.resolve()))
        self.assertEqual(config["reference_files_path"], str(args.reference_files.resolve()))

    def test_workstation_reference_paths_are_preserved(self):
        pipeline_root = Path(__file__).parents[1]
        args = type("Args", (), {})()
        args.base_config = pipeline_root / "config/config.yaml"
        args.site_config = pipeline_root / "config/site_configs/site_config_workstation.yaml"
        args.reference_files = pipeline_root / "config/reference_files/reference_files_workstation.yaml"
        args.library_root = None
        args.pipeline_root = pipeline_root

        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = PREPARE.make_runtime_config(args, Path(tmpdir))
            config = yaml.safe_load(runtime.read_text())

        self.assertEqual(config["site"], "workstation")
        self.assertEqual(
            config["references"]["hg38"]["clinvar"],
            "/data/references/clinvar.vcf.gz",
        )
        self.assertEqual(
            config["references"]["hg38"]["mmi"],
            "/data/references/GRCh38_no_alt_analysis_set.mmi",
        )
        self.assertTrue(config["vep"]["skip"])

    def test_dorado_bam_pass_barcode_map_creates_hydra_tables(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "run" / "bam_pass"
            (root / "barcode01").mkdir(parents=True)
            (root / "barcode02").mkdir(parents=True)
            (root / "unclassified").mkdir(parents=True)
            bam_one = root / "barcode01" / "reads_a.bam"
            bam_two = root / "barcode02" / "reads_b.bam"
            bam_three = root / "barcode02" / "reads_c.bam"
            for bam in (bam_one, bam_two, bam_three, root / "unclassified" / "reads.bam"):
                bam.touch()
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text("barcode,sample,type\n01,SAMPLE_A,N\nbarcode02,SAMPLE_B,T\n")

            args = type("Args", (), {})()
            args.input_dir = root
            args.barcode_map = mapping
            args.sample_type = "N"
            args.input_type = "unmapped_bam"
            args.sequence_id = "run-1"
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "hac"

            samples, units = PREPARE.make_bam_tables_from_barcode_map(args)

        self.assertEqual(samples, [{"sample": "SAMPLE_A"}, {"sample": "SAMPLE_B"}])
        self.assertEqual(len(units), 3)
        self.assertEqual(units[0]["barcode"], "BC")
        self.assertEqual(units[0]["barcode_name"], "barcode01")
        self.assertEqual(units[0]["sample"], "SAMPLE_A")
        self.assertEqual(units[1]["sample"], "SAMPLE_B")
        self.assertEqual(units[1]["processing_unit"], "barcode02_001")
        self.assertTrue(all(row["input_type"] == "unmapped_bam" for row in units))

    def test_dorado_demux_flat_barcode_filenames_are_discovered(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "demux"
            root.mkdir()
            (root / "barcode01.bam").touch()
            (root / "barcode02.bam").touch()

            found = PREPARE.discover_bams_by_barcode(root)

        self.assertEqual(sorted(found), ["barcode01", "barcode02"])

    def test_minknow_tables_defer_dorado_outputs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "minknow"
            root.mkdir()
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text("barcode,sample,type\nbarcode05,SAMPLE_A,N\n")

            args = type("Args", (), {})()
            args.input_dir = root
            args.barcode_map = mapping
            args.sample_type = "N"
            args.input_type = "minknow"
            args.sequence_id = "run-1"
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "hac"

            samples, units = PREPARE.make_minknow_tables(args)
            PREPARE.validate_tables(samples, units, allow_deferred_inputs=True)

        self.assertEqual(samples, [{"sample": "SAMPLE_A"}])
        self.assertEqual(units[0]["input_type"], "dorado_bam")
        self.assertEqual(units[0]["bam"], "basecalling/dorado_demux/barcode05.bam")

    def test_minknow_runtime_config_keeps_dorado_container_settings(self):
        pipeline_root = Path(__file__).parents[1]
        args = type("Args", (), {})()
        args.base_config = pipeline_root / "config/config.yaml"
        args.site_config = pipeline_root / "config/site_configs/site_config_local.yaml"
        args.reference_files = pipeline_root / "config/reference_files/reference_files_local.yaml"
        args.library_root = None
        args.pipeline_root = pipeline_root
        args.input_type = "minknow"
        args.dorado_input = Path("/tmp/run/pod5")
        args.input_dir = args.dorado_input
        args.dorado_model = "hac"
        args.basecalling_model = "hac"
        args.dorado_kit = "SQK-RBK114-24"
        args.dorado_device = "cuda:all"
        args.dorado_models_directory = None

        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = PREPARE.make_runtime_config(args, Path(tmpdir))
            config = yaml.safe_load(runtime.read_text())

        self.assertEqual(
            config["dorado"]["container"],
            "docker://ontresearch/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e",
        )
        self.assertEqual(config["dorado"]["kit"], "SQK-RBK114-24")
        self.assertEqual(config["dorado"]["device"], "cuda:all")

    def test_minknow_runtime_config_supports_local_binary(self):
        pipeline_root = Path(__file__).parents[1]
        args = type("Args", (), {})()
        args.base_config = pipeline_root / "config/config.yaml"
        args.site_config = pipeline_root / "config/site_configs/site_config_local.yaml"
        args.reference_files = pipeline_root / "config/reference_files/reference_files_local.yaml"
        args.library_root = None
        args.pipeline_root = pipeline_root
        args.input_type = "minknow"
        args.dorado_input = Path("/tmp/run/pod5")
        args.input_dir = args.dorado_input
        args.dorado_model = "hac"
        args.basecalling_model = "hac"
        args.dorado_kit = "SQK-RBK114-24"
        args.dorado_device = "cuda:all"
        args.dorado_models_directory = None
        args.dorado_container = ""
        args.dorado_binary = "/usr/local/bin/dorado"

        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = PREPARE.make_runtime_config(args, Path(tmpdir))
            config = yaml.safe_load(runtime.read_text())

        self.assertEqual(config["dorado"]["binary"], "/usr/local/bin/dorado")
        self.assertEqual(config["dorado"]["container"], "none")


if __name__ == "__main__":
    unittest.main()
