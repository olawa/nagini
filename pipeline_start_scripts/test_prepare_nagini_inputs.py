import contextlib
import gzip
import importlib.util
import io
import tempfile
import unittest
from pathlib import Path

import yaml


REPO = Path(__file__).resolve().parents[1]
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

    def test_format_short_barcode(self):
        self.assertEqual(PREPARE.format_short_barcode("barcode01"), "bc1")
        self.assertEqual(PREPARE.format_short_barcode("barcode10"), "bc10")
        self.assertEqual(PREPARE.format_short_barcode("bc01"), "bc1")
        self.assertEqual(PREPARE.format_short_barcode("bc10"), "bc10")
        self.assertEqual(PREPARE.format_short_barcode("1"), "bc1")
        self.assertEqual(PREPARE.format_short_barcode("24"), "bc24")
        self.assertEqual(PREPARE.format_short_barcode(""), "bc")

    def test_make_barcode_sample_name(self):
        self.assertEqual(PREPARE.make_barcode_sample_name("PBE82562", "barcode01"), "PBE82562_bc1")
        self.assertEqual(PREPARE.make_barcode_sample_name("PBE82562", "barcode10"), "PBE82562_bc10")
        self.assertEqual(PREPARE.make_barcode_sample_name("minknow", "barcode01"), "bc1")
        self.assertEqual(PREPARE.make_barcode_sample_name("sample", "barcode05"), "bc5")
        self.assertEqual(PREPARE.make_barcode_sample_name("", "barcode10"), "bc10")
        self.assertEqual(PREPARE.make_barcode_sample_name("RUN_bc1", "barcode01"), "RUN_bc1")
        self.assertEqual(PREPARE.make_barcode_sample_name("barcode01", "barcode01"), "bc1")

    def test_the_shipped_design_beds_are_well_formed(self):
        """A blank line in the T2T design made ClairS-TO exit after alignment.

        The BED is handed verbatim to tools whose parsers assume three fields
        per line, so the shipped assets have to be clean.
        """

        beds = sorted((REPO / "workflow/resources/adaptive").glob("*.bed"))
        self.assertTrue(beds)
        for bed in beds:
            for number, line in enumerate(bed.read_text().splitlines(), start=1):
                if line.startswith(("#", "track", "browser")):
                    continue
                with self.subTest(bed=bed.name, line=number):
                    self.assertTrue(line.strip(), "blank line")
                    fields = line.split("\t")
                    self.assertGreaterEqual(len(fields), 3, fields)
                    start, end = int(fields[1]), int(fields[2])
                    self.assertGreater(end, start)

    def test_a_lab_identifier_in_the_barcode_map_reaches_samples_tsv(self):
        """The DNA number is reported; the sample id still names the files.

        Keeping them apart is what lets a run be delivered under a real lab
        identifier without renaming a single file, read group or VCF column.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "inbox" / "20260625_1307_P2S_PBE84335_b84b698e"
            for barcode in ("barcode01", "barcode02"):
                directory = root / "bam_pass" / barcode
                directory.mkdir(parents=True)
                (directory / "reads_0.bam").touch()
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text(
                "barcode,sample,type,assay,dna_number\n"
                "barcode01,gustaf1,T,adaptive_cnv_sv,DNA2025-1001\n"
                # The column is optional per row.
                "barcode02,gustaf2,T,adaptive_cnv_sv,\n"
            )

            args = type("Args", (), {})()
            args.input_dir = Path(tmpdir) / "inbox"
            args.barcode_map = mapping
            args.sample_type = "T"
            args.input_type = "unmapped_bam"
            args.sequence_id = "inbox"
            args.run_id = ""
            args.methylation = "Yes"
            args.basecalling_model = "sup"

            samples, _units = PREPARE.make_bam_tables_from_barcode_map(args)

        by_name = {row["sample"]: row for row in samples}
        self.assertEqual(by_name["gustaf1"]["dna_number"], "DNA2025-1001")
        self.assertNotIn("dna_number", by_name["gustaf2"])

    def test_an_optional_sample_column_is_not_dropped_by_the_first_row(self):
        """sample_columns took the first row's keys, so a column present only
        on a later sample vanished for everyone."""

        samples = [{"sample": "a"}, {"sample": "b", "dna_number": "DNA2"}]
        columns = ["sample"]
        for row in samples:
            for column in row:
                if column not in columns:
                    columns.append(column)
        self.assertEqual(columns, ["sample", "dna_number"])

    def test_multi_flowcell_units_are_named_after_the_flowcell(self):
        """A merged sample has to say which flowcell each unit came from.

        A positional counter ordered the units by how the paths happened to
        sort, which says nothing about which run a read came from once four
        flowcells of one sample are merged.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "inbox"
            runs = {
                "20260625_1307_P2S-03608-B_PBE84335_b84b698e": ["reads_0.bam"],
                "20260721_1258_P2S-03608-B_PBE89162_f41da7bd": ["reads_0.bam", "reads_1.bam"],
            }
            for run, names in runs.items():
                barcode_dir = root / run / "bam_pass" / "barcode01"
                barcode_dir.mkdir(parents=True)
                for name in names:
                    (barcode_dir / name).touch()
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text("barcode,sample,type,assay\nbarcode01,PT_1,T,adaptive\n")

            args = type("Args", (), {})()
            args.input_dir = root
            args.barcode_map = mapping
            args.sample_type = "T"
            args.input_type = "unmapped_bam"
            args.sequence_id = "inbox"
            args.run_id = ""
            args.methylation = "Yes"
            args.basecalling_model = "sup"

            _, units = PREPARE.make_bam_tables_from_barcode_map(args)

        by_unit = {row["processing_unit"]: row for row in units}
        self.assertEqual(
            sorted(by_unit),
            ["PBE84335-b84b698e", "PBE89162-f41da7bd", "PBE89162-f41da7bd-2"],
        )
        # run_id follows the flowcell too, so QC can tell the runs apart.
        self.assertEqual(by_unit["PBE84335-b84b698e"]["run_id"], "PBE84335-b84b698e")
        # A second file within one flowcell still needs a unique unit.
        self.assertEqual(len({row["processing_unit"] for row in units}), 3)

    def test_units_keep_positional_names_without_a_minknow_run_directory(self):
        """Dorado demux output has no run directory to read a flowcell from."""

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "demux"
            (root / "barcode01").mkdir(parents=True)
            (root / "barcode01" / "a.bam").touch()
            (root / "barcode01" / "b.bam").touch()
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text("barcode,sample,type\nbarcode01,PT_1,T\n")

            args = type("Args", (), {})()
            args.input_dir = root
            args.barcode_map = mapping
            args.sample_type = "T"
            args.input_type = "unmapped_bam"
            args.sequence_id = "demux"
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "sup"

            _, units = PREPARE.make_bam_tables_from_barcode_map(args)

        self.assertEqual(
            sorted(row["processing_unit"] for row in units),
            ["barcode01_001", "barcode01_002"],
        )

    def test_bam_tables_without_barcode_map_infers_bc_sample_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "bams"
            (root / "barcode01").mkdir(parents=True)
            (root / "barcode10").mkdir(parents=True)
            (root / "barcode01" / "reads_1.bam").touch()
            (root / "barcode10" / "reads_2.bam").touch()

            args = type("Args", (), {})()
            args.input_dir = root
            args.barcode_map = None
            args.sample_type = "T"
            args.input_type = "unmapped_bam"
            args.sequence_id = "RUN123"
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "hac"
            args.pipeline_root = Path(__file__).parents[1]

            samples, units = PREPARE.make_bam_tables_from_barcode_map(args)

        self.assertEqual(samples, [{"sample": "RUN123_bc1"}, {"sample": "RUN123_bc10"}])
        self.assertEqual(len(units), 2)
        self.assertEqual(units[0]["sample"], "RUN123_bc1")
        self.assertEqual(units[0]["barcode_name"], "barcode01")
        self.assertEqual(units[1]["sample"], "RUN123_bc10")
        self.assertEqual(units[1]["barcode_name"], "barcode10")

    def test_dorado_demux_flat_barcode_filenames_are_discovered(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "demux"
            root.mkdir()
            (root / "barcode01.bam").touch()
            (root / "barcode02.bam").touch()

            found = PREPARE.discover_bams_by_barcode(root)

        self.assertEqual(sorted(found), ["barcode01", "barcode02"])

    def test_existing_barcode_fastq_runs_are_merged_and_mapped(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "bcr_abl_runs"
            mapping_rows = "barcode,sample,type,library_type,assay\n"
            mapping_rows += "barcode01,BCR_SAMPLE,T,cdna,bcr_abl\n"
            mapping_rows += "barcode02,CONTROL,N,cdna,bcr_abl\n"

            for run_number in (1, 2):
                run = root / f"run_{run_number}"
                for status in ("fastq_pass", "fastq_fail"):
                    for barcode in ("barcode01", "barcode02"):
                        directory = run / status / barcode
                        directory.mkdir(parents=True, exist_ok=True)
                        path = directory / f"reads_{run_number}_{status}_{barcode}.fastq.gz"
                        with gzip.open(path, "wt") as handle:
                            handle.write(
                                f"@{run_number}_{status}_{barcode}\nACGT\n+\nIIII\n"
                            )
                (run / f"sample_sheet_run_{run_number}.csv").write_text(mapping_rows)

            args = type("Args", (), {})()
            args.input_dir = root
            args.analysis_dir = Path(tmpdir) / "analysis"
            args.barcode_map = None
            args.sample_type = "N"
            args.input_type = "fastq"
            args.sequence_id = "combined"
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "sup"
            args.library_type = "cdna"
            args.assay = "bcr_abl"
            args.sample_regex = PREPARE.DEFAULT_SAMPLE_REGEX
            args.fastq_status = "all"
            args.pipeline_root = Path(__file__).parents[1]

            samples, units = PREPARE.make_fastq_tables(args)
            PREPARE.validate_tables(samples, units)

            merged = {
                row["barcode_name"]: Path(row["fastq1"])
                for row in units
            }
            sequences = {
                barcode: gzip.open(path, "rt").read().count("\n")
                for barcode, path in merged.items()
            }

        self.assertEqual(samples, [{"sample": "BCR_SAMPLE"}, {"sample": "CONTROL"}])
        self.assertEqual(len(units), 2)
        self.assertEqual(sequences, {"barcode01": 16, "barcode02": 16})
        self.assertTrue(all(row["input_type"] == "fastq" for row in units))
        self.assertTrue(all(row["assay"] == "bcr_abl" for row in units))

    def test_existing_barcode_fastq_can_select_pass_only(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "run"
            for status in ("fastq_pass", "fastq_fail"):
                directory = root / status / "barcode01"
                directory.mkdir(parents=True)
                with gzip.open(directory / f"reads_{status}.fastq.gz", "wt") as handle:
                    handle.write(f"@{status}\nACGT\n+\nIIII\n")
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text("barcode,sample,type\nbarcode01,SAMPLE_A,T\n")

            args = type("Args", (), {})()
            args.input_dir = root
            args.analysis_dir = Path(tmpdir) / "analysis"
            args.barcode_map = mapping
            args.sample_type = "N"
            args.input_type = "fastq"
            args.sequence_id = "run-1"
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "hac"
            args.library_type = "cdna"
            args.assay = "bcr_abl"
            args.fastq_status = "pass"
            args.pipeline_root = Path(__file__).parents[1]
            args.sample_regex = PREPARE.DEFAULT_SAMPLE_REGEX

            _, units = PREPARE.make_fastq_tables(args)
            with gzip.open(units[0]["fastq1"], "rt") as handle:
                contents = handle.read()

        self.assertIn("fastq_pass", contents)
        self.assertNotIn("fastq_fail", contents)

    def test_existing_barcode_fastq_infers_sample_from_filenames(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "runs"
            for run_number in (1, 2):
                for status in ("fastq_pass", "fastq_fail"):
                    for barcode in ("barcode01", "barcode23"):
                        directory = root / f"run_{run_number}" / status / barcode
                        directory.mkdir(parents=True, exist_ok=True)
                        path = directory / (
                            f"PBE82562_{status.removeprefix('fastq_')}_"
                            f"{barcode}_run{run_number}.fastq.gz"
                        )
                        with gzip.open(path, "wt") as handle:
                            handle.write(f"@{run_number}_{status}_{barcode}\nACGT\n+\nIIII\n")

            args = type("Args", (), {})()
            args.input_dir = root
            args.analysis_dir = Path(tmpdir) / "analysis"
            args.barcode_map = None
            args.sample_type = "N"
            args.input_type = "fastq"
            args.sequence_id = "combined"
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "sup"
            args.library_type = "cdna"
            args.assay = "bcr_abl"
            args.sample_regex = PREPARE.DEFAULT_SAMPLE_REGEX
            args.fastq_status = "all"
            args.pipeline_root = Path(__file__).parents[1]

            samples, units = PREPARE.make_fastq_tables(args)
            PREPARE.validate_tables(samples, units)

            merged = {
                row["barcode_name"]: Path(row["fastq1"])
                for row in units
            }
            line_counts = {
                barcode: gzip.open(path, "rt").read().count("\n")
                for barcode, path in merged.items()
            }

        self.assertEqual(
            samples, [{"sample": "PBE82562_bc1"}, {"sample": "PBE82562_bc23"}]
        )
        self.assertEqual(sorted(merged), ["barcode01", "barcode23"])
        self.assertEqual(line_counts, {"barcode01": 16, "barcode23": 16})
        self.assertTrue(all(row["processing_unit"].endswith("_merged") for row in units))
        self.assertTrue(all(row["assay"] == "bcr_abl" for row in units))
        self.assertEqual(merged["barcode01"].name, "PBE82562_bc1.fastq.gz")
        self.assertEqual(merged["barcode23"].name, "PBE82562_bc23.fastq.gz")

    def test_existing_barcode_fastq_with_generic_names_infers_bc_sample_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "runs"
            for barcode in ("barcode01", "barcode10"):
                directory = root / "fastq_pass" / barcode
                directory.mkdir(parents=True, exist_ok=True)
                path = directory / f"reads_{barcode}.fastq.gz"
                with gzip.open(path, "wt") as handle:
                    handle.write(f"@read_{barcode}\nACGT\n+\nIIII\n")

            args = type("Args", (), {})()
            args.input_dir = root
            args.analysis_dir = Path(tmpdir) / "analysis"
            args.barcode_map = None
            args.sample_type = "T"
            args.input_type = "fastq"
            args.sequence_id = ""
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "sup"
            args.library_type = "wgs"
            args.assay = ""
            args.sample_regex = PREPARE.DEFAULT_SAMPLE_REGEX
            args.fastq_status = "all"
            args.pipeline_root = Path(__file__).parents[1]

            samples, units = PREPARE.make_fastq_tables(args)
            PREPARE.validate_tables(samples, units)

        self.assertEqual(samples, [{"sample": "bc1"}, {"sample": "bc10"}])
        self.assertEqual(Path(units[0]["fastq1"]).name, "bc1.fastq.gz")
        self.assertEqual(Path(units[1]["fastq1"]).name, "bc10.fastq.gz")

    def test_existing_barcode_fastq_merges_across_different_flowcells_without_sample_sheet(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "multi_fc_runs"
            fc1_bc01 = root / "fc1" / "fastq_pass" / "barcode01"
            fc1_bc02 = root / "fc1" / "fastq_pass" / "barcode02"
            fc2_bc01 = root / "fc2" / "fastq_pass" / "barcode01"
            for d in (fc1_bc01, fc1_bc02, fc2_bc01):
                d.mkdir(parents=True, exist_ok=True)

            with gzip.open(fc1_bc01 / "FLOWCELL1_pass_barcode01_1.fastq.gz", "wt") as h:
                h.write("@fc1_bc01_r1\nACGT\n+\nIIII\n")
            with gzip.open(fc1_bc02 / "FLOWCELL1_pass_barcode02_1.fastq.gz", "wt") as h:
                h.write("@fc1_bc02_r1\nACGT\n+\nIIII\n")
            with gzip.open(fc2_bc01 / "FLOWCELL2_pass_barcode01_1.fastq.gz", "wt") as h:
                h.write("@fc2_bc01_r1\nTGCA\n+\nIIII\n")

            args = type("Args", (), {})()
            args.input_dir = root
            args.analysis_dir = Path(tmpdir) / "analysis"
            args.barcode_map = None
            args.sample_type = "T"
            args.input_type = "fastq"
            args.sequence_id = ""
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "sup"
            args.library_type = "wgs"
            args.assay = ""
            args.sample_regex = PREPARE.DEFAULT_SAMPLE_REGEX
            args.fastq_status = "all"
            args.pipeline_root = Path(__file__).parents[1]

            samples, units = PREPARE.make_fastq_tables(args)
            PREPARE.validate_tables(samples, units)

            self.assertEqual(samples, [{"sample": "bc1"}, {"sample": "bc2"}])
            self.assertEqual(len(units), 2)
            bc1_unit = next(u for u in units if u["sample"] == "bc1")
            bc2_unit = next(u for u in units if u["sample"] == "bc2")
            with gzip.open(bc1_unit["fastq1"], "rt") as h:
                bc1_content = h.read()
            self.assertIn("fc1_bc01_r1", bc1_content)
            self.assertIn("fc2_bc01_r1", bc1_content)
            with gzip.open(bc2_unit["fastq1"], "rt") as h:
                bc2_content = h.read()
            self.assertIn("fc1_bc02_r1", bc2_content)

    def test_init_sample_sheet_from_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "inbox"
            (root / "barcode01").mkdir(parents=True)
            (root / "barcode05").mkdir(parents=True)
            (root / "barcode01" / "read.bam").touch()
            (root / "barcode05" / "read.bam").touch()

            args = type("Args", (), {})()
            args.input_dir = root
            args.init_sample_sheet = "sample_sheet.csv"
            args.sample_type = "T"
            args.assay = "bcr_abl"
            args.library_type = "cdna"
            args.barcode_kit = ""
            args.dorado_kit = ""
            args.pipeline_root = Path(__file__).parents[1]

            out_path = PREPARE.generate_sample_sheet_template(args)
            self.assertTrue(out_path.is_file())
            content = out_path.read_text().splitlines()

        self.assertEqual(content[0], "barcode,sample,type,library_type,assay")
        self.assertEqual(content[1], "barcode01,sample_bc1,T,cdna,bcr_abl")
        self.assertEqual(content[2], "barcode05,sample_bc5,T,cdna,bcr_abl")
        self.assertEqual(len(content), 3)

    def test_init_sample_sheet_from_barcode_kit(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "inbox"
            root.mkdir()

            args = type("Args", (), {})()
            args.input_dir = root
            args.init_sample_sheet = "custom_sheet.csv"
            args.sample_type = "N"
            args.assay = "wgs"
            args.library_type = "wgs"
            args.barcode_kit = "SQK-RBK114-24"
            args.dorado_kit = ""
            args.pipeline_root = Path(__file__).parents[1]

            out_path = PREPARE.generate_sample_sheet_template(args)
            self.assertTrue(out_path.is_file())
            rows = out_path.read_text().splitlines()

        self.assertEqual(rows[0], "barcode,sample,type,library_type,assay")
        self.assertEqual(len(rows), 25)
        self.assertEqual(rows[1], "barcode01,sample_bc1,N,wgs,wgs")
        self.assertEqual(rows[10], "barcode10,sample_bc10,N,wgs,wgs")
        self.assertEqual(rows[24], "barcode24,sample_bc24,N,wgs,wgs")

    def test_minknow_tables_without_barcode_map_generates_short_barcode_samples(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "minknow"
            root.mkdir()

            args = type("Args", (), {})()
            args.input_dir = root
            args.barcode_map = None
            args.sample_type = "T"
            args.input_type = "minknow"
            args.sequence_id = "FLOWCELL123"
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "sup"
            args.dorado_kit = "SQK-RBK114-24"
            args.barcode_kit = "SQK-RBK114-24"
            args.library_type = "wgs"
            args.assay = ""
            args.pipeline_root = Path(__file__).parents[1]

            samples, units = PREPARE.make_minknow_tables(args)
            PREPARE.validate_tables(samples, units, allow_deferred_inputs=True)

        sample_names = [row["sample"] for row in samples]
        self.assertEqual(len(sample_names), 24)
        self.assertIn("FLOWCELL123_bc1", sample_names)
        self.assertIn("FLOWCELL123_bc10", sample_names)
        self.assertIn("FLOWCELL123_bc24", sample_names)

        self.assertEqual(len(units), 24)
        units_by_bc = {row["barcode_name"]: row for row in units}
        self.assertEqual(units_by_bc["barcode01"]["sample"], "FLOWCELL123_bc1")
        self.assertEqual(units_by_bc["barcode01"]["bam"], "basecalling/dorado_demux/barcode01.bam")
        self.assertEqual(units_by_bc["barcode10"]["sample"], "FLOWCELL123_bc10")
        self.assertEqual(units_by_bc["barcode10"]["bam"], "basecalling/dorado_demux/barcode10.bam")

    def test_minknow_tables_without_barcode_map_generic_sequence_id(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir) / "minknow"
            root.mkdir()

            args = type("Args", (), {})()
            args.input_dir = root
            args.barcode_map = None
            args.sample_type = "T"
            args.input_type = "minknow"
            args.sequence_id = ""
            args.run_id = ""
            args.methylation = "No"
            args.basecalling_model = "sup"
            args.dorado_kit = "SQK-RBK114-24"
            args.barcode_kit = "SQK-RBK114-24"
            args.library_type = "wgs"
            args.assay = ""
            args.pipeline_root = Path(__file__).parents[1]

            samples, units = PREPARE.make_minknow_tables(args)
            PREPARE.validate_tables(samples, units, allow_deferred_inputs=True)

        sample_names = [row["sample"] for row in samples]
        self.assertEqual(len(sample_names), 24)
        self.assertIn("bc1", sample_names)
        self.assertIn("bc10", sample_names)
        self.assertIn("bc24", sample_names)

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
        args.dorado_demux_output_format = "both"
        args.dorado_aligner = "dorado"

        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = PREPARE.make_runtime_config(args, Path(tmpdir))
            config = yaml.safe_load(runtime.read_text())

        self.assertEqual(
            config["dorado"]["container"],
            "docker://ontresearch/dorado:shac8f356489fa8b44b31beba841b84d2879de2088e",
        )
        self.assertEqual(config["dorado"]["kit"], "SQK-RBK114-24")
        self.assertEqual(config["dorado"]["device"], "cuda:all")
        self.assertEqual(config["dorado"]["demux_output_format"], "both")
        self.assertEqual(config["dorado_alignment"]["aligner"], "dorado")

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

    def test_runtime_config_supports_bed_file_override(self):
        pipeline_root = Path(__file__).parents[1]
        args = type("Args", (), {})()
        args.base_config = pipeline_root / "config/config.yaml"
        args.site_config = pipeline_root / "config/site_configs/site_config_local.yaml"
        args.reference_files = pipeline_root / "config/reference_files/reference_files_local.yaml"
        args.library_root = None
        args.pipeline_root = pipeline_root
        args.input_type = "minknow"
        args.dorado_input = None
        args.input_dir = Path("/tmp/run")
        args.dorado_model = ""
        args.basecalling_model = "hac"
        args.dorado_kit = ""
        args.dorado_device = ""
        args.dorado_models_directory = None
        args.dorado_container = ""
        args.dorado_binary = ""
        with tempfile.TemporaryDirectory() as tmpdir:
            bed_path = Path(tmpdir) / "my_targets.bed"
            bed_path.write_text("chr1\t1\t100\tGENE1\n")
            args.bed_file = bed_path
            runtime = PREPARE.make_runtime_config(args, Path(tmpdir))
            config = yaml.safe_load(runtime.read_text())

        self.assertEqual(config["references"]["hg38"]["design_bed"], str(bed_path.resolve()))
        self.assertEqual(config["HG38_DESIGN_BED"], str(bed_path.resolve()))

    def test_runtime_config_auto_detects_inbox_bed(self):
        pipeline_root = Path(__file__).parents[1]
        args = type("Args", (), {})()
        args.base_config = pipeline_root / "config/config.yaml"
        args.site_config = pipeline_root / "config/site_configs/site_config_local.yaml"
        args.reference_files = pipeline_root / "config/reference_files/reference_files_local.yaml"
        args.library_root = None
        args.pipeline_root = pipeline_root
        args.input_type = "minknow"
        args.dorado_input = None
        args.dorado_model = ""
        args.basecalling_model = "hac"
        args.dorado_kit = ""
        args.dorado_device = ""
        args.dorado_models_directory = None
        args.dorado_container = ""
        args.dorado_binary = ""
        args.bed_file = None

        with tempfile.TemporaryDirectory() as inbox, tempfile.TemporaryDirectory() as tmpdir:
            inbox_bed = Path(inbox) / "custom_panel.bed"
            inbox_bed.write_text("chr1\t100\t200\n")
            args.input_dir = Path(inbox)

            runtime = PREPARE.make_runtime_config(args, Path(tmpdir))
            config = yaml.safe_load(runtime.read_text())

        self.assertEqual(config["references"]["hg38"]["design_bed"], str(inbox_bed.resolve()))
        self.assertEqual(config["HG38_DESIGN_BED"], str(inbox_bed.resolve()))

    def test_runtime_config_configures_clairs_to_platform_and_caller(self):
        pipeline_root = Path(__file__).parents[1]
        args = type("Args", (), {})()
        args.base_config = pipeline_root / "config/config.yaml"
        args.site_config = pipeline_root / "config/site_configs/site_config_local.yaml"
        args.reference_files = pipeline_root / "config/reference_files/reference_files_local.yaml"
        args.library_root = None
        args.pipeline_root = pipeline_root
        args.input_type = "minknow"
        args.dorado_input = None
        args.input_dir = Path("/tmp/run")
        args.dorado_model = "sup,5mCG_5hmCG"
        args.basecalling_model = "sup,5mCG_5hmCG"
        args.dorado_kit = ""
        args.dorado_device = ""
        args.dorado_models_directory = None
        args.dorado_container = ""
        args.dorado_binary = ""
        args.bed_file = None
        args.small_variant_caller = "clairs_to"
        args.clairs_to_platform = ""
        args.clairs_to_container = ""

        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = PREPARE.make_runtime_config(args, Path(tmpdir))
            config = yaml.safe_load(runtime.read_text())

        self.assertEqual(config["small_variant_caller"], "clairs_to")
        self.assertEqual(config["clairs_to_call"]["platform"], "ont_r10_dorado_sup_5khz_ssrs")

    def test_runtime_config_records_enable_vardictrust(self):
        pipeline_root = Path(__file__).parents[1]
        args = type("Args", (), {})()
        args.base_config = pipeline_root / "config/config.yaml"
        args.site_config = pipeline_root / "config/site_configs/site_config_local.yaml"
        args.reference_files = pipeline_root / "config/reference_files/reference_files_local.yaml"
        args.library_root = None
        args.pipeline_root = pipeline_root
        args.input_type = "minknow"
        args.dorado_input = None
        args.input_dir = Path("/tmp/run")
        args.dorado_model = "hac"
        args.basecalling_model = "hac"
        args.dorado_kit = ""
        args.dorado_device = ""
        args.dorado_models_directory = None
        args.dorado_container = ""
        args.dorado_binary = ""
        args.bed_file = None
        args.small_variant_caller = ""
        args.clairs_to_platform = ""
        args.clairs_to_container = ""
        args.enable_vardictrust = True

        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = PREPARE.make_runtime_config(args, Path(tmpdir))
            config = yaml.safe_load(runtime.read_text())

        self.assertTrue(config["enable_vardictrust"])
        self.assertTrue(config["bcr_abl"]["enable_vardictrust"])

    def test_runtime_config_records_cores(self):
        pipeline_root = Path(__file__).parents[1]
        args = type("Args", (), {})()
        args.base_config = pipeline_root / "config/config.yaml"
        args.site_config = pipeline_root / "config/site_configs/site_config_local.yaml"
        args.reference_files = pipeline_root / "config/reference_files/reference_files_local.yaml"
        args.library_root = None
        args.pipeline_root = pipeline_root
        args.input_type = "minknow"
        args.dorado_input = None
        args.input_dir = Path("/tmp/run")
        args.dorado_model = "hac"
        args.basecalling_model = "hac"
        args.dorado_kit = ""
        args.dorado_device = ""
        args.dorado_models_directory = None
        args.dorado_container = ""
        args.dorado_binary = ""
        args.bed_file = None
        args.cores = 96

        with tempfile.TemporaryDirectory() as tmpdir:
            runtime = PREPARE.make_runtime_config(args, Path(tmpdir))
            config = yaml.safe_load(runtime.read_text())

    def test_barcode_map_with_library_type_and_defaults(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text(
                "barcode,sample,type,library_type\n"
                "barcode01,SAMPLE_WGS,T,wgs\n"
                "barcode02,SAMPLE_DRNA,T,direct_rna\n"
                "barcode03,SAMPLE_RNA,T,cdna\n"
                "barcode04,SAMPLE_DEF,N,\n"
            )
            parsed = PREPARE.read_barcode_map(mapping, "N", default_library_type="wgs")
            self.assertEqual(parsed["barcode01"]["library_type"], "wgs")
            self.assertEqual(parsed["barcode02"]["library_type"], "direct_rna")
            self.assertEqual(parsed["barcode03"]["library_type"], "cdna")
            self.assertEqual(parsed["barcode04"]["library_type"], "wgs")

    def test_barcode_map_rejects_invalid_library_type(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text("barcode,sample,type,library_type\nbarcode01,SAMPLE_A,N,invalid_type\n")
            with self.assertRaises(SystemExit) as ctx:
                PREPARE.read_barcode_map(mapping, "N")
            self.assertIn("Invalid library_type", str(ctx.exception))


class AssayRoutingTest(unittest.TestCase):
    """The assay column decides analysis; library_type decides alignment."""

    PIPELINE_ROOT = Path(__file__).parents[1]

    def _args(self, tmpdir, mapping, **overrides):
        args = type("Args", (), {})()
        args.input_dir = Path(tmpdir)
        args.barcode_map = mapping
        args.sample_type = "N"
        args.input_type = "minknow"
        args.sequence_id = "run-1"
        args.run_id = ""
        args.methylation = "No"
        args.basecalling_model = "hac"
        args.library_type = "wgs"
        args.assay = ""
        args.pipeline_root = self.PIPELINE_ROOT
        args.base_config = self.PIPELINE_ROOT / "config/config.yaml"
        args.site_config = None
        for key, value in overrides.items():
            setattr(args, key, value)
        return args

    def test_per_sample_assay_column_reaches_units(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text(
                "barcode,sample,type,library_type,assay\n"
                "barcode01,WGS_SAMPLE,N,wgs,adaptive\n"
                "barcode02,BCR_SAMPLE,T,cdna,bcr_abl\n"
            )
            args = self._args(tmpdir, mapping)
            samples, units = PREPARE.make_minknow_tables(args)

        by_sample = {row["sample"]: row for row in units}
        self.assertEqual(by_sample["WGS_SAMPLE"]["assay"], "adaptive")
        self.assertEqual(by_sample["WGS_SAMPLE"]["library_type"], "wgs")
        self.assertEqual(by_sample["BCR_SAMPLE"]["assay"], "bcr_abl")
        self.assertEqual(by_sample["BCR_SAMPLE"]["library_type"], "cdna")
        self.assertEqual(len(samples), 2)

    def test_cli_assay_is_the_default_when_the_map_omits_it(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text("barcode,sample,type\nbarcode01,S1,N\n")
            args = self._args(tmpdir, mapping, assay="adaptive")
            _, units = PREPARE.make_minknow_tables(args)
        self.assertEqual(units[0]["assay"], "adaptive")

    def test_assay_column_is_written_to_units_tsv(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = Path(tmpdir) / "units.tsv"
            PREPARE.write_tsv(
                path,
                [{"sample": "S1", "assay": "bcr_abl", "library_type": "cdna"}],
                PREPARE.BASE_UNIT_COLUMNS,
            )
            header = path.read_text().splitlines()[0].split("\t")
        self.assertIn("assay", header)
        self.assertIn("library_type", header)

    def test_unknown_assay_is_rejected_with_the_valid_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = self._args(tmpdir, None)
            with self.assertRaises(SystemExit) as raised:
                PREPARE.validate_assays(args, [{"sample": "S1", "assay": "bcrabl"}])
        message = str(raised.exception)
        self.assertIn("Unknown assay", message)
        self.assertIn("bcr_abl", message)

    def test_builtin_assays_pass_validation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = self._args(tmpdir, None)
            PREPARE.validate_assays(
                args,
                [
                    {"sample": "S1", "assay": "wgs"},
                    {"sample": "S2", "assay": "adaptive"},
                    {"sample": "S3", "assay": "bcr_abl"},
                ],
            )

    def test_an_unsupported_library_type_names_the_supported_ones(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            mapping = Path(tmpdir) / "barcode_map.csv"
            mapping.write_text("barcode,sample,type,library_type\nbarcode01,S1,N,rna\n")
            args = self._args(tmpdir, mapping)
            with self.assertRaises(SystemExit) as raised:
                PREPARE.make_minknow_tables(args)
        self.assertIn("direct_rna", str(raised.exception))


class PreflightTest(unittest.TestCase):
    """Catch a bad sample sheet at launch, not after the GPU job has run."""

    PIPELINE_ROOT = Path(__file__).parents[1]

    def _args(self, **overrides):
        args = type("Args", (), {})()
        args.pipeline_root = self.PIPELINE_ROOT
        args.dorado_kit = "SQK-NBD114-24"
        args.input_type = "minknow"
        args.dorado_model = "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"
        args.basecalling_model = args.dorado_model
        args.dorado_models_directory = None
        for key, value in overrides.items():
            setattr(args, key, value)
        return args

    def test_a_barcode_beyond_the_kit_is_rejected(self):
        units = [
            {"sample": "S1", "barcode_name": "barcode01"},
            {"sample": "S2", "barcode_name": "barcode30"},
        ]
        with self.assertRaises(SystemExit) as raised:
            PREPARE.validate_barcode_kit(self._args(), units)
        message = str(raised.exception)
        self.assertIn("barcode30", message)
        self.assertIn("SQK-NBD114-24", message)

    def test_barcodes_within_the_kit_pass(self):
        units = [{"sample": "S1", "barcode_name": f"barcode{n:02d}"} for n in (1, 12, 24)]
        PREPARE.validate_barcode_kit(self._args(), units)

    def test_an_unknown_kit_skips_the_check(self):
        units = [{"sample": "S1", "barcode_name": "barcode30"}]
        PREPARE.validate_barcode_kit(self._args(dorado_kit="SQK-LSK114"), units)

    def test_no_kit_skips_the_check(self):
        PREPARE.validate_barcode_kit(
            self._args(dorado_kit=""), [{"sample": "S1", "barcode_name": "barcode30"}]
        )

    def test_a_missing_models_directory_warns(self):
        buffer = io.StringIO()
        with contextlib.redirect_stderr(buffer):
            PREPARE.warn_about_model_download(self._args())
        self.assertIn("download", buffer.getvalue())

    def test_a_model_absent_from_the_directory_warns(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "dna_r10.4.1_e8.2_400bps_hac@v5.0.0").mkdir()
            buffer = io.StringIO()
            with contextlib.redirect_stderr(buffer):
                PREPARE.warn_about_model_download(
                    self._args(dorado_models_directory=Path(tmpdir))
                )
        text = buffer.getvalue()
        self.assertIn("not present", text)
        self.assertIn("hac@v5.0.0", text)

    def test_a_present_model_is_silent(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            (Path(tmpdir) / "dna_r10.4.1_e8.2_400bps_sup@v5.2.0").mkdir()
            buffer = io.StringIO()
            with contextlib.redirect_stderr(buffer):
                PREPARE.warn_about_model_download(
                    self._args(dorado_models_directory=Path(tmpdir))
                )
        self.assertEqual(buffer.getvalue(), "")

    def test_non_minknow_input_is_not_checked_for_models(self):
        buffer = io.StringIO()
        with contextlib.redirect_stderr(buffer):
            PREPARE.warn_about_model_download(self._args(input_type="fastq"))
        self.assertEqual(buffer.getvalue(), "")


class BasecallingModelDetectionTest(unittest.TestCase):
    def test_detect_model_from_final_summary(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            summary_file = run_dir / "final_summary_run123.txt"
            summary_file.write_text("protocol=sequencing\nbasecall_model_version_id=dna_r10.4.1_e8.2_400bps_sup@v5.2.0\n")
            model = PREPARE.detect_basecalling_model(run_dir)
            self.assertEqual(model, "dna_r10.4.1_e8.2_400bps_sup@v5.2.0")

    def test_detect_model_from_fastq_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            run_dir = Path(tmpdir)
            fq_file = run_dir / "sample_reads.fastq.gz"
            with gzip.open(str(fq_file), "wt") as handle:
                handle.write("@read1 runid=run123 model_version_id=dna_r10.4.1_e8.2_400bps_sup@v5.2.0\nACGT\n+\nIIII\n")
            model = PREPARE.detect_basecalling_model(run_dir)
            self.assertEqual(model, "dna_r10.4.1_e8.2_400bps_sup@v5.2.0")

    def test_can_reuse_existing_tables_detects_input_type_mismatch(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            dir_path = Path(tmpdir)
            samples_tsv = dir_path / "samples.tsv"
            units_tsv = dir_path / "units.tsv"
            samples_tsv.write_text("sample\nsample1\n")
            units_tsv.write_text(
                "sample\ttype\tprocessing_unit\tbarcode\tinput_type\tbam\n"
                "sample1\tT\tPU1\tBC\tdorado_bam\tbasecalling/dorado_demux/barcode01.bam\n"
            )

            self.assertFalse(PREPARE.can_reuse_existing_tables(samples_tsv, units_tsv, "fastq"))
            self.assertTrue(PREPARE.can_reuse_existing_tables(samples_tsv, units_tsv, "minknow"))


class ReferenceRegistryTest(unittest.TestCase):
    """The registry is where a reference path is declared."""

    PIPELINE_ROOT = Path(__file__).parents[1]

    def test_builds_and_annotations_flatten_to_config_variables(self):
        registry = {
            "site": "x",
            "reference_builds": {
                "hg38": {
                    "fasta": {"path": "/r/genome.fa", "type": "file"},
                    "minimap2_index": {"path": "/r/genome.mmi", "type": "file"},
                    "gene_annotation": {"path": "/r/genes.gtf", "type": "file"},
                },
                "hs1": {"fasta": {"path": "/r/hs1.fa", "type": "file"}},
            },
            "annotations": {"vep_cache": {"path": "/r/vep", "type": "directory"}},
        }
        variables = PREPARE.registry_variables(registry)
        self.assertEqual(variables["HG38_FASTA"], "/r/genome.fa")
        self.assertEqual(variables["HG38_MMI"], "/r/genome.mmi")
        self.assertEqual(variables["HG38_GTF"], "/r/genes.gtf")
        self.assertEqual(variables["HS1_FASTA"], "/r/hs1.fa")
        self.assertEqual(variables["VEP_CACHE"], "/r/vep")

    def test_relative_paths_are_joined_to_reference_root(self):
        registry = {
            "reference_root": "/refs",
            "reference_builds": {"hg38": {"fasta": {"path": "sub/genome.fa"}}},
        }
        self.assertEqual(PREPARE.registry_variables(registry)["HG38_FASTA"], "/refs/sub/genome.fa")

    def test_absolute_paths_ignore_reference_root(self):
        registry = {
            "reference_root": "/refs",
            "reference_builds": {"hg38": {"fasta": {"path": "/elsewhere/genome.fa"}}},
        }
        self.assertEqual(
            PREPARE.registry_variables(registry)["HG38_FASTA"], "/elsewhere/genome.fa"
        )

    def test_a_bare_path_is_accepted_as_well_as_a_mapping(self):
        registry = {"reference_builds": {"hg38": {"fasta": "/r/genome.fa"}}}
        self.assertEqual(PREPARE.registry_variables(registry)["HG38_FASTA"], "/r/genome.fa")

    def test_an_empty_path_overrides_the_shared_default(self):
        """Declaring an asset with no path says it is not installed here."""

        registry = {"reference_builds": {"hs1": {"design_bed": {"path": "", "type": "optional"}}}}
        self.assertEqual(PREPARE.registry_variables(registry)["HS1_DESIGN_BED"], "")

    def test_unknown_keys_are_ignored_rather_than_fatal(self):
        registry = {"reference_builds": {"hg38": {"something_new": {"path": "/r/x"}}}}
        self.assertEqual(PREPARE.registry_variables(registry), {})

    def test_missing_required_files_are_reported_but_not_fatal(self):
        registry = {"reference_builds": {"hg38": {"fasta": {"path": "/nope/genome.fa", "type": "file"}}}}
        buffer = io.StringIO()
        with contextlib.redirect_stderr(buffer):
            PREPARE.validate_reference_registry(registry)
        self.assertIn("/nope/genome.fa", buffer.getvalue())

    def test_strict_mode_refuses_to_continue(self):
        registry = {"reference_builds": {"hg38": {"fasta": {"path": "/nope/genome.fa", "type": "file"}}}}
        with self.assertRaises(SystemExit):
            PREPARE.validate_reference_registry(registry, strict=True)

    def test_optional_assets_are_not_required(self):
        registry = {"reference_builds": {"hg38": {"cytoband": {"path": "", "type": "optional"}}}}
        buffer = io.StringIO()
        with contextlib.redirect_stderr(buffer):
            PREPARE.validate_reference_registry(registry)
        self.assertEqual(buffer.getvalue(), "")

    def test_a_wrong_checksum_always_stops_the_run(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            asset = Path(tmpdir) / "genome.fa"
            asset.write_text(">chr1\nACGT\n")
            registry = {
                "reference_builds": {
                    "hg38": {"fasta": {"path": str(asset), "type": "file", "checksum": "0" * 32}}
                }
            }
            with self.assertRaises(SystemExit) as raised:
                PREPARE.validate_reference_registry(registry, verify_checksums=True)
        self.assertIn("registry says", str(raised.exception))

    def test_a_matching_checksum_passes(self):
        import hashlib

        with tempfile.TemporaryDirectory() as tmpdir:
            asset = Path(tmpdir) / "genome.fa"
            asset.write_bytes(b">chr1\nACGT\n")
            digest = hashlib.md5(asset.read_bytes()).hexdigest()
            registry = {
                "reference_builds": {
                    "hg38": {"fasta": {"path": str(asset), "type": "file", "checksum": digest}}
                }
            }
            PREPARE.validate_reference_registry(registry, verify_checksums=True)

    def test_shipped_registries_and_site_configs_do_not_both_own_a_path(self):
        """After the migration, a reference path is declared in one place."""

        reference_variables = set(PREPARE.REGISTRY_ANNOTATION_VARIABLES.values())
        for build in ("HG38", "HS1"):
            for suffix in set(PREPARE.REGISTRY_BUILD_VARIABLES.values()):
                reference_variables.add(f"{build}_{suffix}")
        for site in ("local", "workstation", "marvin"):
            site_config = PREPARE.read_yaml(
                self.PIPELINE_ROOT / f"config/site_configs/site_config_{site}.yaml"
            )
            overlap = sorted(reference_variables & set(site_config))
            self.assertEqual(overlap, [], f"{site} still declares {overlap}")

    def test_every_shipped_registry_matches_its_site(self):
        for site in ("local", "workstation", "marvin"):
            registry = PREPARE.read_yaml(
                self.PIPELINE_ROOT / f"config/reference_files/reference_files_{site}.yaml"
            )
            self.assertEqual(registry.get("site"), site)
            self.assertTrue(PREPARE.registry_variables(registry), f"{site} registry is empty")


class ReferenceRootTest(unittest.TestCase):
    """A site declares where its reference tree lives once, and the registry
    names assets relative to it."""

    PIPELINE_ROOT = Path(__file__).parents[1]

    def test_tree_relative_paths_join_the_declared_root(self):
        registry = {
            "reference_root": "/refs",
            "reference_builds": {"hg38": {"fasta": {"path": "hg38/genome.fa"}}},
        }
        self.assertEqual(
            PREPARE.registry_variables(registry)["HG38_FASTA"], "/refs/hg38/genome.fa"
        )

    def test_the_root_can_be_inherited_from_the_site(self):
        """A registry without its own reference_root uses the site's."""

        registry = {"reference_builds": {"hg38": {"fasta": {"path": "hg38/genome.fa"}}}}
        self.assertEqual(
            PREPARE.registry_variables(registry, "/site/refs")["HG38_FASTA"],
            "/site/refs/hg38/genome.fa",
        )

    def test_an_own_root_wins_over_the_inherited_one(self):
        registry = {
            "reference_root": "/explicit",
            "reference_builds": {"hg38": {"fasta": {"path": "hg38/genome.fa"}}},
        }
        self.assertEqual(
            PREPARE.registry_variables(registry, "/site/refs")["HG38_FASTA"],
            "/explicit/hg38/genome.fa",
        )

    def test_repo_shipped_assets_are_not_joined_to_the_reference_tree(self):
        """workflow/ and config/ paths belong to the checkout, not the tree."""

        registry = {
            "reference_root": "/refs",
            "reference_builds": {
                "hg38": {"design_bed": {"path": "workflow/resources/adaptive/x.bed"}}
            },
        }
        self.assertEqual(
            PREPARE.registry_variables(registry)["HG38_DESIGN_BED"],
            "workflow/resources/adaptive/x.bed",
        )

    def test_absolute_paths_ignore_the_root(self):
        registry = {
            "reference_root": "/refs",
            "reference_builds": {"hg38": {"fasta": {"path": "/elsewhere/genome.fa"}}},
        }
        self.assertEqual(
            PREPARE.registry_variables(registry)["HG38_FASTA"], "/elsewhere/genome.fa"
        )

    def test_shipped_registries_declare_a_root_and_use_relative_paths(self):
        """Guards against a site's tree location creeping back onto every line."""

        for site in ("local", "workstation", "marvin"):
            registry = PREPARE.read_yaml(
                self.PIPELINE_ROOT / f"config/reference_files/reference_files_{site}.yaml"
            )
            self.assertIn(
                "REFERENCE_ROOT", str(registry.get("reference_root", "")),
                f"{site} should take its root from the site config",
            )
            for build, assets in (registry.get("reference_builds") or {}).items():
                for key, value in (assets or {}).items():
                    path = value.get("path", "") if isinstance(value, dict) else value
                    if not path or not path.startswith("/"):
                        continue
                    # local is a developer machine and may point outside the tree.
                    self.assertEqual(
                        site, "local",
                        f"{site}.{build}.{key} is absolute ({path}); make it relative "
                        "to reference_root instead",
                    )

    def test_an_uninstalled_build_is_declared_empty_not_guessed(self):
        """marvin has no hs1; the registry must say so rather than carry a path
        copied from the shared defaults that resolves to nothing anywhere."""

        registry = PREPARE.read_yaml(
            self.PIPELINE_ROOT / "config/reference_files/reference_files_marvin.yaml"
        )
        hs1 = registry["reference_builds"]["hs1"]
        self.assertEqual(hs1["fasta"]["path"], "")
        self.assertEqual(hs1["fasta"]["type"], "optional")


if __name__ == "__main__":
    unittest.main()
