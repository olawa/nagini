import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("rename_nagini_samples.py")
SPEC = importlib.util.spec_from_file_location("rename_nagini_samples", SCRIPT)
RENAME = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(RENAME)


def write_map(directory, text):
    path = Path(directory) / "ids.csv"
    path.write_text(text)
    return path


class IdMapTest(unittest.TestCase):
    def test_a_header_row_is_optional(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            with_header = write_map(tmpdir, "old_id,new_id\ngustaf1,DNA12345\n")
            self.assertEqual(RENAME.read_id_map(with_header), {"gustaf1": "DNA12345"})
            without = Path(tmpdir) / "plain.csv"
            without.write_text("gustaf1,DNA12345\n")
            self.assertEqual(RENAME.read_id_map(without), {"gustaf1": "DNA12345"})

    def test_the_genome_is_appended_to_every_new_id(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_map(tmpdir, "gustaf1,DNA12345\ngustaf2,DNA12346\n")
            self.assertEqual(
                RENAME.read_id_map(path, genome="hg38"),
                {"gustaf1": "DNA12345_hg38", "gustaf2": "DNA12346_hg38"},
            )

    def test_two_samples_cannot_become_the_same_id(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_map(tmpdir, "gustaf1,DNA12345\ngustaf2,DNA12345\n")
            with self.assertRaises(SystemExit):
                RENAME.read_id_map(path)

    def test_an_id_that_is_both_old_and_new_is_rejected(self):
        """Swapping two ids depends on the order the rewrite happens in.

        gustaf1 -> gustaf2 followed by gustaf2 -> DNA1 would put both samples
        under one name, silently.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            path = write_map(tmpdir, "gustaf1,gustaf2\ngustaf2,DNA12345\n")
            with self.assertRaises(SystemExit):
                RENAME.read_id_map(path)

    def test_a_new_id_that_would_break_a_path_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            for bad in ("DNA 12345", "sub/dir", ""):
                path = write_map(tmpdir, f"gustaf1,{bad}\n")
                with self.assertRaises(SystemExit):
                    RENAME.read_id_map(path)


class SubstitutionTest(unittest.TestCase):
    def setUp(self):
        self.mapping = {"gustaf1": "DNA12345", "gustaf10": "DNA99999"}
        self.pattern = RENAME.build_pattern(self.mapping)

    def test_a_shorter_id_does_not_match_inside_a_longer_one(self):
        """gustaf1 must not rewrite the first seven characters of gustaf10."""

        self.assertEqual(
            RENAME.substitute("gustaf10_T.bam", self.pattern, self.mapping),
            "DNA99999_T.bam",
        )
        self.assertEqual(
            RENAME.substitute("gustaf1_T.bam", self.pattern, self.mapping),
            "DNA12345_T.bam",
        )

    def test_the_id_is_rewritten_wherever_it_appears_as_a_token(self):
        text = (
            "alignment/merged/gustaf1_T.bam\n"
            "@RG\tID:gustaf1_T_PBE84335-b84b698e_BC\tSM:gustaf1_T\n"
            "sample\ttype\ngustaf1\tT\n"
        )
        updated = RENAME.substitute(text, self.pattern, self.mapping)
        self.assertNotIn("gustaf1", updated)
        self.assertEqual(updated.count("DNA12345"), 4)


class AnalysisDirectoryTest(unittest.TestCase):
    """End to end over a directory shaped like a finished run."""

    def _build(self, root):
        root = Path(root)
        (root / "alignment/merged").mkdir(parents=True)
        (root / "reports/sv_summary").mkdir(parents=True)
        (root / "qc/cramino").mkdir(parents=True)
        (root / "cnv/rindels").mkdir(parents=True)

        (root / "samples.tsv").write_text("sample\ttumor_content\ngustaf1\t1.0\n")
        (root / "units.tsv").write_text(
            "sample\ttype\tprocessing_unit\ngustaf1\tT\tPBE84335-b84b698e\n"
        )
        (root / "alignment/merged/gustaf1_T.bam").write_bytes(b"BAM\x01")
        (root / "alignment/merged/gustaf1_T.bam.bai").write_bytes(b"BAI")
        (root / "qc/cramino/gustaf1_T.txt").write_text("File name\tgustaf1_T.bam\n")
        (root / "reports/sv_summary/gustaf1_T.tsv").write_text(
            "sample\ttype\tcaller\ngustaf1\tT\tSniffles2\n"
        )
        (root / "reports/gustaf1_T.nagini_report.xlsx").write_bytes(b"PK\x03\x04gustaf1")
        with gzip.open(root / "cnv/rindels/gustaf1_T.vcf.gz", "wt") as handle:
            handle.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tgustaf1\n")
        (root / "cnv/rindels/gustaf1_T.vcf.gz.tbi").write_bytes(b"TBI")
        return root

    def _run(self, root, mapping, **kwargs):
        pattern = RENAME.build_pattern(mapping)
        for path in sorted(Path(root).rglob("*")):
            if not path.is_file():
                continue
            kind = RENAME.classify(path)
            if kind == "text":
                RENAME.rewrite_text(path, pattern, mapping, True)
            elif kind == "gzip":
                RENAME.rewrite_gzip(path, pattern, mapping, True)
        return RENAME.rename_paths(Path(root), pattern, mapping, True)

    def test_paths_and_contents_are_renamed_together(self):
        mapping = {"gustaf1": "DNA12345_hg38"}
        with tempfile.TemporaryDirectory() as tmpdir:
            root = self._build(tmpdir)
            self._run(root, mapping)

            self.assertTrue((root / "alignment/merged/DNA12345_hg38_T.bam").exists())
            self.assertTrue((root / "alignment/merged/DNA12345_hg38_T.bam.bai").exists())
            self.assertFalse((root / "alignment/merged/gustaf1_T.bam").exists())
            self.assertIn("DNA12345_hg38", (root / "samples.tsv").read_text())
            self.assertIn(
                "DNA12345_hg38\tT\tSniffles2",
                (root / "reports/sv_summary/DNA12345_hg38_T.tsv").read_text(),
            )
            with gzip.open(root / "cnv/rindels/DNA12345_hg38_T.vcf.gz", "rt") as handle:
                self.assertIn("DNA12345_hg38", handle.read())

    def test_a_stale_index_is_removed_when_its_file_is_rewritten(self):
        """Rewriting a BGZF file moves every offset the .tbi records."""

        mapping = {"gustaf1": "DNA12345"}
        with tempfile.TemporaryDirectory() as tmpdir:
            root = self._build(tmpdir)
            self._run(root, mapping)
            # Either re-created by tabix or removed, but never left stale.
            index = root / "cnv/rindels/DNA12345_T.vcf.gz.tbi"
            if index.exists():
                self.assertNotEqual(index.read_bytes(), b"TBI")

    def test_binary_files_are_not_rewritten_as_text(self):
        mapping = {"gustaf1": "DNA12345"}
        with tempfile.TemporaryDirectory() as tmpdir:
            root = self._build(tmpdir)
            self._run(root, mapping)
            # The BAM is renamed but its bytes are untouched without --bams.
            self.assertEqual(
                (root / "alignment/merged/DNA12345_T.bam").read_bytes(), b"BAM\x01"
            )
            self.assertEqual(RENAME.classify(root / "alignment/merged/DNA12345_T.bam.bai"), "skip")

    def test_an_existing_target_path_stops_the_rename(self):
        mapping = {"gustaf1": "gustaf2"}
        with tempfile.TemporaryDirectory() as tmpdir:
            root = self._build(tmpdir)
            (root / "alignment/merged/gustaf2_T.bam").write_bytes(b"other sample")
            with self.assertRaises(SystemExit):
                RENAME.rename_paths(
                    root, RENAME.build_pattern(mapping), mapping, apply_changes=True
                )


class GenomeFromConfigTest(unittest.TestCase):
    def test_the_build_is_read_from_the_run_config(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            (root / "config.nagini.yaml").write_text("reference:\n  build: hs1\n")
            self.assertEqual(RENAME.genome_from_config(root), "hs1")

    def test_a_missing_config_is_not_an_error_here(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            self.assertEqual(RENAME.genome_from_config(Path(tmpdir)), "")


if __name__ == "__main__":
    unittest.main()
