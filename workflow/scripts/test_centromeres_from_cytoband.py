import importlib.util, io, tempfile, unittest
from pathlib import Path
SCRIPT = Path(__file__).with_name("centromeres_from_cytoband.py")
SPEC = importlib.util.spec_from_file_location("cent", SCRIPT)
CENT = importlib.util.module_from_spec(SPEC); SPEC.loader.exec_module(CENT)

class CentromeresFromCytobandTest(unittest.TestCase):
    """Rindels falls back to a built-in hg38 centromere set when given none."""

    CYTO = (
        "chr9\t0\t2200000\tp24.3\tgneg\n"
        "chr9\t43000000\t45500000\tp11.1\tacen\n"
        "chr9\t45500000\t47000000\tq11\tacen\n"
        "chr9\t47000000\t138394717\tq12\tgpos50\n"
        "chr20\t26400000\t28100000\tp11.1\tacen\n"
    )
    def test_acen_rows_merge_into_one_interval_per_chromosome(self):
        rows = CENT.centromeres(io.StringIO(self.CYTO))
        self.assertEqual(rows, [("chr9", 43000000, 47000000), ("chr20", 26400000, 28100000)])
    def test_non_centromere_bands_are_ignored(self):
        rows = CENT.centromeres(io.StringIO("chr1\t0\t100\tp36\tgneg\n"))
        self.assertEqual(rows, [])
    def test_a_missing_cytoband_yields_an_empty_bed_not_a_crash(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "c.bed"
            self.assertEqual(CENT.main(["--cytoband", str(Path(d) / "nope.bed"), "--output", str(out)]), 0)
            self.assertEqual(out.read_text(), "")
    def test_a_gzipped_cytoband_is_read(self):
        import gzip
        with tempfile.TemporaryDirectory() as d:
            src = Path(d) / "c.txt.gz"
            with gzip.open(src, "wt") as h: h.write(self.CYTO)
            out = Path(d) / "c.bed"
            CENT.main(["--cytoband", str(src), "--output", str(out)])
            self.assertIn("chr9\t43000000\t47000000\tcentromere", out.read_text())
class ScriptDirectiveCompatibilityTest(unittest.TestCase):
    def test_no_script_directive_file_has_a_future_import(self):
        """Snakemake prepends a preamble to a `script:` file.

        A `from __future__` import is then no longer the first statement and
        the script fails to parse -- which only shows up when the rule runs,
        after everything it depends on has already been computed.
        """

        import re

        rules = Path(__file__).resolve().parents[1] / "rules"
        scripts = set()
        for rule_file in rules.glob("*.smk"):
            scripts.update(re.findall(r'"\.\./scripts/([a-z_0-9]+\.py)"', rule_file.read_text()))
        self.assertTrue(scripts, "no script: directives found; the pattern may be stale")

        # Parsed, not grepped: a comment explaining the rule is not a violation
        # of it, and this file's own explanation would otherwise trip the test.
        import ast

        offenders = []
        for name in sorted(scripts):
            tree = ast.parse((Path(__file__).with_name(name)).read_text())
            if any(
                isinstance(node, ast.ImportFrom) and node.module == "__future__"
                for node in ast.walk(tree)
            ):
                offenders.append(name)
        self.assertEqual(offenders, [])


if __name__ == "__main__":
    unittest.main()
