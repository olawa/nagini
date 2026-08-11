import importlib.util
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).with_name("sv_summary.py")
SPEC = importlib.util.spec_from_file_location("sv_summary", SCRIPT)
SUMMARY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(SUMMARY)


class SvSummaryTest(unittest.TestCase):
    def test_cytoband_lookup_and_iscn_notation(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            cytoband = Path(tmpdir) / "cytoBand.txt"
            cytoband.write_text(
                "chr1\t0\t100\tp36.33\tgneg\n"
                "chr2\t0\t100\tp25.3\tgneg\n"
            )
            lookup = SUMMARY.CytobandLookup(str(cytoband))
            self.assertEqual(lookup.get_band("1", 10), "p36.33")
            self.assertEqual(
                SUMMARY.cytogenetic_notation(
                    "BND", "chr1", 10, 10, "chr2", 20, lookup
                ),
                "t(1;2)(p36.33;p25.3)",
            )


if __name__ == "__main__":
    unittest.main()
