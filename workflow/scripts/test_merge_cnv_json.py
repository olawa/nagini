import gzip
import io
import os
import tempfile
import unittest

from merge_cnv_json import (
    annotation_parser,
    normalize_chrom,
    parse_cytobands,
    parse_fai,
)


class TestMergeCnvJson(unittest.TestCase):

    def setUp(self):
        self.colors = {
            "gneg": "#e3e3e3",
            "gpos25": "#555555",
            "gpos50": "#393939",
            "acen": "#963232",
        }

    def test_normalize_chrom(self):
        self.assertEqual(normalize_chrom("1"), "chr1")
        self.assertEqual(normalize_chrom("chr1"), "chr1")
        self.assertEqual(normalize_chrom("chrchr1"), "chr1")
        self.assertEqual(normalize_chrom("chrX"), "chrX")
        self.assertEqual(normalize_chrom("chrM"), "chrM")

    def test_parse_cytobands_standard_5_col(self):
        content = (
            "chr1\t0\t2300000\tp36.33\tgneg\n"
            "chr1\t2300000\t5300000\tp36.32\tgpos25\n"
            "chr1\t121500000\t125000000\tp11.1\tacen\n"
            "chr1\t125000000\t128900000\tq11\tacen\n"
        )
        with tempfile.NamedTemporaryFile("w", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            cytobands = parse_cytobands(temp_path, self.colors)
            self.assertIn("chr1", cytobands)
            self.assertEqual(len(cytobands["chr1"]), 4)
            # Centromeres: p11.1 right, q11 left
            self.assertEqual(cytobands["chr1"][2]["direction"], "right")
            self.assertEqual(cytobands["chr1"][3]["direction"], "left")
            self.assertEqual(cytobands["chr1"][0]["color"], "#e3e3e3")
            self.assertEqual(cytobands["chr1"][1]["color"], "#555555")
        finally:
            os.unlink(temp_path)

    def test_parse_cytobands_4_col_and_headers(self):
        content = (
            "#chrom\tchromStart\tchromEnd\tname\n"
            "track name=cytobands\n"
            "chr1\t0\t2300000\tp36.33\n"
            "chr1 2300000 5300000 p36.32\n"
            "\n"
            "chr2\t0\t1000\tp1\tgpos50\n"
        )
        with tempfile.NamedTemporaryFile("w", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            cytobands = parse_cytobands(temp_path, self.colors)
            self.assertIn("chr1", cytobands)
            self.assertEqual(len(cytobands["chr1"]), 2)
            # Default giemsa is "gneg"
            self.assertEqual(cytobands["chr1"][0]["giemsa"], "gneg")
            self.assertEqual(cytobands["chr1"][0]["color"], "#e3e3e3")
            self.assertEqual(cytobands["chr1"][1]["giemsa"], "gneg")

            self.assertIn("chr2", cytobands)
            self.assertEqual(cytobands["chr2"][0]["giemsa"], "gpos50")
            self.assertEqual(cytobands["chr2"][0]["color"], "#393939")
        finally:
            os.unlink(temp_path)

    def test_parse_cytobands_gzipped(self):
        content = b"chr1\t0\t1000\tp1\tgneg\n"
        with tempfile.NamedTemporaryFile("wb", suffix=".gz", delete=False) as f:
            f.write(gzip.compress(content))
            temp_path = f.name

        try:
            cytobands = parse_cytobands(temp_path, self.colors)
            self.assertIn("chr1", cytobands)
            self.assertEqual(len(cytobands["chr1"]), 1)
            self.assertEqual(cytobands["chr1"][0]["name"], "p1")
        finally:
            os.unlink(temp_path)

    def test_parse_cytobands_single_centromere_does_not_crash(self):
        content = (
            "chr1\t0\t1000\tp1\tgneg\n"
            "chr1\t1000\t2000\tcen\tacen\n"
        )
        with tempfile.NamedTemporaryFile("w", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            # Should not sys.exit(1), but warn and keep direction="none"
            cytobands = parse_cytobands(temp_path, self.colors)
            self.assertIn("chr1", cytobands)
            self.assertEqual(cytobands["chr1"][1]["direction"], "none")
        finally:
            os.unlink(temp_path)

    def test_parse_fai(self):
        content = (
            "# header\n"
            "chr1\t248956422\t112\t70\t71\n"
            "chr2\t242193529\t252513197\t70\t71\n"
        )
        with tempfile.NamedTemporaryFile("w", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            fai = list(parse_fai(temp_path))
            self.assertEqual(fai, [("chr1", 248956422), ("chr2", 242193529)])
        finally:
            os.unlink(temp_path)

    def test_annotation_parser(self):
        content = (
            "# header\n"
            "chr1\t100\t200\tGENE1\n"
            "chr1\t100\t200\tGENE1\n"  # Duplicate
            "chr2\t300\t400\tGENE2\textra\n"
        )
        with tempfile.NamedTemporaryFile("w", delete=False) as f:
            f.write(content)
            temp_path = f.name

        try:
            parser = annotation_parser()
            annotations = list(parser(temp_path))
            self.assertEqual(annotations, [("chr1", 100, 200, "GENE1"), ("chr2", 300, 400, "GENE2")])
        finally:
            os.unlink(temp_path)


if __name__ == "__main__":
    unittest.main()
