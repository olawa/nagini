"""Tests for the barcode yield table."""

import importlib.util
import tempfile
import unittest
from pathlib import Path


HERE = Path(__file__).parent
_spec = importlib.util.spec_from_file_location(
    "dorado_demux_summary", HERE / "dorado_demux_summary.py"
)
summary = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(summary)


def write_summary(path, rows, header="read_id\tbarcode"):
    path.write_text(header + "\n" + "\n".join(rows) + "\n")


def write_units(path, pairs):
    lines = ["sample\tbarcode_name"]
    lines += [f"{sample}\t{barcode}" for sample, barcode in pairs]
    path.write_text("\n".join(lines) + "\n")


def read_table(path):
    rows = []
    header = None
    for line in path.read_text().splitlines():
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if header is None:
            header = fields
            continue
        rows.append(dict(zip(header, fields)))
    return rows


class BarcodeNormalisationTest(unittest.TestCase):
    def test_barcode_spellings_collapse(self):
        for value in ("barcode05", "BC5", "bc05", "5", "05"):
            self.assertEqual(summary.normalize_barcode(value), "barcode05")

    def test_unclassified_is_recognised(self):
        self.assertEqual(summary.normalize_barcode("unclassified"), "unclassified")
        self.assertEqual(summary.normalize_barcode("unknown"), "unclassified")

    def test_empty_values_are_dropped(self):
        for value in ("", "  ", "NA", "none", "-"):
            self.assertEqual(summary.normalize_barcode(value), "")


class CountingTest(unittest.TestCase):
    def test_reads_are_counted_per_barcode(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "s.tsv"
            write_summary(path, [f"r{i}\tbarcode01" for i in range(3)]
                                + ["r9\tbarcode02", "r10\tunclassified"])
            counts = summary.count_barcodes(path)
        self.assertEqual(counts["barcode01"], 3)
        self.assertEqual(counts["barcode02"], 1)
        self.assertEqual(counts["unclassified"], 1)

    def test_alternative_column_names_are_accepted(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "s.tsv"
            write_summary(path, ["r1\tx\tbarcode03"],
                          header="read_id\tfilename\tbarcode_arrangement")
            counts = summary.count_barcodes(path)
        self.assertEqual(counts["barcode03"], 1)

    def test_a_summary_without_a_barcode_column_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "s.tsv"
            write_summary(path, ["r1\t100"], header="read_id\tlength")
            with self.assertRaises(SystemExit):
                summary.count_barcodes(path)


class TableTest(unittest.TestCase):
    def _build(self, tmp, summary_rows, unit_pairs):
        root = Path(tmp)
        summary_path = root / "s.tsv"
        write_summary(summary_path, summary_rows)
        units_path = root / "units.tsv"
        write_units(units_path, unit_pairs)
        output = root / "qc" / "dorado_demux_mqc.tsv"
        counts = summary.count_barcodes(summary_path)
        rows = summary.build_rows(counts, summary.read_barcode_samples(units_path))
        summary.write_mqc(output, rows, "SQK-NBD114-24")
        return output

    def test_expected_barcode_with_no_reads_is_still_a_row(self):
        """A missing row and a zero row mean different things."""

        with tempfile.TemporaryDirectory() as tmp:
            output = self._build(
                tmp,
                ["r1\tbarcode01", "r2\tbarcode01"],
                [("PAT01", "barcode01"), ("PAT02", "barcode02")],
            )
            rows = {row["Barcode"]: row for row in read_table(output)}
        self.assertEqual(rows["barcode01"]["Reads"], "2")
        self.assertEqual(rows["barcode01"]["Status"], "OK")
        self.assertEqual(rows["barcode02"]["Reads"], "0")
        self.assertEqual(rows["barcode02"]["Status"], "NO READS")
        self.assertEqual(rows["barcode02"]["Sample"], "PAT02")

    def test_percentages_separate_classified_from_the_whole_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = self._build(
                tmp,
                ["r1\tbarcode01", "r2\tbarcode02", "r3\tunclassified", "r4\tunclassified"],
                [("PAT01", "barcode01"), ("PAT02", "barcode02")],
            )
            rows = {row["Barcode"]: row for row in read_table(output)}
        # Two classified reads, one each; four reads in the run.
        self.assertEqual(rows["barcode01"]["Percent of classified"], "50.0")
        self.assertEqual(rows["barcode01"]["Percent of run"], "25.0")
        self.assertEqual(rows["unclassified"]["Percent of run"], "50.0")

    def test_output_is_multiqc_custom_content(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = self._build(tmp, ["r1\tbarcode01"], [("PAT01", "barcode01")])
            text = output.read_text()
        self.assertTrue(output.name.endswith("_mqc.tsv"))
        self.assertIn("# id: 'dorado_demux'", text)
        self.assertIn("# plot_type: 'table'", text)
        self.assertIn("SQK-NBD114-24", text)

    def test_barcodes_without_a_sample_still_appear(self):
        with tempfile.TemporaryDirectory() as tmp:
            output = self._build(tmp, ["r1\tbarcode07"], [("PAT01", "barcode01")])
            rows = {row["Barcode"]: row for row in read_table(output)}
        self.assertIn("barcode07", rows)
        self.assertEqual(rows["barcode07"]["Sample"], "")


if __name__ == "__main__":
    unittest.main()
