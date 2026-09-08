import gzip
import importlib.util
import tempfile
import unittest
from pathlib import Path

SCRIPT = Path(__file__).with_name("cnv_summary.py")
SPEC = importlib.util.spec_from_file_location("cnv_summary", SCRIPT)
CNV_SUMMARY = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(CNV_SUMMARY)


class CnvSummaryTest(unittest.TestCase):
    def test_a_gzipped_cytoband_file_is_read(self):
        """UCSC ships cytoBand gzipped, and plain open() yields binary noise.

        The failure mode is silent: no exception, just no band annotation, so
        this asserts the bands actually come back.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            plain = Path(tmpdir) / "cytoBand.txt"
            rows = "chr9\t0\t2200000\tp24.3\tgneg\nchr9\t2200000\t4600000\tp24.2\tgpos25\n"
            plain.write_text(rows)
            gzipped = Path(tmpdir) / "cytoBand.txt.gz"
            with gzip.open(gzipped, "wt") as handle:
                handle.write(rows)

            for path in (plain, gzipped):
                lookup = CNV_SUMMARY.CytobandLookup(str(path))
                self.assertEqual(lookup.get_band("chr9", 1000), "p24.3", path.name)
                self.assertEqual(lookup.get_band("chr9", 3000000), "p24.2", path.name)

    def test_a_gzipped_target_bed_is_read(self):
        """A design or SV-annotation BED can be bgzipped, as on marvin."""

        with tempfile.TemporaryDirectory() as tmpdir:
            bed = Path(tmpdir) / "targets.bed.gz"
            with gzip.open(bed, "wt") as handle:
                handle.write("chr9\t100\t200\tABL1\n")
            intervals = CNV_SUMMARY.read_target_bed(str(bed))
            self.assertEqual(intervals["chr9"], [(100, 200, "ABL1")])

    def test_classify_event(self):
        self.assertEqual(CNV_SUMMARY.classify_event(1, 2, "Amp", 0.5), "GAIN")
        self.assertEqual(CNV_SUMMARY.classify_event(2, 3, "Amp", 1.2), "AMPLIFICATION")
        self.assertEqual(CNV_SUMMARY.classify_event(0, 1, "Del", -0.8), "LOSS")
        self.assertEqual(CNV_SUMMARY.classify_event(1, 1, "Diploid", 0.0), "DIPLOID")
        self.assertEqual(CNV_SUMMARY.classify_event(0, 2, "LOH", 0.0), "LOH")

    def test_parse_ascat_bed_and_process(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            ascat_path = Path(tmpdir) / "sample.ascat.bed"
            ascat_path.write_text(
                "# ASCAT-style Heuristic Copy Number Segments (exploratory, internal markers)\n"
                "# Heuristic Purity: 0.950\n"
                "# Heuristic Ploidy: 2.100\n"
                "# Columns: chrom, start, end, CN_A, CN_B, BAF, logR, n_sites, state, n_het, n_hom, n_total\n"
                "chr9\t38000000\t45000000\t0\t1\t0.000\t-0.900\t15\tDel\t0\t15\t15\n"
                "chr21\t13000000\t46000000\t2\t3\t0.333\t1.250\t80\tAmp\t50\t30\t80\n"
                "chr1\t10000\t50000000\t1\t1\t0.500\t0.010\t100\tDiploid\t50\t50\t100\n"
            )

            target_bed_path = Path(tmpdir) / "targets.bed"
            target_bed_path.write_text(
                "chr9\t37976085\t44500000\tPAX5\n"
                "chr21\t12991260\t46699124\tRUNX1\n"
            )

            rows, summary = CNV_SUMMARY.process_cnv_segments(
                ascat_bed=str(ascat_path),
                target_bed=str(target_bed_path),
                cytoband_file="",
                sample="TEST_SAMPLE",
                sample_type="T",
            )

            self.assertEqual(summary["purity"], 0.950)
            self.assertEqual(summary["ploidy"], 2.100)
            self.assertEqual(summary["total_segments"], 3)
            self.assertEqual(summary["gains"], 1)
            self.assertEqual(summary["losses"], 1)
            self.assertEqual(summary["target_overlapping_cnvs"], 2)

            # chr9 deletion
            chr9_row = next(r for r in rows if r["chrom"] == "chr9")
            self.assertEqual(chr9_row["event_type"], "LOSS")
            self.assertEqual(chr9_row["copy_number"], 1)
            self.assertTrue(chr9_row["is_target_region"])
            self.assertEqual(chr9_row["target_genes"], "PAX5")

            # chr21 amplification
            chr21_row = next(r for r in rows if r["chrom"] == "chr21")
            self.assertEqual(chr21_row["event_type"], "AMPLIFICATION")
            self.assertEqual(chr21_row["copy_number"], 5)
            self.assertTrue(chr21_row["is_target_region"])
            self.assertEqual(chr21_row["target_genes"], "RUNX1")

            # chr1 diploid
            chr1_row = next(r for r in rows if r["chrom"] == "chr1")
            self.assertEqual(chr1_row["event_type"], "DIPLOID")
            self.assertEqual(chr1_row["copy_number"], 2)
            self.assertFalse(chr1_row["is_target_region"])


if __name__ == "__main__":
    unittest.main()
