import gzip
import importlib.util
import json
import tempfile
import unittest
from pathlib import Path

import pysam

HERE = Path(__file__).parent


def load_module(name: str):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FUSION_COUNTER = load_module("bcr_abl_fusion_counter")
BRANFORD_CALLER = load_module("bcr_abl_branford_caller")
PHASER = load_module("bcr_abl_phaser")
REPORT = load_module("nagini_report")


def create_synthetic_bam(bam_path: Path, reads: list[dict], references: list[tuple[str, int]]):
    """Create a synthetic BAM file for testing."""
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": name, "LN": length} for name, length in references]}
    with pysam.AlignmentFile(str(bam_path), "wb", header=header) as out_bam:
        for r_idx, r in enumerate(reads):
            aln = pysam.AlignedSegment()
            aln.query_name = r.get("name", f"read_{r_idx}")
            aln.query_sequence = r["seq"]
            aln.query_qualities = [30] * len(r["seq"])
            aln.reference_id = [ref[0] for ref in references].index(r["chrom"])
            aln.reference_start = r["start"]  # 0-based
            aln.cigar = r.get("cigar", [(0, len(r["seq"]))])  # default match
            aln.mapping_quality = r.get("mapq", 60)
            aln.is_supplementary = r.get("supplementary", False)
            aln.is_reverse = r.get("is_reverse", False)
            out_bam.write(aln)
    pysam.index(str(bam_path))


class TestBcrAblBranfordCaller(unittest.TestCase):
    def test_short_aa_notation(self):
        self.assertEqual(BRANFORD_CALLER.short_aa_notation("p.Thr315Ile"), "T315I")
        self.assertEqual(BRANFORD_CALLER.short_aa_notation("p.Glu255Lys"), "E255K")
        self.assertEqual(BRANFORD_CALLER.short_aa_notation("p.Tyr253His"), "Y253H")
        self.assertEqual(BRANFORD_CALLER.short_aa_notation("p.Val299Leu"), "V299L")

    def test_branford_calling_on_synthetic_bam(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "test.bam"
            bed_path = tmp / "branford.bed"
            all_sites_tsv = tmp / "all_sites.tsv"
            muts_tsv = tmp / "mutations.tsv"
            vcf_path = tmp / "mutations.vcf"

            # Create test BED: T315 (pos 100) and E255 (pos 50) on test chr9
            bed_path.write_text(
                "chr9\t49\t50\tNM_1:c.(763G>A)\tG\tA\tp.Glu255Lys\n"
                "chr9\t99\t100\tNM_1:c.(944C>T)\tC\tT\tp.Thr315Ile\n"
            )

            # Reads spanning pos 50 (index 49) and pos 100 (index 99)
            # Reference: 150 bp chr9
            # 80 reads with WT (pos 50 = G, pos 100 = C)
            # 20 reads with T315I mutation (pos 100 = T)
            reads = []
            for i in range(80):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "C"
                reads.append({"name": f"wt_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})
            for i in range(20):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "T"
                reads.append({"name": f"mut_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})

            create_synthetic_bam(bam_path, reads, [("chr9", 200)])

            # Run caller main
            BRANFORD_CALLER.sys.argv = [
                "bcr_abl_branford_caller.py",
                "--bam", str(bam_path),
                "--bed", str(bed_path),
                "--sample", "TEST_SAMPLE",
                "--output-all", str(all_sites_tsv),
                "--output-mutations", str(muts_tsv),
                "--output-vcf", str(vcf_path),
                "--min-af", "0.05",
                "--min-depth", "50",
            ]
            BRANFORD_CALLER.main()

            self.assertTrue(all_sites_tsv.is_file())
            self.assertTrue(muts_tsv.is_file())
            self.assertTrue(vcf_path.is_file())

            mut_content = muts_tsv.read_text().splitlines()
            records = [line for line in mut_content if not line.startswith("#") and line.strip()]
            self.assertEqual(len(records), 2)  # Header + 1 detected mutation
            self.assertIn("T315I", records[1])
            self.assertIn("DETECTED", records[1])
            self.assertIn("20.0", records[1])  # 20% VAF


class TestBcrAblPhaser(unittest.TestCase):
    def test_single_molecule_cis_phasing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "test_phasing.bam"
            muts_tsv = tmp / "mutations.tsv"
            phasing_tsv = tmp / "phasing.tsv"
            phasing_json = tmp / "phasing.json"

            # 2 mutations: E255K at pos 50 (G>A), T315I at pos 100 (C>T)
            muts_tsv.write_text(
                "sample\tchrom\tpos\tref\talt\tmutation\taa_change\tname\tdepth\tad_ref\tad_alt\tvaf_percent\tstatus\n"
                "SAMPLE_1\tchr9\t50\tG\tA\tE255K\tp.Glu255Lys\tNM_1\t100\t80\t20\t20.0\tDETECTED\n"
                "SAMPLE_1\tchr9\t100\tC\tT\tT315I\tp.Thr315Ile\tNM_1\t100\t80\t20\t20.0\tDETECTED\n"
            )

            # Create 20 reads with BOTH E255K (A at 49) and T315I (T at 99) -> CIS
            # and 80 WT reads (G at 49, C at 99)
            reads = []
            for i in range(80):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "C"
                reads.append({"name": f"wt_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})
            for i in range(20):
                seq = list("N" * 150)
                seq[49] = "A"
                seq[99] = "T"
                reads.append({"name": f"cis_mut_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})

            create_synthetic_bam(bam_path, reads, [("chr9", 200)])

            PHASER.sys.argv = [
                "bcr_abl_phaser.py",
                "--bam", str(bam_path),
                "--mutations-tsv", str(muts_tsv),
                "--sample", "SAMPLE_1",
                "--output-tsv", str(phasing_tsv),
                "--output-json", str(phasing_json),
                "--min-spanning-reads", "10",
            ]
            PHASER.main()

            self.assertTrue(phasing_tsv.is_file())
            content = phasing_tsv.read_text().splitlines()
            records = [line for line in content if not line.startswith("#") and line.strip()]
            self.assertEqual(len(records), 2)  # Header + 1 pair
            self.assertIn("CIS (Compound)", records[1])

            with phasing_json.open() as h:
                data = json.load(h)
            self.assertEqual(data["phased_pairs_count"], 1)
            pair = data["phased_pairs"][0]
            self.assertEqual(pair["cis_reads_11"], 20)
            self.assertEqual(pair["trans_reads_10"], 0)
            self.assertEqual(pair["trans_reads_01"], 0)
            self.assertEqual(pair["phase_state"], "CIS (Compound)")

    def test_single_molecule_trans_phasing(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            bam_path = tmp / "test_trans.bam"
            muts_tsv = tmp / "mutations.tsv"
            phasing_tsv = tmp / "phasing.tsv"
            phasing_json = tmp / "phasing.json"

            muts_tsv.write_text(
                "sample\tchrom\tpos\tref\talt\tmutation\taa_change\tname\tdepth\tad_ref\tad_alt\tvaf_percent\tstatus\n"
                "SAMPLE_2\tchr9\t50\tG\tA\tE255K\tp.Glu255Lys\tNM_1\t100\t90\t10\t10.0\tDETECTED\n"
                "SAMPLE_2\tchr9\t100\tC\tT\tT315I\tp.Thr315Ile\tNM_1\t100\t90\t10\t10.0\tDETECTED\n"
            )

            # 10 reads with E255K alone (A at 49, C at 99)
            # 10 reads with T315I alone (G at 49, T at 99)
            # 80 reads with WT (G at 49, C at 99) -> TRANS
            reads = []
            for i in range(80):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "C"
                reads.append({"name": f"wt_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})
            for i in range(10):
                seq = list("N" * 150)
                seq[49] = "A"
                seq[99] = "C"
                reads.append({"name": f"mut1_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})
            for i in range(10):
                seq = list("N" * 150)
                seq[49] = "G"
                seq[99] = "T"
                reads.append({"name": f"mut2_read_{i}", "chrom": "chr9", "start": 0, "seq": "".join(seq)})

            create_synthetic_bam(bam_path, reads, [("chr9", 200)])

            PHASER.sys.argv = [
                "bcr_abl_phaser.py",
                "--bam", str(bam_path),
                "--mutations-tsv", str(muts_tsv),
                "--sample", "SAMPLE_2",
                "--output-tsv", str(phasing_tsv),
                "--output-json", str(phasing_json),
            ]
            PHASER.main()

            with phasing_json.open() as h:
                data = json.load(h)
            pair = data["phased_pairs"][0]
            self.assertEqual(pair["cis_reads_11"], 0)
            self.assertEqual(pair["trans_reads_10"], 10)
            self.assertEqual(pair["trans_reads_01"], 10)
            self.assertEqual(pair["phase_state"], "TRANS (Polyclonal)")


class TestBcrAblFusionCounter(unittest.TestCase):
    def test_assign_exon(self):
        exon_map = {"e1": 100, "e13": 500, "e14": 1000}
        self.assertEqual(FUSION_COUNTER.assign_exon(105, exon_map, tolerance=50), "e1")
        self.assertEqual(FUSION_COUNTER.assign_exon(510, exon_map, tolerance=50), "e13")
        self.assertIsNone(FUSION_COUNTER.assign_exon(750, exon_map, tolerance=50))


class TestBcrAblReportIntegration(unittest.TestCase):
    def test_generates_workbook_with_bcr_abl_tab(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            out_xlsx = tmp / "report.xlsx"
            fusion_tsv = tmp / "fusion.tsv"
            muts_tsv = tmp / "muts.tsv"
            phasing_tsv = tmp / "phasing.tsv"

            fusion_tsv.write_text(
                "sample\tfusion_isoform\talt_name\tprotein_form\tcluster\tbcr_exon\tabl_exon\tread_count\tfraction_percent\tstatus\tdescription\n"
                "SAMPLE_1\te14a2\tb3a2\tp210\tM-BCR\te14\ta2\t1500\t95.5\tPRIMARY\tMajor CML fusion\n"
            )
            muts_tsv.write_text(
                "sample\tchrom\tpos\tref\talt\tmutation\taa_change\tname\tdepth\tad_ref\tad_alt\tvaf_percent\tstatus\timatinib\tdasatinib\tnilotinib\tbosutinib\tponatinib\tasciminib\n"
                "SAMPLE_1\tchr9\t130872896\tC\tT\tT315I\tp.Thr315Ile\tNM_1\t5000\t4000\t1000\t20.0\tDETECTED\tResistant\tResistant\tResistant\tResistant\tSensitive\tSensitive\n"
            )
            phasing_tsv.write_text(
                "sample\tmutation_1\tmutation_2\tphase_state\tcis_fraction_percent\tcis_reads_11\ttrans_reads_10\ttrans_reads_01\twt_reads_00\tother_reads\ttotal_spanning_reads\tpos_1\tpos_2\tref_1\talt_1\tref_2\talt_2\tclinical_interpretation\n"
            )

            REPORT.generate_report(
                output=str(out_xlsx),
                sample="SAMPLE_1",
                sample_type="T",
                reference_build="hg38",
                annotation="none",
                snv_vcf="",
                sv_tsv="",
                translocations_tsv="",
                wgs_coverage="",
                target_coverage="",
                cramino="",
                bcr_abl_fusion=str(fusion_tsv),
                bcr_abl_mutations=str(muts_tsv),
                bcr_abl_phasing=str(phasing_tsv),
            )

            self.assertTrue(out_xlsx.is_file())
            self.assertGreater(out_xlsx.stat().st_size, 1000)


if __name__ == "__main__":
    unittest.main()
