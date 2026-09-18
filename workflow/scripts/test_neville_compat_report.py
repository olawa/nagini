import sys
import tempfile
import unittest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from neville_compat_report import (
    SNV_COLUMNS,
    SV_COLUMNS,
    build_neville_rows,
    parse_vcf_header,
    write_tsv,
)

CSQ = (
    "Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|HGVSc|HGVSp|"
    "Existing_variation|CLIN_SIG|MAX_AF|MAX_AF_POPS|VARIANT_CLASS"
)
HEADER = (
    "##fileformat=VCFv4.2\n"
    '##VEP="v111" cache="/data/ref_genomes/VEP"\n'
    f'##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: {CSQ}">\n'
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
)


def csq(allele="T", cons="missense_variant", sym="TP53", hgvsc="ENST0:c.818G>A",
        hgvsp="ENSP0:p.Arg273His", existing="", clin="", maxaf="", pops=""):
    return "|".join([allele, cons, "MODERATE", sym, "ENSG0", "Transcript", "ENST0",
                     "protein_coding", "8/11", hgvsc, hgvsp, existing, clin, maxaf, pops, "SNV"])


def write_vcf(records):
    with tempfile.NamedTemporaryFile("w", suffix=".vcf", delete=False) as tf:
        tf.write(HEADER)
        for r in records:
            tf.write("\t".join(r) + "\n")
        return tf.name


class NevilleCompatTest(unittest.TestCase):
    def test_column_order_matches_neville(self):
        """The whole point is that these diff against neville's sheets."""
        self.assertEqual(SNV_COLUMNS, [
            "Filter", "RunID", "DNAnr", "Gene", "Chr", "Pos", "Ref", "Alt", "AF", "DP",
            "Transcript", "Mutation cds", "ENSP", "Consequence", "COSMIC ids on pos",
            "Clinical Significance", "dbSNP", "Max Pop AF", "Max Pop", "Callers",
        ])
        # SVs is the same list with "SV length" inserted after Alt.
        self.assertEqual(SV_COLUMNS[8], "SV length")
        self.assertEqual([c for c in SV_COLUMNS if c != "SV length"], SNV_COLUMNS)

    def test_parses_csq_layout_from_header(self):
        path = write_vcf([])
        fields, vep_line = parse_vcf_header(path)
        self.assertIn("SYMBOL", fields)
        self.assertIn("MAX_AF_POPS", fields)
        self.assertIn("v111", vep_line)
        Path(path).unlink()

    def test_af_is_a_fraction_not_a_percentage(self):
        """neville writes AF as 0-1; a percentage here would be 100x off."""
        path = write_vcf([(
            "chr17", "7673802", ".", "C", "T", "60", "PASS",
            "CSQ=" + csq(), "GT:DP:AD:AF", "0/1:800:700,100",
        )])
        snv, _ = build_neville_rows(path, run_id="R1", sample="S1")
        self.assertEqual(snv[0][SNV_COLUMNS.index("AF")], 0.125)
        Path(path).unlink()

    def test_af_falls_back_to_ad_over_dp(self):
        path = write_vcf([(
            "chr17", "7673802", ".", "C", "T", "60", "PASS",
            "CSQ=" + csq(), "GT:DP:AD", "0/1:800:700,100",
        )])
        snv, _ = build_neville_rows(path, run_id="R1", sample="S1")
        self.assertEqual(snv[0][SNV_COLUMNS.index("AF")], 0.125)
        Path(path).unlink()

    def test_split_at_fifty_bp(self):
        """An ITD over 50 bp belongs in neville's SVs sheet, not SNV_Indels."""
        long_alt = "A" + "G" * 78
        path = write_vcf([
            ("chr5", "171410539", ".", "C", "CTCTG", "60", "PASS",
             "CSQ=" + csq(allele="TCTG", sym="NPM1"), "GT:DP:AD:AF", "0/1:1200:744,456"),
            ("chr13", "28034100", ".", "A", long_alt, "60", "PASS",
             "CSQ=" + csq(allele="G" * 78, sym="FLT3"), "GT:DP:AD:AF", "0/1:1000:700,300"),
        ])
        snv, sv = build_neville_rows(path, run_id="R1", sample="S1")
        self.assertEqual([r[SNV_COLUMNS.index("Gene")] for r in snv], ["NPM1"])
        self.assertEqual([r[SV_COLUMNS.index("Gene")] for r in sv], ["FLT3"])
        self.assertEqual(sv[0][SV_COLUMNS.index("SV length")], 79)
        Path(path).unlink()

    def test_soft_filtered_variants_are_kept_with_their_flag(self):
        """neville soft-filters; the flag goes in the Filter column, nothing is dropped."""
        path = write_vcf([(
            "chr17", "7674850", ".", "A", "G", "30", "VAF;Biotype",
            "CSQ=" + csq(), "GT:DP:AD:AF", "0/1:900:882,18",
        )])
        snv, _ = build_neville_rows(path, run_id="R1", sample="S1")
        self.assertEqual(len(snv), 1)
        self.assertEqual(snv[0][SNV_COLUMNS.index("Filter")], "VAF;Biotype")
        Path(path).unlink()

    def test_below_min_af_is_dropped(self):
        path = write_vcf([(
            "chr17", "7674860", ".", "A", "G", "30", "PASS",
            "CSQ=" + csq(), "GT:DP:AD:AF", "0/1:1000:998,2",
        )])
        snv, sv = build_neville_rows(path, run_id="R1", sample="S1")
        self.assertEqual((snv, sv), ([], []))
        Path(path).unlink()

    def test_existing_variation_splits_into_cosmic_and_dbsnp(self):
        path = write_vcf([(
            "chr17", "7673802", ".", "C", "T", "60", "PASS",
            "CSQ=" + csq(existing="COSV53360203&rs28934576", clin="pathogenic",
                         maxaf="0.0001", pops="gnomADe"),
            "GT:DP:AD:AF", "0/1:800:700,100",
        )])
        snv, _ = build_neville_rows(path, run_id="R1", sample="S1")
        row = snv[0]
        self.assertEqual(row[SNV_COLUMNS.index("COSMIC ids on pos")], "COSV53360203")
        self.assertEqual(row[SNV_COLUMNS.index("dbSNP")], "rs28934576")
        self.assertEqual(row[SNV_COLUMNS.index("Clinical Significance")], "pathogenic")
        self.assertEqual(row[SNV_COLUMNS.index("Max Pop AF")], "0.0001")
        Path(path).unlink()

    def test_hgvsc_splits_into_transcript_and_cds(self):
        path = write_vcf([(
            "chr17", "7673802", ".", "C", "T", "60", "PASS",
            "CSQ=" + csq(), "GT:DP:AD:AF", "0/1:800:700,100",
        )])
        snv, _ = build_neville_rows(path, run_id="R1", sample="S1")
        self.assertEqual(snv[0][SNV_COLUMNS.index("Transcript")], "ENST0")
        self.assertEqual(snv[0][SNV_COLUMNS.index("Mutation cds")], "c.818G>A")
        Path(path).unlink()

    def test_run_id_and_sample_populate_neville_columns(self):
        path = write_vcf([(
            "chr17", "7673802", ".", "C", "T", "60", "PASS",
            "CSQ=" + csq(), "GT:DP:AD:AF", "0/1:800:700,100",
        )])
        snv, _ = build_neville_rows(path, run_id="RUN123", sample="D25-test007")
        self.assertEqual(snv[0][SNV_COLUMNS.index("RunID")], "RUN123")
        self.assertEqual(snv[0][SNV_COLUMNS.index("DNAnr")], "D25-test007")
        self.assertEqual(snv[0][SNV_COLUMNS.index("Callers")], "rindels")
        Path(path).unlink()

    def test_tsv_round_trip_has_the_header_row(self):
        path = write_vcf([(
            "chr17", "7673802", ".", "C", "T", "60", "PASS",
            "CSQ=" + csq(), "GT:DP:AD:AF", "0/1:800:700,100",
        )])
        snv, _ = build_neville_rows(path, run_id="R1", sample="S1")
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as tf:
            out = tf.name
        write_tsv(snv, SNV_COLUMNS, out)
        lines = Path(out).read_text().splitlines()
        self.assertEqual(lines[0].split("\t"), SNV_COLUMNS)
        self.assertEqual(len(lines), 2)
        Path(out).unlink()
        Path(path).unlink()


if __name__ == "__main__":
    unittest.main()
