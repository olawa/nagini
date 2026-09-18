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


class _Namespace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)


class SvSummaryMainTest(unittest.TestCase):
    """main() reads a module-level `snakemake`, so the test supplies one."""

    def _run(self, vcf_text, tmp, caller=None):
        vcf = tmp / "in.vcf"
        vcf.write_text(vcf_text)
        out = {name: tmp / f"out.{name}" for name in ("tsv", "translocations", "json")}
        params = _Namespace(caller=caller) if caller else _Namespace()
        SUMMARY.snakemake = _Namespace(
            input=_Namespace(vcf=str(vcf)),
            output=_Namespace(**{k: str(v) for k, v in out.items()}),
            params=params,
            wildcards=_Namespace(sample="gustaf1", type="T"),
        )
        try:
            SUMMARY.main()
        finally:
            del SUMMARY.snakemake
        rows = out["tsv"].read_text().splitlines()
        return rows[0].split("\t"), [line.split("\t") for line in rows[1:]]

    VCF = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr9\t100\t.\tA\tT\t30\tPASS\tDP=40\n"
        "chr9\t200\t.\tAC\tA\t30\tPASS\tDP=40\n"
        "chr9\t300\t.\tN\t<DEL>\t30\tPASS\tSVTYPE=DEL;END=900;SVLEN=-600\n"
        "chr9\t400\t.\tN\tN[chr21:500[\t30\tPASS\tDP=12\n"
    )

    def test_non_structural_records_are_skipped(self):
        """rindels --sv-lite writes SNVs and indels into the same VCF.

        Without the filter every one of those became an "OTHER" row, which
        turns a 200 MB VCF into a report table nobody can open.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            columns, rows = self._run(self.VCF, Path(tmpdir))

        self.assertEqual([row[columns.index("svtype")] for row in rows], ["DEL", "BND"])

    def test_a_breakend_carries_the_partner_position(self):
        """Naming chr20 without saying where on it is not a breakpoint.

        `end` equals `pos` for a breakend, so it could not stand in for the
        partner coordinate; the value was already being computed for the
        cytoband lookup and the ISCN notation.
        """

        with tempfile.TemporaryDirectory() as tmpdir:
            columns, rows = self._run(self.VCF, Path(tmpdir))

        by_type = {row[columns.index("svtype")]: row for row in rows}
        bnd = by_type["BND"]
        self.assertEqual(bnd[columns.index("chrom")], "chr9")
        self.assertEqual(bnd[columns.index("pos")], "400")
        self.assertEqual(bnd[columns.index("chrom2")], "chr21")
        self.assertEqual(bnd[columns.index("pos2")], "500")

    def test_an_intrachromosomal_sv_has_no_partner_position(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            columns, rows = self._run(self.VCF, Path(tmpdir))
        by_type = {row[columns.index("svtype")]: row for row in rows}
        # A DEL's second coordinate is `end` on the same chromosome.
        self.assertEqual(by_type["DEL"][columns.index("pos2")], "0")
        self.assertEqual(by_type["DEL"][columns.index("end")], "900")

    def test_a_symbolic_translocation_falls_back_to_info_end(self):
        """A TRA written with CHR2/END has no bracket notation to parse."""

        vcf = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr9\t1000\t.\tN\t<TRA>\t60\tPASS\tSVTYPE=TRA;CHR2=chr20;END=2000\n"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            columns, rows = self._run(vcf, Path(tmpdir))
        self.assertEqual(rows[0][columns.index("chrom2")], "chr20")
        self.assertEqual(rows[0][columns.index("pos2")], "2000")

    def test_the_caller_column_tags_each_row(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            columns, rows = self._run(self.VCF, Path(tmpdir), caller="Rindels")
        caller_at = columns.index("caller")
        self.assertEqual({row[caller_at] for row in rows}, {"Rindels"})

    def test_the_caller_defaults_to_sniffles2(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            columns, rows = self._run(self.VCF, Path(tmpdir))
        self.assertEqual({row[columns.index("caller")] for row in rows}, {"Sniffles2"})


if __name__ == "__main__":
    unittest.main()
