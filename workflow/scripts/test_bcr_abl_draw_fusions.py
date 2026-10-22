import argparse
import os
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace

import bcr_abl_draw_fusions as draw


def completed(returncode=0):
    return SimpleNamespace(returncode=returncode)


class Fake:
    """A stand-in for subprocess.run that creates the files each tool would."""

    def __init__(self, r_ok=True, r_writes=True, fallback_ok=True, fallback_writes=True, pdftoppm_png=b"page"):
        self.calls = []
        self.r_ok, self.r_writes = r_ok, r_writes
        self.fallback_ok, self.fallback_writes = fallback_ok, fallback_writes
        self.pdftoppm_png = pdftoppm_png

    def __call__(self, command):
        self.calls.append(list(command))
        head = os.path.basename(command[0])
        if head == "pdftoppm":
            prefix = command[-1]
            Path(f"{prefix}-1.png").write_bytes(self.pdftoppm_png)
            return completed(0)
        if "draw_bcr_abl_fusion.py" in " ".join(command):
            out = command[command.index("--output") + 1]
            if self.fallback_writes:
                Path(out).write_bytes(b"python-diagram")
            return completed(0 if self.fallback_ok else 1)
        # draw_fusions.R
        out = next(a.split("=", 1)[1] for a in command if a.startswith("--output="))
        if self.r_writes:
            Path(out).write_bytes(b"r-plot")
        return completed(0 if self.r_ok else 1)

    def of(self, needle):
        return [c for c in self.calls if needle in " ".join(c)]


class DiscoveryTest(unittest.TestCase):
    def test_configured_file_wins_when_it_exists(self):
        self.assertEqual(draw.first_existing("a", ["b"], isfile=lambda p: p in {"a", "b"}), "a")

    def test_falls_back_to_the_first_existing_candidate(self):
        self.assertEqual(draw.first_existing("a", ["x", "y", "z"], isfile=lambda p: p in {"y", "z"}), "y")

    def test_candidates_follow_the_genome_build(self):
        """The Arriba image holds both builds; hg19 used to win for hg38 runs
        only because it was listed first."""

        everything = lambda p: True  # noqa: E731 -- both builds present, as in the image
        for build, want, never in (("hg38", "hg38", "hg19"), ("hg19", "hg19", "hg38")):
            cyto = draw.first_existing("", draw.candidates_for_build(draw.CYTOBAND_CANDIDATES, build), everything)
            dom = draw.first_existing("", draw.candidates_for_build(draw.PROTEIN_DOMAIN_CANDIDATES, build), everything)
            for chosen in (cyto, dom):
                self.assertIn(want, chosen.lower(), f"{build}: {chosen}")
                self.assertNotIn(never, chosen.lower(), f"{build}: {chosen}")

    def test_unknown_build_gets_no_reference_file(self):
        self.assertEqual(draw.candidates_for_build(draw.CYTOBAND_CANDIDATES, "t2t"), [])

    def test_configured_file_still_wins_regardless_of_build(self):
        chosen = draw.first_existing(
            "/my/cytobands.tsv", draw.candidates_for_build(draw.CYTOBAND_CANDIDATES, "hg38"), lambda p: True
        )
        self.assertEqual(chosen, "/my/cytobands.tsv")

    def test_nothing_found_is_empty(self):
        self.assertEqual(draw.first_existing("", ["x"], isfile=lambda p: False), "")

    def test_bundled_script_is_preferred_and_run_with_rscript(self):
        cmd = draw.find_draw_command(Path("/b/draw_fusions.r"), isfile=lambda p: True, which=lambda n: "/usr/bin/x")
        self.assertEqual(cmd, ["Rscript", "/b/draw_fusions.r"])

    def test_path_then_known_install_locations(self):
        self.assertEqual(
            draw.find_draw_command(Path("/none"), isfile=lambda p: False, which=lambda n: "/usr/bin/draw_fusions.R"),
            ["draw_fusions.R"],
        )
        self.assertEqual(
            draw.find_draw_command(Path("/none"), isfile=lambda p: p == "/usr/local/bin/draw_fusions.R", which=lambda n: None),
            ["/usr/local/bin/draw_fusions.R"],
        )
        self.assertEqual(draw.find_draw_command(Path("/none"), isfile=lambda p: False, which=lambda n: None), [])


class ArgumentTest(unittest.TestCase):
    def test_panels_follow_the_available_reference_files(self):
        self.assertEqual(draw.choose_panels("c", "p"), "fusion,circos,domains,readcounts")
        self.assertEqual(draw.choose_panels("c", ""), "fusion,circos,readcounts")
        self.assertEqual(draw.choose_panels("", "p"), "fusion,domains,readcounts")
        self.assertEqual(draw.choose_panels("", ""), "fusion,readcounts")

    def test_draw_arguments(self):
        args = draw.build_draw_args(
            ["Rscript", "d.r"], fusions="f.tsv", bam="a.bam", output="o.pdf", gtf="g.gtf",
            sample="S1", cytobands="cy.tsv", protein_domains="", extra="--minConfidenceForCircosPlot=low --foo 1",
        )
        self.assertEqual(args[:2], ["Rscript", "d.r"])
        self.assertIn("--sampleName=S1", args)
        self.assertIn("--transcriptSelection=canonical", args)
        self.assertIn("--plotPanels=fusion,circos,readcounts", args)
        self.assertIn("--cytobands=cy.tsv", args)
        self.assertFalse(any(a.startswith("--proteinDomains") for a in args))
        self.assertEqual(args[-3:], ["--minConfidenceForCircosPlot=low", "--foo", "1"])

    def test_primary_fusion_table_keeps_header_and_first_row(self):
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = Path(tmp) / "f.tsv", Path(tmp) / "p.tsv"
            src.write_text("h\nrow1\nrow2\n")
            draw.primary_fusion_table(str(src), str(dst))
            self.assertEqual(dst.read_text(), "h\nrow1\n")


class DrawTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.root = Path(self.tmp.name)
        (self.root / "fusions.tsv").write_text("#gene1\tgene2\nBCR\tABL1\nBCR\tABL1\n")
        (self.root / "ref.gtf").write_text("x")
        self.args = argparse.Namespace(
            sample="S1", fusions=str(self.root / "fusions.tsv"),
            fusion_summary="summary.tsv", mutations="m.tsv", phasing="p.tsv", exon_coverage="e.tsv",
            bam="a.bam", gtf=str(self.root / "ref.gtf"), cytobands="", protein_domains="",
            build="hg19", extra="",
            output_pdf=str(self.root / "plot.pdf"), output_png=str(self.root / "plot.png"),
        )
        self.messages = []

    def tearDown(self):
        self.tmp.cleanup()

    def go(self, fake, draw_available=True, pdftoppm=False):
        bundled = self.root / "draw_fusions.r"
        if draw_available:
            bundled.write_text("# r")
        return draw.draw_fusion_plots(
            self.args, run=fake,
            which=lambda n: "/usr/bin/pdftoppm" if (n == "pdftoppm" and pdftoppm) else None,
            bundled=bundled, fallback_script=self.root / "draw_bcr_abl_fusion.py",
            python="python3", log=self.messages.append,
        )

    def test_arriba_plots_are_used_when_they_work(self):
        fake = Fake()
        self.assertEqual(self.go(fake), 0)
        self.assertEqual(len(fake.of("Rscript")), 2)
        self.assertEqual(fake.of("draw_bcr_abl_fusion.py"), [])
        self.assertEqual(Path(self.args.output_pdf).read_bytes(), b"r-plot")

    def test_png_is_drawn_from_the_primary_fusion_only(self):
        fake = Fake()
        self.go(fake)
        pdf_call, png_call = fake.of("Rscript")
        self.assertIn(f"--fusions={self.args.fusions}", pdf_call)
        primary = next(a for a in png_call if a.startswith("--fusions=")).split("=", 1)[1]
        self.assertNotEqual(primary, self.args.fusions)

    def test_first_pdf_page_replaces_the_png_when_pdftoppm_exists(self):
        fake = Fake()
        self.go(fake, pdftoppm=True)
        self.assertEqual(Path(self.args.output_png).read_bytes(), b"page")

    def test_no_draw_script_goes_straight_to_the_python_diagram(self):
        fake = Fake()
        self.assertEqual(self.go(fake, draw_available=False), 0)
        self.assertEqual(fake.of("Rscript"), [])
        self.assertEqual(len(fake.of("draw_bcr_abl_fusion.py")), 2)
        self.assertEqual(Path(self.args.output_pdf).read_bytes(), b"python-diagram")

    def test_missing_gtf_skips_arriba(self):
        self.args.gtf = str(self.root / "absent.gtf")
        fake = Fake()
        self.assertEqual(self.go(fake), 0)
        self.assertEqual(fake.of("Rscript"), [])

    def test_failed_r_run_with_empty_pdf_falls_back(self):
        fake = Fake(r_ok=False, r_writes=False)
        self.assertEqual(self.go(fake), 0)
        self.assertEqual(len(fake.of("Rscript")), 2)
        self.assertEqual(len(fake.of("draw_bcr_abl_fusion.py")), 2)

    def test_only_the_empty_deliverable_is_redrawn(self):
        Path(self.args.output_pdf).write_bytes(b"already there")
        fake = Fake(r_writes=False)
        self.go(fake)
        fallback = fake.of("draw_bcr_abl_fusion.py")
        self.assertEqual(len(fallback), 1)
        self.assertIn(self.args.output_png, fallback[0])

    def test_failing_fallback_for_the_pdf_is_an_error(self):
        fake = Fake(r_writes=False, fallback_ok=False, fallback_writes=False)
        with self.assertRaises(RuntimeError):
            self.go(fake, draw_available=False)

    def test_missing_png_at_the_end_fails_the_step(self):
        class PdfOnly(Fake):
            def __call__(self, command):
                result = super().__call__(command)
                if "draw_bcr_abl_fusion.py" in " ".join(command) and command[command.index("--output") + 1].endswith(".png"):
                    Path(command[command.index("--output") + 1]).unlink(missing_ok=True)
                    return completed(1)
                return result

        self.assertEqual(self.go(PdfOnly(), draw_available=False), 1)
        self.assertTrue(any("No plot was produced" in m for m in self.messages))


if __name__ == "__main__":
    unittest.main()
