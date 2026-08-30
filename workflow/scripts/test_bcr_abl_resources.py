"""Regression tests for the BCR-ABL1 reference resources.

These lock the Branford BED files and the fusion-junction coordinate table
against pickett_bcr_abl_pipeline's own annotation files, so the two cannot
silently drift apart again the way they did the first time: 71 rows instead
of 73, a fabricated variant, 9 mislabelled amino-acid changes, and a
BCR-exon-13/14 swap in the fusion counter that made every e14a2 read report
as e13a2 and vice versa.

Tests that need a file outside this repository (pickett's checkout, a
reference FASTA) skip cleanly when it is not present, rather than failing --
this suite is expected to run in CI where neither exists.
"""

import importlib.util
import re
import subprocess
import unittest
from pathlib import Path


HERE = Path(__file__).parent
RESOURCES = HERE.parent / "resources" / "bcr_abl"
PICKETT_BED = (
    Path.home()
    / "dev/cgu_pipelines/wp2/pickett_bcr_abl_pipeline/config/annotations_hg19"
    / "branford_SNVs_230118.bed"
)

_spec = importlib.util.spec_from_file_location(
    "bcr_abl_fusion_counter", HERE / "bcr_abl_fusion_counter.py"
)
fusion_counter = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(fusion_counter)


def read_bed(path):
    rows = []
    for line in path.read_text().splitlines():
        if line.strip():
            rows.append(line.split("\t"))
    return rows


class BranfordBedTest(unittest.TestCase):
    def test_hg19_and_hg38_have_73_rows(self):
        for name in ("branford_positions_hg19.bed", "branford_positions_hg38.bed"):
            rows = read_bed(RESOURCES / name)
            self.assertEqual(len(rows), 73, f"{name} has {len(rows)} rows, expected 73")

    def test_hg38_is_hg19_shifted_by_the_fixed_offset(self):
        hg19 = read_bed(RESOURCES / "branford_positions_hg19.bed")
        hg38 = read_bed(RESOURCES / "branford_positions_hg38.bed")
        self.assertEqual(len(hg19), len(hg38))
        for row19, row38 in zip(hg19, hg38):
            self.assertEqual(int(row19[1]) - int(row38[1]), 2875387)
            self.assertEqual(int(row19[2]) - int(row38[2]), 2875387)
            # Everything except the coordinates must be identical between builds.
            self.assertEqual(row19[3:], row38[3:])

    def test_t315i_is_at_the_published_hg38_position(self):
        """chr9:130872896 C>T, quoted identically by two independent sources
        (rindels' hg38 resistance-mutation list and the nanoMonitor clinical
        mutation table) before this pipeline's own list existed."""

        hg38 = read_bed(RESOURCES / "branford_positions_hg38.bed")
        t315i = [row for row in hg38 if row[6] == "p.Thr315Ile"]
        self.assertEqual(len(t315i), 1)
        self.assertEqual(t315i[0][0], "chr9")
        self.assertEqual(int(t315i[0][2]), 130872896)
        self.assertEqual((t315i[0][4], t315i[0][5]), ("C", "T"))

    # Inherited as-is from pickett_bcr_abl_pipeline's own production Branford
    # list (chr9:133748407-133748408, NM_005157.6:c.(1069A>G)): the label
    # K375E implies codon 375, but cDNA position 1069 is codon 357. The
    # genomic position, REF and ALT are independently confirmed correct
    # (validated against a real hg19 FASTA and hg38 elsewhere in this file);
    # only the human-readable codon number in pickett's own label disagrees
    # with its own cDNA numbering. Not something to silently correct without
    # the source Branford publication -- flagged here rather than fixed, and
    # worth a question to whoever maintains pickett's list.
    _KNOWN_UPSTREAM_CODON_LABEL_MISMATCH = {"NM_005157.6:c.(1069A>G)"}

    def test_amino_acid_column_is_internally_consistent(self):
        """Catches the AA1<->AA3 conversion errors the hand-typed list had:
        wrong residue, wrong position number, or ref/alt letters swapped."""

        pattern = re.compile(r"p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})")
        hgvs_pattern = re.compile(r"c\.\((\d+)([ACGT])>([ACGT])\)")
        for name in ("branford_positions_hg19.bed", "branford_positions_hg38.bed"):
            for row in read_bed(RESOURCES / name):
                hgvs, aa = row[3], row[6]
                if hgvs in self._KNOWN_UPSTREAM_CODON_LABEL_MISMATCH:
                    continue
                aa_match = pattern.fullmatch(aa)
                self.assertIsNotNone(aa_match, f"{name}: unparseable AA change {aa!r}")
                cdna_match = hgvs_pattern.search(hgvs)
                self.assertIsNotNone(cdna_match, f"{name}: unparseable HGVS {hgvs!r}")
                # The codon number implied by cDNA position must match the AA number.
                cdna_pos = int(cdna_match.group(1))
                aa_number = int(aa_match.group(2))
                self.assertEqual(
                    (cdna_pos - 1) // 3 + 1, aa_number,
                    f"{name}: {hgvs} does not correspond to codon {aa_number} in {aa}",
                )

    @unittest.skipUnless(PICKETT_BED.is_file(), "pickett checkout not available on this machine")
    def test_matches_picketts_authoritative_list_exactly(self):
        pickett_rows = {
            (fields[0], int(fields[2]), fields[4], fields[5])
            for fields in read_bed(PICKETT_BED)
        }
        ours = {
            (fields[0], int(fields[2]), fields[4], fields[5])
            for fields in read_bed(RESOURCES / "branford_positions_hg19.bed")
        }
        self.assertEqual(pickett_rows, ours)

    def test_ref_alleles_match_a_real_hg19_reference_if_available(self):
        candidates = [
            Path.home() / "ref/genomes/hg19_pickett/ref_genome.fa",
        ]
        fasta = next((path for path in candidates if path.is_file()), None)
        if fasta is None:
            self.skipTest("no local hg19 FASTA available")
        script = (
            HERE.parent.parent.parent
            / "bcr_abl_af_eval/workflow/scripts/branford_bed_to_vcf.py"
        )
        if not script.is_file():
            self.skipTest("bcr_abl_af_eval checkout not available")
        result = subprocess.run(
            [
                "python3", str(script),
                "--bed", str(RESOURCES / "branford_positions_hg19.bed"),
                "--reference", str(fasta),
                "--out-prefix", "/tmp/test_bcr_abl_ref_check",
            ],
            capture_output=True, text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("All 73 REF alleles match", result.stdout)


class FusionJunctionTableTest(unittest.TestCase):
    """Locks the coordinate table bcr_abl_fusion_counter.py actually reads at
    runtime -- bcr_abl_features.json is not consulted for these values, so a
    fix there alone would not have changed any fusion call."""

    def test_e13_and_e14_are_not_swapped(self):
        """e13 must be upstream of (smaller genomic coordinate than) e14 for
        both BCR (+ strand): this is the exact swap that made every e14a2
        read report as e13a2."""

        for build in ("hg19", "hg38"):
            exons = fusion_counter.KNOWN_JUNCTIONS[build]["bcr_exons"]
            self.assertLess(
                exons["e13"], exons["e14"],
                f"{build}: e13 ({exons['e13']}) is not upstream of e14 ({exons['e14']})",
            )

    def test_abl_exon2_is_upstream_of_exon3(self):
        for build in ("hg19", "hg38"):
            exons = fusion_counter.KNOWN_JUNCTIONS[build]["abl_exons"]
            self.assertLess(exons["a2"], exons["a3"])

    def test_junction_positions_fall_inside_the_matching_breakpoint_cluster(self):
        """Sanity bound: a donor/acceptor coordinate that drifted by 140 kb
        (the previous a2 error) lands outside every known breakpoint region."""

        # (chrom, start, end) from rindels' independently curated clusters,
        # cross-checked earlier against pickett's exon BED.
        clusters_hg19 = {
            "e1": ("chr22", 23_523_145, 23_596_187),
            "e13": ("chr22", 23_631_678, 23_634_851),
            "e14": ("chr22", 23_631_678, 23_634_851),
            "a2": ("chr9", 133_574_969, 133_730_488),
            "a3": ("chr9", 133_574_969, 133_730_488),
        }
        hg19 = fusion_counter.KNOWN_JUNCTIONS["hg19"]
        for exon, (chrom, start, end) in clusters_hg19.items():
            pos = hg19["bcr_exons"].get(exon) or hg19["abl_exons"].get(exon)
            self.assertIsNotNone(pos)
            self.assertTrue(
                start <= pos <= end,
                f"{exon} at {pos} falls outside its breakpoint cluster {chrom}:{start}-{end}",
            )

    def test_hg38_offsets_match_the_branford_bed_offset(self):
        """The chr9 offset used here must be the same one that produced the
        hg38 Branford BED, or ABL1-side fusion calls and Branford calls would
        disagree about where the gene is."""

        hg19_a2 = fusion_counter.KNOWN_JUNCTIONS["hg19"]["abl_exons"]["a2"]
        hg38_a2 = fusion_counter.KNOWN_JUNCTIONS["hg38"]["abl_exons"]["a2"]
        self.assertEqual(hg19_a2 - hg38_a2, 2875387)

        hg19_e14 = fusion_counter.KNOWN_JUNCTIONS["hg19"]["bcr_exons"]["e14"]
        hg38_e14 = fusion_counter.KNOWN_JUNCTIONS["hg38"]["bcr_exons"]["e14"]
        self.assertEqual(hg19_e14 - hg38_e14, 342187)


if __name__ == "__main__":
    unittest.main()
