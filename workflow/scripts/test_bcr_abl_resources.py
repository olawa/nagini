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

    def test_matches_picketts_authoritative_list_exactly(self):
        try:
            if not PICKETT_BED.is_file():
                self.skipTest("pickett checkout not available on this machine")
            pickett_rows = {
                (fields[0], int(fields[2]), fields[4], fields[5])
                for fields in read_bed(PICKETT_BED)
            }
        except (PermissionError, OSError):
            self.skipTest("pickett checkout not accessible")
        ours = {
            (fields[0], int(fields[2]), fields[4], fields[5])
            for fields in read_bed(RESOURCES / "branford_positions_hg19.bed")
        }
        self.assertEqual(pickett_rows, ours)

    def test_ref_alleles_match_a_real_hg19_reference_if_available(self):
        try:
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
            if result.returncode != 0 and "Operation not permitted" in result.stderr:
                self.skipTest("sandbox permission prevented running external script")
            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("All 73 REF alleles match", result.stdout)
        except (PermissionError, OSError):
            self.skipTest("sandbox permission prevented running external script")



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

    # 1-based ends of BCR exons in Ensembl GRCh37 (ENST00000305877), 2026-10-26.
    ENSEMBL_BCR_DONORS = {"e1": 23524426, "e6": 23613779, "e13": 23631808, "e14": 23632600, "e19": 23654023}

    def test_bcr_donors_are_the_ensembl_exon_ends(self):
        """e6 stood in intron 1, 18 kb from exon 6, and e19 was 36 bp short."""
        self.assertEqual(fusion_counter.KNOWN_JUNCTIONS["hg19"]["bcr_exons"], self.ENSEMBL_BCR_DONORS)

    def test_every_bcr_exon_region_contains_its_donor(self):
        for exon, (start, end) in fusion_counter._BCR_EXON_REGIONS_HG19.items():
            donor = fusion_counter.KNOWN_JUNCTIONS["hg19"]["bcr_exons"][exon]
            self.assertTrue(start <= donor <= end, f"{exon}: donor {donor} outside {start}-{end}")

    def test_features_json_bcr_donors_agree_with_the_counter(self):
        import json
        features = json.loads((RESOURCES / "bcr_abl_features.json").read_text())
        seen = 0

        def walk(node):
            nonlocal seen
            if isinstance(node, dict):
                exon = node.get("bcr_exon")
                for build in ("hg19", "hg38"):
                    key = f"bcr_ref_pos_{build}"
                    if exon and key in node:
                        seen += 1
                        self.assertEqual(node[key], fusion_counter.KNOWN_JUNCTIONS[build]["bcr_exons"][exon], f"{exon} {key}")
                for value in node.values():
                    walk(value)
            elif isinstance(node, list):
                for value in node:
                    walk(value)

        walk(features)
        self.assertGreater(seen, 0)

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


class ExonBedTest(unittest.TestCase):

    def test_hg19_and_hg38_exons_have_same_number_of_rows(self):
        hg19 = read_bed(RESOURCES / "bcr_abl1_exons_hg19.bed")
        hg38 = read_bed(RESOURCES / "bcr_abl1_exons_hg38.bed")
        self.assertEqual(len(hg19), len(hg38))
        self.assertEqual(len(hg19), 15)

    def test_hg38_exons_are_shifted_by_exact_offsets(self):
        hg19 = read_bed(RESOURCES / "bcr_abl1_exons_hg19.bed")
        hg38 = read_bed(RESOURCES / "bcr_abl1_exons_hg38.bed")
        for row19, row38 in zip(hg19, hg38):
            chrom = row19[0]
            offset = 2875387 if chrom == "chr9" else 342187
            self.assertEqual(int(row19[1]) - int(row38[1]), offset)
            self.assertEqual(int(row19[2]) - int(row38[2]), offset)
            self.assertEqual(row19[3], row38[3])

    # Verified against the hg19 sequence on 2026-10-22 (UCSC and Ensembl
    # GRCh37 agree): AG immediately before each of these exon starts, ATG at
    # the 1b and BCR exon 1 starts. The file came from pickett with 1-based
    # starts in the 0-based column, so every exon began one base late.
    VERIFIED_STARTS_0BASED = {
        "ABL1_Exon1_NM_007313.2": 133589706,   # ATG
        "ABL1_Exon2_NM_005157.6": 133729450,   # AG at 133729448-449
        "ABL1_Exon5_NM_005157.6": 133747515,   # AG at 133747513-514; c.823
        "BCR_Exon1_NM_004327.4": 23523147,     # ATG
        "BCR_Exon13_NM_004327.4": 23631703,
        "BCR_Exon14_NM_004327.4": 23632525,
    }

    def test_exon_starts_are_zero_based(self):
        starts = {row[3]: int(row[1]) for row in read_bed(RESOURCES / "bcr_abl1_exons_hg19.bed")}
        for name, start in self.VERIFIED_STARTS_0BASED.items():
            self.assertEqual(starts[name], start, f"{name}: BED start must be 0-based")

    def test_codon_276_phase_matches_the_exon_5_boundary(self):
        """c.826 (first base of codon 276) is 133747518 0-based in the Branford
        BED. Exon 4 ends at c.822, so exon 5 begins at c.823 -- three bases
        earlier, 133747515. A one-base error in either file breaks this."""

        exon5 = next(row for row in read_bed(RESOURCES / "bcr_abl1_exons_hg19.bed")
                     if row[3] == "ABL1_Exon5_NM_005157.6")
        d276n = next(row for row in read_bed(RESOURCES / "branford_positions_hg19.bed")
                     if row[3] == "NM_005157.6:c.(826G>A)")
        self.assertEqual(int(d276n[1]) - int(exon5[1]), 826 - 823)


class ExonBoundarySequenceTest(unittest.TestCase):
    """The check that found the one-base error, run against a local reference
    where there is one (the workstation has pickett's hg19)."""

    CANDIDATES = (Path.home() / "ref/genomes/hg19_pickett/ref_genome.fa",)

    def test_exon_boundaries_sit_on_canonical_splice_sites(self):
        fasta_path = next((p for p in self.CANDIDATES if p.is_file()), None)
        if fasta_path is None:
            self.skipTest("no local hg19 FASTA available")
        try:
            import pysam
        except ImportError:
            self.skipTest("pysam is required to read the reference")
        first = {"ABL1_Exon1_NM_007313.2", "ABL1_Exon1_NM_005157.6", "BCR_Exon1_NM_004327.4"}
        last = {"ABL1_Exon11_NM_005157.6"}
        with pysam.FastaFile(str(fasta_path)) as fasta:
            contig = lambda c: c if c in fasta.references else c.removeprefix("chr")  # noqa: E731
            for chrom, start, end, name in (r[:4] for r in read_bed(RESOURCES / "bcr_abl1_exons_hg19.bed")):
                start, end, c = int(start), int(end), contig(chrom)
                if name not in first:
                    self.assertEqual(fasta.fetch(c, start - 2, start).upper(), "AG", f"{name}: no acceptor before {start}")
                if name not in last:
                    self.assertIn(fasta.fetch(c, end, end + 2).upper(), ("GT", "GC"), f"{name}: no donor after {end}")


class JunctionSequenceTest(unittest.TestCase):
    """Every fusion junction the counter knows sits on a canonical splice site
    in the reference -- the check that found e6 in an intron."""

    def test_donors_and_acceptors_are_canonical(self):
        fasta_path = next((p for p in ExonBoundarySequenceTest.CANDIDATES if p.is_file()), None)
        if fasta_path is None:
            self.skipTest("no local hg19 FASTA available")
        try:
            import pysam
        except ImportError:
            self.skipTest("pysam is required to read the reference")
        junctions = fusion_counter.KNOWN_JUNCTIONS["hg19"]
        with pysam.FastaFile(str(fasta_path)) as fasta:
            c22 = "chr22" if "chr22" in fasta.references else "22"
            c9 = "chr9" if "chr9" in fasta.references else "9"
            for exon, end in junctions["bcr_exons"].items():           # 1-based last exon base
                self.assertIn(fasta.fetch(c22, end, end + 2).upper(), ("GT", "GC"), f"BCR {exon}: no donor after {end}")
            for exon, start in junctions["abl_exons"].items():         # 1-based first exon base
                self.assertEqual(fasta.fetch(c9, start - 3, start - 1).upper(), "AG", f"ABL1 {exon}: no acceptor before {start}")


class GtfAndJunctionConsistencyTest(unittest.TestCase):
    """The GTF, the fusion counter and the features JSON all restate the exon
    BED. They were each one base late at every exon start because each was
    derived from it; lock them to the corrected BED instead."""

    @staticmethod
    def _gtf_exons(build):
        rows = {}
        for line in (RESOURCES / f"bcr_abl_{build}.gtf").read_text().splitlines():
            if line.startswith("#") or not line.strip():
                continue
            f = line.split("\t")
            if f[2] == "exon":
                rows[(f[0], int(f[4]))] = int(f[3])
        return rows

    def test_gtf_exon_starts_are_the_bed_starts_plus_one(self):
        for build in ("hg19", "hg38"):
            gtf = self._gtf_exons(build)
            bed = read_bed(RESOURCES / f"bcr_abl1_exons_{build}.bed")
            matched = 0
            for chrom, start, end, name in (r[:4] for r in bed):
                if (chrom, int(end)) in gtf:
                    matched += 1
                    self.assertEqual(gtf[(chrom, int(end))], int(start) + 1, f"{build} {name}")
            self.assertEqual(matched, 14, build)

    def test_fusion_acceptors_are_the_first_base_of_the_abl1_exon(self):
        bed = {r[3]: (int(r[1]), int(r[2])) for r in read_bed(RESOURCES / "bcr_abl1_exons_hg19.bed")}
        abl = fusion_counter.KNOWN_JUNCTIONS["hg19"]["abl_exons"]
        self.assertEqual(abl["a2"], bed["ABL1_Exon2_NM_005157.6"][0] + 1)
        self.assertEqual(abl["a3"], bed["ABL1_Exon3_NM_005157.6"][0] + 1)
        bcr = fusion_counter.KNOWN_JUNCTIONS["hg19"]["bcr_exons"]
        self.assertEqual(bcr["e13"], bed["BCR_Exon13_NM_004327.4"][1])
        self.assertEqual(bcr["e14"], bed["BCR_Exon14_NM_004327.4"][1])

    def test_features_json_acceptors_agree_with_the_counter(self):
        import json
        features = json.loads((RESOURCES / "bcr_abl_features.json").read_text())
        expected = {
            build: fusion_counter.KNOWN_JUNCTIONS[build]["abl_exons"] for build in ("hg19", "hg38")
        }
        seen = 0

        def walk(node):
            nonlocal seen
            if isinstance(node, dict):
                exon = node.get("abl_exon")
                for build in ("hg19", "hg38"):
                    key = f"abl_ref_pos_{build}"
                    if exon in ("a2", "a3") and key in node:
                        seen += 1
                        self.assertEqual(node[key], expected[build][exon], f"{exon} {key}")
                for value in node.values():
                    walk(value)
            elif isinstance(node, list):
                for value in node:
                    walk(value)

        walk(features)
        self.assertGreater(seen, 0)


if __name__ == "__main__":
    unittest.main()

