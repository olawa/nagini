import unittest
from pathlib import Path
import pandas as pd


class LibraryRoutingTest(unittest.TestCase):
    def setUp(self):
        self.output_spec = {
            "directory": ".",
            "files": [
                {
                    "name": "aligned BAM",
                    "output": "alignment/nagini_bam/{sample}_{type}.bam",
                },
                {
                    "name": "Cramino QC",
                    "output": "qc/cramino/{sample}_{type}.txt",
                },
                {
                    "name": "whole genome coverage summary",
                    "output": "qc/mosdepth/{sample}_{type}.mosdepth.summary.txt",
                },
                {
                    "name": "target region coverage summary",
                    "output": "qc/mosdepth_bed/{sample}_{type}.mosdepth.summary.txt",
                    "when": "has_design_bed",
                },
                {
                    "name": "Sniffles2 SV VCF",
                    "output": "cnv_sv/sniffles2_call/{sample}_{type}.vcf.gz",
                    "when": "sv_calling",
                },
                {
                    "name": "ClairS-TO normalized VCF",
                    "output": "snv_indels/clairs_to_postprocess/{sample}_{type}.filtered.norm.vcf.gz",
                    "when": "clairs_to_enabled",
                },
                {
                    "name": "Nagini Excel report",
                    "output": "reports/{sample}_{type}.nagini_report.xlsx",
                },
            ],
        }

        self.samples = pd.DataFrame([
            {"sample": "SAMPLE_WGS"},
            {"sample": "SAMPLE_AMP"},
            {"sample": "SAMPLE_RNA"},
        ]).set_index("sample", drop=False)

        self.units = pd.DataFrame([
            {
                "sample": "SAMPLE_WGS",
                "type": "T",
                "library_type": "wgs",
                "processing_unit": "PU001",
                "barcode": "BC",
            },
            {
                "sample": "SAMPLE_AMP",
                "type": "T",
                "library_type": "amplicon",
                "processing_unit": "PU002",
                "barcode": "BC",
            },
            {
                "sample": "SAMPLE_RNA",
                "type": "T",
                "library_type": "cdna",
                "processing_unit": "PU003",
                "barcode": "BC",
            },
        ])

    def test_routing_functions(self):
        units = self.units

        def get_library_type(sample):
            col = "library_type"
            if col not in units.columns:
                return "wgs"
            types = sorted(
                {
                    str(val).strip().lower()
                    for val in units.loc[units["sample"] == sample, col].tolist()
                    if str(val).strip() and str(val).strip().lower() not in ("nan", "none", "")
                }
            )
            return types[0] if types else "wgs"

        def sample_needs_dna_variant_calling(sample):
            return get_library_type(sample) in ("wgs", "amplicon")

        def sample_needs_sv_calling(sample):
            return get_library_type(sample) == "wgs"

        def sample_needs_rna_analysis(sample):
            return get_library_type(sample) == "cdna"

        self.assertEqual(get_library_type("SAMPLE_WGS"), "wgs")
        self.assertEqual(get_library_type("SAMPLE_AMP"), "amplicon")
        self.assertEqual(get_library_type("SAMPLE_RNA"), "cdna")

        self.assertTrue(sample_needs_dna_variant_calling("SAMPLE_WGS"))
        self.assertTrue(sample_needs_dna_variant_calling("SAMPLE_AMP"))
        self.assertFalse(sample_needs_dna_variant_calling("SAMPLE_RNA"))

        self.assertTrue(sample_needs_sv_calling("SAMPLE_WGS"))
        self.assertFalse(sample_needs_sv_calling("SAMPLE_AMP"))
        self.assertFalse(sample_needs_sv_calling("SAMPLE_RNA"))

        self.assertFalse(sample_needs_rna_analysis("SAMPLE_WGS"))
        self.assertFalse(sample_needs_rna_analysis("SAMPLE_AMP"))
        self.assertTrue(sample_needs_rna_analysis("SAMPLE_RNA"))

    def test_output_file_compilation_per_sample(self):
        units = self.units
        samples = self.samples
        output_spec = self.output_spec

        def get_library_type(sample):
            col = "library_type"
            if col not in units.columns:
                return "wgs"
            types = sorted(
                {
                    str(val).strip().lower()
                    for val in units.loc[units["sample"] == sample, col].tolist()
                    if str(val).strip() and str(val).strip().lower() not in ("nan", "none", "")
                }
            )
            return types[0] if types else "wgs"

        def sample_needs_dna_variant_calling(sample):
            return get_library_type(sample) in ("wgs", "amplicon")

        def sample_needs_sv_calling(sample):
            return get_library_type(sample) == "wgs"

        def sample_needs_rna_analysis(sample):
            return get_library_type(sample) == "cdna"

        def has_design_bed():
            return False

        def clairs_to_enabled():
            return True

        def get_unit_types_for_sample(sample):
            return sorted(set(units.loc[units["sample"] == sample, "type"]))

        outdir = Path(".")
        output_files = []
        for file_definition in output_spec["files"]:
            condition = file_definition.get("when")
            allowed_types = file_definition.get("types")

            if condition == "clairs_to_enabled" and not clairs_to_enabled():
                continue
            elif condition == "has_design_bed" and not has_design_bed():
                continue

            for sample in samples.index:
                clean_sample = str(sample).strip()

                if condition in ("clairs_to_enabled", "sample_needs_dna_variant_calling"):
                    if not sample_needs_dna_variant_calling(clean_sample):
                        continue
                elif condition in ("sv_calling", "sample_needs_sv_calling"):
                    if not sample_needs_sv_calling(clean_sample):
                        continue
                elif condition == "sample_needs_rna_analysis":
                    if not sample_needs_rna_analysis(clean_sample):
                        continue

                for unit_type in get_unit_types_for_sample(sample):
                    if allowed_types and unit_type not in allowed_types:
                        continue
                    clean_type = str(unit_type).strip()
                    template = str(file_definition["output"]).strip()
                    output_files.append(
                        str(
                            outdir
                            / template.format(
                                sample=clean_sample,
                                type=clean_type,
                            )
                        )
                    )

        # WGS sample has all standard outputs
        self.assertIn("alignment/nagini_bam/SAMPLE_WGS_T.bam", output_files)
        self.assertIn("qc/cramino/SAMPLE_WGS_T.txt", output_files)
        self.assertIn("qc/mosdepth/SAMPLE_WGS_T.mosdepth.summary.txt", output_files)
        self.assertIn("cnv_sv/sniffles2_call/SAMPLE_WGS_T.vcf.gz", output_files)
        self.assertIn("snv_indels/clairs_to_postprocess/SAMPLE_WGS_T.filtered.norm.vcf.gz", output_files)
        self.assertIn("reports/SAMPLE_WGS_T.nagini_report.xlsx", output_files)

        # Amplicon sample has SNVs, BAM, Cramino, report, but NO SVs
        self.assertIn("alignment/nagini_bam/SAMPLE_AMP_T.bam", output_files)
        self.assertIn("qc/cramino/SAMPLE_AMP_T.txt", output_files)
        self.assertIn("qc/mosdepth/SAMPLE_AMP_T.mosdepth.summary.txt", output_files)
        self.assertNotIn("cnv_sv/sniffles2_call/SAMPLE_AMP_T.vcf.gz", output_files)
        self.assertIn("snv_indels/clairs_to_postprocess/SAMPLE_AMP_T.filtered.norm.vcf.gz", output_files)
        self.assertIn("reports/SAMPLE_AMP_T.nagini_report.xlsx", output_files)

        # cDNA sample has BAM, Cramino, report, but NO SNVs and NO SVs
        self.assertIn("alignment/nagini_bam/SAMPLE_RNA_T.bam", output_files)
        self.assertIn("qc/cramino/SAMPLE_RNA_T.txt", output_files)
        self.assertIn("qc/mosdepth/SAMPLE_RNA_T.mosdepth.summary.txt", output_files)
        self.assertNotIn("cnv_sv/sniffles2_call/SAMPLE_RNA_T.vcf.gz", output_files)
        self.assertNotIn("snv_indels/clairs_to_postprocess/SAMPLE_RNA_T.filtered.norm.vcf.gz", output_files)
        self.assertIn("reports/SAMPLE_RNA_T.nagini_report.xlsx", output_files)

        # target coverage excluded because has_design_bed is False
        self.assertNotIn("qc/mosdepth_bed/SAMPLE_WGS_T.mosdepth.summary.txt", output_files)


if __name__ == "__main__":
    unittest.main()
