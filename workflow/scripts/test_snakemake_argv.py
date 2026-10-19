"""Each script's snakemake_argv must produce arguments its own parser accepts,
carrying the rule's values to the right options."""

import importlib.util
import unittest
from pathlib import Path
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent


def load(name):
    spec = importlib.util.spec_from_file_location(name, HERE / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def smk(**sections):
    base = dict(wildcards=SimpleNamespace(sample="S1", type="T"), threads=4)
    base.update({k: SimpleNamespace(**v) if isinstance(v, dict) else v for k, v in sections.items()})
    return SimpleNamespace(**base)


class BranfordTest(unittest.TestCase):
    def setUp(self):
        self.mod = load("bcr_abl_branford_caller")
        self.params = dict(
            min_af=0.01, subclonal_af=0.002, min_depth=50, min_alt_reads=5, max_depth=1000000,
            min_base_quality=12, min_mapping_quality=1, min_snr=3.0, min_minor_strand_ratio=0.05,
            min_mean_baseq=14.0, min_kmer_rate=0.6, max_sor=3.0, kmer_window=3, extra="",
        )

    def job(self, ref="ref.fa", known=(), extra=""):
        params = dict(self.params, extra=extra)
        return smk(
            input=dict(bam="a.bam", bed="p.bed", features="f.json", ref_fasta=ref, known_variants=list(known)),
            output=dict(all_sites="all.tsv", mutations="mut.tsv", vcf="m.vcf"),
            params=params,
        )

    def test_arguments_parse_and_carry_the_values(self):
        args = self.mod.parse_args(self.mod.snakemake_argv(self.job(known=["kv.vcf"])))
        self.assertEqual(str(args.bam), "a.bam")
        self.assertEqual(args.sample, "S1")
        self.assertEqual(args.min_af, 0.01)
        self.assertEqual(args.min_mapping_quality, 1)
        self.assertEqual(args.max_sor, 3.0)
        self.assertEqual(args.threads, 4)
        self.assertEqual(str(args.known_variants), "kv.vcf")
        self.assertEqual(str(args.reference_fasta), "ref.fa")

    def test_optional_inputs_are_left_out_when_absent(self):
        argv = self.mod.snakemake_argv(self.job(ref="", known=[]))
        self.assertNotIn("--reference-fasta", argv)
        self.assertNotIn("--known-variants", argv)
        self.mod.parse_args(argv)

    def test_extra_is_split_like_a_shell_would(self):
        argv = self.mod.snakemake_argv(self.job(extra="--min-depth 99"))
        self.assertEqual(argv[argv.index("--kmer-window") + 2:argv.index("--threads")], ["--min-depth", "99"])


class PhaserTest(unittest.TestCase):
    def test_arguments(self):
        mod = load("bcr_abl_phaser")
        argv = mod.snakemake_argv(smk(
            input=dict(bam="a.bam", mutations="m.tsv"),
            output=dict(tsv="p.tsv", json="p.json"),
            params=dict(min_spanning_reads=10, min_mapq=1, min_base_quality=12),
        ))
        args = mod.parse_args(argv)
        self.assertEqual((args.sample, args.min_spanning_reads, args.min_mapq, args.min_base_quality, args.threads),
                         ("S1", 10, 1, 12, 4))
        self.assertEqual(str(args.mutations_tsv), "m.tsv")


class ExonCoverageTest(unittest.TestCase):
    def test_arguments(self):
        mod = load("bcr_abl_exon_coverage")
        args = mod.parse_args(mod.snakemake_argv(smk(
            input=dict(bam="a.bam", bed="e.bed"), output=dict(tsv="c.tsv"),
            params=dict(min_mapq=1, max_depth=1000000),
        )))
        self.assertEqual((args.sample, args.min_mapq, args.max_depth, args.threads), ("S1", 1, 1000000, 4))


class FusionCounterTest(unittest.TestCase):
    def test_arguments(self):
        mod = load("bcr_abl_fusion_counter")
        args = mod.parse_args(mod.snakemake_argv(smk(
            input=dict(bam="a.bam", features="f.json"),
            output=dict(tsv="f.tsv", json="f.json", arriba="arriba.tsv"),
            params=dict(tolerance=150, min_mapq=1, genome_build="hg19"),
        )))
        self.assertEqual((args.sample, args.genome_build, args.tolerance, args.min_mapq), ("S1", "hg19", 150, 1))
        self.assertEqual(str(args.arriba_output), "arriba.tsv")


class DemuxQcTest(unittest.TestCase):
    def test_arguments(self):
        mod = load("dorado_demux_summary")
        args = mod.parse_args(mod.snakemake_argv(smk(
            input=dict(summary="s.tsv", units="u.tsv"), output=dict(mqc="q.tsv"), params=dict(kit="SQK-NBD114-24"),
        )))
        self.assertEqual(args.kit, "SQK-NBD114-24")


if __name__ == "__main__":
    unittest.main()
