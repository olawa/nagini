"""Tests for sample routing.

These import the module the workflow actually uses, so a change in
nagini_routing.py is what they check -- not a copy of it.
"""

import importlib.util
import unittest
from pathlib import Path

import pandas as pd
import yaml


HERE = Path(__file__).parent
REPO = HERE.parents[1]

_spec = importlib.util.spec_from_file_location("nagini_routing", HERE / "nagini_routing.py")
routing = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(routing)


def units_frame(rows):
    return pd.DataFrame(rows)


class AssayRegistryTest(unittest.TestCase):
    def test_builtin_assays_resolve_with_defaults_filled_in(self):
        assays = routing.resolve_assays({})
        self.assertLessEqual(
            {"wgs", "adaptive", "amplicon", "cdna", "bcr_abl"}, set(assays)
        )
        for name, definition in assays.items():
            self.assertEqual(
                set(definition), set(routing.ASSAY_DEFAULTS), f"{name} has odd keys"
            )

    def test_bcr_abl_is_spliced_and_variant_called(self):
        """The combination the old single-column routing could not express."""

        assays = routing.resolve_assays({})
        bcr_abl = assays["bcr_abl"]
        self.assertEqual(bcr_abl["library_type"], "cdna")
        self.assertTrue(bcr_abl["snv_calling"])
        self.assertFalse(bcr_abl["sv_calling"])
        self.assertEqual(routing.minimap2_preset(bcr_abl["library_type"]), "splice")

    def test_adaptive_is_wgs_confined_to_the_design(self):
        assays = routing.resolve_assays({})
        adaptive = assays["adaptive"]
        self.assertEqual(adaptive["library_type"], "wgs")
        self.assertTrue(adaptive["snv_calling"])
        self.assertTrue(adaptive["sv_calling"])
        self.assertTrue(adaptive["target_restricted"])
        self.assertEqual(routing.minimap2_preset(adaptive["library_type"]), "map-ont")

    def test_a_new_assay_is_a_config_edit(self):
        assays = routing.resolve_assays(
            {
                "assays": {
                    "my_panel": {
                        "library_type": "amplicon",
                        "snv_calling": True,
                        "target_restricted": True,
                        "modules": ["my_panel"],
                    }
                }
            }
        )
        self.assertIn("my_panel", assays)
        self.assertFalse(assays["my_panel"]["sv_calling"])  # default applied
        self.assertEqual(assays["my_panel"]["modules"], ["my_panel"])

    def test_config_can_override_a_builtin_without_replacing_it(self):
        assays = routing.resolve_assays({"assays": {"bcr_abl": {"sv_calling": True}}})
        self.assertTrue(assays["bcr_abl"]["sv_calling"])
        self.assertEqual(assays["bcr_abl"]["library_type"], "cdna")
        self.assertEqual(assays["bcr_abl"]["modules"], ["bcr_abl"])

    def test_unsupported_keys_and_library_types_are_rejected(self):
        with self.assertRaises(routing.RoutingError) as raised:
            routing.resolve_assays({"assays": {"x": {"snv_callng": True}}})
        self.assertIn("unsupported key", str(raised.exception))

        with self.assertRaises(routing.RoutingError) as raised:
            routing.resolve_assays({"assays": {"x": {"library_type": "rna"}}})
        self.assertIn("library_type", str(raised.exception))


class SampleRoutingTest(unittest.TestCase):
    def setUp(self):
        self.assays = routing.resolve_assays({})

    def test_assay_column_decides_analysis(self):
        units = units_frame([
            {"sample": "S_WGS", "assay": "wgs"},
            {"sample": "S_ADAPT", "assay": "adaptive"},
            {"sample": "S_BCR", "assay": "bcr_abl"},
        ])
        self.assertTrue(routing.sample_needs_sv_calling(units, "S_WGS", self.assays))
        self.assertTrue(routing.sample_needs_sv_calling(units, "S_ADAPT", self.assays))
        self.assertFalse(routing.sample_needs_sv_calling(units, "S_BCR", self.assays))
        self.assertTrue(routing.sample_needs_snv_calling(units, "S_BCR", self.assays))
        self.assertEqual(
            routing.get_library_type(units, "S_BCR", self.assays), "cdna"
        )

    def test_modules_gate_assay_specific_analysis(self):
        units = units_frame([
            {"sample": "S_BCR", "assay": "bcr_abl"},
            {"sample": "S_WGS", "assay": "wgs"},
        ])
        self.assertTrue(routing.sample_has_module(units, "S_BCR", self.assays, "bcr_abl"))
        self.assertFalse(routing.sample_has_module(units, "S_WGS", self.assays, "bcr_abl"))
        self.assertEqual(routing.sample_modules(units, "S_WGS", self.assays), [])

    def test_units_library_type_overrides_the_assay_default(self):
        """A BCR-ABL panel run on gDNA must not be spliced."""

        units = units_frame([
            {"sample": "S", "assay": "bcr_abl", "library_type": "amplicon"},
        ])
        self.assertEqual(routing.get_library_type(units, "S", self.assays), "amplicon")
        self.assertEqual(
            routing.minimap2_preset(routing.get_library_type(units, "S", self.assays)),
            "map-ont",
        )
        self.assertTrue(routing.sample_needs_snv_calling(units, "S", self.assays))
        warnings = routing.routing_warnings(units, ["S"], self.assays)
        self.assertEqual(len(warnings), 1)
        self.assertIn("library_type=amplicon", warnings[0])

    def test_unknown_assay_names_the_valid_ones(self):
        units = units_frame([{"sample": "S", "assay": "bcrabl"}])
        with self.assertRaises(routing.RoutingError) as raised:
            routing.get_assay(units, "S", self.assays)
        message = str(raised.exception)
        self.assertIn("not defined", message)
        self.assertIn("bcr_abl", message)

    def test_conflicting_values_within_a_sample_are_rejected(self):
        units = units_frame([
            {"sample": "S", "assay": "wgs"},
            {"sample": "S", "assay": "bcr_abl"},
        ])
        with self.assertRaises(routing.RoutingError):
            routing.get_assay(units, "S", self.assays)


class BackwardCompatibilityTest(unittest.TestCase):
    """Sample sheets written before the assay column must route unchanged."""

    def setUp(self):
        self.assays = routing.resolve_assays({})

    def test_no_columns_at_all_defaults_to_wgs(self):
        units = units_frame([{"sample": "S", "type": "T"}])
        self.assertEqual(routing.get_assay(units, "S", self.assays), "wgs")
        self.assertTrue(routing.sample_needs_snv_calling(units, "S", self.assays))
        self.assertTrue(routing.sample_needs_sv_calling(units, "S", self.assays))

    def test_legacy_library_type_values_keep_their_old_routing(self):
        units = units_frame([
            {"sample": "S_WGS", "library_type": "wgs"},
            {"sample": "S_AMP", "library_type": "amplicon"},
            {"sample": "S_CDNA", "library_type": "cdna"},
        ])
        # wgs: SNV + SV, unspliced
        self.assertTrue(routing.sample_needs_snv_calling(units, "S_WGS", self.assays))
        self.assertTrue(routing.sample_needs_sv_calling(units, "S_WGS", self.assays))
        # amplicon: SNV only
        self.assertTrue(routing.sample_needs_snv_calling(units, "S_AMP", self.assays))
        self.assertFalse(routing.sample_needs_sv_calling(units, "S_AMP", self.assays))
        # cdna without an assay: spliced, and no calling, exactly as before
        self.assertFalse(routing.sample_needs_snv_calling(units, "S_CDNA", self.assays))
        self.assertFalse(routing.sample_needs_sv_calling(units, "S_CDNA", self.assays))
        self.assertEqual(
            routing.minimap2_preset(
                routing.get_library_type(units, "S_CDNA", self.assays)
            ),
            "splice",
        )

    def test_empty_and_missing_cells_are_treated_as_absent(self):
        units = units_frame([
            {"sample": "S1", "assay": "", "library_type": ""},
            {"sample": "S2", "assay": "NA", "library_type": "nan"},
        ])
        for sample in ("S1", "S2"):
            self.assertEqual(routing.get_assay(units, sample, self.assays), "wgs")


class PresetTest(unittest.TestCase):
    def test_only_rna_libraries_are_spliced(self):
        self.assertEqual(routing.minimap2_preset("cdna"), "splice")
        self.assertEqual(routing.minimap2_preset("direct_rna"), "splice")
        self.assertEqual(routing.minimap2_preset("wgs"), "map-ont")
        self.assertEqual(routing.minimap2_preset("amplicon"), "map-ont")

    def test_the_configured_default_is_used_for_unspliced_libraries(self):
        self.assertEqual(routing.minimap2_preset("wgs", "lr:hq"), "lr:hq")
        self.assertEqual(routing.minimap2_preset("cdna", "lr:hq"), "splice")


class OutputConditionTest(unittest.TestCase):
    def setUp(self):
        self.assays = routing.resolve_assays({})

    def test_every_condition_in_the_shipped_output_spec_is_known(self):
        """Guards against a typo in `when:` silently dropping an output."""

        output_spec = yaml.safe_load((REPO / "config/output_files.yaml").read_text())
        self.assertEqual(routing.unknown_conditions(output_spec, self.assays), [])

    def test_a_typo_is_reported(self):
        spec = {"files": [{"name": "x", "output": "x", "when": "sv_callng"}]}
        self.assertEqual(routing.unknown_conditions(spec, self.assays), ["sv_callng"])

    def test_assay_scoped_conditions_are_accepted(self):
        self.assertTrue(routing.is_known_condition("assay_module:bcr_abl", self.assays))
        self.assertTrue(routing.is_known_condition("assay:bcr_abl", self.assays))
        self.assertFalse(routing.is_known_condition("assay:nonexistent", self.assays))
        self.assertFalse(routing.is_known_condition("assay_module:", self.assays))


class SummaryTest(unittest.TestCase):
    def test_summary_reports_what_will_run(self):
        assays = routing.resolve_assays({})
        units = units_frame([
            {"sample": "S_BCR", "assay": "bcr_abl"},
            {"sample": "S_WGS", "assay": "wgs"},
        ])
        summary = {row["sample"]: row for row in routing.routing_summary(units, ["S_BCR", "S_WGS"], assays)}
        self.assertEqual(summary["S_BCR"]["preset"], "splice")
        self.assertEqual(summary["S_BCR"]["modules"], ["bcr_abl"])
        self.assertTrue(summary["S_WGS"]["sv_calling"])
        self.assertEqual(summary["S_WGS"]["preset"], "map-ont")


class AssaySettingTest(unittest.TestCase):
    def setUp(self):
        self.assays = routing.resolve_assays({})

    def test_targeted_assays_lower_the_allele_fraction_floor(self):
        """A 20 000x amplicon and a 30x genome cannot share a floor."""

        units = units_frame([
            {"sample": "S_WGS", "assay": "wgs"},
            {"sample": "S_BCR", "assay": "bcr_abl"},
        ])
        self.assertEqual(
            routing.assay_setting(units, "S_WGS", self.assays, "clairs_to_call", "snv_min_af", 0.05),
            0.05,
        )
        self.assertEqual(
            routing.assay_setting(units, "S_BCR", self.assays, "clairs_to_call", "snv_min_af", 0.05),
            0.005,
        )

    def test_a_setting_without_an_override_falls_back_to_the_default(self):
        units = units_frame([{"sample": "S", "assay": "bcr_abl"}])
        self.assertEqual(
            routing.assay_setting(units, "S", self.assays, "clairs_to_call", "extra", "-x"), "-x"
        )
        self.assertIsNone(
            routing.assay_setting(units, "S", self.assays, "nope", "nope")
        )

    def test_config_can_add_overrides_to_a_new_assay(self):
        assays = routing.resolve_assays(
            {"assays": {"panel": {"library_type": "amplicon",
                                  "config_overrides": {"clairs_to_call": {"snv_min_af": 0.02}}}}}
        )
        units = units_frame([{"sample": "S", "assay": "panel"}])
        self.assertEqual(
            routing.assay_setting(units, "S", assays, "clairs_to_call", "snv_min_af", 0.05), 0.02
        )

    def test_a_malformed_override_is_rejected(self):
        with self.assertRaises(routing.RoutingError):
            routing.resolve_assays({"assays": {"x": {"config_overrides": {"clairs_to_call": 0.01}}}})


class ClairsToPlatformTest(unittest.TestCase):
    def test_known_models_map_to_their_error_model(self):
        for model, expected in [
            ("dna_r10.4.1_e8.2_400bps_sup@v5.2.0", "ont_r10_dorado_sup_5khz"),
            ("dna_r10.4.1_e8.2_400bps_hac@v5.0.0", "ont_r10_dorado_hac_5khz"),
            ("dna_r10.4.1_e8.2_400bps_sup@v4.1.0", "ont_r10_dorado_sup_4khz"),
            ("dna_r10.4.1_e8.2_400bps_hac@v4.1.0", "ont_r10_dorado_hac_4khz"),
        ]:
            platform, warning = routing.clairs_to_platform(model)
            self.assertEqual(platform, expected, model)
            self.assertIsNone(warning, model)

    def test_an_explicit_sampling_token_wins_over_the_version(self):
        platform, warning = routing.clairs_to_platform("dna_r10.4.1_e8.2_5khz_400bps_sup@v4.0.0")
        self.assertEqual(platform, "ont_r10_dorado_sup_5khz")
        self.assertIsNone(warning)

    def test_an_unknown_model_warns_and_falls_back_to_sup_5khz(self):
        for model in ("", "sup", "hac", "dna_r9.4.1_e8_sup@v3.6.0", "something_new@v9.9.9"):
            platform, warning = routing.clairs_to_platform(model)
            self.assertEqual(platform, "ont_r10_dorado_sup_5khz", model)
            self.assertIsNotNone(warning, model)
            self.assertIn("ont_r10_dorado_sup_5khz", warning)

    def test_the_fallback_is_configurable(self):
        platform, warning = routing.clairs_to_platform("mystery", default="ont_r10_dorado_hac_4khz")
        self.assertEqual(platform, "ont_r10_dorado_hac_4khz")
        self.assertIn("ont_r10_dorado_hac_4khz", warning)

    def test_the_table_can_be_extended_from_config(self):
        platform, warning = routing.clairs_to_platform(
            "dna_r9.4.1_e8_sup@v3.6.0", {"r9,sup,4khz": "ont_r9_guppy"}
        )
        self.assertEqual(platform, "ont_r9_guppy")
        self.assertIsNone(warning)

    def test_model_parsing_reports_what_it_read(self):
        parsed = routing.parse_basecalling_model("dna_r10.4.1_e8.2_400bps_sup@v5.2.0")
        self.assertEqual(parsed["chemistry"], "r10")
        self.assertEqual(parsed["accuracy"], "sup")
        self.assertEqual(parsed["sampling"], "5khz")
        self.assertEqual(parsed["version"], "5.2.0")

    def test_the_model_is_read_from_units(self):
        units = units_frame([
            {"sample": "S", "basecalling_model": "dna_r10.4.1_e8.2_400bps_sup@v5.2.0"},
        ])
        self.assertEqual(
            routing.get_basecalling_model(units, "S"),
            "dna_r10.4.1_e8.2_400bps_sup@v5.2.0",
        )
        self.assertEqual(routing.get_basecalling_model(units_frame([{"sample": "S"}]), "S"), "")


class BarcodeKitTest(unittest.TestCase):
    def test_the_count_is_read_from_the_kit_name(self):
        self.assertEqual(routing.barcode_kit_size("SQK-NBD114-24"), 24)
        self.assertEqual(routing.barcode_kit_size("SQK-RBK114-96"), 96)
        self.assertEqual(routing.barcode_kit_size("sqk-nbd114-24"), 24)

    def test_kits_that_do_not_follow_the_convention_come_from_the_table(self):
        self.assertEqual(routing.barcode_kit_size("EXP-NBD104"), 12)
        self.assertEqual(routing.barcode_kit_size("SQK-RBK004"), 12)

    def test_a_chemistry_version_is_not_read_as_a_count(self):
        self.assertIsNone(routing.barcode_kit_size("SQK-LSK114"))
        self.assertIsNone(routing.barcode_kit_size("SQK-NBD114"))
        self.assertIsNone(routing.barcode_kit_size(""))

    def test_barcodes_beyond_the_kit_are_reported(self):
        outside = routing.barcodes_outside_kit(
            ["barcode01", "barcode24", "barcode30"], "SQK-NBD114-24"
        )
        self.assertEqual([name for name, _ in outside], ["barcode30"])

    def test_an_unknown_kit_checks_nothing(self):
        self.assertEqual(routing.barcodes_outside_kit(["barcode30"], "SQK-LSK114"), [])

    def test_the_table_can_be_extended_from_config(self):
        self.assertEqual(routing.barcode_kit_size("MY-KIT", {"MY-KIT": 8}), 8)
        self.assertEqual(
            [name for name, _ in routing.barcodes_outside_kit(["barcode09"], "MY-KIT", {"MY-KIT": 8})],
            ["barcode09"],
        )


class CallerResolutionTest(unittest.TestCase):
    def setUp(self):
        self.assays = routing.resolve_assays({})

    def test_auto_routes_normal_to_deepvariant_and_tumor_to_clairs_to(self):
        units = units_frame([
            {"sample": "S_NORM", "assay": "wgs", "type": "N"},
            {"sample": "S_TUMOR", "assay": "wgs", "type": "T"},
        ])
        self.assertEqual(
            routing.get_sample_caller(units, "S_NORM", self.assays, "auto", "N"),
            "deepvariant",
        )
        self.assertEqual(
            routing.get_sample_caller(units, "S_TUMOR", self.assays, "auto", "T"),
            "clairs_to",
        )

    def test_amplicon_and_bcr_abl_route_to_clairs_to(self):
        units = units_frame([
            {"sample": "S_AMP", "assay": "amplicon"},
            {"sample": "S_BCR", "assay": "bcr_abl"},
        ])
        self.assertEqual(
            routing.get_sample_caller(units, "S_AMP", self.assays, "auto"),
            "clairs_to",
        )
        self.assertEqual(
            routing.get_sample_caller(units, "S_BCR", self.assays, "auto"),
            "clairs_to",
        )

    def test_explicit_global_caller_overrides_auto(self):
        units = units_frame([
            {"sample": "S_WGS", "assay": "wgs"},
        ])
        self.assertEqual(
            routing.get_sample_caller(units, "S_WGS", self.assays, "deepvariant", "T"),
            "deepvariant",
        )
        self.assertEqual(
            routing.get_sample_caller(units, "S_WGS", self.assays, "clairs_to", "N"),
            "clairs_to",
        )

    def test_assay_override_takes_precedence(self):
        assays = routing.resolve_assays({
            "assays": {
                "custom_wgs": {
                    "library_type": "wgs",
                    "snv_calling": True,
                    "config_overrides": {"snv_indels": {"small_variant_caller": "mutect2"}},
                }
            }
        })
        units = units_frame([{"sample": "S_CUST", "assay": "custom_wgs"}])
        self.assertEqual(
            routing.get_sample_caller(units, "S_CUST", assays, "auto", "N"),
            "mutect2",
        )


if __name__ == "__main__":
    unittest.main()


