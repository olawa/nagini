"""Check the clinical corrections and explicit resource override contract."""

import json
import tempfile
import unittest
from pathlib import Path

from bcr_abl_annotations import load_features, load_tki_profiles


class AnnotationResourceTest(unittest.TestCase):
    def test_eln_2025_asciminib_resistance_has_source_evidence(self):
        features = load_features()
        profiles = load_tki_profiles()
        for mutation in ("F359V", "M244V", "Q252H"):
            self.assertEqual(profiles[mutation]["Asciminib"], "Resistant")
            source = features["tki_profile_evidence"][mutation]["Asciminib"]["source"]
            self.assertEqual(source, "https://doi.org/10.1038/s41375-025-02664-w")

    def test_partial_override_keeps_other_drugs_and_does_not_mutate_defaults(self):
        original = load_tki_profiles()
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "custom.json"
            path.write_text(json.dumps({"tki_profiles": {"T315I": {"Asciminib": "Unknown"}}}))
            custom = load_tki_profiles(path)
        self.assertEqual(custom["T315I"]["Asciminib"], "Unknown")
        self.assertEqual(custom["T315I"]["Ponatinib"], original["T315I"]["Ponatinib"])
        self.assertEqual(load_tki_profiles(), original)

    def test_broken_or_invalid_explicit_override_is_not_silently_ignored(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "custom.json"
            with self.assertRaises(FileNotFoundError):
                load_tki_profiles(path)
            path.write_text("not JSON")
            with self.assertRaises(ValueError):
                load_tki_profiles(path)
            path.write_text(json.dumps({"tki_profiles": {"T315I": {"Asciminib": "Sensitiv"}}}))
            with self.assertRaises(ValueError):
                load_tki_profiles(path)


if __name__ == "__main__":
    unittest.main()
