"""Load BCR-ABL1 annotations from the versioned bundled knowledge resource.

Clinical profiles are maintained in JSON, not duplicated in the caller.
Explicit overrides must be readable and valid; a broken override is not a
reason to silently fall back to a different set of clinical annotations.
"""

from __future__ import annotations

import json
from pathlib import Path


DEFAULT_FEATURES_PATH = Path(__file__).resolve().parents[1] / "resources" / "bcr_abl" / "bcr_abl_features.json"
TKI_NAMES = ("Imatinib", "Dasatinib", "Nilotinib", "Bosutinib", "Ponatinib", "Asciminib")
RESPONSE_VALUES = {"Sensitive", "Intermediate", "Reduced", "Resistant", "Unknown"}


def load_features(path=None):
    source = Path(path) if path is not None else DEFAULT_FEATURES_PATH
    with source.open() as handle:
        data = json.load(handle)
    if not isinstance(data, dict):
        raise ValueError(f"Expected an object in BCR-ABL1 features: {source}")
    profiles = data.get("tki_profiles", {})
    if not isinstance(profiles, dict):
        raise ValueError(f"tki_profiles must be an object: {source}")
    for mutation, responses in profiles.items():
        if not isinstance(responses, dict):
            raise ValueError(f"Invalid TKI profile for {mutation}: {source}")
        for drug, response in responses.items():
            if drug not in TKI_NAMES or response not in RESPONSE_VALUES:
                raise ValueError(f"Invalid TKI response {mutation}/{drug}={response!r}: {source}")
    return data


def load_tki_profiles(path=None):
    """Return independent profiles, merging an optional override per drug."""
    profiles = {key: dict(value) for key, value in load_features()["tki_profiles"].items()}
    if path is not None and Path(path).resolve() != DEFAULT_FEATURES_PATH:
        for mutation, responses in load_features(path).get("tki_profiles", {}).items():
            profiles.setdefault(mutation, {}).update(responses)
    return profiles


DEFAULT_TKI_PROFILES = load_tki_profiles()
