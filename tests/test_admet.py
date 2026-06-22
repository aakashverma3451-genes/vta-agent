"""ADMET node tests — hermetic (fake model, no torch load)."""
from __future__ import annotations

import pytest

import vta.nodes.admet as admet
from vta.state import new_state


class _FakeModel:
    def predict(self, smiles):
        import pandas as pd
        return pd.DataFrame([
            {"hERG": 0.02, "Bioavailability_Ma": 0.70, "Solubility_AqSolDB": -0.80}
            for _ in smiles
        ])


def _state():
    st = new_state("x", "x")
    st["lead_candidates"] = [
        {"ligand": "Ribavirin", "smiles": "C1=NC...", "score": 0.6, "positive_control": True},
    ]
    return st


def test_admet_annotates_leads(monkeypatch):
    monkeypatch.setattr(admet, "_admet_available", lambda: True)
    admet_ai = pytest.importorskip("admet_ai")
    monkeypatch.setattr(admet_ai, "ADMETModel", _FakeModel)

    out = admet.admet_node(_state())
    a = out["lead_candidates"][0]["admet"]
    assert a["herg"] == 0.02 and a["oral"] == 0.70 and a["solubility"] == -0.80
    assert out["lead_candidates"][0]["admet_flag"] == ""        # 0.02 < 0.5 → no risk
    assert any("ADMET:" in l for l in out["audit_trail"])


def test_admet_flags_herg_risk(monkeypatch):
    monkeypatch.setattr(admet, "_admet_available", lambda: True)
    admet_ai = pytest.importorskip("admet_ai")

    class Risky(_FakeModel):
        def predict(self, smiles):
            df = super().predict(smiles); df["hERG"] = 0.9; return df
    monkeypatch.setattr(admet_ai, "ADMETModel", Risky)

    out = admet.admet_node(_state())
    assert out["lead_candidates"][0]["admet_flag"] == "hERG risk"


def test_admet_skips_gracefully_when_unavailable(monkeypatch):
    monkeypatch.setattr(admet, "_admet_available", lambda: False)
    out = admet.admet_node(_state())
    assert out["lead_candidates"][0].get("admet") is None
    assert any("[skip] ADMET" in l for l in out["audit_trail"])
