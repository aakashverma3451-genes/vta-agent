"""Runtime-detected ADMET provider contracts."""
from __future__ import annotations

from abc import ABC, abstractmethod


class ADMETProvider(ABC):
    name: str

    @abstractmethod
    def available(self) -> bool:
        ...

    @abstractmethod
    def predict(self, smiles: list[str]) -> list[dict]:
        ...


class ADMETAIProvider(ADMETProvider):
    name = "ADMET-AI"

    def available(self) -> bool:
        try:
            import admet_ai  # noqa: F401
            return True
        except Exception:
            return False

    def predict(self, smiles: list[str]) -> list[dict]:
        from admet_ai import ADMETModel
        model = ADMETModel()
        preds = model.predict(smiles)
        rows = []
        for i, smi in enumerate(smiles):
            rows.append({
                "herg": _row_value(preds, i, "hERG"),
                "oral": _row_value(preds, i, "Bioavailability_Ma"),
                "solubility": _row_value(preds, i, "Solubility_AqSolDB"),
                "applicability_domain": "unknown",
                "provider": self.name,
                "smiles": smi,
            })
        return rows


def _row_value(preds, i: int, key: str):
    try:
        import pandas as pd
        if isinstance(preds, pd.DataFrame):
            value = preds.iloc[i][key]
            return round(float(value), 3)
    except Exception:
        pass
    if isinstance(preds, dict) and key in preds:
        return round(float(preds[key]), 3)
    return None


def available_providers() -> list[ADMETProvider]:
    return [p for p in [ADMETAIProvider()] if p.available()]
