"""Provider interface for ADMET prediction backends."""
from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass


@dataclass(frozen=True)
class ADMETPrediction:
    smiles: str
    endpoints: dict[str, float | None]
    provider: str
    warnings: tuple[str, ...] = ()


class ADMETProvider(ABC):
    """Abstract base class for ADMET providers."""

    name: str

    @abstractmethod
    def available(self) -> bool:
        """Return True when the provider can run in this environment."""

    @abstractmethod
    def predict(self, smiles: list[str]) -> list[ADMETPrediction]:
        """Predict ADMET endpoints for a list of SMILES strings."""


class ADMETAIProvider(ADMETProvider):
    """Adapter contract for the optional admet-ai package."""

    name = "admet-ai"

    def available(self) -> bool:
        try:
            import admet_ai  # noqa: F401
            return True
        except Exception:
            return False

    def predict(self, smiles: list[str]) -> list[ADMETPrediction]:
        if not self.available():
            return [
                ADMETPrediction(s, {}, self.name, ("admet-ai not installed",))
                for s in smiles
            ]
        from admet_ai import ADMETModel

        model = ADMETModel()
        raw = model.predict(smiles)
        out: list[ADMETPrediction] = []
        for i, value in enumerate(smiles):
            endpoints: dict[str, float | None] = {}
            for endpoint in ("hERG", "Bioavailability_Ma", "Solubility_AqSolDB"):
                try:
                    endpoints[endpoint] = float(raw.iloc[i][endpoint])
                except Exception:
                    endpoints[endpoint] = None
            out.append(ADMETPrediction(value, endpoints, self.name))
        return out
