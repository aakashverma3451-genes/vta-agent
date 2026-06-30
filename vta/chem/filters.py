"""Annotation-only medicinal chemistry filters."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any


@dataclass(frozen=True)
class ChemistryFlags:
    pains: list[str]
    brenk: list[str]
    aggregator: bool
    ro5_violations: int | None
    beyond_ro5: bool | None
    synthetic_accessibility: float | None
    method: str

    def as_dict(self) -> dict[str, Any]:
        return {
            "pains": self.pains,
            "brenk": self.brenk,
            "aggregator": self.aggregator,
            "ro5_violations": self.ro5_violations,
            "beyond_ro5": self.beyond_ro5,
            "synthetic_accessibility": self.synthetic_accessibility,
            "method": self.method,
        }


def _fallback_pains(smiles: str) -> list[str]:
    alerts = []
    if "C(=S)S" in smiles or "NC(=S)S" in smiles:
        alerts.append("rhodanine-like")
    if smiles.count("O") >= 2 and "c1" in smiles and "O)" in smiles:
        alerts.append("polyphenol/catechol-like")
    return alerts


def _fallback_sa(heavy_atoms: int, rings: int) -> float:
    return round(min(10.0, 1.5 + heavy_atoms / 12.0 + rings * 0.35), 2)


def chemistry_flags(smiles: str | None) -> dict[str, Any]:
    """Return annotation-only chemistry flags for a SMILES string."""
    if not smiles:
        return ChemistryFlags([], [], False, None, None, None, "no-smiles").as_dict()
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski
        from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    except Exception:
        return ChemistryFlags(
            _fallback_pains(smiles), [], False, None, None, None,
            "fallback-no-rdkit").as_dict()

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ChemistryFlags([], [], False, None, None, None, "invalid-smiles").as_dict()

    pains: list[str] = _fallback_pains(smiles)
    brenk: list[str] = []
    try:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
        for match in FilterCatalog(params).GetMatches(mol):
            pains.append(match.GetDescription())
    except Exception:
        pass
    try:
        params = FilterCatalogParams()
        params.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)
        for match in FilterCatalog(params).GetMatches(mol):
            brenk.append(match.GetDescription())
    except Exception:
        pass

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    rot = Lipinski.NumRotatableBonds(mol)
    heavy = mol.GetNumHeavyAtoms()
    rings = Lipinski.RingCount(mol)
    ro5 = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    beyond_ro5 = mw > 700 or logp > 7 or hbd > 8 or hba > 15 or rot > 15
    aggregator = logp > 5.5 and mw > 450 and rings >= 3
    sa = _fallback_sa(heavy, rings)
    return ChemistryFlags(
        sorted(set(pains)), sorted(set(brenk)), aggregator, int(ro5),
        bool(beyond_ro5), sa, "rdkit+heuristics").as_dict()
