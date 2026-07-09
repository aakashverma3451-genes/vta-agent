"""Covalent-warhead classification for SARS-CoV-2 Mpro (Cys145) inhibitors.

Mpro has both **covalent** (electrophilic warhead reacting with the Cys145 thiol —
nirmatrelvir-class) and **non-covalent** inhibitors. AutoDock Vina scores non-covalent
binding only, so a covalent binder's affinity is not meaningful in a non-covalent dock.
This module classifies binding mode so the benchmark can keep the two classes separate —
the same binding-mode discipline that stopped HCV NI and NNI actives being pooled.

Precedence (experiment first):
  1. An explicit covalent annotation (e.g. COVID Moonshot ``covalent_warhead`` /
     ``Covalent Fragment``) is authoritative — trust it either way.
  2. Otherwise, substructure rules over the Cys145 electrophilic warheads decide; a match
     is conservatively called covalent (keeps reactive chemotypes out of the non-covalent
     headline set).
  3. Unparseable SMILES with no annotation → ``ambiguous`` (never silently pooled).

Pure + dependency-light: RDKit is used when present, but the module imports without it and
annotation-only classification still works offline.
"""
from __future__ import annotations

from typing import Optional

# SMARTS for the electrophilic warheads that react with the Mpro catalytic cysteine.
# Names are reported back so a reviewer can see *why* a compound was called covalent.
_WARHEAD_SMARTS: dict[str, str] = {
    "nitrile": "[CX2]#[NX1]",                       # reversible covalent (nirmatrelvir-like)
    "aldehyde": "[CX3H1]=O",
    "alpha_ketoamide": "[CX3](=O)[CX3](=O)[NX3]",
    "michael_acceptor": "[CX3]=[CX3][CX3]=[OX1]",   # acrylamide / vinyl ketone
    "haloacetamide": "[Cl,Br,I][CH2][CX3](=O)[NX3]",
    "vinyl_sulfone": "[CX3]=[CX3][SX4](=O)=O",
    "epoxide": "[OX2r3]1[CX4r3][CX4r3]1",
    "boronic_acid": "[BX3]([OX2H])[OX2H]",
}

# Truthy / falsey spellings seen in Moonshot + ChEMBL annotation columns.
_TRUE = {"true", "1", "yes", "y", "covalent", "t"}
_FALSE = {"false", "0", "no", "n", "non-covalent", "noncovalent", "f", ""}


def _annotation_to_bool(annotation: object) -> Optional[bool]:
    """Coerce a free-form covalent annotation to True/False, or None if uninformative."""
    if annotation is None:
        return None
    if isinstance(annotation, bool):
        return annotation
    text = str(annotation).strip().lower()
    if text in _TRUE:
        return True
    if text in _FALSE:
        return None if text == "" else False
    return None


def matched_warheads(smiles: str | None) -> list[str]:
    """Names of covalent warheads found in `smiles` (RDKit SMARTS), or []."""
    if not smiles:
        return []
    try:
        from rdkit import Chem
    except Exception:
        return []
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return []
    found = []
    for name, smarts in _WARHEAD_SMARTS.items():
        patt = Chem.MolFromSmarts(smarts)
        if patt is not None and mol.HasSubstructMatch(patt):
            found.append(name)
    return found


def classify_binding_mode(smiles: str | None, annotation: object = None) -> dict:
    """Classify a ligand as covalent / non_covalent / ambiguous.

    `annotation` is an optional experimental covalent flag (Moonshot ``covalent_warhead``
    or ``Covalent Fragment``); when informative it wins over substructure inference.

    Returns ``{"mode", "basis", "warheads", "annotation"}`` so the decision is auditable.
    """
    ann = _annotation_to_bool(annotation)
    if ann is True:
        return {"mode": "covalent", "basis": "annotation",
                "warheads": matched_warheads(smiles), "annotation": True}
    if ann is False:
        # Authoritative experimental non-covalent call — trust it (do not override on a
        # substructure match), but surface any warhead so a reviewer can see the tension.
        return {"mode": "non_covalent", "basis": "annotation",
                "warheads": matched_warheads(smiles), "annotation": False}

    warheads = matched_warheads(smiles)
    if warheads:
        return {"mode": "covalent", "basis": "substructure",
                "warheads": warheads, "annotation": None}
    # No annotation, no warhead match: non-covalent if it parsed, else ambiguous.
    try:
        from rdkit import Chem
        parsed = smiles is not None and Chem.MolFromSmiles(smiles) is not None
    except Exception:
        parsed = bool(smiles)
    if parsed:
        return {"mode": "non_covalent", "basis": "substructure",
                "warheads": [], "annotation": None}
    return {"mode": "ambiguous", "basis": "unparseable",
            "warheads": [], "annotation": None}
