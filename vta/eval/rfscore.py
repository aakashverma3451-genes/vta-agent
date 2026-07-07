"""RF-Score-style learned rescorer (Phase S) — a transparent, in-environment learned scorer.

Off-the-shelf pretrained rescorers (GNINA, RTMScore, ODDT) are un-installable on this arm64
stack (Docker/CDN, dgl wheels, OpenBabel-3.x incompatibility). This is the honest alternative:
RF-Score's own, published method (Ballester & Mitchell 2010, Bioinformatics 26:1169) built on
the stack that works here (openbabel read + numpy + scikit-learn) so we can still ask — with a
learned scorer — whether it beats the trivial 2D-QSAR on activity cliffs.

RF-Score v1 features: for a docked protein–ligand pose, count contacts between each
(protein element, ligand element) pair within a distance cutoff (12 Å). Protein ∈ {C,N,O,S},
ligand ∈ {C,N,O,F,P,S,Cl,Br,I} → 36 features. A RandomForest maps those to pIC50.

CRUCIAL leakage control: this model is TARGET-SPECIFIC (trained on Moonshot Mpro), so it is
evaluated ONLY out-of-fold under **scaffold-clustered** leave-out CV (compounds sharing a
Bemis–Murcko scaffold never straddle the train/test split), and always reported against the
same 2D-QSAR baseline on the identical cliff pairs. A win is a *target-specific* result, not a
general-transfer claim (that needs a PDBbind-pretrained scorer on a native box).
"""
from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

PROT_ELEMS = ("C", "N", "O", "S")
LIG_ELEMS = ("C", "N", "O", "F", "P", "S", "Cl", "Br", "I")
CUTOFF = 12.0
FEATURE_NAMES = [f"{pe}-{le}" for pe in PROT_ELEMS for le in LIG_ELEMS]

# AutoDock (PDBQT) atom type → element. H and pseudo/glue atoms (G, G0, CG… — macrocycle
# closure markers openbabel chokes on) map to None and are skipped. This direct parse is what
# fixes the coverage bug where openbabel dropped whole poses on a single 'G' record.
_AD2ELEM = {"C": "C", "A": "C", "N": "N", "NA": "N", "NS": "N", "OA": "O", "OS": "O", "O": "O",
            "SA": "S", "S": "S", "F": "F", "P": "P", "CL": "Cl", "BR": "Br", "I": "I",
            "HD": None, "H": None, "HS": None, "G": None, "GA": None, "G0": None, "G1": None,
            "CG": None, "CG0": None, "CG1": None, "Z": None, "W": None}
_ELEM_KEEP = set(PROT_ELEMS) | set(LIG_ELEMS)


def _elem_from_pdb(line: str) -> Optional[str]:
    el = line[76:78].strip().capitalize()          # standard PDB element column
    if not el:
        name = line[12:16].strip()                 # fallback: infer from atom name
        el = (name[:2].capitalize() if name[:2] in ("CL", "BR") else name[:1].upper()) if name else ""
    return el if el in _ELEM_KEEP else None


def _elem_from_pdbqt(line: str) -> Optional[str]:
    ad = line.rstrip().split()[-1].upper()          # AutoDock type is the last token
    return _AD2ELEM.get(ad)


def read_atoms(path: str, fmt: str) -> Dict[str, np.ndarray]:
    """Parse a PDB/PDBQT → {element: (N,3) coords} for the RF-Score element sets.

    Direct fixed-column parse (no openbabel): robust to AutoDock pseudo-atoms and H, so a single
    unusual record never drops the whole structure. Elements outside the RF-Score sets are ignored.
    """
    pick = _elem_from_pdbqt if fmt == "pdbqt" else _elem_from_pdb
    by_elem: Dict[str, List[List[float]]] = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            el = pick(line)
            if el is None:
                continue
            try:
                xyz = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            except ValueError:
                continue
            by_elem.setdefault(el, []).append(xyz)
    return {el: np.asarray(v, dtype=float) for el, v in by_elem.items()}


def rfscore_features(receptor_by_elem: Dict[str, np.ndarray], ligand_path: str,
                     ligand_fmt: str = "pdbqt", cutoff: float = CUTOFF) -> Optional[np.ndarray]:
    """36-D RF-Score v1 contact-count feature vector for one docked pose (None if unreadable)."""
    try:
        lig = read_atoms(ligand_path, ligand_fmt)
    except Exception:
        return None
    feats = np.zeros(len(PROT_ELEMS) * len(LIG_ELEMS), dtype=float)
    k = 0
    for pe in PROT_ELEMS:
        pcoords = receptor_by_elem.get(pe)
        for le in LIG_ELEMS:
            lcoords = lig.get(le)
            if pcoords is not None and pcoords.size and lcoords is not None and lcoords.size:
                # pairwise distances protein(pe) × ligand(le); count within cutoff
                d = np.sqrt(((pcoords[:, None, :] - lcoords[None, :, :]) ** 2).sum(-1))
                feats[k] = float((d <= cutoff).sum())
            k += 1
    return feats


def murcko_scaffold(smiles: str) -> str:
    """Generic Bemis–Murcko scaffold SMILES (for leave-scaffold-out grouping); '' on failure."""
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
        m = Chem.MolFromSmiles(smiles)
        if m is None:
            return ""
        gen = MurckoScaffold.MakeScaffoldGeneric(MurckoScaffold.GetScaffoldForMol(m))
        return Chem.MolToSmiles(gen)
    except Exception:
        return ""


def scaffold_cv_predict(X: np.ndarray, y: np.ndarray, groups: List[str], *,
                        n_splits: int = 5, seed: int = 0) -> np.ndarray:
    """Out-of-fold RandomForest pIC50 predictions under scaffold-grouped K-fold CV.

    Compounds sharing a scaffold group are never split across train/test, so a cliff pair's
    members are predicted by a model that never saw their scaffold — no analog leakage. Returns
    an out-of-fold prediction per row (NaN if a row was never in a usable test fold).
    """
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.model_selection import GroupKFold

    uniq = sorted(set(groups))
    n_splits = max(2, min(n_splits, len(uniq)))
    preds = np.full(len(y), np.nan, dtype=float)
    gkf = GroupKFold(n_splits=n_splits)
    for tr, te in gkf.split(X, y, groups):
        rf = RandomForestRegressor(n_estimators=200, random_state=seed, n_jobs=-1)
        rf.fit(X[tr], y[tr])
        preds[te] = rf.predict(X[te])
    return preds
