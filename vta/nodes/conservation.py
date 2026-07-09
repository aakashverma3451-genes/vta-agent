"""vta.nodes.conservation — real per-pocket conservation (retire the 0.5 placeholder).

SPEC #1. `pockets.py` sets every pocket's `conservation` to a neutral 0.5 because real
conservation needs a multiple-sequence alignment. This node computes it for real, so
`rank.py`'s EXISTING 10% conservation term carries signal: a pocket whose residues are
evolutionarily conserved across viral homologs is harder to escape by mutation — a
genuine antiviral druggability axis.

Metric: Jensen–Shannon divergence per MSA column vs a background AA distribution
(Capra & Singh 2007, Bioinformatics 23(15):1875, doi:10.1093/bioinformatics/btm270) —
the best-performing conservation measure for functional residues and natively bounded
to [0,1]. A column dominated by one residue diverges far from background (→1); a column
that looks like background diverges little (→0). Columns with many gaps are down-weighted.

Pocket residues come by PROXIMITY: the experimental site is only a `center`, so we take
residues with any atom within R=6 Å of it. To avoid the numbering hazard (experimental
chains like 8PSO are not guaranteed 1-indexed/contiguous), we map by OBSERVED RESIDUE
ORDER in the structure, not by raw residue numbers, and derive the target sequence from
the structure itself (the MSA's row 0). No MSA → keep 0.5 and label it honestly.
"""
from __future__ import annotations

import math

from vta.state import VTAState

R_PROXIMITY = 6.0                       # Å — pocket residue = any atom within this of center
_CONSERVATION_FALLBACK = 0.5            # the neutral placeholder we replace on the real path

_AA = "ACDEFGHIKLMNPQRSTVWY"

# Background AA frequencies (Robinson & Robinson 1991, as used by Capra & Singh's ref impl).
_BACKGROUND = {
    "A": 0.078, "R": 0.051, "N": 0.045, "D": 0.054, "C": 0.019, "Q": 0.043,
    "E": 0.063, "G": 0.074, "H": 0.022, "I": 0.051, "L": 0.090, "K": 0.057,
    "M": 0.022, "F": 0.039, "P": 0.052, "S": 0.071, "T": 0.058, "W": 0.013,
    "Y": 0.032, "V": 0.064,
}

_THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


# ── injectable seam (monkeypatched in tests / stubbed to None in CI) ──────────
def fetch_homolog_msa(sequence: str, taxon: str | None = None) -> list[str] | None:
    """Gapped homolog MSA rows (row 0 = target `sequence`), or None if unavailable.

    Real path: a committed alignment built offline via homolog search + MAFFT
    (`vta.data.msa`), so the gate is reproducible without network. Returning None makes
    the node fall back to the labelled 0.5 placeholder (graceful-fallback discipline).
    """
    from vta.data.msa import load_msa_for
    return load_msa_for(sequence, taxon)


# ── pure: JSD conservation (Capra & Singh 2007) ──────────────────────────────
def _norm_background() -> list[float]:
    tot = sum(_BACKGROUND[a] for a in _AA)
    return [_BACKGROUND[a] / tot for a in _AA]


def _kl(p: list[float], q: list[float]) -> float:
    """KL(p||q) in bits; 0*log0 := 0, and q has no zeros (background-derived mixture)."""
    s = 0.0
    for pi, qi in zip(p, q):
        if pi > 0.0:
            s += pi * math.log2(pi / qi)
    return s


def jensen_shannon_conservation(columns: list[str]) -> list[float]:
    """Per-column conservation in [0,1] = JSD(column AA dist, background), gap-weighted.

    `columns[i]` is the i-th alignment column as a string across all rows (may contain
    '-' / 'X'). Identical column → ~1.0; background-like column → ~0.0; gappy columns
    are down-weighted by their non-gap fraction.
    """
    q = _norm_background()
    out: list[float] = []
    for col in columns:
        observed = [c for c in col if c in _BACKGROUND]
        if not observed:
            out.append(0.0)
            continue
        counts = {a: 0 for a in _AA}
        for c in observed:
            counts[c] += 1
        p = [counts[a] / len(observed) for a in _AA]
        m = [(pi + qi) / 2.0 for pi, qi in zip(p, q)]
        jsd = 0.5 * _kl(p, m) + 0.5 * _kl(q, m)          # log2 → bounded by 1.0
        gap_weight = len(observed) / len(col)            # Capra gap down-weighting
        out.append(round(max(0.0, min(1.0, jsd * gap_weight)), 4))
    return out


# ── pure: structure parsing + proximity mapping ──────────────────────────────
def parse_residues(pdb_text: str) -> list[dict]:
    """Ordered residues observed in the PDB: [{key, aa, atoms:[(x,y,z)]}], first model.

    Residues are returned in observed order (NOT by residue number) keyed by
    (chain, resseq, icode), so downstream mapping is index-based and immune to gapped
    or non-1-indexed numbering.
    """
    residues: list[dict] = []
    index: dict[tuple, dict] = {}
    for line in pdb_text.splitlines():
        if line.startswith("ENDMDL"):
            break                                        # first model only
        if not line.startswith(("ATOM", "HETATM")):
            continue
        resname = line[17:20].strip()
        aa = _THREE_TO_ONE.get(resname)
        if aa is None:
            continue                                     # skip ligands/waters/non-standard
        key = (line[21], line[22:26].strip(), line[26].strip())   # chain, resseq, icode
        try:
            xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
        except ValueError:
            continue
        rec = index.get(key)
        if rec is None:
            rec = {"key": key, "aa": aa, "atoms": []}
            index[key] = rec
            residues.append(rec)
        rec["atoms"].append(xyz)
    return residues


def structure_sequence(residues: list[dict]) -> str:
    """One-letter target sequence in observed order (the MSA's row-0 reference)."""
    return "".join(r["aa"] for r in residues)


def pocket_residue_indices(residues: list[dict], center, radius: float = R_PROXIMITY) -> list[int]:
    """Observed-order indices of residues with any atom within `radius` Å of `center`."""
    cx, cy, cz = center
    r2 = radius * radius
    hits: list[int] = []
    for i, res in enumerate(residues):
        for (x, y, z) in res["atoms"]:
            if (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= r2:
                hits.append(i)
                break
    return hits


def _row0_column_map(row0: str) -> list[int]:
    """Map observed-residue order-index → alignment column, via row 0's non-gap chars."""
    return [col for col, ch in enumerate(row0) if ch not in ("-", ".")]


# ── node ─────────────────────────────────────────────────────────────────────
def conservation_node(state: VTAState) -> VTAState:
    structures = state.get("structures") or {}
    pockets = state.get("pockets") or {}
    taxon = (state.get("classification") or {}).get("taxon")

    for name, pocket_list in pockets.items():
        struct = structures.get(name) or {}
        pdb_path = struct.get("pdb_path")
        if not pdb_path:
            continue
        try:
            with open(pdb_path) as fh:
                residues = parse_residues(fh.read())
        except OSError:
            residues = []
        if not residues:
            state["audit_trail"].append(
                f"[skip] conservation[{name}]: no parseable residues; kept 0.5")
            continue

        target = structure_sequence(residues)
        msa = fetch_homolog_msa(target, taxon)
        if not msa or msa[0].replace("-", "").replace(".", "") != target:
            reason = "no MSA" if not msa else "MSA row0 != structure sequence"
            state["audit_trail"].append(
                f"[skip] conservation[{name}]: {reason}; kept 0.5 placeholder")
            continue

        ncols = len(msa[0])
        cols = ["".join(row[c] for row in msa if c < len(row)) for c in range(ncols)]
        jsd = jensen_shannon_conservation(cols)
        col_of = _row0_column_map(msa[0])               # order-index → alignment column

        # Persist the per-residue JSD keyed by resseq (SPEC #4): conservation_contacts
        # reads this to make conservation ligand-specific. Keyed by resseq so a docked
        # pose's receptor (same structure PDB) can look residues up by number.
        res_map = {
            residues[i]["key"][1]: jsd[col_of[i]]
            for i in range(min(len(residues), len(col_of)))
        }
        state.setdefault("residue_conservation", {})[name] = res_map

        for pocket in pocket_list:
            idxs = pocket_residue_indices(residues, pocket["center"])
            scores = [jsd[col_of[i]] for i in idxs if i < len(col_of)]
            if scores:
                pocket["conservation"] = round(sum(scores) / len(scores), 3)
                pocket.setdefault("detectors", [])
                state["audit_trail"].append(
                    f"conservation[{name}] pocket {pocket.get('id')}: "
                    f"{pocket['conservation']} (JSD over {len(scores)} residues "
                    f"within {R_PROXIMITY}Å, {len(msa)} homologs)")
            else:
                state["audit_trail"].append(
                    f"[skip] conservation[{name}] pocket {pocket.get('id')}: "
                    f"no residues within {R_PROXIMITY}Å; kept 0.5")

    state["versions"]["conservation"] = (
        "JSD/Capra&Singh-2007 vs Robinson background, proximity 6Å, MAFFT MSA")
    return state
