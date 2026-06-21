"""structure_node (REAL) — produce a 3D structure for each extracted protein.

Plan Task 2.1. Decision order per protein, in priority order:

  1. EXPERIMENTAL-FIRST. If the subunit has a solved cryo-EM chain, fetch that PDB
     and extract THAT CHAIN (not the whole complex — see snag #2). A real 2.4 Å
     structure beats any prediction. Experimental structures get mean_plddt = None,
     meaning "trusted ground truth", NOT "low/missing" — Phase-2 docking must treat
     None as trusted (snag #3).

  2. ESMFold FALLBACK. For a protein with no experimental chain AND length within the
     public API's ceiling, fold via api.esmatlas.com and compute mean pLDDT.

  3. REFUSE (honestly). The api.esmatlas.com endpoint has a HARD 400-residue ceiling
     (probed: 413 "Sequence is longer than 400" at 500 aa; 504 timeout at 400 aa).
     PB1 (~500) / PB2 (~450) exceed it. Rather than crash, we record pdb_path=None
     with a clear reason so the graph proceeds and Phase 2 (local GPU ESMFold or
     chunking) picks it up.

Chain-mapping note (from the actual 8PSO/8PT2 COMPND records): chain A = PA-like,
chain B = "putative PB1", chain C = RdRp. There is NO chain explicitly labelled
PB2, so PB2 has no confident experimental structure here and falls to the ESMFold
path (which then refuses, being >400 aa). This nomenclature gap is real and flagged,
not hidden.

Network I/O goes through `fetch_rcsb_pdb` / `fold_esmfold` so tests can monkeypatch
them and stay offline.
"""
from __future__ import annotations

import os

import requests

from vta.state import VTAState

# Public ESMFold API ceiling — hard limit confirmed by probe (413 above 400).
ESMFOLD_API_MAX_AA = 400

# Confident subunit → (PDB id, chain) from the 8PSO COMPND records. Only entries we
# can defend biologically are listed; PB2 is intentionally absent (no labelled chain).
EXPERIMENTAL_PDB = {
    "PA":  ("8PSO", "A"),   # COMPND: "POLYMERASE ACIDIC PROTEIN (PA-LIKE)" — 316 res
    "PB1": ("8PSO", "B"),   # COMPND: "PUTATIVE PB1" — 515 res (matches ~500 aa ORF)
}

_OUT_DIR = "structures"


# ── injectable network seams (monkeypatched in tests) ────────────────────────
def fetch_rcsb_pdb(pdb_id: str) -> str:
    """Download a PDB file from RCSB. Raises on HTTP error."""
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60)
    r.raise_for_status()
    return r.text


def fold_esmfold(sequence: str) -> str:
    """Fold a single sequence via the public ESMFold API. Raises on HTTP error."""
    r = requests.post(
        "https://api.esmatlas.com/foldSequence/v1/pdb/",
        data=sequence,
        headers={"Content-Type": "text/plain"},
        timeout=300,
    )
    r.raise_for_status()
    return r.text


# ── pure helpers ─────────────────────────────────────────────────────────────
def extract_chain(pdb_text: str, chain: str) -> str:
    """Keep only the atom records for `chain` — the fix for whole-complex reuse."""
    keep = []
    for line in pdb_text.splitlines():
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM", "TER", "ANISOU") and len(line) > 21 and line[21] == chain:
            keep.append(line)
    keep.append("END")
    return "\n".join(keep) + "\n"


def mean_plddt(pdb_text: str) -> float:
    """ESMFold writes per-residue pLDDT into the B-factor column (CA atoms)."""
    vals = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                pass
    return round(sum(vals) / len(vals), 1) if vals else 0.0


# ── node ─────────────────────────────────────────────────────────────────────
def structure_node(state: VTAState) -> VTAState:
    os.makedirs(_OUT_DIR, exist_ok=True)
    structures: dict[str, dict] = {}

    for name, prot in (state.get("extracted_proteins") or {}).items():
        seq = prot.get("sequence", "")
        n = len(seq)
        path = os.path.join(_OUT_DIR, f"{state['run_id']}_{name}.pdb")
        rec = {"pdb_path": None, "mean_plddt": None, "method": None, "source": None}

        if name in EXPERIMENTAL_PDB:
            pdb_id, chain = EXPERIMENTAL_PDB[name]
            try:
                chain_pdb = extract_chain(fetch_rcsb_pdb(pdb_id), chain)
                with open(path, "w") as fh:
                    fh.write(chain_pdb)
                # experimental = trusted ground truth → mean_plddt stays None (snag #3)
                rec.update(pdb_path=path, method="experimental",
                           source=f"{pdb_id}:{chain}")
                state["audit_trail"].append(
                    f"Structure[{name}]: EXPERIMENTAL {pdb_id} chain {chain} "
                    f"(pLDDT=None = trusted ground truth)")
            except Exception as e:  # network/HTTP failure — degrade, don't crash
                rec["method"] = "experimental_failed"
                state["audit_trail"].append(
                    f"Structure[{name}]: WARNING experimental fetch failed ({e}); "
                    f"no structure")

        elif n <= ESMFOLD_API_MAX_AA:
            try:
                pdb_text = fold_esmfold(seq)
                plddt = mean_plddt(pdb_text)
                with open(path, "w") as fh:
                    fh.write(pdb_text)
                flag = " (LOW CONFIDENCE)" if plddt < 70 else ""
                rec.update(pdb_path=path, mean_plddt=plddt, method="esmfold",
                           source="api.esmatlas.com")
                state["audit_trail"].append(
                    f"Structure[{name}]: ESMFold, mean pLDDT {plddt}{flag}")
            except Exception as e:
                rec["method"] = "esmfold_failed"
                state["audit_trail"].append(
                    f"Structure[{name}]: WARNING ESMFold failed ({e}); no structure")

        else:
            # > API ceiling and no experimental chain → refuse honestly.
            rec["method"] = "refused_too_long"
            state["audit_trail"].append(
                f"Structure[{name}]: REFUSED — {n} aa > ESMFold API ceiling "
                f"({ESMFOLD_API_MAX_AA}); needs local ESMFold/chunking (Phase 2)")

        structures[name] = rec

    state["structures"] = structures
    state["versions"]["esmfold"] = "esm2 (api.esmatlas.com, <=400aa) + experimental PDB"
    return state
