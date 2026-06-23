"""Shared test fixtures — keep the structure node offline.

The real `structure_node` fetches experimental PDBs from RCSB and folds via the
ESMFold API. Tests monkeypatch those two network seams (`fetch_rcsb_pdb`,
`fold_esmfold`) so the graph runs end-to-end without a network.
"""
from __future__ import annotations

import pytest


def _atom(serial: int, chain: str, resseq: int, bfac: float = 50.0) -> str:
    """A column-correct PDB CA ATOM line (chain at col 22, B-factor at 61-66)."""
    return (
        f"ATOM  {serial:>5}  CA  ALA {chain}{resseq:>4}    "
        f"{0.0:8.3f}{0.0:8.3f}{0.0:8.3f}  1.00{bfac:6.2f}           C"
    )


def fake_complex_pdb(chains=("A", "B"), n_res: int = 20) -> str:
    """A tiny multi-chain PDB standing in for an RCSB complex (e.g. 8PSO A/B)."""
    lines, serial = ["HEADER    FAKE COMPLEX"], 1
    for ch in chains:
        for r in range(1, n_res + 1):
            lines.append(_atom(serial, ch, r))
            serial += 1
        lines.append(f"TER   {serial:>5}      ALA {ch}{n_res:>4}")
        serial += 1
    lines.append("END")
    return "\n".join(lines) + "\n"


def fake_esmfold_pdb(plddt: float = 85.0, n_res: int = 50) -> str:
    """A single-chain PDB with pLDDT in the B-factor column, like ESMFold output."""
    lines = [_atom(i, "A", i, bfac=plddt) for i in range(1, n_res + 1)]
    return "\n".join(lines) + "\nEND\n"


@pytest.fixture
def mock_structure_net(monkeypatch):
    """Patch the structure network seams AND force pockets/docking to the mock path.

    Tool auto-discovery (vta.toolconfig) now finds the locally-built fpocket/vina, so
    without this a graph test would run real FPocket + real Vina (slow, networked,
    non-deterministic). Forcing the binaries to None keeps the suite hermetic; the
    real path is exercised separately by scripts/validate_controls.py.
    """
    import vta.nodes.admet as admet
    import vta.nodes.boltzina as boltzina
    import vta.nodes.docking as docking
    import vta.nodes.pockets as pockets
    import vta.nodes.structure as structure

    monkeypatch.setattr(structure, "fetch_rcsb_pdb", lambda pdb_id: fake_complex_pdb())
    monkeypatch.setattr(structure, "fold_esmfold", lambda seq: fake_esmfold_pdb())
    monkeypatch.setattr(structure, "_boltz2_bin", lambda: None)   # no Boltz-2 in CI
    monkeypatch.setattr(pockets, "_fpocket_bin", lambda: None)
    monkeypatch.setattr(pockets, "_p2rank_bin", lambda: None)         # no Java in CI
    monkeypatch.setattr(docking, "_vina_bin", lambda: None)
    monkeypatch.setattr(boltzina, "_boltzina_available", lambda: False)  # no GPU in CI
    monkeypatch.setattr(admet, "_admet_available", lambda: False)        # no torch in CI
    return structure
