"""Conservation node tests (SPEC #1) — hermetic (MSA seam monkeypatched)."""
from __future__ import annotations

import vta.nodes.conservation as C
from vta.state import new_state


# ── PDB helpers ──────────────────────────────────────────────────────────────
_ONE_TO_THREE = {v: k for k, v in C._THREE_TO_ONE.items()}


def _atom(serial, atom, aa, chain, resseq, x, y, z, icode=" "):
    res = _ONE_TO_THREE[aa]
    return (f"ATOM  {serial:>5} {atom:<4} {res} {chain}{resseq:>4}{icode}   "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {atom[0]:>2}")


def _residue_lines(serial0, aa, chain, resseq, base, icode=" "):
    """Two atoms (N, CA) for one residue near `base` xyz."""
    x, y, z = base
    return [
        _atom(serial0, "N", aa, chain, resseq, x, y, z, icode),
        _atom(serial0 + 1, "CA", aa, chain, resseq, x + 0.5, y, z, icode),
    ]


# ── JSD (Capra & Singh) ──────────────────────────────────────────────────────
def test_jsd_identical_column_is_high():
    # A perfectly conserved column → far from background → ~1.0.
    [score] = C.jensen_shannon_conservation(["WWWWWW"])
    assert score > 0.9


def test_jsd_background_like_column_is_low():
    # A column drawn to match the background distribution → ~0.
    col = "".join(C._AA)                       # one of each AA ≈ flat, close-ish to bg
    [score] = C.jensen_shannon_conservation([col])
    assert score < 0.25


def test_jsd_gappy_column_is_downweighted():
    full = C.jensen_shannon_conservation(["WWWW"])[0]
    gappy = C.jensen_shannon_conservation(["WW--"])[0]   # same residue, half gaps
    assert gappy < full
    assert abs(gappy - full * 0.5) < 1e-6                # linear gap down-weighting


def test_jsd_all_gap_column_is_zero():
    assert C.jensen_shannon_conservation(["----"]) == [0.0]


# ── structure parsing + proximity ────────────────────────────────────────────
def test_parse_and_sequence_observed_order():
    lines = (_residue_lines(1, "A", "B", 10, (0, 0, 0))
             + _residue_lines(3, "C", "B", 11, (10, 0, 0))
             + _residue_lines(5, "D", "B", 12, (20, 0, 0)))
    res = C.parse_residues("\n".join(lines) + "\n")
    assert C.structure_sequence(res) == "ACD"
    assert [r["key"][1] for r in res] == ["10", "11", "12"]


def test_proximity_is_index_based_under_nonstandard_numbering():
    # NUMBERING HAZARD: resseq starts at 100 and is non-contiguous (100,101,105).
    # Proximity must still resolve by OBSERVED ORDER (indices 0,1,2), not residue number.
    lines = (_residue_lines(1, "A", "B", 100, (0, 0, 0))
             + _residue_lines(3, "C", "B", 101, (50, 0, 0))
             + _residue_lines(5, "D", "B", 105, (100, 0, 0)))   # 3rd observed = index 2
    res = C.parse_residues("\n".join(lines) + "\n")
    assert C.structure_sequence(res) == "ACD"                   # contiguous despite gaps
    near = C.pocket_residue_indices(res, center=[100, 0, 0], radius=6.0)
    assert near == [2]                                          # only the 3rd residue


def test_row0_column_map_skips_gaps():
    assert C._row0_column_map("A-CD--E") == [0, 2, 3, 6]


# ── node end-to-end ──────────────────────────────────────────────────────────
def _state_with_structure(tmp_path, pdb_text, center):
    pdb = tmp_path / "PB1.pdb"
    pdb.write_text(pdb_text)
    st = new_state("cons-test", "cons-test")
    st["structures"] = {"PB1": {"pdb_path": str(pdb)}}
    st["pockets"] = {"PB1": [{"id": 1, "center": center, "conservation": 0.5}]}
    return st


def test_node_overwrites_placeholder_with_real_score(monkeypatch, tmp_path):
    # Structure "ACDEF"; pocket sits on residues 0-1. Homolog MSA: col0 fully conserved
    # (A), col1 variable -> mean JSD should be a real value in (0,1), not 0.5.
    lines = (_residue_lines(1, "A", "B", 1, (0, 0, 0))
             + _residue_lines(3, "C", "B", 2, (1, 0, 0))
             + _residue_lines(5, "D", "B", 3, (40, 0, 0))
             + _residue_lines(7, "E", "B", 4, (60, 0, 0))
             + _residue_lines(9, "F", "B", 5, (80, 0, 0)))
    st = _state_with_structure(tmp_path, "\n".join(lines) + "\n", center=[0.5, 0, 0])

    msa = ["ACDEF", "ACDEF", "AGDEF", "AKDEF"]      # col0 = AAAA (conserved), col1 varies
    monkeypatch.setattr(C, "fetch_homolog_msa", lambda seq, taxon=None: msa)
    out = C.conservation_node(st)
    cons = out["pockets"]["PB1"][0]["conservation"]
    assert cons != 0.5 and 0.0 <= cons <= 1.0
    assert any("conservation[PB1] pocket 1" in l for l in out["audit_trail"])
    assert "conservation" in out["versions"]


def test_node_keeps_placeholder_when_no_msa(monkeypatch, tmp_path):
    lines = _residue_lines(1, "A", "B", 1, (0, 0, 0)) + _residue_lines(3, "C", "B", 2, (1, 0, 0))
    st = _state_with_structure(tmp_path, "\n".join(lines) + "\n", center=[0.5, 0, 0])
    monkeypatch.setattr(C, "fetch_homolog_msa", lambda seq, taxon=None: None)
    out = C.conservation_node(st)
    assert out["pockets"]["PB1"][0]["conservation"] == 0.5
    assert any("[skip] conservation[PB1]: no MSA" in l for l in out["audit_trail"])


def test_node_skips_when_msa_row0_mismatches(monkeypatch, tmp_path):
    lines = _residue_lines(1, "A", "B", 1, (0, 0, 0)) + _residue_lines(3, "C", "B", 2, (1, 0, 0))
    st = _state_with_structure(tmp_path, "\n".join(lines) + "\n", center=[0.5, 0, 0])
    monkeypatch.setattr(C, "fetch_homolog_msa", lambda seq, taxon=None: ["WXYZ", "WXYZ"])
    out = C.conservation_node(st)
    assert out["pockets"]["PB1"][0]["conservation"] == 0.5
    assert any("row0 != structure sequence" in l for l in out["audit_trail"])
