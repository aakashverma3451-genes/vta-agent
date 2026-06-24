"""conservation_contacts tests (SPEC #4) — hermetic (poses written to tmp_path)."""
from __future__ import annotations

import vta.nodes.conservation_contacts as CC
from vta.nodes.conservation import _THREE_TO_ONE
from vta.state import new_state

_ONE_TO_THREE = {v: k for k, v in _THREE_TO_ONE.items()}


def _atom(serial, atom, aa, chain, resseq, x, y, z, record="ATOM"):
    res = _ONE_TO_THREE[aa]
    pad = f"{record:<6}"
    return (f"{pad}{serial:>5} {atom:<4} {res} {chain}{resseq:>4}    "
            f"{x:>8.3f}{y:>8.3f}{z:>8.3f}  1.00  0.00          {atom[0]:>2}")


def _receptor_pdbqt(residues):
    """residues: [(aa, resseq, (x,y,z))] → a minimal receptor .pdbqt (one CA each)."""
    lines = [_atom(i + 1, "CA", aa, "B", rs, *xyz) for i, (aa, rs, xyz) in enumerate(residues)]
    return "\n".join(lines) + "\nEND\n"


def _pose_pdbqt(coords, second_model=None):
    """coords: [(x,y,z)] ligand atoms in pose/model 1 (HETATM); optional 2nd model after."""
    lines = ["MODEL 1"]
    lines += [_atom(i + 1, "C", "G", "L", 1, *c, record="HETATM") for i, c in enumerate(coords)]
    lines.append("ENDMDL")
    if second_model is not None:
        lines.append("MODEL 2")
        lines += [_atom(i + 1, "C", "G", "L", 1, *c, record="HETATM")
                  for i, c in enumerate(second_model)]
        lines.append("ENDMDL")
    return "\n".join(lines) + "\n"


def _write(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


# Three residues: 10 (high JSD) at origin, 11 (low JSD) at x=20, 12 (mid) at x=40.
_RESIDUES = [("W", 10, (0, 0, 0)), ("A", 11, (20, 0, 0)), ("D", 12, (40, 0, 0))]
_RES_CONS = {"PB1": {"10": 0.95, "11": 0.10, "12": 0.50}}


def _state_with_pose(tmp_path, ligand_xyz, pocket_cons=0.5, name="lig"):
    receptor = _write(tmp_path, "rec.pdbqt", _receptor_pdbqt(_RESIDUES))
    pose = _write(tmp_path, f"{name}.pdbqt", _pose_pdbqt(ligand_xyz))
    st = new_state("cc", "cc")
    st["residue_conservation"] = {k: dict(v) for k, v in _RES_CONS.items()}
    st["docking_results"] = [{
        "ligand": name, "protein": "PB1", "pocket": 1, "conservation": pocket_cons,
        "pose_path": pose, "receptor_path": receptor,
    }]
    return st


def test_ligand_near_conserved_residue_scores_high(tmp_path):
    # Ligand sits on residue 10 (JSD 0.95) → conservation ≈ 0.95, not the 0.5 placeholder.
    st = _state_with_pose(tmp_path, [(0, 0, 0)])
    out = CC.conservation_contacts_node(st)
    assert out["docking_results"][0]["conservation"] == 0.95
    assert any("contact-conservation: ligand-weighted 1" in l for l in out["audit_trail"])


def test_ligand_near_variable_residue_scores_low(tmp_path):
    # Ligand on residue 11 (JSD 0.10) → low conservation.
    st = _state_with_pose(tmp_path, [(20, 0, 0)])
    out = CC.conservation_contacts_node(st)
    assert out["docking_results"][0]["conservation"] == 0.10


def test_two_ligands_same_pocket_get_different_scores(tmp_path):
    # THE POINT of SPEC #4: two poses in one pocket touching different residues must
    # no longer share one pocket value.
    receptor = _write(tmp_path, "rec.pdbqt", _receptor_pdbqt(_RESIDUES))
    pose_hi = _write(tmp_path, "hi.pdbqt", _pose_pdbqt([(0, 0, 0)]))      # residue 10
    pose_lo = _write(tmp_path, "lo.pdbqt", _pose_pdbqt([(20, 0, 0)]))    # residue 11
    st = new_state("cc", "cc")
    st["residue_conservation"] = {k: dict(v) for k, v in _RES_CONS.items()}
    st["docking_results"] = [
        {"ligand": "Hi", "protein": "PB1", "pocket": 1, "conservation": 0.5,
         "pose_path": pose_hi, "receptor_path": receptor},
        {"ligand": "Lo", "protein": "PB1", "pocket": 1, "conservation": 0.5,
         "pose_path": pose_lo, "receptor_path": receptor},
    ]
    out = CC.conservation_contacts_node(st)
    hi, lo = out["docking_results"]
    assert hi["conservation"] != lo["conservation"]
    assert hi["conservation"] > lo["conservation"]


def test_no_pose_keeps_pocket_value(tmp_path):
    # Mock-docking record (no pose_path) → v1 pocket value untouched + labelled skip.
    st = new_state("cc", "cc")
    st["residue_conservation"] = {"PB1": {"10": 0.95}}
    st["docking_results"] = [{"ligand": "m", "protein": "PB1", "pocket": 1,
                              "conservation": 0.5}]
    out = CC.conservation_contacts_node(st)
    assert out["docking_results"][0]["conservation"] == 0.5
    assert any("[skip] contact-conservation: no pose" in l for l in out["audit_trail"])


def test_no_residue_map_keeps_pocket_value(tmp_path):
    # Real pose but no MSA-derived residue map → keep v1 value + labelled skip.
    st = _state_with_pose(tmp_path, [(0, 0, 0)])
    st["residue_conservation"] = {}
    out = CC.conservation_contacts_node(st)
    assert out["docking_results"][0]["conservation"] == 0.5
    assert any("no residue map" in l for l in out["audit_trail"])


def test_only_first_pose_model_is_used(tmp_path):
    # A second model far away must NOT pull in extra residues — parse model 1 only.
    receptor = _write(tmp_path, "rec.pdbqt", _receptor_pdbqt(_RESIDUES))
    pose = _write(tmp_path, "p.pdbqt", _pose_pdbqt([(0, 0, 0)], second_model=[(20, 0, 0)]))
    st = new_state("cc", "cc")
    st["residue_conservation"] = {k: dict(v) for k, v in _RES_CONS.items()}
    st["docking_results"] = [{"ligand": "x", "protein": "PB1", "pocket": 1,
                              "conservation": 0.5, "pose_path": pose,
                              "receptor_path": receptor}]
    out = CC.conservation_contacts_node(st)
    assert out["docking_results"][0]["conservation"] == 0.95   # only residue 10, not 11


def test_resseq_keys_match_between_map_and_receptor(tmp_path):
    # Key-match guard: resseq parsed from the receptor must hit the residue_conservation map.
    from vta.nodes.conservation import parse_residues
    residues = parse_residues(_receptor_pdbqt(_RESIDUES))
    parsed_keys = {r["key"][1] for r in residues}
    assert parsed_keys == set(_RES_CONS["PB1"])      # {"10","11","12"}


def test_noop_when_no_docking_results():
    st = new_state("cc", "cc")
    out = CC.conservation_contacts_node(st)
    assert out.get("docking_results") in (None, [])
