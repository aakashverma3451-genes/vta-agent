"""Hermetic tests for the Mpro benchmark harness + active-site wiring (no real Vina)."""
from __future__ import annotations

import vta.nodes.structure as structure
from scripts.run_mpro_benchmark import control_draw_spread, select_library, to_ligand
from vta.nodes.pockets import EXPERIMENTAL_ACTIVE_SITE, pockets_node
from vta.nodes.structure import EXPERIMENTAL_PDB, structure_node
from vta.state import new_state


def _entry(name, score, active):
    return {"name": name, "score": score, "positive_control": active}


def test_control_draw_spread_deterministic_and_ordered():
    actives = [_entry(f"a{i}", 10 - i, True) for i in range(8)]
    inactives = [_entry(f"d{i}", -i, False) for i in range(20)]
    out1 = control_draw_spread(actives + inactives, draws=6, seed=7)
    out2 = control_draw_spread(actives + inactives, draws=6, seed=7)
    assert out1 == out2  # seeded determinism
    assert out1["draws"] == 6
    for metric, vals in out1["spread"].items():
        assert vals["min"] <= vals["median"] <= vals["max"]


def test_control_draw_spread_needs_both_classes():
    out = control_draw_spread([_entry("a", 1, True)], draws=3)
    assert out["draws"] == 0


def _fake_strat(n_moon_act=60, n_chembl_act=40, n_moon_inact=60):
    return {
        "non_covalent_actives": (
            [{"id": f"M{i}", "smiles": "CCO", "source": "moonshot"} for i in range(n_moon_act)]
            + [{"id": f"C{i}", "smiles": "CCN", "source": "chembl"} for i in range(n_chembl_act)]
        ),
        "inactives": [{"id": f"I{i}", "smiles": "CCC", "source": "moonshot"}
                      for i in range(n_moon_inact)],
    }


def test_select_library_prefers_moonshot_and_samples_requested_counts():
    a, i = select_library(_fake_strat(), n_actives=50, n_inactives=50, seed=17)
    assert len(a) == 50 and len(i) == 50
    assert all(r["source"] == "moonshot" for r in a)  # enough Moonshot -> no ChEMBL needed
    assert all(r["source"] == "moonshot" for r in i)


def test_select_library_tops_up_actives_from_chembl_when_short():
    strat = _fake_strat(n_moon_act=10, n_chembl_act=80, n_moon_inact=60)
    a, _ = select_library(strat, n_actives=50, n_inactives=50, seed=17)
    assert len(a) == 50
    assert any(r["source"] == "chembl" for r in a)  # topped up


def test_to_ligand_sanitizes_id_and_sets_label():
    lig = to_ligand({"id": "MAT-POS-/ 1", "smiles": "CCO", "source": "moonshot"}, True, 3)
    assert lig["positive_control"] is True
    assert "/" not in lig["chembl_id"] and " " not in lig["chembl_id"]


def test_mpro_is_wired_into_structure_and_active_site():
    assert EXPERIMENTAL_PDB.get("MPRO") == ("7L11", "A")
    site = EXPERIMENTAL_ACTIVE_SITE.get("MPRO")
    assert site and site["center"] == [-21.815, -4.216, -27.984]
    assert "7L11" in site["source"]


def test_structure_and_pockets_use_experimental_mpro(monkeypatch):
    # Tiny fake 7L11 chain-A PDB so structure_node stays offline.
    fake_pdb = (
        "ATOM      1  N   SER A   1      11.104  13.207  10.000  1.00 20.00           N\n"
        "ATOM      2  CA  SER A   1      12.560  13.207  10.000  1.00 20.00           C\n"
        "ATOM      3  C   HIS A  41      10.000   1.000   0.000  1.00 20.00           C\n"
        "END\n"
    )
    monkeypatch.setattr(structure, "fetch_rcsb_pdb", lambda pdb_id: fake_pdb)
    st = new_state("mpro_test", "mpro_test")
    st["extracted_proteins"] = {"MPRO": {"sequence": "M" * 306, "length_aa": 306, "plddt": None}}
    st = structure_node(st)
    assert st["structures"]["MPRO"]["pdb_path"]
    assert st["structures"]["MPRO"]["source"] == "7L11:A"
    st = pockets_node(st)
    pocket = st["pockets"]["MPRO"][0]
    assert pocket["center"] == [-21.815, -4.216, -27.984]
    assert pocket["method"] == "experimental_active_site"
