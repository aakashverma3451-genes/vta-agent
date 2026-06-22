"""CLI tests — the autonomous entry point (hermetic: mocked classify + tools off)."""
from __future__ import annotations

import pytest

pytest.importorskip("langgraph", reason="requires langgraph")

import taxonagent  # noqa: E402
from vta import cli  # noqa: E402


def _contract(conf: float = 98.0) -> dict:
    return {
        "strain_id": "TiLV", "genus": "Tilapinevirus", "species": "Tilapia lake virus",
        "confidence_pct": conf, "confidence_basis": "aa identity, not calibrated",
        "ictv_rule_match": [],
        "extracted_proteins": {"PB1": {"sequence": "M" * 500, "length_aa": 500, "plddt": None}},
        "reference_alignment": {}, "audit_trail": [], "timestamp": "t",
        "taxonagent_version": "0.1.0", "kg_version": "ICTV-VMR-MSL40",
    }


def _fasta(tmp_path):
    f = tmp_path / "genome.fasta"
    f.write_text(">x\nACGTACGT\n")
    return str(f)


def test_cli_run_end_to_end(monkeypatch, tmp_path, mock_structure_net):
    monkeypatch.chdir(tmp_path)                               # keep outputs/ in tmp
    monkeypatch.setattr(taxonagent, "classify_genome", lambda f: _contract())
    state = cli.run(_fasta(tmp_path))
    assert state["route"] == "proceed"
    assert state["lead_candidates"]                          # produced leads autonomously
    assert any(l.startswith("Report:") for l in state["audit_trail"])  # emitted a report


def test_cli_main_run_returns_zero(monkeypatch, tmp_path, mock_structure_net):
    monkeypatch.chdir(tmp_path)
    monkeypatch.setattr(taxonagent, "classify_genome", lambda f: _contract())
    assert cli.main(["run", _fasta(tmp_path)]) == 0


def test_cli_no_command_prints_help(capsys):
    assert cli.main([]) == 1
    assert "Autonomous viral target assessment" in capsys.readouterr().out
