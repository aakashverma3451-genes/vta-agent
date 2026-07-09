"""Hermetic tests for the OpenAI-compatible adapter (Open WebUI integration).

Uses FastAPI's TestClient — no server, no network, no docking. Verifies the endpoints Open WebUI
depends on and that the honest triage verdict + non-removable disclaimer flow through.
"""
from __future__ import annotations

import json

import pytest
from fastapi.testclient import TestClient

from vta.service.openai_adapter import app
from vta.service.triage_chat import run_triage

client = TestClient(app)


def test_models_endpoint_advertises_the_model():
    r = client.get("/v1/models")
    assert r.status_code == 200
    ids = [m["id"] for m in r.json()["data"]]
    assert "vta-agent-triage" in ids


def test_health():
    assert client.get("/health").json()["status"] == "ok"


def _chat(text, stream=False):
    return client.post("/v1/chat/completions", json={
        "model": "vta-agent-triage", "stream": stream,
        "messages": [{"role": "user", "content": text}]})


def test_mpro_chat_returns_honest_downgrade_verdict():
    r = _chat("What about SARS-CoV-2 Mpro?")
    assert r.status_code == 200
    content = r.json()["choices"][0]["message"]["content"]
    assert "MPRO" in content
    assert "does NOT beat" in content                  # honest 2D-baseline verdict
    assert "hypotheses" in content.lower()             # non-removable disclaimer
    assert "Pinned benchmark" in content               # pinned frozen benchmark


def test_nucleotide_target_routes_to_annotate_only():
    r = _chat("triage HCV NS5B for me")
    content = r.json()["choices"][0]["message"]["content"]
    assert "annotate_only" in content


def test_raw_sequence_is_deferred_not_docked():
    seq = "M" + "ACDEFGHIKLMNPQRSTVWY" * 5
    content = _chat(seq).json()["choices"][0]["message"]["content"]
    assert "full folding" in content and "hypotheses" in content


def test_streaming_is_sse_and_terminates():
    r = _chat("Mpro", stream=True)
    assert r.status_code == 200
    body = r.text
    assert body.startswith("data: ") and body.rstrip().endswith("data: [DONE]")
    # every non-DONE data line is a valid chat.completion.chunk
    chunks = [ln[6:] for ln in body.splitlines() if ln.startswith("data: ") and "[DONE]" not in ln]
    for c in chunks:
        assert json.loads(c)["object"] == "chat.completion.chunk"


def test_api_key_enforced_when_set(monkeypatch):
    monkeypatch.setenv("VTA_API_KEY", "secret")
    assert client.get("/v1/models").status_code == 401
    assert client.get("/v1/models", headers={"Authorization": "Bearer secret"}).status_code == 200


def test_run_triage_is_offline_and_grounded():
    out = run_triage("main protease")
    assert out["target"] == "MPRO" and out["decision"] == "full_dock"
    assert "disclaimer" not in out["markdown"].lower() or "hypotheses" in out["markdown"]
