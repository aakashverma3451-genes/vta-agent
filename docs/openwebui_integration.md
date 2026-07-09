# Open WebUI integration — chat with VTA-Agent

VTA-Agent ships a small **OpenAI-compatible adapter** so [Open WebUI](https://github.com/open-webui/open-webui)
(or any OpenAI client) can drive the agent as a chat "model". No fork of Open WebUI is needed —
it connects to the adapter over the standard OpenAI API.

```
Open WebUI  ──OpenAI /v1──▶  vta.service.openai_adapter  ──▶  dossier → triage → verification
   (chat UI)                    (FastAPI, this repo)              + honesty envelope (offline)
```

## What a chat turn does

You name a target; the agent runs the **real reasoning layer** (the same committed nodes the
pipeline uses) and returns an honest verdict grounded in the frozen benchmark:

- **routing decision** — `full_dock | annotate_only | defer | refuse` (with rationale)
- **benchmark grade + 95% CI** (BEDROC, ROC-AUC) from `outputs/phase9/locked_benchmark.json`
- **docking vs the trivial 2-D baseline** — the paired-test verdict (Mpro: *does NOT beat*)
- **pose reliability**, scoring caveats, and — for Mpro — the powered activity-cliff result
- the **non-removable disclaimer**: ranked, uncertainty-bearing hypotheses, not efficacy claims

It is **offline and fast** — no docking, no network folding — so an interactive turn can't hang.
A raw genome/FASTA is politely deferred (that needs the full folding + docking pipeline via
`build_app`, out of scope for a chat turn).

## 1. Start the adapter

```bash
# from the repo root, using the project venv
PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python -m vta.service
# serves http://0.0.0.0:8000  (endpoints: /v1/models, /v1/chat/completions, /health)
```

Environment: `VTA_HOST` (default `0.0.0.0`), `VTA_PORT` (default `8000`),
`VTA_API_KEY` (optional — if set, clients must send `Authorization: Bearer <key>`).

Verify:

```bash
curl -s http://localhost:8000/v1/models
curl -s http://localhost:8000/v1/chat/completions -H 'Content-Type: application/json' \
  -d '{"model":"vta-agent-triage","messages":[{"role":"user","content":"triage SARS-CoV-2 Mpro"}]}'
```

## 2. Run Open WebUI and point it at the adapter

Open WebUI is a separate app — run it however you deploy it (Docker is the usual path):

```bash
docker run -d -p 3000:8080 --add-host=host.docker.internal:host-gateway \
  -v open-webui:/app/backend/data --name open-webui \
  ghcr.io/open-webui/open-webui:main
```

Then in Open WebUI: **Settings → Connections → OpenAI API → add a connection**
- **API Base URL:** `http://host.docker.internal:8000/v1` (if Open WebUI is in Docker and the
  adapter runs on the host) — or `http://localhost:8000/v1` if both are on the host
- **API Key:** anything (or the `VTA_API_KEY` you set)

Pick the **`vta-agent-triage`** model in a new chat and ask, e.g. *"triage HCV NS5B"* or
*"does docking beat 2-D on Mpro?"*.

## Deployability — honest scope

- **The adapter is production-shaped and deployable now:** a small stateless FastAPI service,
  no GPU, no network, hermetically tested (`tests/test_openai_adapter.py`). Put it behind a
  reverse proxy + set `VTA_API_KEY` for anything exposed.
- **It serves the honest-triage verdict, not a full screening run.** Running real folding +
  docking interactively is out of scope; that is the deferred deployment work
  (`VTA-Agent_Deployment_Implementation_Plan.md`, D1–D10 — queue, workers, containers,
  governance) and would run as a background job, not a chat turn.
- Everything the chat says is grounded in committed benchmark artifacts and carries the
  non-removable honesty envelope — the same discipline as the rest of the pipeline.
