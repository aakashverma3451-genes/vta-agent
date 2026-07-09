# Cloudflare deployment

Two independent pieces — deploy either or both.

## 1. Figure set (Cloudflare Pages) — static, no code changes

The research-article figure set and pitch-deck diagrams are plain static files.

**Via dashboard:** Cloudflare dashboard → Workers & Pages → Create → Pages → Connect to Git →
select this repo/branch (`session-a`) →
- **Build command:** (leave empty)
- **Build output directory:** `docs/figures`

**Via CLI** (needs `npm i -g wrangler` and `wrangler login`):
```bash
wrangler pages deploy docs/figures --project-name=vta-agent-figures
```

Serves `docs/figures/index.html` as the landing page, linking to
`vta_figure_set.html` and the four `pitch/*.png` diagrams.

## 2. API (Cloudflare Worker, proxying to Vercel)

Cloudflare's compute runtime is JS/WASM-first. Its Python support is an experimental
Pyodide-based runtime that does not reliably run FastAPI + pydantic, so the actual triage
logic (`vta.service.openai_adapter`) stays on Vercel — already deployed and tested there.
`cloudflare/worker/` is a thin proxy Worker: it forwards `/health` and `/v1/*` to that
Vercel deployment, so the same API is reachable on a Cloudflare domain/edge (custom domain,
WAF, caching) without duplicating the reasoning layer in JS.

**Before deploying:** edit `cloudflare/worker/wrangler.toml` and set `UPSTREAM_URL` to your
real Vercel deployment URL (no trailing slash), e.g. `https://vta-agent.vercel.app`.

```bash
cd cloudflare/worker
wrangler login          # once
wrangler deploy
```

Verify:
```bash
curl -s https://vta-agent-api.<your-subdomain>.workers.dev/health
curl -s https://vta-agent-api.<your-subdomain>.workers.dev/v1/chat/completions \
  -H 'Content-Type: application/json' \
  -d '{"model":"vta-agent-triage","messages":[{"role":"user","content":"triage SARS-CoV-2 Mpro"}]}'
```

If you'd rather not depend on Vercel at all, the alternative is rewriting
`vta/service/openai_adapter.py` and `vta/service/triage_chat.py` as a native JS/TS Worker
(the reasoning logic itself — dossier/triage/verification — is simple enough to port; it
just isn't done here). Ask if you want that path instead of the proxy.
