/**
 * VTA-Agent API proxy — Cloudflare Worker.
 *
 * Cloudflare's compute runtime is JS/WASM-first; its Python support (Pyodide-based) does not
 * reliably run FastAPI + pydantic, so the real triage logic (vta.service.openai_adapter) stays
 * on Vercel, where it's already deployed and tested. This Worker just forwards the OpenAI-
 * compatible routes to that upstream, so the API is also reachable on a Cloudflare domain/edge
 * (custom domain, WAF, caching, etc.) without duplicating or reimplementing the reasoning layer.
 *
 * Set the upstream with `wrangler secret put UPSTREAM_URL` (or the `vars` in wrangler.toml for
 * a non-secret value), e.g. https://your-app.vercel.app — no trailing slash.
 */

const ALLOWED_PREFIXES = ["/health", "/v1/"];

function corsHeaders() {
  return {
    "Access-Control-Allow-Origin": "*",
    "Access-Control-Allow-Methods": "GET, POST, OPTIONS",
    "Access-Control-Allow-Headers": "Content-Type, Authorization",
  };
}

export default {
  async fetch(request, env) {
    const url = new URL(request.url);

    if (request.method === "OPTIONS") {
      return new Response(null, { status: 204, headers: corsHeaders() });
    }

    const allowed = ALLOWED_PREFIXES.some((p) => url.pathname.startsWith(p));
    if (!allowed) {
      return new Response(
        JSON.stringify({ error: "not found", allowed_paths: ALLOWED_PREFIXES }),
        { status: 404, headers: { "Content-Type": "application/json", ...corsHeaders() } },
      );
    }

    if (!env.UPSTREAM_URL) {
      return new Response(
        JSON.stringify({ error: "UPSTREAM_URL is not configured on this Worker" }),
        { status: 500, headers: { "Content-Type": "application/json", ...corsHeaders() } },
      );
    }

    const upstream = new URL(env.UPSTREAM_URL);
    upstream.pathname = url.pathname;
    upstream.search = url.search;

    const upstreamRequest = new Request(upstream.toString(), {
      method: request.method,
      headers: request.headers,
      body: ["GET", "HEAD"].includes(request.method) ? undefined : request.body,
    });

    const upstreamResponse = await fetch(upstreamRequest);
    const response = new Response(upstreamResponse.body, upstreamResponse);
    for (const [k, v] of Object.entries(corsHeaders())) response.headers.set(k, v);
    return response;
  },
};
