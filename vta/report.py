"""vta.report — render a final VTAState into a self-contained HTML report.

One file, inline CSS, no JS framework (the audit trail uses native <details>). Built
for a 30-second read: classification + confidence badge, how each protein was
structured, the ranked leads, and a collapsible audit trail.

`render_report(state) -> str` is pure (unit-testable, no I/O); `write_report` is the
thin wrapper that writes outputs/{run_id}_report.html — the same split discipline as
classify()/build_contract().
"""
from __future__ import annotations

import html
import os
from typing import Any

from vta.state import VTAState

_CSS = """
:root { --green:#1e8449; --amber:#b9770e; --red:#a93226; --ink:#1c2833; --mute:#5d6d7e;
        --line:#e5e8eb; --bg:#f7f9fa; --card:#fff; --ctrl:#eafaf1; }
* { box-sizing:border-box; }
body { font:15px/1.5 -apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
       color:var(--ink); background:var(--bg); margin:0; padding:32px; }
.wrap { max-width:960px; margin:0 auto; }
h1 { font-size:24px; margin:0 0 2px; }
h2 { font-size:15px; text-transform:uppercase; letter-spacing:.05em; color:var(--mute);
     margin:28px 0 12px; }
.sub { color:var(--mute); font-style:italic; margin:0 0 16px; }
.card { background:var(--card); border:1px solid var(--line); border-radius:10px;
        padding:18px 20px; margin-bottom:12px; }
.badge { display:inline-block; padding:4px 12px; border-radius:999px; color:#fff;
         font-weight:600; font-size:13px; }
.g{background:var(--green);} .a{background:var(--amber);} .r{background:var(--red);}
.grid { display:grid; grid-template-columns:repeat(auto-fill,minmax(220px,1fr)); gap:12px; }
.prot { font-weight:600; } .meta { color:var(--mute); font-size:13px; margin-top:4px; }
table { width:100%; border-collapse:collapse; background:var(--card);
        border:1px solid var(--line); border-radius:10px; overflow:hidden; }
th,td { padding:9px 12px; text-align:left; border-bottom:1px solid var(--line); font-size:14px; }
th { background:#eef2f4; color:var(--mute); font-weight:600; text-transform:uppercase;
     font-size:12px; letter-spacing:.03em; }
td.num { text-align:right; font-variant-numeric:tabular-nums; }
tr.ctrl { background:var(--ctrl); }
.tag { font-size:11px; font-weight:700; color:var(--green); }
.scorecell { display:flex; align-items:center; justify-content:flex-end; gap:8px; }
.bar { background:var(--line); border-radius:4px; height:8px; width:64px; flex:0 0 auto; }
.bar > span { display:block; height:100%; background:var(--green); border-radius:4px; }
.note { color:var(--mute); font-size:13px; margin:10px 2px 0; }
.note b { color:var(--ink); }
details { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:0 18px; }
summary { cursor:pointer; padding:14px 0; font-weight:600; }
pre { font:13px/1.55 ui-monospace,SFMono-Regular,Menlo,monospace; white-space:pre-wrap;
      color:#2e3b46; margin:0 0 16px; }
.foot { color:var(--mute); font-size:12px; margin-top:24px; border-top:1px solid var(--line);
        padding-top:12px; }
.empty { color:var(--mute); font-style:italic; }
"""


def _esc(v: Any) -> str:
    return html.escape("" if v is None else str(v))


def _badge(conf: float | None, route: str | None) -> str:
    if route == "defer" or conf is None or conf < 85:
        return f'<span class="badge r">DEFERRED — {_esc(conf)}%</span>'
    cls, label = ("g", "HIGH") if conf >= 95 else ("a", "FLAG")
    return f'<span class="badge {cls}">{label} CONFIDENCE — {_esc(conf)}%</span>'


def _header(state: VTAState) -> str:
    tr = state.get("taxon_result") or {}
    genus = _esc(tr.get("genus") or "—")
    species = _esc(tr.get("species") or "—")
    basis = _esc(tr.get("confidence_basis") or "")
    conf = state.get("classification_confidence")
    badge = _badge(conf, state.get("route"))
    return (
        f'<div class="card"><h1>{species}</h1>'
        f'<p class="sub">genus <b>{genus}</b> &nbsp;·&nbsp; strain '
        f'{_esc(tr.get("strain_id") or "—")}</p>{badge}'
        f'<p class="meta">{basis}</p></div>'
    )


def _structures(state: VTAState) -> str:
    structs = state.get("structures") or {}
    if not structs:
        return '<p class="empty">No structures (run deferred or no proteins).</p>'
    cards = []
    for name, s in structs.items():
        method = s.get("method") or "—"
        if method == "experimental":
            detail = f'experimental PDB {_esc(s.get("source"))} · trusted ground truth'
        elif method == "esmfold":
            detail = f'ESMFold · mean pLDDT {_esc(s.get("mean_plddt"))}'
        elif method == "refused_too_long":
            detail = 'refused — exceeds ESMFold API ceiling (Phase 2: local fold)'
        else:
            detail = _esc(method)
        cards.append(f'<div class="card"><span class="prot">{_esc(name)}</span>'
                     f'<div class="meta">{detail}</div></div>')
    return f'<div class="grid">{"".join(cards)}</div>'


def _leads(state: VTAState) -> str:
    leads = state.get("lead_candidates") or []
    if not leads:
        return '<p class="empty">No lead candidates — run was deferred to a human expert.</p>'
    scores = [r.get("score") or 0 for r in leads]
    smax = max(scores) or 1
    head = ("<tr><th>#</th><th>Ligand</th><th>ChEMBL</th><th>Protein</th>"
            "<th class='num'>ΔG</th><th class='num'>LE</th>"
            "<th class='num'>cons</th><th class='num'>score</th></tr>")
    rows = []
    for i, r in enumerate(leads, 1):
        ctrl = r.get("positive_control")
        tag = ' <span class="tag">★ control</span>' if ctrl else ""
        pct = round(100 * (r.get("score") or 0) / smax)
        bar = (f'<div class="scorecell">{_esc(r.get("score"))}'
               f'<div class="bar"><span style="width:{pct}%"></span></div></div>')
        rows.append(
            f'<tr class="{"ctrl" if ctrl else ""}"><td class="num">{i}</td>'
            f'<td>{_esc(r.get("ligand"))}{tag}</td>'
            f'<td>{_esc(r.get("ligand_id"))}</td>'
            f'<td>{_esc(r.get("protein"))}</td>'
            f'<td class="num">{_esc(r.get("dG"))}</td>'
            f'<td class="num">{_esc(r.get("le"))}</td>'
            f'<td class="num">{_esc(r.get("conservation"))}</td>'
            f'<td class="num">{bar}</td></tr>'
        )
    notes = [
        '<p class="note"><b>Ranking:</b> leads are scored by <b>ligand efficiency</b> '
        '(binding energy per heavy atom), not raw ΔG. A compact, efficient binder can '
        'therefore outrank a larger molecule with a stronger absolute ΔG — and large '
        'prodrugs (whose active metabolite is smaller) rank conservatively. This is the '
        'standard medicinal-chemistry correction for docking\'s size bias.</p>'
    ]
    if all((r.get("conservation") == 0.5) for r in leads):
        notes.append('<p class="note"><b>Conservation:</b> shown as a 0.5 placeholder — '
                     'MSA-based pocket conservation is planned (a 10% term in the score).</p>')
    return f"<table>{head}{''.join(rows)}</table>{''.join(notes)}"


def _audit(state: VTAState) -> str:
    lines = state.get("audit_trail") or []
    body = _esc("\n".join(lines)) or "(empty)"
    return f"<details><summary>Audit trail ({len(lines)} decisions)</summary><pre>{body}</pre></details>"


def _footer(state: VTAState) -> str:
    v = state.get("versions") or {}
    tr = state.get("taxon_result") or {}
    parts = [f"{_esc(k)}={_esc(val)}" for k, val in v.items()]
    parts.append(f"taxonagent={_esc(tr.get('taxonagent_version'))}")
    parts.append(f"kg={_esc(tr.get('kg_version'))}")
    return (f'<p class="foot">run {_esc(state.get("run_id"))} &nbsp;·&nbsp; '
            f'{" · ".join(parts)}</p>')


def render_report(state: VTAState) -> str:
    """Render a VTAState into a single self-contained HTML string (pure, no I/O)."""
    tr = state.get("taxon_result") or {}
    title = _esc(tr.get("species") or "VTA-Agent report")
    return (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>VTA-Agent — {title}</title><style>{_CSS}</style></head><body><div class='wrap'>"
        "<h2>Viral Target Assessment</h2>"
        f"{_header(state)}"
        "<h2>Structures</h2>"
        f"{_structures(state)}"
        "<h2>Lead candidates</h2>"
        f"{_leads(state)}"
        "<h2>Provenance</h2>"
        f"{_audit(state)}"
        f"{_footer(state)}"
        "</div></body></html>"
    )


def write_report(state: VTAState, out_dir: str = "outputs") -> str:
    """Render and write outputs/{run_id}_report.html; return the path."""
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"{state.get('run_id', 'run')}_report.html")
    with open(path, "w") as fh:
        fh.write(render_report(state))
    return path
