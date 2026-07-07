"""Render a final VTAState into a self-contained HTML report."""
from __future__ import annotations

import html
import os
from typing import Any

from vta.report_envelope import build_envelope
from vta.report_extra import footer_provenance, scientific_status
from vta.state import VTAState
from vta.nodes.verification import build_verdict

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
td.risk { color:var(--red); font-weight:600; }
.stable{color:var(--green);font-weight:600;}
.moderate{color:var(--amber);font-weight:600;}
.unstable{color:var(--red);font-weight:600;}
details { background:var(--card); border:1px solid var(--line); border-radius:10px; padding:0 18px; }
summary { cursor:pointer; padding:14px 0; font-weight:600; }
pre { font:13px/1.55 ui-monospace,SFMono-Regular,Menlo,monospace; white-space:pre-wrap;
      color:#2e3b46; margin:0 0 16px; }
.foot { color:var(--mute); font-size:12px; margin-top:24px; border-top:1px solid var(--line);
        padding-top:12px; }
.empty { color:var(--mute); font-style:italic; }
.envelope { border:2px solid var(--amber); background:#fffdf5; }
.envelope h3 { margin:0 0 6px; font-size:15px; color:var(--ink); }
.envelope .disc { font-weight:600; color:var(--ink); margin:0 0 10px; }
.envelope .pin { color:var(--mute); font-size:12px; margin:0 0 10px;
                 font-variant-numeric:tabular-nums; }
.envelope .verdict-no { color:var(--red); font-weight:600; }
.envelope .verdict-yes { color:var(--green); font-weight:600; }
.envelope .ood { color:var(--red); font-weight:600; }
.envelope ul { margin:8px 0 0; padding-left:18px; }
.envelope li { color:var(--mute); font-size:13px; margin:3px 0; }
.vbanner { border-radius:10px; padding:14px 18px; margin-bottom:12px; font-weight:600; }
.vbanner.down { background:#fdecea; border:2px solid var(--red); color:var(--red); }
.vbanner.pass { background:var(--ctrl); border:1px solid var(--green); color:var(--green); }
.vbanner.neutral { background:#eef2f4; border:1px solid var(--line); color:var(--mute); }
.vbanner .sub { color:var(--ink); font-weight:400; font-style:normal; margin-top:6px;
                font-size:13px; }
.annot .cav { color:var(--red); font-size:13px; margin:4px 0 0; }
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
    has_admet = any(r.get("admet") for r in leads)
    has_chem = any(r.get("active_species") or r.get("chemistry_flags") for r in leads)
    admet_head = ("<th class='num'>hERG</th><th class='num'>oral</th>"
                  "<th class='num'>solub</th>") if has_admet else ""
    chem_head = "<th>Active species</th><th>Flags</th>" if has_chem else ""
    head = (f"<tr><th>#</th><th>Ligand</th><th>ChEMBL</th><th>Protein</th>"
            f"<th class='num'>ΔG</th><th class='num'>LE</th>"
            f"<th class='num'>cons</th><th class='num'>score</th>{chem_head}{admet_head}</tr>")
    rows = []
    for i, r in enumerate(leads, 1):
        ctrl = r.get("positive_control")
        tag = ' <span class="tag">★ control</span>' if ctrl else ""
        pct = round(100 * (r.get("score") or 0) / smax)
        bar = (f'<div class="scorecell">{_esc(r.get("score"))}'
               f'<div class="bar"><span style="width:{pct}%"></span></div></div>')
        admet_cells = ""
        chem_cells = ""
        if has_chem:
            species = r.get("active_species") or {}
            flags = r.get("chemistry_flags") or {}
            source = r.get("species_source")
            active_label = species.get("active_form")
            if source and source != "parent":
                active_label = f"{active_label} ({source})"
            flag_names = []
            if flags.get("pains"):
                flag_names.append("PAINS")
            if flags.get("brenk"):
                flag_names.append("Brenk")
            if flags.get("aggregator"):
                flag_names.append("aggregator")
            if flags.get("beyond_ro5"):
                flag_names.append("bRo5")
            chem_cells = (
                f'<td>{_esc(active_label)}</td>'
                f'<td>{_esc(", ".join(flag_names) or "—")}</td>'
            )
        if has_admet:
            a = r.get("admet") or {}
            herg = a.get("herg")
            herg_cls = "num risk" if (herg is not None and herg > 0.5) else "num"
            admet_cells = (f'<td class="{herg_cls}">{_esc(herg)}</td>'
                           f'<td class="num">{_esc(a.get("oral"))}</td>'
                           f'<td class="num">{_esc(a.get("solubility"))}</td>')
        rows.append(
            f'<tr class="{"ctrl" if ctrl else ""}"><td class="num">{i}</td>'
            f'<td>{_esc(r.get("ligand"))}{tag}</td>'
            f'<td>{_esc(r.get("ligand_id"))}</td>'
            f'<td>{_esc(r.get("protein"))}</td>'
            f'<td class="num">{_esc(r.get("dG"))}</td>'
            f'<td class="num">{_esc(r.get("le"))}</td>'
            f'<td class="num">{_esc(r.get("conservation"))}</td>'
            f'<td class="num">{bar}</td>{chem_cells}{admet_cells}</tr>'
        )
    notes = [
        '<p class="note"><b>Ranking:</b> leads are ranked by <b>AutoDock Vina affinity '
        '(ΔG)</b> — the primary and only ranking term. Ligand efficiency (LE) and '
        'conservation are shown per compound but carry <b>weight 0</b>: LE was demoted to a '
        'reported annotation by the Phase 11 WI-6 paired-test gate (the LE-led composite did '
        'not beat ΔG-only on the powered Mpro benchmark; Kenny 2019). See the honesty '
        'envelope above for the benchmark this ranking rests on.</p>'
    ]
    if all((r.get("conservation") == 0.5) for r in leads):
        notes.append('<p class="note"><b>Conservation:</b> shown as a 0.5 placeholder — '
                     'no homolog MSA was available, so the real JSD score fell back to '
                     'neutral for this run. Annotation only (weight 0 in the ranking).</p>')
    else:
        notes.append('<p class="note"><b>Conservation:</b> real per-pocket Jensen–Shannon '
                     'divergence vs background (Capra &amp; Singh 2007) over a homolog MSA; '
                     'higher = more evolutionarily conserved. Annotation only (weight 0 in '
                     'the ranking).</p>')
    if has_admet:
        notes.append('<p class="note"><b>ADMET</b> (ADMET-AI): hERG = cardiotoxicity '
                     'probability (lower better; red = >0.5 risk), oral = predicted oral '
                     'bioavailability, solub = aqueous solubility (log mol/L). Annotation '
                     'only — does not affect ranking.</p>')
    if has_chem:
        notes.append('<p class="note"><b>Chemistry flags</b>: active-species/prodrug '
                     'annotations and PAINS/Brenk/aggregator/Ro5-style flags are '
                     'annotation-only. Parent prodrugs are not mechanistic proof of RdRp '
                     'inhibition unless the active species is modelled.</p>')
    return f"<table>{head}{''.join(rows)}</table>{''.join(notes)}"


def _md_section(state: VTAState) -> str:
    leads = state.get("md_validated_leads") or []
    if not leads:
        return ""
    head = ("<tr><th>#</th><th>Ligand</th><th>MD Verdict</th>"
            "<th class='num'>RMSD (Å)</th><th class='num'>MM-GBSA</th>"
            "<th>Key contacts</th><th class='num'>MD score</th></tr>")
    rows = []
    for i, r in enumerate(leads, 1):
        det = r.get("md_details") or {}
        rmsd = (det.get("rmsd") or {}).get("mean_rmsd")
        mmgbsa = ((det.get("mmgbsa") or {}).get("mean_binding_energy"))
        contacts_raw = (det.get("contacts") or {}).get("top_contacts") or []
        contacts = ", ".join(k for k, _ in contacts_raw[:4]) or "—"
        badge = r.get("md_badge", "—")
        if "STABLE" in badge:
            badge_html = f'<span class="stable">✓ {_esc(badge)}</span>'
        elif "MODERATE" in badge:
            badge_html = f'<span class="moderate">~ {_esc(badge)}</span>'
        elif "UNSTABLE" in badge:
            badge_html = f'<span class="unstable">✗ {_esc(badge)}</span>'
        else:
            badge_html = _esc(badge)
        rows.append(
            f'<tr><td class="num">{i}</td><td>{_esc(r.get("ligand"))}</td>'
            f'<td>{badge_html}</td>'
            f'<td class="num">{_esc(rmsd) if rmsd is not None else "—"}</td>'
            f'<td class="num">{_esc(mmgbsa) if mmgbsa is not None else "pending"}</td>'
            f'<td>{_esc(contacts)}</td>'
            f'<td class="num">{_esc(r.get("md_score"))}</td></tr>'
        )
    note = ('<p class="note"><b>MD validation</b> (OpenMM, 100 ns NPT, AMBER ff14SB + '
            'OpenFF 2.0): RMSD is mean ligand RMSD after backbone alignment. '
            'STABLE &lt;2 Å, MODERATE 2–4 Å, UNSTABLE &gt;4 Å (Yamaotsu &amp; '
            'Hirono 2016). MM-GBSA is pending AmberTools integration.</p>')
    return (f"<h2>MD-validated leads</h2>"
            f"<table>{head}{''.join(rows)}</table>{note}")


def _fep_section(state: VTAState) -> str:
    leads = state.get("fep_validated_leads") or []
    if not leads:
        return ""
    head = ("<tr><th>#</th><th>Ligand</th><th>Status</th>"
            "<th class='num'>ΔG</th><th class='num'>Error</th><th>Method</th></tr>")
    results = state.get("fep_results") or {}
    rows = []
    for i, lead in enumerate(leads, 1):
        result = results.get(lead.get("ligand"), {})
        rows.append(
            f'<tr><td class="num">{i}</td><td>{_esc(lead.get("ligand"))}</td>'
            f'<td>{_esc(lead.get("fep_badge") or result.get("status"))}</td>'
            f'<td class="num">{_esc(lead.get("fep_delta_g"))}</td>'
            f'<td class="num">{_esc(lead.get("fep_error"))}</td>'
            f'<td>{_esc(result.get("method") or result.get("reason"))}</td></tr>'
        )
    note = ('<p class="note"><b>FEP / ABFE</b>: optional final validation over the '
            'top MD-vetted leads. Values are shown only when an external ABFE runner '
            'is available; skipped rows are explicit and do not affect ranking.</p>')
    return f"<h2>FEP validation</h2><table>{head}{''.join(rows)}</table>{note}"


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
            f'{" · ".join(parts)}</p>{footer_provenance(state, _esc)}')


def _envelope(state: VTAState) -> str:
    """The non-removable honesty envelope (D0.5), rendered as the first card of every report.

    Always emitted — deferred, empty, or full run — so no output shows a ranking without its
    disclaimer, pinned benchmark, trivial-baseline verdict, pose reliability, and caveats.
    """
    env = state.get("honesty_envelope") or build_envelope(state)
    pin = env.get("benchmark_pin") or {}
    if pin.get("available"):
        pin_line = (f'benchmark: <b>{_esc(pin.get("artifact"))}</b> · frozen '
                    f'{_esc(pin.get("frozen_at"))} · hash {_esc(pin.get("content_hash"))} · '
                    f'gate: {_esc(pin.get("gate_decision"))}')
    else:
        pin_line = ('benchmark: <span class="ood">no pinned validation basis — '
                    'ranking is exploratory</span>')

    basis_rows = []
    for b in env.get("ranking_basis") or []:
        tb = b.get("trivial_baseline") or {}
        pr = b.get("pose_reliability") or {}
        beats = tb.get("beats_trivial_2d_baseline")
        if tb.get("available") and beats is not None:
            vcls = "verdict-yes" if beats else "verdict-no"
            vtxt = "beats 2D-sim" if beats else "does NOT beat 2D-sim"
            verdict = f'<span class="{vcls}">{vtxt}</span>'
        else:
            verdict = '<span class="mute">untested</span>'
        mci = b.get("metrics_ci") or {}
        roc = (mci.get("ROC_AUC") or {})
        roc_txt = (f'{_esc(roc.get("median"))} {_esc(roc.get("ci95"))}'
                   if roc else "—")
        basis_rows.append(
            f'<tr><td>{_esc(b.get("run_protein"))}</td>'
            f'<td>{_esc(b.get("benchmark_target"))}</td>'
            f'<td>{_esc(b.get("grade"))}</td>'
            f'<td class="num">{roc_txt}</td>'
            f'<td>{verdict}</td>'
            f'<td>{_esc(pr.get("status"))}</td></tr>'
        )
    basis_table = ""
    if basis_rows:
        basis_table = (
            '<table><tr><th>Run target</th><th>Benchmark</th><th>Grade</th>'
            '<th class="num">ROC-AUC (95% CI)</th><th>Docking vs trivial baseline</th>'
            f'<th>Pose</th></tr>{"".join(basis_rows)}</table>')

    ood = env.get("out_of_validated_domain") or []
    ood_html = ""
    if ood:
        ood_html = (f'<p class="ood">⚠ Out-of-validated-domain: {_esc(", ".join(ood))} — '
                    'no frozen benchmark covers this target; its ranking is exploratory and '
                    'is NOT backed by a validation grade.</p>')

    caveats = env.get("scoring_caveats") or []
    caveats_html = ""
    if caveats:
        items = "".join(f"<li>{_esc(c)}</li>" for c in caveats)
        caveats_html = f'<p class="note"><b>Scoring caveats:</b></p><ul>{items}</ul>'

    return (
        '<div class="card envelope"><h3>Honesty envelope — read before the leads</h3>'
        f'<p class="disc">{_esc(env.get("disclaimer"))}</p>'
        f'<p class="pin">{pin_line}</p>'
        f'{basis_table}{ood_html}'
        f'<p class="note">{_esc(env.get("ranking_note"))}</p>'
        f'{caveats_html}</div>'
    )


def _verification_banner(state: VTAState) -> str:
    """R4 verdict banner — structural, like the envelope: always rendered on a proceed run.

    A `downgrade` verdict tells the reader plainly that the docking ranking below did not clear
    the physics+statistics gate and is to be read as a ligand-based annotation, not enrichment.
    """
    if not (state.get("lead_candidates") or state.get("triage_decision")
            or state.get("annotation_rankings")):
        return ""   # deferred/empty run: nothing was ranked, nothing to gate
    verdict = state.get("verification_verdict") or build_verdict(state)
    v = verdict.get("verdict")
    if v == "downgrade":
        downgraded = sorted(p for p, t in (verdict.get("per_target") or {}).items()
                            if t.get("verdict") == "downgrade")
        return (
            '<div class="vbanner down">⚠ Verification: structure-based enrichment NOT '
            'demonstrated — docking ranking DOWNGRADED to a ligand-based annotation'
            f'{(" for " + ", ".join(downgraded)) if downgraded else ""}.'
            '<div class="sub">The ranking below did not beat a trivial 2D-similarity baseline '
            'by a paired-bootstrap test (or its pose reliability is unestablished). Treat it as '
            'an annotation, not validated enrichment.</div></div>')
    if v == "pass":
        return ('<div class="vbanner pass">✓ Verification: docking ranking cleared the hard '
                'gate (beats the 2D baseline by a paired test; pose-reliable).</div>')
    if v in {"defer", "refuse"}:
        return (f'<div class="vbanner neutral">Verification: run {v.upper()} — no docking '
                'enrichment claim is made.</div>')
    return ""


def _annotations(state: VTAState) -> str:
    """R5 labelled ligand-based annotations for annotate_only targets."""
    annots = state.get("annotation_rankings") or {}
    if not annots:
        return ""
    blocks = []
    for protein, a in annots.items():
        cavs = "".join(f'<p class="cav">⚠ {_esc(c)}</p>' for c in (a.get("caveats") or []))
        top = a.get("ranking") or []
        if top:
            rows = "".join(
                f'<tr><td class="num">{i}</td><td>{_esc(r.get("ligand"))}</td>'
                f'<td class="num">{_esc(r.get("similarity_to_actives"))}</td></tr>'
                for i, r in enumerate(top[:10], 1))
            body = (f'<p class="meta">{_esc(a.get("method"))}</p>'
                    f'<table><tr><th>#</th><th>Ligand</th>'
                    f'<th class="num">2D-sim to actives</th></tr>{rows}</table>')
        else:
            body = f'<p class="empty">{_esc(a.get("note"))}</p>'
        blocks.append(f'<div class="card annot"><span class="prot">{_esc(protein)}</span> '
                      f'<span class="meta">— {_esc(a.get("status"))}</span>{cavs}{body}</div>')
    return f'<h2>Ligand-based annotations (docking out-of-domain)</h2>{"".join(blocks)}'


def render_report(state: VTAState) -> str:
    """Render a VTAState into a single self-contained HTML string (pure, no I/O)."""
    tr = state.get("taxon_result") or {}
    title = _esc(tr.get("species") or "VTA-Agent report")
    return (
        "<!doctype html><html lang='en'><head><meta charset='utf-8'>"
        "<meta name='viewport' content='width=device-width,initial-scale=1'>"
        f"<title>VTA-Agent — {title}</title><style>{_CSS}</style></head><body><div class='wrap'>"
        "<h2>Viral Target Assessment</h2>"
        f"{_envelope(state)}"
        f"{_verification_banner(state)}"
        f"{_header(state)}"
        "<h2>Structures</h2>"
        f"{_structures(state)}"
        "<h2>Lead candidates</h2>"
        f"{_leads(state)}"
        f"{_annotations(state)}"
        f"{scientific_status(state, _esc)}"
        f"{_md_section(state)}"
        f"{_fep_section(state)}"
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
