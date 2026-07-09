"""Phase 10: Ensemble docking on Mpro — quantify ranking change vs single-structure.

Compares a 3-conformer receptor ensemble (7L11:A holo primary + 6Y2E:A apo +
7K3T:A holo secondary) against the Phase 9 single-structure benchmark.

Design:
- Phase 9 primary scores (7L11:A) are loaded from the existing JSON (no re-docking).
- 6Y2E:A (apo Mpro, Zhang 2020, 1.75 Å) and 7K3T:A (holo, Moonshot inhibitor)
  are fetched from RCSB and prepped via OpenBabel (same fallback as Phase 9).
- Each of the 100 Phase 9 compounds is docked against the 2 new conformers using
  cached ligand PDBQTs from structures/.ligand_pdbqt/.
- Box center is His41/Cys145 CA midpoint computed per conformer — different crystal
  forms have different absolute coordinates; the 7L11 centroid cannot be reused.
- Ensemble score per compound = minimum dG across all 3 (best predicted binding).
- Enrichment computed with bootstrap CIs and compared to Phase 9 single-receptor.
- Ranking delta: for each top-10 active in Phase 9, record its ensemble rank vs
  single-structure rank (positive = worsened, negative = improved).
- Significance is a PAIRED bootstrap of Δ(ensemble−single) per compound (WI-3), not CI overlap.
- Honest about failures: any compound that fails all 3 conformers gets dG=None and
  is excluded from enrichment (labelled skip, same policy as Phase 9).
"""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

from vta.eval.metrics import bootstrap_enrichment_report, enrichment_report
from vta.nodes.docking import _LIG_CACHE, _run_vina
from vta.nodes.structure import extract_chain, fetch_rcsb_pdb
from vta.toolconfig import find_tool

# ── conformers ────────────────────────────────────────────────────────────────
# Primary: 7L11:A (Phase 9, holo, non-covalent Moonshot inhibitor XF1)
# Ensemble 1: 6Y2E:A (apo Mpro, Zhang et al. 2020, 1.75 Å — open/apo conformation)
# Ensemble 2: 7K3T:A (holo Mpro, Moonshot non-covalent inhibitor)
ENSEMBLE_CONFORMERS = [
    {"label": "6Y2E_apo",  "pdb_id": "6Y2E", "chain": "A",
     "note": "apo Mpro, Zhang 2020, 1.75 Å — open S4-subpocket loop"},
    {"label": "7K3T_holo", "pdb_id": "7K3T", "chain": "A",
     "note": "Mpro + non-covalent Moonshot inhibitor, holo alternative"},
]

# Phase 9 box center (XF1 centroid in 7L11:A) — fallback only; per-conformer centers
# are computed from His41/Cys145 CA midpoint to handle different crystal frames.
MPRO_CENTER = [-21.815, -4.216, -27.984]
VINA_EXHAUSTIVENESS = 8
VINA_SEED = 42
VINA_MODES = 5

PHASE9_JSON = Path("outputs/phase9/mpro_noncovalent_benchmark.json")
OUT_DIR = Path("outputs/phase10")
OUT_JSON = OUT_DIR / "ensemble_benchmark.json"
OUT_MD = OUT_DIR / "ensemble_benchmark.md"
STRUCTURES_DIR = Path("structures")


# ── network seam (monkeypatched in tests) ────────────────────────────────────
def _fetch_rcsb(pdb_id: str) -> str:
    return fetch_rcsb_pdb(pdb_id)


# ── box-center helper ─────────────────────────────────────────────────────────
def _dyad_center(pdb_path: str, chain: str = "A") -> list[float] | None:
    """His41 + Cys145 CA midpoint — per-conformer active-site center for Mpro.

    Each PDB crystal structure lives in its own coordinate frame, so the 7L11
    XF1 centroid cannot be reused for 6Y2E or 7K3T. His41 and Cys145 are
    present in every Mpro structure; their CA midpoint reliably places the Vina
    search box on the catalytic dyad regardless of crystal form.
    """
    his41 = cys145 = None
    for ln in open(pdb_path):
        if not ln.startswith("ATOM"):
            continue
        if ln[21] != chain:
            continue
        if ln[12:16].strip() != "CA":
            continue
        resnum = int(ln[22:26].strip())
        xyz = [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
        if resnum == 41:
            his41 = xyz
        elif resnum == 145:
            cys145 = xyz
    if his41 and cys145:
        return [round((his41[i] + cys145[i]) / 2, 3) for i in range(3)]
    return None


# ── receptor prep ─────────────────────────────────────────────────────────────
def _prep_receptor_obabel(pdb_path: str, pdbqt: str) -> str | None:
    """OpenBabel receptor prep — same logic as docking._prep_receptor_obabel."""
    if os.path.exists(pdbqt):
        return pdbqt
    try:
        from openbabel import pybel
    except Exception:
        return None
    protein = [ln for ln in open(pdb_path) if ln.startswith(("ATOM", "TER"))]
    if not protein:
        return None
    tmp_pdb = f"{pdbqt}.receptor.pdb"
    with open(tmp_pdb, "w") as fh:
        fh.write("".join(protein))
        fh.write("END\n")
    try:
        mol = next(pybel.readfile("pdb", tmp_pdb))
    except Exception:
        return None
    try:
        mol.addh()
        mol.write("pdbqt", pdbqt, opt={"r": True}, overwrite=True)
    except Exception:
        return None
    if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) == 0:
        return None
    body = open(pdbqt).read()
    with open(pdbqt, "w") as fh:
        fh.write("REMARK VTA Phase 10 receptor prep: OpenBabel fallback\n")
        fh.write(body)
    return pdbqt


def fetch_and_prep_conformer(conf: dict) -> tuple[str | None, list[float] | None]:
    """Fetch PDB, extract chain, prep PDBQT; return (receptor_path, box_center).

    Box center is the His41/Cys145 CA midpoint in this structure's coordinate
    frame — essential because different crystal structures have different absolute
    coordinates (6Y2E active site differs from 7L11 by ~35 Å in raw PDB coords).
    """
    pdb_id, chain, label = conf["pdb_id"], conf["chain"], conf["label"]
    STRUCTURES_DIR.mkdir(exist_ok=True)
    pdb_path = STRUCTURES_DIR / f"phase10_{pdb_id}_{chain}.pdb"
    pdbqt_path = STRUCTURES_DIR / f"phase10_{pdb_id}_{chain}_rec.pdbqt"

    if not pdb_path.exists():
        try:
            raw = _fetch_rcsb(pdb_id)
            chain_pdb = extract_chain(raw, chain)
            if not chain_pdb.strip():
                print(f"  [warn] no ATOM records for {pdb_id} chain {chain}")
                return None, None
            pdb_path.write_text(chain_pdb)
        except Exception as e:
            print(f"  [warn] RCSB fetch failed for {pdb_id}: {e}")
            return None, None

    result = _prep_receptor_obabel(str(pdb_path), str(pdbqt_path))
    center = _dyad_center(str(pdb_path), chain)
    if not center:
        center = MPRO_CENTER
        print(f"  [warn] His41/Cys145 not found in {pdb_id}:{chain} — using MPRO_CENTER fallback")
    if result:
        print(f"  receptor ready: {pdbqt_path}  box_center: {center}")
    else:
        print(f"  [warn] receptor prep failed for {label} ({pdb_id}:{chain})")
        return None, None
    return result, center


# ── per-conformer docking ─────────────────────────────────────────────────────
def dock_compounds_against_receptor(
    compounds: list[dict],
    receptor_pdbqt: str,
    conformer_label: str,
    vina_bin: str,
    center: list[float] | None = None,
) -> dict[str, float | None]:
    """Dock all compounds against one receptor; return {compound_id: best_dG}."""
    box_center = center if center is not None else MPRO_CENTER
    scores: dict[str, float | None] = {}
    n_ok = n_fail = 0
    for c in compounds:
        cid = c["ligand_id"]
        lig_pdbqt = os.path.join(_LIG_CACHE, f"{cid}.pdbqt")
        if not os.path.exists(lig_pdbqt):
            scores[cid] = None
            n_fail += 1
            continue
        out_pose = str(STRUCTURES_DIR / f"phase10_{cid}_{conformer_label}.pdbqt")
        dG = _run_vina(vina_bin, receptor_pdbqt, lig_pdbqt, box_center, out_pose)
        scores[cid] = dG
        if dG is not None:
            n_ok += 1
        else:
            n_fail += 1
    print(f"  [{conformer_label}] docked {n_ok}/{n_ok+n_fail}, {n_fail} failed")
    return scores


# ── ensemble merge ────────────────────────────────────────────────────────────
def merge_ensemble(
    primary_rows: list[dict],
    new_scores: dict[str, dict[str, float | None]],
) -> list[dict]:
    """For each compound, take best dG across primary + all new conformers."""
    merged = []
    for row in primary_rows:
        cid = row["ligand_id"]
        primary_dG = row.get("dG")
        candidates: list[tuple[float, str]] = []
        if primary_dG is not None:
            candidates.append((primary_dG, "7L11_primary"))
        for label, score_dict in new_scores.items():
            dG_new = score_dict.get(cid)
            if dG_new is not None:
                candidates.append((dG_new, label))
        if not candidates:
            continue  # all conformers failed → labelled skip
        best_dG, best_label = min(candidates, key=lambda x: x[0])
        merged.append({
            **row,
            "ensemble_dG": best_dG,
            "ensemble_conformer": best_label,
            "primary_dG": primary_dG,
            "conformer_scores": {
                "7L11_primary": primary_dG,
                **{lbl: score_dict.get(cid) for lbl, score_dict in new_scores.items()},
            },
        })
    return merged


# ── enrichment helpers ────────────────────────────────────────────────────────
def _entries_from_rows(rows: list[dict], score_key: str = "ensemble_dG") -> list[dict]:
    out = []
    for r in rows:
        dG = r.get(score_key)
        if dG is None:
            continue
        out.append({
            "name": r["ligand_id"],
            "score": round(-float(dG), 4),
            "positive_control": bool(r.get("positive_control")),
            "dG": dG,
        })
    return out


def _ranking_delta(phase9_rows: list[dict], ensemble_rows: list[dict]) -> list[dict]:
    """For each top-10 Phase 9 active, report ensemble rank vs single-structure rank."""
    def rank_list(rows: list[dict], score_key: str) -> dict[str, int]:
        sorted_rows = sorted(
            [r for r in rows if r.get(score_key) is not None],
            key=lambda r: r[score_key]
        )
        return {r["ligand_id"]: i + 1 for i, r in enumerate(sorted_rows)}

    phase9_ranks = rank_list(phase9_rows, "dG")
    ensemble_ranks = rank_list(ensemble_rows, "ensemble_dG")

    top10_actives = sorted(
        [(cid, rank) for cid, rank in phase9_ranks.items()
         if any(r["ligand_id"] == cid and r.get("positive_control") for r in phase9_rows)],
        key=lambda x: x[1]
    )[:10]

    delta = []
    for cid, p9_rank in top10_actives:
        ens_rank = ensemble_ranks.get(cid)
        delta.append({
            "compound_id": cid,
            "phase9_rank": p9_rank,
            "ensemble_rank": ens_rank,
            "delta_rank": (ens_rank - p9_rank) if ens_rank is not None else None,
            "improved": (ens_rank is not None and ens_rank < p9_rank),
        })
    return delta


# ── main ──────────────────────────────────────────────────────────────────────
def run(n_resamples: int = 2000) -> dict:
    phase9 = json.loads(PHASE9_JSON.read_text())
    phase9_rows = phase9.get("rows", [])
    if not phase9_rows:
        raise RuntimeError(f"No Phase 9 rows found in {PHASE9_JSON}")

    vina = find_tool("vina", "VINA_BIN")
    if not vina:
        raise RuntimeError("AutoDock Vina not found — set VINA_BIN or install vina")

    compounds = [{"ligand_id": r["ligand_id"], "positive_control": r.get("positive_control"),
                  "smiles": r.get("smiles")} for r in phase9_rows]

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    STRUCTURES_DIR.mkdir(exist_ok=True)

    new_scores: dict[str, dict[str, float | None]] = {}
    conformer_status: dict[str, str] = {}
    conformer_centers: dict[str, list[float]] = {}
    for conf in ENSEMBLE_CONFORMERS:
        label = conf["label"]
        print(f"\nPreparing conformer {label} ({conf['pdb_id']}:{conf['chain']})...")
        receptor, center = fetch_and_prep_conformer(conf)
        if not receptor:
            conformer_status[label] = "prep_failed"
            print(f"  [skip] {label} receptor prep failed — excluded from ensemble")
            continue
        conformer_status[label] = "ok"
        conformer_centers[label] = center
        print(f"  Docking 100 compounds against {label} (center={center})...")
        new_scores[label] = dock_compounds_against_receptor(
            compounds, receptor, label, vina, center=center
        )

    if not new_scores:
        return _write_blocked(phase9, conformer_status)

    ensemble_rows = merge_ensemble(phase9_rows, new_scores)
    n_excluded = len(phase9_rows) - len(ensemble_rows)

    phase9_bench = phase9.get("benchmark", {})
    phase9_boot = (phase9.get("bootstrap") or {}).get("bootstrap") or {}

    ens_entries = _entries_from_rows(ensemble_rows, "ensemble_dG")
    ens_report = enrichment_report(ens_entries)
    ens_bootstrap = bootstrap_enrichment_report(ens_entries, n_resamples=n_resamples, seed=42)

    ranking_delta = _ranking_delta(phase9_rows, ensemble_rows)

    # Phase 11 WI-3: significance by PAIRED bootstrap on the shared per-compound scores
    # (single vs ensemble ΔG), NOT by comparing marginal CI ranges. CIs are still reported.
    from vta.eval.significance import paired_bootstrap_delta
    _labels = [1 if r.get("positive_control") else 0 for r in ensemble_rows]
    _ens_s = [-float(r["ensemble_dG"]) for r in ensemble_rows]
    _sng_s = [-float(r.get("primary_dG", r["dG"])) for r in ensemble_rows]
    paired = {m: paired_bootstrap_delta(_ens_s, _sng_s, _labels, m, n=n_resamples, seed=42)
              for m in ("bedroc", "roc_auc")}

    p9_bedroc_ci = (phase9_boot.get("bedroc") or {}).get("ci95")
    ens_bedroc_ci = (ens_bootstrap.get("bootstrap", {}).get("bedroc") or {}).get("ci95")
    p9_roc_ci = (phase9_boot.get("roc_auc") or {}).get("ci95")
    ens_roc_ci = (ens_bootstrap.get("bootstrap", {}).get("roc_auc") or {}).get("ci95")

    payload = {
        "phase": "Phase 10 — Ensemble Docking (Mpro)",
        "phase9_structure": "7L11 chain A (primary, holo non-covalent)",
        "ensemble_conformers": [
            {"label": c["label"], "pdb_id": c["pdb_id"], "chain": c["chain"],
             "note": c["note"], "status": conformer_status.get(c["label"], "not_run"),
             "box_center": conformer_centers.get(c["label"])}
            for c in ENSEMBLE_CONFORMERS
        ],
        "n_compounds": len(phase9_rows),
        "n_excluded_all_conformers_failed": n_excluded,
        "n_ensemble_docked": len(ensemble_rows),
        "new_conformers_included": list(new_scores.keys()),
        "primary_center_7L11": MPRO_CENTER,
        "vina_seed": VINA_SEED,
        "phase9_benchmark": {
            "bedroc": phase9_bench.get("bedroc"),
            "roc_auc": phase9_bench.get("roc_auc"),
            "log_auc": phase9_bench.get("log_auc"),
            "EF1%": (phase9_bench.get("ef") or {}).get("EF1%"),
            "ci_bedroc": p9_bedroc_ci,
            "ci_roc_auc": p9_roc_ci,
        },
        "ensemble_benchmark": {
            "bedroc": ens_report.get("bedroc"),
            "roc_auc": ens_report.get("roc_auc"),
            "log_auc": ens_report.get("log_auc"),
            "EF1%": (ens_report.get("ef") or {}).get("EF1%"),
            "ci_bedroc": ens_bedroc_ci,
            "ci_roc_auc": ens_roc_ci,
        },
        "delta": {
            "bedroc": (round(ens_report["bedroc"] - phase9_bench["bedroc"], 4)
                       if ens_report.get("bedroc") and phase9_bench.get("bedroc") else None),
            "roc_auc": (round(ens_report["roc_auc"] - phase9_bench["roc_auc"], 4)
                        if ens_report.get("roc_auc") and phase9_bench.get("roc_auc") else None),
            "paired_bedroc": paired["bedroc"],      # paired-bootstrap Δ(ensemble−single)
            "paired_roc_auc": paired["roc_auc"],
        },
        "ranking_delta_top10_actives": ranking_delta,
        "ensemble_bootstrap": ens_bootstrap,
        "ensemble_rows": ensemble_rows,
        "interpretation": _interpret(phase9_bench, ens_report, ranking_delta, paired["bedroc"]),
    }
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    _write_md(payload)
    print(f"\nPhase 10 complete → {OUT_JSON}")
    return payload


def _interpret(p9: dict, ens: dict, delta: list, paired_bedroc: dict) -> str:
    bedroc_delta = (ens.get("bedroc", 0) or 0) - (p9.get("bedroc", 0) or 0)
    improved_leads = sum(1 for d in delta if d.get("improved"))
    # Phase 11 WI-3: verdict from the PAIRED bootstrap Δ (ensemble − single), not CI overlap.
    sig = paired_bedroc.get("significant")
    if sig is None:
        msg = "paired test unavailable."
    elif sig and paired_bedroc.get("favours_a"):
        msg = (f"paired BEDROC Δ={paired_bedroc['median_delta']} (95% CI {paired_bedroc['ci95']}) "
               "is significantly > 0 — the ensemble improves enrichment.")
    elif sig:
        msg = (f"paired BEDROC Δ={paired_bedroc['median_delta']} (95% CI {paired_bedroc['ci95']}) "
               "is significantly < 0 — the ensemble worsens enrichment.")
    else:
        msg = (f"paired BEDROC Δ={paired_bedroc['median_delta']} (95% CI {paired_bedroc['ci95']}) "
               "includes 0 — no significant ensemble effect; single-structure signal is stable "
               "under conformational sampling (paired test, not CI overlap).")
    return (
        f"Ensemble (3 conformers) BEDROC {ens.get('bedroc')} vs Phase 9 single-structure "
        f"{p9.get('bedroc')} (delta {round(bedroc_delta, 4):+.4f}). "
        f"{msg} "
        f"Top-10 Phase 9 actives: {improved_leads}/10 improved rank in ensemble. "
        "Gate unchanged: DO NOT PROMOTE without a leakage-controlled DL improvement whose "
        "paired-bootstrap Δ vs Vina has a 95% CI strictly > 0."
    )


def _write_blocked(phase9: dict, conformer_status: dict) -> dict:
    payload = {
        "phase": "Phase 10 — Ensemble Docking (BLOCKED)",
        "status": "blocked",
        "reason": "All ensemble conformer receptor preps failed — check OpenBabel installation.",
        "conformer_status": conformer_status,
        "phase9_benchmark": phase9.get("benchmark"),
    }
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_JSON.write_text(json.dumps(payload, indent=2))
    print(f"BLOCKED — {payload['reason']}")
    return payload


def _write_md(payload: dict) -> None:
    p9 = payload["phase9_benchmark"]
    ens = payload["ensemble_benchmark"]
    delta = payload["delta"]
    ranking = payload.get("ranking_delta_top10_actives", [])

    def _ci(ci: list | None) -> str:
        return f"[{ci[0]:.3f}, {ci[1]:.3f}]" if ci else "n/a"

    lines = [
        "# Phase 10: Ensemble Docking — Mpro",
        "",
        f"Ensemble: 7L11:A (primary holo) + {', '.join(c['label'] for c in payload['ensemble_conformers'])}.",
        f"N compounds: {payload['n_compounds']}. "
        f"Excluded (all conformers failed): {payload['n_excluded_all_conformers_failed']}.",
        "Box centers: per-conformer His41/Cys145 CA midpoint (crystal frames differ).",
        "",
        "## Enrichment comparison",
        "",
        "| Metric | Phase 9 (single) | Phase 10 (ensemble) | Delta | paired Δ 95% CI (sig?) |",
        "|--------|-----------------|---------------------|-------|-------------|",
        f"| BEDROC(α=20) | {p9.get('bedroc')} {_ci(p9.get('ci_bedroc'))} | "
        f"{ens.get('bedroc')} {_ci(ens.get('ci_bedroc'))} | "
        f"{delta.get('bedroc'):+.4f} | {(delta.get('paired_bedroc') or {}).get('ci95')} "
        f"({(delta.get('paired_bedroc') or {}).get('significant')}) |",
        f"| ROC-AUC | {p9.get('roc_auc')} {_ci(p9.get('ci_roc_auc'))} | "
        f"{ens.get('roc_auc')} {_ci(ens.get('ci_roc_auc'))} | "
        f"{delta.get('roc_auc'):+.4f} | {(delta.get('paired_roc_auc') or {}).get('ci95')} "
        f"({(delta.get('paired_roc_auc') or {}).get('significant')}) |",
        f"| logAUC | {p9.get('log_auc')} | {ens.get('log_auc')} | "
        f"{round((ens.get('log_auc') or 0) - (p9.get('log_auc') or 0), 4):+.4f} | — |",
        f"| EF1% | {p9.get('EF1%')} | {ens.get('EF1%')} | "
        f"{round((ens.get('EF1%') or 0) - (p9.get('EF1%') or 0), 4):+.4f} | — |",
        "",
        "## Ranking delta — top 10 Phase 9 actives",
        "",
        "| Compound | Phase 9 rank | Ensemble rank | Delta | Improved? |",
        "|----------|-------------|---------------|-------|-----------|",
    ]
    for d in ranking:
        ens_r = d.get("ensemble_rank") or "—"
        dlt = f"{d['delta_rank']:+d}" if d.get("delta_rank") is not None else "—"
        lines.append(
            f"| {d['compound_id']} | {d['phase9_rank']} | {ens_r} | {dlt} | {d.get('improved')} |"
        )
    lines += [
        "",
        "## Interpretation",
        "",
        payload.get("interpretation", ""),
        "",
        "Gate: **DO NOT PROMOTE** — no DL rescore wired; modest Mpro docking signal unchanged.",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--resamples", type=int, default=2000)
    args = parser.parse_args()
    out = run(args.resamples)
    if out.get("status") == "blocked":
        return
    p9 = out["phase9_benchmark"]
    ens = out["ensemble_benchmark"]
    d = out["delta"]
    print(f"\nPhase 9  BEDROC={p9['bedroc']} ROC-AUC={p9['roc_auc']}")
    print(f"Ensemble BEDROC={ens['bedroc']} ROC-AUC={ens['roc_auc']}")
    print(f"Delta    BEDROC={d['bedroc']:+.4f} ROC-AUC={d['roc_auc']:+.4f}")
    print(f"paired-Δ BEDROC CI={(d.get('paired_bedroc') or {}).get('ci95')} "
          f"sig={(d.get('paired_bedroc') or {}).get('significant')}")
    print(f"Top-10 actives improved: "
          f"{sum(1 for x in out['ranking_delta_top10_actives'] if x.get('improved'))}/10")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
