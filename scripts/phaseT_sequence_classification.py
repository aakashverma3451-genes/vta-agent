"""Phase T — per-record viral sequence classification benchmark.

Classifies every contig in `~/taxonagent/data/raw/ViralSequences/*.fna` individually and
scores the predicted family against the ground-truth family encoded in the filename.

Two methodological points this script exists to get right:

1. **Per-record, not per-file.** `classify_genome()` (and the `run_pipeline()` it wraps)
   concatenates every record in a FASTA into one pseudo-genome — correct for its documented
   contract (one FASTA = one genome), but wrong here: these files each hold many *unrelated*
   metagenomic contigs (`k141_*`, distinct species and hosts). Concatenating 72 Picornaviridae
   contigs and predicting once would give 20 predictions, not 262. We split and classify each.

2. **`diamond` is a hard requirement, not an optional degradation.** It is the only step that
   assigns taxonomy. Absent, every record returns `genus=''`, `species=''`, `confidence=0.0`
   and an accuracy plot would read 0% for a purely environmental reason. We fail loudly here
   rather than silently produce a zero.

Ground truth is the ICTV **family** in the filename; the prediction is DIAMOND's
`taxonomy_summary["best_hit_family"]` (best-bitscore hit against the VMR exemplar DB).

Outcomes are recorded in three classes, never collapsed:
  correct       — predicted family == true family
  wrong         — predicted a different family
  unclassified  — no DIAMOND hit at all (usually a partial contig with no recognisable ORF)

Accuracy is therefore reported twice, and both numbers belong in any figure:
  strict accuracy   = correct / all records          (unclassified counts against)
  conditional acc.  = correct / classified records   (coverage-adjusted)

Usage:
    PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python \
        scripts/phaseT_sequence_classification.py [--workers 5] [--limit N]
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import shutil
import sys
import tempfile
import time
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

SEQ_DIR = Path(os.environ.get(
    "VTA_VIRAL_SEQ_DIR", str(Path.home() / "taxonagent/data/raw/ViralSequences")))
OUT = Path("outputs/phaseT/sequence_classification.json")


# ── FASTA splitting ────────────────────────────────────────────────────────
def parse_fasta(path: Path) -> List[Tuple[str, str]]:
    """Return [(header, sequence)] for a (possibly multi-record) FASTA."""
    records: List[Tuple[str, str]] = []
    header: Optional[str] = None
    chunks: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            elif header is not None:
                chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(chunks)))
    return records


def true_family(fna: Path) -> str:
    """`Coronaviridae_self.fna` -> `Coronaviridae`."""
    return fna.stem.replace("_self", "")


# ── one record ─────────────────────────────────────────────────────────────
def classify_record(task: Dict[str, Any]) -> Dict[str, Any]:
    """Classify a single contig. Runs in a worker process."""
    from taxonagent.pipeline.sequence_pipeline import run_pipeline

    header, seq, truth = task["header"], task["seq"], task["true_family"]
    contig_id = header.split()[0]
    t0 = time.time()

    tmpdir = tempfile.mkdtemp(prefix="phaseT_")
    try:
        fa = Path(tmpdir) / "record.fna"
        fa.write_text(f">{header}\n{seq}\n")

        # run_pipeline is chatty; keep worker stdout clean.
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf), contextlib.redirect_stderr(buf):
                result = run_pipeline(fa)
        except Exception as e:  # a crash on one contig must not kill the sweep
            return {"contig_id": contig_id, "header": header, "true_family": truth,
                    "predicted_family": None, "status": "error", "error": repr(e),
                    "seq_length": len(seq), "elapsed_s": round(time.time() - t0, 2)}

        summary = result.taxonomy_summary or {}
        pred = summary.get("best_hit_family") or None

        if pred is None:
            status = "unclassified"
        elif pred == truth:
            status = "correct"
        else:
            status = "wrong"

        rdrp = (result.rdrp_hits or [{}])[0] if result.rdrp_hits else {}
        return {
            "contig_id": contig_id,
            "header": header,
            "true_family": truth,
            "predicted_family": pred,
            "status": status,
            "best_hit_species": summary.get("best_hit_species"),
            "best_hit_pident": summary.get("best_hit_pident"),
            "best_hit_evalue": summary.get("best_hit_evalue"),
            "family_votes": summary.get("family_votes"),
            "rdrp_profile": rdrp.get("profile_name"),
            "rdrp_evalue": rdrp.get("evalue"),
            "num_orfs": result.num_orfs,
            "seq_length": result.seq_length,
            "n_diamond_hits": len(result.diamond_hits or []),
            "pipeline_errors": list(result.errors or []),
            "elapsed_s": round(time.time() - t0, 2),
        }
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


# ── aggregate ──────────────────────────────────────────────────────────────
def summarize(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    per_family: Dict[str, Dict[str, Any]] = {}
    by_fam: Dict[str, List[Dict]] = defaultdict(list)
    for r in records:
        by_fam[r["true_family"]].append(r)

    for fam, rows in sorted(by_fam.items()):
        n = len(rows)
        counts = Counter(r["status"] for r in rows)
        correct = counts["correct"]
        classified = counts["correct"] + counts["wrong"]
        confusion = Counter(r["predicted_family"] for r in rows
                            if r["status"] == "wrong" and r["predicted_family"])
        per_family[fam] = {
            "n_records": n,
            "correct": correct,
            "wrong": counts["wrong"],
            "unclassified": counts["unclassified"],
            "error": counts["error"],
            "strict_accuracy": round(correct / n, 4) if n else None,
            "conditional_accuracy": round(correct / classified, 4) if classified else None,
            "n_classified": classified,
            "confused_with": dict(confusion.most_common()),
            "median_pident": _median([r["best_hit_pident"] for r in rows
                                      if r.get("best_hit_pident")]),
        }

    n = len(records)
    counts = Counter(r["status"] for r in records)
    classified = counts["correct"] + counts["wrong"]
    return {
        "overall": {
            "n_records": n,
            "n_families": len(by_fam),
            "correct": counts["correct"],
            "wrong": counts["wrong"],
            "unclassified": counts["unclassified"],
            "error": counts["error"],
            "strict_accuracy": round(counts["correct"] / n, 4) if n else None,
            "conditional_accuracy": (round(counts["correct"] / classified, 4)
                                     if classified else None),
            "macro_strict_accuracy": (
                round(sum(v["strict_accuracy"] for v in per_family.values())
                      / len(per_family), 4) if per_family else None),
        },
        "per_family": per_family,
    }


def _median(xs: List[float]) -> Optional[float]:
    xs = sorted(x for x in xs if x is not None)
    if not xs:
        return None
    m = len(xs) // 2
    return round(xs[m] if len(xs) % 2 else (xs[m - 1] + xs[m]) / 2, 2)


# ── what VTA-Agent (not TaxonAgent) would DO with these calls ───────────────
def router_analysis(records: List[Dict[str, Any]]) -> Dict[str, Any]:
    """Apply the real confidence router to each contig and score what it acts on.

    TaxonAgent's `confidence_pct` IS `best_hit_pident` (an aa-identity band — see
    vta/nodes/router.py's honesty note), so the router's decision is fully determined by
    the identity we already recorded. Thresholds are imported, never restated, so this
    analysis cannot drift from the shipped router.

    The claim this supports is narrow and must stay narrow: the router does NOT *detect*
    misclassifications. It thresholds identity, and on this dataset every wrong call
    happens to be a low-identity call. That is a correlation observed over 11 errors, not
    an error-detection mechanism. The cost — correct calls also deferred — is reported
    beside it, because a router that defers everything would trivially "catch" all errors.
    """
    from vta.nodes.router import FLAG_THRESHOLD, PROCEED_THRESHOLD

    def decide(rec: Dict[str, Any]) -> str:
        # unclassified/error -> TaxonAgent emits confidence 0.0 -> defer
        conf = rec.get("best_hit_pident") if rec["status"] in ("correct", "wrong") else None
        conf = conf or 0.0
        if conf >= PROCEED_THRESHOLD:
            return "proceed"
        if conf >= FLAG_THRESHOLD:
            return "flag"
        return "defer"

    by_route: Dict[str, Counter] = defaultdict(Counter)
    for r in records:
        by_route[decide(r)][r["status"]] += 1

    n = len(records)
    acted_on = by_route["proceed"] + by_route["flag"]
    n_acted = sum(acted_on.values())
    n_wrong_total = sum(1 for r in records if r["status"] == "wrong")
    n_correct_total = sum(1 for r in records if r["status"] == "correct")

    return {
        "thresholds": {"proceed_gte": PROCEED_THRESHOLD, "flag_gte": FLAG_THRESHOLD,
                       "source": "vta.nodes.router (imported, not restated)"},
        "confidence_is": "TaxonAgent confidence_pct == best-hit aa identity; a similarity "
                         "heuristic, NOT a calibrated posterior",
        "counts": {rt: dict(c) for rt, c in sorted(by_route.items())},
        "n_acted_on": n_acted,
        "n_deferred": n - n_acted,
        "frac_acted_on": round(n_acted / n, 4) if n else None,
        "precision_on_acted": (round(acted_on["correct"] / n_acted, 4) if n_acted else None),
        "wrong_calls_deferred": f"{by_route['defer']['wrong']}/{n_wrong_total}",
        "correct_calls_deferred": f"{by_route['defer']['correct']}/{n_correct_total}",
        "interpretation": (
            "The router acts on a small, high-identity minority and defers the rest. On this "
            "set every misclassification falls below the defer threshold, so nothing wrong "
            "reaches docking — but so do most correct calls. Report both numbers together: "
            "precision on the acted-on set is meaningless without the deferral rate beside it."),
        "scope": (
            "These are divergent metagenomic contigs (median identity ~60%); on well-"
            "characterised reference genomes identities run high and far more would proceed. "
            "Do not generalise this deferral rate to other input distributions."),
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=5,
                    help="parallel workers (DIAMOND uses 4 threads each)")
    ap.add_argument("--limit", type=int, default=None, help="only first N records (smoke test)")
    ap.add_argument("--reanalyze", action="store_true",
                    help="recompute summary + router analysis from the existing JSON records "
                         "(no re-classification; the sweep takes ~18 min)")
    args = ap.parse_args()

    if args.reanalyze:
        if not OUT.exists():
            print(f"FATAL: {OUT} does not exist — run the sweep first.", file=sys.stderr)
            return 2
        payload = json.loads(OUT.read_text())
        records = payload["records"]
        payload["summary"] = summarize(records)
        payload["router_analysis"] = router_analysis(records)
        OUT.write_text(json.dumps(payload, indent=2))
        _print_report(payload)
        print(f"  -> {OUT} (reanalyzed, {len(records)} records; not re-classified)")
        return 0

    if not shutil.which("diamond"):
        print("FATAL: `diamond` is not on PATH. It is the ONLY step that assigns taxonomy;\n"
              "without it every record returns an empty genus/species at 0.0% confidence and\n"
              "this benchmark would report 0% accuracy for an environmental reason.\n"
              "Install it (`brew install diamond`) and re-run.", file=sys.stderr)
        return 2

    if not SEQ_DIR.is_dir():
        print(f"FATAL: sequence dir not found: {SEQ_DIR}", file=sys.stderr)
        return 2

    tasks: List[Dict[str, Any]] = []
    for fna in sorted(SEQ_DIR.glob("*.fna")):
        truth = true_family(fna)
        for header, seq in parse_fasta(fna):
            tasks.append({"header": header, "seq": seq, "true_family": truth})
    if args.limit:
        tasks = tasks[: args.limit]

    print(f"[phaseT] {len(tasks)} records across "
          f"{len({t['true_family'] for t in tasks})} families; {args.workers} workers")

    t0 = time.time()
    records: List[Dict[str, Any]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        futs = {pool.submit(classify_record, t): t for t in tasks}
        for i, fut in enumerate(as_completed(futs), 1):
            r = fut.result()
            records.append(r)
            mark = {"correct": "✓", "wrong": "✗", "unclassified": "·", "error": "!"}[r["status"]]
            print(f"  [{i:3d}/{len(tasks)}] {mark} {r['true_family']:<18} "
                  f"-> {str(r['predicted_family']):<18} {r['contig_id']}", flush=True)

    records.sort(key=lambda r: (r["true_family"], r["contig_id"]))
    elapsed = time.time() - t0
    summary = summarize(records)

    payload = {
        "artifact": "VTA-Agent Phase T — per-record viral family classification",
        "source_dir": str(SEQ_DIR),
        "ground_truth": "ICTV family encoded in the source filename (<Family>_self.fna)",
        "prediction": "DIAMOND best-bitscore hit vs the ICTV VMR exemplar DB "
                      "(taxonomy_summary.best_hit_family)",
        "unit": "one metagenomic contig per record (records classified INDEPENDENTLY; "
                "run_pipeline concatenates records when given a multi-record FASTA, which "
                "would conflate unrelated contigs)",
        "caveats": [
            "Contigs are metagenomic assemblies and may be partial; a contig with no "
            "recognisable ORF yields no DIAMOND hit and is recorded as `unclassified`, "
            "not as a wrong prediction.",
            "Class support is heavily imbalanced (see per_family.n_records): families with "
            "n=1 can only score 0% or 100%. Always plot n alongside accuracy.",
            "best_hit_pident is amino-acid identity to the nearest exemplar, a similarity "
            "heuristic, NOT a calibrated posterior probability.",
            "This measures TaxonAgent's family assignment only — the upstream classify step "
            "of the VTA pipeline. It says nothing about docking or the Phase-R gate.",
        ],
        "diamond_version": _diamond_version(),
        "elapsed_s": round(elapsed, 1),
        "summary": summary,
        "router_analysis": router_analysis(records),
        "records": records,
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2))

    print(f"\n[phaseT] done in {elapsed/60:.1f} min")
    _print_report(payload)
    print(f"  -> {OUT}")
    return 0


def _print_report(payload: Dict[str, Any]) -> None:
    o = payload["summary"]["overall"]
    print(f"  records      : {o['n_records']} across {o['n_families']} families")
    print(f"  correct      : {o['correct']}")
    print(f"  wrong        : {o['wrong']}")
    print(f"  unclassified : {o['unclassified']}   error: {o['error']}")
    print(f"  strict acc   : {o['strict_accuracy']}")
    print(f"  conditional  : {o['conditional_accuracy']}  (of classified only)")
    print(f"  macro strict : {o['macro_strict_accuracy']}  (mean over families)")
    ra = payload.get("router_analysis")
    if ra:
        print(f"\n  -- VTA-Agent router (not TaxonAgent) --")
        for rt in ("proceed", "flag", "defer"):
            print(f"  {rt:<12} : {ra['counts'].get(rt, {})}")
        print(f"  acts on      : {ra['n_acted_on']}/{o['n_records']} "
              f"({ra['frac_acted_on']:.1%}); precision on those {ra['precision_on_acted']}")
        print(f"  wrong calls deferred   : {ra['wrong_calls_deferred']}")
        print(f"  correct calls deferred : {ra['correct_calls_deferred']}  (the cost)")


def _diamond_version() -> Optional[str]:
    import subprocess
    try:
        out = subprocess.run(["diamond", "--version"], capture_output=True, text=True, timeout=20)
        return out.stdout.strip().splitlines()[0] if out.stdout else None
    except Exception:
        return None


if __name__ == "__main__":
    raise SystemExit(main())
