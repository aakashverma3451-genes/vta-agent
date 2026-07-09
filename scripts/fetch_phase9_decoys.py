"""Build property-matched, scaffold-distinct presumed decoys for the HCV NS5B NI
active-site actives (Phase 9B-3).

Mirrors the DUD-E / DEKOIS philosophy: candidate decoys are physicochemically
matched to each active (MW, logP, HBD, HBA, rotatable bonds, formal charge) but
scaffold-distinct, and are *presumed* inactives — not experimentally verified.

Network access is isolated behind an injectable ``fetch_fn`` seam so tests stay
hermetic (monkeypatch the seam) and a fetch failure degrades to a labelled,
provenance-bearing result rather than crashing.

Honesty contract:
  * Salts/counterions are stripped before any property is computed.
  * Active scaffolds are excluded from the decoy set.
  * If fewer than ``per_active`` decoys per active are recovered, the output is
    explicitly labelled demonstration-grade, not publication-grade.
  * No SMILES, decoy, or property is fabricated; missing data is reported.
"""
from __future__ import annotations

import argparse
import json
import statistics
import time
from pathlib import Path
from typing import Callable, Optional

import requests

from vta.eval.splits import scaffold_key

CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"
ACTIVES_PATH = Path("vta/data/phase9_actives/hcv_ns5b_ni_active_site.json")
OUT_SMI = Path("vta/data/decoys_cache/hcv_ns5b_ni_matched.smi")
OUT_PROVENANCE = Path("outputs/phase9/hcv_ns5b_ni_decoy_provenance.json")
OUT_QUALITY = Path("outputs/phase9/hcv_ns5b_ni_matching_quality.json")

# DUD-E-style physicochemical matching windows. MW is widened progressively in
# query_candidates; the others are fixed local filters.
TOL = {
    "mw": 50.0,
    "alogp": 1.5,
    "hba": 3,
    "hbd": 2,
    "rtb": 4,
    "charge": 1,
}
PUBLICATION_GRADE_PER_ACTIVE = 30

# Network seam type: takes a ChEMBL query-param dict, returns the molecules list.
FetchFn = Callable[[dict], list[dict]]


def largest_fragment_smiles(smiles: str | None) -> str | None:
    """Return the largest (most heavy atoms) dotted fragment — strips salts.

    Falls back to a plain string split when RDKit is unavailable so the seam
    stays usable in minimal environments.
    """
    if not smiles:
        return None
    if "." not in smiles:
        return smiles
    try:
        from rdkit import Chem

        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("unparseable")
        frags = Chem.GetMolFrags(mol, asMols=True, sanitizeFrags=False)
        if not frags:
            return smiles
        best = max(frags, key=lambda m: m.GetNumHeavyAtoms())
        return Chem.MolToSmiles(best)
    except Exception:
        parts = [p for p in smiles.split(".") if p]
        return max(parts, key=len) if parts else smiles


def compute_properties(smiles: str | None) -> Optional[dict]:
    """RDKit physicochemical descriptors of the salt-stripped parent, or None."""
    parent = largest_fragment_smiles(smiles)
    if not parent:
        return None
    try:
        from rdkit import Chem
        from rdkit.Chem import Crippen, Descriptors, Lipinski

        mol = Chem.MolFromSmiles(parent)
        if mol is None:
            return None
        return {
            "mw": round(Descriptors.MolWt(mol), 2),
            "alogp": round(Crippen.MolLogP(mol), 2),
            "hba": Lipinski.NumHAcceptors(mol),
            "hbd": Lipinski.NumHDonors(mol),
            "rtb": Lipinski.NumRotatableBonds(mol),
            "charge": Chem.GetFormalCharge(mol),
            "scaffold": scaffold_key(parent),
            "parent_smiles": parent,
        }
    except Exception:
        return None


def property_match(active_props: dict, cand_props: dict, tol: dict = TOL) -> bool:
    """True if a candidate falls inside every physicochemical window of an active."""
    return (
        abs(cand_props["mw"] - active_props["mw"]) <= tol["mw"]
        and abs(cand_props["alogp"] - active_props["alogp"]) <= tol["alogp"]
        and abs(cand_props["hba"] - active_props["hba"]) <= tol["hba"]
        and abs(cand_props["hbd"] - active_props["hbd"]) <= tol["hbd"]
        and abs(cand_props["rtb"] - active_props["rtb"]) <= tol["rtb"]
        and abs(cand_props["charge"] - active_props["charge"]) <= tol["charge"]
    )


def default_chembl_fetch(params: dict, timeout: int = 45, retries: int = 4) -> list[dict]:
    """Live ChEMBL molecule query with backoff.

    ChEMBL throttles rapid sequential queries (read timeouts / 429s), so transient
    failures are retried with exponential backoff. Only a genuinely exhausted retry
    budget raises requests.RequestException, which the caller degrades honestly.
    """
    last_exc: Optional[Exception] = None
    for attempt in range(retries):
        try:
            response = requests.get(CHEMBL, params={**params, "format": "json"}, timeout=timeout)
            response.raise_for_status()
            return response.json().get("molecules", [])
        except requests.RequestException as exc:
            last_exc = exc
            time.sleep(1.5 * (attempt + 1))
    raise last_exc  # type: ignore[misc]


def query_candidates(
    active_props: dict,
    per_active: int,
    pages: int,
    fetch_fn: FetchFn,
    widths: tuple[int, ...] = (50, 90, 140),
) -> list[dict]:
    """Query ChEMBL by progressively wider MW windows, keep local property matches."""
    page_size = min(100, max(40, per_active * 3))
    rows: list[dict] = []
    seen_ids: set[str] = set()
    for width in widths:
        for page in range(pages):
            params = {
                "limit": page_size,
                "offset": page * page_size,
                "molecule_properties__mw_freebase__gte": max(100, int(active_props["mw"] - width)),
                "molecule_properties__mw_freebase__lte": int(active_props["mw"] + width),
            }
            for mol in fetch_fn(params):
                structures = mol.get("molecule_structures") or {}
                smiles = structures.get("canonical_smiles")
                chembl_id = mol.get("molecule_chembl_id")
                if not smiles or not chembl_id or chembl_id in seen_ids:
                    continue
                cand_props = compute_properties(smiles)
                if cand_props and property_match(active_props, cand_props):
                    seen_ids.add(chembl_id)
                    rows.append({
                        "name": mol.get("pref_name") or chembl_id,
                        "chembl_id": chembl_id,
                        "smiles": smiles,
                        "props": cand_props,
                    })
            if len(rows) >= per_active:
                return rows
        if rows:
            # Got something at this width; don't widen further unless we have nothing.
            if len(rows) >= max(1, per_active // 2):
                return rows
    return rows


def load_actives(path: Path = ACTIVES_PATH) -> list[dict]:
    payload = json.loads(path.read_text())
    return payload.get("actives", [])


def build_decoys(
    actives: list[dict],
    per_active: int = 30,
    pages: int = 1,
    fetch_fn: FetchFn = default_chembl_fetch,
    delay: float = 0.34,
) -> dict:
    """Assemble property-matched, scaffold-distinct decoys; return provenance dict."""
    active_props = [
        {**a, "_props": compute_properties(a.get("canonical_smiles"))} for a in actives
    ]
    active_ids = {a.get("molecule_chembl_id") for a in actives}
    active_scaffolds = {
        ap["_props"]["scaffold"] for ap in active_props if ap["_props"]
    }

    decoys: list[dict] = []
    failures: list[dict] = []
    seen = set(active_ids)
    per_active_counts: dict[str, int] = {}

    for ap in active_props:
        active_id = ap.get("molecule_chembl_id")
        active_name = ap.get("pref_name") or active_id
        per_active_counts[active_id] = 0
        if not ap["_props"]:
            failures.append({
                "matched_active": active_id,
                "error": "active SMILES did not parse; cannot compute matching window",
            })
            continue
        try:
            candidates = query_candidates(ap["_props"], per_active, pages, fetch_fn)
        except requests.RequestException as exc:
            failures.append({
                "matched_active": active_id,
                "error": f"{type(exc).__name__}: {exc}",
            })
            continue
        kept = 0
        for cand in candidates:
            if cand["chembl_id"] in seen:
                continue
            if cand["props"]["scaffold"] in active_scaffolds:
                continue  # scaffold-distinct from every active
            seen.add(cand["chembl_id"])
            decoys.append({
                "name": cand["name"],
                "chembl_id": cand["chembl_id"],
                "smiles": cand["smiles"],
                "matched_active": active_id,
                "property_match_method": "ChEMBL MW window query; local MW/logP/HBA/HBD/RTB/charge match; scaffold-distinct",
                "props": cand["props"],
            })
            kept += 1
            if kept >= per_active:
                break
        per_active_counts[active_id] = kept
        if delay:
            time.sleep(delay)

    n_actives = len(actives)
    n_decoys = len(decoys)
    per_active_actual = (n_decoys / n_actives) if n_actives else 0.0
    publication_grade = (
        n_actives > 0
        and not failures
        and all(c >= PUBLICATION_GRADE_PER_ACTIVE for c in per_active_counts.values())
    )
    grade = "publication-grade-candidate" if publication_grade else "demonstration-grade"

    provenance = {
        "source": "ChEMBL molecule API",
        "method": "property-matched presumed decoys; MW window query + local MW/logP/HBA/HBD/RTB/charge match; scaffold-distinct; salts stripped to largest fragment",
        "target": "HCV NS5B NI active site",
        "actives_source": str(ACTIVES_PATH),
        "n_actives": n_actives,
        "per_active_requested": per_active,
        "pages_requested": pages,
        "n_decoys": n_decoys,
        "decoys_per_active_mean": round(per_active_actual, 3),
        "decoys_per_active_min": min(per_active_counts.values()) if per_active_counts else 0,
        "decoys_per_active_max": max(per_active_counts.values()) if per_active_counts else 0,
        "per_active_counts": per_active_counts,
        "grade": grade,
        "publication_grade": publication_grade,
        "n_failures": len(failures),
        "matching_tolerances": TOL,
        "output_smi": str(OUT_SMI),
        "caveat": (
            "Presumed decoys are not experimentally verified inactives. NI actives are "
            "nucleotide triphosphates / phosphoramidate prodrugs with unusual "
            "physicochemistry, so property-matched ChEMBL decoys are scarce; treat any "
            "low decoys-per-active count as an underpowered demonstration, not a claim."
        ),
        "decoys": decoys,
        "failures": failures,
    }
    return provenance


def summarize_matching_quality(provenance: dict, actives: list[dict]) -> dict:
    """Per-property distribution comparison of actives vs recovered decoys."""
    active_props = [p for p in (compute_properties(a.get("canonical_smiles")) for a in actives) if p]
    decoy_props = [d["props"] for d in provenance["decoys"] if d.get("props")]

    def _stats(values: list[float]) -> dict:
        if not values:
            return {"n": 0, "mean": None, "median": None, "min": None, "max": None}
        return {
            "n": len(values),
            "mean": round(statistics.fmean(values), 3),
            "median": round(statistics.median(values), 3),
            "min": round(min(values), 3),
            "max": round(max(values), 3),
        }

    properties = ["mw", "alogp", "hba", "hbd", "rtb", "charge"]
    distribution = {
        prop: {
            "actives": _stats([p[prop] for p in active_props]),
            "decoys": _stats([p[prop] for p in decoy_props]),
        }
        for prop in properties
    }

    return {
        "target": "HCV NS5B NI active site",
        "actives_source": str(ACTIVES_PATH),
        "n_actives": provenance["n_actives"],
        "n_decoys": provenance["n_decoys"],
        "decoys_per_active_mean": provenance["decoys_per_active_mean"],
        "decoys_per_active_min": provenance["decoys_per_active_min"],
        "decoys_per_active_max": provenance["decoys_per_active_max"],
        "publication_grade_threshold_per_active": PUBLICATION_GRADE_PER_ACTIVE,
        "grade": provenance["grade"],
        "publication_grade": provenance["publication_grade"],
        "scaffold_distinct_from_actives": True,
        "matching_tolerances": TOL,
        "property_distribution": distribution,
        "interpretation": (
            "Decoy property distributions should overlap the actives' windows while "
            "remaining scaffold-distinct. A small n_decoys or low decoys-per-active "
            "means the enrichment benchmark built on this set is underpowered and must "
            "be labelled a demonstration, consistent with the TiLV PB1 precedent."
        ),
    }


def write_outputs(provenance: dict, quality: dict) -> None:
    OUT_SMI.parent.mkdir(parents=True, exist_ok=True)
    OUT_PROVENANCE.parent.mkdir(parents=True, exist_ok=True)
    usable = [d for d in provenance["decoys"] if d.get("smiles") and d.get("chembl_id")]
    OUT_SMI.write_text(
        "".join(
            f"{d['smiles']}\t{d['name']}\t{d['chembl_id']}\tmatched_to={d['matched_active']}\n"
            for d in usable
        )
    )
    OUT_PROVENANCE.write_text(json.dumps(provenance, indent=2))
    OUT_QUALITY.write_text(json.dumps(quality, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--per-active", type=int, default=30)
    parser.add_argument("--pages", type=int, default=2)
    parser.add_argument("--delay", type=float, default=0.34)
    args = parser.parse_args()

    actives = load_actives()
    provenance = build_decoys(actives, args.per_active, pages=args.pages, delay=args.delay)
    quality = summarize_matching_quality(provenance, actives)
    write_outputs(provenance, quality)

    print(f"actives={provenance['n_actives']} decoys={provenance['n_decoys']} "
          f"per_active_mean={provenance['decoys_per_active_mean']} grade={provenance['grade']}")
    print(OUT_SMI)
    print(OUT_PROVENANCE)
    print(OUT_QUALITY)


if __name__ == "__main__":
    main()
