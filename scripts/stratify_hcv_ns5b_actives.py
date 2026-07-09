"""Stratify HCV NS5B actives into NI, NNI, off-target, or ambiguous classes."""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import requests

INFILE = Path("vta/data/phase9_actives/hcv_ns5b.json")
OUTFILE = Path("vta/data/phase9_actives/hcv_ns5b_stratified.json")
MECHANISM = "https://www.ebi.ac.uk/chembl/api/data/mechanism.json"

KNOWN_NI = {
    "SOFOSBUVIR",
    "MERICITABINE",
    "VALOPICITABINE",
    "R1479",
    "PSI-6130",
    "PSI-6206",
    "GS-461203",
}
KNOWN_NNI = {
    "FILIBUVIR",
    "DASABUVIR",
    "BECLABUVIR",
    "TEGOBUVIR",
    "DELEOBUVIR",
    "SETROBUVIR",
    "LOMIBUVIR",
    "NESBUVIR",
    "RADALBUIVIR",
}
KNOWN_OFF_TARGET = {
    "GLECAPREVIR": "NS3/NS4A protease inhibitor",
    "ELBASVIR": "NS5A inhibitor",
}
NNI_ASSAYS = {
    "CHEMBL5348065": "Filibuvir/benzofuran NS5B non-nucleoside inhibitor series",
}
INDIRECT_ASSAYS = {
    "CHEMBL5733646": "PBMC interferon-conditioned HCV replicon assay; indirect antiviral effect",
}
DIRECT_NS5B_ASSAYS = {
    "CHEMBL6128809": "Direct inhibition of Hepatitis C virus NS5B polymerase",
    "CHEMBL5391242": "Direct inhibition of Hepatitis C virus NS5B polymerase",
    "CHEMBL5635900": "Direct inhibition of Hepatitis C virus NS5B polymerase",
}
OFF_TARGET_TERMS = ("NS3/NS4A", "NS3", "NS5A", "SERINE PROTEASE", "PROTEASE")
NS5B_TERMS = ("RNA-DIRECTED RNA POLYMERASE", "RNA-DEPENDENT RNA POLYMERASE", "NS5B")


def _mechanisms(chembl_id: str) -> list[dict]:
    r = requests.get(MECHANISM, params={"molecule_chembl_id": chembl_id, "limit": 20, "format": "json"}, timeout=30)
    r.raise_for_status()
    return r.json().get("mechanisms", [])


def _looks_nucleoside(smiles: str | None) -> bool:
    if not smiles:
        return False
    s = smiles.lower()
    sugar = ("[c@h]1o" in s or "[c@@h]1o" in s) and "co" in s
    base = "n1cn" in s or "n1cc" in s or "n2ccc" in s
    phosphate = "p(=o)" in s
    return (sugar and base) or (phosphate and base)


def classify(active: dict, mechanisms: list[dict]) -> dict:
    name = (active.get("pref_name") or "").upper()
    smiles = active.get("canonical_smiles")
    mech_text = " | ".join((m.get("mechanism_of_action") or "") for m in mechanisms).upper()
    evidence = []

    if name in KNOWN_OFF_TARGET:
        evidence.append(f"known HCV drug mechanism is not NS5B: {KNOWN_OFF_TARGET[name]}")
        return {
            "mechanism_class": "off_target_hcv_drug",
            "binding_site": "not NS5B benchmark active",
            "benchmark_scope": "exclude",
            "confidence": "high",
            "evidence": evidence,
        }
    if name in KNOWN_NI or _looks_nucleoside(smiles):
        evidence.append("known NI or nucleoside/nucleotide-like scaffold")
        return {
            "mechanism_class": "NI",
            "binding_site": "NS5B catalytic active site",
            "benchmark_scope": "active_site",
            "confidence": "medium" if name not in KNOWN_NI else "high",
            "evidence": evidence,
        }
    if name in KNOWN_NNI:
        evidence.append("known NS5B non-nucleoside inhibitor")
        return {
            "mechanism_class": "NNI",
            "binding_site": "NS5B allosteric pocket",
            "benchmark_scope": "allosteric",
            "confidence": "high",
            "evidence": evidence,
        }
    assay_id = active.get("assay_chembl_id")
    if assay_id in NNI_ASSAYS:
        evidence.append(NNI_ASSAYS[assay_id])
        return {
            "mechanism_class": "NNI",
            "binding_site": "NS5B allosteric pocket",
            "benchmark_scope": "allosteric",
            "confidence": "medium",
            "evidence": evidence,
        }
    if assay_id in INDIRECT_ASSAYS:
        evidence.append(INDIRECT_ASSAYS[assay_id])
        return {
            "mechanism_class": "indirect_antiviral_assay",
            "binding_site": "not assignable to NS5B pocket",
            "benchmark_scope": "exclude",
            "confidence": "high",
            "evidence": evidence,
        }
    if any(term in mech_text for term in OFF_TARGET_TERMS):
        evidence.append(f"approved/curated mechanism is not NS5B: {mech_text}")
        return {
            "mechanism_class": "off_target_hcv_drug",
            "binding_site": "not NS5B benchmark active",
            "benchmark_scope": "exclude",
            "confidence": "high",
            "evidence": evidence,
        }
    if any(term in mech_text for term in NS5B_TERMS):
        evidence.append(f"ChEMBL mechanism supports NS5B polymerase: {mech_text}")
        return {
            "mechanism_class": "NNI_or_NI_ambiguous",
            "binding_site": "unknown NS5B pocket",
            "benchmark_scope": "ambiguous",
            "confidence": "low",
            "evidence": evidence,
        }
    if active.get("target_chembl_id") == "CHEMBL5375":
        if assay_id in DIRECT_NS5B_ASSAYS:
            evidence.append(DIRECT_NS5B_ASSAYS[assay_id])
        else:
            evidence.append("ChEMBL activity row maps to CHEMBL5375 but binding site is unknown")
        return {
            "mechanism_class": "NS5B_site_unknown",
            "binding_site": "unknown NS5B pocket",
            "benchmark_scope": "ambiguous",
            "confidence": "low",
            "evidence": evidence,
        }
    return {
        "mechanism_class": "ambiguous",
        "binding_site": "unknown",
        "benchmark_scope": "ambiguous",
        "confidence": "low",
        "evidence": ["no mechanism or binding-site evidence"],
    }


def stratify(path: Path = INFILE, fetch_mechanisms: bool = True) -> dict:
    data = json.loads(path.read_text())
    out = []
    for active in data["actives"]:
        mechanisms = _mechanisms(active["molecule_chembl_id"]) if fetch_mechanisms else []
        row = dict(active)
        row["chembl_mechanisms"] = mechanisms
        row["stratification"] = classify(active, mechanisms)
        out.append(row)
        if fetch_mechanisms:
            time.sleep(0.1)
    counts: dict[str, int] = {}
    for row in out:
        key = row["stratification"]["mechanism_class"]
        counts[key] = counts.get(key, 0) + 1
    payload = {
        "source": str(path),
        "target": data["target"],
        "n_actives": len(out),
        "counts": counts,
        "headline_scope": "NI active-site benchmark only when enough NI actives exist",
        "pocket_mapping": {
            "NI": "NS5B catalytic active site",
            "NNI": "NS5B allosteric thumb/palm pockets; separate benchmark required",
        },
        "actives": out,
    }
    OUTFILE.write_text(json.dumps(payload, indent=2))
    return payload


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--offline", action="store_true")
    args = parser.parse_args()
    out = stratify(fetch_mechanisms=not args.offline)
    print(json.dumps(out["counts"], indent=2))
    print(OUTFILE)


if __name__ == "__main__":
    main()
