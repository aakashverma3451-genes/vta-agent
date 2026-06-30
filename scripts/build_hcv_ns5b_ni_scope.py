"""Build the HCV NS5B NI active-site benchmark scope."""
from __future__ import annotations

import json
from pathlib import Path

from scripts.stratify_hcv_ns5b_actives import KNOWN_NI, _looks_nucleoside

INFILE = Path("vta/data/phase9_actives/hcv_rdrp_mechanism_target.json")
OUTFILE = Path("vta/data/phase9_actives/hcv_ns5b_ni_active_site.json")


def build_scope(path: Path = INFILE) -> dict:
    data = json.loads(path.read_text())
    actives = []
    for row in data["actives"]:
        name = (row.get("pref_name") or "").upper()
        if name in KNOWN_NI or _looks_nucleoside(row.get("canonical_smiles")):
            out = dict(row)
            out["mechanism_class"] = "NI"
            out["binding_site"] = "NS5B catalytic active site"
            out["benchmark_scope"] = "hcv_ns5b_active_site"
            out["stratification_evidence"] = (
                "nucleoside/nucleotide-like scaffold or known NI from ChEMBL "
                "RNA-directed RNA polymerase mechanism target"
            )
            actives.append(out)
    payload = {
        "label": "hcv_ns5b_ni_active_site",
        "source": str(path),
        "source_target": data["target"],
        "n_actives": len(actives),
        "binding_site": "NS5B catalytic active site",
        "mechanism_class": "NI",
        "pocket_mapping": {
            "NI": "dock to NS5B catalytic active site with Mg2+ and active triphosphate species where curated",
        },
        "caveat": (
            "Derived from ChEMBL RNA-directed RNA polymerase mechanism target, not only CHEMBL5375. "
            "Assay-level review is still required before publication claims."
        ),
        "actives": actives,
    }
    OUTFILE.write_text(json.dumps(payload, indent=2))
    return payload


def main() -> None:
    out = build_scope()
    print(f"NI active-site actives={out['n_actives']}")
    print(OUTFILE)


if __name__ == "__main__":
    main()
