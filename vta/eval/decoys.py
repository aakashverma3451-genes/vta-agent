"""Property-matched decoy loading contracts for enrichment benchmarks."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class DecoySet:
    name: str
    path: str
    source: str
    property_matched: bool
    caveat: str
    license: str = "user-provided"


def registered_decoy_sets(root: str = "vta/data/decoys_cache") -> dict[str, DecoySet]:
    return {
        "demo": DecoySet(
            name="demo",
            path=str(Path(root) / "decoys_demo.smi"),
            source="committed demo decoys",
            property_matched=False,
            caveat="Demonstration-class only; not property-matched.",
            license="repository-local",
        ),
        "deepcoy": DecoySet(
            name="deepcoy",
            path=str(Path(root) / "deepcoy.smi"),
            source="DeepCoy property-matched decoys",
            property_matched=True,
            caveat="Use only with target-specific generation metadata and leakage checks.",
        ),
        "dekois2": DecoySet(
            name="dekois2",
            path=str(Path(root) / "dekois2.smi"),
            source="DEKOIS 2.0",
            property_matched=True,
            caveat="Use target-specific benchmark sets; verify redistribution license.",
        ),
        "phase8_chembl": DecoySet(
            name="phase8_chembl",
            path=str(Path(root) / "phase8_tilv_pb1_chembl_matched.smi"),
            source="ChEMBL API physicochemical property matching",
            property_matched=True,
            caveat="Presumed decoys; not experimentally verified inactives; underpowered.",
            license="ChEMBL CC BY-SA 3.0",
        ),
        "hcv_ns5b_ni": DecoySet(
            name="hcv_ns5b_ni",
            path=str(Path(root) / "hcv_ns5b_ni_matched.smi"),
            source="ChEMBL API physicochemical property matching (Phase 9B-3)",
            property_matched=True,
            caveat="Presumed decoys for HCV NS5B NI active-site actives; not experimentally "
            "verified inactives. Check decoys-per-active in the matching-quality report "
            "before treating any benchmark built on this set as powered.",
            license="ChEMBL CC BY-SA 3.0",
        ),
    }


def load_decoy_smiles(name: str = "demo", root: str = "vta/data/decoys_cache") -> tuple[list[dict], DecoySet]:
    sets = registered_decoy_sets(root)
    if name not in sets:
        raise KeyError(f"unknown decoy set: {name}")
    meta = sets[name]
    path = Path(meta.path)
    if not path.exists():
        raise FileNotFoundError(f"decoy set '{name}' not found at {path}")
    rows = []
    for i, line in enumerate(path.read_text().splitlines(), 1):
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        rows.append({"name": parts[1] if len(parts) > 1 else f"{name}_{i}", "smiles": parts[0]})
    return rows, meta
