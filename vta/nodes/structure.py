"""structure_node (REAL) — produce a 3D structure for each extracted protein.

Decision order per protein, in priority order:

  1. EXPERIMENTAL-FIRST. If the subunit has a solved cryo-EM chain, fetch that PDB
     and extract THAT CHAIN. A real 2.4 Å structure beats any prediction.
     Experimental structures get mean_plddt = None ("trusted ground truth").

  2. ESMFold (≤400 aa). For proteins within the public API ceiling, fold via
     api.esmatlas.com and compute mean pLDDT.

  3. Boltz-2 (>400 aa). For large orphan proteins that ESMFold refuses (e.g. PB2
     at 450 aa, TiLV segments 5–10), try the `boltz` CLI (Passaro et al. 2025,
     bioRxiv:2025.06.14.659707). Needs the `boltz` package + a GPU for production;
     degrades to step 4 if absent. Output CIF is converted to PDB via gemmi so the
     rest of the pipeline (FPocket, Vina) sees a standard PDB file.

  4. REFUSE (honestly). No experimental chain, no folding tool available → record
     pdb_path=None with a clear reason so the graph proceeds without crashing.

Chain-mapping note (8PSO COMPND): chain A = PA-like, chain B = putative PB1,
chain C = RdRp. No chain explicitly labelled PB2 → PB2 falls to ESMFold (refuses
at 450 aa) then Boltz-2 if installed.

Network I/O goes through `fetch_rcsb_pdb` / `fold_esmfold` / `fold_boltz2` so
tests can monkeypatch them and stay offline.
"""
from __future__ import annotations

import os
import subprocess
import tempfile

import requests

from vta.state import VTAState

# Public ESMFold API ceiling — hard limit confirmed by probe (413 above 400).
ESMFOLD_API_MAX_AA = 400

# Confident subunit → (PDB id, chain) from the 8PSO COMPND records. Only entries we
# can defend biologically are listed; PB2 is intentionally absent (no labelled chain).
EXPERIMENTAL_PDB = {
    "PA":  ("8PSO", "A"),   # COMPND: "POLYMERASE ACIDIC PROTEIN (PA-LIKE)" — 316 res
    "PB1": ("8PSO", "B"),   # COMPND: "PUTATIVE PB1" — 515 res (matches ~500 aa ORF)
}

_OUT_DIR = "structures"


# ── Boltz-2 helpers ──────────────────────────────────────────────────────────
def _boltz2_bin() -> str | None:
    from vta.toolconfig import find_tool
    return find_tool("boltz", "BOLTZ2_BIN")


def _cif_to_pdb(cif_path: str, pdb_path: str) -> bool:
    """Convert Boltz-2 CIF output → PDB using gemmi (already a dep via meeko)."""
    try:
        import gemmi
        st = gemmi.read_structure(cif_path)
        st.write_pdb(pdb_path)
        return os.path.exists(pdb_path)
    except Exception:
        return False


def fold_boltz2(sequence: str, name: str, out_dir: str) -> tuple[str, float] | None:
    """Fold a protein with Boltz-2 CLI; return (pdb_path, mean_plddt) or None.

    Boltz-2 input is a YAML file specifying the sequence. Output is a CIF in
    <out_dir>/predictions/<name>_model_0.cif which we convert to PDB via gemmi.
    Requires the `boltz` CLI (pip install boltz) and ideally a CUDA GPU; will
    run on CPU but is very slow (minutes per protein).
    """
    boltz = _boltz2_bin()
    if not boltz:
        return None
    with tempfile.TemporaryDirectory() as tmpdir:
        # Write YAML input — protein-only (no ligand for structure prediction).
        yaml_path = os.path.join(tmpdir, f"{name}.yaml")
        with open(yaml_path, "w") as fh:
            fh.write(f"sequences:\n  - protein:\n      id: A\n      sequence: {sequence}\n")
        result = subprocess.run(
            [boltz, "predict", yaml_path, "--out_dir", out_dir,
             "--cache", os.path.join(out_dir, ".boltz_cache")],
            capture_output=True, text=True, timeout=1800,  # 30 min ceiling
        )
        if result.returncode != 0:
            return None
        # Boltz-2 writes: <out_dir>/predictions/<yaml_stem>/<name>_model_0.cif
        import glob
        cifs = glob.glob(os.path.join(out_dir, "predictions", "**", "*.cif"),
                         recursive=True)
        if not cifs:
            return None
        cif_path = cifs[0]
        pdb_path = os.path.join(out_dir, f"{name}_boltz2.pdb")
        if not _cif_to_pdb(cif_path, pdb_path):
            return None
        # Boltz-2 writes confidence (pLDDT) in B-factor column, same as ESMFold.
        plddt = mean_plddt(open(pdb_path).read())
        return pdb_path, plddt


# ── injectable network seams (monkeypatched in tests) ────────────────────────
def fetch_rcsb_pdb(pdb_id: str) -> str:
    """Download a PDB file from RCSB. Raises on HTTP error."""
    r = requests.get(f"https://files.rcsb.org/download/{pdb_id}.pdb", timeout=60)
    r.raise_for_status()
    return r.text


def fold_esmfold(sequence: str) -> str:
    """Fold a single sequence via the public ESMFold API. Raises on HTTP error."""
    r = requests.post(
        "https://api.esmatlas.com/foldSequence/v1/pdb/",
        data=sequence,
        headers={"Content-Type": "text/plain"},
        timeout=300,
    )
    r.raise_for_status()
    return r.text


def fetch_alphafold(uniprot: str) -> str:
    """Fetch a predicted structure from AlphaFold DB by UniProt accession.

    AlphaFold DB has NO length ceiling (unlike the public ESMFold API's 400-aa cap)
    and returns a curated model with no local GPU needed, so it's the natural fill
    for large subunits like PB1/PB2 that ESMFold refuses. Like ESMFold/Boltz-2 it
    stores per-residue pLDDT in the B-factor column, so `mean_plddt` applies
    unchanged. Raises on HTTP error (404 = AlphaFold has no model for this
    accession). Catalogued in `vta.data.databases` as key "alphafold".
    """
    r = requests.get(
        f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb",
        timeout=60,
    )
    r.raise_for_status()
    return r.text


# ── pure helpers ─────────────────────────────────────────────────────────────
def _uniprot_accession(prot: dict) -> str | None:
    """The protein's UniProt accession under any of the keys upstream may use."""
    return prot.get("uniprot") or prot.get("uniprot_id") or prot.get("accession")
def extract_chain(pdb_text: str, chain: str) -> str:
    """Keep only the atom records for `chain` — the fix for whole-complex reuse."""
    keep = []
    for line in pdb_text.splitlines():
        rec = line[:6].strip()
        if rec in ("ATOM", "HETATM", "TER", "ANISOU") and len(line) > 21 and line[21] == chain:
            keep.append(line)
    keep.append("END")
    return "\n".join(keep) + "\n"


def mean_plddt(pdb_text: str) -> float:
    """ESMFold writes per-residue pLDDT into the B-factor column (CA atoms)."""
    vals = []
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA":
            try:
                vals.append(float(line[60:66]))
            except ValueError:
                pass
    return round(sum(vals) / len(vals), 1) if vals else 0.0


# ── node ─────────────────────────────────────────────────────────────────────
def structure_node(state: VTAState) -> VTAState:
    os.makedirs(_OUT_DIR, exist_ok=True)
    structures: dict[str, dict] = {}

    for name, prot in (state.get("extracted_proteins") or {}).items():
        seq = prot.get("sequence", "")
        n = len(seq)
        path = os.path.join(_OUT_DIR, f"{state['run_id']}_{name}.pdb")
        rec = {"pdb_path": None, "mean_plddt": None, "method": None, "source": None}

        if name in EXPERIMENTAL_PDB:
            pdb_id, chain = EXPERIMENTAL_PDB[name]
            try:
                chain_pdb = extract_chain(fetch_rcsb_pdb(pdb_id), chain)
                with open(path, "w") as fh:
                    fh.write(chain_pdb)
                # experimental = trusted ground truth → mean_plddt stays None (snag #3)
                rec.update(pdb_path=path, method="experimental",
                           source=f"{pdb_id}:{chain}")
                state["audit_trail"].append(
                    f"Structure[{name}]: EXPERIMENTAL {pdb_id} chain {chain} "
                    f"(pLDDT=None = trusted ground truth)")
            except Exception as e:  # network/HTTP failure — degrade, don't crash
                rec["method"] = "experimental_failed"
                state["audit_trail"].append(
                    f"Structure[{name}]: WARNING experimental fetch failed ({e}); "
                    f"no structure")

        elif n <= ESMFOLD_API_MAX_AA:
            try:
                pdb_text = fold_esmfold(seq)
                plddt = mean_plddt(pdb_text)
                with open(path, "w") as fh:
                    fh.write(pdb_text)
                flag = " (LOW CONFIDENCE)" if plddt < 70 else ""
                rec.update(pdb_path=path, mean_plddt=plddt, method="esmfold",
                           source="api.esmatlas.com")
                state["audit_trail"].append(
                    f"Structure[{name}]: ESMFold, mean pLDDT {plddt}{flag}")
            except Exception as e:
                rec["method"] = "esmfold_failed"
                state["audit_trail"].append(
                    f"Structure[{name}]: WARNING ESMFold failed ({e}); no structure")

        elif _boltz2_bin():
            # > ESMFold ceiling but Boltz-2 is available — try it.
            try:
                result = fold_boltz2(seq, name, _OUT_DIR)
                if result:
                    b2_path, plddt = result
                    flag = " (LOW CONFIDENCE)" if plddt < 70 else ""
                    rec.update(pdb_path=b2_path, mean_plddt=plddt,
                               method="boltz2", source="boltz predict (local)")
                    state["audit_trail"].append(
                        f"Structure[{name}]: Boltz-2, {n} aa, "
                        f"mean pLDDT {plddt}{flag}")
                else:
                    rec["method"] = "boltz2_failed"
                    state["audit_trail"].append(
                        f"Structure[{name}]: WARNING Boltz-2 produced no output; "
                        f"no structure")
            except Exception as e:
                rec["method"] = "boltz2_failed"
                state["audit_trail"].append(
                    f"Structure[{name}]: WARNING Boltz-2 failed ({e}); no structure")

        else:
            # > API ceiling, no Boltz-2 → refuse honestly.
            rec["method"] = "refused_too_long"
            state["audit_trail"].append(
                f"Structure[{name}]: REFUSED — {n} aa > ESMFold API ceiling "
                f"({ESMFOLD_API_MAX_AA}); install boltz for large-protein folding")

        structures[name] = rec

    state["structures"] = structures
    used = "esm2 (api.esmatlas.com, ≤400aa) + experimental PDB"
    if _boltz2_bin():
        used += " + Boltz-2 (>400aa)"
    state["versions"]["esmfold"] = used
    return state
