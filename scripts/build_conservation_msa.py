"""Build the committed TiLV-PB1 homolog MSA used by conservation scoring (SPEC #1).

Real homolog search + alignment, run OFFLINE and committed so the validation gate is
reproducible without network:

    1. target = 8PSO chain B (TiLV PB1), derived with the SAME parser the node uses,
       so the committed alignment's row-0 ungapped sequence matches at run time.
    2. blastp vs nr (viruses) → real viral-polymerase homologs.
    3. efetch full protein sequences for the top diverse hits.
    4. MAFFT-align (target first → row 0), write vta/data/msa/TiLV_PB1.afa.

Run from the worktree:
    PYTHONPATH=.:../taxonagent/src ../taxonagent/venv/bin/python \
        scripts/build_conservation_msa.py
"""
from __future__ import annotations

import io
import os
import subprocess
import sys
import time

from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML

from vta.nodes.structure import extract_chain, fetch_rcsb_pdb
from vta.nodes.conservation import parse_residues, structure_sequence

Entrez.email = os.environ.get("NCBI_EMAIL", "aakashverma3451@gmail.com")
_OUT = os.path.join(os.path.dirname(__file__), "..", "vta", "data", "msa", "TiLV_PB1.afa")
_MAX_HOMOLOGS = 30


def target_sequence() -> str:
    return structure_sequence(parse_residues(extract_chain(fetch_rcsb_pdb("8PSO"), "B")))


def blast_accessions(seq: str) -> list[str]:
    print("blastp vs nr (viruses) — this takes several minutes…", flush=True)
    handle = NCBIWWW.qblast(
        "blastp", "nr", seq, hitlist_size=80, expect=1e-3,
        entrez_query="viruses[organism]",
    )
    record = NCBIXML.read(handle)
    accs: list[str] = []
    for aln in record.alignments:
        acc = aln.accession
        if acc and acc not in accs:
            accs.append(acc)
    print(f"  {len(accs)} unique hit accessions", flush=True)
    return accs


def fetch_sequences(accs: list[str], target_len: int) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    for i in range(0, len(accs), 20):
        batch = accs[i:i + 20]
        h = Entrez.efetch(db="protein", id=",".join(batch), rettype="fasta", retmode="text")
        for rec in SeqIO.parse(io.StringIO(h.read()), "fasta"):
            s = str(rec.seq).replace("X", "")
            if 0.5 * target_len <= len(s) <= 1.6 * target_len:   # drop fragments/fusions
                out.append((rec.id, str(rec.seq)))
        time.sleep(0.4)
        if len(out) >= _MAX_HOMOLOGS:
            break
    return out[:_MAX_HOMOLOGS]


def main() -> None:
    target = target_sequence()
    print(f"target 8PSO:B = {len(target)} aa", flush=True)
    homologs = fetch_sequences(blast_accessions(target), len(target))
    print(f"kept {len(homologs)} homologs after length filter", flush=True)
    if len(homologs) < 3:
        sys.exit("too few homologs for a meaningful MSA — aborting (keep 0.5 fallback)")

    unaligned = f">target_8PSO_B\n{target}\n" + "".join(
        f">{acc}\n{seq}\n" for acc, seq in homologs)
    proc = subprocess.run(["mafft", "--auto", "--anysymbol", "-"],
                          input=unaligned, capture_output=True, text=True, check=True)

    os.makedirs(os.path.dirname(_OUT), exist_ok=True)
    with open(_OUT, "w") as fh:
        fh.write(proc.stdout)

    rows = [str(r.seq) for r in SeqIO.parse(io.StringIO(proc.stdout), "fasta")]
    row0 = rows[0].replace("-", "").upper()
    ok = row0 == target.upper()
    print(f"wrote {_OUT}: {len(rows)} rows x {len(rows[0])} cols; "
          f"row0==target: {ok}")
    if not ok:
        sys.exit("row0 ungapped != target — MSA reordered; fix before committing")


if __name__ == "__main__":
    main()
