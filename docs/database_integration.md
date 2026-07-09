# Database Integration Guide

This document explains how database integration is done in VTA-Agent and how to add
new sources without breaking the pipeline's reproducibility model.

## Integration Model

VTA-Agent uses two integration styles:

- Dedicated modules for sources used directly by the current graph.
- External adapters for sources that are available to future nodes, benchmarks, or
  manual review but are not part of the default lead-ranking path.

Dedicated modules live close to the code that consumes them. Examples:

- `vta.nodes.structure`: RCSB PDB, ESMFold, AlphaFold DB.
- `vta.data.uniprot`: UniProt accession and sequence lookup.
- `vta.data.ligands`: ChEMBL ligand cache and refresh.
- `vta.data.pubchem`: PubChem SMILES fallback.
- `vta.data.msa`: committed MSA lookup for conservation.

External adapters live in `vta/data/external_databases.py`. Each adapter records:

- `key`: stable registry key matching `vta/data/databases.py`.
- `name`: display name.
- `access`: access mode such as `public_json`, `public_web`, `public_download`, or
  `manual_download`.
- `base_url`: human or API landing URL.
- `description`: what the source contributes.
- `build_url`: deterministic lookup URL builder.
- `json_api`: whether `fetch_json()` can call the source directly.
- `restricted`: whether automated fetching should not run by default.

## Registry Contract

`vta/data/databases.py` is the machine-readable catalog. Each entry should include:

- Stable `key`.
- Human `name`.
- `Category`.
- One-line `provides` description.
- Landing `url`.
- `Status.INTEGRATED`.
- Optional `api`.
- `wired_in` module path.

`wired_in` must import successfully. Tests enforce this so the docs cannot claim an
integration that has no code.

## How To Add A New Source

1. Decide whether the source is used by the graph now.
2. If yes, create a dedicated helper module or node seam.
3. If no, add an `ExternalAdapter`.
4. Add the source to `DATABASES` with `wired_in` set to the owning module.
5. Add tests for URL construction, adapter coverage, and importability.
6. Update `docs/databases.md`.
7. Run the focused tests.

Example adapter pattern:

```python
"example": ExternalAdapter(
    "example", "Example DB", "public_json",
    "https://example.org/api/search",
    "What the database provides.",
    _param("https://example.org/api/search", "q"), True),
```

Example registry pattern:

```python
Database("example", "Example DB", Category.LIGAND,
         "What the database provides",
         "https://example.org", Status.INTEGRATED,
         api="https://example.org/api/search",
         wired_in="vta.data.external_databases")
```

## Network And Test Rules

Runtime network code must be isolated behind small functions. Unit tests should not
call external services. Test URL construction, adapter metadata, and mocked responses
instead.

Use cached or committed assets for reproducible graph behavior. The default pipeline
should not depend on a remote service unless the node has a labelled fallback.

## What Integration Does Not Mean

Integration does not mean every database is automatically used in the default graph.
Some sources are:

- Reference databases for manual review.
- Benchmark datasets to download before validation.
- Web tools without public batch APIs.
- Licensed or vendor-distributed datasets.

For those, VTA-Agent now has an explicit adapter and registry entry, but a future node
must still decide how the data affects scoring, reporting, or validation.

## Validation Commands

Focused checks:

```bash
python -m pytest tests/test_databases.py tests/test_external_databases.py -q
```

Repository environment used in this workspace:

```bash
python -m pytest tests/test_databases.py tests/test_external_databases.py -q
```
