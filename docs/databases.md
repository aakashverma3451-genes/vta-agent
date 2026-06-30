# Databases for Drug Design

Reference catalog of public databases relevant to VTA-Agent's structure-based
virtual-screening workflow:

```text
target structures -> pocket detection -> docking -> ranking -> ADMET -> validation
```

The machine-readable registry lives in `vta/data/databases.py`. Dedicated pipeline
sources have their own modules. The rest are integrated through
`vta/data/external_databases.py`, which records the access mode and builds stable
source-specific lookup URLs.

Legend: ✅ dedicated pipeline module · 🔌 external adapter

## 1. Protein / Target Structure

| Database | What it provides | Integration |
|----------|------------------|-------------|
| **PDB (RCSB)** | Experimental 3D structures from X-ray, cryo-EM, and NMR | ✅ `vta.nodes.structure` |
| **ESM Atlas / ESMFold** | On-the-fly single-sequence structure prediction | ✅ `vta.nodes.structure` |
| **AlphaFold DB** | AI-predicted structures for proteins with no ESMFold length cap | ✅ `vta.nodes.structure` |
| **UniProt** | Accession resolution for AlphaFold DB fallback | ✅ `vta.data.uniprot` |
| **Homolog MSA** | Viral homolog alignments for per-pocket JSD conservation | ✅ `vta.data.msa` |
| **PDBe** | European PDB mirror with richer annotations | 🔌 `vta.data.external_databases` |
| **SWISS-MODEL Repository** | Homology models of proteins | 🔌 `vta.data.external_databases` |
| **OPM** | Orientations of Proteins in Membranes | 🔌 `vta.data.external_databases` |

## 2. Small-Molecule / Ligand Libraries

| Database | What it provides | Integration |
|----------|------------------|-------------|
| **ChEMBL** | Bioactive molecules with measured activity such as IC50 and Ki | ✅ `vta.data.ligands` |
| **PubChem** | Name-to-SMILES fallback when ChEMBL misses a compound | ✅ `vta.data.pubchem` |
| **DrugBank** | Approved and experimental drugs with targets and pharmacology | 🔌 `vta.data.external_databases` |
| **ZINC22** | Purchasable compounds for virtual screening | 🔌 `vta.data.external_databases` |
| **Enamine REAL** | Make-on-demand compound library | 🔌 `vta.data.external_databases` |
| **BindingDB** | Measured protein-ligand binding affinities | 🔌 `vta.data.external_databases` |

## 3. Protein-Ligand Complex / Binding Data

| Database | What it provides | Integration |
|----------|------------------|-------------|
| **PDBbind** | Curated protein-ligand complexes with affinities for scoring and ML | 🔌 `vta.data.external_databases` |
| **BindingMOAD** | High-quality protein-ligand binding data | 🔌 `vta.data.external_databases` |
| **PDBe-KB** | Aggregated functional and binding-site annotations | 🔌 `vta.data.external_databases` |

## 4. Binding Site / Pocket

| Database | What it provides | Integration |
|----------|------------------|-------------|
| **scPDB** | Annotated druggable binding sites and pocket benchmarks | 🔌 `vta.data.external_databases` |
| **CASTp** | Computed pocket geometry and topology | 🔌 `vta.data.external_databases` |

## 5. Target / Pathway / Disease

| Database | What it provides | Integration |
|----------|------------------|-------------|
| **Open Targets** | Target-disease associations with evidence scoring | 🔌 `vta.data.external_databases` |
| **Therapeutic Target Database (TTD)** | Known and explored therapeutic targets | 🔌 `vta.data.external_databases` |
| **KEGG** | Pathways, diseases, and drug interactions | 🔌 `vta.data.external_databases` |
| **Guide to Pharmacology (IUPHAR/BPS)** | Drug targets and ligand pharmacology | 🔌 `vta.data.external_databases` |

## 6. ADMET / Toxicity / Pharmacokinetics

| Database | What it provides | Integration |
|----------|------------------|-------------|
| **SwissADME** | ADME and drug-likeness web predictions | 🔌 `vta.data.external_databases` |
| **ADMETlab / admetSAR** | ADMET and toxicity predictions | 🔌 `vta.data.external_databases` |
| **Tox21 / ToxCast** | High-throughput toxicity screening data | 🔌 `vta.data.external_databases` |

VTA-Agent's production ADMET annotation still uses the local `admet_ai` package when
installed. Public ADMET sites are integrated as reference adapters unless a future
node deliberately consumes them.

## 7. Bioactivity / Assay / Screening Benchmark

| Database | What it provides | Integration |
|----------|------------------|-------------|
| **DUD-E** | Decoy sets for docking and virtual-screening validation | 🔌 `vta.data.external_databases` |
| **LIT-PCBA** | Unbiased experimental benchmark for virtual screening and ML | 🔌 `vta.data.external_databases` |

## 8. Repurposing / Approved Drugs

| Database | What it provides | Integration |
|----------|------------------|-------------|
| **Drug Repurposing Hub (Broad)** | Curated repurposing candidates with annotations | 🔌 `vta.data.external_databases` |
| **ClinicalTrials.gov** | Clinical-trial status of compounds | 🔌 `vta.data.external_databases` |

## Integration Status

All databases in the registry now have code-level integration:

- Dedicated modules are used by the current VTA graph or data layer.
- External adapters expose source metadata, access mode, URL construction, and JSON
  fetching where a public JSON endpoint exists.
- Some external sources are web portals, licensed downloads, or manual benchmark
  downloads. Their adapters make the integration explicit, but the default graph does
  not automatically consume those datasets.

## Maintenance Rule

When adding or changing a source:

- Update `vta/data/databases.py`.
- Add or update an adapter in `vta/data/external_databases.py` unless the source has a
  dedicated pipeline module.
- Mark the registry `wired_in` field with the owning module.
- Update this document and `docs/database_integration.md`.
- Run `python -m pytest tests/test_databases.py tests/test_external_databases.py -q`.
