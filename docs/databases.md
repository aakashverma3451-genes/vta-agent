# Databases for Drug Design

Reference catalog of the public databases that matter for a structure-based virtual
screening pipeline like VTA-Agent (target structures → pocket detection → docking →
MD → ADMET → ranking).

A machine-readable mirror of this list lives in
[`vta/data/databases.py`](../vta/data/databases.py) — the single source of truth the
pipeline reads. The `wired` column below reflects which sources are **actually
integrated today** vs. catalogued for future work; keep the two in sync.

---

## 1. Protein / Target Structure
| Database | What it provides | Wired |
|----------|------------------|:-----:|
| **PDB (RCSB)** | Experimental 3D structures (X-ray, cryo-EM, NMR) | ✅ `structure.py` |
| **PDBe** | European PDB mirror with richer annotations | — |
| **AlphaFold DB** | AI-predicted structures (200M+ proteins) | ◐ planned |
| **ESM Atlas / ESMFold** | On-the-fly single-sequence folding | ✅ `structure.py` |
| **UniProt** | Sequences, function, variants, PTMs | ◐ planned |
| **SWISS-MODEL Repository** | Homology models | — |
| **OPM** | Orientations of Proteins in Membranes | — |
| **CATH / SCOP** | Fold / domain classification | — |

## 2. Small-Molecule / Ligand Libraries
| Database | What it provides | Wired |
|----------|------------------|:-----:|
| **ChEMBL** | Bioactive molecules + measured activity (IC50, Ki) | ✅ `data/ligands.py` |
| **PubChem** | 100M+ compounds, bioassays, properties | ◐ planned |
| **DrugBank** | Approved/experimental drugs + targets + pharmacology | — |
| **ZINC / ZINC20 / ZINC22** | Purchasable compounds for VS (billions) | — |
| **Enamine REAL** | Make-on-demand library (40B+) | — |
| **BindingDB** | Measured protein–ligand binding affinities | — |
| **COCONUT / NPASS** | Natural products | — |

## 3. Protein–Ligand Complexes / Binding Data
| Database | What it provides | Wired |
|----------|------------------|:-----:|
| **PDBbind** | Curated complexes + affinities (scoring/ML training) | — |
| **BindingMOAD** | High-quality protein–ligand binding data | — |
| **PDBe-KB** | Aggregated functional/binding-site annotations | — |
| **scPDB** | Druggable binding sites (pocket benchmark) | — |

## 4. Binding Site / Pocket
*(relevant to the P2Rank / fpocket consensus step)*
| Database | What it provides | Wired |
|----------|------------------|:-----:|
| **scPDB** | Annotated druggable pockets | — |
| **CASTp** | Computed pocket geometry/topology | — |
| **ProBiS / PDBeFold** | Binding-site comparison | — |

## 5. Target / Pathway / Disease
| Database | What it provides | Wired |
|----------|------------------|:-----:|
| **Open Targets** | Target–disease associations w/ evidence | — |
| **TTD** | Known/explored therapeutic targets | — |
| **KEGG** | Pathways, diseases, drug interactions | — |
| **Reactome** | Biological pathways | — |
| **DisGeNET** | Gene–disease associations | — |
| **GtoPdb (IUPHAR/BPS)** | Drug targets + ligand pharmacology | — |

## 6. ADMET / Toxicity / Pharmacokinetics
*(relevant to the ADMET annotation node)*
| Database | What it provides | Wired |
|----------|------------------|:-----:|
| **SwissADME** | ADME + drug-likeness predictions | ◐ planned |
| **admetSAR / ADMETlab** | ADMET + toxicity predictions | — |
| **Tox21 / ToxCast** | High-throughput toxicity screening data | — |
| **CTD** | Chemical–gene–disease interactions | — |
| **SIDER** | Drug side effects from labels | — |

## 7. Bioactivity / Assay / Screening Benchmarks
| Database | What it provides | Wired |
|----------|------------------|:-----:|
| **PubChem BioAssay** | HTS assay results | — |
| **DUD-E / DEKOIS / LIT-PCBA** | Decoy sets for docking/VS validation | — |

## 8. Repurposing / Approved Drugs
| Database | What it provides | Wired |
|----------|------------------|:-----:|
| **DrugBank** | Approved/experimental drug encyclopedia | — |
| **Drug Repurposing Hub (Broad)** | Repurposing candidates | — |
| **ClinicalTrials.gov** | Trial status of compounds | — |
| **FDA Orange Book** | Approved drug products | — |

---

## ⭐ Must-have core for this pipeline
1. **PDB + AlphaFold DB** → receptor structures
2. **UniProt** → target sequence/function
3. **ZINC / Enamine REAL** → screening library
4. **ChEMBL + BindingDB** → known actives & affinities (validation / ML)
5. **PDBbind** → scoring-function / GNINA-style training & benchmarking
6. **DrugBank** → repurposing & known drugs
7. **SwissADME / ADMETlab** → the ADMET annotation step
8. **DUD-E / LIT-PCBA** → benchmark docking accuracy

> Legend: ✅ integrated · ◐ planned / partial · — catalogued only
