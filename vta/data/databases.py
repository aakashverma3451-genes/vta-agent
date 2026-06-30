"""vta.data.databases — registry of the drug-design databases VTA-Agent draws on.

The human-readable version of this catalog is `docs/databases.md`; THIS module is
the machine-readable single source of truth the pipeline reads, mirroring the
self-describing discipline of `vta.toolconfig` (which catalogs external *tools*).

Each entry is honest about whether the source is actually wired into the pipeline:

    Status.INTEGRATED — code has a real module/adapter for this source
    Status.PLANNED    — a concrete integration seam is intended, not built yet
    Status.CATALOGUED — known-relevant, listed for completeness, no plan yet

Keeping this beside the code means a node can ask "where does my receptor / ligand /
ADMET data come from?" programmatically, and `tests/test_databases.py` can assert the
catalog stays consistent with `docs/databases.md`.
"""
from __future__ import annotations

import enum
from dataclasses import dataclass


class Category(str, enum.Enum):
    STRUCTURE = "Protein / Target Structure"
    LIGAND = "Small-Molecule / Ligand Library"
    COMPLEX = "Protein–Ligand Complex / Binding Data"
    POCKET = "Binding Site / Pocket"
    TARGET = "Target / Pathway / Disease"
    ADMET = "ADMET / Toxicity / Pharmacokinetics"
    BENCHMARK = "Bioactivity / Assay / Screening Benchmark"
    REPURPOSING = "Repurposing / Approved Drugs"


class Status(str, enum.Enum):
    INTEGRATED = "integrated"   # a node really pulls from it today
    PLANNED = "planned"         # a concrete seam is intended, not built yet
    CATALOGUED = "catalogued"   # relevant, listed for completeness


@dataclass(frozen=True)
class Database:
    """One drug-design data source."""

    key: str                    # short stable id, e.g. "rcsb_pdb"
    name: str                   # display name, e.g. "PDB (RCSB)"
    category: Category
    provides: str               # one line: what you get from it
    url: str                    # human landing page
    status: Status
    api: str | None = None      # programmatic endpoint, if one is used/intended
    wired_in: str | None = None # module that integrates it (when INTEGRATED)


# ── the catalog ──────────────────────────────────────────────────────────────
# Ordered by category to match docs/databases.md. `wired_in` names the module that
# owns the integration. Default graph sources have dedicated modules; future-facing
# sources use vta.data.external_databases adapters.
DATABASES: tuple[Database, ...] = (
    # 1. Structure ------------------------------------------------------------
    Database("rcsb_pdb", "PDB (RCSB)", Category.STRUCTURE,
             "Experimental 3D structures (X-ray, cryo-EM, NMR)",
             "https://www.rcsb.org", Status.INTEGRATED,
             api="https://files.rcsb.org/download/{pdb_id}.pdb",
             wired_in="vta.nodes.structure"),
    Database("esm_atlas", "ESM Atlas / ESMFold", Category.STRUCTURE,
             "On-the-fly single-sequence structure prediction",
             "https://esmatlas.com", Status.INTEGRATED,
             api="https://api.esmatlas.com/foldSequence/v1/pdb/",
             wired_in="vta.nodes.structure"),
    Database("alphafold", "AlphaFold DB", Category.STRUCTURE,
             "AI-predicted structures for 200M+ proteins (no length cap)",
             "https://alphafold.ebi.ac.uk", Status.INTEGRATED,
             api="https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_v4.pdb",
             wired_in="vta.nodes.structure"),
    Database("uniprot", "UniProt", Category.STRUCTURE,
             "Protein sequences + accession resolution (feeds AlphaFold DB)",
             "https://www.uniprot.org", Status.INTEGRATED,
             api="https://rest.uniprot.org/uniprotkb/search",
             wired_in="vta.data.uniprot"),
    Database("pdbe", "PDBe", Category.STRUCTURE,
             "European PDB mirror with richer annotations",
             "https://www.ebi.ac.uk/pdbe", Status.INTEGRATED,
             api="https://www.ebi.ac.uk/pdbe/api/",
             wired_in="vta.data.external_databases"),
    Database("swissmodel", "SWISS-MODEL Repository", Category.STRUCTURE,
             "Homology models of proteins",
             "https://swissmodel.expasy.org", Status.INTEGRATED,
             api="https://swissmodel.expasy.org/repository/uniprot/{accession}",
             wired_in="vta.data.external_databases"),
    Database("opm", "OPM", Category.STRUCTURE,
             "Orientations of Proteins in Membranes",
             "https://opm.phar.umich.edu", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),
    Database("homolog_msa", "Homolog MSA (NCBI/UniProt + MAFFT)", Category.STRUCTURE,
             "Viral homolog alignments → per-pocket JSD conservation (SPEC #1)",
             "https://mafft.cbrc.jp/alignment/software/", Status.INTEGRATED,
             api="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
             wired_in="vta.data.msa"),

    # 2. Ligand libraries -----------------------------------------------------
    Database("chembl", "ChEMBL", Category.LIGAND,
             "Bioactive molecules with measured activity (IC50, Ki)",
             "https://www.ebi.ac.uk/chembl", Status.INTEGRATED,
             api="https://www.ebi.ac.uk/chembl/api/data/molecule",
             wired_in="vta.data.ligands"),
    Database("pubchem", "PubChem", Category.LIGAND,
             "100M+ compounds; name→SMILES fallback when ChEMBL misses",
             "https://pubchem.ncbi.nlm.nih.gov", Status.INTEGRATED,
             api="https://pubchem.ncbi.nlm.nih.gov/rest/pug",
             wired_in="vta.data.pubchem"),
    Database("drugbank", "DrugBank", Category.LIGAND,
             "Approved/experimental drugs + targets + pharmacology",
             "https://www.drugbank.com", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),
    Database("zinc", "ZINC22", Category.LIGAND,
             "Purchasable compounds for virtual screening (billions)",
             "https://zinc.docking.org", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),
    Database("enamine_real", "Enamine REAL", Category.LIGAND,
             "Make-on-demand library (40B+ compounds)",
             "https://enamine.net/compound-collections/real-compounds",
             Status.INTEGRATED, wired_in="vta.data.external_databases"),
    Database("bindingdb", "BindingDB", Category.LIGAND,
             "Measured protein–ligand binding affinities",
             "https://www.bindingdb.org", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),

    # 3. Complexes / binding data --------------------------------------------
    Database("pdbbind", "PDBbind", Category.COMPLEX,
             "Curated protein–ligand complexes + affinities (scoring/ML)",
             "http://www.pdbbind.org.cn", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),
    Database("binding_moad", "BindingMOAD", Category.COMPLEX,
             "High-quality protein–ligand binding data",
             "http://www.bindingmoad.org", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),
    Database("pdbe_kb", "PDBe-KB", Category.COMPLEX,
             "Aggregated functional / binding-site annotations",
             "https://www.ebi.ac.uk/pdbe/pdbe-kb", Status.INTEGRATED,
             api="https://www.ebi.ac.uk/pdbe/graph-api/",
             wired_in="vta.data.external_databases"),

    # 4. Pockets --------------------------------------------------------------
    Database("scpdb", "scPDB", Category.POCKET,
             "Annotated druggable binding sites (pocket benchmark)",
             "http://bioinfo-pharma.u-strasbg.fr/scPDB", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),
    Database("castp", "CASTp", Category.POCKET,
             "Computed pocket geometry / topology",
             "http://sts.bioe.uic.edu/castp", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),

    # 5. Target / pathway / disease ------------------------------------------
    Database("open_targets", "Open Targets", Category.TARGET,
             "Target–disease associations with evidence scoring",
             "https://www.opentargets.org", Status.INTEGRATED,
             api="https://api.platform.opentargets.org/api/v4/graphql",
             wired_in="vta.data.external_databases"),
    Database("ttd", "Therapeutic Target Database", Category.TARGET,
             "Known and explored therapeutic targets",
             "https://db.idrblab.net/ttd", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),
    Database("kegg", "KEGG", Category.TARGET,
             "Pathways, diseases, drug interactions",
             "https://www.kegg.jp", Status.INTEGRATED,
             api="https://rest.kegg.jp", wired_in="vta.data.external_databases"),
    Database("gtopdb", "Guide to Pharmacology (IUPHAR/BPS)", Category.TARGET,
             "Drug targets + ligand pharmacology",
             "https://www.guidetopharmacology.org", Status.INTEGRATED,
             api="https://www.guidetopharmacology.org/services",
             wired_in="vta.data.external_databases"),

    # 6. ADMET ----------------------------------------------------------------
    # SwissADME is a web form with no public batch API; the pipeline does real ADMET
    # locally via admet_ai, while the external adapter records the reference source.
    Database("swissadme", "SwissADME", Category.ADMET,
             "ADME + drug-likeness (web form, no public API)",
             "http://www.swissadme.ch", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),
    Database("admetlab", "ADMETlab / admetSAR", Category.ADMET,
             "ADMET + toxicity predictions",
             "https://admetmesh.scbdd.com", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),
    Database("tox21", "Tox21 / ToxCast", Category.ADMET,
             "High-throughput toxicity screening data",
             "https://tripod.nih.gov/tox21", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),

    # 7. Benchmarks -----------------------------------------------------------
    Database("dude", "DUD-E", Category.BENCHMARK,
             "Decoy sets for docking / virtual-screening validation",
             "http://dude.docking.org", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),
    Database("lit_pcba", "LIT-PCBA", Category.BENCHMARK,
             "Unbiased experimental benchmark for VS / ML",
             "https://drugdesign.unistra.fr/LIT-PCBA", Status.INTEGRATED,
             wired_in="vta.data.external_databases"),

    # 8. Repurposing ----------------------------------------------------------
    Database("repurposing_hub", "Drug Repurposing Hub (Broad)", Category.REPURPOSING,
             "Curated repurposing candidates with annotations",
             "https://www.broadinstitute.org/drug-repurposing-hub",
             Status.INTEGRATED, wired_in="vta.data.external_databases"),
    Database("clinicaltrials", "ClinicalTrials.gov", Category.REPURPOSING,
             "Clinical-trial status of compounds",
             "https://clinicaltrials.gov", Status.INTEGRATED,
             api="https://clinicaltrials.gov/api/v2/studies",
             wired_in="vta.data.external_databases"),
)

# Fast lookup by stable key.
_BY_KEY: dict[str, Database] = {db.key: db for db in DATABASES}


def get(key: str) -> Database:
    """Return the database with this key, or raise KeyError."""
    return _BY_KEY[key]


def by_category(category: Category) -> list[Database]:
    """All databases in a category, catalog order preserved."""
    return [db for db in DATABASES if db.category is category]


def by_status(status: Status) -> list[Database]:
    """All databases at a given integration status."""
    return [db for db in DATABASES if db.status is status]


def integrated() -> list[Database]:
    """The data sources a node actually pulls from today."""
    return by_status(Status.INTEGRATED)


def summary() -> str:
    """One-line-per-source overview, grouped by category (for logs / reports)."""
    mark = {Status.INTEGRATED: "✅", Status.PLANNED: "◐", Status.CATALOGUED: "—"}
    lines: list[str] = []
    for cat in Category:
        members = by_category(cat)
        if not members:
            continue
        lines.append(f"\n{cat.value}:")
        for db in members:
            wired = f"  [{db.wired_in}]" if db.wired_in else ""
            lines.append(f"  {mark[db.status]} {db.name} — {db.provides}{wired}")
    return "\n".join(lines).lstrip("\n")
