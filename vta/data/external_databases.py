"""Adapters for drug-design databases not on the default VTA path.

The core pipeline has dedicated modules for RCSB, ESMFold/AlphaFold, UniProt, MSA,
ChEMBL, and PubChem. This module integrates the remaining catalogued sources as
explicit adapters so future nodes can discover how to query them without hard-coding
URLs throughout the codebase.

Some sources are public APIs; others are web portals or licensed downloads. Those are
still represented with deterministic URL builders, but `fetch_json` only runs for
sources that expose a JSON endpoint.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Callable

import requests


@dataclass(frozen=True)
class ExternalAdapter:
    key: str
    name: str
    access: str
    base_url: str
    description: str
    build_url: Callable[[str], str]
    json_api: bool = False
    restricted: bool = False


def _quote(value: str) -> str:
    return requests.utils.quote((value or "").strip())


def _identity(base: str) -> Callable[[str], str]:
    return lambda _query: base


def _path(base: str) -> Callable[[str], str]:
    return lambda query: f"{base.rstrip('/')}/{_quote(query)}"


def _param(base: str, name: str) -> Callable[[str], str]:
    return lambda query: f"{base}?{name}={_quote(query)}"


ADAPTERS: dict[str, ExternalAdapter] = {
    "pdbe": ExternalAdapter(
        "pdbe", "PDBe", "public_json",
        "https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary",
        "PDB mirror and annotation APIs keyed by PDB id.",
        _path("https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary"), True),
    "swissmodel": ExternalAdapter(
        "swissmodel", "SWISS-MODEL Repository", "public_json",
        "https://swissmodel.expasy.org/repository/uniprot",
        "Homology-model lookup keyed by UniProt accession.",
        _path("https://swissmodel.expasy.org/repository/uniprot"), True),
    "opm": ExternalAdapter(
        "opm", "OPM", "public_web",
        "https://opm.phar.umich.edu/proteins",
        "Membrane-protein orientation records keyed by OPM/PDB search terms.",
        _param("https://opm.phar.umich.edu/proteins", "search")),
    "drugbank": ExternalAdapter(
        "drugbank", "DrugBank", "licensed_download",
        "https://go.drugbank.com",
        "Approved and investigational drug metadata; full data requires a license.",
        _param("https://go.drugbank.com/unearth/q", "searcher"), False, True),
    "zinc": ExternalAdapter(
        "zinc", "ZINC22", "public_web",
        "https://zinc.docking.org",
        "Purchasable compound search and tranche downloads.",
        _param("https://zinc.docking.org/substances/home", "q")),
    "enamine_real": ExternalAdapter(
        "enamine_real", "Enamine REAL", "download_or_vendor",
        "https://enamine.net/compound-collections/real-compounds",
        "Make-on-demand library; production use usually consumes vendor files.",
        _identity("https://enamine.net/compound-collections/real-compounds")),
    "bindingdb": ExternalAdapter(
        "bindingdb", "BindingDB", "public_webservice",
        "https://www.bindingdb.org/rwd/bind/chemsearch/marvin/SDFdownload.jsp",
        "Measured binding data export; query shape depends on ligand or target mode.",
        _param("https://www.bindingdb.org/rwd/bind/chemsearch/marvin/SDFdownload.jsp",
               "monomerid")),
    "pdbbind": ExternalAdapter(
        "pdbbind", "PDBbind", "manual_download",
        "http://www.pdbbind.org.cn",
        "Curated complexes and affinities; datasets are downloaded manually.",
        _identity("http://www.pdbbind.org.cn")),
    "binding_moad": ExternalAdapter(
        "binding_moad", "BindingMOAD", "public_download",
        "http://www.bindingmoad.org",
        "Protein-ligand complex annotations distributed as downloadable files.",
        _identity("http://www.bindingmoad.org")),
    "pdbe_kb": ExternalAdapter(
        "pdbe_kb", "PDBe-KB", "public_json",
        "https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages",
        "Aggregated functional and binding-site annotations keyed by PDB id.",
        _path("https://www.ebi.ac.uk/pdbe/graph-api/pdbe_pages"), True),
    "scpdb": ExternalAdapter(
        "scpdb", "scPDB", "public_download",
        "http://bioinfo-pharma.u-strasbg.fr/scPDB",
        "Annotated druggable pockets; consumed as downloaded benchmark files.",
        _identity("http://bioinfo-pharma.u-strasbg.fr/scPDB")),
    "castp": ExternalAdapter(
        "castp", "CASTp", "public_web",
        "http://sts.bioe.uic.edu/castp",
        "Computed pocket geometry and topology from submitted structures.",
        _identity("http://sts.bioe.uic.edu/castp")),
    "open_targets": ExternalAdapter(
        "open_targets", "Open Targets", "public_graphql",
        "https://api.platform.opentargets.org/api/v4/graphql",
        "Target-disease association GraphQL API.",
        _identity("https://api.platform.opentargets.org/api/v4/graphql"), True),
    "ttd": ExternalAdapter(
        "ttd", "Therapeutic Target Database", "public_web",
        "https://db.idrblab.net/ttd",
        "Known and explored therapeutic targets.",
        _param("https://db.idrblab.net/ttd/search/ttd/target", "search_api_fulltext")),
    "kegg": ExternalAdapter(
        "kegg", "KEGG", "public_text",
        "https://rest.kegg.jp/find/drug",
        "Pathway, disease, and drug REST records.",
        _path("https://rest.kegg.jp/find/drug")),
    "gtopdb": ExternalAdapter(
        "gtopdb", "Guide to Pharmacology (IUPHAR/BPS)", "public_json",
        "https://www.guidetopharmacology.org/services",
        "Drug target and ligand pharmacology service endpoints.",
        _path("https://www.guidetopharmacology.org/services/ligands"), True),
    "swissadme": ExternalAdapter(
        "swissadme", "SwissADME", "web_only",
        "http://www.swissadme.ch",
        "ADME and drug-likeness web tool; no public batch API.",
        _identity("http://www.swissadme.ch")),
    "admetlab": ExternalAdapter(
        "admetlab", "ADMETlab / admetSAR", "web_or_package",
        "https://admetmesh.scbdd.com",
        "ADMET and toxicity web predictions; VTA uses local admet_ai for automation.",
        _identity("https://admetmesh.scbdd.com")),
    "tox21": ExternalAdapter(
        "tox21", "Tox21 / ToxCast", "public_portal",
        "https://tripod.nih.gov/tox21",
        "High-throughput toxicity screening portal and downloads.",
        _identity("https://tripod.nih.gov/tox21")),
    "dude": ExternalAdapter(
        "dude", "DUD-E", "public_download",
        "http://dude.docking.org",
        "Docking decoy sets for retrospective validation.",
        _identity("http://dude.docking.org")),
    "lit_pcba": ExternalAdapter(
        "lit_pcba", "LIT-PCBA", "public_download",
        "https://drugdesign.unistra.fr/LIT-PCBA",
        "Less-biased virtual-screening benchmark datasets.",
        _identity("https://drugdesign.unistra.fr/LIT-PCBA")),
    "repurposing_hub": ExternalAdapter(
        "repurposing_hub", "Drug Repurposing Hub (Broad)", "public_download",
        "https://www.broadinstitute.org/drug-repurposing-hub",
        "Curated repurposing candidates and annotations.",
        _identity("https://www.broadinstitute.org/drug-repurposing-hub")),
    "clinicaltrials": ExternalAdapter(
        "clinicaltrials", "ClinicalTrials.gov", "public_json",
        "https://clinicaltrials.gov/api/v2/studies",
        "Clinical-trial status for compounds and interventions.",
        _param("https://clinicaltrials.gov/api/v2/studies", "query.term"), True),
}


def get_adapter(key: str) -> ExternalAdapter:
    """Return the adapter for a catalog key, or raise KeyError."""
    return ADAPTERS[key]


def build_url(key: str, query: str = "") -> str:
    """Build the source-specific lookup URL for `query`."""
    return get_adapter(key).build_url(query)


def fetch_json(key: str, query: str = "", *, timeout: int = 30) -> dict[str, Any] | None:
    """Fetch JSON for adapters that expose a JSON endpoint; return None otherwise."""
    adapter = get_adapter(key)
    if not adapter.json_api or adapter.restricted:
        return None
    r = requests.get(adapter.build_url(query), timeout=timeout)
    if r.status_code == 404:
        return None
    r.raise_for_status()
    return r.json()


def adapter_summary() -> list[dict[str, str]]:
    """Serializable overview for docs, reports, or health checks."""
    return [
        {"key": a.key, "name": a.name, "access": a.access, "base_url": a.base_url}
        for a in ADAPTERS.values()
    ]
