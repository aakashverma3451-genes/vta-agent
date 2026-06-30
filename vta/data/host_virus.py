"""Host-virus adapter metadata; real ingestion remains opt-in."""
from __future__ import annotations

ADAPTERS = {
    "VirHostNet": {"license": "source-specific", "access": "manual_download"},
    "viruses.STRING": {"license": "STRING terms", "access": "api_or_download"},
    "Reactome": {"license": "CC BY 4.0", "access": "open"},
    "SIGNOR": {"license": "CC BY 4.0", "access": "open"},
    "HuRI": {"license": "source-specific", "access": "open_with_terms"},
    "DepMap": {"license": "CC BY 4.0", "access": "open"},
}


def host_virus_adapters() -> dict:
    return ADAPTERS
