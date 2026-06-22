"""vta.toolconfig — self-configuring discovery of external tool binaries.

For VTA-Agent to run drug design AUTONOMOUSLY (one command, no human setup), the
nodes must find their tools themselves instead of relying on a human to export
FPOCKET_BIN / VINA_BIN. Resolution order for each tool:

    1. explicit env var (override / CI)         e.g. VINA_BIN
    2. the tool on PATH                          shutil.which(name)
    3. known local build locations               <repo>/../tools/...

Returns None if the tool genuinely isn't present, so callers degrade to their
labelled fallback (the graceful-fallback discipline used across the pipeline).
"""
from __future__ import annotations

import os
import shutil
from pathlib import Path

# vta/toolconfig.py -> vta/ -> vta-agent/ -> vtaagent/  (tools live beside the repo)
_REPO = Path(__file__).resolve().parents[1]          # .../vta-agent
_WS = _REPO.parent                                    # .../vtaagent

# Known local build outputs (created during this project's tool builds).
_LOCAL = {
    "fpocket": [
        _WS / "tools" / "fpocket-src" / "bin" / "fpocket",
        _REPO / "tools" / "fpocket-src" / "bin" / "fpocket",
    ],
    "vina": [
        _WS / "tools" / "vina",
        _REPO / "tools" / "vina",
    ],
}


def find_tool(name: str, env_var: str | None = None) -> str | None:
    """Resolve an external tool binary path, or None if unavailable."""
    if env_var and os.environ.get(env_var):
        return os.environ[env_var]
    on_path = shutil.which(name)
    if on_path:
        return on_path
    for cand in _LOCAL.get(name, []):
        if cand.exists() and os.access(cand, os.X_OK):
            return str(cand)
    return None
