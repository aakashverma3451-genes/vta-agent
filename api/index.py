"""Vercel serverless entry — exposes the VTA-Agent OpenAI-compatible adapter as an ASGI app.

Vercel's Python runtime serves the module-level `app`. We anchor sys.path and chdir to the repo
root so the committed benchmark artifacts (read via relative paths in report_envelope /
triage_chat) resolve on the serverless filesystem — those files are bundled by
`vercel.json` → functions.includeFiles = "outputs/**".

The triage path is intentionally light (no rdkit/numpy/langgraph), so the whole function fits
well under Vercel's size limit and cold-starts fast.
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
os.chdir(ROOT)

from vta.service.openai_adapter import app  # noqa: E402  (must follow sys.path/chdir setup)

# `app` is what Vercel serves.
__all__ = ["app"]
