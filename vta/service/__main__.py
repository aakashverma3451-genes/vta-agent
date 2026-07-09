"""Launch the OpenAI-compatible adapter:  python -m vta.service  (host/port via env).

  VTA_HOST (default 0.0.0.0) · VTA_PORT (default 8000) · VTA_API_KEY (optional; require a key)

Then point Open WebUI at  http://<host>:8000/v1  (any API key) and pick model `vta-agent-triage`.
"""
from __future__ import annotations

import os

import uvicorn


def main() -> None:
    uvicorn.run("vta.service.openai_adapter:app",
                host=os.environ.get("VTA_HOST", "0.0.0.0"),
                port=int(os.environ.get("VTA_PORT", "8000")), log_level="info")


if __name__ == "__main__":
    main()
