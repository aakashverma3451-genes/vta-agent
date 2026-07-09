"""OpenAI-compatible adapter — lets Open WebUI (or any OpenAI client) drive VTA-Agent as a model.

Exposes the two endpoints Open WebUI needs to treat this as an OpenAI provider:
  GET  /v1/models             — advertises the `vta-agent-triage` model
  POST /v1/chat/completions   — runs the honest triage (streaming or not)
plus GET /health.

Auth is intentionally permissive: Open WebUI always sends an `Authorization: Bearer <key>`, and
we accept any value (set VTA_API_KEY to require a specific one). No key is validated by default
so a local deployment "just works"; put it behind a reverse proxy for anything exposed.

The turn is offline and fast — it runs the reasoning layer + honesty envelope (see
`triage_chat.run_triage`), never docking or network folding, so it can't hang a chat.
"""
from __future__ import annotations

import json
import os
import time
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, Header, HTTPException
from fastapi.responses import StreamingResponse
from pydantic import BaseModel

from vta.service.triage_chat import run_triage

MODEL_ID = "vta-agent-triage"
app = FastAPI(title="VTA-Agent OpenAI-compatible adapter", version="1.0")


class Message(BaseModel):
    role: str
    content: Any = ""


class ChatRequest(BaseModel):
    model: str = MODEL_ID
    messages: List[Message] = []
    stream: bool = False


def _check_auth(authorization: Optional[str]) -> None:
    required = os.environ.get("VTA_API_KEY")
    if required and (authorization or "").removeprefix("Bearer ").strip() != required:
        raise HTTPException(status_code=401, detail="invalid api key")


def _last_user_text(messages: List[Message]) -> str:
    for m in reversed(messages):
        if m.role == "user":
            c = m.content
            if isinstance(c, list):  # OpenAI content-parts form
                return " ".join(p.get("text", "") for p in c if isinstance(p, dict))
            return str(c or "")
    return ""


@app.get("/health")
def health() -> Dict[str, str]:
    return {"status": "ok", "model": MODEL_ID}


@app.get("/v1/models")
def list_models(authorization: Optional[str] = Header(default=None)) -> Dict[str, Any]:
    _check_auth(authorization)
    return {"object": "list", "data": [{
        "id": MODEL_ID, "object": "model", "created": 0, "owned_by": "vta-agent",
    }]}


def _chunk(delta: Dict[str, Any], finish: Optional[str] = None) -> str:
    payload = {"id": "vta-triage", "object": "chat.completion.chunk", "created": int(time.time()),
               "model": MODEL_ID, "choices": [{"index": 0, "delta": delta, "finish_reason": finish}]}
    return f"data: {json.dumps(payload)}\n\n"


@app.post("/v1/chat/completions")
def chat_completions(req: ChatRequest,
                     authorization: Optional[str] = Header(default=None)):
    _check_auth(authorization)
    result = run_triage(_last_user_text(req.messages))
    # A short, honest "thinking" preamble (the pipeline's own audit lines) then the verdict.
    steps = result.get("steps") or []
    preamble = "".join(f"> {s}\n" for s in steps)
    body = (preamble + "\n" if preamble else "") + result["markdown"]

    if not req.stream:
        return {
            "id": "vta-triage", "object": "chat.completion", "created": int(time.time()),
            "model": MODEL_ID,
            "choices": [{"index": 0, "message": {"role": "assistant", "content": body},
                         "finish_reason": "stop"}],
            "usage": {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0},
        }

    def stream():
        yield _chunk({"role": "assistant"})
        # stream the thinking lines first, then the verdict paragraph by paragraph
        for s in steps:
            yield _chunk({"content": f"> {s}\n"})
        yield _chunk({"content": "\n"})
        for block in result["markdown"].split("\n\n"):
            yield _chunk({"content": block + "\n\n"})
        yield _chunk({}, finish="stop")
        yield "data: [DONE]\n\n"

    return StreamingResponse(stream(), media_type="text/event-stream")
