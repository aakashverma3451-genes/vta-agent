"""report_node — final graph node: emit the self-contained HTML report.

Runs at the end of BOTH paths (proceed and defer), so every run produces a report
(a deferred run's report just shows the classification + defer reason, no leads).
Thin wrapper over vta.report.write_report; the rendering logic lives there.
"""
from __future__ import annotations

from vta.report import write_report
from vta.report_envelope import build_envelope
from vta.state import VTAState


def report_node(state: VTAState) -> VTAState:
    # Build the honesty envelope (D0.5) BEFORE writing the report and stash it on the state, so
    # the same object is available to the (future) API response and run manifest — not just the
    # HTML. render_report renders it unconditionally, so a run can never emit a naked ranking.
    state["honesty_envelope"] = build_envelope(state)
    path = write_report(state)
    state["audit_trail"].append(f"Report: wrote {path}")
    return state
