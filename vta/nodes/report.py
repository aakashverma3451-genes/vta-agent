"""report_node — final graph node: emit the self-contained HTML report.

Runs at the end of BOTH paths (proceed and defer), so every run produces a report
(a deferred run's report just shows the classification + defer reason, no leads).
Thin wrapper over vta.report.write_report; the rendering logic lives there.
"""
from __future__ import annotations

from vta.report import write_report
from vta.state import VTAState


def report_node(state: VTAState) -> VTAState:
    path = write_report(state)
    state["audit_trail"].append(f"Report: wrote {path}")
    return state
