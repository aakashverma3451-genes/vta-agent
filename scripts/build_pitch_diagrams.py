"""Build the four competition-deck diagrams (BGI Research AI+X finals pitch).

Pure matplotlib (patches + patheffects), no seaborn/plotly. Values match the committed
benchmark artifacts already cited in docs/figures/vta_figure_set.html (Figs 1-5) and
docs/PROJECT_DOCUMENT.md — this script only re-renders them in the deck's brand palette,
it introduces no new numbers.

Usage:
    PYTHONPATH=. ../taxonagent/venv/bin/python scripts/build_pitch_diagrams.py
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle, Rectangle
from matplotlib.lines import Line2D

OUT = Path("docs/figures/pitch")
OUT.mkdir(parents=True, exist_ok=True)

# ── brand palette ──────────────────────────────────────────────────────────
PRIMARY = "#0099CC"
ACCENT = "#00AEEF"
BG = "#FFFFFF"
INK = "#1A1A2E"
GOOD = "#2E7D32"
GOOD_T = "#C8E6C9"
BAD = "#B71C1C"
BAD_T = "#FFCDD2"
AMBER = "#F57F17"
AMBER_T = "#FFF9C4"
BORDER = "#E0E0E0"
LABEL = "#546E7A"
GREY = "#B0BEC5"

plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["DejaVu Sans", "Helvetica", "Arial"]


def _fig(w_px: int, h_px: int, dpi: int = 100):
    fig = plt.figure(figsize=(w_px / dpi, h_px / dpi), dpi=dpi)
    fig.patch.set_facecolor(BG)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, w_px)
    ax.set_ylim(0, h_px)
    ax.invert_yaxis()  # work in "screen" coords: y=0 at top
    ax.axis("off")
    return fig, ax


def _rrect(ax, x, y, w, h, fc, ec=None, lw=1.2, radius=12, zorder=2):
    p = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad=0,rounding_size={radius}",
                        linewidth=lw, edgecolor=ec or "none", facecolor=fc, zorder=zorder,
                        mutation_aspect=1)
    ax.add_patch(p)
    return p


def _save(fig, name):
    path = OUT / name
    fig.savefig(path, dpi=300, facecolor=BG)
    plt.close(fig)
    print(f"wrote {path}")


# ─────────────────────────────────────────────────────────────────────────
# DIAGRAM 1 — problem / innovation two-column overview
# ─────────────────────────────────────────────────────────────────────────
def diagram1():
    W, H = 1920, 1080
    fig, ax = _fig(W, H)

    divider_x = W / 2
    ax.add_line(Line2D([divider_x, divider_x], [60, H - 100], color=ACCENT, linewidth=2.5, zorder=1))

    problems = [
        ("Fragmented Workflows", "Genome to candidates requires 6+ disconnected manual steps."),
        ("Inflated Benchmarks", "Most pipelines never test whether they beat a trivial 2D baseline."),
        ("False Confidence", "Docking scores reported as rankings with no uncertainty quantification."),
        ("No Honest Validation", "Negative results hidden; overfitting to easy decoy sets."),
    ]
    innovations = [
        ("Autonomous End-to-End", "Genome FASTA → ranked hypotheses in one agent loop."),
        ("Self-Auditing Triage", "Routes each target to the right method — or refuses."),
        ("Honest Validation", "Every metric: median + 95% CI, vs a 2D-similarity baseline."),
        ("Non-Removable Honesty Envelope", "Every report pinned to a frozen benchmark hash."),
    ]

    def _col(x0, title, items, dot_color):
        ax.text(x0, 90, title, fontsize=26, fontweight="bold", color=dot_color,
                 ha="left", va="center")
        y = 190
        card_h = 190
        for label, desc in items:
            Circle_ = Circle((x0 + 26, y + 30), 26, facecolor=dot_color, edgecolor="none", zorder=3)
            ax.add_patch(Circle_)
            ax.text(x0 + 26, y + 30, "!" if dot_color == BAD else "✓", fontsize=20,
                     color="white", fontweight="bold", ha="center", va="center", zorder=4)
            ax.text(x0 + 74, y + 14, label, fontsize=18.5, fontweight="bold", color=INK,
                     ha="left", va="center")
            ax.text(x0 + 74, y + 50, desc, fontsize=13.5, color=LABEL, ha="left", va="center",
                     wrap=True)
            ax.add_line(Line2D([x0, x0 + 780], [y + card_h - 30, y + card_h - 30],
                                color=BORDER, linewidth=1))
            y += card_h

    _col(90, "The Field's Problem", problems, BAD)
    _col(divider_x + 90, "What VTA-Agent Does", innovations, PRIMARY)

    # bottom bar
    _rrect(ax, 0, H - 90, W, 90, PRIMARY, radius=0)
    ax.text(W / 2, H - 45, "False-confidence rate:  100% (raw pipeline)   →   0% (full architecture)",
            fontsize=22, color="white", fontweight="bold", ha="center", va="center")

    _save(fig, "diagram1_problem_innovation.png")


# ─────────────────────────────────────────────────────────────────────────
# DIAGRAM 2 — the 6-node reasoning architecture pipeline
# ─────────────────────────────────────────────────────────────────────────
def diagram2():
    W, H = 1920, 1080
    fig, ax = _fig(W, H)

    # subtitle bar
    _rrect(ax, 0, 0, W, 70, PRIMARY, radius=0)
    ax.text(W / 2, 35, "Phase-R Reasoning Architecture  ·  "
                        "Autonomous routing — never changes docking weights",
            fontsize=19, color="white", fontweight="bold", ha="center", va="center")

    nodes = [
        ("INPUT\nCASCADE", "Genome → TaxonAgent → Structure\n(PDB/ESMFold/AlphaFold/Boltz-2)\n→ Pocket",
         "#E3F2FD", PRIMARY, INK, False),
        ("R1\nDOSSIER", "Provenance · pLDDT\nMetal dependence\nLigand class · Benchmarkability",
         "#E8F5E9", GOOD, INK, False),
        ("R2 TRIAGE\nROUTER", "◆ full_dock / annotate_only\n/ defer / refuse",
         PRIMARY, "white", "white", True),
        ("EXECUTOR", "Docking (Vina 1.2.5) OR\nAnnotation ranker +\n2D-sim baseline",
         "#FFF8E1", AMBER, INK, False),
        ("R4 VERIFICATION\nGATE", "Gate A: Redock RMSD < 2Å\nGate B: Beats 2D baseline\nGate C: AD check",
         "#FCE4EC", BAD, INK, False),
        ("HONESTY\nREPORT", "Pinned benchmark hash\n95% CI · disclaimer\nAudit trail",
         "#F3E5F5", "#6A1B9A", INK, False),
    ]
    principles = ["Know the source", "Know the target", "Choose the method",
                  "Execute honestly", "Gate every claim", "Ship with evidence"]

    n = len(nodes)
    node_w, node_h = 232, 240
    spine_y = 330
    margin = 70
    span = W - 2 * margin
    xs = [margin + span * (i + 0.5) / n for i in range(n)]

    for i in range(n - 1):
        x1 = xs[i] + node_w / 2
        x2 = xs[i + 1] - node_w / 2
        arr = FancyArrowPatch((x1, spine_y), (x2, spine_y), arrowstyle="-|>", mutation_scale=22,
                               linewidth=3, color=ACCENT, zorder=1,
                               shrinkA=2, shrinkB=2)
        ax.add_patch(arr)
        if i == 1:
            ax.text((x1 + x2) / 2, spine_y - 26, "decision", fontsize=10.5, style="italic",
                     color=LABEL, ha="center", zorder=5,
                     bbox=dict(facecolor=BG, edgecolor="none", pad=1))
        if i == 4:
            ax.text((x1 + x2) / 2, spine_y - 26, "verified\nclaim only", fontsize=10, style="italic",
                     color=LABEL, ha="center", zorder=5, linespacing=1.1,
                     bbox=dict(facecolor=BG, edgecolor="none", pad=1))

    for i, (label, sub, fc, hc, bc, solid) in enumerate(nodes):
        x = xs[i] - node_w / 2
        y = spine_y - node_h / 2
        _rrect(ax, x, y, node_w, node_h, fc, ec=(hc if not solid else "none"), lw=2, radius=16, zorder=3)
        ax.text(xs[i], y + 44, label, fontsize=15.5, fontweight="bold", color=hc,
                ha="center", va="center", zorder=4, linespacing=1.25)
        ax.text(xs[i], y + node_h / 2 + 28, sub, fontsize=10.8, color=bc,
                ha="center", va="center", zorder=4, linespacing=1.6)
        ax.text(xs[i], spine_y + node_h / 2 + 46, principles[i], fontsize=13, style="italic",
                color=LABEL, ha="center", va="center")

    # branches under node 3 (R2 triage router, index 2)
    bx = xs[2]
    by0 = spine_y + node_h / 2
    branches = [
        (bx - 230, "✓ full_dock", "High-confidence target", GOOD_T, GOOD),
        (bx, "⚠ annotate_only", "Out-of-domain", AMBER_T, AMBER),
        (bx + 230, "✗ defer / refuse", "Low-conf. / impossible", BAD_T, BAD),
    ]
    branch_y = by0 + 210
    for bxi, label, sub, fc, ec in branches:
        ax.add_line(Line2D([bx, bxi], [by0 + 100, branch_y - 30], color=LABEL, linewidth=1.3,
                            linestyle="--", zorder=1))
        _rrect(ax, bxi - 108, branch_y - 30, 216, 68, fc, ec=ec, lw=1.6, radius=10, zorder=3)
        ax.text(bxi, branch_y - 8, label, fontsize=13, fontweight="bold", color=ec,
                ha="center", va="center", zorder=4)
        ax.text(bxi, branch_y + 16, sub, fontsize=9.8, color=LABEL, ha="center", va="center", zorder=4)

    _save(fig, "diagram2_reasoning_architecture.png")


# ─────────────────────────────────────────────────────────────────────────
# DIAGRAM 3 — component ablation (bar chart + heatmap)
# ─────────────────────────────────────────────────────────────────────────
def diagram3():
    fig = plt.figure(figsize=(19.2, 9.0), dpi=100)
    fig.patch.set_facecolor(BG)

    levels = ["L0\nraw", "L1\n+dossier", "L2\n+triage", "L3\n+verify", "L4\nfull"]
    acc = [0, 0, 60, 100, 100]
    fc = [100, 100, 40, 0, 0]

    axA = fig.add_axes([0.06, 0.14, 0.40, 0.72])
    axA.set_facecolor(BG)
    xs = range(len(levels))
    width = 0.34
    b1 = axA.bar([x - width / 2 for x in xs], acc, width, color=GOOD, label="decision accuracy %", zorder=3)
    b2 = axA.bar([x + width / 2 for x in xs], fc, width, color=BAD, label="false-confidence %", zorder=3)
    axA.axhline(100, color=GOOD, linewidth=1.2, linestyle="--", zorder=2)
    axA.text(4.35, 102, "target (100%)", fontsize=10, color=GOOD, ha="right")
    axA.set_ylim(0, 112)
    axA.set_xticks(list(xs))
    axA.set_xticklabels(levels, fontsize=12)
    axA.set_ylabel("%", fontsize=13)
    axA.grid(axis="y", color="#F0F0F0", zorder=0)
    for spine in ("top", "right"):
        axA.spines[spine].set_visible(False)
    for bars in (b1, b2):
        for bar in bars:
            h = bar.get_height()
            color = "white" if h >= 25 else INK
            va = "center" if h >= 25 else "bottom"
            y = h / 2 if h >= 25 else h + 2
            axA.text(bar.get_x() + bar.get_width() / 2, y, f"{int(h)}", ha="center", va=va,
                      fontsize=10.5, fontweight="bold", color=color, zorder=4)
    axA.legend(loc="upper left", bbox_to_anchor=(0.0, 1.14), ncol=1, frameon=False, fontsize=11)
    axA.set_title("a  decision accuracy vs false-confidence, per level", fontsize=17,
                   fontweight="bold", color=INK, loc="left", pad=14)

    # ── B: heatmap ──
    axB = fig.add_axes([0.55, 0.14, 0.40, 0.72])
    axB.set_facecolor(BG)
    targets = ["NS5B", "GPX*", "GPY*", "Mpro", "PB1"]
    cols = ["L0", "L1", "L2", "L3", "L4"]
    grid = [
        ["F", "F", "C", "C", "C"],
        ["F", "F", "C", "C", "C"],
        ["F", "F", "C", "C", "C"],
        ["F", "F", "F", "C", "C"],
        ["F", "F", "F", "C", "C"],
    ]
    ncols, nrows = len(cols), len(targets)
    cw, ch = 1.0, 1.0
    for r in range(nrows):
        for c in range(ncols):
            val = grid[r][c]
            fcolor = GOOD_T if val == "C" else BAD_T
            ecolor = GOOD if val == "C" else BAD
            rect = Rectangle((c * cw, (nrows - 1 - r) * ch), cw, ch, facecolor=fcolor,
                              edgecolor=BORDER, linewidth=1.2)
            axB.add_patch(rect)
            axB.text(c * cw + cw / 2, (nrows - 1 - r) * ch + ch / 2, val, ha="center", va="center",
                      fontsize=13, fontweight="bold", color=ecolor)
    axB.set_xlim(0, ncols)
    axB.set_ylim(0, nrows)
    axB.set_xticks([c + 0.5 for c in range(ncols)])
    axB.set_xticklabels(cols, fontsize=12)
    axB.set_yticks([nrows - 1 - r + 0.5 for r in range(nrows)])
    axB.set_yticklabels(targets, fontsize=12)
    axB.tick_params(length=0)
    for spine in axB.spines.values():
        spine.set_visible(False)
    axB.set_title("b  which layer fixes which held-out target (disjoint sets)", fontsize=17,
                   fontweight="bold", color=INK, loc="left", pad=14)
    axB.text(0, -0.45, "*low-pLDDT predicted structures (refuse/defer)", fontsize=10.5,
              color=LABEL, ha="left", va="top")

    fig.text(0.5, 0.035,
              "Router fixes 3 targets (NS5B, GPX*, GPY*)  ·  Gate fixes 2 (Mpro, PB1)  ·  "
              "Neither alone is sufficient  ·  Run-to-run consistency: 1.0",
              fontsize=13, style="italic", color=LABEL, ha="center")

    _save(fig, "diagram3_ablation.png")


# ─────────────────────────────────────────────────────────────────────────
# DIAGRAM 4 — benchmark results (cliffs + BEDROC + ladder/table)
# ─────────────────────────────────────────────────────────────────────────
def diagram4():
    fig = plt.figure(figsize=(19.2, 9.0), dpi=100)
    fig.patch.set_facecolor(BG)

    fam, famv, famc = ["2-D", "Vina", "RF"], [0.67, 0.42, 0.39], [PRIMARY, BAD, BAD]

    def _mini_arena(ax, title=None):
        ax.set_facecolor(BG)
        for i, (name, v, c) in enumerate(zip(fam, famv, famc)):
            yy = 2 - i
            ax.plot([0.30, v], [yy, yy], color=BORDER, linewidth=1, zorder=1)
            ax.scatter([v], [yy], color=c, s=60, zorder=3)
            ax.text(v + 0.012, yy, f"{v:.2f}", fontsize=9.5, color=c, va="center", fontweight="bold")
            ax.text(0.285, yy, name, fontsize=9.5, va="center", ha="right", color=c)
        ax.axvline(0.50, color=GREY, linestyle="--", linewidth=1, zorder=1)
        ax.set_xlim(0.20, 0.78)
        ax.set_ylim(-0.7, 2.7)
        ax.axis("off")
        if title:
            ax.set_title(title, fontsize=10.5, color=INK, loc="left", pad=6)

    # ── A: cliff ranking accuracy dot-whisker (top-left) ──
    axA = fig.add_axes([0.075, 0.60, 0.255, 0.28])
    axA.set_facecolor(BG)
    methods = ["2-D kNN\nQSAR", "AutoDock\nVina", "RF-Score\n(learned)"]
    med = [0.666, 0.418, 0.392]
    lo = [0.64, 0.39, 0.36]
    hi = [0.70, 0.45, 0.42]
    colors = [PRIMARY, LABEL, LABEL]
    ys = [2, 1, 0]
    for y, m, l, h, c in zip(ys, med, lo, hi, colors):
        axA.errorbar(m, y, xerr=[[m - l], [h - m]], fmt="o", color=c, ecolor=c,
                      elinewidth=2, capsize=6, markersize=9, zorder=3)
    axA.axvline(0.50, color=BAD, linestyle="--", linewidth=1.3, zorder=1)
    axA.set_ylim(-0.6, 2.85)
    axA.text(0.50, 2.65, "chance 0.50", color=BAD, fontsize=9.5, ha="center", va="bottom")
    axA.set_yticks(ys)
    axA.set_yticklabels(methods, fontsize=10.5)
    axA.set_xlim(0.30, 0.72)
    axA.tick_params(axis="x", labelsize=9.5)
    axA.grid(axis="x", color="#F0F0F0", zorder=0)
    for spine in ("top", "right"):
        axA.spines[spine].set_visible(False)
    fig.text(0.075, 0.905, "a  Cliff-pair ranking accuracy, 1,193 pairs",
              fontsize=14, fontweight="bold", color=INK, ha="left")
    fig.text(0.075, 0.875, "(median, 95% CI)", fontsize=11, color=LABEL, ha="left")
    fig.text(0.075, 0.542, "ranking accuracy — is the more-potent analog ranked higher?",
              fontsize=9.5, color=LABEL, ha="left")

    # A mini panel (below)
    axA2 = fig.add_axes([0.075, 0.315, 0.255, 0.17])
    _mini_arena(axA2)
    fig.text(0.075, 0.495, "Fair arena (activity cliffs), ranking acc.", fontsize=10.5,
              color=INK, ha="left")
    fig.text(0.075, 0.29, "→ both docking methods below chance on activity cliffs",
              fontsize=10, style="italic", color=BAD, ha="left")

    # ── B: Mpro BEDROC bars (top-centre) ──
    axB = fig.add_axes([0.40, 0.60, 0.255, 0.28])
    axB.set_facecolor(BG)
    bmethods = ["2-D-sim", "Vina", "random"]
    bvals = [0.92, 0.68, 0.50]
    bcolors = [AMBER, PRIMARY, GREY]
    ybs = [2, 1, 0]
    axB.barh(ybs, bvals, height=0.5, color=bcolors, zorder=3)
    for y, v, c in zip(ybs, bvals, bcolors):
        axB.text(v + 0.02, y, f"{v:.2f}", va="center", fontsize=10.5, fontweight="bold", color=c)
    axB.set_yticks(ybs)
    axB.set_yticklabels(bmethods, fontsize=10.5)
    axB.set_xlim(0, 1.05)
    axB.tick_params(axis="x", labelsize=9.5)
    axB.grid(axis="x", color="#F0F0F0", zorder=0)
    for spine in ("top", "right"):
        axB.spines[spine].set_visible(False)
    fig.text(0.40, 0.905, "b  Enrichment (Moonshot Mpro),", fontsize=14, fontweight="bold",
              color=INK, ha="left")
    fig.text(0.40, 0.875, "BEDROC(α=20)", fontsize=11, color=LABEL, ha="left")
    fig.text(0.40, 0.542, "→ 2-D wins (expected on an analog-clustered set)", fontsize=9.5,
              style="italic", color=LABEL, ha="left")

    axB2 = fig.add_axes([0.40, 0.315, 0.255, 0.17])
    _mini_arena(axB2)
    fig.text(0.40, 0.495, "Fair arena (activity cliffs), ranking acc.", fontsize=10.5,
              color=INK, ha="left")

    # ── C: table (top-right) ──
    axC = fig.add_axes([0.695, 0.55, 0.285, 0.36])
    axC.set_facecolor(BG)
    axC.axis("off")
    axC.set_xlim(0, 10)
    axC.set_ylim(0, 10)

    rows = [
        ("SARS-CoV-2 Mpro", "50 measured Moonshot\ninactives (real)", "0.68\n[0.36, 0.89]", "powered", GOOD_T, GOOD),
        ("TiLV PB1", "ChEMBL property-\nmatched (4 actives)", "0.41\n[0.00, 1.00]", "underpowered", AMBER_T, AMBER),
        ("HCV NS5B NI", "24/43 recover\n0 matched decoys", "—", "un-benchmarkable", BAD_T, BAD),
    ]
    axC.add_patch(Rectangle((0, 9.0), 10, 0.9, facecolor=PRIMARY))
    heads = ["Target", "Control set", "BEDROC", "Grade"]
    hx = [0.15, 3.5, 6.35, 8.15]
    for hxi, h in zip(hx, heads):
        axC.text(hxi, 9.45, h, fontsize=10, fontweight="bold", color="white", va="center")
    ry = 8.55
    for name, ctrl, bedroc, grade, gfc, gec in rows:
        rh = 2.75
        axC.text(0.15, ry - rh / 2 + 0.55, name, fontsize=10, fontweight="bold", color=INK, va="center")
        axC.text(3.5, ry - rh / 2 + 0.3, ctrl, fontsize=8.4, color=LABEL, va="center", linespacing=1.6)
        axC.text(6.35, ry - rh / 2 + 0.3, bedroc, fontsize=8.6, color=INK, va="center", linespacing=1.6)
        axC.add_patch(FancyBboxPatch((7.95, ry - rh / 2 + 0.12), 2.0, 0.42,
                                      boxstyle="round,pad=0,rounding_size=0.08",
                                      facecolor=gfc, edgecolor="none"))
        axC.text(8.95, ry - rh / 2 + 0.33, grade, fontsize=7.6, fontweight="bold", color=gec,
                  ha="center", va="center")
        axC.add_line(Line2D([0, 10], [ry - rh, ry - rh], color=BORDER, linewidth=1))
        ry -= rh

    # ── ladder (bottom-right, own row — no overlap with table) ──
    axD = fig.add_axes([0.695, 0.20, 0.285, 0.28])
    axD.set_facecolor(BG)
    axD.axis("off")
    axD.set_xlim(0, 10)
    axD.set_ylim(0, 5)
    axD.text(0, 4.75, "the validation ladder", fontsize=13, fontweight="bold", color=INK, va="top")

    stages = ["3/4\ntop-5", "1/4\n(matched)", "CIs\n=[0,1]", "powered\nMpro", "Vina\n<2-D", "redock\n1.65 Å"]
    stage_fc = [BAD_T, BAD_T, AMBER_T, GOOD_T, GOOD_T, GOOD_T]
    stage_ec = [BAD, BAD, AMBER, GOOD, GOOD, GOOD]
    n = len(stages)
    lad_x0 = 0.15
    step_x = 1.55
    step_y = 0.28
    box_w, box_h = 1.3, 0.95
    lad_y0 = 1.9
    for i, (s, fc_, ec_) in enumerate(zip(stages, stage_fc, stage_ec)):
        x = lad_x0 + i * step_x
        y = lad_y0 + i * step_y
        axD.add_patch(FancyBboxPatch((x, y), box_w, box_h, boxstyle="round,pad=0,rounding_size=0.07",
                                      facecolor=fc_, edgecolor=ec_, linewidth=1.3))
        axD.text(x + box_w / 2, y + box_h / 2, s, fontsize=7.8, fontweight="bold", color=ec_,
                  ha="center", va="center", linespacing=1.15)
        if i < n - 1:
            x2 = lad_x0 + (i + 1) * step_x
            y2 = lad_y0 + (i + 1) * step_y
            arr = FancyArrowPatch((x + box_w, y + box_h / 2), (x2, y2 + box_h / 2),
                                   arrowstyle="-|>", mutation_scale=13, linewidth=1.6, color=ACCENT)
            axD.add_patch(arr)
    axD.text(0, 0.35,
              "inflated gate → matched decoys → CIs → real inactives\n"
              "→ honest negative → pose reliable", fontsize=8, color=LABEL, va="top", linespacing=1.5)

    _save(fig, "diagram4_benchmark.png")


def main():
    diagram1()
    diagram2()
    diagram3()
    diagram4()

    from PIL import Image
    print()
    for f in ["diagram1_problem_innovation.png", "diagram2_reasoning_architecture.png",
              "diagram3_ablation.png", "diagram4_benchmark.png"]:
        p = OUT / f
        assert p.exists() and p.stat().st_size > 0, f"missing or empty: {p}"
        img = Image.open(p)
        print(f, img.size)


if __name__ == "__main__":
    main()
