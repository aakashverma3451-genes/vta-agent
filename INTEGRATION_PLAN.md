# VTA-Agent — Accuracy Integration Plan

**What this is:** a sequenced, reasoned plan for adding the accuracy-improving tools
from *"Research Papers That Can Improve Accuracy"* and *"Molecular Dynamics
Implementation Guide"* into VTA-Agent — grounded in (a) what the pipeline actually is
today, and (b) the real environment (Apple-Silicon arm64, no GPU, where every external
binary so far has needed a build saga). Each item says **what it is, why it helps, how
it plugs in, the benefit, the effort/risk, and how we'll verify it.**

> Status today: `classify → router → structure → pockets → dock(Vina) → rank(LE-led)
> → report`. Validation gate **PASSES** (controls recover). 16 tests, all real tools
> except a few placeholders (conservation).

---

## 1. The guiding idea — a validated discovery funnel

The papers all plug into one mental model: a funnel that gets **narrower and more
expensive** at each stage. Cheap, broad methods screen; expensive, accurate methods
validate. VTA-Agent already implements the top of the funnel; these tools extend it
downward and sideways.

```
        ChEMBL library            ──►  (broad, cheap)
   [pockets] → [Vina dock]        ──►  we are here
   [LE-led rank] → top-20         ──►  we are here
        │
        ├─►  ADMET filter         (Phase 3a) drug-likeness, per lead     ⚡ cheap
        ├─►  DL re-score          (Phase 3b) GNINA/Boltz re-rank top-N   🔧 GPU-ish
        ├─►  pocket consensus     (Phase 3c) P2Rank cross-check          🟢 cheap
        │
        ▼
   top-5 leads
        │
        └─►  MD validation        (Phase 4) 100 ns stability + MM-GBSA   🏗️ GPU-days
                │
                ▼
   top-3 leads
        └─►  FEP                   (Phase 5) absolute ΔG, paper-grade     🏗️ expert
```

**Design principle (unchanged):** every new capability is a node that **auto-detects
its dependency and degrades gracefully to a labelled fallback** — exactly how
`structure`/`pockets`/`docking` already behave. This keeps CI hermetic (no GPU, no
external binaries needed for tests) while the real path runs wherever the tool exists.
The §4.1.1-style contract means each addition is additive: it writes new fields, never
breaks existing ones.

---

## 2. The reusable integration pattern (how every item below lands)

We've now done this four times (ESMFold, FPocket, Vina, ChEMBL). Each integration is
the same five moves — naming them once so the per-item sections stay short:

1. **Probe first.** Confirm the tool installs/runs on *this* box before building
   against it (the ESMFold ceiling, the fpocket arm64 build, the Vina Rosetta binary
   were all found this way). No building against an assumed dependency.
2. **Isolate the dependency** behind a `_tool_available()` / `_BIN` seam so tests
   monkeypatch it and stay offline.
3. **New node, additive state.** Add a node that writes *new* `VTAState` fields; never
   mutate the meaning of existing ones.
4. **Graceful fallback**, labelled in the audit trail and `versions` (`real` vs
   `MOCK`/`skipped`), so a missing GPU/binary never crashes the graph.
5. **Verify**: a hermetic unit test (mock path) + one real end-to-end check + a
   validation-gate re-run when it touches ranking.

---

## 3. Phase 3a — ADMET-AI (the quick win) ⚡

**Idea.** After ranking, predict each lead's drug-likeness/safety profile from its
SMILES with ADMET-AI (a Chemprop-RDKit GNN, #1 on the TDC ADMET leaderboard) and
surface it. We already carry real canonical SMILES on every lead, so this is a clean
bolt-on.

**Why it helps.** Binding affinity alone is not drug discovery. A compound that binds
tightly but is hERG-toxic or non-bioavailable is a dead end. ADMET adds an orthogonal
quality axis the current pipeline is blind to — and it's the single most *visible*
upgrade in a demo (new safety columns next to the leads).

**How it plugs in.**
- New node `admet_node` after `rank`, before `report`:
  `… → rank → admet → report`.
- State: each `lead_candidate` gains `admet: {herg, bioavailability, solubility, …}`;
  optionally an `admet_flag` (e.g. "hERG risk") for quick scanning.
- Report: `vta/report.py` `_leads()` gains columns (hERG / oral-avail / solubility),
  with at-risk cells tinted. A short methodology note like the LE one.
- Dependency seam: `_admet_available()` (try `import admet_ai`); fallback writes
  `admet: None` + audit `"[skip] ADMET: admet-ai not installed"`.

**Benefit.** Turns "docking pipeline" into "drug-discovery platform" — a whole new
decision dimension, near-zero architectural risk, big presentation payoff.

**Effort / deps.** ⚡ Drop-in. `pip install admet-ai` (pulls chemprop + torch — CPU
fine for ~20 leads). The only risk is the PyPI/TLS flakiness we've hit; torch is a
large download.

**Environment fit.** ✅ CPU-only, arm64-friendly. No binary build.

**Verify.** Unit test with a stubbed predictor (hermetic). One real call on the 9
cached ligands → confirm sensible values (e.g. ribavirin high solubility). Report
renders the new columns. No gate re-run needed (ADMET doesn't change ranking unless we
choose to weight it — keep it a *filter/annotation* first, not a score term).

**Citation status.** Tool real & verifiable (Swanson et al., *Bioinformatics* 2024,
`pip install admet-ai`).

---

## 4. Phase 3b — Deep-learning re-scoring (GNINA / Boltz-2 / Boltzina) 🔧

**Idea.** Re-score the existing Vina poses with a learned scoring function and blend
that into the ranking. Three candidates, same slot:
- **GNINA 1.3** — Vina + CNN rescoring, same PDBQT input.
- **Boltz-2 / Boltzina** — DL affinity prediction / docking-guided re-score.

**Why it helps.** Scoring is VTA-Agent's most sensitive layer — we already proved that
(ΔG-led → LE-led flipped the validation gate). A learned scorer that correlates better
with experiment than Vina's empirical function is the highest-ceiling accuracy gain in
the docs. The coronavirus-Mpro benchmark (§1.4) is the read that tells us which is
actually best for *viral* targets before we commit.

**How it plugs in.**
- Either swap the engine in `docking.py` (GNINA: same PDBQT, change `_vina_bin` →
  `_gnina_bin`, parse CNN score) **or** add a `rescore_node` after `dock` that takes
  the Vina poses and adds a `dl_score` per record.
- State: docking records gain `cnn_score` / `dl_affinity`; `rank.py` adds a weighted
  term (and we **re-run the validation gate** to recalibrate weights, exactly as we did
  for LE).
- Fallback: no GNINA/GPU → keep Vina score, label `versions[docking]` accordingly.

**Benefit.** Potentially the biggest correlation-with-truth improvement available;
keeps our funnel intact (re-scores, doesn't replace).

**Effort / deps.** 🔧 Medium-to-saga **on this box.** GNINA needs the `gnina` binary
(CUDA CNN — really wants a GPU); Boltz-2/Boltzina need GPU + model weights. The docs
call GNINA "drop-in," but that assumes a Linux/GPU host. **On arm64-no-GPU this is not
drop-in.**

**Environment fit.** ⚠️ Poor locally. Best done on a GPU host/Colab, or deferred.
Recommendation: **read §1.4 benchmark now** (free), build the `rescore_node` *seam*
now (so it's ready), but run the real DL scorer on a GPU later.

**Verify.** Benchmark-informed choice → gate re-run with the DL term → confirm controls
still pass and ideally separate *better* than LE alone.

**Citation status.** GNINA, Boltz-2 real. Boltzina (arXiv 2025) and the 2026 benchmark
PMCs — verify before citing.

---

## 5. Phase 3c — P2Rank pocket consensus 🟢

**Idea.** Run P2Rank (ML pocket predictor) alongside FPocket; treat pockets both tools
agree on as high-confidence, flag disagreements.

**Why it helps (and the honest caveat).** For *predicted* structures (TiLV orphan
segments) pocket quality is uncertain, and a consensus reduces false pockets. **But**
for our *validated* target we already bypass blind pocket-finding with the
experimental active site (8PSO bound-CTP), so the marginal value **right now** is low.
It matters once we screen proteins without a known site.

**How it plugs in.** `pockets_node` already has the `EXPERIMENTAL_ACTIVE_SITE` →
FPocket → mock cascade; add P2Rank as a parallel detector and a consensus step before
FPocket. State: pockets gain `consensus: bool` / `detectors: [...]`.

**Effort / deps.** ⚡ Java, runs on PDBs like FPocket; arm64-OK. Low risk.

**Environment fit.** ✅ Fine (needs a JRE).

**Verify.** On 8PSO, confirm P2Rank's top pocket is near the known NTP site (a second
independent vote for our active-site choice).

**Citation status.** P2Rank real (Krivák & Hoksza 2018).

---

## 6. Phase 4 — Molecular dynamics validation (Level 2) 🏗️

**Idea.** Take the top-5 leads, run ~100 ns explicit-solvent MD per complex (OpenMM),
and extract three signals — **pose stability (RMSD), MM-GBSA binding energy, persistent
contacts** — then re-rank with that evidence. This is the second doc in full.

**Why it helps.** Vina is a photograph; MD is a video. A great-looking static pose can
fall out of the pocket in 10 ns; MD distinguishes real binders from docking artifacts.
The prior TiLV paper (Sumon et al. 2023) did MD — a reviewer *will* expect it, so this
is the line between "hackathon project" and "publishable."

**How it plugs in.** Exactly the second doc's design — a **separate, opt-in Phase B**
so it never blocks the fast pipeline:
```
… rank → md_select → md_simulate → md_analyze → md_rerank → report
build_app(include_md=True)
```
- New nodes: `md_select` (top-5), `md_simulate` (OpenMM 100 ns), `md_analyze`
  (MDAnalysis RMSD/contacts + MM-GBSA), `md_rerank` (STABLE boost / UNSTABLE penalty).
- State: `md_candidates`, `md_results`, `md_analysis`, `md_validated_leads`.
- Report: a new "MD-validated leads" section (RMSD verdict, MM-GBSA, key contacts).
- Fallback: `_check_openmm()` already in the doc's skeleton → labelled `[MOCK]`.

**Benefit.** The credibility capstone: "Ribavirin stayed bound 100 ns (RMSD < 2 Å),
MM-GBSA −28 kcal/mol, H-bonds to the RdRp catalytic residues." That's paper-grade.

**Effort / deps.** 🏗️ Major — ~400 lines / 4 nodes (skeletons already drafted in the
guide and they match our architecture), **plus 24–48 GPU-hours** for 5 × 100 ns. CPU
can run a 10 ns Level-1 smoke test only.

**Environment fit.** ❌ Not on this box for production (no GPU). **Build the nodes now
(opt-in, mock-tested), run for real on a GPU host later.** A Level-1 (10 ns) CPU run is
feasible as a demo-of-the-mechanism.

**Verify.** Mock-path test (OpenMM absent → labelled skip, graph still completes). One
Level-1 10 ns CPU run on the top control to prove the wiring. Production Level-2 on GPU.

**Citation status.** OpenMM/OpenMMDL/MDAnalysis/MM-GBSA all real. Match Sumon et al.
2023 (verify the exact reference).

---

## 7. Phase 5 — FEP / absolute binding free energy 🏗️ (paper only)

**Idea.** For the final top-3, compute absolute binding free energies (ABFE) via
alchemical MD (OpenFE / Boltz-ABFE). RMS error ~1–2 kcal/mol vs experiment.

**Why / fit.** Publication-grade quantitative claims only. Requires alchemical-MD
expertise and serious compute. **Out of scope until post-submission** — listed for
completeness so the funnel is honest end-to-end.

---

## 8. Sequencing & recommendation

| # | Phase | Improves | Effort *here* | When | Gate re-run? |
|---|-------|----------|---------------|------|--------------|
| 1 | **3a ADMET-AI** | new safety/PK axis on leads | ⚡ CPU drop-in | **now** | no |
| 2 | 3c P2Rank | pocket confidence (orphans) | ⚡ JRE | opportunistic | no |
| 3 | 1.4 benchmark read | which DL scorer to pick | ⚡ reading | now (free) | n/a |
| 4 | 3b DL re-score seam | scoring ceiling | 🔧 seam now, GPU later | seam now | **yes (later)** |
| 5 | 4 MD nodes (mock) | pose validation wiring | 🔧 build | this week | no (until run) |
| 6 | 4 MD production | publishable validation | 🏗️ GPU-days | GPU host | rerank |
| 7 | 5 FEP | quantitative ΔG | 🏗️ expert | paper revision | n/a |

**Bottom line:** do **ADMET-AI now** (real win, real fit, zero GPU). Build the
**DL-rescore and MD seams** so the architecture is ready, but run their heavy compute
on a GPU host. P2Rank is a cheap opportunistic add. FEP is post-submission.

**Hackathon vs paper split:**
- *Hackathon (this week):* ADMET-AI columns + (optional) P2Rank consensus + a Level-1
  10 ns MD smoke run as a "validation preview." All demoable on this laptop.
- *Paper:* GPU DL-rescoring, full Level-2 MD on top-5, FEP on top-3, conservation calc.

---

## 9. Risks & honesty notes

- **Environment is the real constraint, not the code.** GNINA/Boltz/MD all assume
  Linux+GPU; the docs' "drop-in" labels are host-dependent. We've already paid this
  tax three times (fpocket build, Vina Rosetta, ChEMBL throttle) — plan for it.
- **Citations:** the load-bearing tools are real and verifiable; several 2025–2026
  blog/PMC items in the source docs are not yet verifiable — treat as leads, verify
  before any paper cites them.
- **Don't let new score terms drift the gate silently.** Anything touching `rank.py`
  (DL score, MD rerank, ADMET-as-score) must re-run `scripts/validate_controls.py` so
  we recalibrate against the controls — the same discipline that caught the size bias.
- **Keep ADMET an annotation first, a score term second.** Adding it straight into the
  composite would re-open the calibration question; surface it as columns/flags first,
  decide on weighting deliberately.
