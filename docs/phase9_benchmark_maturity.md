# Phase 9 Benchmark Maturity Status

## Current Conclusion

Phase 9 confirms the Phase 8 interpretation: the TiLV PB1 benchmark is scientifically
useful but not statistically mature. The active count, not the decoy count, is now the
main blocker.

## 9C: Confidence Intervals

Bootstrap confidence intervals were added to `vta/eval/metrics.py` and applied to the
current TiLV PB1 matched benchmark.

Artifact:

- `outputs/phase8/benchmark_ci_tilv_pb1.md`
- `outputs/phase8/matched_benchmark_tilv_pb1.json`

Current TiLV PB1 matched benchmark:

- N = 54: 4 actives / 50 ChEMBL property-matched presumed decoys.
- BEDROC(alpha=20): point 0.4062; bootstrap median 0.4062, 95% CI [0.0009, 1.0].
- EF1%: point 13.5; bootstrap median 9.0, 95% CI [0.0, 39.825].
- logAUC: point 0.3883; bootstrap median 0.3877, 95% CI [0.1327, 1.0].
- ROC-AUC: point 0.785; bootstrap median 0.7921, 95% CI [0.5976, 1.0].

Interpretation: the intervals are too wide for a publication-grade enrichment claim.
This is expected with only four actives.

## 9A: Expanded Actives

New actives fetcher:

- `scripts/fetch_phase9_actives.py`

Generated actives artifacts:

- `vta/data/phase9_actives/hcv_ns5b.json`
  - Target: ChEMBL `CHEMBL5375`, Hepatitis C virus NS5B RNA-dependent RNA polymerase.
  - Unique ChEMBL actives: 52.
  - Status after stratification: not directly usable as one pooled benchmark.
- `vta/data/phase9_actives/influenza_pa_pb1.json`
  - Target: ChEMBL `CHEMBL3137263`, influenza PA/PB1 interaction.
  - Unique ChEMBL actives: 2.
  - Status: insufficient for a mature orthomyxo/PB1-class benchmark.
- `vta/data/phase9_actives/hcv_ns5b_stratified.json`
  - Stratifies CHEMBL5375 rows by mechanism/binding-site evidence.
  - Counts: 1 NI, 6 NNI, 8 direct NS5B site-unknown, 35 indirect antiviral assay,
    2 off-target HCV drugs.
  - Status: confirms that pooled CHEMBL5375 actives must not be docked into one pocket.
- `vta/data/phase9_actives/hcv_ns5b_ni_active_site.json`
  - Extracted from ChEMBL `CHEMBL4296320`, the RNA-directed RNA polymerase mechanism
    target used by sofosbuvir/filibuvir-style mechanisms.
  - NI / nucleoside-like active-site actives: 43.
  - Status: practical powered active-site scope for HCV NS5B after assay-level review.

The HCV NS5B NI active-site file is the practical Phase 9 route to statistical power
and the second target. The raw CHEMBL5375 file is a stratification input, not a pooled
benchmark set. The TiLV PB1 benchmark remains active-limited unless curated literature
adds substantially more PB1/RdRp actives.

## 9B: HCV NS5B Stratification Result

The HCV NS5B benchmark must not pool all CHEMBL5375 actives into one pocket. The raw
target export contains:

- direct active-site NI evidence: too small in CHEMBL5375 alone,
- allosteric NNI series: separate allosteric-pocket benchmark only,
- indirect HCV replicon/interferon assay rows: excluded from pocket benchmark,
- approved HCV drugs with off-target mechanisms such as NS3/NS4A or NS5A: excluded.

Headline scope selected for the second target:

- Class: NI / nucleoside-like inhibitors.
- Pocket: NS5B catalytic active site.
- Active count: 43.
- Decoys: generated in 9B-3 below (133 presumed decoys, demonstration-grade — see that
  section before running any active-site benchmark).

## 9B-3: HCV NS5B NI Active-Site Decoys

Property-matched, scaffold-distinct presumed decoys were generated for the 43 NI
active-site actives.

Builder and artifacts:

- `scripts/fetch_phase9_decoys.py` (injectable ChEMBL seam, salts stripped to largest
  fragment, hermetic tests in `tests/test_phase9_decoys.py`).
- `vta/data/decoys_cache/hcv_ns5b_ni_matched.smi`
- `outputs/phase9/hcv_ns5b_ni_decoy_provenance.json`
- `outputs/phase9/hcv_ns5b_ni_matching_quality.json`
- Registered as decoy set `hcv_ns5b_ni` in `vta/eval/decoys.py`.

Matching method: ChEMBL MW-window query, then local match on MW, logP, HBA, HBD,
rotatable bonds, and formal charge; scaffold-distinct from every active; deduplicated
and excluding active ChEMBL ids.

Result:

- Actives used: 43.
- Decoys generated: 133.
- Decoys per active: mean 3.09, min 0, max 27.
- Actives with at least one matched decoy: 19 / 43.
- Actives with zero matched decoys: 24 / 43.
- Network failures during fetch: 0 (with backoff/retry on the seam).
- **Grade: demonstration-grade, not publication-grade** (threshold is 30 decoys per
  active = 1290 total).

Honest interpretation:

The NI active-site actives are nucleotide triphosphates and phosphoramidate prodrugs.
The triphosphate forms have extreme physicochemistry (logP near -2, many H-bond
donors/acceptors), and ChEMBL contains essentially no scaffold-distinct molecules that
match those windows — hence 24 actives recover zero decoys. The 133 decoys cluster on
the 19 more drug-like (largely prodrug) actives, so the aggregate decoy logP
distribution (median ~2.3) skews more lipophilic than the actives (median ~-0.8). This
is a real limitation of physicochemical decoy matching for nucleotide active forms, not
a fetch error.

Consequence for the benchmark: the HCV NS5B NI active-site decoy set is currently an
underpowered demonstration set, in the same honesty class as TiLV PB1. A powered
active-site benchmark would require either a much larger property-matched decoy pool
(e.g. DeepCoy/DEKOIS generation tuned to nucleotide chemistry) or restricting the active
set to the matchable prodrug subset and saying so explicitly. The HCV benchmark (9B-4)
should not be run as a publication claim on this decoy set; it may be run only as a
labelled demonstration.

## 9D: Mechanistic Input Caveat — and the executed parent-vs-active-form test

The TiLV PB1 benchmark still docks parent compounds for remdesivir, sofosbuvir, and
molnupiravir because curated triphosphate active-form structures are not committed; that
specific re-dock remains a **labelled skip** (we will not fabricate the triphosphate SMILES).

But the underlying question — *how much of poor nucleotide recovery is a prodrug-input
artifact vs. a genuine scoring limitation?* — was answered directly on a set that needs no
fabrication. The HCV NS5B NI active-site set contains **both** real parent/prodrug forms and
ChEMBL-curated triphosphate active forms, so both groups were docked into the real NS5B
catalytic site (**2XI3 chain A, with the experimental catalytic Mg²⁺ retained**, box centred
on the GDD aspartate-220/318/319 centroid) with AutoDock Vina
(`scripts/run_phase9d_active_form.py`, `outputs/phase9/phase9d_active_form_comparison.json`):

| Form | n docked | median ΔG | best ΔG | worst ΔG |
|---|---|---|---|---|
| Triphosphate active form | 26 | −7.57 | −8.81 | −6.93 |
| Parent / prodrug | 17 | −7.47 | −7.90 | −6.71 |

**Median ΔG difference (triphosphate − parent) = −0.10 kcal/mol** — well inside Vina's
~1–2 kcal/mol error. The matched pair makes the point sharply: **sofosbuvir (parent) docked
at −7.90, *stronger* than its own triphosphate GS-461203 at −7.35.**

Honest conclusion: at the rigid-Vina level, the active triphosphate form is **not** scored
meaningfully better than the parent prodrug — Vina does not capture the triphosphate's
catalytic-metal coordination advantage. So the poor parent recovery is **not primarily a
"wrong species docked" artifact**; it reflects a genuine scoring limitation (Vina cannot
reward the mechanistic feature that defines these inhibitors). This is consistent with, and
mechanistically explains, the Phase 9B-3 finding that the nucleotide class resists
matched-decoy validation. (ΔG here is a within-set score comparison, not enrichment — HCV
NS5B NI is un-benchmarkable by matched decoys.)

## 9E: Frozen validation artifact

The matured, CI-bearing benchmark is frozen as a single hash-stamped artifact,
`outputs/phase9/locked_benchmark.json` (`scripts/freeze_benchmark.py`): the three targets
with their metrics + 95% CIs, the **exact Mpro compound splits** (active/inactive ChEMBL/CID
lists + selection seed), decoy/inactive provenance, the 9D active-form comparison, and the
gate decision, with a content hash so drift is detectable. Policy is immutable: a changed
benchmark gets a new dated artifact, not an edit. The gate
(`outputs/phase9/gate_decision_mpro.md`) stays **DO NOT PROMOTE** — even a powered benchmark
does not license a ranking change without a leakage-controlled DL improvement at
non-overlapping CIs, and no DL rescore is wired.

## Gate Status

No ranking term is promoted. The benchmark is now CI-bearing, but it is still:

- single-target for scored docking evidence,
- active-limited for TiLV PB1,
- based on presumed rather than experimentally verified decoys,
- not yet active-form/triphosphate corrected.

## Phase 9C: Two-Track Resolution — Nucleotide Meta-Finding + Powered Mpro Target

Phase 9B-3 forced a fork. The HCV NS5B NI set cannot be validated by matched-decoy
enrichment (above), so two tracks run in parallel: keep the nucleotide set as a labelled
mechanistic case **and** write the impossibility itself up as a result (Track A), and move
the powered enrichment number to a drug-like target with real measured inactives (Track B).

### Meta-finding (first-class result): matched-decoy validation cannot assess nucleotide-analog antivirals

Matched-decoy validation **with a purchasable-library, scaffold-distinct, property-matched
decoy recipe (DUD-E / DEKOIS style) fails for nucleotide-analog antivirals.** Their
mechanistically-active species are nucleoside **triphosphates**, whose physicochemistry
(formal charge near the phosphates, logP ≈ −2, 6–8 H-bond donors, 9–16 acceptors) has almost
no scaffold-distinct, property-matched neighbours reachable by that recipe. The evidence is
direct: under this recipe, for the 43 HCV NS5B NI active-site actives **24/43 recovered zero**
matched decoys and the whole set averaged only 3.09 decoys/active against a 30/active
threshold (`outputs/phase9/hcv_ns5b_ni_matching_quality.json`) — a property of the recipe ×
chemical space, not a fetch failure. This is tied to the recipe, not asserted as impossible in
general: **property-unmatched / charge-extrema (DUDE-Z), generative (DeepCoy), or curated real
inactive nucleotide decoys are the untested alternatives** (WI-4/WI-5). Whether any of those
yields above-chance enrichment on this class is an open, honest question pending an extended
NS5B dock; DeepCoy is not installed in this environment.

The corollary is a concrete warning about the field: a pipeline that appears to "validate"
a nucleotide antiviral by retrospective enrichment is almost always **docking the parent
prodrug**, not the active triphosphate — exactly the caveat already flagged for TiLV PB1,
where remdesivir/sofosbuvir/molnupiravir are docked as parents because curated triphosphate
forms are not committed (§9D). Nucleotide antivirals must be validated through the
prospective experimental loop (Phase 11), not retrospective matched-decoy enrichment.

### Track B: SARS-CoV-2 Mpro (non-covalent) — the powered target

Mpro (3CLpro) was chosen for the powered number: a single well-defined active site, a
different enzyme class from RdRp (supporting a real generalization claim), and — decisively —
COVID Moonshot provides **experimentally measured actives AND inactives**, removing the
presumed-decoy weakness entirely.

Assembled and committed (no docking required):

- Dataset (`vta/data/mpro/mpro_dataset.json`): 2848 measured actives + 930 measured
  inactives from Moonshot (`covid_submissions_all_info.csv`, primary) + ChEMBL
  `CHEMBL4523582` (top-up), each with measured activity value and source.
- Binding-mode stratification (`vta/data/mpro/mpro_stratified.json`): **1945 non-covalent
  actives** (770 Moonshot) for the headline, **903 covalent actives held out** (Vina cannot
  score covalent binders — same binding-mode discipline that separated HCV NI from NNI),
  via the covalent-warhead classifier `vta/chem/warheads.py`.
- Active site wired experimentally: 7L11 chain A, box at the bound non-covalent inhibitor
  XF1 centroid (3.1 Å from Cys145 Sγ, confirming non-covalent engagement) —
  `EXPERIMENTAL_PDB`/`EXPERIMENTAL_ACTIVE_SITE` in the structure/pocket nodes.
- Benchmark harness `scripts/run_mpro_benchmark.py`: ranks by Vina affinity, reports every
  metric as median + 95% bootstrap CI plus a multi-draw control-spread.

**Result (real AutoDock Vina, 50 non-covalent actives / 50 measured Moonshot inactives):**

| Metric | Point | Bootstrap median | 95% CI |
|---|---|---|---|
| BEDROC(α=20) | 0.677 | 0.682 | [0.357, 0.893] |
| logAUC | 0.199 | 0.206 | [0.137, 0.301] |
| ROC-AUC | 0.580 | 0.580 | [0.467, 0.690] |
| EF1% | 2.0 | 1.96 | [0.0, 2.44] |

Receptor prep was unblocked with an **OpenBabel fallback** (`vta/nodes/docking._prep_receptor`
tries Meeko first — so the TiLV path is untouched — then falls back to a rigid-receptor
OpenBabel PDBQT, which succeeds where Meeko 0.7.1 deterministically fails on the Mpro chain).
The docking is real (`dG_provenance.tool == "AutoDock Vina"`); no mock score is used.

The headline is **methodological, and honest about its two layers**:

1. *Benchmark design is publication-grade.* With experimentally measured inactives, binding-
   mode stratification, real Vina, and bootstrap + multi-draw CIs, the intervals are
   **informative** — BEDROC [0.357, 0.893], ROC-AUC [0.467, 0.690] — not the [0,1] of TiLV.
   The pipeline *can* produce a powered, CI-bearing enrichment number; that was the point.
2. *The docking signal itself is modest.* ROC-AUC ≈ 0.58 with a CI that **crosses 0.5**, so
   rigid-receptor Vina only weakly separates non-covalent Mpro actives from measured
   inactives; BEDROC 0.68 (CI lower bound 0.36) indicates some early enrichment. EF1% is
   ceiling-limited at the 50/50 ratio (max 2.0); a larger inactive:active ratio would sharpen
   the early-enrichment metrics. This modest number is reported as-is, not inflated.

### Multi-target table + gate

`outputs/phase9/multitarget_ci_table.md` places all three targets side by side with CIs and
an explicit grade/status: TiLV PB1 (underpowered demonstration, CIs ≈ [0,1]), HCV NS5B NI
(un-benchmarkable meta-finding), SARS-CoV-2 Mpro non-covalent (powered benchmark, informative
CIs, modest docking signal). The gate (`outputs/phase9/gate_decision_mpro.md`) is **DO NOT
PROMOTE**: even with a powered benchmark, promotion requires a DL/consensus signal that beats
the Vina baseline with non-overlapping CIs, and `consensus_node`'s `cnn_affinity`/
`boltzina_score` are not populated (no DL rescore wired). Consensus/DL outputs stay
annotation-only and ranking weights are unchanged.
