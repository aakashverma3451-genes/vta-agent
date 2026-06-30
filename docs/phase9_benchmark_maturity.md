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

## 9D: Mechanistic Input Caveat

The TiLV PB1 benchmark still docks parent compounds for remdesivir, sofosbuvir, and
molnupiravir because curated triphosphate active-form structures are not committed.
Therefore the current 1/4 top-5 recovery should not be treated as the final
mechanistically-correct number.

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

Property-matched-decoy validation (DUD-E / DEKOIS style) **systematically cannot assess
nucleotide-analog antivirals.** Their mechanistically-active species are nucleoside
**triphosphates**, whose physicochemistry (formal charge near the phosphates, logP ≈ −2,
6–8 H-bond donors, 9–16 acceptors) has no scaffold-distinct, property-matched neighbours in
chemical databases. The evidence is direct: for the 43 HCV NS5B NI active-site actives,
**24/43 recovered zero** matched decoys and the whole set averaged only 3.09 decoys/active
against a 30/active publication threshold — not a fetch failure, a structural property of
chemical space (`outputs/phase9/hcv_ns5b_ni_matching_quality.json`).

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

**Status: BLOCKED on receptor prep, no number claimed.** AutoDock Vina is installed and the
TiLV 8PSO receptor preps successfully, but Meeko 0.7.1 fails deterministically on the Mpro
chain (`update_H_positions: Updated 1 H positions but deleted N`) across 11 structures tried
(apo and holo), independent of hydrogens, HETATM, chain extraction, termini, or blunt-ends;
no fallback prep tool (OpenBabel/reduce/ADFR) is installed. The harness therefore writes a
labelled blocked artifact (`outputs/phase9/mpro_noncovalent_benchmark.json`) rather than a
mock/fabricated enrichment. Unblocking needs only a working receptor-prep path; the dataset,
stratification, and wiring are complete.

### Multi-target table + gate

`outputs/phase9/multitarget_ci_table.md` places all three targets side by side with CIs and
an explicit grade/status: TiLV PB1 (underpowered demonstration), HCV NS5B NI
(un-benchmarkable meta-finding), SARS-CoV-2 Mpro non-covalent (powered target, blocked on
receptor prep). The gate (`outputs/phase9/gate_decision_mpro.md`) is **DO NOT PROMOTE**:
with no powered number yet, and no DL/consensus rescore wired to beat the Vina baseline with
non-overlapping CIs, consensus/DL outputs stay annotation-only and ranking weights are
unchanged.
