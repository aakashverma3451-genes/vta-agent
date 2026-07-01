# Phase 9 Multi-Target Enrichment Summary

Every metric is reported as median + 95% bootstrap CI, never a bare point. Each
target carries an explicit grade/status. Two of the three results are negative or
blocked — and the negatives are the scientific contribution, not a failure to hide.

| Target | Class / pocket | Actives | Control set | BEDROC(α20) med [CI] | EF1% med [CI] | logAUC med [CI] | ROC-AUC med [CI] | Grade / status |
|---|---|---|---|---|---|---|---|---|
| TiLV PB1 | RdRp NTP site (8PSO:B) | 4 | ChEMBL property-matched presumed decoys | 0.4062 [0.0009, 1.0] | 9.0 [0.0, 39.825] | 0.3877 [0.1327, 1.0] | 0.7921 [0.5976, 1.0] | underpowered demonstration |
| HCV NS5B NI | catalytic active site (triphosphate) | 43 | property-matched decoys IMPOSSIBLE (24/43 actives recovered 0 decoys; 3.093/active) | n/a | n/a | n/a | n/a | **un-benchmarkable** (meta-finding) |
| SARS-CoV-2 Mpro (non-covalent) | 3CLpro active site (7L11:A) | 50 | measured Moonshot inactives | 0.6823 [0.3574, 0.8933] | 1.9608 [0.0, 2.439] | 0.2056 [0.137, 0.3012] | 0.5795 [0.4669, 0.6903] | powered benchmark — real measured inactives, informative CIs (not [0,1]); docking signal MODEST — ROC-AUC CI [0.4669, 0.6903] includes 0.5, EF1% ceiling-limited at 50/50 |

## Reading the table

- **TiLV PB1**: a real docking run, but four actives make every CI span almost the
  whole [0,1] range — a labelled demonstration, not an enrichment claim.
- **HCV NS5B NI**: not assessable by matched-decoy enrichment at all (see meta-finding).
- **SARS-CoV-2 Mpro**: the intended powered target with real measured inactives; the
  dataset and active-site wiring are complete, but the headline number is blocked on a
  Meeko receptor-prep bug. No mock or fabricated number is reported in its place.
