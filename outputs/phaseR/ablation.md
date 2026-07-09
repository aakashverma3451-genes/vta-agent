# Phase R / R6 — component ablation (decision-level correctness)

Held-out targets spanning classes; each level toggles one reasoning layer on. Scored on the emitted CLAIM vs a curated ground-truth disposition — NOT enrichment.

## Decision accuracy / false-confidence per level

| Level | decision accuracy | false-confidence rate | OOD-flag recall |
|---|---|---|---|
| L0_raw | 0.0 | 1.0 | 0.0 |
| L1_dossier | 0.0 | 1.0 | 0.0 |
| L2_triage | 0.6 | 0.4 | 1.0 |
| L3_verification | 1.0 | 0.0 | 1.0 |
| L4_full_agent | 1.0 | 0.0 | 1.0 |

## Per-target emitted claim (L0 raw → L4 full agent)

| Target | class | ground truth | L0_raw | L1_dossier | L2_triage | L3_verification | L4_full_agent |
|---|---|---|---|---|---|---|---|
| MPRO | protease (in-domain, powered) | dock_but_downgraded | dock_and_claim_enrichment | dock_and_claim_enrichment | dock_and_claim_enrichment | dock_but_downgraded | dock_but_downgraded |
| NS5B | polymerase / nucleotide (out-of-domain) | annotate_only | dock_and_claim_enrichment | dock_and_claim_enrichment | annotate_only | annotate_only | annotate_only |
| PB1 | polymerase (underpowered demonstration) | dock_but_downgraded | dock_and_claim_enrichment | dock_and_claim_enrichment | dock_and_claim_enrichment | dock_but_downgraded | dock_but_downgraded |
| GPX | predicted, very low binding-site pLDDT | refuse | dock_and_claim_enrichment | dock_and_claim_enrichment | refuse | refuse | refuse |
| GPY | predicted, borderline binding-site pLDDT | defer | dock_and_claim_enrichment | dock_and_claim_enrichment | defer | defer | defer |

## Which layer fixed which target

- **MPRO** — corrected at **L3_verification**
- **NS5B** — corrected at **L2_triage**
- **PB1** — corrected at **L3_verification**
- **GPX** — corrected at **L2_triage**
- **GPY** — corrected at **L2_triage**

Run-to-run consistency: 1.0

**Routing-type-dependent, as HemaGuide predicts: the raw pipeline (L0) makes an unwarranted structure-based enrichment claim on every target; the ROUTER fixes the out-of-domain targets (annotate/defer/refuse) and the VERIFICATION gate fixes the in-domain-but-loses-to-2D targets (downgrade). No single layer is sufficient — the architecture, not the scorer, is what makes the decisions honest.**
