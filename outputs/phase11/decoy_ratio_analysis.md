# Phase 11 WI-4 — Decoy/inactive ratio analysis (Mpro)

Available: 1945 non-covalent actives, 930 measured inactives.
Current benchmark: 50/50 (1:1) — ratio-degenerate (EF1% ceiling 2.0), kept labelled.

## Achievable ≥30:1 WITHOUT fabrication
- Use **30 actives × 900 measured Moonshot inactives = 30:1** (900 ≤ 930 available). All real, no DeepCoy, no fabrication.
- Cost: dock ~850 additional real inactives with production Vina (multi-hour run) — NOT executed this pass.

## Synthetic decoys
- DeepCoy: **not installed** → generative decoys unavailable here.
- DUDE-Z property-unmatched decoys: generatable from ChEMBL, but still require docking (same multi-hour cost). Every synthetic decoy would log tool+version+seed; none generated → none added (no unlabelled padding).

**Status:** analysis only. 30:1 is reachable with real inactives via an extended dock; current headline remains 1:1, labelled ratio-degenerate.
