# Phase S — in-house RF-Score rescorer on activity cliffs: DONE (learned rescoring also fails)

**Date:** 2026-07-07 · **Owner:** session-a (opus)

**Why in-house:** every off-the-shelf pretrained rescorer is un-installable on this arm64 stack —
GNINA (10 GB amd64 image + unreliable Docker Hub CDN, 4+ failed pulls), RTMScore (`dgl` has no
working arm64/torch-2.12 wheel), ODDT/RF-Score (needs OpenBabel-2.x `OBElementTable`, removed in
OpenBabel 3.x). So RF-Score's own published method was rebuilt on the working stack (openbabel
read + numpy + scikit-learn) to still ask the question with a *learned* scorer.

## Method
RF-Score v1 contact-count features (protein∈{C,N,O,S} × ligand∈{C,N,O,F,P,S,Cl,Br,I} within 12 Å
= 36 features) from every committed Vina pose → RandomForest → **scaffold-clustered leave-out CV**
(compounds sharing a Bemis–Murcko scaffold never straddle train/test → no analog leakage) →
out-of-fold pIC50 predictions, scored on the SAME 1,193 cliff pairs and 2D-kNN baseline as the
powered cliff benchmark, through the identical paired-bootstrap gate.

## Result (`outputs/phaseS/rfscore_cliff.{json,md}`) — 1,109 poses, 425 scaffold groups

| Method | cliff-ranking accuracy (median [95% CI]); chance = 0.50 |
|---|---|
| **RF-Score (learned)** | **0.392 [0.365, 0.420]** — significantly **BELOW** chance |
| 2D-kNN QSAR | **0.666 [0.639, 0.692]** — significantly above (matches the committed benchmark) |
| Vina (−ΔG) | 0.418 [0.390, 0.446] (matches the committed benchmark) |

- Paired **RF − 2D-kNN = −0.273 [−0.313, −0.234]** — RF loses to the trivial QSAR, decisively.
- Paired **RF − Vina = −0.026 [−0.063, +0.011]** — RF is **statistically indistinguishable from
  raw Vina** (CI includes 0).

## Verdict — learned rescoring does NOT rescue docking on cliffs
The learned RF-Score rescorer is *also* significantly below chance — the same
anti-correlation-with-potency failure as Vina — and performs no better than raw Vina. Both
pose-based scorers (physics and learned) lose to a trivial ligand-based QSAR on activity cliffs.
This strengthens the powered fair-arena negative: it is not merely that *AutoDock Vina* is weak
here — a learned scorer trained on this very target fails too.

## Honest scope (what this does and does NOT settle)
- **TARGET-SPECIFIC**, not general transfer: RF-Score is trained on Moonshot Mpro itself
  (out-of-fold, scaffold-CV — no leakage), so this is not the "does a PDBbind-pretrained scorer
  transfer" question. That still needs GNINA/RTMScore on a native Linux/GPU box.
- **RF-Score v1 is coarse** (contact counts) and weaker than a 3D-CNN. A null here is *suggestive*
  that cliffs are hard for pose-based scorers, but does NOT prove GNINA (finer 3D read) would also
  fail — that remains the one open follow-up, deferred to a native environment.
- No ranking change; the R4 gate keeps docking annotation-grade.
