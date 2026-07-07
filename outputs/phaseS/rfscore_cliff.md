# Phase S — in-house RF-Score learned rescorer on activity cliffs

RF-Score v1 contact-count features + RandomForest; scaffold-clustered leave-out CV

- compounds featurized: 1109 across 425 scaffold groups
- cliff pairs scored (RF+Vina+2D all available): 1193

## Cliff-pair ranking accuracy (median [95% CI]); chance = 0.50

| Method | accuracy |
|---|---|
| **RF-Score (learned)** | 0.3923 [0.3646, 0.4199] |
| 2D-kNN QSAR | 0.6655 [0.6387, 0.6924] |
| Vina (−ΔG) | 0.4183 [0.3898, 0.4459] |

Paired RF−2D: Δ -0.2733 [-0.3127, -0.2339]; RF−Vina: Δ -0.026 [-0.0629, 0.0109].

**The learned RF-Score rescorer does NOT beat both the 2D-kNN QSAR and chance on activity cliffs (scaffold-CV, out-of-fold) — learned rescoring did not demonstrate cliff-ranking skill over the trivial baseline.**

## Caveats
- TARGET-SPECIFIC: RF-Score is trained on Moonshot Mpro itself (out-of-fold under scaffold-clustered CV — no analog leakage), so this is NOT a general-transfer claim. A PDBbind-pretrained scorer (GNINA/RTMScore) on a native box answers that separately.
- RF-Score v1 (coarse contact counts) is weaker than modern CNN/GNN scorers; a null here is suggestive, not proof that GNINA would also fail on these cliffs.
- Poses are the committed Vina docks; the same 1,193-pair set and 2D-kNN baseline as the powered cliff benchmark, through the identical paired gate.
