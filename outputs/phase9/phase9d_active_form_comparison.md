# Phase 9D: Parent prodrug vs active triphosphate (HCV NS5B catalytic site)

Structure: 2XI3 chain A (GTP + catalytic Mg2+). Box center [9.451, 7.204, 11.876] (centroid of catalytic GDD aspartates 220/318/319 CA). Metal: Mg2+ retained from experimental 2XI3 HETATM (curated, not fabricated).
Scoring: AutoDock Vina ΔG (kcal/mol); ligands docked AS GIVEN (no species resolution).

## ΔG by form (more negative = stronger predicted binding)

| Form | n docked | median ΔG | best ΔG | worst ΔG |
|---|---|---|---|---|
| Triphosphate active form | 26 | -7.573 | -8.81 | -6.925 |
| Parent / prodrug | 17 | -7.469 | -7.897 | -6.708 |

Median ΔG (triphosphate − parent): **-0.104** kcal/mol (negative = triphosphate binds stronger).

## Matched pair

- Sofosbuvir (parent): ΔG -7.897
- GS-461203 (sofosbuvir triphosphate): ΔG -7.346

## Interpretation

More negative ΔG = stronger predicted binding. If the triphosphate group docks markedly stronger than the parent/prodrug group at the metal-containing catalytic site, the poor parent recovery is substantially a prodrug-input artifact, not purely a scoring limitation. ΔG is not enrichment: HCV NS5B NI is un-benchmarkable by matched decoys (Phase 9B-3 meta-finding), so this is a within-set score comparison, reported as such.

## TiLV active-form

- labelled skip: TiLV PB1 controls have no curated triphosphate active-form SMILES (vta/chem/species.py active_form_smiles=None); not fabricated.
