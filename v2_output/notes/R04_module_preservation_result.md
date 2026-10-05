# Module preservation, both directions (R04)

Date: 2026-10-01. Script: `Code/Revision_v2/R04_module_preservation.R`.
Tables: `Table_S19` (all modules), `Table_S20` / `S20b` (designated modules).
Raw: `intermediate/R04_preservation.rds`, `R04_preservation_raw.rds`.
200 permutations, both directions, common gene universe `U`.

## Designated translation modules

| reference | module | n genes | tested in | Zsummary | medianRank | BH FDR | class |
|---|---|---|---|---|---|---|---|
| HF | purple | 153 | HCC | 11.03 | 3 | 4.2e-41 | strong (>10) |
| HCC | magenta | 146 | HF | 7.09 | 5 | 8.0e-15 | weak-moderate (2-10) |

HF→HCC: transmission module ranks **2 of 23** named modules.
HCC→HF: ranks **6 of 20**.

Direction-level totals (gold pseudo-module and grey unassigned bin dropped):

| direction | named modules | strong >10 | weak-mod 2-10 | none <2 |
|---|---|---|---|---|
| HF → HCC | 23 | 3 | 9 | 11 |
| HCC → HF | 20 | 2 | 12 | 6 |

## The plan's contingency is NOT triggered, but the framing changes

The plan flagged: *"跨平台 modulePreservation 可能显示不保留 —— 预案：先跑同技术配对；
若连配对都不保留，则如实报告并重新拟定标题."* Neither designated module falls below
the Zsummary = 2 floor, so the modules **are** preserved and the contingency in its
literal form does not fire.

**But the result refutes v1's claim, and this is the part that must reach the
manuscript.** v1 used Fisher's exact test on gene overlap to argue that the limited
overlap indicated disease-specific network remodelling. The formal test says the
opposite for the translation module: it is preserved. So:

- The claim "the translation module is disease-specific" **must be removed**.
- The Title survives only in a narrower sense: the shared object is the module
  **architecture**; what differs between diseases is the **direction of its
  association with disease state** (the mirror).
- v2 must present preservation **before** the mirror, so a reader cannot infer a
  structural difference that does not exist.

This is a net gain for the paper's coherence — "same architecture, opposite disease
association" is a sharper claim than "different modules" — but it is a claim v1 did
not make and would not have survived making.

## Asymmetry must be reported, not averaged

HF purple is **strongly** preserved in HCC; HCC magenta is only **weakly-to-moderately**
preserved in HF. Roughly half the modules in each direction show no evidence of
preservation at all. Reporting only the favourable direction would be exactly the
selective reporting R2 and R3 objected to.

## Caveat: maxModuleSize = 1000 subsamples the large modules

`WGCNA::modulePreservation` defaults to `maxModuleSize = 1000`, documented as:
*"maximum module size used for calculations. Modules larger than maxModuleSize will be
reduced by randomly sampling maxModuleSize genes."*

Verified against `Table_S4_module_inventory.csv`:

| module | true size (S4) | moduleSize in S19 |
|---|---|---|
| HF turquoise | 2132 | 1000 |
| HF blue | 1525 | 1000 |
| HF brown | 1516 | 1000 |
| HCC turquoise | 3771 | 1000 |
| HCC blue | 1948 | 1000 |
| HCC brown | 389 | 389 |

So the Zsummary values for turquoise, blue and brown in each direction are computed on
a **random 1,000-gene subsample**, not the whole module. Both designated translation
modules (153, 146 genes) are well under the cap and are unaffected. No conclusion
depends on the three subsampled values. Must be disclosed in the Supporting
Information legend.

## Reruns note

The production artefacts were written at 20:45 by a run started 19:37. An earlier
smoke run (2 permutations) had written to the same production path because the final
`saveRDS` bypassed `tn()`; that was fixed before this run and the production files now
carry the real 200-permutation results (`R04_preservation.rds` mtime 20:45:45).
`R04_preservation_raw_SMOKE.rds` and `SMOKE_Table_S19/S20/S20b` are leftovers from that
smoke run and must be deleted before packaging so they cannot be shipped by accident.
