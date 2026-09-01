# Layer 3 harness correction — 2026-09-01

The first native-R Layer-3 run reported 14 failures with observed PoA values equal to their priors and `BackgroundSSe = NA` in the spatial mapping check.

Cause: `core_one_round()` returned `z$PoASummary[1,]`. In the PoA output contract, row 1 is the Round-0 prior summary. The actual surveillance update is the row with `Round == 1`.

Correction: the harness now selects `z$PoASummary[z$PoASummary$Round == 1L,]` and stops unless exactly one such row exists.

This is a validation-harness correction only. `source/INApestPoA.R` is unchanged.
