# INApest corrected production source release - 2026-08-17

This directory contains the eight current INApest / INApestMeta production functions, reconciled against the GitHub `main` branch on 2026-08-17 and then updated with the cumulative corrections developed and tested during the analytical-companion work.

## Production files

- `INApest.R`
- `INApestParallel.R`
- `INApestParallelInAScene.R`
- `INApestMetaParallel.r`
- `INApestMetaMultipleLandUse.r`
- `INApestMetaParallelMultipleLandUse.r`
- `INApestMetaTransitionMatrix.r`
- `INApestMetaTransitionMatrixParallel.r`

The analytical companion is intentionally not included here; these are the stochastic production functions intended to replace/update the repository copies.

## Current-main functionality retained

The release was built on top of the current-main source rather than replacing it with an older snapshot. In particular, it retains the existing time-varying 3-D SDD/LDD support and the newer transition-matrix/local-dynamics arguments and logic already present on `main`.

## Cumulative corrections/features included

1. **Time-varying habitat initialization in binary INApest**
   - When `EnvEstabProb` is a node x timestep matrix, `BPAM` is initialized using `EnvEstabProb[,1]` before the timestep loop.
   - Applied to `INApest`, `INApestParallel`, and `INApestParallelINAscene`.

2. **Information retention/decay**
   - `InfoRetentionProb` can be a scalar, node vector, or node x timestep matrix.
   - Information can decay after management/spread and can subsequently be refreshed by information transfer or direct detection.

3. **Programmed information stopping**
   - New backward-compatible argument `InfoPersistenceSteps = NA`.
   - Accepts a scalar, node vector, or node x timestep matrix of non-negative whole numbers; `NA` leaves that node/timestep on the `InfoRetentionProb` pathway.
   - Programmed stopping has priority where a finite persistence value is supplied.
   - Direct detection refreshes the local-evidence clock. Binary INApest also refreshes the clock from informed extant infestation; Meta population models use realized management-caused mortality as local evidence without changing the existing mortality draw.

4. **Serial transition-matrix timestep indexing**
   - Current-timestep detection uses `NodeDetectionProb[,,timestep]`, eliminating the stale loop-index use.
   - The current-main mortality calculation already uses the correct `timestep` index and is retained.

5. **Multiple-land-use SEAM dimensional correction**
   - Serial MLU information transfer collapses node x land-use detection to node level with `apply(Detected,1,max)`, matching the intended node x node SEAM calculation.

6. **Multiple-land-use LDD population-share ordering**
   - In the serial MLU local-dynamics helper, source population shares (`Pn`) are applied to LDD propagule production before the multinomial dispersal draw.

7. **`Nperm = 1` dimension preservation**
   - Result extraction and invasion-probability summaries preserve the permutation dimension with `drop = FALSE` and singleton-dimension repairs where required.
   - Applied across serial and parallel implementations.

8. **INAscene information refresh ordering**
   - When information loss is active, SEAM refresh is applied after stop/decay so incoming information can restore a node whose prior information has just expired.
   - The default no-information-loss path remains on the existing INAscene route.

## Validation performed for this reconciled release

- Fresh copies of all eight files were retrieved from GitHub `main` on 2026-08-17 before reconciliation.
- All eight corrected source files parse successfully in WebR / R 4.6.0.
- The three serial production files (`INApest.R`, `INApestMetaMultipleLandUse.r`, and `INApestMetaTransitionMatrix.r`) also source successfully in WebR / R 4.6.0.
- Structural checks confirm the cumulative fixes above are present in the intended functions.
- The underlying correction sets were previously regression-tested across the eight public APIs; parallel bodies were tested with sequential shims where WebR lacked the usual PSOCK/doParallel/INA stack.

This validation does not replace a full installed-package/concurrent-worker integration test in the normal desktop/server R environment after the files are committed.

## Audit files

- `corrected_vs_current_main.patch` - unified diff from the GitHub `main` files retrieved on 2026-08-17 to this corrected release.
- `SOURCE_SHA256SUMS.txt` - hashes of the eight corrected production files.
- `CURRENT_MAIN_BASELINE_SHA256SUMS.txt` - hashes of the eight GitHub `main` files used as the reconciliation base.
- `VALIDATION.txt` - concise validation record.
