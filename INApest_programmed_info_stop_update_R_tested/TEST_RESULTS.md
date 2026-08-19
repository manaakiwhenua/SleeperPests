# Programmed information-stop test results

## Summary

All eight public INApest / INApestMeta implementations passed the programmed-stop model tests on a 30 x 30 (900-node) landscape.

### Scalar programmed stop

With all nodes initially informed, `ManageProb = 1`, and a last-known-presence event at timestep 1:

- `InfoPersistenceSteps = 2`: management profile = `900 -> 900 -> 900 -> 0`.
- `InfoPersistenceSteps = 0`: management profile = `900 -> 0 -> 0 -> 0`.
- Multiple-land-use models have two managed land uses per node, so the corresponding profiles are `1800 -> 1800 -> 1800 -> 0` and `1800 -> 0 -> 0 -> 0`.

### Spatial and temporal parameterisation

- Node vector: half nodes `2`, half nodes `0` -> `900 -> 450 -> 450 -> 0`.
- Node x timestep matrix: policy change for half the nodes at timestep 2 -> `900 -> 900 -> 450 -> 0`.
- Multiple-land-use equivalents are doubled.

### Priority over stochastic decay

With `InfoPersistenceSteps = 2` and `InfoRetentionProb = 0`, all functions retained the programmed-stop profile (`900 -> 900 -> 900 -> 0`, or doubled for MLU) and emitted exactly one priority warning.

A mixed node-vector test with programmed persistence on half the nodes and stochastic retention 0 on the other half produced `900 -> 450 -> 450 -> 0`.

### Local-evidence reset

For Meta models, a management kill at timestep 2 reset the clock and extended management through timestep 3 (`900 -> 900 -> 900 -> 0`). Without a management kill, the corresponding profile was `900 -> 900 -> 0 -> 0`.

For occupancy INApest models, continuing informed extant presence resets the clock each timestep, so management continued through all four timesteps in the comparison test.

### Information from other nodes

A dedicated donor-recipient test used `InfoPersistenceSteps = 0` at the recipient. With no SEAM link, recipient management was `1 -> 0 -> 0 -> 0` (or `2 -> 0 -> 0 -> 0` for MLU). With a deterministic donor-to-recipient SEAM link, incoming information refreshed the recipient after each local stop and management was `1 -> 1 -> 1 -> 1` (or `2 -> 2 -> 2 -> 2` for MLU) in every implementation.

### Backward compatibility and regressions

- Every implementation was exactly backward compatible with the previous info-decay version when `InfoPersistenceSteps` was omitted.
- All eight stochastic `InfoRetentionProb` regression suites passed.
- Serial temporal SDD/LDD regression: INApest, INApestMetaMultipleLandUse and INApestMetaTransitionMatrix all passed 2D compatibility, dynamic SDD and dynamic LDD tests.
- Parallel-body temporal SDD/LDD regression: all five parallel implementations passed through the sequential test shims.
- All eight final production files parse successfully in R 4.6.0.
- Source-invariant checks pass for all eight files.

## Runtime caveat

Some WebR Node wrapper invocations report an `UnwindProtectException` after the R script has completed. In those cases every model PASS marker and result file had already been written. The same tests split into smaller direct runs return normally where the wrapper permits. No R assertion failed.
