# INApest temporal connectivity update

Prepared against the current `main` versions in:

`INApest INApestMeta core function code`

Reference implementation for shared recent changes: `INApestMetaTransitionMatrix.r` (starting-copy SHA-256 `fc30ee4f655ad6d7d0371d21c6906898697eb066ce60c632a97ef3430bb6ea8b`).

## Files updated

- `INApest.R`
- `INApestParallel.R`
- `INApestParallelInAScene.R`
- `INApestMetaParallel.r`
- `INApestMetaMultipleLandUse.r`
- `INApestMetaParallelMultipleLandUse.r`
- `INApestMetaTransitionMatrix.r`
- `INApestMetaTransitionMatrixParallel.r`

## Changes

1. `SDDprob` and `LDDprob` may now be supplied either as the existing 2D matrices or as 3D arrays with dimensions `nodes x nodes x Ntimesteps`.
2. For 3D inputs, the simulation selects the current `[, , timestep]` slice immediately before existing dispersal calculations. Existing dispersal kernels remain matrix-based.
3. Mixed inputs are supported: one of `SDDprob`/`LDDprob` can be static (2D) while the other is temporal (3D).
4. Existing 2D behaviour is retained; no unrelated parallel-core scheduling changes are included.
5. The safer initial-invasion setup used by the transition implementation was migrated to sibling functions where the same initialisation pattern existed. This avoids failures when optional initial-invasion inputs are left at their defaults and handles matrix `InvasionRisk` consistently using the first timestep.
6. In serial `INApestMetaTransitionMatrix.r`, management mortality now indexes `NodeMortalityProb[, s, timestep]`, matching the parallel implementation, rather than the stale pre-generation loop index `t`.

Transition-specific biology (stage weighting, blocked-transition mortality, density-dependent dispersal, seedbank reporting, etc.) was not migrated into models where it does not apply.

## Test status

Model-level R tests were run after an R runtime became available. On a 30 x 30 landscape, direct serial tests passed for `INApest`, `INApestMetaMultipleLandUse`, and `INApestMetaTransitionMatrix`, including 2D compatibility, temporal SDD, temporal LDD, mixed 2D/3D inputs, dimensional validation, and guarded initialisation cases.

The five parallel implementation bodies also passed the 2D/repeated-3D, dynamic-SDD, and dynamic-LDD suites using test-only sequential shims because the uploaded WebR environment lacks the normal parallel/package dependencies. The `INApestParallelINAscene` test used an `INAscene` test double. See `TEST_RESULTS.md` for the exact scope and limitations.

All eight production R files parse in R 4.6.0. The source-invariant and independent 900-node temporal-connectivity regression tests also pass.

The model-level tests use `Nperm = 2`. A first run with `Nperm = 1` exposed a pre-existing output-summary dimension-drop issue; it was not changed as part of this update.

## Apply the patch

From inside the repository directory containing the eight R files:

```bash
patch -p1 < /path/to/INApest_temporal_connectivity.patch
```

## Tests

Environment-independent checks:

```bash
python3 tests/test_temporal_connectivity_30x30.py
python3 tests/test_source_invariants.py
```

Direct model-level R test in a normal R environment:

```bash
Rscript tests/test_models_temporal_connectivity_30x30.R
```

Parallel tests in that script require `abind` and `doParallel`; the INAScene test also requires `INA`.

Additional WebR test scripts used during verification are included as `tests/test_serial_edges_30x30.R` and `tests/test_parallel_shim_30x30.R`. The latter is explicitly a test shim and does not modify the production R source.
