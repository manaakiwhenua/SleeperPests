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
5. The safer initial-invasion setup used by the transition implementation was migrated to sibling functions where the same initialisation pattern existed. This also avoids failures when optional initial-invasion inputs are left at their defaults and handles matrix `InvasionRisk` consistently using the first timestep.
6. In serial `INApestMetaTransitionMatrix.r`, management mortality now indexes `NodeMortalityProb[, s, timestep]`, matching the parallel implementation, rather than the stale pre-generation loop index `t`.

Transition-specific biology (stage weighting, blocked-transition mortality, density-dependent dispersal, seedbank reporting, etc.) was not migrated into models where it does not apply.

## Tests run here

The following passed in the build environment:

```text
PASS: 30x30 temporal-connectivity regression tests
nodes=900, timesteps=4
static combined non-zero links=4380
barrier timestep non-zero links=4320
changed-LDD timestep non-zero links=4380
verified: 2D == repeated 3D, mixed 2D/3D, SDD temporal change, LDD temporal change, validation, establishment weighting

PASS: source invariants and delimiter balance for 8 modified R files
PASS: R model-test script delimiter balance
PASS: patch applies cleanly and reproduces all 8 modified files
```

The build environment did **not** contain GNU R / `Rscript`, so the supplied model-level R regression harness could not be executed here. It is included as `tests/test_models_temporal_connectivity_30x30.R` for execution in an R environment with the model dependencies installed.

## Apply the patch

From inside the repository directory containing the eight R files:

```bash
patch -p1 < /path/to/INApest_temporal_connectivity.patch
```

The patch was tested by applying it to clean copies of all eight starting files and byte-comparing the result with the supplied `modified/` files.


The two environment-independent checks can also be rerun from the bundle root:

```bash
python3 tests/test_temporal_connectivity_30x30.py
python3 tests/test_source_invariants.py
```

The connectivity test requires NumPy.

## Run the model-level R test

From the root of this bundle:

```bash
Rscript tests/test_models_temporal_connectivity_30x30.R
```

The script uses a 30 x 30 (900-node) landscape and checks static 2D versus repeated 3D connectivity, dynamic SDD, and dynamic LDD across the public implementations. Parallel tests require `abind` and `doParallel`; the INAScene test also requires `INA` and is skipped if unavailable.
