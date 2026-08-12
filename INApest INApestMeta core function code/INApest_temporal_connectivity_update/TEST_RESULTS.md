# Test results: temporal SDD/LDD update

Tests were executed on 12 August 2026 using the uploaded WebR 0.6.0 runtime, which ran R 4.6.0. The landscape was 30 x 30 (900 nodes), with three simulation timesteps and two realisations for model-level tests.

## Production source status

No production R source file required a further change as a result of the model-level tests. The `modified/` files are byte-identical to those in the original update bundle supplied before WebR was uploaded.

## Direct model-level tests in R

The following serial implementations were executed directly from the production source files:

| Function | 2D = repeated 3D | dynamic SDD | dynamic LDD | Result |
|---|---:|---:|---:|---|
| `INApest` | PASS | PASS | PASS | PASS |
| `INApestMetaMultipleLandUse` | PASS | PASS | PASS | PASS |
| `INApestMetaTransitionMatrix` | PASS | PASS | PASS | PASS |

Additional direct R edge tests passed for these serial implementations:

- mixed 3D SDD + 2D LDD;
- mixed 2D SDD + 3D LDD;
- rejection of an SDD array with the wrong number of timestep slices;
- rejection of an LDD array with the wrong number of timestep slices;
- guarded no-initial-invasion smoke tests for `INApest` and `INApestMetaMultipleLandUse`.

## Parallel implementation body tests

The uploaded WebR runtime does not provide the external R packages `abind`, `doParallel`, or `INA`, and its Node environment cannot create a normal PSOCK cluster without an additional WebSocket server dependency. To test the new simulation logic without modifying production source, a test-only sequential backend was used.

The original model bodies were loaded from the production source with only parallel/package plumbing replaced in the test environment. For `INApestParallelINAscene`, a small `INAscene` test double was used to expose the supplied biophysical adjacency matrix and return the nested result structure expected by the wrapper.

| Function | 2D = repeated 3D | dynamic SDD | dynamic LDD | Result |
|---|---:|---:|---:|---|
| `INApestParallel` | PASS | PASS | PASS | PASS |
| `INApestParallelINAscene` | PASS | PASS | PASS | PASS |
| `INApestMetaParallel` | PASS | PASS | PASS | PASS |
| `INApestMetaParallelMultipleLandUse` | PASS | PASS | PASS | PASS |
| `INApestMetaTransitionMatrixParallel` | PASS | PASS | PASS | PASS |

These tests validate the temporal-connectivity branches and simulation bodies, but they are not a substitute for a final concurrency/integration run with the real `doParallel`/PSOCK/`INA` dependencies in a normal R installation.

## Source and environment-independent checks

All eight production R files parse successfully in R 4.6.0.

The source-invariant test passed for all eight modified files, including delimiter balance and the expected temporal SDD/LDD branches.

The independent 900-node connectivity regression also passed, verifying:

- legacy 2D connectivity equals repeated 3D connectivity;
- mixed 2D/3D connectivity;
- timestep-specific SDD change;
- timestep-specific LDD change;
- dimensional validation;
- establishment weighting and probability bounds.

## Existing Nperm = 1 behaviour

The first direct run used `Nperm = 1` and exposed an existing result-summary issue: expressions such as `InvasionResults[, timestep, ]` drop to a vector, after which `rowSums()` errors. This is unrelated to the temporal-connectivity changes. The model-level regression therefore uses `Nperm = 2`. The production functions were not changed for this separate issue in order to keep this update tightly scoped.
