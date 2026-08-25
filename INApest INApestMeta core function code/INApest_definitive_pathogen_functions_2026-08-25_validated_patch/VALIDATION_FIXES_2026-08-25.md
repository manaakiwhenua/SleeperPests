# Fresh webR regression and validation fixes — 25 August 2026

A fresh runtime regression was performed with webR using R 4.6.0 (`wasm32-unknown-emscripten`).

## What happened on the first fresh run

The as-uploaded consolidated bundle passed 6 of 10 regression blocks. Four regression blocks exposed implementation/API issues rather than stochastic disagreements:

1. `INApestMeta` / `INApestMetaParallel` did not retain the original pathogen specification required later by pathogen-detection bookkeeping (`PathogenOriginal`).
2. `INApestPathogenTransitionMatrix.R` lost matrix dimensions when calculating blocked stage-transition counts; after preserving dimensions the count matrix also needed integer storage for the stochastic count code.
3. `INApestMetaParallel` calculated pathogen detections within workers but omitted `PathogenDetected` from the worker return object.
4. `INApestMetaPointParallel` and `INApestPointTransitionMatrixParallel` accepted pathogen settings indirectly through `...` but did not expose `Pathogen` explicitly in the public function signature.

A test-fixture ambiguity was also corrected in `test_INApestPathogen_Meta_MLU.R`: when the model has exactly two nodes and two land uses, `Beta = c(0,0)` is intentionally ambiguous. The test now supplies an explicit 2 x 2 node-by-land-use matrix. The resolver itself was not weakened.

## Result after minimal fixes

The unchanged 10-block consolidated regression suite then passed **10/10**:

- binary pathogen mechanisms;
- Meta/MultipleLandUse mechanisms and integration;
- point pathogen mechanisms;
- point-transition pathogen mechanisms;
- node transition-matrix pathogen mechanisms;
- pathogen detection -> information triggering;
- serial parent-function pathogen APIs;
- parallel-wrapper pathogen APIs;
- current-release `Pathogen = NULL` equivalence;
- deterministic `INApestMeta` serial/parallel-wrapper parity.

All 15 corrected source files also passed the delimiter/string/comment-aware structural audit.

## Evidence boundaries

- The `Pathogen = NULL` test compares explicit `Pathogen = NULL` with omission **within this current release**. It is not a byte-for-byte or runtime comparison against a historical pre-pathogen INApest source file.
- webR reported no usable multi-core count, so the parallel parity test exercised the single-worker wrapper/fallback path. This validates wrapper logic but is **not** a true multi-worker PSOCK validation. A native-R multi-worker check is still recommended.
- The longer behavioural/report demonstration script is included but was not completed at its high replication counts in webR because the WebAssembly run was prohibitively slow. It is intended for native R or for a deliberately lower-rep exploratory run.

See `validation_results/webr_R4.6.0/` for the fresh regression table, summary and session information.
