# Validation status — validated patch, 25 August 2026

## Fresh runtime validation completed

The corrected pathogen-capable source set was executed under webR with **R 4.6.0** on `wasm32-unknown-emscripten` using the consolidated regression runner in `tests/run_consolidated_pathogen_regression_2026-08-25.R`.

Result: **10/10 regression blocks PASS; 0 FAIL; 0 SKIP.**

The regression blocks cover:

1. binary pathogen occupancy, directed transmission, clearance, pathogen introduction and pathogen-associated host extinction;
2. SIS/SIR/SEIR state accounting, recruitment/thinning reconciliation, progression, recovery, waning immunity, pathogen mortality, pathogen introduction, land-use/time resolution and Meta/MLU integration;
3. point pathogen state, contact/transmission, recruits, progression, recovery, waning, mortality and introduction;
4. independence of demographic and pathogen state in the point transition-matrix model;
5. demographic-stage x pathogen-state accounting, cross-stage transmission and stage-transition movement in the node transition-matrix model;
6. pathogen detection triggering INApest information only when requested;
7. serial parent-function pathogen APIs;
8. parallel-wrapper pathogen APIs;
9. current-release `Pathogen = NULL` equivalence;
10. deterministic serial versus parallel-wrapper `INApestMeta` parity.

Full results and session information are in `validation_results/webr_R4.6.0/`.

## Structural audit

All **15/15** corrected source files passed the delimiter/string/comment-aware source audit. The fresh corrected-source audit is stored at `validation_results/static_audit_corrected_sources.txt`.

## Fixes exposed by the fresh run

The first run of the as-uploaded consolidated bundle passed 6/10 regression blocks and exposed four small implementation/API defects plus one ambiguous MLU test fixture. The minimal corrections and their rationale are documented in `VALIDATION_FIXES_2026-08-25.md`.

## Important limits

### Historical backward compatibility

The current regression verifies that explicitly using `Pathogen = NULL` is equivalent to omitting `Pathogen` under the same current source and seed. The supplied older `test_Pathogen_NULL_and_parallel.R` cannot provide a historical pre-pathogen comparison because its referenced `baseline/INApestMeta.r` and some development filenames are absent from the release. A true historical-source regression remains a separate task.

### Multi-worker parallel execution

webR did not expose a usable `parallel::detectCores()` value. The deterministic serial/parallel-wrapper comparison therefore validates the wrapper/fallback path, **not true multi-worker PSOCK execution**. A native-R multi-worker regression remains recommended before making a definitive multi-worker claim.

### Behavioural parameter-response demonstrations

The high-replication report demonstration script is supplied as `tests/report_behavioural_validation_highrep_2026-08-25.R`. Its full high-replication run was not completed in webR because the WebAssembly runtime was prohibitively slow. This is separate from the completed 10/10 regression suite.

## Earlier supporting disease-behaviour evidence

Earlier development validation also exercised generic epidemiological behaviours such as epidemic thresholds/stochastic fade-out, infectious-period effects, contact formulation, movement, seasonal effects, waning immunity and qualitative migration/flyway/reseeding behaviour. These are supporting evidence for the architecture, but should remain clearly distinguished from the fresh regression of this corrected consolidated build.

## PoF integration patch

A subsequent narrow proof-of-freedom integration patch added `INApestPathogenPoF.R`, serial transition-matrix pathogen-detection persistence, and exposure of the existing point resolver. All core source files pass the structural audit after these changes. Deterministic binary, SEIR, stage/land-use aggregation, and point-surveillance mathematical benchmarks pass. `tests/test_pathogen_pof_integration_patch.R` is included for fresh R execution; the current Node/webR launcher fails before R initialisation in this environment, so that new R regression is not claimed as executed here.
