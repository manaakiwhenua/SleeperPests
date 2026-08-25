# Validation status

## Completed in this release build

- All 15 files in `src/` passed the supplied delimiter/string/comment-aware structural audit.
- Required public model/helper symbols are present; see `SYMBOL_CHECK.txt`.
- SHA-256 checksums are supplied for every release file in `SHA256SUMS.txt`.
- Exact current GitHub blob SHAs were recorded for the relevant upstream baselines.
- The package contains focused R regression scripts for binary, Meta/MLU, point, point-transition and node transition-matrix pathogen behaviour.

## Earlier executed evidence carried into the development line

The `LocalDynamicsArgs` mechanism used by the transition-matrix work previously passed an R/webR serial/parallel regression suite, including static arguments, timestep slicing, resolver functions, collision checks and serial/parallel behaviour.

The standalone vertebrate point engine used as the interaction-capable parent for the point-transition facade was previously exercised in the vertebrate validation work, including state-specific movement and interaction behaviour.

These earlier results support the underlying mechanisms, but are not substitutes for a fresh runtime regression of this exact consolidated release.

## Not completed in the present session

A fresh end-to-end R execution of **this exact 15-file release** was not possible in the available runtime. Native `R`/`Rscript` was unavailable. Attempts to initialize the available webR runtime failed in its launcher/worker setup before R executed the regression scripts.

Accordingly:

- structural audit: **PASS**;
- current source completeness: **PASS**;
- fresh full R regression for this exact release: **NOT EXECUTED**;
- no unexecuted R test is labelled as PASS.

## Included R regression scripts

`tests/` contains:

- `test_binary_pathogen.R`
- `test_INApestPathogen_Meta_MLU.R`
- `test_Pathogen_NULL_and_parallel.R`
- `test_INApestPointPathogen.R`
- `test_INApestPointTransitionPathogen.R`
- `test_INApestPathogenTransitionMatrix.R`

These should be the first tests run in a native R environment before treating the release as runtime-validated for production use.

## Pathogen detection can trigger information

All pathogen-capable models now use the common `INApestPathogen()` fields `DetectionProb` and `DetectionTriggersInfo`. The default `DetectionTriggersInfo = FALSE` preserves earlier behaviour. When TRUE, detection of infectious pathogen state creates local information in the parent INApest model. In abundance models a per-infected-host detection probability is aggregated as `1 - (1 - p)^I`; binary occupancy uses `p` when pathogen is present; point models draw detection per infectious point. Initial pathogen detection can seed information before timestep 1, while later detections affect management from the next timestep under the existing INApest information/management timing.

### Pathogen detection -> information trigger update

Static source audit: PASS across all 15 source files after the update.
A dedicated regression script is included at `tests/test_pathogen_detection_triggers_information.R`.
The exact consolidated update has not been executed under native R in this environment because native R is unavailable; do not interpret the static audit as an end-to-end runtime PASS.
