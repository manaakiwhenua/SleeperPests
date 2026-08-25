# INApest definitive pathogen-capable source bundle — 25 August 2026

This bundle is the consolidated source release produced from the current INApest development line plus the generic pathogen component developed in this work.

## What is in `src/`

Every file in `src/` is a complete R source file. There are no patch scripts or apply-to-GitHub steps in the source directory.

The generic pathogen architecture covers:

- binary `INApest`: binary host occupancy plus binary pathogen occupancy;
- `INApestMeta`: aggregate SIS/SIR/SEIR-style pathogen state alongside total host abundance;
- `INApestMetaMultipleLandUse`: the same aggregate state with land-use structure;
- `INApestMetaTransitionMatrix`: demographic-stage abundance plus a separate demographic-stage × pathogen-state product state;
- `INApestMetaPoint`: persistent pathogen state on explicit points;
- `INApestPointTransitionMatrix`: demographic stage and pathogen state retained independently on each point;
- serial and parallel counterparts.

The common specification is `INApestPathogen()` in `INApestPathogen.R`. Point models additionally use `INApestPointPathogenInteraction()` from `INApestPointPathogen.R`. Node transition-matrix models use helpers in `INApestPathogenTransitionMatrix.R`.

## Recommended source order

For abundance Meta/MLU models, each model source currently embeds the common pathogen helpers, so it may be sourced directly. For a common explicit setup, source `INApestPathogen.R` first and then the desired model file.

For point models:

1. `INApestPathogen.R`
2. `INApestPointPathogen.R`
3. `INApestMetaPoint.R` or `INApestPointTransitionMatrix.R`
4. the corresponding parallel wrapper if required.

For node transition-matrix models:

1. `INApestPathogen.R`
2. `INApestPathogenTransitionMatrix.R`
3. `INApestMetaTransitionMatrix.r`
4. optionally `INApestMetaTransitionMatrixParallel.r`.

## Backward-compatibility principle

The pathogen component is optional. The intended contract is that absence of a pathogen specification follows the ordinary host-model path. Existing host demography, dispersal, capacity, information and management remain owned by the INApest core. The pathogen component maintains pathogen state and only changes host abundance through an explicitly specified pathogen-associated host-loss process.

## Validation status

See `VALIDATION_STATUS.md`. All 15 source files passed the supplied structural source audit. The full current pathogen family has **not** been executed end-to-end in R in this session because a working native R runtime was unavailable and the available webR launcher failed before R initialization. The included R regression tests are therefore supplied for execution, not reported as fresh PASS results.

## Provenance

See `PROVENANCE.md` for exact GitHub blob SHAs and notes on which release files are direct extensions versus derived/refactored implementations.

## Pathogen detection can trigger information

All pathogen-capable models now use the common `INApestPathogen()` fields `DetectionProb` and `DetectionTriggersInfo`. The default `DetectionTriggersInfo = FALSE` preserves earlier behaviour. When TRUE, detection of infectious pathogen state creates local information in the parent INApest model. In abundance models a per-infected-host detection probability is aggregated as `1 - (1 - p)^I`; binary occupancy uses `p` when pathogen is present; point models draw detection per infectious point. Initial pathogen detection can seed information before timestep 1, while later detections affect management from the next timestep under the existing INApest information/management timing.
