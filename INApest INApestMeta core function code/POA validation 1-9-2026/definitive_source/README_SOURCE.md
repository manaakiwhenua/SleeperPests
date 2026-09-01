# Definitive frozen source

This folder contains the source state frozen after the completed PoA validation programme.

## Biological simulation engines

- `INApest.R`
- `INApestParallel.R`
- `INApestMeta.r`
- `INApestMetaParallel.r`
- `INApestMetaMultipleLandUse.r`
- `INApestMetaParallelMultipleLandUse.r`
- `INApestMetaTransitionMatrix.r`
- `INApestMetaTransitionMatrixParallel.r`
- `INApestMetaPoint.R`
- `INApestMetaPointParallel.R`
- `INApestPointTransitionMatrix.R`
- `INApestPointTransitionMatrixParallel.R`

## Shared simulation helpers retained with the definitive core set

- `INApestPathogen.R`
- `INApestPathogenTransitionMatrix.R`
- `INApestPointPathogen.R`

These helpers are not required for ordinary non-pathogen PoA, but are part of the definitive pathogen-enabled INApest simulation source on which the patched engines were based.

## PoA wrapper and general PoA helpers

- `INApestPoA.R`

The common Bayesian PoA core, model-family adapters, observation-model helpers, Background/InfoTriggered likelihood handling, compatible-history logic and user-facing PoA wrappers are contained in this file. There is no additional separate general PoA helper file in the frozen implementation.

See `SOURCE_ORDER.md` for sourcing guidance and `OBSERVATION_ARCHITECTURE.md` for the surveillance contract.
