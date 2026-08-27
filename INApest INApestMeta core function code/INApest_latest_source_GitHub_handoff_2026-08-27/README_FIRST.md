# INApest latest source GitHub handoff — 2026-08-27

This bundle is for manual upload to:

`manaakiwhenua/SleeperPests/INApest INApestMeta core function code/`

## What to upload

Upload the **files inside**:

`github_upload/INApest INApestMeta core function code/`

Those files are the latest consolidated source set available at the end of the 2026-08-27 analytical/validation work. Later validated fixes override earlier source snapshots.

The folder `supporting_analytical_modules_NOT_REQUIRED_FOR_RUNTIME/` contains modular development sources used to build the final unified `INApestAnalytical.R`. They are included for provenance/maintainability but are **not required for runtime if the unified analytical source is uploaded**.

## Important latest fixes included

- Binary `InformationAcquisition = host/pathogen/both` support in `INApest.R`.
- Meta information-acquisition changes in `INApestMeta.r`.
- MLU information-acquisition changes in `INApestMetaMultipleLandUse.r`.
- Windows PSOCK LocalDynamicsArgs resolver capture in:
  - `INApestMetaParallel.r`
  - `INApestMetaParallelMultipleLandUse.r`
  - `INApestMetaTransitionMatrixParallel.r`
- MetaPoint parallel `parLapply` argument fix in `INApestMetaPointParallel.R` (`model_fun`, not conflicting `fun`).
- Vertebrate Point custom-Birth zero-offspring edge-case fix in `INApestPointTransitionMatrix.R`.
- Unified `INApestAnalytical.R` containing the latest analytical family through vertebrate host + pathogen methods.

## Files that are current but unchanged relative to the prior definitive source handoff

Some files are included even where no new source edit was required, so this bundle can be treated as a complete current source set rather than a patch-only bundle.

## Validation status relevant to the latest additions

Native Windows R 4.4.1 validation completed successfully for the major latest analytical layers:

- MLU pathogen: all validation blocks passed after targeted test corrections.
- Transition-Matrix pathogen: 7/7 PASS, including PSOCK Cores=2.
- Point pathogen (MetaPoint + PointTransitionMatrix): 9/9 PASS after FIX1, including PSOCK Cores=2.
- Vertebrate host analytical: 8/8 PASS after zero-offspring source hotfix.
- Vertebrate + pathogen analytical: 8/8 PASS.

See `SOURCE_MANIFEST.csv` for source provenance and upload action.
