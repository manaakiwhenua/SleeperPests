# Source provenance

Repository: `manaakiwhenua/SleeperPests`

Repository directory: `INApest INApestMeta core function code`

GitHub baselines were inspected on 25 August 2026. Where a current repository source existed, its blob SHA was recorded so the baseline is reproducible even if `main` later changes.

| Release component | Current GitHub baseline / derivation | GitHub blob SHA | Release note |
|---|---|---|---|
| `INApest.R` | `INApest.R` | `b60bb344fd637077fed7546ad13a2cf9aa3ef604` | Definitive clean/refactored binary implementation aligned with current core semantics and adding optional binary pathogen occupancy. Not byte-identical to the baseline. |
| `INApestParallel.R` | `INApestParallel.R` | `f341134cb9b10be4928c6eb6dcde9b9aed18a673` | Definitive parallel wrapper around the release serial binary engine. |
| `INApestMeta.r` | Derived serial Meta source from the current Meta code line; no top-level serial `INApestMeta.r` was available when this line was rebased | n/a | Complete serial source, including `LocalDynamicsArgs` and pathogen support. |
| `INApestMetaParallel.r` | `INApestMetaParallel.r` | `d5fac3331953a1732b5cbbbc8033bfc2e9f504e3` | Current-line parallel Meta extended with `LocalDynamicsArgs` and aggregate pathogen state. |
| `INApestMetaMultipleLandUse.r` | `INApestMetaMultipleLandUse.r` | `00002a506f907726ed676309497560983b969ba4` | Current-line MLU extended with aggregate pathogen state. |
| `INApestMetaParallelMultipleLandUse.r` | `INApestMetaParallelMultipleLandUse.r` | `26cbfd4af85f6e5ba76678d89d86d813258af91c` | Current-line parallel MLU extended with aggregate pathogen state. |
| `INApestMetaTransitionMatrix.r` | `INApestMetaTransitionMatrix.r` | `655f1a113d48158e55ee6eae9a09197673c47203` | Current-line transition matrix extended using a separate demographic-stage × pathogen-state product state. |
| `INApestMetaTransitionMatrixParallel.r` | `INApestMetaTransitionMatrixParallel.r` | `ac2ddb83646f1babc37fa3e227f29856d457566b` | Parallel counterpart of the product-state transition implementation. |
| `INApestMetaPoint.R` | `INApestMetaPoint.R` | `6a624c37df818b3c54507adf006c195674f59eac` | Complete compact refactor aligned to the current point-engine contract, with persistent `pathogen_state`. Not byte-identical to baseline. |
| `INApestMetaPointParallel.R` | `INApestMetaPointParallel.R` | `0a4a27f14acbe8f880dd6d085bee7a45be4404f4` | Current thin parallel wrapper extended to retain/combine pathogen events. |
| `INApestPointTransitionMatrix.R` | `INApestPointTransitionMatrix.R`; release facade uses the standalone validated vertebrate point-transition engine as its parent engine | `439f516edffe69a3b42ab2544e9ec0e08ac8c4b3` | Complete source. `Pathogen` is passed through the generic interaction hook. Host biology is designed to preserve parent behaviour when specialist interaction is absent; output schema may contain additional engine fields. |
| `INApestPointTransitionMatrixParallel.R` | `INApestPointTransitionMatrixParallel.R` | `0cf40e6b836267a4a9e3d43d09d4e0a74f639f94` | Complete thin parallel wrapper retaining `PathogenEvents`. |
| `INApestPathogen.R` | New common component | n/a | Shared SIS/SIR/SEIR/Binary pathogen specification and aggregate-state engine. |
| `INApestPathogenTransitionMatrix.R` | New adapter | n/a | Product-state helper layer for node transition-matrix models. |
| `INApestPointPathogen.R` | New adapter | n/a | Generic per-point pathogen interaction engine. |

## Meaning of “definitive” in this bundle

“Definitive” means these are the consolidated complete source files selected for the 25 August 2026 pathogen-capable development release: they can be sourced directly and no patch step is needed. It does **not** mean every release file is a byte-for-byte copy of GitHub plus a minimal textual diff. The binary and ordinary point serial engines were deliberately cleaned/refactored while retaining the targeted model semantics, and the point-transition facade uses the validated standalone interaction-enabled engine.

That distinction is recorded here so subsequent work can choose either this consolidated release as the new development baseline or perform a stricter line-by-line upstream merge if required.
