INApest latest source-file bundle
Prepared 2026-08-19

FILES
-----
INApest.R
  Binary-presence serial INApest from the current GitHub main branch.

INApestParallel.R
  Binary-presence parallel INApest using the harmonised base-R PSOCK/parLapply architecture.

INApestParallelInAScene.R
  InAScene binary parallel variant using the harmonised base-R PSOCK/parLapply architecture.

INApestMeta.r
  Ordinary abundance/Meta serial model.
  Includes FecundityReduction and a user-definable LocalDynamics hook.
  The current default local.dynamics function is defined in this same source file.

INApestMetaParallel.r
  Parallel counterpart of INApestMeta.r using the harmonised PSOCK/parLapply architecture.
  Includes FecundityReduction and user-definable LocalDynamics.

INApestMetaMultipleLandUse.r
  Multiple-land-use serial Meta model.
  Includes land-use-aware FecundityReduction and user-definable LocalDynamics.
  The current default local.dynamicsLU function is defined in this same source file.

INApestMetaParallelMultipleLandUse.r
  Parallel counterpart of the multiple-land-use Meta model using the harmonised PSOCK/parLapply architecture.

INApestMetaTransitionMatrix.r
  Stage-structured transition-matrix serial Meta model.
  Includes FecundityReduction; user-definable LocalDynamics; and optional movement during stage transitions using separate TransitionSDDprob and TransitionLDDprob inputs, including time-varying matrices/arrays.
  The current default local.dynamics.transition.matrix function is defined in this same source file.

INApestMetaTransitionMatrixParallel.r
  Parallel counterpart of the transition-matrix Meta model using the harmonised PSOCK/parLapply architecture.

ParallelSetup.r
  Updated parallel setup helper consistent with the base-R PSOCK implementation.

NOTES
-----
* FecundityReduction = 0 retains the previous Meta fecundity behaviour.
* Transition movement is opt-in; omitting the transition-movement arguments retains local stage transitions.
* Custom LocalDynamics functions are supported in all Meta families.
* Multi-process parallel execution uses base R parallel/PSOCK rather than foreach/doParallel.
