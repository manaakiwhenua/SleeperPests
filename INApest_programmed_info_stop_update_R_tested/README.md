# INApest programmed information-stop update

This bundle contains the eight updated INApest / INApestMeta production functions plus regression tests.

## New parameter

`InfoPersistenceSteps = NA`

Accepted forms:

- single non-negative whole number;
- vector of length `nodes`;
- matrix `nodes x Ntimesteps`;
- `NA` at any node/timestep means use `InfoRetentionProb` there instead.

`InfoPersistenceSteps = NA` everywhere is the backward-compatible default.

A finite value gives the number of *additional timesteps after the most recent known local presence* for which information remains available for management. For example, if last known local presence occurs at timestep 1 and `InfoPersistenceSteps = 2`, information is available for management in timesteps 1, 2 and 3 and is stopped after timestep 3 unless refreshed.

When both programmed persistence and stochastic decay are supplied, programmed persistence takes priority for node/timestep combinations where `InfoPersistenceSteps` is finite. A single warning is emitted at function entry. `InfoRetentionProb` applies only where `InfoPersistenceSteps` is `NA`.

## What resets local known presence

- `INApest`, `INApestParallel`, `INApestParallelINAscene`: an informed node with an extant infestation.
- INApestMeta population models: a realised management-caused kill in the node.
- Direct pest detection also records local known presence in all models.
- Information arriving from SEAM or external information sources refreshes `HaveInfo`, but does not itself create a local-presence record.

For the non-transition Meta models, natural and management mortality are combined in the existing population draw. To avoid rewriting that draw, the update conditions on the realised deaths and samples whether at least one was management-caused. When management is the only possible cause, the reset is deterministic.

## Ordering within a timestep

1. existing information is available for management;
2. management and population/spread operations run;
3. last-known-local-presence is updated from extant known occupancy (INApest) or management kills (Meta);
4. programmed information stopping is applied;
5. stochastic `InfoRetentionProb` decay is applied only where no programmed stop is supplied;
6. SEAM, external information and direct detection can refresh information;
7. direct detection also resets the local-presence clock.

## Other correctness fixes found during testing

- `INApestMetaTransitionMatrix.r`: current-timestep detection now uses `NodeDetectionProb[,,timestep]` rather than the stale earlier loop index `t`.
- `INApestMetaMultipleLandUse.r`: SEAM communication now collapses the node x land-use `Detected` matrix to node level with `apply(Detected,1,max)`, matching the parallel multiple-land-use implementation.
- `INApestParallelInAScene.R`: when information loss is active, SEAM refresh is applied after the stop/decay step so incoming information can restore a node whose information just expired. The default no-loss path remains unchanged.

## Tests

The main programmed-stop tests use a 30 x 30 landscape (900 nodes), four timesteps and two stochastic realisations. They cover all eight public implementations and verify:

- exact backward compatibility when `InfoPersistenceSteps` is omitted;
- scalar persistence;
- zero-step stop;
- node-vector persistence;
- node x timestep persistence;
- priority over stochastic decay with exactly one warning;
- mixed programmed-stop / stochastic-decay nodes;
- local evidence resetting the clock;
- Meta management-kill versus no-kill behaviour;
- invalid dimension and non-whole-value rejection;
- cross-node SEAM refresh after local information expiry.

The earlier information-decay and temporal-connectivity regressions were rerun on the final code and still pass.

Serial models were executed directly in WebR/R 4.6.0. Parallel model bodies were executed with sequential test shims because the WebR environment does not provide the normal PSOCK/doParallel/INA package stack. This validates the model logic but is not a full concurrent-worker or real INA integration test.
