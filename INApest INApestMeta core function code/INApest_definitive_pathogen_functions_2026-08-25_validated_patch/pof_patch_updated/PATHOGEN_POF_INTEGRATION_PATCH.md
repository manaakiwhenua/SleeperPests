# Pathogen PoF integration patch — 2026-08-25

This patch aligns the validated pathogen-capable INApest source release with the pathogen proof-of-freedom companion.

## Changes

1. `INApestMetaTransitionMatrix.r` now saves `PathogenDetectedLargeOut.rds` in serial runs, matching the parallel output contract and retaining realised pathogen-detection histories for future conditioned-history replay.
2. `INApestPointPathogenInteraction()` now exposes its existing internal point parameter resolver as `ResolvePoint`. Existing behaviour is unchanged; this lets PoF use exactly the same detection-probability resolution semantics as the simulator.
3. `INApestPathogenPoF.R` now accepts:
   - an in-memory model result object;
   - a saved point-model result `.rds` object;
   - the common filename stem for Meta-family pathogen output arrays; or
   - `list(OutputDir=..., ModelName=...)` for Meta-family outputs.
4. Point-model observation likelihoods use `PathogenInteraction$ResolvePoint` when available, with a backwards-compatible fallback for older interaction objects.

## Biological/statistical semantics unchanged

- Freedom is derived from latent pathogen state, never from detection state.
- SIS/SIR are free when `I == 0` everywhere.
- SEIR is free when `E + I == 0` everywhere.
- Surveillance likelihood uses infectious hosts/points under the current default detection model.
- Dynamic compatible-history replay remains the next development step when pathogen detection changes information and future management.

## Validation status

The supplied validated pathogen release had already passed 10/10 webR R 4.6.0 regressions before this narrow patch. A new R regression script, `tests/test_pathogen_pof_integration_patch.R`, covers the new contracts. In the current execution environment webR still fails to initialise under Node before R starts (`__dirname is not defined`), so this new R regression script is included but has not been freshly executed here. Source-level checks and deterministic mathematical benchmarks have been rerun.
