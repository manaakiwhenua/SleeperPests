# Test and validation summary

## Production code corrections

The cumulative corrected source contains the following fixes identified during this development sequence:

1. Serial transition management mortality used stale loop index `t`; corrected to `timestep`.
2. Serial transition detection used stale loop index `t`; corrected to `timestep`.
3. Serial multiple-land-use SEAM multiplied a node x node matrix by node x land-use detection without first collapsing detection to node level; corrected to use node-level detection.
4. Time-varying `EnvEstabProb` could be used before `BPAM` was initialized in binary INApest; corrected in `INApest`, `INApestParallel`, and `INApestParallelINAscene` by initializing from timestep 1.
5. Serial multiple-land-use LDD calculated population-share weighting after the multinomial dispersal draw; corrected so `Pn` enters `Qout` before `Qin` is generated.
6. `Nperm = 1` could drop the realisation dimension. Output aggregation now uses `drop=FALSE`, and the older foreach parallel implementations explicitly restore/reshape singleton permutation output.

The incremental patch in this release starts from the preceding programmed-info-stop bundle, so items 1–3 were already present in that baseline and are not repeated in the patch.

## Direct bug tests

WebR / R 4.6.0 tests completed with the following results:

- `Nperm = 1` completed for all eight public function versions.
- Saved invasion RDS dimensions were explicitly checked:
  - standard models: nodes x timesteps x 1;
  - MLU models: nodes x land uses x timesteps x 1.
- Time-varying node x timestep `EnvEstabProb` completed in all three binary INApest implementations.
- Corrected serial MLU LDD diagnostic: simulated destination population = 18.553, versus analytical corrected expectation ~18.85. The pre-fix simulator was ~36.0 in the earlier diagnostic.

Raw results: `results/bug_fix_results.csv`.

The combined WebR parallel-shim harness raises a WebR non-local-control exception after all PASS markers and CSV output have been written. This is the same runtime teardown behaviour observed in earlier parallel tests and is not an R model assertion failure. The saved completion log is included as `results/bug_progress.txt`.

## User-facing analytical interface

`tests/test_user_facing.R` passed for:

- binary INApest baseline;
- explicit binary outside-edge escape;
- INApestMeta;
- INApestMetaTransitionMatrix including intrinsic local lambda;
- corrected INApestMetaMultipleLandUse semantics;
- dynamic information smoke test;
- time-varying connectivity smoke test.

Raw summary: `results/user_facing_summary.csv`.

## Branching escape upgrade

The final analytical release uses lineage branching/no-escape recursion as the primary escape probability and keeps the cumulative first-moment Poisson value as a comparator.

Eight-step analytical examples:

| Case | Branching | First-moment Poisson |
|---|---:|---:|
| Meta baseline | 0.1440 | 0.1595 |
| Meta managed | 0.0709 | 0.0766 |
| Transition baseline | 0.2219 | 0.2484 |
| Transition managed | 0.1212 | 0.1323 |
| MLU baseline | 0.1366 | 0.1595 |
| MLU managed | 0.0583 | 0.0629 |

Raw results: `results/branching_escape_comparison.csv`.

The Meta branching implementation was specifically corrected during final testing so integer dispersers are allocated over the reconstructed internal + external multinomial kernel. Separately truncating only the outside fraction can incorrectly turn small but real export pathways into zero.

## Earlier simulation validation retained in this release

### Binary INApest

- baseline extinction: simulation 0.555; branching 0.551;
- all-informed management: 0.690 vs 0.667;
- explicit outer-sink escape at timestep 8: simulation 0.0825, branching 0.1018, first-moment Poisson 0.1438;
- managed outer-sink escape: 0.0300, 0.0367, 0.0482.

### INApestMeta

- baseline multiplier 0.8617; simulation final mean 0.537; analytical 0.420; simulation extinction 0.706; branching 0.688;
- all-informed management multiplier 0.7142; simulation final mean 0.157; analytical 0.134; extinction 0.854 vs branching 0.873;
- focused 1,000-realisation information test: strong SEAM reduced final mean from 0.226 to 0.172 and increased extinction from 0.799 to 0.840.

### Transition matrix

- baseline intrinsic lambda 1.4225; landscape multiplier 1.0576; simulation extinction 0.475; branching 0.500;
- all-informed management intrinsic lambda 1.2173 but landscape multiplier 0.9034; extinction 0.713 vs branching 0.720.

### Multiple land use

- baseline multiplier 0.8540; simulation final mean 0.5583; analytical 0.5610;
- corrected LDD diagnostic confirms the population-share fix as above.

Raw validation CSVs are stored under `results/phase1/` and `results/phase2/`.

## Remaining test limitation

WebR does not reproduce the normal installed `abind`/`doParallel`/`INA` package stack and PSOCK worker environment. Parallel model bodies were therefore exercised with sequential test shims. A conventional-R parallel smoke test remains recommended before merging to a production branch.

## Programmed information stopping: analytical age-state extension (14 August 2026)

`INApestAnalytical()` now represents finite `InfoPersistenceSteps` values with explicit time-since-local-evidence states rather than reducing the programmed stop to a memoryless retention probability.

WebR / R 4.6.0 checks passed for:

- exact finite-state clock timing for a hand-calculated one-type case;
- continual local-evidence reset for an informed extant binary-INApest infestation;
- expiry of non-local/pre-existing information in the absence of new local evidence;
- a full INApestMeta example with detection followed by a two-step programmed stopping window;
- all four analytical model families;
- scalar, node-specific and node x timestep `InfoPersistenceSteps` inputs;
- priority of finite `InfoPersistenceSteps` over `InfoRetentionProb`;
- exact backward compatibility with the previous analytical release when `InfoPersistenceSteps = NA` across all four model families.

The deterministic clock itself is therefore explicit in the rare-state operator. The remaining approximation is shared node-level evidence at finite prevalence: in Meta/transition/MLU, a management kill can leave information at a node after the killed focal individual is gone, and informed-but-uninvaded nodes can precondition later recolonisation. Those effects require a richer node-level finite-prevalence mean-field or stochastic simulation for exact treatment.
