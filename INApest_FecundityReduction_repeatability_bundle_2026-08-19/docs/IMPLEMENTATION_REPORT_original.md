# INApest fecundity-under-management and parallelisation update

## Scope

This working-copy update implements the four requested tasks against the current top-level INApest source family used for the change:

1. Standardise parallel processing across the active parallel variants.
2. Add `FecundityReduction` to the abundance/Meta and transition-matrix model families.
3. Test serial/parallel agreement and backwards compatibility.
4. Run trial simulations demonstrating fecundity reduction across mortality levels and land uses, while exercising every agreed parameterisation form.

The update has **not** been pushed to GitHub. The supplied source files and patch are ready for review/application.

## Step 1 — Standardised parallel processing

The following parallel functions now use the same base-R architecture:

- `INApestParallel`
- `INApestParallelInAScene`
- `INApestMetaParallel`
- `INApestMetaParallelMultipleLandUse`
- `INApestMetaTransitionMatrixParallel`

The common pattern is:

- determine the available core count;
- define a one-permutation worker closure inside the model function;
- use `lapply()` when only one worker is available;
- otherwise create a base-R PSOCK cluster with `parallel::makeCluster()`;
- initialise independent worker RNG streams with `parallel::clusterSetRNGStream()`;
- execute permutations with `parallel::parLapply()`;
- stop the cluster explicitly, with `on.exit()` protection.

This removes the mixed `foreach/%dopar%` versus `parLapply()` implementation. The previous `doParallel`, `foreach` and `abind` dependencies are no longer needed by the core parallel functions. Result assembly formerly using `abind` is handled with base R `simplify2array()` while retaining the permutation dimension when `Nperm = 1`.

The transition-matrix parallel implementation no longer maintains a long manual `clusterExport()` list. Because the worker is a closure defined inside the model call, model inputs and local helper functions—including `FecundityReduction`—are part of the worker environment. This substantially reduces the risk of a new model parameter being omitted from worker setup.

`INApestParallelInAScene` calls `INA::INAscene()` explicitly inside its worker so the worker call does not rely on package attachment state.

`ParallelSetup.r` has also been simplified: base R supplies `parallel`; `doParallel`, `foreach` and `abind` are no longer core setup requirements. `INA` remains the external dependency for the InAScene variant.

## Step 2 — `FecundityReduction`

### General meaning

`FecundityReduction` is placed after the mortality arguments and before `SpreadReduction` in the relevant function signatures. Values must be finite and between 0 and 1:

- `0` = no management effect on per-capita fecundity;
- `1` = complete suppression of fecundity for managed individuals.

The order of effects is deliberate. Management mortality is applied first. `FecundityReduction` then reduces the reproductive output of individuals that survived mortality, and only where management is active. This preserves the distinction between killing individuals and reducing the fecundity of survivors.

There is no `FecundityReductionSD` parameter.

### Ordinary `INApestMetaParallel`

Supported forms are:

- scalar;
- vector of length number of nodes;
- matrix `nodes x Ntimesteps`.

The operative calculation is equivalent to:

`EffectiveReproductivePopulation = N0 * (1 - NodeFecundityReduction * Managing)`

before the existing propagule-production draw.

### Multiple-land-use Meta models

Both serial and parallel versions support:

- scalar;
- vector of length `Nlanduses`;
- matrix `nodes x Nlanduses`;
- array `nodes x Nlanduses x Ntimesteps`.

Each land use's reproductive contribution is reduced before node-level propagules are summed. The same reduced reproductive contributions are used when assigning LDD propagules back to source land uses. Consequently, a node containing several land uses can have different fecundity responses to management within the same timestep.

The serial version retains compatibility with older custom `LocalDynamics` functions when fecundity reduction is inactive. If an active non-zero `FecundityReduction` is requested, a custom dynamics function must accept `nodefecundityreduction` (or `...`), otherwise the model gives a clear error rather than silently ignoring the new management effect.

### Transition-matrix Meta models

Both serial and parallel versions support:

- scalar;
- vector of length number of nodes;
- vector of length `Nstages`;
- matrix `nodes x Ntimesteps`;
- array `nodes x Nstages x Ntimesteps`.

Fecundity is reduced stage-by-stage when the first-row fecundity contributions in the transition matrix are converted to expected propagules. The supplied transition matrix itself is not modified.

`FecundityReduction[1]` is accepted so the vector aligns naturally with stage numbering, but stage 1 is not used in the fecundity calculation because transition-matrix fecundities are the entries for stages `2:Nstages`.

If the number of nodes equals `Nstages`, an unnamed vector of that shared length is ambiguous. The model deliberately errors and asks for a dimensioned node-by-time matrix or node-by-stage-by-time array instead of guessing whether the vector means nodes or stages.

The same custom-`LocalDynamics` compatibility rule used in the multiple-land-use serial model is applied to both transition-matrix implementations.

## Step 3 — Validation

### Parsing and dependency audit

All updated source files, `ParallelSetup.r`, and both validation scripts parse successfully under WebR's R 4.6.0 parser. Static inspection confirms that the active parallel implementations no longer contain `%dopar%`, `registerDoParallel`, `library(doParallel)`, `library(abind)`, the old fixed `set.seed(12345 + i_perm)` pattern, or `closeAllConnections()`.

### Serial/parallel agreement in the executable environment

WebR reports no native cores and cannot create PSOCK worker processes. Therefore, the parallel functions select their one-core fallback in this environment. This is still a useful test because it executes the **parallel function's actual one-permutation worker closure and result-combination path**, with the same RNG stream as the corresponding serial model.

Across the tested binary INApest, multiple-land-use Meta and transition-matrix Meta models, all 12 compared saved outputs were exactly equal: every `max_abs_difference` was 0. This included population/stage-population, invasion, information/management and detection outputs as applicable.

For every supported multiple-land-use and transition-matrix `FecundityReduction` shape, matched serial and parallel-function runs were also exactly equal on this path. All nine shape-specific comparisons had `max_abs_difference = 0`.

### Backwards compatibility

The pre-change serial multiple-land-use and transition-matrix source files were retained as references. With `FecundityReduction = 0`, the new models reproduced the old population outputs exactly (`max_abs_difference = 0`). This confirms that the default value retains existing model behaviour in the tested scenarios.

### Parameterisation smoke tests

All agreed forms executed successfully:

- ordinary Meta: scalar, node vector, node-by-timestep matrix;
- multiple land use: scalar, land-use vector, node-by-land-use matrix, node-by-land-use-by-timestep array;
- transition matrix: scalar, node vector, stage vector, node-by-timestep matrix, node-by-stage-by-timestep array.

Two edge-case checks also passed: changing only transition-matrix stage 1 reduction produced exactly zero output difference, and the deliberate node/stage vector ambiguity produced the intended error.

### True PSOCK limitation and supplied native test

The supplied WebR runtime cannot instantiate native PSOCK subprocesses (`parallel::makeCluster()` is unsupported there), so I do **not** claim to have executed the multi-process branch in this environment.

`native_psock_validation.R` is included for a normal native R installation. It requires at least three detected cores so the functions actually choose more than one worker, runs binary INApest, multiple-land-use Meta and transition-matrix Meta serial/parallel comparisons, and checks that parallel and serial final-population summaries are within four combined Monte Carlo standard errors. This is the remaining runtime check before treating the multi-process execution as fully validated on the target machine.

The InAScene file parses and has the same PSOCK architecture, but its runtime path was not executed here because the WebR environment does not have the external `INA` package.

## Step 4 — Trial simulations

### Fecundity reduction by management mortality

A simple ordinary Meta scenario was run for six timesteps with 60 stochastic realisations per combination. Management mortality was 0, 0.3 or 0.6, crossed with `FecundityReduction` of 0, 0.25, 0.5, 0.75 and 1.

At mortality 0, mean final total population fell from **678.0** at no fecundity reduction to **470.2**, **305.5**, **185.3** and **100.0** as fecundity reduction increased. The final value of 100 at complete fecundity suppression is the initial total population in this no-mortality setup: survivors persist but produce no new propagules.

At mortality 0.3, the corresponding means were **101.1**, **60.6**, **37.3**, **22.5** and **11.4**. At mortality 0.6 they were **3.15**, **2.47**, **1.43**, **0.90** and **0.52**. Thus mortality and fecundity suppression act together as intended rather than one replacing the other.

The full 5th–95th percentile results are in `mortality_by_fecundity_results.csv` and visualised in `mortality_fecundity_effect.png`.

### Land-use-specific fecundity reduction

The land-use demonstration used three land uses plus one mixed node, again over six timesteps and 60 realisations. A gradient scenario assigned fecundity reductions of **0, 0.5 and 0.9** to land uses 1–3.

At mortality 0.3, final mean population in the single-land-use nodes was **20.5** for land use 1 (`FR=0`), **10.5** for land use 2 (`FR=0.5`) and **3.9** for land use 3 (`FR=0.9`). This is the intended within-model ordering and demonstrates that the land-use vector is being applied to reproductive contributions rather than collapsed to one node-level value.

For comparison, total mean final population at mortality 0.3 was **114.8** with no fecundity reduction, **42.7** with a uniform 0.5 reduction, **49.2** with the 0/0.5/0.9 land-use gradient, and **16.4** with a uniform 0.9 reduction.

The complete results are in `landuse_mortality_fecundity_results.csv`; the land-use gradient is illustrated in `landuse_fecundity_gradient_effect.png`.

## Files supplied

### Updated source files

- `INApestParallel.R`
- `INApestParallelInAScene.R`
- `INApestMetaParallel.r`
- `INApestMetaMultipleLandUse.r`
- `INApestMetaParallelMultipleLandUse.r`
- `INApestMetaTransitionMatrix.r`
- `INApestMetaTransitionMatrixParallel.r`
- `ParallelSetup.r`

`INApest.R` itself was not changed; it was used as the serial reference for the binary-presence parallelisation check.

### Validation and evidence

- `validation_fecundity_reduction.R` — executed WebR validation/demonstration suite.
- `native_psock_validation.R` — remaining true-multiworker native-R validation.
- `equivalence_results.csv`
- `parameterisation_smoke_results.csv`
- `parameterisation_parallel_equivalence.csv`
- `mortality_by_fecundity_results.csv`
- `landuse_mortality_fecundity_results.csv`
- `validation_environment.csv`
- `mortality_fecundity_effect.png`
- `landuse_fecundity_gradient_effect.png`
- `inapest_fecundity_parallel.patch` — unified patch against the downloaded pre-change top-level files.

## Recommended final check before merging

Run `source("native_psock_validation.R")` from the updated source directory on a native R installation with at least three detected cores. If that passes, the code has both exact one-core serial/parallel equivalence evidence and an actual multi-process PSOCK smoke/statistical-equivalence check.
