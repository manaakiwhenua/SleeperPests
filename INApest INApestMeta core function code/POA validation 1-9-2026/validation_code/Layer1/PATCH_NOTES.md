# INApest PoA observation-engine patch — 31 August 2026

## Purpose

This patch implements the PoF-style separation between biological state, surveillance observations, information state and management response for INApest proof-of-absence (PoA).

The patch deliberately does **not** calculate PoA inside the biological engines. Engines simulate and retain the information needed by the PoA companion; `INApestPoA.R` calculates observation likelihoods and performs Bayesian conditioning.

## Canonical surveillance streams

1. **Background** surveillance: routine/passive surveillance that operates regardless of prior information.
2. **InfoTriggered** surveillance: additional targeted surveillance/search that operates only where information existed **before the current host-surveillance round**.

`HaveInfo` remains a response state. It is not evidence by itself.

## Timing rule

For each host-surveillance round the engine freezes the pre-surveillance information state before either stream is evaluated. A new Background detection can create information for later response, but it cannot retroactively activate InfoTriggered surveillance in that same round.

## New public arguments

Node/binary models append:

- `InfoTriggeredDetectionProb = 0`
- `InfoTriggeredDetectionSD = NULL`

Point models additionally append:

- `InfoTriggeredDetectionSpatial = NULL`

`INApestMeta`, MLU and transition-matrix serial/parallel models additionally append:

- `ReturnResults = FALSE`

The default targeted-surveillance pathway is off. Existing public argument order is preserved; all new arguments are appended.

## New engine outputs

Node-family engines retain:

- `BackgroundDetectedResults`
- `InfoTriggeredDetectedResults`
- `InformationStateBeforeSurveillanceResults`
- `HaveInfoResults`
- `BackgroundDetectionProbabilityResults`
- `InfoTriggeredDetectionProbabilityResults`

The probability outputs are **realised per-unit/per-individual detection probabilities**, including any `DetectionSD` uncertainty. They are not surveillance-system sensitivity and they are not PoA `q` values.

For MLU the event arrays are node x timestep x permutation, while realised detection probabilities remain node x land-use x timestep x permutation.

For transition-matrix models the event arrays are node x timestep x permutation, while realised detection probabilities remain node x stage x timestep x permutation.

Point engines retain corresponding event counts in `Summary` and, for every point snapshot:

- `have_info_before_surveillance`
- `background_detection_prob`
- `info_triggered_detection_prob`

They also retain stream-specific events in `EventLog`.

## Legacy outputs

`DetectedResults` is retained unchanged as the legacy persistent known-present state (`HaveInfo * Invaded`) in node models. The PoA companion no longer treats it as a reliable new-detection event for sequential compatible-history inference.

Existing standard saved outputs are retained. New output files are additive.

## PoA companion changes

`INApestPoA.R` now:

- uses `Background` and `InfoTriggered` as canonical evidence source names;
- accepts `Surveillance` and `Information` only as deprecated aliases;
- distinguishes event history from `HaveInfo` state;
- supports all six biological architectures and their parallel wrappers;
- restores a serial `INApestMetaPoA()` wrapper now that Meta has `ReturnResults`;
- can calculate likelihood from a user-supplied observation model;
- automatically calculates no-detection likelihood from **recorded engine detection probabilities** when patched model results are supplied;
- preserves legacy stored no-detection-probability support;
- allows hybrid recorded/stored observation sources for frozen analyses;
- prevents sequential propagation from relying on persistent legacy `DetectedResults`.

The automatic recorded-probability likelihood is:

`q = product((1 - p)^N)`

across the relevant observation units, with `InfoTriggered` probabilities multiplied by the pre-surveillance information gate. For point models the equivalent product is evaluated directly over the points present at the surveillance round.

## Backward-compatibility intent

The patch is additive. With `InfoTriggeredDetectionProb = 0` and `InfoTriggeredDetectionSD = NULL`, the new targeted stream is not sampled. New probability/event bookkeeping is designed not to consume additional random draws on the default pathway.

This source-level claim still requires the supplied native-R regression run before the patch should be described as runtime-validated.

## Validation state in this bundle

Completed here:

- static delimiter/source-contract check across all 12 serial/parallel core files plus the integrated PoA companion: PASS;
- public-formal compatibility audit against `INApest_latest_source_GitHub_handoff_2026-08-27`: PASS for every directly changed public biological function;
- unchanged point-parallel wrapper hashes equal the authoritative baseline.

Not completed in this environment:

- native R execution of the supplied runtime suites;
- true multi-worker PSOCK comparison;
- full seeded baseline-versus-patched stochastic equality checks.

Run `RUN_VALIDATION.ps1` in native R before promoting this patch to the definitive repository state.

## Native-R MLU adapter correction — 31 Aug 2026
The first native-R architecture run reached 26/29 PASS. B05, R04 and D03 all traced to one MLU compatibility-shape issue: `DetectedResults` is a node x timestep x permutation observation/state array in the patched MLU engine, while the MLU adapter's legacy compatibility line still tried to aggregate it using the biological node x land-use x timestep x permutation dimensions. The adapter now routes legacy `DetectedResults` through `.INApestPoAAggregateDetectionEvent()`, the same shape-aware helper used for the explicit Background and InfoTriggered event streams. No biological engine code or likelihood mathematics changed in this correction.

## Native-R correction after engine validation (2026-08-31)

A second native-R pass exposed two implementation issues that were not visible in static checks:

1. `INApest.R` used `identical(dim(x), c(...))` for some matrix/array validation. A valid nodes x timesteps matrix could therefore be rejected solely because `Ntimesteps` was numeric rather than integer. Dimension validation now compares dimension values rather than R storage mode.
2. Plain `INApestMeta` and `INApestMetaParallel` packed multiple one-node/one-timestep worker matrices through nested `simplify2array()` calls. R could collapse singleton dimensions before unpacking. Workers now return named lists and results are bound explicitly to nodes x timesteps x permutations arrays, matching the robust pathogen-result pattern and preserving singleton dimensions.

These changes do not alter surveillance probabilities, biological dynamics, or random-number draws.
