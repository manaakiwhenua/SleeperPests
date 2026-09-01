# Validation code

This folder consolidates the corrected code used for the three-layer PoA validation programme.

- `Layer1/` — architecture, engine runtime and backward-compatibility checks.
- `Layer2/` — exact internal Bayesian / finite-state behavioural benchmarks.
- `Layer3/` — external published benchmarks of increasing complexity.
- `RUN_ALL_VALIDATION.ps1` — convenience runner for all three layers.

## Consolidated layout note

The original validated layer packages each carried a local `source/` copy. In this frozen bundle there is one authoritative copy in `../definitive_source/` instead. The R validation scripts were changed only in their source-file discovery lines so that they can resolve this shared folder. The benchmark calculations, assertions, tolerances and test definitions were not changed. See `VALIDATION_LAYOUT_ONLY_DIFF.patch`.

Layer 1's `validation/baseline_2026-08-27/` is intentionally retained because those files are test fixtures used to prove default-path backward compatibility. They are not definitive source.
