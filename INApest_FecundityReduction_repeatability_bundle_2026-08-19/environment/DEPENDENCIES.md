# Runtime dependencies

## Core fecundity-validation workflow

The ordinary Meta, multiple-land-use Meta and transition-matrix Meta functions used in this bundle rely on base R only. The standardised parallel implementations use the base-R `parallel` package and PSOCK clusters. No `foreach`, `doParallel` or `abind` packages are required by these current core functions.

The archived and current one-core validation evidence was run under WebR reporting **R 4.6.0 (2026-04-24)**. WebR cannot start native PSOCK worker processes, so that environment exercises the parallel functions' one-core fallback, worker body and result-combination logic but not multi-process execution.

Run `scripts/05_native_psock_validation_current.R` in a normal native R installation with at least three detected cores to force the actual PSOCK branch.

## Other source files in the current snapshot

`INApestParallelInAScene.R` requires the external `INA` package because it calls `INA::INAscene()`. It is retained in the current source snapshot because parallelisation was standardised across the broader source family, but it is not used by the FecundityReduction experiments.

## Report reproduction

The final HTML report is self-contained. `report/FecundityReduction_repeatability_report.qmd` is supplied as an editable report source. Rendering the QMD requires Quarto; its analysis tables and figures are read from the saved CSV/PNG outputs in this bundle.
