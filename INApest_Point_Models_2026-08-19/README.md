# INApest point-model development package — 2026-08-19

This package contains the current point-based invasion modelling engines, a native-R animated demonstration, regression tests and two HTML reports.

## Start here

1. Open **`INApest_Point_Model_Executive_Summary.html`** for the concise overview and embedded animation.
2. Open **`INApest_Point_Model_Technical_Report.html`** for the expert-modeller report.
3. Run **`QUICKSTART.R`** from the extracted package root in R/RStudio.

## Current model names

- `INApestMetaPoint()` — canonical point-based analogue of INApestMeta. `INApestPoint()` remains an alias for compatibility.
- `INApestPointTransitionMatrix()` — stage-structured point engine with optional movement during adjacent stage progression.

Both serial engines support native base-R spatial grids and spatial schedules for habitat and management. GIS libraries are optional upstream tools for preparing numeric grids; they are not required by the core simulation.

## Folder guide

- `models/` — serial engines and parallel wrappers
- `demos/` — native-R two-stage animation script and GIF
- `outputs/` — validation and illustrative simulation CSVs
- `tests/` — regression scripts and available logs
- `docs/` — compatibility/status notes

## Validation status

Serial engines, base-R habitat logic, management spatial modifiers, spatial schedules, transition-matrix mechanisms and the `INApestMetaPoint` rename have passed the included R 4.6 WebAssembly checks. Genuine multi-process PSOCK execution cannot be launched in that environment; the included parallel test performs the native two-worker check when run in ordinary R with at least two cores.
