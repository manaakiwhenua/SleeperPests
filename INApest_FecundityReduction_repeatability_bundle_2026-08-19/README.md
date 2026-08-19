# INApest FecundityReduction repeatability bundle

**Bundle date:** 19 August 2026  
**Purpose:** freeze the implementation, validation evidence and demonstration simulations for the INApestMeta `FecundityReduction` upgrade, while also showing that the same evidence is reproduced by the current later source defaults.

Start with:

- `report/FecundityReduction_repeatability_report.html` — detailed self-contained technical report for modellers who do not already know INApest.
- `inputs/scenario_manifest.csv` — compact description of the validation and demonstration scenarios.
- `inputs/fecundity_parameterisation_contract.csv` — allowed input shapes and model interpretation.

## What is frozen here

Two source snapshots are retained deliberately:

1. `source/fecundity_validation_snapshot/` is the exact code state used for the original fecundity validation evidence.
2. `source/latest/` is the current source set after later compatible work, including restored user-supplied `LocalDynamics` hooks and transition movement in the transition-matrix model.

The original validation workflow was rerun against the latest source defaults under the recorded WebR environment. All six archived CSV evidence files were reproduced byte-for-byte. A separate serial-only current-source reproduction also regenerated the three scientific parameterisation/effect tables exactly.

## Main reproduction routes

### Native R — scientific results

From the bundle root:

```r
Rscript scripts/00_run_all.R
```

On a normal multi-core machine this reproduces the scientific parameterisation and demonstration tables using current **serial** Meta functions, avoiding differences caused solely by parallel random-number streams. It also regenerates the figures.

### Native R — actual multi-worker PSOCK

```r
Rscript scripts/05_native_psock_validation_current.R
```

This requires at least three detected cores. It forces the parallel functions to use multiple PSOCK workers and assesses serial/parallel Monte Carlo agreement rather than exact realisation-by-realisation equality.

### Exact archived one-core validation environment

```text
python environment/run_webr_exact_reproduction.py
```

The working WebR runtime used for the repeatability check is retained under `environment/webr_runtime/`. The harness additionally requires Python Playwright and a usable Chromium/Chrome browser. On bundle assembly it reproduced all six archived evidence CSVs byte-for-byte.

## Evidence directories

- `results/archived/` — frozen original evidence.
- `results/latest_rerun/` — archived validation workflow rerun against latest source defaults.
- `results/native_serial_reproduction/` — serial current-source reproduction of scientific tables.
- `results/webr_exact_rerun/` — repeat run using the WebR environment retained in this bundle.

## Important interpretation

`FecundityReduction` is a proportional reduction in per-capita fecundity **after survival/management mortality and only where management is active**. It is separate from mortality and from `SpreadReduction`.

- ordinary Meta: scalar, node vector or node × timestep matrix;
- multiple-land-use Meta: scalar, land-use vector, node × land-use matrix or node × land-use × timestep array;
- transition-matrix Meta: scalar, node vector, stage vector, node × timestep matrix or node × stage × timestep array.

The default value is zero and preserves previous behaviour in the tested regression cases.

## Provenance and checksums

`docs/source_snapshot_manifest.csv` records SHA-256 hashes for each file in both source snapshots. `SHA256SUMS.txt` at the bundle root records checksums for the complete handover bundle, excluding the checksum file itself.

The WebR runtime is third-party software retained for reproducibility; see `environment/THIRD_PARTY_NOTE.txt` before redistributing it outside the project.
