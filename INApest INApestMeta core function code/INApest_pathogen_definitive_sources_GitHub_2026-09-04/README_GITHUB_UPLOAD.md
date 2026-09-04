# GitHub handoff: definitive commented INApest pathogen simulation engines

Date: 2026-09-04

This is a clean source handoff for the `manaakiwhenua/SleeperPests` repository. The mirrored repository path is:

`github_upload/INApest INApestMeta core function code/`

The live repository uses this core directory for the main INApest simulation-engine source files. This handoff contains **18 destination source files** selected from the fully verified architecture-specific commented validation bundle.

## What to upload

Copy the contents of:

`github_upload/INApest INApestMeta core function code/`

into the repository directory:

`INApest INApestMeta core function code/`

and replace the corresponding source files. Do not copy the CSV/Markdown provenance files into the core source directory unless you want them versioned there as documentation.

For a local clone on Windows, `APPLY_TO_LOCAL_REPO.ps1` performs the overlay after you provide the repository root.

## Important INApestPathogen.R selection

The architecture-validation archive contains two executable snapshots named `INApestPathogen.R`. They differ only in an abundance-model initialisation line that is outside the Binary pathogen branch. The Binary-local snapshot retains the pre-repair `pmin()` clipping behaviour; validation of Meta exposed that silent-clipping defect and the later generic helper removed it.

The repository handoff therefore uses the **corrected generic `tm/INApestPathogen.R` snapshot** as the single root `INApestPathogen.R`. Its Binary branch is unchanged relative to the Binary-local helper, while its abundance branch contains the validated Meta/MLU/TM repair. The architecture-specific Binary-local copy remains preserved in the full definitive commented validation bundle, but it is not the appropriate root repository helper.

## Validation and freeze status

All 18 upload files originate from the pinned commented source set. Their executable code is inherited from the frozen validation programme: 1520 / 1520 architecture-level assertions PASS, with a separate 77 / 77 Vertebrate Node development-gate record.

`SOURCE_MAPPING.csv` records the exact source selection and both the frozen executable and commented SHA-256 hashes. `COMMENT_ONLY_VERIFICATION.csv` retains the complete 27-file architecture-specific comment-only verification record.

## Suggested commit message

`Pin definitive commented pathogen simulation-engine sources (2026-09-04)`
