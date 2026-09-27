# INApest current source handoff — 28 September 2026

This bundle consolidates the accepted production source files identified for the SleeperPests main script folder.

## Included directly

- 19 accepted simulation-engine R source files from `INApest_SimEngine_GitHub_Source_20260924`.
- Frozen biocontrol/multi-target PoF source `INApestProofOfFreedom.R` from the 24-Sep PoF handoff.
- Provenance and native PoF validation status files.

## Canonical analytical source

The accepted core analytical source is pinned to immutable SleeperPests commit `6e4b9f032bec9001bf7241a9d59625894c776f9b`, Git blob `2a217a157e036be20b581517ceacfe066a686719`, SHA-256 `49238e62a9f4c99445072d08356ca33081c60b37f1a4b268b328794990d163c5`.

The raw project-library copy could not be materialised into this container. To avoid substituting a different file, run `GET_PINNED_INApestAnalytical.ps1` (Windows) or `GET_PINNED_INApestAnalytical.sh` to fetch the exact immutable file into `R/` and verify its SHA.

## Biocontrol analytical checkpoint

The latest report line reviewed is based on frozen parent `BCAN_v11_FROZEN_20260924.zip`, SHA-256 `6a87b84d72ce927d565518878a970367bb9e92b1e17d558862c12d5fcf165c3d`. This is an analytical/scientific checkpoint rather than a replacement engine source file.

## GitHub comparison

Current `main/INApest INApestMeta core function code/` is behind this accepted source state: several newer biocontrol, vertebrate-parallel and PoF files are absent, while some existing core engine files differ. The two pathogen source files were exact matches during the comparison.
