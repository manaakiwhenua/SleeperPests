# POA_validation_source

Frozen INApest proof-of-absence (PoA) source and validation bundle.
Freeze date: 2026-09-01.

## Folder layout

- `definitive_source/` — the frozen INApest simulation core source set plus the definitive PoA wrapper/core/adapters.
- `validation_code/` — corrected validation code for Layers 1, 2 and 3, plus independent benchmark/oracle inputs and the 27-Aug pre-patch baseline fixtures required for backward-compatibility checks.

The PoA wrapper in `definitive_source/INApestPoA.R` is the same frozen file used by the successful Layer 2 and Layer 3 validations (SHA-256 recorded in `FROZEN_SOURCE_SHA256.csv`).

## Run all validation layers on Windows

From this bundle root in PowerShell:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\validation_code\RUN_ALL_VALIDATION.ps1 -Rscript "C:\Program Files\R\R-4.4.1\bin\x64\Rscript.exe"
```

Each layer can also be run separately from its own folder.

## Validation evidence already obtained

See `VALIDATION_STATUS.md`.

## Important freeze rule

Treat `definitive_source/` as frozen. If any source file is edited, regenerate the SHA-256 manifest and rerun Layers 1-3 before calling the new state validated.
