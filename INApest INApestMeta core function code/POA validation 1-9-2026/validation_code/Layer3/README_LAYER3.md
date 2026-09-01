# INApest PoA Layer 3 external validation

This bundle is the next validation layer after the completed Layer 1 software/runtime checks and Layer 2 exact internal Bayesian benchmarks.

## Run on Windows / native R
From PowerShell in this extracted folder:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\RUN_LAYER3_VALIDATION.ps1 -Rscript "C:\Program Files\R\R-4.4.1\bin\x64\Rscript.exe"
```

Outputs are written to `validation_output`:
- `layer3_test_results.csv`
- `layer3_numeric_comparisons.csv`
- `layer3_summary.txt`
- `layer3_console.log`
- `sessionInfo.txt`

## How to read the statuses
`PASS` means an externally published/reproducible target was met. `SUPPORTED` is deliberately weaker: a simpler public-data calculation is consistent with a richer published result. `REVIEW` records the nutria source-data limitation. The runner exits non-zero only for `FAIL`.

The bundle does **not** modify the validated PoA source. `source/INApestPoA.R` is a frozen copy of the Layer 2 source used for these external tests.
