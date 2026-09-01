# INApest PoA Layer 2 validation bundle

This bundle tests the **inference mathematics**, following the Layer-1 runtime validation of the PoA observation-engine patch.

## Run on the validated Windows R environment

Open PowerShell in this folder and run:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
.\RUN_LAYER2_VALIDATION.ps1 -Rscript "C:\Program Files\R\R-4.4.1\bin\x64\Rscript.exe"
```

Successful completion ends with:

```text
PASS: Layer 2 exact internal validation completed successfully.
```

Outputs are written to `validation_output`.

## What is independent of the PoA R code?

`expected/build_layer2_oracle.py` calculates the closed-form targets directly.  
`expected/independent_finite_state_check.py` independently verifies the trajectory fixtures for reinvasion, extinction/reinvasion, management feedback and unsupported-class propagation.

These Python checks were executed while the bundle was built. They do not call `INApestPoA.R`.

## Runtime status at packaging

- Independent closed-form oracle: PASS.
- Independent finite-state trajectory check: PASS.
- R source delimiter/contract preflight: PASS.
- Native R Layer-2 suite: requires execution in the supplied Windows R 4.4.1 environment.
