# Frozen validation status

Native validation environment reported by the completed runs:

- R version: R 4.4.1 (2024-06-14 ucrt)
- Platform: x86_64-w64-mingw32
- Windows / PowerShell execution

## Layer 1 — structural and runtime validation

**PASS.** The final corrected Layer-1 runner completed all supplied native-R validation scripts successfully.

This covers the common observation architecture, engine-level Background/InfoTriggered bookkeeping across the six biological architectures and parallel implementations, causal timing, realised detection probabilities, singleton-dimension regressions, and seeded default-path backward compatibility against the 27-Aug baseline.

The common observation architecture block was explicitly reported as **29/29 PASS** before the final full Layer-1 run.

## Layer 2 — exact internal behavioural validation

**10/10 PASS; 0 FAIL.** Maximum asserted numerical error: **0.0001269672**.

The suite covers closed-form Bayes, repeated negative surveillance, reinvasion, extinction plus reinvasion, two surveillance pathways, compatible-history management feedback, high-PoA convergence, unsupported-class safeguards, simulation-mode agreement and positive detection.

## Layer 3 — published external validation

**21 PASS; 3 SUPPORTED; 1 REVIEW; 0 FAIL.**

- PASS = exact or published-rounding reproduction.
- SUPPORTED = simplified public-data reconstruction behaves consistently with a more detailed published model, but is not an exact raw-data replication.
- REVIEW = the source-data boundary prevents full independent replication.

The REVIEW classification is the nutria full spatial trajectory, for which the underlying USDA surveillance data are not openly available. No synthetic substitute is presented as validation.

## Frozen interpretation

The source state in `definitive_source/` is the validated PoA implementation corresponding to these results. Any later source change creates a new validation state and should not inherit these results without rerunning the validation ladder.
