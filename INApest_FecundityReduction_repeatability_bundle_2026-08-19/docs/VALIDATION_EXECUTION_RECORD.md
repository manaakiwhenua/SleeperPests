# Validation execution record

Bundle assembly date: 19 August 2026.

## Executed during bundle assembly

1. **Original FecundityReduction validation workflow against current latest source defaults**
   - Runtime: WebR, R 4.6.0 (2026-04-24).
   - Result: completed successfully.
   - All six core archived CSV evidence files were reproduced byte-for-byte.

2. **Current ordinary Meta serial/parallel supplemental checks**
   - Scalar, node-vector and node-by-timestep `FecundityReduction` forms.
   - Result: all three current serial/parallel comparisons had maximum absolute difference 0 on the one-core path.
   - Backward custom-`LocalDynamics` contract: old custom function works when `FecundityReduction = 0`; active fecundity reduction produces the intended interface error if the new argument is not accepted.

3. **Current-source serial scientific reproduction**
   - Parameterisation smoke table.
   - Mortality × fecundity demonstration table.
   - Land-use demonstration table.
   - Result: all three tables matched the archived versions exactly.

4. **Portable exact WebR reproduction harness**
   - `environment/run_webr_exact_reproduction.py` was executed using the WebR runtime retained inside this bundle.
   - Result: all six archived evidence CSVs reproduced byte-for-byte.

5. **R script parse check**
   - All `.R` scripts under `scripts/` parsed successfully with the WebR R parser.

6. **HTML report QA**
   - The final HTML was rendered in headless Chromium and inspected at the executive-summary, multiple-land-use, validation and reproduction sections.

## Not executed in this environment

**True multi-process PSOCK validation.** WebR cannot create native R worker subprocesses. `scripts/05_native_psock_validation_current.R` is supplied for execution in a normal native-R environment with at least three detected cores. It forces the multi-worker branch and tests serial/parallel Monte Carlo agreement rather than exact RNG identity.
