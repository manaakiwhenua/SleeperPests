# INApest PoA Layer 2: exact internal behavioural validation

**Date:** 31 August 2026  
**Purpose:** test whether the PoA inference behaves correctly on deliberately small problems whose answers can be calculated independently.

Layer 1 established that the patched software runs and exposes the intended observation contract. Layer 2 asks a different question: **when the correct inputs are supplied, does the PoA calculation return the correct probability?**

The suite uses synthetic INApest outputs rather than complex ecological runs wherever possible. This isolates the Bayesian and compatible-history logic from spread, population, habitat and management implementation details.

## Common one-round calculation

Let

- \(A\) = pest is absent at the surveillance round;
- \(P\) = pest is present;
- \(\pi=P(A)\) = prior probability of absence;
- \(q=P(D=0\mid P)\) = probability of no detection if the pest is present;
- \(P(D=0\mid A)=1\) when false positives are not modelled.

After observing zero detections,

\[
P(A\mid D=0)=\frac{\pi}{\pi+(1-\pi)q}.
\]

This is the reference calculation used repeatedly below. Detection probability, surveillance-system sensitivity (SSe), and posterior PoA remain distinct quantities.

## Benchmarks

### L2-01 — One surveillance round

Prior PoA is 0.8 and detection probability conditional on presence is 0.8, so \(q=0.2\):

\[
P(A\mid D=0)=\frac{0.8}{0.8+0.2(0.2)}=0.9523809524.
\]

Expected Background SSe is 0.8.

### L2-02 — Repeated zero-detection rounds

With the same hidden state and independent surveillance each round, after \(k\) zero-detection rounds:

\[
P(A\mid D_1=0,\ldots,D_k=0)
=\frac{0.8}{0.8+0.2(0.2)^k}.
\]

Expected posteriors are 0.9523809524, 0.9900990099 and 0.9980039920 for rounds 1–3.

### L2-03 — Reinvasion between rounds

After round 1 the posterior is 0.9523809524. Let \(r=0.1\) be the probability that an absent system is reinvaded before round 2. With no extinction of already-present populations,

\[
P(A_2)=P(A_1\mid D_1=0)(1-r)=0.8571428571.
\]

After another zero-detection round with \(q=0.2\), expected PoA is 0.9677419355.

### L2-04 — Extinction plus reinvasion

Let \(r=0.2\) be reinvasion probability from the absent state and \(e=0.3\) be extinction probability from the present state. If \(a_1\) is round-1 posterior PoA, then

\[
P(A_2)=a_1(1-r)+(1-a_1)e.
\]

With prior PoA 0.7 and detection probability 0.6, the expected values are:

- round-1 posterior: 0.8536585366;
- round-2 prior: 0.7268292683;
- round-2 posterior: 0.8693115519.

### L2-05 — Background plus information-triggered surveillance

There are two present hidden states with equal prior weight. Both receive Background surveillance with detection probability 0.6. Only one was informed before surveillance and therefore receives targeted surveillance with detection probability 0.5.

For the two present states:

\[
q^{BG}=0.4,
\]

\[
q^{INFO}\in\{1,0.5\},
\]

and, conditional on the hidden state, the two pathways are independent:

\[
q^{ALL}=q^{BG}q^{INFO}.
\]

Averaging over present states gives Background SSe = 0.6, InfoTriggered SSe = 0.25, Combined SSe = 0.7, and posterior PoA = 0.7692307692 after zero detections from both sources.

### L2-06 — Detection changes later management

This is the key compatible-history benchmark. At round 1, a present trajectory that is actually detected is successfully managed and becomes absent by round 2. Another present trajectory is not detected and remains present.

The observed record is **zero detections**. Therefore the detected-and-managed trajectory is incompatible with the observation and cannot be used to represent future biology.

The correct compatible-history calculation gives:

- round-1 posterior PoA = \(2/3\);
- round-2 prior PoA = \(2/3\);
- round-2 posterior PoA = 0.8.

A deliberately incorrect calculation that likelihood-weights round 1 but then propagates the detected management history anyway gives round-2 prior PoA = \(5/6\) and posterior PoA = 0.9090909091. The test therefore discriminates sharply between correct and incorrect sequential propagation.

### L2-07 — High PoA without false certainty

Prior PoA is 0.999. Conditional on presence, abundance is 1 with probability

\[
\alpha=\frac{\sqrt{2}}{2},
\]

and abundance is 20 otherwise. Individual detection probability is 0.2. Therefore the exact mean no-detection probability conditional on presence is

\[
\bar q=\alpha(0.8)+(1-\alpha)(0.8)^{20}=0.5690622539.
\]

The exact posterior is

\[
P(A\mid D=0)=\frac{0.999}{0.999+0.001\bar q}=0.9994306924.
\]

Finite particle approximations are evaluated at 1,000, 10,000, 50,000 and 200,000 particles. The 200,000-particle result must be within \(5\times10^{-6}\) of the exact value and must remain below 1.

### L2-08 — Unsupported posterior class

If likelihood weighting leaves positive posterior mass on the present class but every realised present trajectory conflicts with the observation, there is no valid trajectory available to propagate that class. The correct behaviour is to **stop with an explicit error**, not silently convert inadequate simulation support into 100% PoA.

### L2-09 — Simulation-mode benchmark

A 100,000-particle realised-detection benchmark independently checks the older simulation-conditioning pathway. With prior PoA 0.8 and detection probability 0.8, the simulation-mode estimate must lie within 0.003 of the exact posterior 0.9523809524.

### L2-10 — Positive detection

For a true positive binary detection with no false-positive model,

\[
P(D>0\mid A)=0.
\]

Therefore posterior exact absence must be zero. With prior PoA 0.8 and detection probability 0.8, the marginal probability of the positive observation is \(0.2\times0.8=0.16\).

## Evidence produced

The R suite writes:

- `layer2_test_results.csv` — test-level PASS/FAIL;
- `layer2_numeric_comparisons.csv` — observed versus independently calculated targets;
- `layer2_high_poa_convergence.csv` — high-PoA particle sequence;
- `layer2_summary.txt` — concise runtime summary;
- `sessionInfo.txt` — R/runtime provenance.

The `expected/` directory contains a Python-built closed-form oracle and an independent finite-state trajectory check. Those files do not call the R PoA implementation.
