# INApest PoA Layer 3 — external published validation methods

## Purpose
Layer 1 established that the software architecture runs correctly. Layer 2 established that the Bayesian and compatible-history mathematics agrees with independently calculable internal benchmarks. Layer 3 asks a different question: **does the PoA framework reproduce behaviour and numerical results reported independently in the published eradication-surveillance literature?**

Validation deliberately increases in complexity.

## Evidence classes
- **PASS** — an independently published numerical target or published equation is reproduced exactly or within the precision with which it was published.
- **SUPPORTED** — a deliberately simpler reconstruction using public summary inputs is quantitatively consistent with a more detailed published model. It is useful evidence, but is not described as an exact replication.
- **REVIEW** — a meaningful external target exists, but the public evidence is insufficient for independent rerunning. This is an evidence boundary, not a software failure.
- **FAIL** — a target that should be reproducible is missed.

## Stage 3A — simple published Bayes
With perfect specificity and no detections,

\[
PoA=\frac{Prior}{1-SSe(1-Prior)}.
\]

Ramsey et al. (2023) give two illustrative combinations that both exceed 0.90: Prior=0.5 with SSe=0.9, and Prior=0.8 with SSe=0.6. Both give 0.9090909 exactly. Required surveillance sensitivity for a target is

\[
SSe_{req}=\frac{PoA_{target}-Prior}{PoA_{target}(1-Prior)}.
\]

The INApest test evaluates the real `INApestPoACore()` against these targets using a synthetic two-state biological carrier. The target calculation is independent of INApest.

## Stage 3B — repeated published surveys
Ward et al. (2016) report four successive surveillance sensitivities and posterior probabilities of eradication for Argentine ants. The test reproduces every table row from the published SSe and the effective prior implied by that rounded row. This is an equation-level external reproduction, not a claim that the unpublished/random spatial realisations have been reconstructed.

A less circular check is also included. The published initial prior is PERT(mode=0.25, range=0–0.75). Under the standard beta-PERT parameterisation (lambda=4), its median combined with the published March SSe=0.149 gives a PoE close to the published March median 0.312. The small effective prior discounts between later surveys are also checked for consistency with the paper's stated very-low re-introduction assumption.

## Stage 3C — spatial detection
Ward et al. use a half-normal spatial detection function

\[
p(d)=g_0\exp\left(-\frac{d^2}{2\sigma^2}\right).
\]

For overlapping devices/search points along a path, the probability of at least one detection for a nest on the path is

\[
p_{line}=1-\prod_j[1-p(d_j)].
\]

The test independently reconstructs the same-cell/on-line probabilities for bait vials, visual search and sniffer dogs from the published `g0`, `sigma` and spacing values. The independently calculated spatial SSe is then passed to INApest's observation likelihood unchanged.

## Stage 3D — broadscale and re-introduction logic
Anderson et al. (2017) aggregate surveillance-unit detection into management-zone sensitivity:

\[
Se_j=1-(1-PdAve_j\,Prp_j)^{P_u^*}.
\]

They then apply the same Bayesian zero-detection update as above and carry the posterior into the next surveillance period, optionally discounting it for re-introduction risk. Tests reproduce their one-year 0.95 calibration and the published 10-management-zone compounding example, then exercise the published repeated-update and re-introduction structure through the INApest core.

## Stage 3E — realistic Delmarva nutria application
The nutria analysis is a much richer external case: four management zones, 100-m cells, active and public surveillance, spatial risk, increasing design prevalence, and annual Bayesian updating. Publicly reproducible checks include:

1. the published occupied-cell growth rule from 2015 to 2022;
2. the low-intensity public surveillance component, using
   \[
   SSe_{public}=1-(1-0.03)^k,
   \]
   where `k` is the number of occupied cells;
3. comparison of that simple component with published 2022 zone sensitivities;
4. the equivalent cumulative SSe required to move a prior of 0.01 to the published mean PoA of 0.75.

The public-only approximation is expected to be close for Virginia, Maryland and Delaware, which depended heavily on public surveillance. It should be much lower than the Blackwater total because active surveillance was concentrated in Blackwater.

### Important replication boundary
The article states that the underlying USDA surveillance data are not openly available. Therefore the exact spatial 2015–2020 analysis that produced mean PoA=0.75 cannot be independently rerun from public inputs. The validation suite marks this as **REVIEW** rather than substituting invented surveillance data.

## Interpretation
A clean Layer 3 result should not be summarized as “all external studies were exactly replicated.” The correct interpretation is:
- simple published Bayes: exact reproduction;
- repeated published surveys: table/equation reproduction plus prior-distribution consistency;
- spatial detection: direct kernel reproduction;
- broadscale logic: direct reproduction of published calculations and recurrence;
- realistic nutria case: several published components reproduced or supported, with full source-data replication explicitly unavailable.
