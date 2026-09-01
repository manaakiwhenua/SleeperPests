# Layer 3 literature sources

## Stage 3A — simple published Bayesian examples
Ramsey, D. S. L., Anderson, D. P. & Gormley, A. M. (2023). **Invasive species eradication: How do we declare success?** *Cambridge Prisms: Extinction*, 1, e4. DOI: 10.1017/ext.2023.1.

Published targets used:
- Prior 0.5 and SSe 0.9 gives PoA > 0.90 (exact Bayes value 0.9090909).
- Prior 0.8 and SSe 0.6 gives the same PoA > 0.90.
- The stopping-rule calculation for Prior 0.9 gives required SSe approximately 0.53 for target PoA 0.95 and approximately 0.91 for target PoA 0.99.

## Stages 3B–3C — repeated and spatial surveillance
Ward, D. F., Anderson, D. P. & Barron, M. C. (2016). **Using spatially explicit surveillance models to provide confidence in the eradication of an invasive ant.** *Scientific Reports*, 6, 34953. DOI: 10.1038/srep34953.

Published parameter-set-2 medians used:

| Survey | SSe | PoE |
|---|---:|---:|
| March 2013 | 0.149 | 0.312 |
| October 2013 | 0.713 | 0.611 |
| November 2013 | 0.734 | 0.855 |
| February 2014 | 0.736 | 0.957 |

Spatial detection inputs used:
- baited vial: g0=0.548, sigma=1.331, spacing 2 m; published same-cell probability about 0.70;
- visual search: g0=0.733, sigma=0.4, spacing 1 m; published same-cell probability about 0.75;
- sniffer dog: g0=0.750, sigma=1.65, spacing 2 m; published same-cell probability about 0.90.

The associated 2014 Landcare Research / Envirolink report additionally records the prior as a PERT distribution with mode 0.25 and range 0–0.75 and annual re-introduction risk as a PERT distribution with mode 0.01 and range 0–0.012.

## Stage 3D — broadscale and re-introduction logic
Anderson, D. P., Gormley, A. M., Ramsey, D. S. L. et al. (2017). **Bio-economic optimisation of surveillance to confirm broadscale eradications of invasive pests and diseases.** *Biological Invasions*, 19, 2869–2884. DOI: 10.1007/s10530-017-1490-5.

Published targets / equations used:
- MZ sensitivity: Se = 1 - (1 - PdAve × Prp)^Pu*.
- Posterior freedom after negative surveillance: Prior / [1 - Se(1-Prior)].
- Posterior becomes the next prior, optionally discounted for re-introduction.
- Their one-year 0.95 calibration uses Prior=0.70, Pd=0.90, Prp=0.98 and Pu*=1.
- Ten independent management zones each with freedom 0.95 imply about 0.40 probability that at least one still contains a residual population.

## Stage 3E — realistic multi-zone application
Anderson, D. P., Pepper, M. A., Travers, S. et al. (2022). **Confirming the broadscale eradication success of nutria (Myocastor coypus) from the Delmarva Peninsula, USA.** *Biological Invasions*, 24, 3509–3521. DOI: 10.1007/s10530-022-02855-x.

Published targets / parameters used:
- prior PoA approximately 0.01 in May 2015;
- mean PoA 0.75 at December 2020 (95% CI 0.30–0.93);
- forecast target PoA 0.95 by June 2022 under continued negative surveillance;
- public annual detection probability mean 0.03;
- 2015 occupied-cell design values: Blackwater 1, Eastern shore Virginia 4, Maryland 9, Delaware 11;
- annual increase: one occupied cell per year;
- 2022 values: 8, 11, 16, 18 respectively;
- 2022 mean zone sensitivities: 0.90, 0.33, 0.40, 0.43 respectively.

**Public-data boundary:** the article states that the USDA Wildlife Services surveillance data are not openly available. The Python model code is available, but the source surveillance data required for a complete independent 2015–2020 spatial rerun are not. Layer 3 therefore distinguishes exact/publicly reproducible checks from supporting consistency and explicitly marks full raw-data replication as REVIEW rather than pretending synthetic data are an external validation.
