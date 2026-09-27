# INApest simulation-engine source handoff

**Source status:** latest accepted source set used by the frozen PoF biocontrol / multi-target v0.1 checkpoint  
**Date:** 24 September 2026

This directory contains the complete 19-file simulation-engine source set used in the accepted multi-target PoF validation.

- 18 files are byte-identical to the accepted validated biocontrol source line.
- `INApestVertebratePointParallel.R` contains the narrow, validated PoF/PSOCK facade patch.
- The patch only exports point-biocontrol helper symbols to workers and preserves/recombines `BiocontrolHistory` and `BiocontrolPointEvents`.
- No pest, pathogen, biocontrol, movement, attack, recruitment, or management dynamics were changed by that patch.

`SOURCE_PROVENANCE.csv` records the canonical source SHA-256 values and identifies the one patched facade.

The later spatial-search / explicit-point-agent v0.1.4 work is intentionally not included here because it remained a candidate rather than an accepted frozen engine release.
