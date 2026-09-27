# INApest proof-of-freedom source handoff

**Checkpoint:** PoF biocontrol / multi-target v0.1  
**Frozen:** 24 September 2026  
**Native validation:** Windows R 4.4.1, 4/4 blocks PASS

## GitHub source

`R/INApestProofOfFreedom.R` is the canonical frozen source. It preserves the 1 September 2026 definitive pathogen PoF implementation and adds:

- biocontrol-agent proof of freedom;
- simultaneous pathogen + biocontrol inference;
- particle-wise joint proof of freedom;
- biocontrol surveillance likelihoods with stage selection; and
- state extraction for all eight supported biocontrol architectures.

The joint probability is computed from the same posterior particles and is not the product of marginal PoFs.

## Tests

`tests/` contains the exact four-block native validation suite used for the accepted checkpoint. The final run under R 4.4.1 passed all four blocks.

## Companion engine requirement

For Vertebrate Point PSOCK runs with active biocontrol, use the companion simulation-engine handoff dated 24 September 2026. Its `INApestVertebratePointParallel.R` contains the validated facade-only worker-export/output-preservation patch required by the PoF wrapper.
