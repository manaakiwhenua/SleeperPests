# Pathogen detection -> information update

Added a common pathogen-surveillance option across the pathogen-capable INApest family.

Use in the common pathogen specification:

```r
Pathogen <- INApestPathogen(
  ...,
  DetectionProb = 0.2,
  DetectionTriggersInfo = TRUE
)
```

`DetectionTriggersInfo = FALSE` is the backward-compatible default.

Semantics:
- Binary INApest: pathogen detection probability applies when pathogen occupancy is present.
- Meta abundance models: `DetectionProb` is per infectious host; node detection is `1-(1-p)^I`.
- Multiple-land-use models: per-infectious-host detection is aggregated across all land-use classes to the node information state.
- Node transition-matrix models: per-infectious-host detection is aggregated across demographic stages to the node information state.
- Point models: detection is drawn per infectious point and recorded in `PathogenEvents$pathogen_detected`.
- Initial pathogen detection may seed information before timestep 1. Later pathogen detections occur after the current timestep's management decision and therefore affect management from the next timestep, matching the existing host-detection information timing.

New node-model output: `PathogenDetectedLargeOut.rds` where applicable.
