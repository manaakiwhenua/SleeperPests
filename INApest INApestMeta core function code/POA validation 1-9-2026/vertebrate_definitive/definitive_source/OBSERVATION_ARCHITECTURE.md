# PoA observation architecture v2

## Separation of responsibilities

**Biological engine**

- simulates infestation/population state;
- simulates Background and InfoTriggered detection events;
- retains the information state before surveillance;
- retains the realised per-unit detection probabilities used in each simulated history;
- applies information/management feedback to future biological states.

**Observation adapter**

- maps each architecture to a common hidden-state representation;
- calculates `P(no detection | hidden state, realised detectability)`;
- keeps Background and InfoTriggered likelihoods separate.

**PoA inference core**

- applies prior PoA if supplied;
- likelihood-weights hidden states using observed surveillance;
- calculates stream-specific and combined surveillance-system sensitivity;
- propagates only realised histories compatible with the observed event history when later management depends on detections;
- reports posterior exact absence and optional functional freedom separately.

## Why pre-surveillance information is retained

InfoTriggered surveillance is conditional on the knowledge available when the targeted search was planned or initiated. Therefore the engine records information immediately before the host-surveillance round. A Background detection made during that round may update information for future response, but it cannot be used to claim that targeted search also occurred earlier in the same round.

## Why realised detection probabilities are retained

INApest permits uncertainty in detection probability through `DetectionSD`, and point/stage/land-use models can make detectability heterogeneous. The probability used by a particular simulated history therefore can differ among particles. Retaining that realised probability allows the observation layer to use the same surveillance model that generated the event and management history without putting Bayesian inference inside the engine.

## Combined no-detection likelihood

Where the two surveillance streams are conditionally independent given the simulated hidden state and their realised detection probabilities:

`q_combined = q_background * q_info_triggered`

This is an observation-model statement, not an assumption that the two programmes are operationally independent. A custom observation model can be supplied where dependence must be represented explicitly.
