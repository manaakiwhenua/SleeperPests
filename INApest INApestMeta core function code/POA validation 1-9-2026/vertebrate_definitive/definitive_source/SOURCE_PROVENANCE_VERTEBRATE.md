# Vertebrate PoA source provenance

- `INApestPoA.R` — copied byte-for-byte from the previously frozen `POA_validation_source` bundle; SHA-256 `603f93cbfb3331c98abfb808d76f9c77e4857a5e9ddd56e9f3e2d0721f88b5df`.
- `INApestPoAVertebrate.R` — new additive PoA adapter/wrapper companion created 1 September 2026. It does not redefine `INApestPoACore` or the existing six model-family adapters/wrappers.
- `INApestPointTransitionMatrix.R` — based on the previously frozen PoA source. The public `INApestVertebratePoint()` arguments and biological event ordering are unchanged; the only vertebrate-PoA addition is `ControlObservationHistory`, which retains the pre-control detection/kill opportunity and realised control-death state needed to evaluate the observation likelihood correctly.
- `INApestVertebrateNode.R` — based on the standalone vertebrate-node development source in the user Library (`file_000000001ad881fa839aa7761b4e40b3`, 21-Aug vertebrate source family). The PoA patch adds the canonical Background/InfoTriggered observation outputs, retains pre-control abundance/detectability, and adds suppressible output/progress side effects. The vertebrate Birth, HomeRange, Control and Interaction biological mechanisms are otherwise unchanged.

The current public GitHub root contains `INApestVertebratePoint()` inside `INApestPointTransitionMatrix.R`, but a standalone runtime `INApestVertebrateNode.R` was not found at the root during this extension work; the separately validated standalone node development source was therefore used for the node extension.
