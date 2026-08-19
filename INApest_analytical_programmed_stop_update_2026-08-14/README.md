# INApest analytical programmed-stop update

Date: 14 August 2026

This update extends the user-facing `INApestAnalytical()` screening function so finite `InfoPersistenceSteps` values are represented with explicit time-since-local-evidence states rather than a memoryless approximation.

Key points:

- finite programmed persistence takes priority over `InfoRetentionProb`, matching the simulation functions;
- scalar, node-specific and node x timestep persistence inputs are supported;
- where `InfoPersistenceSteps` is `NA`, the previous stochastic information-retention pathway is retained;
- binary INApest treats an informed extant infestation as continuing local evidence;
- SEAM/non-local information without new local evidence is kept distinct from a local clock reset;
- the deterministic clock is explicit, but shared node-level evidence and informed-but-uninvaded states remain approximations away from rarity in Meta/transition/MLU models;
- `InfoPersistenceSteps = NA` is exactly backward compatible with the preceding analytical release in the four-family regression suite.

The accompanying QMD and HTML have only targeted changes to the approved executive summary, information-state methods, limitations, outputs and future-development notes.

No production INApest simulation functions are changed by this analytical update.
