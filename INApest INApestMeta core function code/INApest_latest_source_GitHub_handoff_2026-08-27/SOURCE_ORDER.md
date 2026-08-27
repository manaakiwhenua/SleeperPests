# Suggested source order

For a clean R session, a conservative source order is:

1. `INApestPathogen.R`
2. `INApestPointPathogen.R`
3. `INApestPathogenTransitionMatrix.R`
4. Core serial engines (`INApest.R`, `INApestMeta.r`, MLU, TM, MetaPoint, PointTransitionMatrix)
5. Parallel wrappers
6. `INApestAnalytical.R`

`INApestAnalytical.R` is the consolidated analytical interface and includes the latest binary, Meta, MLU, Transition-Matrix, Point, vertebrate-host and vertebrate-pathogen analytical methods.
