# Suggested source order

For ordinary non-pathogen PoA work:

1. Source the biological engine(s) required for the analysis.
2. Source `INApestPoA.R`.

For pathogen-enabled simulations, source the relevant unchanged pathogen helper before running the biological model:

- `INApestPathogen.R`
- `INApestPathogenTransitionMatrix.R` for node transition-matrix pathogen work
- `INApestPointPathogen.R` for point pathogen work

Parallel wrappers require their matching serial engine to be sourced first where the wrapper delegates to it (binary and point wrappers).

The PoA companion does not require pathogen helpers unless the underlying simulation does.
