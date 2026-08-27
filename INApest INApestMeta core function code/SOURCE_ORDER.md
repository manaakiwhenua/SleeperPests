# Recommended source combinations

The core files are designed to be sourced by model family. The validation scripts use isolated environments to avoid accidental masking among complete standalone files.

## Binary

```r
source("INApestPathogen.R")
source("INApest.R")
source("INApestParallel.R")
```

## Meta abundance

```r
source("INApestMeta.r")
source("INApestMetaParallel.r")
```

## Multiple land use

```r
source("INApestMetaMultipleLandUse.r")
source("INApestMetaParallelMultipleLandUse.r")
```

## Demographic transition matrix

```r
source("INApestPathogen.R")
source("INApestPathogenTransitionMatrix.R")
source("INApestMetaTransitionMatrix.r")
source("INApestMetaTransitionMatrixParallel.r")
```

## Point abundance

```r
source("INApestPathogen.R")
source("INApestPointPathogen.R")
source("INApestMetaPoint.R")
source("INApestMetaPointParallel.R")
```

## Point transition matrix

```r
source("INApestPathogen.R")
source("INApestPointPathogen.R")
source("INApestPointTransitionMatrix.R")
source("INApestPointTransitionMatrixParallel.R")
```

For validation, prefer the family-isolated loading pattern used in the supplied scripts rather than sourcing every complete standalone source file into the same global environment.
