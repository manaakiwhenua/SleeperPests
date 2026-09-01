# Source order for vertebrate PoA

The common PoA core remains `INApestPoA.R`. Vertebrate adapters/wrappers are additive in `INApestPoAVertebrate.R`.

## Vertebrate node PoA

```r
source("INApestVertebrateNode.R")
source("INApestPoA.R")
source("INApestPoAVertebrate.R")
```

## Vertebrate point PoA

```r
source("INApestPointTransitionMatrix.R")
source("INApestPoA.R")
source("INApestPoAVertebrate.R")
```

## Both vertebrate architectures

Source the node engine first, then the point-transition engine, then the PoA files. The shared vertebrate helper definitions used by the node engine (`.iv_home_range_node`, `.iv_node_control`, `.iv_node_birth`, `.iv_node_interaction`, `.iv_node_area_default`, `.iv_node_device_default`) were statically checked and are identical in the two supplied engine files.

```r
source("INApestVertebrateNode.R")
source("INApestPointTransitionMatrix.R")
source("INApestPoA.R")
source("INApestPoAVertebrate.R")
```
