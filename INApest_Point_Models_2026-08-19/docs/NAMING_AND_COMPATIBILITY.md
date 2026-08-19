# Naming and compatibility

`INApestMetaPoint()` is the canonical name for the point-based engine that inherits the management/environment architecture of `INApestMeta` while replacing predefined spatial nodes with explicit individual coordinates.

For compatibility with development scripts, sourcing `models/INApestMetaPoint.R` also defines:

```r
INApestPoint <- INApestMetaPoint
```

Returned objects inherit both `INApestMetaPoint` and `INApestPoint`, so older class checks continue to work.
