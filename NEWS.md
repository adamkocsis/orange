# Change log of the R package 'orange'

# orange 0.1.0-11 - 2026-06-24

## Added

- The `kentsamples` data for testing and demonstration.
- The `sc_shape` function to calculate points on small circles 
- The `sc_center` function to calculate a small circle based on points on the sphere
- The `sc_in` function to calculate whether a bunch of points are inside a small circle or not (based on planes, rather than distances, about 10x faster).

## Changed

- Omitted irrelevant functions from `biodome` (based on an initial screening).
- Omitted `ranges_` prefix from range calculation functions. 
- The function `cenrad` was renamed to `radius`
- Refactored functions and standardized interface.

## Known Missing 

- Hull-based functions 
- `plot=TRUE` for `data.frame` methods when a metric is iterated with `tax` 
- `mstlength`, `maxdist` with `tax` entry and either `full=TRUE` or `dm!=NULL`, because the tracing of coordates in the original dataset is not yet implemented

* * *

## Initialization

- The package was forked from the 'biodome' project (v0.2.0-6) to efficiently focus on geographic range estimation. The corresponding functions were removed from that package.


* * *

