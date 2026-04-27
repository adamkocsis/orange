# Change log of the R package 'orange'

# orange 0.1.0-6 - 2026-04-27

## Added

- The `kentsamples` data for testing and demonstration.
- The `circumcenter` function to calculate the small circle 
- The `smallcircles` function to calculate points on small circles

## Changed

- Omitted most likely irrelevant functions (based on an initial screening.
- Omitted `ranges_` prefix from range calculation functions. 
- Refactored functions and standardized interface.

## Known Missing 

- `plot=TRUE` for `data.frame` methods when a metric is iterated with `tax` 
- `maxdist` with `tax` entry and either `full=TRUE` or `dm!=NULL` 

## Initialization

- The package was forked from the 'biodome' project (v0.2.0-6) to efficiently focus on geographic range estimation. The corresponding functions were removed from that package.


* * *

