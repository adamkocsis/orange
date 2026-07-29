
# orange <img src="man/figures/logo.png" align="right" />

[![](https://img.shields.io/badge/devel%20version-0.1.0--15-green.svg)](https://github.com/adamkocsis/orange)
[![](https://www.r-pkg.org/badges/version/orange?color=blue)](https://cran.r-project.org/package=orange)
[![](http://cranlogs.r-pkg.org/badges/grand-total/orange?color=yellow)](https://cran.r-project.org/package=orange)
[![CRAN
checks](https://badges.cranchecks.info/summary/orange.svg)](https://cran.r-project.org/web/checks/check_results_orange.html)

## Spherical Descriptors of Geographic Distributions

This package is a collation of tools and metrics to characterize
distribution data on the surface of a sphere. The primary group of these
metrics is those that describe the extent of a distribution,
i.e. geographic ranges/range sizes (depending on who you talk to). The
calculation of geographic ranges can be executed using point coordinate
data, cells on a discretized sphere, vector polygons and rasters.
Besides ensuring the use a geometrically correct implementations, the
package offers the exploration of partial results for visual diagnostics
and further post-processing.

------------------------------------------------------------------------

## Dependencies

The functions of the package are dependent on geographic discretization
structures. Due to its utility, the only hard dependency is `icosa`,
which implements near equal area gridding of the sphere. Additional
spatial structures include vector-spatial and raster objects implmented
in the `sf` and `terra` packages, respectively.

------------------------------------------------------------------------

## Example use

``` r
library(orange)

# background map
library(chronosphere)
ne <- fetch("NaturalEarth", verbose=FALSE)

# some example raw data (Pinna nobilis from OBIS)
data(pinna)
nobilis <- pinna[pinna$species=="Pinna nobilis", ]

# only coords
coords <- unique(nobilis[, c("decimalLongitude", "decimalLatitude")]) 
# rename to standard
colnames(coords) <- c("long", "lat") 

# omit issue localities
coords <- coords[-c(127,155), ]

# basic plot
plot(ne$geometry, reset=FALSE, xlim=c(-15, 40), ylim=c(30, 50),
    col="#ffe7ba", border=NA)
points(coords, pch=3, col="#5d7327", cex=2)

# Maximum great circle distance
mgcd <- mgcd(coords, plot=TRUE)

# Occupancy on a hexagonal grid
hex <- hexagrid(deg=2, sf=TRUE) # from package icosa
occ <- occupancy(coords, s=hex, plot=TRUE)

# Minimum spanning tree length
mst <- mstlength(coords, plot=TRUE, plot.args=list(lwd=3))


# centroid-circle hull
cc <- shull(coords, method="centroidcircle")
plot(cc, add=TRUE)
```

<img src="man/figures/range_example.png" alt="Equirectangular projection of the globe, showing some range calculation methods using the Pinna example dataset.">

The Maximum great circle distance:

``` r
# Maximum great circle distance
mgcd
```

    ## [1] 3452.84

Count of occupied cells:

``` r
# Occupancy count 
occ
```

    ## [1] 25

------------------------------------------------------------------------

This is a pre-alpha version, and like R, comes with absolutely no
warranty.
