# Identify whether a given set of points is in, on or out of a small circle

This method calculates the plane of the small circle and uses the normal
to evaluate whether points are in or outside the small circle.

## Usage

``` r
sc_in(x, center, r, origin = c(0, 0, 0))
```

## Arguments

- x:

  A set of points to be determined a matrix of longtidues and latitudes.

- center:

  The center of a small circle

- r:

  The surface radius of a small circle.

- origin:

  The origin of the sphere.

## Value

A numeric vector indicating whether the points are in (1), on (0) or
outside of the small circle.

## Examples

``` r
# set.seed(1)
# generate a random small circle
set.seed(1)
pol <- icosa::rpsphere(1000, output="polar")
# use the first as a small circle's center
center <- pol[1, , drop=FALSE]
radius <- 5000 # define radius in km

# draw a small circle
circle <- sc_shape(x=center, r=radius, breaks=200)

# visualize
plot(NULL, NULL, xlim=c(-180, 180), ylim=c(-90, 90))
points(pol, col="red", pch=3, cex=2)
arcs(circle, col="blue", lwd=2)

# find points in small circle (about 10x faster than distance evaluation)
inside <- sc_in(x=pol, center=center, r=radius)
points(pol[inside==1, ], col="green", pch=3, cex=2)

```
