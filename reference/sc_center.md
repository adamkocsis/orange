# Identify the center of a small circle based on three points of a sphere

Identify the center of a small circle based on three points of a sphere

## Usage

``` r
sc_center(x, origin = c(0, 0, 0), output = "polar")
```

## Arguments

- x:

  A matrix of 3 points, either polar longitude and latitude coordinates,
  or XYZ Cartesian coordinates.

- origin:

  The origin of the sphere.

- output:

  If set to "polar" then the function will return the longitude and
  latitude of the small circle center.

## Value

A list with two elements, the center o the small circle and the surface
radius of the small circle.

## Examples

``` r
# generate 3 points on a sphere
set.seed(2)
ps <- icosa::rpsphere(3, output="polar")
small <- sc_center(x=ps)
circle<- sc_shape(x=small$center, r=small$r, output="polar")

plot(NULL, NULL, xlim=c(-180, 180), ylim=c(-90, 90))
points(ps, col="red", pch=3, cex=4)
points(small$center, col="blue", pch=16)
points(circle, col="green", pch=16)
```
