# Spherical hull geometries and hull area calculations

Calculation of hull geometries in spherical space

## Usage

``` r
shull(x, s, ...)

# S4 method for class 'matrix,missing'
shull(
  x,
  method,
  long = NULL,
  lat = NULL,
  duplicates = FALSE,
  plot = FALSE,
  plot.args = NULL,
  sphererad = authRadius
)

# S4 method for class 'data.frame,missing'
shull(x, long = "long", lat = "lat", ...)

# S4 method for class 'matrix,trigrid'
shull(
  x,
  s,
  method,
  long = NULL,
  lat = NULL,
  duplicates = FALSE,
  plot = FALSE,
  plot.args = NULL,
  sphererad = authRadius,
  drop = TRUE
)

shullarea(x, ...)

# Default S3 method
shullarea(
  x,
  method,
  metric = "area",
  full = FALSE,
  sphererad = authRadius,
  ...
)

# S3 method for class 'shull'
shullarea(x, metric = "area", s = NULL, full = FALSE, ...)
```

## Arguments

- x:

  Either a 2-column numeric matrix with two columns: longitudes and
  latitudes, or a `data.frame` with these columns.

- s:

  An external spatial structure to resolve the hull (e.g. a `hexagrid`
  object.)

- ...:

  Additional arguments passed to class-specific methods.

- method:

  `character` Hull construction method. Currently available:
  `"centroidcircle"`

- long:

  `character`, column name of the longitudes.

- lat:

  `character`, column name of the latitudes.

- duplicates:

  `logical`, should identical coordinates be included in the calculation
  (default is `FALSE`).

- plot:

  Logical, should the result be plotted? Will plot over active plot (as
  in `add=TRUE`), if there is any.

- plot.args:

  List arguments passed to the plotting function: `lines`.

- sphererad:

  `numeric` The radius of the sphere used in the calculations, defaults
  to the authalic radius of Earth in km (6370.997).

- drop:

  `logical` In case `s` is provided, should the returned object be a
  hull-class object, or only the identifiers of `s` (`drop=FALSE`)?

- metric:

  `character` What metric should be used for area of the hull? `area`
  returns the values in square of the unit of `sphererad` (defaults to
  square km), `prop` returns the area as a proportion of the sphere,
  `count` returns the area as the number of components it occupies in
  `s`, when applicable.

- full:

  Logical switch indicating whether only the estimate should be shown
  (`FALSE`), or the hull itself should be returned as well?.

## Value

A `hull`-class object, a vector of identifiers (if applicable and
`drop=FALSE`) or the area of te hull as a numeric.

## Details

Constructing and enclosing geometry around a distribution on the surface
of a sphere represent unique challenges. The `"centroidcircle"` method
is using the centroid of the distribution as the center of small circle
that defines a spherical cap. This is likely not the smallest possible
cap to cover the distribution and takes the heterogeneous density of the
distribution into account.

## Examples

``` r
# example data
data(kentsamples)
p <- kentsamples$dateline_m

# plotting
plot(NULL, NULL, xlim=c(-180, 180), ylim=c(-90,90), xlab="", ylab="")
points(p, pch=3, col="red")

# basic spherical hull definition (native)
cc <- shull(p, method="centroidcircle")
plot(cc, add=TRUE)

# spherical hull area calculation
shullarea(cc)
#> [1] 137466948

# basic hull definition (icosahedral grid)
hex <- icosa::hexagrid(spacing=8, sf=TRUE)
#> Selecting hexagrid with tessellation vector: c(3, 3).
#> Spacing: 7.686 degrees.
# dropping hull container
cc_hex <- shull(p, s=hex, method="centroidcircle")
#> Unfinished implementation! - Proper lookup of small circles required!
plot(hex, cc_hex, col="#0000DD33", border="#FFFFFF77", add=TRUE)



# with shull container
cc_hex2 <- shull(p, s=hex, method="centroidcircle", drop=FALSE)
#> Unfinished implementation! - Proper lookup of small circles required!
shullarea(cc_hex2, s=hex,  metric="prop")
#> [1] 0.2684729
```
