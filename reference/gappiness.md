# The gappiness a of shape

Gappiness refers to the internal gaps that defined by a set of
discretized cells or a point set that covers some cells in a
discretization structre.

## Usage

``` r
gappiness(x, s, ...)

# S4 method for class 'character,trigrid'
gappiness(x, s, exclude = NULL, full = FALSE)

# S4 method for class 'matrix,trigrid'
gappiness(
  x,
  s,
  long = NULL,
  lat = NULL,
  duplicates = FALSE,
  plot = FALSE,
  plot.args = NULL,
  full = FALSE,
  exclude = NULL
)
```

## Arguments

- x:

  The list of faces that are part of the shape.

- s:

  Spatial discretization strucutre (currently only `trigrid`-class
  icosahedral grid).

- ...:

  Arguments passed to class-specific methods.

- exclude:

  The list of faces that is to be excluded from the calculation

- full:

  `logical`, should only the estimate (`FALSE`) be returned, or
  additional data as well?(`TRUE`).

- long:

  `character`, column name of the longitudes.

- lat:

  `character`, column name of the latitudes.

- duplicates:

  `logical`, should identical coordinates be included in the calculation
  (default is `FALSE`)

- plot:

  Logical, should the result be plotted? Will plot over active plot (as
  in `add=TRUE`), if here is any.

- plot.args:

  List arguments passed to the plotting function: `plot`.

## Value

A proportion that gives the ratio of hole-cells compared to all occupied
cells.

## Examples

``` r
# create a grid
hex <- hexagrid(2, sf=TRUE)

# an example shape
shape <- paste0("F", c(4, 5, 11, 13, 15, 21, 24, 26, 32, 33, 34, 35, 36))

# the gappiness
gappiness(shape, hex)
#> [1] 0.2777778
```
