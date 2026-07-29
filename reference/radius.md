# Calculate ranges with the radius-group of methods

This family of metrics rely on estimating the distance of an extent from
a fixed point.

## Usage

``` r
radius(x, s, ...)

# S4 method for class 'matrix,missing'
radius(
  x,
  p = NULL,
  long = NULL,
  lat = NULL,
  duplicates = FALSE,
  plot = FALSE,
  plot.args = NULL,
  full = FALSE,
  q = 1
)

# S4 method for class 'data.frame,missing'
radius(
  x,
  p = NULL,
  long = "long",
  lat = "lat",
  tax = NULL,
  duplicates = FALSE,
  q = 1,
  plot = FALSE,
  plot.args = NULL,
  full = FALSE
)
```

## Arguments

- x:

  Either a 2D numeric `matrix` with two columns: longitudes and
  latitudes, a `data.frame` with the same information - or a character
  vector of cell identifiers.

- s:

  Structure to substitute the points, either missing (using coordinate
  pairs) or a `trigrid` (icosahedral grid from the package icosa).

- ...:

  Additional arguments passed to class-specific methods.

- p:

  A single point of reference (longitude/latitude). If missing, the
  point of reference will be the centroid given by `centroid`.

- long:

  `character`, column name of the longitudes.

- lat:

  `character`, column name of the latitudes.

- duplicates:

  `logical`, should identical coordinates be included in the calculation
  (default is `FALSE`)

- plot:

  Logical, should the result be plotted? Will plot over active plot (as
  in `add=TRUE`).

- plot.args:

  List arguments passed to the plotting function:
  [`sf::plot`](https://r-spatial.github.io/sf/reference/plot.html).

- full:

  `logical`, should only the estimate (`FALSE`) be returned, or
  additional data as well?(`TRUE`).

- q:

  `numeric`, a value between 0 and 1, the quantile.

- tax:

  `character`, used only in the `data.frame` method. Column name of
  groups (e.g. taxa) that allows the iteration of the method for
  multiple groups.

## Value

A list with an estimate an two indices the rows of the input matrix that
represent the length of the tree (or one of them).

## Details

This group of methods rely on calculating single distances.

## Examples

``` r
# 1. Records
data(pinna)
# Subset to Pinna nobilis
nobilis <- pinna[pinna$species=="Pinna nobilis", ]
plot(nobilis[c("decimalLongitude", "decimalLatitude")], pch=16, col="#00BBAA66")

# 2. calculate and visualize
rad <- radius(nobilis, long="decimalLongitude", lat="decimalLatitude", plot=TRUE, full=TRUE)
```
