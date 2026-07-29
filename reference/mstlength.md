# Calculate ranges with the minimum spanning tree length method

This family of metrics rely on constructing a minimum spanning tree from
the distances between the points, and use its length to describe
spreading.

## Usage

``` r
mstlength(x, s, ...)

# S4 method for class 'matrix,missing'
mstlength(
  x,
  dm = NULL,
  long = NULL,
  lat = NULL,
  duplicates = FALSE,
  plot = FALSE,
  plot.args = NULL,
  full = FALSE,
  q = 1
)

# S4 method for class 'data.frame,missing'
mstlength(
  x,
  long = "long",
  lat = "lat",
  tax = NULL,
  dm = NULL,
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
  latitudes, a `data.frame` with the same information.

- s:

  Structure to replace the points, either missing (coordinate pairs) or
  a `trigrid` (icosahedral grid from the package icosa).

- ...:

  Additional arguments passed to class-specific methods.

- dm:

  If there is a pre-made distance matrix, it can be plugged in here. If
  this is provided, the default coordinates will not be used.

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

This metrics includes maximum great circle distance and similar methods.

## Examples

``` r
# 1. Records
data(pinna)
# Subset to Pinna nobilis
nobilis <- pinna[pinna$species=="Pinna nobilis", ]
plot(nobilis[c("decimalLongitude", "decimalLatitude")], pch=16, col="#00BBAA66")

# 2. calculate and visualize
mst <- mstlength(nobilis, long="decimalLongitude", lat="decimalLatitude", plot=TRUE, full=TRUE)

mst$estimate
#> [1] 11163.76
```
