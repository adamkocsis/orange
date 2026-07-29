# Calculate latitudinal ranges

Calculate latitudinal ranges

## Usage

``` r
latrange(x, ...)

# S4 method for class 'matrix'
latrange(
  x,
  long = NULL,
  lat = NULL,
  q = 1,
  duplicates = FALSE,
  plot = FALSE,
  plot.args = NULL,
  full = FALSE
)

# S4 method for class 'data.frame'
latrange(
  x,
  tax = NULL,
  q = 1,
  long = "long",
  lat = "lat",
  duplicates = FALSE,
  plot = FALSE,
  plot.args = NULL,
  full = FALSE
)
```

## Arguments

- x:

  Either a 2D numeric `matrix` with two columns: longitudes and
  latitudes, a `data.frame` with the same information.

- ...:

  Arguments passed to class-specific methods.

- long:

  `character`, column name of the longitudes.

- lat:

  `character`, column name of the latitudes.

- q:

  `numeric`, a value between 0 and 1, the quantile.

- duplicates:

  `logical`, should identical coordinates be included in the calculation
  (default is `FALSE`)

- plot:

  Logical, should the result be plotted? Will plot over active plot (as
  in `add=TRUE`), if here is any.

- plot.args:

  List arguments passed to the plotting function: `points`.

- full:

  `logical`, should only the estimate (`FALSE`) be returned, or
  additional data as well?(`TRUE`).

- tax:

  `character`, used only in the `data.frame` method. Column name of
  groups (e.g. taxa) that allows the iteration of the method for
  multiple groups.

## Value

Either a single numeric or a list with an estimate and other
information. If iterated using \`tax\`, then either a named vector or
list of lists.

## Examples

``` r
# 1. Records
data(pinna)
# Subset to Pinna nobilis
nobilis <- pinna[pinna$species=="Pinna nobilis", ]
plot(nobilis[c("decimalLongitude", "decimalLatitude")], pch=16, col="#00BBAA66")

# Number of unique coordinate pairs
lr <- latrange(nobilis, long="decimalLongitude", lat="decimalLatitude", full=TRUE)

abline(h=lr$range, col="darkred", lty=3, lwd=4)
```
