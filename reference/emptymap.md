# Empty world map (Plate Carée)

Empty world map (Plate Carée)

## Usage

``` r
emptymap(
  xlim = c(-180, 180),
  ylim = c(-90, 90),
  grat = c(30, 30),
  col = "#DDDDFE",
  grat.col = "#FFFFFF",
  lwd = 2,
  thick = lwd * 2,
  border = grat.col
)
```

## Arguments

- xlim:

  Standard `xlim` argument of par in degress longitude.

- ylim:

  Standard `ylim` argument of par in degrees latitude.

- grat:

  Graticule resolution. Use `grat=NULL`, if you don't want any.

- col:

  Plot background, RGBA possible.

- grat.col:

  Graticule color RGBA possible.

- thick:

  The 0 coordinate circles are highlighted with thicker lines, given
  here

- border:

  Border color of rectange on the outside of the plot.

## Value

The function has no return value.

## Examples

``` r
data(kentsamples)
emptymap(l)
#> Error: object 'l' not found
points(kentsamples$central_s, pch=16, col=1)
#> Error in plot.xy(xy.coords(x, y), type = type, ...): plot.new has not been called yet
points(kentsamples$arctic_m, pch=16, col=2)
#> Error in plot.xy(xy.coords(x, y), type = type, ...): plot.new has not been called yet
points(kentsamples$antarctic_l, pch=16, col=3)
#> Error in plot.xy(xy.coords(x, y), type = type, ...): plot.new has not been called yet
points(kentsamples$dateline_l, pch=16, col=4)
#> Error in plot.xy(xy.coords(x, y), type = type, ...): plot.new has not been called yet
```
