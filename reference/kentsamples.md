# Samples generated from various Kent distributions

Geographic test and demonstration sample set

## Usage

``` r
data(kentsamples)
```

## Format

A `liest` with 16 matrices (1000 x 2) with columns being equal to
longitudes and latitudes.

- `scientificname`:

  Taxon name.

- `central_s`:

  Central position, small spread.

- `central_m`:

  Central position, medium spread.

- `central_l`:

  Central position, large spread.

- `central_xl`:

  Central position, extra large spread.

- `arctic_s`:

  Arctic position, small spread.

- `arctic_m`:

  Arctic position, medium spread.

- `arctic_l`:

  Arctic position, large spread.

- `arctic_xl`:

  Arctic position, extra large spread.

- `antarctic_s`:

  Antarctic position, small spread.

- `antarctic_m`:

  Antarctic position, medium spread.

- `antarctic_l`:

  Antarctic position, large spread.

- `antarctic_xl`:

  Antarctic position, extra large spread.

- `dateline_s`:

  Dateline position, small spread.

- `dateline_m`:

  Dateline position, medium spread.

- `dateline_l`:

  Dateline position, large spread.

- `dateline_xl`:

  Dateline position, extra large spread.

## Details

A set of coordinates generated from Kent distributions with the sample
number n = 1000. The samples are elements of a list, 16 in total with
the combination of location and spread (i.e. the kappa parameter or the
Kent distribution). There are four locations: `central` (0°N 0°E),
`arctic` (85°N 5°E), `antarctic` (85°S 5°E) and `dateline` (170°E 5°N),
and there are four sizes: small (`s`, kappa = 50), medium (`m`, kappa =
15), large (`l`, kappa = 5), and extra-large (`xl`, kappa = 2).
