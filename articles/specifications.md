# Specifications

## I. General

### Overview

#### About this document

This is the software requirement specification and design log of the
`orange` package. As there is a large combination of input variables for
every implemented approach, this document is intended to be coordinating
the development and and ensure a homogeneous interface and structure
through the package.

This document provides: - Define boundaries for the software package. -
Define rules for consistent interface development - Provide rationale
for implementation decisions - Forms the basis for the scientific
article describing the effort.

------------------------------------------------------------------------

#### Purpose of the `orange` package

The boundaries of this package can be drawn at the simple task: using
**spherical shape information** (e.g. taxon distribution data), it
allows the calculation of **numeric descriptors** (shape) that
characterize the sphes using their **intrinsic, geometric
characteristics**, and how they are modified by providing some
information that modify simple spherical geometry (e.g. a mask that
limits habitable space).

------------------------------------------------------------------------

#### Correctness

Although these calculations can be represented with ‘simple’
abstractions, the input data can vary extremely and the implementation
of these methods can be executed with varying awareness of the
difficulties imposed by spherical geometry, that may be erroneous in
some cases (such as convex hull calculations). Therefore it is a core
philosophy of the `orange` package, that implementations need to be
correct, regardless whether the input data are globally distributed or
confined to a simple region.

------------------------------------------------------------------------

#### Transparency

Since the metrics implemented here often reduce complex information to a
single number, users must get access to intermediate data objects to
accurately assess the correctness of what they are calculating - as well
as being able to use the intermediate objects in other, unexpected
calculations.

------------------------------------------------------------------------

### Input specifications

*Repetition goes to the index page.*

Although many metrics defined here have been available in simple
application scenarios, they can be extended almost trivially to other
combinations of input data, which may have better uses in some cases.
This section reviews the types of input data used as the basis of
geographic shape analysis (especially for taxon distribution data),
which are later used to systematically review and develop specific
methods.

------------------------------------------------------------------------

#### Primary input: distribution sample

The functions of metrics are S4 methods that are implemented for
multiple primary data input classes. The primary input is (following R
tradition) parametrized as `x`. The `x` parameter can be any of the
following input data structures:

- *coordmat*: A two-column tabular structure that inherits from
  `matrix`. The simplest input format, implies long-lat for the two
  columns.
- *coorddat*: A `data.frame` with geographic coordinates. The column of
  coordinates need to be specified, they default to `long="long"` and
  `lat="lat"`, but these can be specified. May contain subsets (implied:
  different taxa), for which most functions support iteration ( `tax` is
  the index for this).
- *sfpoints*: A multipoint object, either of classes `sf` or `sfc`.
- *sf-polygons*: A multipolygon class object `sf` or `sfc`.
- *icell*: Vector of cell identifiers of an icosahedral grid
  ([`icosa::hexagrid`](https://icosa-grid.github.io/R-icosa/reference/hexagrid-class.html)
  or
  [`icosa::trigrid`](https://icosa-grid.github.io/R-icosa/reference/trigrid-class.html)),
  a `character` vector. Note that these can only be used together with
  an `s=icosa` type icosahedral grid.
- *rast*: A `SpatRaster` class object defined in the `terra` package.

Note that the input types *coordmat* -\> *coorddat* -\> *sfpoints*
represent a natural order of promotion. If a method is described for a
lower rank type, it is automatically assumed that without the addition
of other parameters, the methods for the higher-rank types are
implemented as well. Arguments based on sf are not yet supported, they
will be included later.

In case of some input types, specific components can be specified with
the following parameters. - `long`: Column name indicating longitudes
(`character`) - `lat`: Column name indicating latitudes (`character`) -
`tax`: Column name indicating taxa, or other functional units that can
be used to split a dataset and the calculation of a metric to be
iterated with it.

------------------------------------------------------------------------

#### Secondary input: spatial structures

Many methods either rely or can incorporate additional spatial
structures to interact with the original distribution data. These
pre-specified structures are denoted with `s`. This parameter might have
a pivotal influence on how a metric behaves. For example `occupancy` is
defined in most cases as how much of `s` is occupied by `x`. In other
cases, `s` can be used as a discretization substitute for the original
data `x`, for example, by replacing the original point distribution data
with the centers of cells in an icosahedral grid. The different types of
`s` include:

- `missing`: If `s` is an parameter of a function, then there might be
  cases when spatial structuring is implied.
- `character`: Used when a `x=coorddat`, a column name that included
  implied structure, such as locality names.
- `icosa`: An icosahedral grid object defined with either
  [`icosa::trigrid`](https://icosa-grid.github.io/R-icosa/reference/trigrid-class.html)
  or the derived class
  [`icosa::hexagrid`](https://icosa-grid.github.io/R-icosa/reference/hexagrid-class.html).  
- `sf-polygons`: an sf/sfc class object that includes polygons.
- `rast`: A SpatRaster-type object.

------------------------------------------------------------------------

#### Additional data parameters

This is a list of additional parameters that may be necessary for some
metrics.

- `p`: Geographic coordinates (long, lat) of a single point given as a
  two-column matrix.
- `mask`: Another spatial structure that allows the exclusion of a
  certain area from the calculation of a metric.
- `q`: A single real number in the range of `[0,1]`, expressing
  confidence or critical interval, `numeric`.
- `prop`: In some cases this switch can be used to trigger proportional
  results.
- `dm`: A distance matrix that may replace a default distance matrix
  (e.g. a great circle distance matrix). Frequently used to refine the
  spatial environment, e.g. consider distances only on water, only on
  land, etc.
- `unique`: A logical switch indicating whether duplicate coordinates be
  omitted before the calculations, which can have effects on the results
  in some cases (e.g. centroid calculation). For the sake of efficiency
  and transparency, **this is set to to `TRUE` by default, which means
  that every coordinate pair is considered only once in the
  calculations**. The effects can be explored after a deliberate
  inclusion of duplicates.

------------------------------------------------------------------------

#### Behavioral modifiers

This list of parameters that change the behavior of functions. - `full`:
Boolean switch indicating whether a single estimate or
diagnostic/partial data should be returned from an object. - `plot`:
Boolean switch indicating whether the standard visualization of the
method should be executed or not. Defaults to `FALSE`. - `iter`: Natural
number indicating the number of iterations to be executed. Usually used
for metrics using Monte Carlo estimation.

------------------------------------------------------------------------

### Output specification

The default behavior of estimators is to return a single `numeric` value
for every distribution, which is the result of the argumentation
`full=FALSE`. This can be switched to `TRUE`, which will make range and
shape metrics output objects that inherit from the `orange` class. Every
`orange` class object must have at least two elements: - `$estimate`:
The numeric estimate calculated with the function. - *partials*:
Additional list elements that can be used to inspect the accuracy of the
method, as well as plotting additional ftod

------------------------------------------------------------------------

### Metrics and methods

This sections deals with some implementation details.

#### Function call stack notes

There are three tiers of functions: - 1. **Generics** that define the
names of the methods. - 2. **Methods** that invoke an internal logic
applicable for as set of input classes. - 3. **Internals** that
implement a logic.

**Methods** are defined using the rules of S4 method dispatch, based on
the primary arguments `x` and `s`,

------------------------------------------------------------------------

## II. Range/extent esimators

### Occupancy

The metrics are implemented with the `occupancy` generic function.
Occupancy is the number of elements in a spatial structure (`s`) that is
indicated to occupied by the distribution object (`x`). These spatial
structures can be implied (coordinates), named localities, or irregular
and regular (grids) sets of spherical polygons.

#### Overview of occupancy methods

| Status | Metric | Method | Args |
|----|----|----|----|
| ✅ | [No. of occupied coordinate pairs](#number-of-occupied-coordinate-pairs) | `x`: *coordmat* `s`:`missing` |  |
| ✅ | [No. of occupied coordinate pairs, iterated for taxa](#number-of-occupied-coordinate-pairs-iterated-for-taxa) | `x`: *coorddat* `s`:`missing` | `tax` `long`, `lat` |
| ✅ | [No. of occupied named localities](#number-of-occupied-named-localities) | `x`: *coorddat*, `s`:`char` |  |
| ✅ | [No. of occupied named localities, iterated for taxa](#number-of-occupied-named-localities-iterated-for-taxa) | `x`: *coorddat*, `s`:`char` | `tax` |
| ✅ | [No. of occupied grid cells in an icosahedral grid based on point data](#occupied-grid-cells-in-an-icosahedral-grid-based-on-point-data) | `x`: *coordmat*, `s`:`icosa` | `q=1`, `long`, `lat` |
| ✅ | [No. of occupied grid cells in an icosahedral grid based on point data, iterated for taxa](#occupied-grid-cells-in-an-icosahedral-grid-based-on-point-data-iterated-for-taxa) | `x`: *coorddat*, `s`:`icosa` | `tax`, `q=1`, `long`, `lat` |
| ❓ | [No. of occupied grid cells in an icosahedral grid, with a confidence cut-off](#occupied-grid-cells-in-an-icosahedral-grid-with-a-confidence-cut-off) | `x`: *coordmat*, `s`:`icosa` | `q<1`, `method`=“min” |
| ❓ | [No. of occupied grid cells in an icosahedral grid, with a confidence cut-off, iterated for taxa](#occupied-grid-cells-in-an-icosahedral-grid-with-a-confidence-cut-off-iterated-for-taxa) | `x`: *coorddat*, `s`:`icosa` | `tax`, `q<1`, `method`=“min” |
| ❌ | [No. of occupied polygons in an sf or sfc polygon/multipolygon object](#occupied-polygons-in-an-sf-or-sfc-polygonmultipolygon-object) | `x`: *coordmat*, `s`:`sf-polygons` | `q=1` |
| ❌ | [No. of occupied polygons in an sf or sfc polygon/multipolygon object, iterated for taxa](#occupied-polygons-in-an-sf-or-sfc-polygonmultipolygon-object-iterated-for-taxa) | `x`: *coorddat*, `s`:`sf-polygons` | `tax`, `q=1` |
| ❌ | [No. of occupied cells in a SpatRaster object](#occupied-cells-in-a-spatraster-object) | `x`: *rast* | `threshold`=0 |
| ❌ | [No. of occupied cells in an icosahedral grid based on a SpatRaster object](#occupied-cells-in-an-icosahedral-grid-based-on-a-spatraster-object) | `x`: *rast*, `s`:`icosa` | `threshold`=0 |
| ❌ | [No. of occupied polygons based on a SpatRaster object](#occupied-polygons-based-on-a-spatraster-object) | `x`: *rast*, `s`:`sf-polygons` | `threshold`=0 |

Occupancy methods

#### Number of occupied coordinate pairs

Returns the number of occupied coordinate pairs.

| Field        | Value                                                     |
|--------------|-----------------------------------------------------------|
| Usage        | `occupancy(x=coordmat)`                                   |
| Internal     | `occupancy_coords`                                        |
| Basic return | The number of occupied coordinate pairs.                  |
| Full return  | `$estimate`: The number of occupied coordinate pairs.     |
|              | `$occupied`: Unique coordiate pairs occupied by the data. |

[Back to table](#overview-of-occupancy-methods)

#### Number of occupied coordinate pairs iterated for taxa

Returns the number of occupied coordinate pairs - can be repeated for
every taxon.

| Field | Value |
|----|----|
| Usage | `occupancy(x=coorddat, tax)` |
| Internal | `occupancy_coords` |
| Basic return | The number of occupied coordinate pairs: a vector of occupancies for every taxon. |
| Full return | A list-array of taxa, every element referring to one taxon. The elements of these entries are: |
|  | `$estimate`: The number of occupied coordinate pairs. |
|  | `$occupied`: Unique coordiate pairs occupied by the taxon. |

[Back to table](#overview-of-occupancy-methods)

#### Number of occupied named localities

Returns the number of different locality entries.

| Field | Value |
|----|----|
| Usage | `occupancy(x=coorddat, s)` |
| Internal | `occupancy_coords` |
| Basic return | The number of localities. |
| Full return | `$estimate`: The number of localities occupie. |
|  | `$occupied`: Vector of unique locality entries occupied by the data. |

[Back to table](#overview-of-occupancy-methods)

#### Number of occupied named localities, iterated for taxa

Returns the number of different locality entries - can be repeated for
every taxon.

| Field | Value |
|----|----|
| Usage | `occupancy(x=coorddat, s, tax)` |
| Internal | `occupancy_coords` |
| Basic return | The number of occupied localities: a named vector of occupancies for every taxon. |
| Full return | A list-array of taxa, every element referring to one taxon. The elements of these entries are: |
|  | `$estimate`: The number of occupied localities. |
|  | `$occupied`: Vector of unique locality entries occupied by the taxon. |

[Back to table](#overview-of-occupancy-methods)

#### Occupied grid cells in an icosahedral grid based on point data

Returns the number of occupied grid cells in an icosahedral grid.

| Field | Value |
|----|----|
| Usage | `occupancy(x=coordmat, s=icosa)` |
| Internal | `occupancy_coords` |
| Basic return | The number of grid cells occupied. |
| Full return | `$estimate`: The number of grid cells occupied. |
|  | `$occupied`: Vector of unique grid cells identifiers occupied by the data. |
| Test case | `pinna` and `hexagrid(deg=5)` |

[Back to table](#overview-of-occupancy-methods)

#### Occupied grid cells in an icosahedral grid based on point data, iterated for taxa

Returns the number of occupied grid cells in an icosahedral grid.

| Field | Value |
|----|----|
| Usage | `occupancy(x=coorddat, s=icosa, tax)` |
| Internal | `occupancy_coords` |
| Basic return | The number of occupied grid cells: a named vector of occupancies for every taxon. |
| Full return | A list-array of taxa, every element referring to one taxon. The elements of these entries are: |
|  | `$estimate`: The number of grid cells occupied by the taxon. |
|  | `$occupied`: Vector of unique grid cells occupied by the taxon. |
| Test case | `pinna` and `hexagrid(deg=5)` |

[Back to table](#overview-of-occupancy-methods)

#### Occupied grid cells in an icosahedral grid, with a confidence cut-off

Returns the expected number of occupied grid cells that represent `q`\*
100% of the overal number of records.

Multiple implementations possible. These can be: - minimizing the range
`method="min"`: tabulating the record distribution on cells, sorting
them in decreasing order of record counts, and omitting the rarest cells
until only the `q` proportion of records are considered - the
expectation for the range, `method="mean", iter=100`. Records are
randomly drawn until `q` proportion of the overall cover is reached.

| Field | Value |
|----|----|
| Usage | `occupancy(x=coordmat, s=icosa, q=0.95)` |
| Internal | `occupancy_coords_icosa` |
| Return | `$estimate`: the number of grid cells occupied, `$cells`: The IDs of occupied cells |
| Test case | `pinna` and `hexagrid(deg=5)` |

[Back to table](#overview-of-occupancy-methods)

#### Occupied grid cells in an icosahedral grid, with a confidence cut-off, iterated for taxa

Returns the expected number of occupied grid cells that represent `q`\*
100% of the overal number of records.

Multiple implementations possible. These can be: - minimizing the range
`method="min"`: tabulating the record distribution on cells, sorting
them in decreasing order of record counts, and omitting the rarest cells
until only the `q` proportion of records are considered - the
expectation for the range, `method="mean", iter=100`. Records are
randomly drawn until `q` proportion of the overall cover is reached.

| Field | Value |
|----|----|
| Usage | `occupancy(x=coorddat, s=icosa, tax, q=0.95)` |
| Internal | `occupancy_coords_icosa` |
| Return | `$estimate`: the number of grid cells occupied, `$cells`: The IDs of occupied cells |
| Test case | `pinna` and `hexagrid(deg=5)` |

[Back to table](#overview-of-occupancy-methods)

#### Occupied polygons in an sf or sfc polygon/multipolygon object

The original definition for occupancy. Executes an intersection join,
and calculates the number of occupied polygons.

| Field | Value |
|----|----|
| Usage | `occupancy(x=coordmat, s=sfc)` |
| Internal | `occupancy_coords_sfc` |
| Return | `$estimate`: the number of grid cells occupied, `$ids`: The IDs of occupied polygons |
| Test case | `countries` and some terrestrial data |

The `sf`-method invokes the `sfc` method.

[Back to table](#overview-of-occupancy-methods)

#### Occupied polygons in an sf or sfc polygon/multipolygon object, iterated for taxa

The original definition for occupancy. Executes an intersection join,
and calculates the number of occupied polygons.

| Field | Value |
|----|----|
| Usage | `occupancy(x=coorddat, s=sfc, tax)` |
| Internal | `occupancy_coords_sfc` |
| Return | `$estimate`: the number of grid cells occupied, `$ids`: The IDs of occupied polygons |
| Test case | `countries` and some terrestrial data |

[Back to table](#overview-of-occupancy-methods)

#### Occupied cells in a SpatRaster object

Returs the number of cells where the value is at least the value given
as `threshold`.

| Field    | Value                                           |
|----------|-------------------------------------------------|
| Usage    | `occupancy(x=rast)`                             |
| Internal | `occupancy_coords_rast`                         |
| Return   | `$estimate`: the number of grid cells occupied, |

[Back to table](#overview-of-occupancy-methods)

#### Occupied cells in an icosahedral grid based on a SpatRaster object

Returns the number of cells where at least one raster cell has a value
that is at least the value given as `threshold`.

| Field | Value |
|----|----|
| Usage | `occupancy(x=rast, s=icosa)` |
| Internal | `occupancy_coords_rast` |
| Return | `$estimate`: the number of grid cells occupied, `$cells`: The IDs of occupied cells |

[Back to table](#overview-of-occupancy-methods)

#### Occupied polygons based on a SpatRaster object

Returns the number of polygon s where at least one raster cell has a
value that is at least the value given as `threshold`.

| Field | Value |
|----|----|
| Usage | `occupancy(x=sfc, s=icosa)` |
| Internal | `occupancy_coords_rast` |
| Return | `$estimate`: the number of grid cells occupied, `$ids`: The IDs of occupied polygons |

[Back to table](#overview-of-occupancy-methods)

------------------------------------------------------------------------

### Proportional occupancy

These metrics are also implemented with the `occupancy` function, and
are triggered when the `prop` general argument is set to a non-`NULL`
value. The implementation of these metrics are strongly tied to normal
occupany calculations, and their output matches that of normal
`occupancy` calls, except that these will return proportions instead of
counts as `estimate`s. Two different metrics are defined here, depending
on what a proportion relates to:

1.  Global proportional occupancies (`prop="global"`)  
2.  Relative proportional occupancies (`prop="relative"`)

#### Global proportional occupancies

*Global proportional occupancies* express how many out of the **globally
available** spatial slots are occupied by the distribution data `x`,
where global is given by the spatial structure itself. This can be the
total number of grid cells, the total number of polygons, etc. These
metrics are available only in those cases when `s` is a self-sufficient
spatial structure, such as `icosa`, `sf-polygons` or a `SpatRaster`
object. For cases when `s` is `missing` or `character`, setting
`prop="global"` triggers an error: since in these cases the number
globally available spatial slots is unknown to the function.

#### Relative proportional occupancies

In contrast, *relative proportional occupancies* express how many out of
the **sampled** spatial slots is occupied by the data, where what is
sampled is given by `x` itself. For instance, when a `tax` argument is
used to iterate the calculation of the metric across different taxa, the
relative proportional occupancies will return how many of the spatial
slots out of the totally sampled ones are occupied by a taxon. This can
be particularly useful when comparing data from multiple temporal
horizons, where the spatial coverage of overall sampling varies. When
this method is triggered with `tax=NULL`, the estimates default to `1.0`
(i.e. total occupancy).

------------------------------------------------------------------------

### Maximum Distance

The method is implemented with the `maxdist` function. For convenience,
the alias `mgcd` is also provided.

##### Overview

|  | Internal Name | Method | Args | What it does |
|----|----|----|----|----|
| ✅ | `maxdist_coords` | `x`: *coordmat*, `dm`:NULL, | `q=1` | The maximum great circle distance in a point set. |
| ✅ | `maxdist_coords` | `x`: *coorddat*, `dm`:NULL, | `q=1`, `long`, `lat` | The maximum great circle distance in a point set `data.frame`. |
| ✅ | `maxdist_coords` | `x`: *coordmat*, `dm`:*dm*, | `q=1` | The maximum great circle distance in a point set with a pre-defined distance matrix. |
| ✅ | `maxdist_coords` | `x`: *coorddat*, `dm`:*dm*, | `q=1`, `long`, `lat` | The maximum great circle distance in a point set `data.frame` with a pre-defined distance matrix. |
| ❌ | `maxdist_coords` | `x`: *coordmat*, `dm`:NULL | `q`=*q*, | Quantile great circle distance |
| ❌ | `maxdist_coords` | `x`: *coorddat*, `dm`:NULL | `q`=*q*, `long`, `lat` | Quantile great circle distance |
| ❌ | `maxdist_coords` | `x`: *coordmat*, `dm`:*dm*, | `q`=*q*, | Quantile great circle distance in a point set with a pre-defined distance matrix. |
| ❌ | `maxdist_coords` | `x`: *coorddat*, `dm`:*dm*, | `q`=*q*, `long`, `lat` | Quantile great circle distance in a point set `data.frame` with a pre-defined distance matrix. |
| ❌ | `maxdist_coords` | `x`: *coordmat*, `icosa` | `q=1` | The maximum great circle distance in a point set. |
| ❌ | `maxdist_coords` | `x`: *coorddat*, `icosa` | `q=1`, `long`, `lat` | The maximum great circle distance in a point set `data.frame`. |
| ❌ | `maxdist_coords` | `x`: *coordmat*, `icosa` | `q`=*q*, | Quantile great circle distance |
| ❌ | `maxdist_coords` | `x`: *coorddat*, `icosa` | `q`=*q*, `long`, `lat` | Quantile great circle distance |

##### Maximum Great Circle Distance with a given points

The basis for the method.

| Field | Value |
|----|----|
| Usage | `maxdist(coordmat)`, `mgcd(x=coordmat)` |
| Internal | `maxdist_coords` |
| Return | `$estimate`: the number of grid cells occupied, `$where`: Indices of point pairs. |
| Test case | `pinna` |

When the distance matrix parameter `dm` is not provided, the method will
omit duplicates and calculate the Great Circle Distance Matrix with
[`icosa::arcdistmat`](https://icosa-grid.github.io/R-icosa/reference/arcdistmat.html).

##### Maximum Distance with a given points specifying a distance matrix

This can be used in case topography/habitat-defined routes are more
important than great circle distances, or when the distance matrix would
have to recalculated.

| Field | Value |
|----|----|
| Usage | `maxdist(coords, dm=dm)` |
| Internal | `maxdist_coords` |
| Return | `$estimate`: the maximum distance calculated, `$index`: The row and column index of the relevant point pair. |
| Test case | `pinna` |

##### Quantile distance between points (Maximum distance with a confidence cutoff)

The same behavior as with the methods above, except that a
frequency-distribution of the distances is calculated first, and the `q`
quantile of the frequency distribution becomes the estimate. Works both
with and without a pre-specify distance matrix.

| Field     | Value                                 |
|-----------|---------------------------------------|
| Usage     | `maxdist(coords, dm=dm, q=0.95)`      |
| Internal  | `maxdist_coords`                      |
| Return    | `$estimate`: the distance calculated, |
| Test case | `pinna`                               |

The visualization of the method is based on a histogram.

------------------------------------------------------------------------

#### Latitudinal range - `latrange`

##### Overview

|  | Internal Name | Method | Args | What it does |
|----|----|----|----|----|
| ✅ | `range` | `x`: *coordmat* | `q=1`, `long`, `lat` | Latitudinal range of a long-lat point cloud. |
| ✅ | `range` | `x`: *coorddat*, | `q=1`, `tax` , `long`, `lat` | Latitudinal range of a long-lat point cloud. |
| ❌ | `range` | `x`: *coordmat*, | `q<1`, `long`, `lat` | Quantile latitudinal range of a long-lat point cloud. |
| ❌ | `range` | `x`: *coorddat*, | `q<1`, `tax` , `long`, `lat` | Quantile latitudinal range of a long-lat point cloud. |
| ❌ | `range` | `x`: *SpatRaster*, | `q=1`, `long`, `lat`, `threshold=0` | Latitudinal range of an `SpatRaster`. |
| ❌ | `range` | `x`: *sf-points*, | `q=1`, `long`, `lat` | Latitudinal range of an `sf` point cloud. |

##### Latitudinal range of a long-lat point cloud

Returns the latitudinal range of a single longitude-latitude point
cloud.

| Field | Value |
|----|----|
| Usage | `occupancy(x=coordmat)` |
| Internal | `range` |
| Return | `$estimate`: the number of occupied coordinate pairs, `$range`: minimum and maximum longiude |
| Test case | `pinna`, [`divDyn::corals`](https://rdrr.io/pkg/divDyn/man/corals.html) |

- `duplicates` have no effect when `q=1`

##### Latitudinal range of a long-lat point cloud (`data.frame`)

Returns the latitudinal range of a single or multiple longitude-latitude
point cloud(s).

| Field | Value |
|----|----|
| Usage | `occupancy(x=coorddf)` |
| Internal | `range` |
| Return | `$estimate`: the number of occupied coordinate pairs, `$range`: minimum and maximum longiude |
| Test case | `pinna`, [`divDyn::corals`](https://rdrr.io/pkg/divDyn/man/corals.html) |

- Use `tax` to iterate for multiple taxa.

#### Fixed radius - `fixrad`

The generic of this metric is the `fixrad` function. For convenience,
the alias \`\` is also provided. This metric is similar to the Maximum
Distance, but instead of calculating a full distance matrix, distances
compared to a single point are assessed. This point defaults to the
centroid.

##### Overview

| In | Internal Name | Main Input | What it does |
|----|----|----|----|
| ✅ | `fixrad_coords` | `x`=coords, `p`=centroid(coords) | The centroid radius of the point cloud. |
| ✅ | `fixrad_coords` | `x`=coords, `p`=centroid(coords), `q`=*q* | The centroid radius of the point cloud. |

##### Centroid radius

The method calculates every point’s distance to the centroid and
searches for the maximum.

| Field | Value |
|----|----|
| Usage | `fixrad(coords)` |
| Internal | `fixrad_coords` |
| Return | `$estimate`: the distance calculated, `$index`: which coordinate is the farthest from the centroid. |
| Test case | `pinna` |

##### Centroid quantile radius

The method calculates every point’s distance to the centroid, tabulates
the frequency distribution of the distances, and then returns distance
for a `q` quantile.

| Field | Value |
|----|----|
| Usage | `fixrad(coords, q=0.95)` |
| Internal | `fixrad_coords` |
| Return | `$estimate`: the distance calculated, `$index`: which coordinate is the farthest from the centroid. |
| Test case | `pinna` |

#### Zonal area

#### Minimum Spanning Tree (MST) Length

##### Overview

|  | Internal Name | Method | Args | What it does |
|----|----|----|----|----|
| ✅ | `mstlength_coords` | `x`: *coordmat*, `dm`:NULL, | `q=1` | The length of great circle arcs forming an MST in a point set. |
| ✅ | `mstlength_coords` | `x`: *coorddat*, `dm`:NULL, | `q=1`, `long`, `lat`, `tax` | The length of great circle arcs forming an MST in a point set `data.frame`, iteratated for every `tax`. |
| ✅ | `mstlength_coords` | `x`: *coordmat*, `dm`:*dm*, | `q=1` | The length of great circle arcs forming an MST in a point set with a pre-defined distance matrix. |
| ❌ | `mstlength_coords` | `x`: *coorddat*, `dm`:*dm*, | `q=1`, `long`, `lat`, `tax` | The length of great circle arcs forming an MST in a point set with a pre-defined distance matrix, iteratated for every `tax`. |
| ❌ | `mstlength_coords` | `x`: *coordmat*, `dm`:NULL | `q<1`, | The quantile length of great circle arcs forming an MST in a point set. |
| ❌ | `mstlength_coords` | `x`: *coorddat*, `dm`:NULL | `q<1`, `long`, `lat` , `tax` | The quantile length of great circle arcs forming an MST in a point set, iterated for every `tax`. |
| ❌ | `mstlength_coords` | `x`: *coordmat*, `dm`:*dm*, | `q<1`, | The quantile length of great circle arcs forming an MST in a point set, with a pre-defined distance matrix. |
| ❌ | `mstlength_coords` | `x`: *coorddat*, `dm`:*dm*, | `q<1`, `long`, `lat`, `tax` | The quantile length of great circle arcs forming an MST in a point set, iterated for every `tax`, and with a pre-defined distance matrix. |
| ❌ | `mstlength_icosa` | `x`: *coordmat*, `s`:*icosa* | `q=1` | The length of great circle arcs forming an MST of icosahedral grid cell centers. |
| ❌ | `mstlength_icosa` | `x`: *coorddat*, `s`:*icosa* | `q=1`, `long`, `lat`, `tax` | The length of great circle arcs forming an MST of icosahedral grid cell centers, iterated for every `tax`. |
| ❌ | `mstlength_icosa` | `x`: *coordmat*, `s`:*icosa* | `q<1`, | The quantile length of great circle arcs forming an MST of icosahedral grid cell centers. |
| ❌ | `mstlength_icosa` | `x`: *coorddat*, `s`:*icosa* | `q<1`, `long`, `lat`, `tax` | The quantile length of great circle arcs forming an MST of icosahedral grid cell centers, iterated for every `tax`. |

##### The length of great circle arcs forming an MST in a point set

Takes every coordinate pair and constructs an MST from the great circle
distances between the points. The sum length of this MST becomes the
estimate.

| Field | Value |
|----|----|
| Usage | `mstlength(coordmat)` |
| Internal | `mstlength_coords` |
| Return | `$estimate`: the total length of the MST (in km) |
|  | `$index`: the strucutre of the MST, which coordinate pair is connected with which |
|  | `$show`: Matrix of coordinates indicating what should be visualized |
| Test case | `pinna` |

- When the distance matrix parameter `dm` is not provided, the method
  will omit duplicates and calculate the Great Circle Distance Matrix
  with
  [`icosa::arcdistmat`](https://icosa-grid.github.io/R-icosa/reference/arcdistmat.html).
- `duplicates=TRUE` has no effect on the estimate as the distance
  between matching coordinate pairs is 0.

##### The length of great circle arcs forming an MST in a point set `data.frame`, iteratated for every `tax`.

The same method as above, but iterated for every entry in `tax`.

| Field | Value |
|----|----|
| Usage | `mstlength(coordmat)` |
| Internal | `mstlength_coords` |
| Return | `$estimate`: the total length of the MST (in km) |
|  | `$index`: the strucutre of the MST, which coordinate pair is connected with which |
|  | `$show`: Matrix of coordinates indicating what should be visualized |
| Test case | `pinna` |

##### The length of great circle arcs forming an MST in a point set with a pre-specified distance matrix

This can be used in case topography/habitat-defined routes are more
important than great circle distances, or when the distance matrix would
have to recalculated.

| Field | Value |
|----|----|
| Usage | `mstlength(coords, dm=dm)` |
| Internal | `mstlength_coords` |
| Return | `$estimate`: the total length of the MST (in km) |
|  | `$index`: the strucutre of the MST, which coordinate pair is connected with which |
|  | `$show`: Matrix of coordinates indicating what should be visualized |
| Test case | `pinna` |

##### The quantile length of great circle arcs forming an MST in a point set

- NOT YET DEFINED!
- Candidate method 1: use subsampling (`q` proportion of the points) and
  re-calculate the MST iteratively
- Candidate method 2: mst length \* q?

| Field     | Value                                 |
|-----------|---------------------------------------|
| Usage     | `maxdist(coords, dm=dm, q=0.95)`      |
| Internal  | `maxdist_coords`                      |
| Return    | `$estimate`: the distance calculated, |
| Test case | `pinna`                               |

------------------------------------------------------------------------

#### Hull Area

#### Circle Area

### Additional

#### Centroid

Sets of funtions that calculate the surface projection (geographic
coordinates) of centroid for a given geometry. \#### Overview

|  | Internal Name | Method | Args | What it does |
|----|----|----|----|----|
| ✅ | [`icosa::surfacecentroid`](https://icosa-grid.github.io/R-icosa/reference/surfacecentroid.html) | `x`: *coordmat*/*coordmat* | `long`, `lat` | Centroid coordinates of long-lat points. |
| ✅ | [`icosa::surfacecentroid`](https://icosa-grid.github.io/R-icosa/reference/surfacecentroid.html) | `x`: *coorddat* | `tax`, `long`, `lat` | Centroid coordinates of long-lat points. |

##### Centroid coordinates of long-lat points

Returns the number of occupied coordinate pairs.

| Field    | Value                                    |
|----------|------------------------------------------|
| Usage    | `occupancy(x=coordmat)`                  |
| Internal | `occupancy_coords`                       |
| Return   | The number of occupied coordinate pairs. |

- The `tax` column identifier `character` argument can be used to
  indicate iterated calculation for subsets.
- The `full` argument is not applicable here, there are no intermediate
  data.
- The `duplicates` argument modifies the calculation result. When set to
  `FALSE` (default), the duplicate coordinate entries will count as a
  single point in space. When set to `TRUE` the centroid will
  effectively be weigthed by the number of times the coordinate pair
  occurrs in the data, leading to slightly different results.

## Old

\[x\] \| `mgcd` \| cells, grid, q \| MGCD of cell centroids, or the
great circle distance of a quantile. \|  
\[x\] \| `mgcd` \| coords, q \| The maximum great circle distance of the
points \|  
\[x\] \| `cenrad` \| coords, q \| Centroid radius of points, maximum or
quantile distance. \|  
\[x\] \| `cenrad` \| cells, q \| Centroid radius of points, maximum or
quantile distance. \|  
\[x\] \| `latrange` \| coords \| The latitudinal range of a point cloud.
\|  
\[ \] \| `zoneaarea` \| coords \| The spherical surface area of a zone
defined by the latitudinal range of a point cloud. \|  
\[ \] \| `zoneaarea` \| coords/cells, grid \| The latitudinal range of a
point cloud. \|  
\[x\] \| `mstlength` \| coords \| The length of a minimum spanning tree
based on a set of points \|  
\[ \] \| `mstlength` \| cells, grid \| The length of a minimum spanning
tree based on the centers of cells. \|  
\[ \] \| `chull_cell` \| coords, grid \| The number of cells that a
convex hull occupies. \|  
\[ \] \| `chull_cell` \| cells, grid \| The number of cells that a
convex hull occupies. \|  
\[ \] \| `ahull_cell` \| cells, grid \| The number of cells covered by
an alpha hull on the surface of a sphere. \|  
\[ \] \| `ahull_cell` \| coords, grid \| The number of cells covered by
an alpha hull on the surface of a sphere. \|  
\[ \] \| `mincircle` \| coords, q \| The minimum small circle area,
radius and position that covers q proportion of the points \|  
\[ \] \| `minellipse` \| coords, q \| The covering spherical ellipse of
a set of points, shorter and longer axis . \|  
\[ \] \| `density` \| coords, q \| The densitiy of occurrence records \|

## Other esimators

| In | Function | Input | What it does |
|----|----|----|----|
| \[ \] | `ahull` | coords | The area of an alpha hull on the surface of a sphere defined by the points. |
| \[ \] | `chull` | coords | The spherical convex hull area of a set of points. |

## Shape esimators

| In | Function | Input | What it does |
|----|----|----|----|
| \[x\] | `patches` | cells, grid | The number of patches from a given set of cells. |
| \[ \] | `patches` | coords, grid | The number of patches from a given set of cells. |
| \[ \] | `holes` | coords, grid | The number of holes in the patches defined a given set of cells. |
| \[ \] | `holes` | cells, grid | The number of holes from a given set of cells. |
| \[ \] | `chull_cell_filling` | cells, grid | The proportion of `ranges_occupancy`/`ranges_chull_cell` , \[0,1\] |
| \[ \] | `chull_cell_filling` | coords, grid | The proportion of `ranges_occupancy`/`ranges_chull_cell` , \[0,1\] |
| \[ \] | `eccentricity` | coords, q | The eccentricity of a minimum ellipse |

| In    | Function   | Input  | What it does                       |
|-------|------------|--------|------------------------------------|
| \[x\] | `centroid` | coords | Latitude, longitude of point cloud |

## X. Geometry functions

### Small circles

This group of functions implement operations using small circles on a
sphere.
