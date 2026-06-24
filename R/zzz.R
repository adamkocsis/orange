#' Orange
#'
#' Spherical Ranges and Parametrization of Geographic Shapes
#'
#' A collation of tools and metrics to characterize distribution data on the surface of a sphere or an ellipsoid. The primary group of these metrics is those that describe the extent of a distribution (geographic ranges). The calculation of geographic ranges can be executed using point coordinates data, vector polygons, as well as cells on a discretized sphere. Besides ensuring the use a geometrically correct implementations, the package offers the exploration of partial results for visual diagnostics.
#' This is still the pre-alpha version. Notes about found bugs and suggestions are more than welcome!
#'
#' @author Adam T. Kocsis (adam.t.kocsis@gmail.com)
#' @docType package
#' @name orange
"_PACKAGE"

#' @importFrom igraph V
#' @importFrom igraph random_walk
#' @importFrom igraph induced_subgraph
#' @importFrom utils flush.console
#' @importFrom utils head
#' @importFrom grDevices dev.cur
#' @importFrom graphics abline
#' @importFrom graphics points
#' @importFrom stats as.dist
#' @importFrom stats quantile
#' @import icosa
NULL
