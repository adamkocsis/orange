#' Occurrences of Pinna nobilis from OBIS
#'
#' Occurrence records downloaded from OBIS on 2024-10-23
#'
#' This is an example occurrence record \code{data.frame}.
#'
#' @format A \code{data.frame} with 1466 observations and 6 variables:
#' 	\describe{
#' 		\item{\code{scientificname}}{Taxon name.}
#' 		\item{\code{basisofrecord}}{The basis of the record.}
#' 		\item{\code{date_year}}{Year.}
#' 		\item{\code{depth}}{Depth in meters.}
#' 		\item{\code{decimallongiude}}{Longitude in decimal degrees.}
#' 		\item{\code{decimallatitude}}{Latitude in decimal degrees.}
#' 	}
#'
#' @usage data(pinna)
"pinna"

#' Samples generated from various Kent distributions
#'
#' Geographic test and demonstration sample set
#'
#' A set of coordinates generated from Kent distributions with the sample number n = 1000. The samples are elements of a list, 16 in total with the combination of location and spread (i.e. the kappa parameter or the Kent distribution). There are four locations: \code{central} (0°N 0°E), \code{arctic} (85°N 5°E), \code{antarctic} (85°S 5°E) and \code{dateline} (170°E 5°N), and there are four sizes: small (\code{s}, kappa = 50), medium (\code{m}, kappa = 15), large (\code{l}, kappa = 5), and extra-large (\code{xl}, kappa = 2).
#'
#' @format A \code{liest} with 16 matrices (1000 x 2) with columns being equal to longitudes and latitudes.
#' 	\describe{
#' 		\item{\code{scientificname}}{Taxon name.}
#'		\item{\code{central_s}}{Central position, small spread. }
#'		\item{\code{central_m}}{Central position, medium spread. }
#'		\item{\code{central_l}}{Central position, large spread. }
#'		\item{\code{central_xl}}{Central position, extra large spread. }
#'		\item{\code{arctic_s}}{Arctic position, small spread. }
#'		\item{\code{arctic_m}}{Arctic position, medium spread. }
#'		\item{\code{arctic_l}}{Arctic position, large spread. }
#'		\item{\code{arctic_xl}}{Arctic position, extra large spread. }
#'		\item{\code{antarctic_s}}{Antarctic position, small spread. }
#'		\item{\code{antarctic_m}}{Antarctic position, medium spread. }
#'		\item{\code{antarctic_l}}{Antarctic position, large spread. }
#'		\item{\code{antarctic_xl}}{Antarctic position, extra large spread. }
#'		\item{\code{dateline_s}}{Dateline position, small spread. }
#'		\item{\code{dateline_m}}{Dateline position, medium spread. }
#'		\item{\code{dateline_l}}{Dateline position, large spread. }
#'		\item{\code{dateline_xl}}{Dateline position, extra large spread. }
#' 	}
#'
#' @usage data(kentsamples)
"kentsamples"
