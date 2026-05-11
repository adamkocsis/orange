


#' Function to process various coordinate input
#'
#' Note: this will be buried in the package namespace!
#'
#' @return A matrix of longitude and latitude coordinates
#' @param x A data.frame with occurrence records
#' @param long Column name of longitudes
#' @param lat Column name of latitudes
#' @return A two-column numeric matrix, longitudes and latitudes.
#' @export
SimpleCoordinates <- function(x, long="long", lat="lat"){
	object <- as.matrix(x[, c(long, lat)])
	colnames(object) <- c("long", "lat")
	return(object)
}

################################################################################
# Convex hull (planar)
################################################################################

#' Calculate planar convex hull
#'
#' Only included for the sake of comparison, it usage is in general not recommended.
#'
#' @param coordMat 2D numeric matrix with two columns: longitudes and latitudes.
#' @param gcbound Logical, should the bounds of the convex hull be great circles?
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{sf::plot}.
#' @return A list with an estimate an sf-object that represents the convex hull itself.
#' @export
#' @examples
#' # 1. Canvas
#' hex <- hexagrid(deg=5, sf=TRUE)
#' plot(hex, reset=FALSE, xlim=c(-15, 40), ylim=c(25, 63))
#'
#' # 2. Records
#' data(pinna)
#'
#' # just the coordinates
#' coordMat <- SimpleCoordinates(pinna, long="decimallongitude", lat="decimallatitude")
#' points(coordMat)
#'
#' # 3. calculate and visualize
#' planarChull <- chullplane(coordMat, gcbound=FALSE, plot=TRUE)
#' # plot(planarChull$sf, col="#55000033", add=TRUE)
#'
#' planarChullGC <- chullplane(coordMat, gcbound=TRUE, plot=TRUE)
#' # plot(planarChullGC$sf, col="#00550033", add=TRUE)
chullplane <- function(coordMat, gcbound=FALSE, plot=FALSE, plot.args=NULL){

	if(!requireNamespace("sf", quietly=TRUE)) stop("This function requires the 'sf' extension package. ")

	# the index of the convex hull
	hullIndex <- chull(coordMat)

	# has to be looped
	hullLast <- c(hullIndex, hullIndex[1])

	# the coordinates of the hull
	hullCoords <- coordMat[hullLast, ]

	#
	if(!gcbound){
		polyArc <- sf::st_geometry(sf::st_polygon(list(hullCoords)))
		sf::st_crs(polyArc) <- "WGS84"

	}else{

		# first coordinates
		arcCoords <- hullCoords[1,,drop=FALSE]

		for(i in 2:nrow(hullCoords)){
			# one arc between two points
			oneArc <- icosa::arcpoints(hullCoords[i-1,,drop=FALSE],hullCoords[i,,drop=FALSE], onlyNew=TRUE, output="polar", breaks=100)

			# the previous results, the arc and the point
			arcCoords <- rbind(arcCoords, oneArc, hullCoords[i,,drop=FALSE])
		}

		# create an sf version
		polyArc <- sf::st_geometry(sf::st_polygon(list(arcCoords)))
		sf::st_crs(polyArc) <- "WGS84"
	}

	res <- list(
		estimate= sf::st_area(polyArc),
		index=hullIndex,
		sf=polyArc
	)

	if(plot){
		if(is.null(plot.args)) plot.args <- list(col="#55000033", border="#550000")
		arguments <- c(list(x=polyArc, add=TRUE), plot.args)
		do.call("plot", arguments)
	}

	return(res)

}

################################################################################
# Convex hull (spherical)
################################################################################

#' Calculate spherical convex hull
#'
#' Wrapper around icosa::chullsphere. Still in development
#'
#' @param coordMat 2D numeric matrix with two columns: longitudes and latitudes.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{sf::plot}.
#' @return A list with an estimate an sf-object that represents the convex hull itself.
#' @export
#' @examples
#' # 1. Canvas
#' hex <- hexagrid(deg=5, sf=TRUE)
#' plot(hex, reset=FALSE, xlim=c(-15, 40), ylim=c(25, 63))
#'
#' # 2. Records
#' data(pinna)
#'
#' # just the coordinates
#' coordMat <- SimpleCoordinates(pinna, long="decimallongitude", lat="decimallatitude")
#' points(coordMat)
#'
#' # 3. calculate and visualize
#' sphericalChull <- chullsphere(coordMat, plot=TRUE)
#'
#' # plot(sphericalChull$sf, col="#55000033", add=TRUE)
chullsphere<- function(coordMat, plot=FALSE, plot.args=NULL){

	if(!requireNamespace("sf", quietly=TRUE)) stop("This function requires the 'sf' extension package. ")

	# the index of the convex hull
	hullIndex <- icosa::chullsphere(coordMat)

	# has to be looped
	hullLast <- c(hullIndex, hullIndex[1])

	# the coordinates of the hull
	hullCoords <- coordMat[hullLast, ]

	#
	# first coordinates
	arcCoords <- hullCoords[1,,drop=FALSE]

	for(i in 2:nrow(hullCoords)){
		# one arc between two points
		oneArc <- icosa::arcpoints(hullCoords[i-1,,drop=FALSE],hullCoords[i,,drop=FALSE], onlyNew=TRUE, output="polar", breaks=100)

		# the previous results, the arc and the point
		arcCoords <- rbind(arcCoords, oneArc, hullCoords[i,,drop=FALSE])
	}

	# create an sf version
	polyArc <- sf::st_geometry(sf::st_polygon(list(arcCoords)))
	sf::st_crs(polyArc) <- "WGS84"
	polyArc <- sf::st_wrap_dateline(polyArc)


	res <- list(
		estimate= sf::st_area(polyArc),
		index=hullIndex,
		sf=polyArc
	)

	if(plot){
		if(is.null(plot.args)) plot.args <- list(col="#55000033", border="#550000")
		arguments <- c(list(x=polyArc, add=TRUE), plot.args)
		do.call("plot", arguments)
	}

	return(res)

}





