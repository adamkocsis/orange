## ranges <- function(x, long="long", lat="lat", method=c("centroid", "mgcd", "occupancy", "mst", "cha"), debug=FALSE,  ...){

## 	# the basic coordinates
## 	coordMat <- SimpleCoordinates(x, long=long, lat=lat)

## 	# pre-calculate the distance matrix for faster results
## 	if(any(c("mst","mgcd")%in%method)){
## 		dm <- icosa::arcdistmat(coordMat)
## 	}

## 	# Additional arguments
## 	addArgs <- list(...)

## 	# the results
## 	results <- list()

## 	if(any("centroid"%in%method)) results <- c(results,list(centroid=ranges_centroid_points(coordMat)))

## 	if(any("mgcd"%in%method)) results <- c(results, list(mgcd=ranges_mgcd(coordMat, dm=dm)))

## 	if(any("occupancy"%in%method)){

## 		# copy over the argumnets
## 		if(any(c("sf", "icosa")%in%names(addArgs))){
## 			newArgs <- addArgs[names(addArgs)%in% c("sf", "icosa")]
## 		}
## 		allArgs <- c(list(coordMat=coordMat), newArgs)

## 		# call with these arguments
## 		results <- c(results, list(occupancy=do.call("ranges_occupancy", allArgs))
## 	}

## 	if(any("mst"%in%method)) results <- c(results, list(mst=ranges_mstlength(coordMat))

## 	if(any("cha"%in%method)) results <- c(results, list(cha=ranges_chullsphere(coordMat))





## }




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
# Centroid
################################################################################

#' Calculate surface centroid of coordinates
#'
#' @param coordMat 2D numeric matrix with two columns: longitudes and latitudes.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{points}.
#' @return A longitude and latitude value.
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
#' cent <- ranges_centroid_points(coordMat, plot=TRUE)
#'
#' # secondary visualization
#' # points(cent[1], cent[2])
ranges_centroid_points <- function(coordMat, plot=FALSE, plot.args=NULL){
	result <- icosa::surfacecentroid(coordMat)
	if(plot){
		if(is.null(plot.args)) plot.args <- list(pch=4, cex=1.5, col="red")
		arguments <- c(list(x=result[1]), list(y=result[2]), plot.args)
		do.call("points", arguments)
	}
	return(result)
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
#' planarChull <- ranges_chullplane(coordMat, gcbound=FALSE, plot=TRUE)
#' # plot(planarChull$sf, col="#55000033", add=TRUE)
#'
#' planarChullGC <- ranges_chullplane(coordMat, gcbound=TRUE, plot=TRUE)
#' # plot(planarChullGC$sf, col="#00550033", add=TRUE)
ranges_chullplane <- function(coordMat, gcbound=FALSE, plot=FALSE, plot.args=NULL){

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
#' sphericalChull <- ranges_chullsphere(coordMat, plot=TRUE)
#'
#' # plot(sphericalChull$sf, col="#55000033", add=TRUE)
ranges_chullsphere<- function(coordMat, plot=FALSE, plot.args=NULL){

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


################################################################################
# Great Circle distances
################################################################################

#' Calculate range with the maximum great-circle method
#'
#' @param coordMat 2D numeric matrix with two columns: longitudes and latitudes.
#' @param dm If there is a pre-made distance matrix, it can be plugged in here.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{sf::plot}.
#' @return A list with an estimate an two indices the rows of the input matrix that represent the longest great circle (or one of them).
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
#' mgcd <- ranges_mgcd(coordMat, plot=TRUE)
#'
#' # single line visualization
#' # lines(coordMat[mgcd$index, ])
ranges_mgcd <- function(coordMat, dm=NULL, plot=FALSE, plot.args=NULL){

	# calculate the distance matrix
	if(is.null(dm)) dm <- icosa::arcdistmat(coordMat)

	# between which points is this observed - technically this can be between more than one pair!
	logMat <- dm==max(dm, na.rm=TRUE)
	where <- which(logMat, arr.ind=TRUE)
	for(i in 1:nrow(where)){
		where[i, ] <- sort(where[i, ])
	}
	where <- unique(where)

	# subset
	result <- list(
		estimate=max(dm),
		index = where
	)

	if(plot){
		if(is.null(plot.args)) plot.args <- list(col="#550000", lwd=3)
		arguments <- c(list(x=coordMat[where, ]), plot.args)
		do.call(icosa::arcs, arguments)
	}

	return(result)

}

#' Calculate range with the centroid-radius method
#'
#' @param coordMat 2D numeric matrix with two columns: longitudes and latitudes.
#' @param centroid The centroid given with a latitude-longitude coordinate pair.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{sf::plot}.
#' @param q The quantile that will be returned as an estimate
#' @return A list with an estimate an two indices the rows of the input matrix that represent the longest great circle (or one of them).
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
#' centDist <- ranges_centroid_radius(coordMat, plot=TRUE)
#'
#' lines(rbind(centDist$centroid,
#' 		coordMat[which(centDist$estimate==centDist$distances)[1],]),
#' 		col="blue", lwd=3)
#
ranges_centroid_radius <- function(coordMat, centroid=ranges_centroid_points(coordMat), plot=FALSE, plot.args=NULL, q=0.95){

	# the centroid
	centroidMat <- matrix(centroid, ncol=2, byrow=TRUE)

	# calculate the distance matrix
	dm <- as.numeric(icosa::arcdistmat(centroidMat, coordMat))
	names(dm) <- rownames(coordMat)

	# calculate the quantile
	estimate <- quantile(dm, q)

	# subset
	result <- list(
		estimate=estimate,
		distances = dm,
		centroid = centroid
	)

	# plotting
	if(plot){

		if(is.null(plot.args)) plot.args <- list(col="#88888866")
		for(i in 1:nrow(coordMat)){
			arguments <- c(list(x=rbind(centroid, coordMat[i, ])), plot.args)
			do.call(icosa::arcs, arguments)
		}

		# find the closest
		absDiff <- abs(dm-estimate)
		index <- which(absDiff==min(absDiff))[1]
		arcs(
			x=rbind(centroid, coordMat[index, ]),
			col="red", lwd=2
		)

		# visualize this as a small circle -# dependent on future icosa addition!!
		small<- smallcircles(x=centroidMat, r=estimate)
		arcs(small, col="gray", lty=2, lwd=2)

	}

	return(result)

}
################################################################################
# Latitudinal range
################################################################################

#' Calculate latitudinal range
#'
#' @param coordMat 2D numeric matrix with two columns: longitudes and latitudes.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{sf::plot}.
#' @return A list with an estimate and two latitudes: the minimum (southernmost) and maximum (northernmost) latitude.
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
#' latrange <- ranges_latrange(coordMat, plot=TRUE)
#'
#' # abline(h=latrange$range, lty=2)
ranges_latrange <- function(coordMat, plot=FALSE, plot.args=NULL){
	# the range
	ran <- range(coordMat[, "lat"], na.rm=TRUE)

	# the final object
	result <- list(
		estimate=diff(ran),
		range=ran
	)
	if(plot){
		if(is.null(plot.args)) plot.args <- list(col="#550000", lty=2)
		arguments <- c(list(h=ran), plot.args)
		do.call(abline, arguments)
	}

	# return
	return(result)

}

################################################################################
# Minimum spanning tree length
################################################################################

#' Calculate Minimum Spanning Tree length
#'
#' @param coordMat 2D numeric matrix with two columns: longitudes and latitudes.
#' @param dm If there is a pre-made distance matrix, it can be plugged in here.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param icosa An icosahedral grid to base the spanning tree on. If \code{NULL} then the points original coordinates will be maintained.
#' @param plot.args List arguments passed to the plotting function: \code{lines}.
#' @return A list with an estimate and the input for lines to show the MST.
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
#' mst <- ranges_mstlength(coordMat, plot=TRUE)
#' # lines(mst$show)
ranges_mstlength <- function(coordMat, dm=NULL, plot=FALSE, plot.args=NULL, icosa=NULL){
	if(!requireNamespace("vegan", quietly=TRUE)) stop("This function requires the 'vegan' extension package. ")

	if(is.null(colnames(coordMat))) colnames(coordMat) <- c("long","lat")

	# get rid of redundant entries
	coordMat <- unique(coordMat)

	# reduce to icosahedral gridpoints
	if(!is.null(icosa)){
		# look up the cells
		cells <- icosa::locate(icosa, coordMat)
		# look up the centers
		cents <- icosa::centers(icosa)
		# unique the cells
		un <- unique(cells)
		# the coordinates of the centers
		coordMat <- cents[un, ]
	}

	# calculate the distance matrix
	if(is.null(dm)) dm <- icosa::arcdistmat(coordMat)

	# calculte the minimum spanning tree
	stre <- vegan::spantree(as.dist(dm))

	# create an index input for the arcs function
	showThis <- matrix(ncol=2, nrow=0)
	index <- matrix(ncol=2, nrow=length(stre$kid))

	for(i in 1:length(stre$kid)){
		# the connected point index
		index[i, ] <- c(i+1, stre$kid[i])

		# point 1
		showThis <- rbind(showThis,
			matrix(c(
				coordMat[ i+1, "long"],
				coordMat[ i+1, "lat"],
				coordMat[stre$kid[i] ,"long"],
				coordMat[ stre$kid[i], "lat"], NA, NA
			), ncol=2, byrow=TRUE)
		)

	}

	if(!is.null(icosa)){
		result <- list(
			estimate=sum(stre$dist),
			centers=coordMat,
			index=index,
			show=showThis
		)
	}else{
		result <- list(
			estimate=sum(stre$dist),
			index=index,
			show=showThis
		)
	}
	if(plot){
		if(is.null(plot.args)) plot.args <- list(col="gray")
		arguments <- c(list(x=showThis), plot.args)

		if(!is.null(icosa)) plot(icosa, un, col="#FF000055", add=TRUE, border="white")
		do.call(icosa::arcs, arguments)
		# do.call(arcs, arguments)
	}

	return(result)
}

################################################################################
# Occupancy-related
################################################################################

#' Calculate ranges with the occupancy method
#'
#' @param coordMat 2D numeric matrix with two columns: longitudes and latitudes.
#' @param icosa An icosahedral grid from the package icosa.
#' @param sf An sf object (not yet!).
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{lines}.
#' @param alpha Minimum occupancy with \code{alpha} proportion of occurrences.
#' @return A list with an estimate and the input for lines to show the MST.
#' @export
#' @examples
#' # 1. Canvas
#' hex <- hexagrid(deg=3, sf=TRUE)
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
#' occ <- ranges_occupancy(coordMat, icosa=hex, plot=TRUE)
#'
#' # plot(hex, occ$cells, add=TRUE, col="green")
ranges_occupancy <- function(coordMat, sf=NULL, icosa=NULL, plot=FALSE, plot.args=NULL, alpha=1){


	if(!is.null(icosa)){
		result <- occupancy_icosa(coordMat, icosa, plot=plot, plot.args=plot.args, alpha=alpha)
	}
	if(!is.null(sf)){
		stop("Not yet!")
	}

	return(result)

}


occupancy_icosa <- function(x, icosa, plot=FALSE, plot.args=NULL, alpha=1){

	# the occupied cells by the points
	cells <- icosa::locate(icosa, x)

	if(alpha!=1){
		# tabulate the cells
		tabulatedCells <- table(cells)

		# in decreasinng order
		decreasingCells <- sort(tabulatedCells, decreasing=TRUE)

		# cumulated to get the total
		cumulated <- cumsum(decreasingCells)

		# the number of occurrences to consider
		nOccs <- nrow(x)*alpha

		# which are below the cutoff
		below <- which(cumulated < nOccs)

		# if there is at least something that needs to be omitted
		if(length(below)>0)	{
			first <- max(below)[1]

			if(first!=length(cumulated)) first <- first+1

		# otherwise
		}else{
			first <- length(decreasingCells)
		}
		occupCells <- names(decreasingCells)[1:first]

		# register these as well
		freq <- data.frame(frequency=as.numeric(decreasingCells), keep=FALSE)
		rownames(freq) <- names(decreasingCells)
		freq$keep[1:first] <- TRUE

		# the result object
		res <- list(
			estimate=length(unique(occupCells)),
			cells=occupCells,
			freq=freq
		)


		if(plot){
			if(is.null(plot.args)) plot.args <- list(col="#55000033")
			arguments <- c(list(x=icosa, y=res$cells, add=TRUE), plot.args)
			do.call(icosa::plot, arguments)
		}

	}else{
		# the occupied cell
		occupCells <- unique(cells)

		# the result object
		res <- list(
			estimate=length(unique(occupCells)),
			cells=occupCells
		)

		if(plot){
			if(is.null(plot.args)) plot.args <- list(col="#55000033")
			arguments <- c(list(x=icosa, y=res$cells, add=TRUE), plot.args)
			do.call(icosa::plot, arguments)
		}

	}

	return(res)

}




#' The gappiness of shape on an icosahedral grid
#'
#' @param x The list of faces that are part of the shape.
#' @param icosa The icosahedral grid.
#' @param ... Arguments passed to class-specific methods.
#' @param exclude The list of faces that is to be excluded from the calculation
#' @examples
#' # create a grid
#' hex <- hexagrid(2, sf=TRUE)
#'
#' # an example shape
#' shape <- paste0("F", c(4, 5, 11, 13, 15, 21, 24, 26, 32, 33, 34, 35, 36))
#'
#' # the gappiness
#' ranges_gappiness(shape, hex)
#' @rdname ranges_gappiness
#' @exportMethod ranges_gappiness
setGeneric(
	name="ranges_gappiness",
	def=function(x,icosa,...){
		standardGeneric("ranges_gappiness")
	}
)

#' @rdname ranges_gappiness
setMethod(
	"ranges_gappiness",
	signature=c(x="character", icosa="trigrid"),
		definition=function(x, icosa, exclude=NULL){

		# the holes of this patch
		theHoles<- holes(x, icosa)

		# gappiness is defined as number of hole cells divided by the holes and shape itself
		gap <- (length(theHoles)-length(exclude))/(length(theHoles) + length(unique(x))-length(exclude))

		# return both the metric and the way to plot it
		theGap <- list(estimate=gap, holes=theHoles)

		# return theappiness
		return(theGap)
	}
)


#' @rdname ranges_gappiness
setMethod(
	"ranges_gappiness",
	signature=c(x="matrix", icosa="trigrid"),
		definition=function(x, icosa, exclude=NULL){

			# get the list of faces occupied
			faceList <- locate(icosa, x)

			# calculate the gappiness based on the faces
			gap <- ranges_gappiness(faceList, icosa, exclude=exclude)

			# return the gappiness
			return(gap)
	}
)
