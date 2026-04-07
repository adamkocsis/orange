################################################################################
# Maximum-distance family of matrix Great Circle distances
################################################################################

#' Calculate ranges with the maximum distance method
#'
#' This family of metrics include great circle distances and similar methods
#'
#' @param x 2D numeric matrix with two columns: longitudes and latitudes.
#' @param dm If there is a pre-made distance matrix, it can be plugged in here.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{sf::plot}.
#' @return A list with an estimate an two indices the rows of the input matrix that represent the longest great circle (or one of them).
#' @export
#' @rdname maxdist 
#' @examples
#' # 1. Canvas to plot on
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
#' mgcd <- mgcd(coordMat, plot=TRUE)
#'
#' # single line visualization
#' # lines(coordMat[mgcd$index, ])
setGeneric(
	name="maxdist",
	package="orange",
	def=function(x, icosa, ...){
		standardGeneric("maxdist")
	}

)

#' @rdname maxdist
#' @export
mgcd <- function(x, ...){
	maxdist(x, q=1, ...)
}

#' @rdname maxdist 
setMethod(
	"maxdist",
	signature=c(x="matrix", icosa="missing"),
	definition=function(x,tax, dm=NULL, long="long", lat="lat", duplicates=FALSE, plot=FALSE, plot.args=NULL, full=FALSE){
		# if locality is given
		if(!duplicates) x <- unique(x)

		# if the columns are given, make sure they are interpreted right!
		if(!is.null(long) | !is.null(lat)) x <- x[, c(long, lat), drop=FALSE]

		# for the q
		if(q!=1 & !qTest) stop("Feature not yet implemented!")

		# if locality is given
		# calculate the distance matrix
		if(is.null(dm)){
			dm <- icosa::arcdistmat(x)
		}else{
			# check whether the supplied distance matrix is the same as th pointselt
		}
		# q not implemented yet
		maxdist_coords(dm=dm, full=full, plot=plot, plot.args=plot.args)

		return(result)
	}
)

# uses the matrix method
#' @rdname maxdist 
setMethod(
	"occupancy",
	signature=c(x="data.frame", icosa="missing"),
	definition=function(x, long="long", lat="lat", tax=NULL, dm=NULL, duplicates=FALSE, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
		# the same as the matrix method
		if(is.null(tax)){
			x <- as.matrix(x[, c(long, lat)])
			result <- maxdist(x, long=long, lat=lat, duplicates=duplicates, dm=dm,
				q=q, plot=plot, plot.args=plot.args, full=full)
		}else{
			result <- tapply(
				INDEX=x[,tax],
				X=x[, c(long, lat)],
				FUN=function(a){

					# get rid of unnecessary recursion
					a <- as.matrix(a)
					result <- maxdist(a, long=long, lat=lat, duplicates=duplicates, dm=dm,
						q=q, plot=plot, plot.args=plot.args, full=full)
					}
				)
		}
		return(result)
	}
)


################################################################################
# Internals
################################################################################

# maximum distance method with distance matrix
maxdist_coords <- function(dm, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
	# between which points is this observed - technically this can be between more than one pair!
	estimate <- max(dm, na.rm=TRUE)
	# do this only if necessary
	if(full | plot ){

		# look up where is the estimate
		logMat <- dm==estimate

		# array indices 
		where <- which(logMat, arr.ind=TRUE)
		# make sure there is only one!
		for(i in 1:nrow(where)){
			where[i, ] <- sort(where[i, ])
		}
		where <- unique(where)

		# subset
		result <- list(
			estimate=estimate,
			index = where
		)

		# do a plot
		if(plot){
			if(is.null(plot.args)) plot.args <- list(col="#550000", lwd=3)
			arguments <- c(list(x=coordMat[where, ]), plot.args)
			do.call(icosa::arcs, arguments)
		}
	}

	if(!full) result <- estimate
	# return the final object
	 
	return(result)
}

