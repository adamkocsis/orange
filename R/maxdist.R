################################################################################
# Maximum-distance family of matrix Great Circle distances
################################################################################

#' Calculate ranges with the maximum distance method
#'
#' This family of metrics rely on the maximum distance within a point cloud
#'
#' This metrics includes maximum great circle distance and similar methods.
#'
#' @param x Either a 2D numeric \code{matrix} with two columns: longitudes and latitudes, a \code{data.frame} with the same information.
#' @param long \code{character}, column name of the longitudes.
#' @param lat \code{character}, column name of the latitudes.
#' @param tax \code{character}, used only in the \code{data.frame} method. Column name of groups (e.g. taxa) that allows the iteration of the method for multiple groups.
#' @param q \code{numeric}, a value between 0 and 1, the quantile. 
#' @param dm If there is a pre-made distance matrix, it can be plugged in here. If this is provided, the default coordinates will not be used.
#' @param full \code{logical}, should only the estimate (\code{FALSE}) be returned, or additional data as well?(\code{TRUE}).
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param icosa An icosahedral grid object the inherits from the \code{trigrid} class. Providing this argumnet reduces the point cloud to the centers of the grid cells.
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
	definition=function(x, dm=NULL, long=NULL, lat=NULL, duplicates=FALSE, plot=FALSE, plot.args=NULL, full=FALSE, q=1){


		# if the columns are given, make sure they are interpreted right!
		if(!is.null(long) | !is.null(lat)) x <- x[, c(long, lat), drop=FALSE]

		# create a copy
		if(full) xOrig <- x

		# if given
		if(!duplicates) x <- unique(x)

		# omit missing
		x <- x[!is.na(x[,1]) & !is.na(x[,2]),, drop=FALSE]
		
		# for the q
		if(q!=1 & !qTest) stop("Feature not yet implemented!")

		# if locality is given
		# calculate the distance matrix
		if(is.null(dm)){
			if(length(x)==0){
				if(full){
					retObj <- list(estimate=NA, index=NA)
					class(retObj) <- "orange"
					return(retObj)
				}else{
					return(NA) 
				}
			}
			dm <- icosa::arcdistmat(x)
		}else{
			# check whether the supplied distance matrix is the same as th pointset
		}
		# q not implemented yet
		result <- maxdist_coords(x=x, dm=dm, full=full, plot=plot, plot.args=plot.args)

		# translate the indices to the original
		if(full){
			# the first point's first match
			result$index <- c(
				which(x[result$index[1,1],1] == xOrig[,1] & x[result$index[1,1],2] == xOrig[,2])[1], 
				which(x[result$index[1,2],1] == xOrig[,1] & x[result$index[1,2],2] == xOrig[,2])[1]
			)
			# get rid of the names
			names(result$index) <- NULL
		}

		return(result)
	}
)

# uses the matrix method
#' @rdname maxdist 
setMethod(
	"maxdist",
	signature=c(x="data.frame", icosa="missing"),
	definition=function(x, long="long", lat="lat", tax=NULL, dm=NULL, duplicates=FALSE, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
		# the same as the matrix method
		if(is.null(tax)){
			x <- as.matrix(x[, c(long, lat)])
			result <- maxdist(x, long=long, lat=lat, duplicates=duplicates, dm=dm,
				q=q, plot=plot, plot.args=plot.args, full=full)
		}else{
			if(full) stop("Full output not yet implemented for tax-wise iteration.")
			if(is.null(dm)){
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
			# the distance matrix needs to be spliced!, and the 
			}else{
				stop("Not yet!")
			}
		}
		return(result)
	}
)


################################################################################
# Internals
################################################################################

# maximum distance method with distance matrix
maxdist_coords <- function(x, dm, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
	# between which points is this observed - technically this can be between more than one pair!
	estimate <- max(dm, na.rm=TRUE)
	# do this only if necessary
	if( full | plot ){

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
			arguments <- c(list(x=x[where, ]), plot.args)
			if(dev.cur()==1) plot(x, pch=16)
			do.call(icosa::arcs, arguments)
		}

		class(result) <- "orange"
	}

	if(!full) result <- estimate
	# return the final object
	 
	return(result)
}

