################################################################################
# Radius family of extent metrics
################################################################################

qTest <- FALSE
#' Calculate ranges with the radius-group of methods
#'
#' This family of metrics rely on estimating the distance of an extent from a fixed point.  
#'
#' This group of methods rely on calculating single distances. 
#'
#' @param x Either a 2D numeric \code{matrix} with two columns: longitudes and latitudes, a \code{data.frame} with the same information - or a character vector of cell identifiers.
#' @param s Structure to substitute the points, either missing (using coordinate pairs) or a \code{trigrid} (icosahedral grid from the package icosa).
#' @param p A single point of reference (longitude/latitude). If missing, the point of reference will be the centroid given by \code{centroid}.
#' @param long \code{character}, column name of the longitudes.
#' @param lat \code{character}, column name of the latitudes.
#' @param tax \code{character}, used only in the \code{data.frame} method. Column name of groups (e.g. taxa) that allows the iteration of the method for multiple groups.
#' @param duplicates \code{logical}, should identical coordinates be included in the calculation (default is \code{FALSE})
#' @param q \code{numeric}, a value between 0 and 1, the quantile. 
#' @param full \code{logical}, should only the estimate (\code{FALSE}) be returned, or additional data as well?(\code{TRUE}).
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{sf::plot}.
#' @return A list with an estimate an two indices the rows of the input matrix that represent the length of the tree (or one of them).
#' @export
#' @rdname radius 
#' @examples
#' # 1. Records
#' data(pinna)
#'
#' # 2. calculate and visualize
#' rad <- radius(pinna, long="decimallongitude", lat="decimallatitude", plot=TRUE, full=TRUE)
#'
setGeneric(
	name="radius",
	package="orange",
	def=function(x, s, ...){
		standardGeneric("radius")
	}

)

#' @rdname radius
setMethod(
	"radius",
	signature=c(x="matrix", s="missing"),
	definition=function(x, p=NULL, long=NULL, lat=NULL, duplicates=FALSE, plot=FALSE, plot.args=NULL, full=FALSE, q=1){

		# turn this off for now
		if(q!=1 & !qTest) stop("Feature not yet implemented!")

		# conserve the names if there are any
		nam <- rownames(x)

		# keep track of indices
		rownames(x) <- 1:nrow(x)

		# if the columns are given, make sure they are interpreted right!
		if(!is.null(long) | !is.null(lat)) x <- x[, c(long, lat), drop=FALSE]

		# if given
		if(!duplicates) x <- unique(x)

		# omit missing
		x <- x[!is.na(x[,1]) & !is.na(x[,2]),, drop=FALSE]
		
		# if the centroid, or fixed point is not 
		# calculate it here 
		if(is.null(p)){
			focus <- centroid(x)
		}		

		# the internal
		result <- qradius_coords(x=x, focus=focus, q=q, full=full, plot=plot, plot.args=plot.args)

		if(full){
			# translate the indices - ugh...
			result$index <- as.integer(rownames(x)[result$index])
		}

		return(result)
	}
)

# uses the matrix method
#' @rdname radius 
setMethod(
	"radius",
	signature=c(x="data.frame", s="missing"),
	definition=function(x, p=NULL, long="long", lat="lat", tax=NULL,  duplicates=FALSE, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
		# the same as the matrix method
		if(is.null(tax)){
			x <- as.matrix(x[, c(long, lat)])
			result <- radius(x, p=p, long=long, lat=lat, duplicates=duplicates,
				q=q, plot=plot, plot.args=plot.args, full=full)
		}else{
			if(full) stop("Full output not yet implemented for tax-wise iteration.")
			if(plot) warning("Multi-taxon plotting is not yet supported.")
			result <- tapply(
				INDEX=x[,tax],
				X=x[, c(long, lat)],
				FUN=function(a){

					# get rid of unnecessary recursion
					a <- as.matrix(a)
					result <- radius(a, p=p, long=long, lat=lat, duplicates=duplicates,
						q=q, plot=FALSE, plot.args=plot.args, full=full)
					}
				)
				resNames <- names(result)
				dim(result) <- NULL
				names(result) <- resNames
		}
		return(result)
	}
)


################################################################################
# Internals - logic and plotting included here
################################################################################
qradius_coords <- function(x, focus, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
	# NULL input -
	if(length(x)==0){
		if(full){
			return(
				list(
					focus=c(long=NA,lat=NA),
					estimate=NA,
					index=NA
				)
			)
		}else{
			return(NA)
		}
	}
	# calculate the distances between the point set and the focal point
	dists <- icosa::arcdistmat(points1=x, points2=matrix(focus, ncol=2))

	# between which points is this observed - technically this can be between more than one pair!
	estimate <- as.numeric(quantile(dists,q,  na.rm=TRUE))

	# do this only if necessary
	if( full | plot ){

		# look up where is the estimate
		logMat <- dists==estimate

		# array indices 
		where <- which(logMat)
		# make sure there is only one!
		where <- sort(where)
		where <- where[1]

		# subset
		result <- list(
			focus=focus, 
			estimate=estimate,
			index = where
		)

		# do a plot
		if(plot){
			# default plotting
			if(is.null(plot.args)) plot.args <- list(col="#550000", lwd=3)
			# the main arg. 
			arguments <- c(list(x=rbind(x[where, ], focus)), plot.args)
			# show the points too if no bg
			if(dev.cur()==1) plot(x, pch=16)
			do.call(icosa::arcs, arguments)

			# draft the small circle
			circ <- smallcircles(x=focus, r=estimate)
			arcs(circ, lty=2)
			
		}

		class(result) <- "orange"
	}

	if(!full) result <- estimate
	# return the final object
	 
	return(result)
}
