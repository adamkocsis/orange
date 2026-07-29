#' The gappiness a of shape
#'
#' Proportion of gap cells
#' 
#' Gappiness refers to proporion of area that are covered by the internal gaps that defined by a set of discretized cells or a point set that covers some cells in a discretization structre.
#'
#' @param x The list of faces that are part of the shape.
#' @param s Spatial discretization strucutre (currently only \code{trigrid}-class icosahedral grid).
#' @param long \code{character}, column name of the longitudes.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}), if here is any.
#' @param duplicates \code{logical}, should identical coordinates be included in the calculation (default is \code{FALSE})
#' @param plot.args List arguments passed to the plotting function: \code{plot}.
#' @param lat \code{character}, column name of the latitudes.
#' @param full \code{logical}, should only the estimate (\code{FALSE}) be returned, or additional data as well?(\code{TRUE}).
#' @param ... Arguments passed to class-specific methods.
#' @param exclude The list of faces that is to be excluded from the calculation
#' @return A proportion that gives the ratio of hole-cells compared to all occupied cells.
#' @examples
#' # 1. Only cells
#' # create a grid
#' hex <- hexagrid(2, sf=TRUE)
#'
#' # an example shape
#' shape <- paste0("F", c(4, 5, 11, 13, 15, 21, 24, 26, 32, 33, 34, 35, 36))
#' plot(hex, border="gray80")
#' plot(hex, shape, col="red", add=TRUE)
#'
#' # the gappiness
#' gap <- gappiness(shape, hex, full=TRUE)
#' gap
#' # plot the first gap
#' plot(hex, names(gap$holes[gap$holes==1]), col="#00555577", add=TRUE)
#' 
#' # points
#' set.seed(7)
#' rand <- icosa::rpsphere(20, output="polar")
#' occ <- locate(hex, rand)
#' plot(hex, border="gray80")
#' plot(hex, occ, col="red", add=TRUE)
#' points(rand, col="green", pch=3, cex=2)
#' 
#' # gappiness
#' gap <- gappiness(occ, hex, full=TRUE)
#' plot(hex, names(gap$holes), col="orange", add=TRUE)
#' @rdname gappiness
#' @exportMethod gappiness
setGeneric(
	name="gappiness",
	def=function(x,s,...){
		standardGeneric("gappiness")
	}
)

#' @rdname gappiness
setMethod(
	"gappiness",
	signature=c(x="character", s="trigrid"),
		definition=function(x, s, exclude=NULL, full=FALSE){
			if(any(!x%in%faces(s))) stop("All face names have to part faces of the grid.")

			# get a unique list of occupied faces
			occup <- unique(x)

			if(all(x%in%exclude)){
				gap <- NA
				theHoles <- NULL
			}else{

				# the holes of this patch
				theHoles<- icosa::holes(s,occup)

				# gappiness is defined as number of hole cells divided by the holes and shape itself
				gap <- (length(theHoles)-length(exclude))/(length(theHoles) + length(occup)-length(exclude))
			}

			# return both the metric and the way to plot it
			if(full){
				result <- list(estimate=gap, holes=theHoles, occupied=occup)
			}else{
				result <- gap
			}

		# return theappiness
		return(result)
	}
)


# coordinate pairs - relies on the character method
#' @rdname gappiness
setMethod(
	"gappiness",
	signature=c(x="matrix", s="trigrid"),
	definition=function(x,s,long=NULL,lat=NULL, duplicates=FALSE, plot=FALSE, plot.args=NULL, full=FALSE, exclude=NULL){

		# if locality is given
		if(!is.null(long) & !is.null(lat)) x <- x[,c(long, lat), drop=FALSE]
		if(!duplicates) x <- unique(x)

		# omit missing
		notMiss <- !is.na(x[,1]) & !is.na(x[,2])
		if(sum(notMiss)!=0){
			x <- x[notMiss, , drop=FALSE]

			# get the list of faces occupied
			faceList <- locate(s, x)

			# calculate the gappiness based on the faces
			# calculate full output
			result <- gappiness(
				x=unique(faceList),
				s=s,
				exclude=exclude,
				full=TRUE)
		}else{
			# construct the structure manually
			result <- list(
				estimate = NA,
				holes= NULL,
				occupied = NULL
			) 
		}
		if(plot & sum(notMiss)!=0){
			if(is.null(plot.args)) plot.args <- list(col="#BB000055", lwd=2, border="white")
			arguments <- c(list(x=s, y=result$occupied), plot.args)
			# if no plots are open yet, make one!
			if(grDevices::dev.cur()<=1) plot(x)
			do.call(plot, arguments)
		}
		# streamline output
		if(!full){
			result <- result$estimate
		}
		
		return(result)
	}
)
