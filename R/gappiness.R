#' The gappiness a of shape
#'
#' Gappiness refers to the internal gaps that defined by a set of discretized cells or a point set that covers some cells in a discretization structre.
#'
#' @param x The list of faces that are part of the shape.
#' @param s Spatial discretization strucutre (currently only \code{trigrid}-class icosahedral grid).
#' @param full \code{logical}, should only the estimate (\code{FALSE}) be returned, or additional data as well?(\code{TRUE}).
#' @param ... Arguments passed to class-specific methods.
#' @param exclude The list of faces that is to be excluded from the calculation
#' @return A proportion that gives the ratio of hole-cells compared to all occupied cells.
#' @examples
#' # create a grid
#' hex <- hexagrid(2, sf=TRUE)
#'
#' # an example shape
#' shape <- paste0("F", c(4, 5, 11, 13, 15, 21, 24, 26, 32, 33, 34, 35, 36))
#'
#' # the gappiness
#' gappiness(shape, hex)
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
				theHoles<- holes(s,occup)

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
#' @rdname latrange
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
			if(dev.cur()<=1) plot(x)
			do.call(plot, arguments)
		}
		# streamline output
		if(!full){
			result <- result$estimate
		}
		
		return(result)
	}
)
