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


			if(all(x%in%exclude)){
				gap <- NA
				theHoles <- NULL
			}else{

				# the holes of this patch
				theHoles<- holes(s,x )

				# gappiness is defined as number of hole cells divided by the holes and shape itself
				gap <- (length(theHoles)-length(exclude))/(length(theHoles) + length(unique(x))-length(exclude))
			}

			# return both the metric and the way to plot it
			if(full){
				result <- list(estimate=gap, holes=theHoles)
			}else{
				result <- gap
			}

		# return theappiness
		return(result)
	}
)


## # uses the characer method
## #' @rdname gappiness
## setMethod(
## 	"gappiness",
## 	signature=c(x="matrix", s="trigrid"),
## 		definition=function(x, s, exclude=NULL){

## 			# get the list of faces occupied
## 			faceList <- locate(s, x)

## 			# calculate the gappiness based on the faces
## 			gap <- gappiness(faceList, s, exclude=exclude)

## 			# return the gappiness
## 			return(gap)
## 	}
## )
