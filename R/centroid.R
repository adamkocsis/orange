################################################################################
# Centroid-related code
################################################################################

qTest <- FALSE

#' Calculate ranges with the occupancy method
#'
#' @param x Eiher a 2-column numeric matrix with two columns: longitudes and latitudes, or a \code{data.frame} with these columns.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{lines}.
#' @param long \code{character}, column name of the longitudes.
#' @param lat \code{character}, column name of the latitudes.
#' @param q Minimum occupancy with \code{q} proportion of occurrences.
#' @param full Logical switch indicating whether only estimate should be shown (\code{FALSE}), or other info as well.
#' @return Either a single numeric or a list with an estimate and other information.
#' @rdname occupancy
#' @export
#' @examples
#' # 1. Canvas
#' hex <- hexagrid(deg=3, sf=TRUE)
#' plot(hex, reset=FALSE, xlim=c(-15, 40), ylim=c(25, 63))
#'
#' # 2. Records
#' data(pinna)
#'
#' # Number of unique coordinate pairs
#' cent <- centroid(pinna, long="decimallongitude", lat="decimallatitude")
#'
#' points(cent, col="darkred", pch=3, lwd=4, cex=4)
setGeneric(
	name="centroid",
	package="orange",
	def=function(x, ...){
		standardGeneric("centroid")
	}

)

# coordinate pairs
#' @rdname centroid
setMethod(
	"centroid",
	signature=c(x="matrix"),
	definition=function(x,long=NULL,lat=NULL, duplicates=FALSE){
		# if locality is given
		y <- x
		if(!is.null(long) & !is.null(lat)) y <- x[,c(long, lat)]
		if(duplicates) y <- unique(y)
		# omit missing
		y <- y[!is.na(y[,1]) & !is.na(y[,2]), ]

		# the result
		cent <- icosa::surfacecentroid(y)
		return(cent)
	}
)

## # coordinate pairs
## #' @rdname centroid
## setMethod(
## 	"centroid",
## 	signature=c(x="data.frame"),
## 	definition=function(x,tax=NULL, long="long", lat="lat", duplicates=FALSE){

## 		if(!all(c(long, lat) %in% colnames(x))) stop("The 'long' and 'lat' parameters must be valid column names.")
## 		y <- x[,c(tax, long, lat)]
## 		if(duplicates) y <- unique(y)

## 		# the result
## 		if(!is.null(tax)){
## 			res <- table(y[, tax])
## 		}else{
## 			res <- nrow(y)
## 		}
## 		resNum <- as.numeric(res)
## 		names(resNum) <- names(res)

## 		return(resNum)
## 	}
## )

## # loc enries
## #' @rdname occupancy
## setMethod(
## 	"occupancy",
## 	signature=c(x="data.frame", s="character"),
## 	definition=function(x,s, tax=NULL){

## 		if(!any(s==colnames(x))) stop("The 'loc' argument must be a column in 'x'.")
## 		y <- x[,c(tax, s)]
## 		y <- unique(y)
## 		# make sure that there are no NAs!
## 		y <- y[!is.na(y[,tax]) & !is.na(y[,s]) , ]

## 		# the result
## 		res <- table(y[, tax])
## 		resNum <- as.numeric(res)
## 		names(resNum) <- names(res)

## 		return(resNum)

## 	}
## )
