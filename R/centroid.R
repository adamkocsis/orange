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
		if(!is.null(long) & !is.null(lat)) x <- x[,c(long, lat), drop=FALSE]
		if(!duplicates) x <- unique(x)
		# omit missing
		notMiss <- !is.na(x[,1]) & !is.na(x[,2])
		if(sum(notMiss)!=0){
			x <- x[notMiss, , drop=FALSE]

			# the result
			cent <- icosa::surfacecentroid(x)
		}else{
			cent <- c(NA, NA)
			names(cent) <- c("long", "lat")
		}
		return(cent)
	}
)

# coordinate pairs
#' @rdname centroid
setMethod(
	"centroid",
	signature=c(x="data.frame"),
	definition=function(x,tax=NULL, long="long", lat="lat", duplicates=FALSE){

		if(!all(c(long, lat) %in% colnames(x))) stop("The 'long' and 'lat' parameters must be valid column names.")
		x <- x[,c(tax, long, lat)]
		if(!duplicates) x <- unique(x)

		# the result
		if(!is.null(tax)){
			# iterate with tapply
			resRaw <- tapply(X=x[, c(long,lat)], INDEX=x[,tax], FUN=function(a){
				centroid(as.matrix(a), duplicates=duplicates)
			})

			#make this palateable...
			res <- matrix(NA, ncol=2, nrow=length(resRaw))
			rownames(res) <- names(resRaw)
			colnames(res) <- names(resRaw[[1]])
			for(i in 1:length(resRaw)) res[i,] <- resRaw[[i]]

		}else{
			# fall back to matrix method
			res <- centroid(as.matrix(x), duplicates=duplicates)
		}

		return(res)
	}
)

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
