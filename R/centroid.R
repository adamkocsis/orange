################################################################################
# Centroid-related code
################################################################################

qTest <- FALSE

#' Calculate ranges with the occupancy method
#'
#' @param x Eiher a 2-column numeric matrix with two columns: longitudes and latitudes, or a \code{data.frame} with these columns.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}), if here is any.
#' @param tax \code{character}, used only in the \code{data.frame} method. Column name of groups (e.g. taxa) that allows the iteration of the method for multiple groups.
#' @param plot.args List arguments passed to the plotting function: \code{points}.
#' @param long \code{character}, column name of the longitudes.
#' @param lat \code{character}, column name of the latitudes.
#' @param q Minimum occupancy with \code{q} proportion of occurrences.
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
	definition=function(x,long=NULL,lat=NULL, duplicates=FALSE, plot=FALSE, plot.args=NULL){
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
		if(plot){
			if(is.null(plot.args)) plot.args <- list(col="#BB0000", pch=16, cex=2)
			arguments <- c(list(x=matrix(cent, ncol=2)), plot.args)
			# if no plots are open yet, make one!
			if(dev.cur()==1) plot(x)
			do.call(points, arguments)
		}
		return(cent)
	}
)

# coordinate pairs
#' @rdname centroid
setMethod(
	"centroid",
	signature=c(x="data.frame"),
	definition=function(x,tax=NULL, long="long", lat="lat", duplicates=FALSE, plot=FALSE, plot.args=NULL){

		if(!all(c(long, lat) %in% colnames(x))) stop("The 'long' and 'lat' parameters must be valid column names.")
		x <- x[,c(tax, long, lat)]
		if(!duplicates) x <- unique(x)

		# the result
		if(!is.null(tax)){
			if(plot) warning("Multi-taxon plotting is not yet supported.")
			# iterate with tapply
			resRaw <- tapply(X=x[, c(long,lat)], INDEX=x[,tax], FUN=function(a){
				centroid(as.matrix(a), duplicates=duplicates, plot=FALSE, plot.args=NULL)
			})

			#make this palateable...
			res <- matrix(NA, ncol=2, nrow=length(resRaw))
			rownames(res) <- names(resRaw)
			colnames(res) <- names(resRaw[[1]])
			for(i in 1:length(resRaw)) res[i,] <- resRaw[[i]]

		}else{
			# fall back to matrix method
			res <- centroid(as.matrix(x), duplicates=duplicates, plot=plot, plot.args=plot.args)
		}

		return(res)
	}
)

