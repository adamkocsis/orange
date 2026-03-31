################################################################################
# Occupancy-related code
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
#' # just the coordinates
#' coordMat <- SimpleCoordinates(pinna, long="decimallongitude", lat="decimallatitude")
#' points(coordMat)
#'
#' # 3. calculate and visualize
#' occ <- occupancy(coordMat, icosa=hex, plot=TRUE)
#'
#' # plot(hex, occ$cells, add=TRUE, col="green")
setGeneric(
	name="occupancy",
	package="orange",
	def=function(x, icosa, sf, ...){
		standardGeneric("occupancy")
	}

)

#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="data.frame", icosa="missing", sf="missing"),
	definition=function(x,tax, loc=NULL, long="long", lat="lat", duplicates=FALSE, ...){
		# if locality is given
		if(!is.null(loc)){
			y <- x[,c(tax, loc)]
			if(!duplicates) y <- unique(y)

		# locality not given
		}else{
			y <- x[,c(tax, long, lat)]
			if(!duplicates) y <- unique(y)
		}
		
		# the result
		res <- table(y[, tax])
		resNum <- as.numeric(res)
		names(resNum) <- names(res)
				
		return(resNum)
	
	}
)


## occupancy <- function(coordMat, sf=NULL, icosa=NULL, plot=FALSE, plot.args=NULL, alpha=1){


## 	if(!is.null(icosa)){
## 		result <- occupancy_icosa(coordMat, icosa, plot=plot, plot.args=plot.args, alpha=alpha)
## 	}
## 	if(!is.null(sf)){
## 		stop("Not yet!")
## 	}

## 	return(result)

## }


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
