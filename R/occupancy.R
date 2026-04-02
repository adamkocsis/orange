################################################################################
# Occupancy-related code
################################################################################

qTest <- FALSE

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
	definition=function(x,tax, loc=NULL, long="long", lat="lat", duplicates=FALSE){
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

#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="matrix", icosa="trigrid", sf="missing"),
	definition=function(x, icosa, long=NULL, lat=NULL, duplicates=FALSE, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
		# if locality is given
		if(!duplicates) x <- unique(x)

		# if the columns are given, make sure they are interpreted right!
		if(!is.null(long) | !is.null(lat)) x <- x[, c(long, lat), drop=FALSE]

		if(q!=1 & !qTest) stop("Feature not yet finalized!")

		# invoke the internal
		result <- occupancy_coords_icosa(x,icosa, q=q,  plot=plot, plot.args=plot.args, full=full)

		# return a result
		return(result)

	}
)

# uses the matrix method
#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="data.frame", icosa="trigrid", sf="missing"),
	definition=function(x, icosa, long="long", lat="lat", tax=NULL, duplicates=FALSE, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
		# the same as the matrix method
		if(is.null(tax)){
			x <- as.matrix(x[, c(long, lat)])
			result <- occupancy(x, icosa=icosa, long=long, lat=lat, duplicates=duplicates,
				q=q, plot=plot, plot.args=plot.args, full=full)
		}else{
			result <- tapply(
				INDEX=x[,tax],
				X=x[, c(long, lat)],
				FUN=function(a){

					# get rid of unnecessary recursion
					a <- as.matrix(a)
					result <- occupancy(a, icosa=icosa, long=long, lat=lat, duplicates=duplicates,
						q=q, plot=plot, plot.args=plot.args, full=full)
					}
				)
		}
		return(result)
	}
)




# internal method: look up coordinates with icosa grids
# x: 2 column matrix
occupancy_coords_icosa <- function(x, icosa, plot=FALSE, plot.args=NULL, q=1, full=FALSE){
	# onit missing values
	x <- x[!(is.na(x[,1]) |  is.na(x[,2])),, drop=FALSE]

	if(nrow(x)==0){
		cells<- NULL
	}else{
		# the occupied cells by the points
		cells <- icosa::locate(icosa, x)
	}

	if(q!=1){
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
			# if no plots are open yet, make one!
			if(dev.cur()==1) arguments$add <- NULL
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
			# if no plots are open yet, make one!
			if(dev.cur()==1) arguments$add <- NULL
			do.call(icosa::plot, arguments)
		}

	}
	class(res) <- "orange"

	# provide only estimate by default
	if(!full) res <- res$estimate
	return(res)

}
