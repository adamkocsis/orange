################################################################################
# Occupancy-related code
################################################################################

qTest <- FALSE

#' Calculate ranges with the occupancy method
#'
#' @param x Eiher a 2-column numeric matrix with two columns: longitudes and latitudes, or a \code{data.frame} with these columns.
#' @param s Structure to be occupied, either \code{NULL} (coordinate pairs), \code{character} (column name indicating locality) or a \code{trigrid} (icosahedral grid from the package icosa).
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
#' cpairs <- occupancy(pinna, long="decimallongitude", lat="decimallatitude")
#'
#' # just the coordinates
#'
#' # 3. calculate and visualize
#' occ <- occupancy(pinna, s=hex, plot=TRUE, long="decimallongitude", lat="decimallatitude")
#'
#' # plot(hex, occ$cells, add=TRUE, col="green")
setGeneric(
	name="occupancy",
	package="orange",
	def=function(x, s, ...){
		standardGeneric("occupancy")
	}

)

# coordinate pairs
#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="matrix", s="missing"),
	definition=function(x,tax=NULL, long=NULL, lat=NULL){
		# if locality is given
		y <- x
		if(!is.null(long) & !is.null(lat)) y <- x[,c(long, lat)]
		 y <- unique(y)
		# the result
		resNum <- nrow(y)
		return(resNum)

	}
)

# coordinate pairs
#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="data.frame", s="missing"),
	definition=function(x,tax=NULL, long="long", lat="lat"){

		if(!all(c(long, lat) %in% colnames(x))) stop("The 'long' and 'lat' parameters must be valid column names.")
		y <- x[,c(tax, long, lat)]
		y <- unique(y)

		# the result
		if(!is.null(tax)){
			res <- table(y[, tax])
		}else{
			res <- nrow(y)
		}
		resNum <- as.numeric(res)
		names(resNum) <- names(res)
				
		return(resNum)
	}
)

# loc enries
#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="data.frame", s="character"),
	definition=function(x,s, tax=NULL){

		if(!any(s==colnames(x))) stop("The 'loc' argument must be a column in 'x'.")
		y <- x[,c(tax, s)]
		y <- unique(y)
		# make sure that there are no NAs!
		y <- y[!is.na(y[,tax]) & !is.na(y[,s]) , ]

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
	signature=c(x="matrix", s="trigrid"),
	definition=function(x, s, long=NULL, lat=NULL, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
		# if locality is given
		x <- unique(x)

		# if the columns are given, make sure they are interpreted right!
		if(!is.null(long) | !is.null(lat)) x <- x[, c(long, lat), drop=FALSE]

		if(q!=1 & !qTest) stop("Feature not yet finalized!")

		# invoke the internal
		result <- occupancy_coords_icosa(x,s, q=q,  plot=plot, plot.args=plot.args, full=full)

		# return a result
		return(result)

	}
)

# uses the matrix method
#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="data.frame", s="trigrid"),
	definition=function(x, s, long="long", lat="lat", tax=NULL, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
		# the same as the matrix method
		if(is.null(tax)){
			x <- as.matrix(x[, c(long, lat)])
			result <- occupancy(x, s=s, long=long, lat=lat,
				q=q, plot=plot, plot.args=plot.args, full=full)
		}else{
			result <- tapply(
				INDEX=x[,tax],
				X=x[, c(long, lat)],
				FUN=function(a){

					# get rid of unnecessary recursion
					a <- as.matrix(a)
					result <- occupancy(a, s=s, long=long, lat=lat,
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
		nOccs <- nrow(x)*q

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
