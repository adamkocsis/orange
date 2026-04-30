################################################################################
# Minimum spanning tree length family
################################################################################

qTest <- FALSE
#' Calculate ranges with the minimum spanning tree length method
#'
#' This family of metrics rely on constructing a minimum spanning tree from the distances between the points, and use its length to describe spreading.
#'
#' This metrics includes maximum great circle distance and similar methods.
#'
#' @param x Either a 2D numeric \code{matrix} with two columns: longitudes and latitudes, a \code{data.frame} with the same information.
#' @param s Structure to replace the points, either missing (coordinate pairs) or a \code{trigrid} (icosahedral grid from the package icosa).
#' @param long \code{character}, column name of the longitudes.
#' @param lat \code{character}, column name of the latitudes.
#' @param tax \code{character}, used only in the \code{data.frame} method. Column name of groups (e.g. taxa) that allows the iteration of the method for multiple groups.
#' @param duplicates \code{logical}, should identical coordinates be included in the calculation (default is \code{FALSE})
#' @param q \code{numeric}, a value between 0 and 1, the quantile. 
#' @param dm If there is a pre-made distance matrix, it can be plugged in here. If this is provided, the default coordinates will not be used.
#' @param full \code{logical}, should only the estimate (\code{FALSE}) be returned, or additional data as well?(\code{TRUE}).
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{sf::plot}.
#' @return A list with an estimate an two indices the rows of the input matrix that represent the length of the tree (or one of them).
#' @export
#' @rdname mstlength
#' @examples
#' # 1. Records
#' data(pinna)
#'
#' # 2. calculate and visualize
#' mst <- mstlength(pinna, long="decimallongitude", lat="decimallatitude", plot=TRUE, full=TRUE)
#'
setGeneric(
	name="mstlength",
	package="orange",
	def=function(x, s, ...){
		standardGeneric("mstlength")
	}

)

#' @rdname mstlength
setMethod(
	"mstlength",
	signature=c(x="matrix", s="missing"),
	definition=function(x, dm=NULL, long=NULL, lat=NULL, duplicates=FALSE, plot=FALSE, plot.args=NULL, full=FALSE, q=1){

		# if the columns are given, make sure they are interpreted right!
		if(!is.null(long) | !is.null(lat)) x <- x[, c(long, lat), drop=FALSE]

		# create a copy
		if(full) xOrig <- x

		# if given
		if(!duplicates) x <- unique(x)

		# omit missing
		x <- x[!is.na(x[,1]) & !is.na(x[,2]),, drop=FALSE]
		
		# for the q
		if(q!=1 & !qTest) stop("Feature not yet implemented!")

		# if locality is given
		# calculate the distance matrix
		if(is.null(dm)){
			if(length(x)==0){
				if(full){
					retObj <- list(estimate=NA, index=NA)
					class(retObj) <- "orange"
					return(retObj)
				}else{
					return(NA) 
				}
			}
			dm <- icosa::arcdistmat(x)
		}else{
			# check whether the supplied distance matrix is the same as th pointset
		}
		# q not implemented yet
		result <- mstlength_coords(x=x, dm=dm, q=q, full=full, plot=plot, plot.args=plot.args)

		## if(full){
		## 	# translate the indices to the original
		## }

		return(result)
	}
)

# uses the matrix method
#' @rdname mstlength
setMethod(
	"mstlength",
	signature=c(x="data.frame", s="missing"),
	definition=function(x, long="long", lat="lat", tax=NULL, dm=NULL, duplicates=FALSE, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
		# the same as the matrix method
		if(is.null(tax)){
			x <- as.matrix(x[, c(long, lat)])
			result <- mstlength(x, long=long, lat=lat, duplicates=duplicates, dm=dm,
				q=q, plot=plot, plot.args=plot.args, full=full)
		}else{
			if(full) stop("Full output not yet implemented for tax-wise iteration.")
			if(plot) warning("Multi-taxon plotting is not yet supported.")
			if(is.null(dm)){
				result <- tapply(
					INDEX=x[,tax],
					X=x[, c(long, lat)],
					FUN=function(a){

						# get rid of unnecessary recursion
						a <- as.matrix(a)
						result <- mstlength(a, long=long, lat=lat, duplicates=duplicates, dm=dm,
							q=q, plot=FALSE, plot.args=plot.args, full=full)
						}
					)
					resNames <- names(result)
					dim(result) <- NULL
					names(result) <- resNames

				
			# the distance matrix needs to be spliced!, and the
			}else{
				stop("Not yet!")
			}
		}
		return(result)
	}
)


################################################################################
# Internals
################################################################################
mstlength_coords <- function(x, dm, q=1, plot=FALSE, plot.args=NULL, full=FALSE){
	if(!requireNamespace("vegan", quietly=TRUE)) stop("This function requires the 'vegan' extension package. ")

	# case when there is just one row
	if(nrow(x)==1){
		estimate <- 0
		if(full){
			result <- list(
				estimate=estimate,
				index=matrix(NA, ncol=2, nrow=1),
				show=matrix(NA, ncol=2, nrow=1)
			)
		}else{
			result <- estimate

		}

	}else{
		# calculte the minimum spanning tree
		stre <- vegan::spantree(as.dist(dm))

		# the estimate
		estimate <- sum(stre$dist)


		# complete output necessary - eiher because of
		if(full | plot){

			# create an index input for the arcs function
			showThis <- matrix(ncol=2, nrow=0)
			index <- matrix(ncol=2, nrow=length(stre$kid))

			for(i in 1:length(stre$kid)){
				# the connected point index
				index[i, ] <- c(i+1, stre$kid[i])

				# point 1
				showThis <- rbind(showThis,
					matrix(c(
						x[ i+1, 1],
						x[ i+1, 2],
						x[stre$kid[i] ,1],
						x[ stre$kid[i], 2], NA, NA
					), ncol=2, byrow=TRUE)
				)
			}
		}
		if(full){
			result <- list(
				estimate=estimate,
				index=index,
				show=showThis
			)
		}else{
			#
			result <- estimate
		}

		if(plot){
			if(is.null(plot.args)) plot.args <- list(col="gray")
			arguments <- c(list(x=showThis), plot.args)
			if(dev.cur()==1) plot(x, pch=16)
			do.call(icosa::arcs, arguments)
		}
	}
	return(result)
}
