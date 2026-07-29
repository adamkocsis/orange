################################################################################
# Occupancy-related code
################################################################################

qTest <- FALSE

#' Calculate ranges with the occupancy method
#'
#' High-level, abstract function to calculate how many or what proporotion of components in a pre-specified spatial structure is occupied by a distribution dataset
#'
#' The function by default returns counts (i.e. natural numbers) of how many discrete units are occupied.
#' However, there are many cases, when proportions are much more useful, which can be set with the \code{prop} argument.
#' Proportions are either global (\code{prop="global"}) or relative (\code{prop="relative"}). Global proportions express how much of the overall spatial structure (\code{s}) is occupied; for example, what is the proportion of cells in a grid that are occupied.
#' This method is only applicable, when \code{s} is defined as a spatial structure that is independent from \code{x}.
#' In contrast, relative proportional occupancies express proportions in the sampled set, i.e. what proportion of the overall sampled grid cells or localities are occupied by a taxon. This is particularly useful when \code{tax!=NULL}.
#' @param x Eiher a 2-column numeric matrix with two columns: longitudes and latitudes, or a \code{data.frame} with these columns.
#' @param s Structure to be occupied, either \code{NULL} (coordinate pairs), \code{character} (column name indicating locality) or a \code{trigrid} (icosahedral grid from the package icosa).
#' @param tax \code{character}, used only in the \code{data.frame} method. Column name of groups (e.g. taxa) that allows the iteration of the method for multiple groups.
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}).
#' @param plot.args List arguments passed to the plotting function: \code{lines}.
#' @param long \code{character}, column name of the longitudes.
#' @param lat \code{character}, column name of the latitudes.
#' @param q Minimum occupancy with \code{q} proportion of occurrences (not yet implemented!).
#' @param full Logical switch indicating whether only the estimate should be shown (\code{FALSE}), or other info (i.e. the list of occupied components in \code{s}) as well.
#' @param listarray If the full traceable output is required, should this be organized with list-array (native output of tapply).
#' @param prop Should counts be returned (\code{prop=NULL}), or proportions? If \code{prop="global"}, then global proportions are returned, if \code{prop="relative"}, relative proporitions are calculated.
#' @param ... Additional arguments passed to class-specific methods.
#' @return For single subsets (\code{tax=NULL}) either a single numeric or an orange list with an estimate and other information. Iterations for multiple taxa result in a named numeric vector or a list.
#' @rdname occupancy
#' @export
#' @examples
#' # I. Single taxon: Pinna nobilsi
#' # 1. Records
#' data(pinna)
#' # Subset to Pinna nobilis
#' nobilis <- pinna[pinna$species=="Pinna nobilis", ]
#'
#' # Number of unique coordinate pairs
#' cpairs <- occupancy(nobilis, long="decimalLongitude", lat="decimalLatitude")
#'
#' # 2. Occupancy in icosahedral grid 
#' hex <- hexagrid(deg=3, sf=TRUE)
#' plot(nobilis[c("decimalLongitude", "decimalLatitude")], pch=16, col="#00BBAA66")
#' plot(hex, border="gray70", add=TRUE)
#'
#' # calculate occupancy
#' occ <- occupancy(nobilis, s=hex, plot=TRUE, long="decimalLongitude", lat="decimalLatitude", full=TRUE)
#'
#' # manual coloring from full output
#' plot(hex, occ$occupied, add=TRUE, col="#00BB0088")
#'
#' # global proportional occupancies - relative to the grid
#' occprop <- occupancy(nobilis, s=hex, long="decimalLongitude", lat="decimalLatitude", prop="global")
#' # same as
#' occ$estimate/length(hex)
#'
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
	definition=function(x,long=NULL, lat=NULL, full=FALSE, prop=NULL){
		if(!is.null(prop)) prop <- match.arg(prop, c("global", "relative"))

		# if locality is given
		y <- x
		if(!is.null(long) & !is.null(lat)) y <- x[,c(long, lat)]
		if(ncol(y)!=2) stop ("You must provide a 2-column matrix!")
		# ensure missing values are omitted
		b1 <- !is.na(y[,1])
		b2 <- !is.na(y[,2])

		if(sum(b1)!=sum(b2)) warning("Some coordinate pairs are missing partially, they are omitted!")
		
		y <- y[b1 & b2, , drop=FALSE]
		 y  <- unique(y)
		# the result
		res <- nrow(y)
		if(!is.null(prop)){
			if(prop=="global") stop("Global proportions are not available for this input.")
			if(prop=="relative"){
				if(nrow(x)!=0){
					res <- 1.0 
				}else{
					res <- NA

				}
			}
		}
		if(full){
			res <- list(
				estimate=res,
				occupied=y
			)
			class(res) <- "orange"

		}
		return(res)
	}
)

# coordinate pairs
#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="data.frame", s="missing"),
	definition=function(x,tax=NULL, long="long", lat="lat", full=FALSE, listarray=TRUE, prop=NULL){
		# defend prop
		if(!is.null(prop)) prop <- match.arg(prop, c("global", "relative"))

		if(!all(c(long, lat) %in% colnames(x))) stop("The 'long' and 'lat' parameters must be valid column names.")
		y <- x[,c(tax, long, lat)]
		y <- unique(y)
		if(nrow(y)==0){
			if(!is.null(prop)){
				if(prop=="global") stop("Global proportions are not available for this input.")
			}
		}
		
		# taxon - iteration
		if(!is.null(tax)){
			# not full output
			if(!full){
				# omit any missing!
				y <- y[!is.na(y[,long]) & !is.na(y[,lat]) & !is.na(y[,tax]), ] 

				# in case of 0 rows
				if(nrow(y)==0) return(numeric())
				res <- table(y[, tax])
				resNum <- as.numeric(res)
				names(resNum) <- names(res)
				# are proportional occupancies required
				if(!is.null(prop)){
					if(prop=="global") stop("Global proportions are not available for this input.")
					if(prop=="relative"){
						# recursive call to get the total occupancy
						totalOccup <- occupancy(x, long=long, lat=lat, tax=NULL)

						# calculate the proportional occupancies
						resNum <-resNum/totalOccup
					}
				}
				return(resNum)
			# full output
			}else{
				# if this is to be a list-style output
				if(listarray){
					if(nrow(y)==0){
						res <-numeric()
						return(res)
					}
					# use a simple tapply to iterate functionally
					res <- tapply(
						INDEX=y[, tax],
						X=y[, c(long, lat)],
						FUN=function(a){
							occupancy(a, long=long, lat=lat, full=full)
						}

					)
					# if proprotions are needed
					if(!is.null(prop)){
						if(prop=="global") stop("Global proportions are not available for this input.")
						if(prop=="relative"){
							# recursive call to get the total occupancy
							totalOccup <- occupancy(x, long=long, lat=lat, tax=NULL)

							# calculate the proportional occupancies
							for(i in 1:length(res)) res[[i]]$estimate <- res[[i]]$estimate/totalOccup

						}
					}
					return(res)
				}else{
					stop("Donkey, snow, wine: not yet so fine!")
				}

			}
		# single taxon/set
		}else{

			# use the same as the matrix method
			res <- occupancy(as.matrix(y[, c(long, lat)]), full=full) 

			# are proportional occupancies required
			if(!is.null(prop)){
				if(prop=="global") stop("Global proportions are not available for this input.")
				if(prop=="relative"){
					val <- 1.0
					if(nrow(x)==0) val <- NA
					if(full){
						res$estimate <- val
					}else{
						res <- val
					}
				}
			}else{
				if(nrow(x)==0){
					if(full){
						res <- list(
							estimate=0,
							occupied=y[, c(long, lat)]
						)
						class(res) <- "orange"
						return(res)
					}else{
						return(0)
					}
				}

			}
			return(res)
		}
	}
)

# loc enries
#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="data.frame", s="character"),
	definition=function(x,s, tax=NULL, full=FALSE, prop=NULL, listarray=TRUE){
		# defend prop
		if(!is.null(prop)) prop <- match.arg(prop, c("global", "relative"))

		if(!any(s==colnames(x))) stop("The 'loc' argument must be a column in 'x'.")
		y <- x[,c(tax, s)]
		y <- unique(y)
		# make sure that there are no NAs!
		if(!is.null(tax)){
			 y <- y[!is.na(y[,tax]) & !is.na(y[,s]) , ]
			# the result
			if(!full){
				res <- table(y[, tax])
				resNum <- as.numeric(res)
				names(resNum) <- names(res)
				# proportional occupancies?
				if(!is.null(prop)){
					if(prop=="global") stop("Global proportions are not available for this input.")
					if(prop=="relative"){
						# recursive call to the total dataset
						totalOccup <- occupancy(x=x, s=s, tax=NULL)
						resNum  <- resNum/totalOccup
					}
				}
				result <- resNum
			}else{
				# if this is to be a list-style output
				if(listarray){
					# use a simple tapply to iterate functionally
					res <- tapply(
						INDEX=y[, tax],
						X=y[, s, drop=FALSE],
						FUN=function(a){
							occupancy(a, s=s, full=full)
						}
					)
					if(!is.null(prop)){
						if(prop=="global") stop("Global proportions are not available for this input.")
						if(prop=="relative"){
							# recursive call to the total dataset, total occupancy
							totalOccup <- occupancy(x=x, s=s, tax=NULL)
							# this has to be iterated separately..
							if(length(res)>0){
								for(i in 1:length(res)){
									res[[i]]$estimate <- res[[i]]$estimate/totalOccup
								}
							}
						}
					}
					result <- res
				}else{
					stop("Donkey, snow, wine: not yet so fine!")
				}


			}
		# no taxon iteration
		}else{
			occupied <- levels(factor(y))
			# promote for output consistency
			res <- as.numeric(length(occupied))
			# if proportional occupancies would be required
			if(!is.null(prop)){
				if(prop=="global") stop("Global proportions are not available for this input.")
				if(prop=="relative"){
					if(length(occupied)==0){
						res <- NA
					}else{
						res <- 1.0
					}
				}
			}
			if(full){
				result <- list(
					estimate=res,
					occupied=occupied
				)
				class(result) <- "orange"

			}else{
				result <- res
			}

		}
		return(result)

	}
)

#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="matrix", s="trigrid"),
	definition=function(x, s, long=NULL, lat=NULL, q=1, plot=FALSE, plot.args=NULL, full=FALSE, prop=NULL){
		# defend prop
		if(!is.null(prop)) prop <- match.arg(prop, c("global", "relative"))

		# if locality is given
		x <- unique(x)

		# if the columns are given, make sure they are interpreted right!
		if(!is.null(long) | !is.null(lat)) x <- x[, c(long, lat), drop=FALSE]

		if(q!=1 & !qTest) stop("Feature not yet finalized!")

		# invoke the internal, provides full results
		result <- occupancy_coords_icosa(x,s, q=q,  plot=plot, plot.args=plot.args)
		if(!is.null(prop)){
			if(prop=="global") result$estimate <- result$estimate/as.numeric(length(s))
			if(prop=="relative"){
				if(nrow(x)==0){
					result$estimate <- NA 
				}else{
					result$estimate <- 1.0 
				}
			}

		}

		# return a result
		if(!full) result <- result$estimate
		return(result)

	}
)

# uses the matrix method
#' @rdname occupancy
setMethod(
	"occupancy",
	signature=c(x="data.frame", s="trigrid"),
	definition=function(x, s, long="long", lat="lat", tax=NULL, q=1, plot=FALSE, plot.args=NULL, full=FALSE, prop=NULL){
		# defend prop
		if(!is.null(prop)) prop <- match.arg(prop, c("global", "relative"))

		# the same as the matrix method
		if(is.null(tax)){
			x <- as.matrix(x[, c(long, lat)])
			result <- occupancy(x, s=s, long=long, lat=lat,
				q=q, plot=plot, plot.args=plot.args, full=TRUE)
			if(!is.null(prop)){
				if(prop=="global") result$estimate <- result$estimate/as.numeric(length(s))
				if(prop=="relative"){
					if(nrow(x)==0){
						result$estimate <- NA 
					}else{
						result$estimate <- 1.0
					}
				}
			}
			# always return full output
			if(!full) result <- result$estimate
		}else{
			# method has to depend on whether there is
			# plotting or not!
			if(!plot){ # use tapply
				result <- tapply(
					INDEX=x[,tax],
					X=x[, c(long, lat)],
					FUN=function(a){

						# get rid of unnecessary recursion
						a <- as.matrix(a)
						result <- occupancy(a, s=s, long=long, lat=lat,
							q=q, full=full)
					}
				)
			}else{ # use a for loop - must be tested!, slower
				taxEntries <- levels(factor(x[,tax]))
				result <- list()
				for(i in 1:length(taxEntries)){
					a <- x[which(x[,tax]==taxEntries[i]), c(long, lat)]
					a <- as.matrix(a)
					# sequential running ensures plotting
					result[[i]] <- occupancy(a, s=s, long=long, lat=lat,
						q=q, full=full, plot=plot, plot.args=plot.args)
				}
				dim(result) <- length(taxEntries)
				names(result) <- taxEntries
				if(!full) result <- sapply(result,function(x) x[1])

			}

			# if the result is estimate-only get rid of dimensions, enforce numeric output
			if(!full){
				na <- names(result)
				result <- as.numeric(result)
				names(result) <- na

			}
			if(!is.null(prop)){
				# this has to be done manually here
				if(full){
					if(prop=="global"){
						# calculate the number of faces in the grid
						gridSize <- as.numeric(length(s))
						# calculate the proportions
						if(length(result)>0){
							for(i in 1:length(result)) result[[i]]$estimate <- result[[i]]$estimate/gridSize
						}
					}
					if(prop=="relative"){
						# recursive call to occupancy, to figure out the occupied faces by the total dataset -> what is up with q???
						totalOccup <- occupancy(x=x, s=s, long=long, lat=lat, tax=NULL, q=1)
						if(length(result)>0){
							for(i in 1:length(result)) result[[i]]$estimate <- result[[i]]$estimate/totalOccup
						}
					}
				}else{
					if(prop=="global"){
						# calculate the number of faces in the grid
						gridSize <- as.numeric(length(s))
						# calculate the proportions
						result <- result/gridSize
					}
					if(prop=="relative"){
						# recursive call to occupancy, to figure out the occupied faces by the total dataset -> what is up with q???
						totalOccup <- occupancy(x=x, s=s, long=long, lat=lat, tax=NULL, q=1)
						result <- result/totalOccup
					}

				}
			}
		}
		return(result)
	}
)


################################################################################
# INTERNAL
################################################################################
 
# internal method: look up coordinates with icosa grids
# x: 2 column matrix
occupancy_coords_icosa <- function(x, icosa, plot=FALSE, plot.args=NULL, q=1){
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
			occupied=occupCells,
			freq=freq
		)


		if(plot){
			if(is.null(plot.args)) plot.args <- list(col="#55000033")
			# if no plots are open yet, make one!
			if(dev.cur()==1) arguments$add <- NULL
			arguments <- c(list(x=icosa, y=res$occupied, add=TRUE), plot.args)
			do.call(icosa::plot, arguments)
		}

	}else{
		# the occupied cell
		occupCells <- unique(cells)

		# the result object
		res <- list(
			estimate=length(unique(occupCells)),
			occupied=occupCells
		)

		if(plot){
			if(is.null(plot.args)) plot.args <- list(col="#55000033")
			arguments <- c(list(x=icosa, y=res$occupied, add=TRUE), plot.args)
			# this will force par()
			# Ensure that sometihng is plotted
			openPlot<- is.null(try(text(0,0, label =""), silent=TRUE))
			if(!openPlot){
				arguments$add <- NULL
			}
			do.call(icosa::plot, arguments)
		}

	}
	class(res) <- "orange"

	# provide only estimate by default
	return(res)

}
