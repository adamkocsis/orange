#' Random shapes defined with random walks
#'
#' @param x The grid object
#' @param start The starting cell. defaults to NULL (random start)
#' @param steps The number of steps that the random walk will take. Directly controls the algorithm.
#' @param size The final size of the shape defined in the number of cells that constitute it. If set,
#' @param history Logical flag. Should the output include the complete history of cell occurrence through the walk, or only the cells that were visited? (default)
#'
#' @return A character vector, face names of \code{x}. If \code{history} is set to \code{FALSE} it includes every face name only once, in the order of visitation. If \code{history=TRUE} face names can be duplicated following the history of the random walk.
#' @export
#' @examples
#' hex <- hexagrid(8)
#' # random shape of 15 faces
#' rand <- rwfaces(hex, size=15)
rwfaces <- function(x, start=NULL, steps=NULL, size=NULL, history=FALSE){
	# if both step and size are given: halt
	if(!is.null(size) & !is.null(steps)) stop("Both the number of steps desired size are given.\nPlease provide only one!")

	# check dependency
	if(!requireNamespace("igraph", quietly=TRUE)) stop("This function requires the 'igraph' package.")

	# check whether the grid has a graph
	if(inherits(x@graph, "igraph")){
		gg <- x@graph
	}else{
		message("The grid's @graph slot is empty. Creating graph representation.")
		gg <- icosa::gridgraph(x)
	}

	# the number of faces to generate
	nFaces <- nrow(icosa::centers(x))

	# if starting face is given
	if(!is.null(start)){
		# Translate start: given in facename to an  index
		startIndex <- which(names(igraph::V(gg))==start)

	# if start is not given, start randomly
	}else{
		startIndex <- sample(nFaces,1)
	}


	# if nothing is given: completely random size - anything is possible
	if(is.null(size) & is.null(steps)) size <- sample(nFaces,1)


	# directly generate with a specified step number
	if(is.null(size)){

		# generate ar
		rwGraph <- igraph::random_walk(gg, start=startIndex, steps=steps)

		# the result of the random walk
		rw <- names(rwGraph)

		if(!history){
			# the output as faces - order is the first meeting
			res <- unique(rw)

		}else{
			res <- rw
		}

	# genearte until size is reached
	}else{
		# corner case
		rw <- NULL

		if(size!=0){
			# defend against idiocy
			# halt the function if value is too big
			if(size>nFaces) stop("Desired size is too big")

			# eyeball it...
			# the number already generated
			sizeGen <- 0

			# generate random sequence from here
			newStart <- startIndex

			# iterate this until the given size is reached
			while(sizeGen<size){
				# generate this chunk
				rwGraph <- igraph::random_walk(gg, start=newStart, steps=nFaces*6)
				rwPart <- names(rwGraph)

				# continue the sequence from the last face
				newStart <- which(names(V(gg)) == rwPart[length(rwPart)])

				# add this to the rest
				rw <- c(rw, rwPart)

				# calculate overall size
				sizeGen <- length(unique(rw))

			}

			# select the desired size
			cellList <- unique(rw)[1:size]

			if(!history){
				res <- cellList
			}else{
				# the first step at which the cell is encountered
				firstEncounter <- which(!duplicated(rw))

				# when is the last cell encountered
				lastStep <- firstEncounter[size]

				# all the steps
				res<- rw[1:lastStep]

			}

		}else{
			res <- rw
		}
	}
	return(res)

}



#' Random density fields defined by random walks
#'
#' @param coarse The coarse-resolution grid object (defines output).
#' @param fine The fine-resolution grid object (random walks operate on this grid).
#' @param start The starting cell on the fine-resolution grid. Defaults to NULL (random start)
#' @param steps The number of steps that the random walk will take, vectorized for multiple walks.  Directly controls the algorithm.
#' @param sizes The final size of the shape defined in the number of cells that constitute it. If set,
#'
#' @return A character vector, face names of \code{x}. If \code{history} is set to \code{FALSE} it includes every face name only once, in the order of visitation. If \code{history=TRUE} face names can be duplicated following the history of the random walk.
#' @export
#' @examples
#
#' coarse <- hexagrid(8, sf=TRUE)
#' dens <- rwdensity(coarse, steps=c(200, 500))
#' current <- dens[dens$shape==1, "count"]
#' names(current) <- dens[dens$shape==1, "face"]
#' plot(coarse, current)
#'
#' fine <- hexagrid(2)
#' # random shape of 15 faces
#' dist <- floor(rexp(15)*1500)
#' rand <- rwdensity(coarse, fine, steps=dist)
rwdensity <- function(coarse, fine=NULL, sizes=NULL, steps=NULL, start=NULL, tabulated=TRUE, df=TRUE){

	if(is.null(fine)){

		if(!is.null(sizes)){
			stop("Not yet!")
		}

		# if at least one steps argument is provided
		if(!is.null(steps)){

			# define containers
			faceNames <- NULL
			faceCount <- NULL
			iteration <- NULL

			for(i in 1:length(steps)){
				one <- rwfaces(coarse, steps=steps[i], start=start, history=TRUE)
				if(tabulated){
					# tabulate
					tab<- table(one)

					# store
					faceNames <- c(faceNames, names(tab))
					faceCount <- c(faceCount, as.numeric(tab))
					iteration <- c(iteration,rep(i, length(tab)))
				}
				cat(i, "\r")
				flush.console()
			}
			if(df){

				# the returned object
				res <- data.frame(
					face= faceNames,
					count = faceCount,
					shape = iteration
				)
				return(res)
			}else{
				return(tab)
			}

		}

	# both a coarse and a fine resolutino grid
	}else{

		# get the centers of the finer, basis grid
		fineCenters <- icosa::centers(fine)
		# locate centers, where they are on the coarse grid
		dict <- icosa::locate(coarse, fineCenters)
		if(any(is.na(dict))) stop("Exact alignment of 'fine' and 'coarse' is to be avoided.\nRotate the fine-scale grid!")

		# create dictionary
		names(dict) <- rownames(fineCenters)

		if(!is.null(sizes)){
			stop("Not yet!")
		}

		if(!is.null(steps)){

			# define containers
			faceNames <- NULL
			faceCount <- NULL
			iteration <- NULL
			fineface <- NULL


			for(i in 1:length(steps)){

				# generate fine resolution
				one <- rwfaces(fine, steps=steps[i], start=start, history=TRUE)

				# translate to other grid
				visits <- dict[one]

				if(tabulated){
					# tabulate
					tab<- table(visits)

					# store
					faceNames <- c(faceNames, names(tab))
					faceCount <- c(faceCount, as.numeric(tab))
					iteration <- c(iteration,rep(i, length(tab)))
				}else{

					faceNames <- c(faceNames, visits)
					iteration <- c(iteration, rep(i, length(visits)))
					fineface <- c(fineface,one)

				}
				cat(i, "\r")
				flush.console()
			}
		}

		if(tabulated){
			res <- data.frame(
				face= faceNames,
				count = faceCount,
				shape = iteration
			)
		}else{
			res <- data.frame(
				face= faceNames,
				shape = iteration,
				fine= fineface
			)


		}


		return(res)
	}
}

#' Random walk on a sphere
#'
#' @param x One pair of longitude and latitude coordinates in a matrix.
#' @param step Distance covered in a single step (km)
#' @param time The number of steps to take (integer)
#' @param output The
#' @return A matrix of coordinates.
#' @examples
#' randomWalk <- rwsphere(x=matrix(c(0,0), ncol=2), step=100, time=10)
#' @export
rwsphere <- function(x, step, time, output="polar"){

	# a couple of random points
	randomPoints <- icosa::rpsphere(time)

	# the point
	results <- matrix(NA, ncol=3, nrow=time+1)

	# add the first point to the results
	results[1,] <- icosa::PolToCar(x)

	# for every time bin
	for(i in 1:time){
		# push the previous point in a random direction with step
		results[i+1, ] <- pushFromPoint(
			sourcePoint = results[i,,drop=FALSE],
			dists=step,
			otherPoint=randomPoints[i,,drop=FALSE],
			output="cart")
	}
	if(output=="polar") {
		results <- icosa::CarToPol(results)[,1:2]
	}
	return(results)
}
