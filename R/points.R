#' Generate presences based on different models
#'
#' @param n The number of presences ot generate
#' @param method The method to generate presences
#' @export
#' @rdname presences
presences <- function(n, method="kent",...){
	if(method=="kent"){
		res <- presences_kent(n=n, ...)
	}

	if(method=="kent-icosa"){
		res <- presences_kent_icosa(n=n, ...)
	}
	return(res)

}


#' Generate random presences based on the Kent distribution
#'
#' The function is a wrapper around the \code{rkent} function in the \code{Directional} package.
#'
#' @param n Either a single integer or a vector with the number of points to be generated in one distribution.
#' @param kappa A vector of kappa parameters of the Kent distribution. kappa = a indicates a uniform distribution, the higher kappa is, the concentrated the presences.
#' @param centers The center points around which the points are to created (matrix). If \code{NULL}, then these will be randomized.
#' @param beta The ovalness parameter of the Kent distribution
#' @param listout For multiple distributions, should the outputs be separated (list) or returned in a single matrix?
#' @return A two-column matrix of point coordinates or a list of matrices.
#' @examples
#' # uniform points on a sphere
#' uniform <- presences_kent(n=400, kappa=0)
#' # points around the 50 lon, 80 lat point
#' conc <- presences_kent(n=100, kappa=50, centers=matrix(c(50,80), ncol=2))
#' plot(NULL, NULL, xlim=c(-180, 180), ylim=c(-90,90))
#' points(uniform, col="red", pch=1)
#' points(conc, col="blue", pch=3)
#' @export
presences_kent <- function(n, kappa, centers=NULL, beta=0, listout=FALSE){
	if(!requireNamespace("Directional", quietly=TRUE)) stop("This function requires the 'Directional' package.")


	# the number of pointsets to simulate
	s <- length(n)

	# the number of species
	if(length(kappa)==1) kappa <- rep(kappa, s)

	# if centers not given, do it randomly
	if(is.null(centers)){
		centers <- icosa::rpsphere(s, output="polar")
	}

	# the centers on a unit sphere
	unitCenters <- icosa::PolToCar(centers, radius=1)

	if(!listout){
		pointSet<- matrix(NA, ncol=2, nrow=sum(n))
		cumN <- cumsum(n)
	}else{
		# generate as man
		pointSet <- list()
	}
	
	for(i in 1:s){

		#this shit is bugy
		go <- TRUE
		# try to do this.... iteratively, gives an error with some random numbers
		while(go){
			try({
				# points along a great circle
				unitPoints <- Directional::rkent(n=n[i], k=kappa[i], m=as.numeric(unitCenters[i,,drop=FALSE]), b=beta)
				go <- FALSE
			}, silent=TRUE)
		}

		if(listout){
			# back to longitude and latitude
			pointSet[[i]]  <- icosa::CarToPol(unitPoints)[,c(1,2)]
		}else{
			if(i==1) offset <- 0 else offset <- cumN[i-1]
			pointSet[1:n[i]+offset,]  <- icosa::CarToPol(unitPoints)[,c(1,2)]

		}
	}

	return(pointSet)

}

#' Generate presences from an icosahedral grid frequency table using Kent distributions
#'
#' The function will generate presences from Kent distributions that are centered around the face centers
#' of the given icosahedral grid. The
#'
#' @param x A table or named numeric, with integers and names as the face names of an icosahedral grid
#' @param icosa An icosahedral grid.
#' @param kappa Either a single kappa value, or a vector. The argument \code{NULL} will use 40000*4 divided
#' by the average vertex radius of the grid (selected arbitrarily).
#' @param times Single integer: a multiplier indicating how many times the should the sampling be repeated.
#' @return A 2-column matrix of points.
#' @export
#' @examples
#' # example grid
#' hex <- hexagrid(deg=5, sf=TRUE)
#'
#' # example data
#' sampledFaces <- paste0("F", 50:70)
#' tab <- rep(10, length(sampledFaces))
#' names(tab) <- sampledFaces
#'
#' # generate points
#' ps <- presences_kent_icosa(x=tab, icosa=hex, kappa=NULL, times=1)
#'
#' # visualize
#' plot(hex, tab, reset=FALSE)
#' points(ps, col="red", pch=3)
presences_kent_icosa <- function(x, icosa, kappa=NULL, times=1){
	if(times%%1!=0 | times<=0) stop("The times argument hast to a positive integer")
	# the centers of the faces
	cents<-icosa::centers(icosa)[names(x), ]

	# calculate kappa
	if(is.null(kappa)){
		# calculate mean vertexradius of the grid
		radius <- mean(icosa::vertexradius(icosa, degree=FALSE), na.rm=TRUE)
		# the kappa value
		kappa <- 40000/radius*4
	}

	# calculate the points
	ps <- presences_kent(centers=cents, n=as.numeric(x)*times, kappa=kappa)

	# return the points
	return(ps)
}
