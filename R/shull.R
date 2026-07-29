# Basic hull definitions - s=NULL
authRadius <- 6370.997
#' Spherical hull geometries and hull area calculations
#'
#' @param x Either a 2-column numeric matrix with two columns: longitudes and latitudes, or a \code{data.frame} with these columns.
#' @param s An external spatial structure to resolve the hull (e.g. a \code{hexagrid} object.)
#' @param plot Logical, should the result be plotted? Will plot over active plot (as in \code{add=TRUE}), if there is any.
#' @param long \code{character}, column name of the longitudes.
#' @param lat \code{character}, column name of the latitudes.
#' @param drop \code{logical} In case \code{s} is provided, should the returned object be a hull-class object, or only the identifiers of \code{s} (\code{drop=FALSE})?
#' @param metric \code{character} What metric should be used for area of the hull? \code{area} returns the values in square of the unit of \code{sphererad} (defaults to square km),
#' \code{prop} returns the area as a proportion of the sphere, \code{count} returns the area as the number of components it occupies in \code{s}, when applicable.
#' @param sphererad \code{numeric} The radius of the sphere used in the calculations, defaults to the authalic radius of Earth in km (6370.997).
#' @param duplicates \code{logical}, should identical coordinates be included in the calculation (default is \code{FALSE}).
#' @param ... Additional arguments passed to class-specific methods.
#' @return A \code{hull}-class object, a vector of identifiers (if applicable and \code{drop=FALSE}) or the area of te hull as a numeric. 
#' @rdname shull
#' @export
#' @examples
#' 
#' # example data
#' data(kentsamples)
#' p <- kentsamples$dateline_m
#' 
#' # plotting
#' plot(NULL, NULL, xlim=c(-180, 180), ylim=c(-90,90), xlab="", ylab="")
#' points(p, pch=3, col="red")
#' 
#' # basic spherical hull definition (native)
#' cc <- shull(p, method="centroidcircle")
#' plot(cc, add=TRUE)
#' 
#' # spherical hull area calculation
#' shullarea(cc)
#' 
#' # basic hull definition (icosahedral grid)
#' hex <- icosa::hexagrid(spacing=8, sf=TRUE)
#' # dropping hull container
#' cc_hex <- shull(p, s=hex, method="centroidcircle")
#' plot(hex, cc_hex, col="#0000DD33", border="#FFFFFF77", add=TRUE)
#' 
#' 
#' # with shull container
#' cc_hex2 <- shull(p, s=hex, method="centroidcircle", drop=FALSE)
#' shullarea(cc_hex2, s=hex,  metric="prop")
setGeneric(
	name="shull",
	package="orange",
	def=function(x, s,  ...){
		standardGeneric("shull")
	}
)

# coordinate pairs
#' @rdname shull 
setMethod(
	"shull",
	signature=c(x="matrix", s="missing"),
	definition=function(x,method, long=NULL,lat=NULL, duplicates=FALSE, plot=FALSE, plot.args=NULL, sphererad=authRadius){
		# if locality is given
		if(!is.null(long) & !is.null(lat)) x <- x[,c(long, lat), drop=FALSE]
		if(!duplicates) x <- unique(x)
		# omit missing
		notMiss <- !is.na(x[,1]) & !is.na(x[,2])
		if(sum(notMiss)!=0){
			x <- x[notMiss, , drop=FALSE]

			if(method=="centroidcircle"){
				hullData <- shull_centroid_circle_coords(x, sphererad=sphererad)
			}
		}else{
		}
		if(plot){
		 	if(is.null(plot.args)) plot.args <- list(col="#BB0000", pch=16, cex=2)
		 	arguments <- c(list(x=hullData, plot.args))
		 	## # if no plots are open yet, make one!
		 	## if(grDevices::dev.cur()<=1) plot(x)
		 	do.call(plot, arguments)
			
		}
		return(hullData)
	}
)


#' @rdname shull 
setMethod(
	"shull",
	signature=c(x="data.frame", s="missing"),
	definition=function(x,long="long",lat="lat",...){
		mat <- as.matrix(x[, c(long, lat)])
		res <- shull(mat, ...)
		return(res)

	}
)

# coordinate pairs with icosa
#' @rdname shull
setMethod(
	"shull",
	signature=c(x="matrix", s="trigrid"),
	definition=function(x,s, method,long=NULL,lat=NULL, duplicates=FALSE, plot=FALSE, plot.args=NULL, sphererad=authRadius, drop=TRUE){
		# calculate the hull geometry
		theHull <- shull(x, method=method, sphererad=sphererad)
		if(theHull$type=="centroidcircle"){
			# find the occupied cells - by distance?
			cent <- centers(s)
			message("Unfinished implementation! - Proper lookup of small circles required!")
			adm <- arcdistmat(points1=cent, points2=theHull$center, radius=sphererad)
			res <- rownames(cent)[adm<=theHull$radius]
			if(!drop){
				res=list(
					type="icosa",
					faces=res
				)
				class(res) <- "shull"

			}
		}

		if(plot){
		 	## if(is.null(plot.args)) plot.args <- list(col="#BB0000", pch=16, cex=2)
		 	## arguments <- c(list(x=hullData, plot.args))
		 	## ## # if no plots are open yet, make one!
		 	## ## if(grDevices::dev.cur()<=1) plot(x)
		 	## do.call(plot, arguments)
			
		}
		return(res)
	}
)



#' @rdname shull
#' @export
shullarea <- function(x, ...){
	UseMethod("shullarea", x)
}

#' @rdname shull
#' @export
shullarea.default <- function(x, method, metric="area", full=FALSE, sphererad=authRadius,... ){
		# construct the hull, pass methods!
		theHull <- shull(x, method=method, sphererad=sphererad,...)

		# recursively call to the hull area function
		res <- shullarea(theHull, full=full, ...)

		return(res)
}


#' @rdname shull
#' @export
shullarea.shull <-  function(x, metric="area", s=NULL, full=FALSE,  ...){
	if(x$type=="icosa"){
		if(metric=="area"){
			if(is.null(s)) stop("This metric requires the grid to be passed as 's'. ")
			result <- sum(icosa::surfacearea(s)[x$faces])
		}
		if(metric=="prop"){
			if(is.null(s)) stop("This metric requires the grid to be passed as 's'. ")
			result <- length(x$faces)/as.numeric(length(s))
		}
		if(metric=="count"){
			result <- length(x$faces)
		}
	}

	if(x$type=="centroidcircle"){
		if(metric=="count") stop("This hull method does not support a 'count' metric.")
		# the half-angle of the cap
		theta <- x$radius/x$sphererad

		if(metric=="area"){
			# spherical cap area
			result <- 2 * pi * x$sphererad^2 * (1- cos(theta))
		}
		if(metric=="prop"){
			# The proportion of the cap
			result <- 0.5 * (1- cos(theta))
			
		}
	}
	return(result)
}







################################################################################
# Internals
################################################################################
# centroid circle
# minimum smallcircle
# minimum ellipse
# convex s.s. 
# alpha hull


shull_centroid_circle_coords <- function(x, sphererad=authRadius){
	# caluclate the centroid
	focus <- centroid(x)
	# get the radius
	dists <- icosa::arcdistmat(points1=x, points2=matrix(focus, ncol=2), radius=sphererad)

	# get the maximum distance
	rad <- max(dists)
	focus <- matrix(focus, nrow=1)
	colnames(focus) <- c("long", "lat")
	# construct a hull object
	hullData <- list(
		type="centroidcircle",
		center=focus, radius=rad, sphererad=sphererad
	)
	class(hullData) <- "shull"
	return(hullData)
}



