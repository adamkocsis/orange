# Lower-level functions related to density estimation based on Kent distributions.
# See densityarea.R for the user interface function.

################################################################################
# dkentgeo: Calculate densities of Kent distributions
################################################################################

#' The density of the Kent distribution reparametrized for geographic domain
#'
#' Return the density values of a Kent distribution for various input
#'
#' The function relies on the Directional package.
#' @param x Geographic domain of the distribution. Either a longitude-latitude coordinate matrix, a SpatRaster object.
#' @param mu The center of the distribution given with a vector, a single-row matrix of longitude-latitude values. 
#' @param kappa The Kappa concentration parameter.
#' @param majorheading An initial heading ([0,360])from the center of the distribution for the direction of the major axis (0: north, 90: east).
#' @param majorpoint A pair of coordinates that define the direction of the major axis from the center of the distribution.
#' @param beta The beta parameter defining the excentricity of the distribution. \code{beta=0} defines an isometric distribution. 
#' @return The values of the density function, either as numeric vector, or an input-compatible data container (\code{SpatRaster}).
#' @export
#' @examples
#' # points
#' x <- rpsphere(100, output="polar")[,1:2]
#' mu <- rpsphere(1, output="polar")[,1:2, drop=FALSE]
#' coordres<- dkentgeo(x=x, mu=matrix(0, ncol=2, nrow=1), kappa=100, beta=0)
#' 
#' # terra
#' if(requireNamespace("terra", quietly=TRUE)){
#' library(terra)
#' # no eccentricity 
#' res<- dkentgeo(x=rast(), mu=matrix(0, ncol=2, nrow=1), kappa=100)
#' plot(res)
#' res<- dkentgeo(x=rast(), mu=matrix(c(0,30), ncol=2, nrow=1), kappa=100)
#' plot(res)
#' res<- dkentgeo(x=rast(), mu=matrix(c(0,60), ncol=2, nrow=1), kappa=100)
#' plot(res)
#' res<- dkentgeo(x=rast(), mu=matrix(c(0, 30), ncol=2, nrow=1),kappa=10)
#' plot(res)
#' res<- dkentgeo(x=rast(), mu=matrix(c(0, 60), ncol=2, nrow=1),kappa=5)
#' plot(res)
#' res<- dkentgeo(x=rast(), mu=matrix(c(30,  60), ncol=2, nrow=1),kappa=5)
#' plot(res)
#' res<- dkentgeo(x=rast(), mu=matrix(c(120,  60), ncol=2, nrow=1),kappa=40)
#' plot(res)
#' 
#' # with some eccentricity 
#' res<- dkentgeo(x=rast(), mu=matrix(c(0,  0), ncol=2, nrow=1),
#' 	majorpoint=matrix(c( 60,0), ncol=2, nrow=1), kappa=20, beta=9.99)
#' plot(res)
#' res<- dkentgeo(x=rast(), mu=matrix(c(0,  0), ncol=2, nrow=1),
#' 	majorpoint=matrix(c(45, 45), ncol=2, nrow=1), kappa=20, beta=5)
#' plot(res)
#' res<- dkentgeo(x=rast(), mu=matrix(c(0,  0), ncol=2, nrow=1),
#' 	majorheading=0, kappa=20, beta=8)
#' plot(res)
#' res<- dkentgeo(x=rast(), mu=matrix(c(100,  -60), ncol=2, nrow=1),
#' 	majorheading=90, kappa=20, beta=8)
#' plot(res)
#' res<- dkentgeo(x=rast(), mu=matrix(c(0,  50), ncol=2, nrow=1),
#' 	majorheading=30, kappa=20, beta=8)
#' plot(res)
#' }
setGeneric(
	name="dkentgeo",
	package="orange",
	def=function(x, ...){
		standardGeneric("dkentgeo")
	}
)



# majdir is an another point that describes the direction of the major axis (this will be replaced with an azimuth)
# azimut is used to generate two points along a great circle arc
#' @rdname dkentgeo
setMethod(
	"dkentgeo",
	signature=c(x="matrix"),
	function(x, mu, kappa, majorheading=NULL, majorpoint=NULL, beta=0){
		if(!requireNamespace("Directional", quietly=TRUE)){
			stop("This function relies on the 'Directional' package. Please install it first. ")
		}
		if(beta!=0) if(beta>=(kappa/2)) stop("The 'beta' parameter must be lower than kappa/2.")
		if(beta<0) stop("The 'beta' parameter must be positive")
		# transform the domain's geographic coordinates to unit sphere
		coords <- icosa::PolToCar(x, radius=1, origin=c(0,0,0))
		# the mean direction
		center <- icosa::PolToCar(mu, radius=1, origin=c(0,0,0))

		# the orthogonal vectors need to be calculated, even when beta =0
		# if it is not given in any form, but there is a beta, stop
		if((is.null(majorheading) & is.null(majorpoint)) & beta!=0)
			stop("You must provide a major axis orientation if the distribution is elliptic.")

		# if a heading is given but not a point
		if(!is.null(majorheading)){
			majorpoint <- CoordFromDistAzimuth(coords=mu, azimuth=-majorheading, distance=100, radius=6371)
		}

		# no direction given, then it does not matter, randomize
		if(is.null(majorheading) & is.null(majorpoint)){
			majorpoint <- icosa::rpsphere(1, radius=1, origin=c(0,0,0)) 
		}
		majPoint <- icosa::PolToCar(majorpoint, radius=1, origin=c(0,0,0))

		# calculate the G matrix for directional
		# find orthogonal vectors to meanDir
		# find another point along the great circle arc between center and majPoint:
		planePoints <- icosa::arcpoints(center, majPoint, breaks=1, origin=c(0,0,0), radius=1)
			
		# calculate the orthogonal to this
		P1 <- planePoints[1,]
		P2 <- planePoints[2,]
		P3 <- planePoints[3,]

		# 1st find the Omega: the origin's projection in the plane
		A <- P2-P1
		B <- P3-P1

		# get the normal
		minorVector <- cross(A,B) # function in sc.R
		minorVector <- minorVector / sqrt(sum(minorVector^2))

		# the third orthogonal with the same method
		majorVector <- cross(center, minorVector)
		
		# the three orthogonal vectors are center, major and minor, they make up the gamma matrix
		gamma <- cbind(as.numeric(center), as.numeric(majorVector), as.numeric(minorVector))

		# the result of the kent density
		dvals <- Directional::dkent(y=coords, G=gamma, param=c(kappa, beta))

		return(dvals)
	}
)

#' @rdname dkentgeo
setMethod(
	"dkentgeo",
	signature=c(x="SpatRaster"),
	function(x, ...){
		if(!requireNamespace("terra", quietly=TRUE)){
			stop("This function relies on the 'terra' package. Please install it first. ")
		}
		# calculate density for this..
		xy <- terra::xyFromCell(x, 1:terra::ncell(x))
		# res<- dgeo(x=xy, mu=matrix(0, ncol=2, nrow=1), majdir=NULL, kappa=0.25, beta=0)
		res<- dkentgeo(x=xy, ...)
		# set values and return
		terra::values(x) <- res
		return(x)

	}
)


## # sf geometry method
## setMethod(
## 	"dkentgeo",
## 	signature=c(x="sfc"),
## 	function(x, mu, kappa, beta){

## 	}
## )

## # sf geometry method
## setMethod(
## 	"dkentgeo",
## 	signature=c(x="trigrid"),
## 	function(x, mu, kappa, beta){

## 	}
## )


################################################################################
# kde: Kent-based Kernel Density Estimation 
################################################################################

#' Geographic kernel density estimation using Kent distributions
#'
#' @param x Point distribution data (matrix)
#' @param s Output spatial structure. Currently only a SpatRaster.
#' @param kappa The concentration parameter of the distribution (the inverse of a 'bandwidth')
#' @param beta Parameter of the Kent distribution controlling eccentricity
#' @param majorheading The heading of the major axes.
#' @param ... Arguments passed to the \code{dkentgeo} function specifying the kernels.
#' @return The values of the density function.
#' @export
#' @rdname kde
#' @examples
#' if(requireNamespace(terra, quietly=TRUE)){
#' # generate 50 random points
#' set.seed(1)
#' x <- icosa::rpsphere(50, output="polar")[,c(1:2)]
#' 
#' # wide kernel
#' a<- kde(x=x, s=rast(res=2), kappa=50)
#' plot(a)
#' points(x)
#' 
#' 
#' # mid size
#' a<- kde(x=x, s=rast(res=2), kappa=10)
#' plot(a)
#' points(x)
#' 
#' # larger
#' a<- kde(x=x, s=rast(res=2), kappa=5)
#' plot(a)
#' points(x)
#' 
#' # anisometric
#' a<- kde(x=x, s=rast(res=2), kappa=5, beta=2.4, majorheading=90)
#' plot(a)
#' points(x)
#' 
#' # a kent sample
#' data(kentsamples)
#' da <- kentsamples$dateline_m[1:50,]
#' a<- kde(x=da, s=rast(res=2), kappa=50)
#' plot(a)
#' points(da, col="#FFFFFF44")
#' 
#' }
setGeneric(
	name="kde",
	package="orange",
	def=function(x,s, ...){
		standardGeneric("kde")
	}
)

#' @rdname kde
setMethod(
	"kde",
	signature=c(x="matrix", s="SpatRaster"),
	function(x, s, kappa, beta=0, majorheading=NULL){
		# stop on missing values!
		if(!requireNamespace("terra", quietly=TRUE)){
			stop("This function relies on the 'terra' package. Please install it first. ")
		}
		# initialize the raster to 0s
		xy <- terra::xyFromCell(s, 1:terra::ncell(s))
		res <- rep(0, nrow(xy))

 		# loop through all points
		for(i in 1:nrow(x)){
			# this particular point
			one <- dkentgeo(x=xy, mu=x[i,, drop=FALSE], majorheading=majorheading, kappa=kappa, beta=beta)
			# aggregate
			res <- res + one
		}
		terra::values(s) <- res

		# still needs to be scaled
		# to make this a proper density function, the volume under this surface must integrate to 1!
		# the area of cells relative to one
		unitArea <- terra::cellSize(s)/sum(terra::values(terra::cellSize(s)))
		
		# scaled density, when multiplied with the unit area, this will sum up to one
		news <- s/sum(terra::values((s * unitArea)))

		# a separate function needs to be written to do the limiting
		return(news)

	}
)


## # trigrid method
## setMethod(
## 	"kde",
## 	signature=c(x="trigrid"),
## 	function(x, mu, kappa, beta){

## 	}
## )


################################################################################
# confregion: define confidence regions in the 
################################################################################

#' Confidence region of a geographic kernel density estimator
#'
#' @param x The result of the a geographic kernel density estimation (kde)
#' @param level A single numeric, confidence level (0-1).
#' @param area A single logical, should areas be returned?
#' @param prop A single logical, should areas returned be the proportion of sphere (\code{TRUE}), or square meters (\code{FALSE})?
#' @return If \code{area} is \code{FALSE}, then the same type as \code{x}, with binary values indicating whether a spatial unit is inside the region or not.
#' Otherwise a \code{list}, with the the region and the area.
#' @rdname confregion
#' @export
#' @examples
#' # set seed
#' set.seed(2)
#' coords <- icosa::rpsphere(50, output="polar")[, c(1,2)]
#' if(!requireNamespace("terra", quietly=TRUE)){
#' library(terra)
#' # estimate density
#' den<- kde(x=coords, s=rast(res=2), kappa=20)
#' 
#' # visualize
#' plot(den)
#' points(coords)
#' 
#' # calculate confidence region
#' cr <- confregion(den, 0.95)
#' cr
#' 
#' # the region
#' plot(den)
#' points(coords)
#' plot(cr$region, add=TRUE, col="#ff000044", legend=FALSE)
#' # the area
#' cr$area
#' }
setGeneric(
	name="confregion",
	package="orange",
	def=function(x, ...){
		standardGeneric("confregion")
	}
)

# conditional definition
#if(!requireNamespace("terra", quietly=TRUE)){
#' @rdname confregion
setMethod(
	"confregion",
	signature=c(x="SpatRaster"),
	function(x, level, area=TRUE, prop=TRUE){
		# calculate the area
		allArea <- terra::cellSize(x)
		unitArea <- allArea/sum(terra::values(allArea))

		# calculate the volume of the 
		volume <- x*unitArea

		# isolate vector values
		volVals <- terra::values(volume)
		sVals <- terra::values(x)

		# Aggregate volume of cells with decreasing density values.
		# the order of the density values in decreasing order
		decOrd <- order(sVals, decreasing=TRUE)

		# the reordered volume contributions
		volValsDec<- volVals[decOrd]
		cumulated <- cumsum(volValsDec)

		# the resulting raster
		conf <- x

		# assign the cumulated confidence levels
		terra::values(conf)[decOrd] <- cumulated

		# the confidence region
		res <- conf<level
		terra::values(res)[!terra::values(res)] <- NA

		# plot(res)

		if(area){
			if(prop){
				# to get the area, multipply logical mask with area
				# and sum up values
				areaVal <- sum(terra::values(res*unitArea), na.rm=TRUE)
			}else{
				areaVal <- sum(terra::values(res*allArea), na.rm=TRUE)
			}
			res <- list(region=res, area=areaVal)
		}

		return(res)
	}
)
#}
