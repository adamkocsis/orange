

# getting a cross product
create3D <- function(x) head(c(x, rep(0, 3)), 3)
cross <- function(x, y, i=1:3) {
  x <- create3D(x)
  y <- create3D(y)
  j <- function(i) (i-1) %% 3+1
  return (x[j(i+1)]*y[j(i+2)] - x[j(i+2)]*y[j(i+1)])
}


# return the z coordinate for an x and y coordinate pair for a plane defined by 3 points in planePoints
# @param x a vector of x coordinates
# @param y a vector of y coordinates
# @planePoints A 3-by-3 matrix of three points in given as xyz coordinates (in columns)
# @return the z coordinate on the plane
planez<- function(x, y, planePoints){

	# the points
	P1 <- planePoints[1,]
	P2 <- planePoints[2,]
	P3 <- planePoints[3,]

	# 1st find the Omega: the origin's projection in the plane
	A <- P2-P1
	B <- P3-P1

	# get the normal
	n <- cross(A,B)
	n <- n / sqrt(sum(n^2))
	k <- 0 - sum(n*P1)

	# the z coordinate on the  plane
	z <- (0-n[1]*x - n[2]*y -k)/n[3]

	# does not give good answers when n[3]=0, i.e. when the plane is parallel to z!!
	return(z)
}

#' Identify the center of a small circle based on three points of a sphere
#'
#' @param x A matrix of 3 points, either polar longitude and latitude coordinates, or XYZ Cartesian coordinates.
#' @param origin The origin of the sphere.
#' @param output If set to "polar" then the function will return the longitude and latitude of the small circle center.
#' @return A list with two elements, the center o the small circle and the surface radius of the small circle.
#' @examples
#' # generate 3 points on a sphere
#' set.seed(2)
#' ps <- icosa::rpsphere(3, output="polar")
#' small <- sc_center(x=ps)
#' circle<- sc_shape(x=small$center, r=small$r, output="polar")
#'
#' plot(NULL, NULL, xlim=c(-180, 180), ylim=c(-90, 90))
#' points(ps, col="red", pch=3, cex=4)
#' points(small$center, col="blue", pch=16)
#' points(circle, col="green", pch=16)
#' @export
sc_center <- function(x, origin=c(0,0,0), output="polar"){

#	origin <- c(0, 0,0)
	# translate the Cartesian
	if(ncol(x)==2){
		three <- icosa::PolToCar(x)
	}else{
		three <- x
	}
	if(nrow(x)!=3) stop("Exactly 3 points are required for this method.")


	# the points
	P1 <- three[1,]
	P2 <- three[2,]
	P3 <- three[3,]

	# 1st find the Omega: the origin's projection in the plane
	A <- P2-P1
	B <- P3-P1

	# get the normal
	n <- cross(A,B)
	n <- n / sqrt(sum(n^2))

	# the omega
	omega <- origin - ((origin-P1)%*%n)[1] * n

	# the center
	surfaceCenter <- omega/sqrt(sum(omega^2))*sqrt(sum(P1^2))

	# make sure that the output is a matrix
	surfaceCenter <- matrix(surfaceCenter, nrow=1)

	# the radius
	dis <- arcdist(P1, surfaceCenter)
	if(output=="polar"){
		surfaceCenter <- icosa::CarToPol(surfaceCenter)[-3]
		surfaceCenter <- matrix(surfaceCenter, nrow=1)
	}

	# the return
	return(list(center=surfaceCenter, r=dis))
}


# create small circle
standardSmallCircle <- function(rad, n, radius=6371.007){
	# the center of the circle
	x <- cos(rad)*radius

	# extrinsic radius
	rex <- sin(rad) * radius

	# the angles of the points around the center
	guide <- seq(0, 2*pi, length.out=n)

	# matrix of points
	result <- cbind(
		x= rep(x, n),
		y=cos(guide)*rex,
		z=sin(guide)*rex
	)

	return(result)
}

#' Generate coordinates of small circles
#'
#' Creating points along small circles 
#'
#' The function generates a single reference (standard) small circle around the (0,0) long-lat coordinates using a given radius.
#' This reference circle is then rotated to match a given center. The function is iterated for multiple radiiand centers (separately!).
#' @param x Numeric matrix or character string. Coordinates of the centers of the small circles
#' @param n Integer or NULL. only used if "x= random"
#' @param breaks Integer. The number of points to create from the small circle.
#' @param r The surface radius of the small circle expressed as great circle distance of the small circle's points from the center of the small circles.
#' @param r.ex The extrinsic radius of the small circles, i.e. the radius of the circle in the plane of the small circle.
#' @param r.rad The angle of the small circle in radians.
#' @param r.deg The angle of the small circle in degrees.
#' @param radius The radius of the sphere, defaults to the authalic radius of Earth - relevant only for Cartesian output.
#' @param origin The center of the sphere (don't touch this, unless you really think you know what you are doing!).
#' @param output Output structure for the function. The value output="polar" will return the polar (longitude-latitude) coordinates of the small circles in an array.
#' The setting output="cartesian" will return the 3D cartesian coordinates of the small circles in an array. The option output="sf" will return an sfc geometry collection.
#' @param drop If there is a single small circle to be generated, should its array wrapper be dropped?
#' @param sf.type The type of sf object to be returned.
#' @param sf.wrap.dateline Argument indicating whether returned sf object should be wrapped around the dateline.
#' @return Either a numeric array or and sfc geometry collection.
#' @export
#' @examples
#' set.seed(1)
#' # generate a random small circle
#' central <- sc_shape(r=6000, breaks=40)
#' plot(NULL, xlim=c(-180, 180), ylim=c(-90,90), xlab="longitude", ylab="latitude")
#' points(central, pch=17, col="gray")
#' # repeat, same radius multiple centers
#' smaller <- sc_shape(x=central, r=1000, breaks=100)
#' for(i in 1:nrow(central)) points(smaller[,,i], col=i, pch=16)
#'
#' # Random circles with different radii
#' diff200 <- sc_shape(r=runif(200, 100, 10000))
#' plot(NULL, xlim=c(-180, 180), ylim=c(-90,90), xlab="", ylab="", axes=F)
#' for(i in 1:(dim(diff200)[3])){
#'  points(diff200[,,i], col=i,pch=16)
#' }
#'
#' # Same center, different radius
#' cent <- central[10, ,drop=FALSE]
#' outward <- sc_shape(cent, r=seq(100, 15000, length.out=20)) 
#' plot(NULL, xlim=c(-180, 180), ylim=c(-90,90), xlab="", ylab="", axes=F)
#' for(i in 1:(dim(outward)[3])){
#'  points(outward[,,i], col=i,pch=16)
#' }
sc_shape <- function(x="random", r=NULL, r.ex=NULL, r.rad=NULL, r.deg=NULL, n=NULL, breaks=100,radius=6371.007, origin=c(0,0,0), output="polar", sf.type="polygon", sf.wrap.dateline=TRUE, drop=TRUE){

## x="random"
## r=NULL
## r.ex=NULL
## r.rad=NULL
## r.deg=NULL
## n=NULL
## breaks=100
## radius=authRadius
## origin=c(0,0,0)
## output="polar"
## sf.type="polygon"
## sf.wrap.dateline=TRUE


	# if the surface radius is given
	if(!is.null(r)) r.rad <- r / radius

	# if the extrinsic radius is given
	if(!is.null(r.ex))  r.rad <- asin(r.ex / radius)

	# if the degree is given
	if(!is.null(r.deg))  r.rad <- r.deg / 180 * pi


	# rotation result format to generate
	if(output=="polar" | output=="sf"){
		gen <- "polar"
		coordNames <- c("long", "lat")
	}else{
		gen <- "cartesian"
		coordNames <- c("x", "y", "z")
	}

	# how is this need to be processed?
	if(inherits(x, "character")){
		if(any(x!="random")) stop("Invalid 'x' argument.")
		if(is.null(n) & length(r.rad)<1) stop("You must provide the number of small circles to generate.")
		if(is.null(n)) n <- length(r.rad)
		if(length(n)!=1 | !is.numeric(n)) stop("Provide a single integer 'n'.")
		if(length(n)<1 | n%%1!=0) stop("Provide a single integer 'n'.")
			


		# to hold the data
		results <- array(NA, dim=c(breaks, length(coordNames), n))
		dimnames(results)[[2]] <- coordNames

		# same circle rotated may times
		if(length(r.rad)==1){

			# 1. create a standard small circle
			stand <- standardSmallCircle(r.rad,  n = breaks)

			# rotate the small circles
			for(i in 1:n){
				results[,,i]<- rotate(stand, angles="random", pivot=origin, output=gen)
			}
		}else{

			# rotate the small circles
			for(i in 1:n){
				# 1. create a standard small circle with a given radius
				stand <- standardSmallCircle(r.rad[i],  n = breaks)

				# and rotate it
				results[,,i]<- rotate(stand, angles="random", pivot=origin, output=gen)
			}

		}

		# drop array
		if(n==1 & drop){
			results <- results[,,1]
		}

	}

	if(inherits(x, "numeric")) if(is.null(dim(x))) x <- matrix(x, nrow=1)
	if(inherits(x, "matrix")){
		# make sure it is longitude latitude
		if(ncol(x)==3) x <- CarToPol(x, norad=TRUE)

		if(length(r.rad)==1 & nrow(x)>=1){
			# 1. create a standard small circle
			stand <- standardSmallCircle(r.rad,  n = breaks)

			# the number of points to generate
			n <- nrow(x)

			# to hold the data
			results <- array(NA, dim=c(breaks, length(coordNames), n))
			dimnames(results)[[2]] <- coordNames

			for(i in 1:n){
				results[,,i]<- rotate(stand, long=x[i,1], lat=x[i,2], reflong=0, pivot=origin, output=gen)
			}
		}
		if(length(r.rad)>1 & nrow(x)==1){
			# the numer of circles to generate
			n<- length(r.rad)
			# to hold the data
			results <- array(NA, dim=c(breaks, length(coordNames), n))
			dimnames(results)[[2]] <- coordNames

			for(i in 1:n){
				# 1. create a standard small circle
				stand <- standardSmallCircle(r.rad[i],  n = breaks)
				results[,,i]<- rotate(stand, long=x[1,1], lat=x[1,2], reflong=0, pivot=origin, output=gen)
			}

		
		}
		# drop array
		if(n==1 & drop){
			results <- results[,,1]
		}

	}

	# if the output is an sf object
	if(output=="sf"){
		warning("This feature is not yet thoroughly tested!")
		if(!requireNamespace("sf", quiet=TRUE)) stop("This output option requires the 'sf' package.")

		polyList <- list()
		for(i in 1:n){
			coordSet <- results[,,i]
			polyList[[i]] <- sf::st_polygon(list(coordSet[c(1:nrow(coordSet),1), ]))
		}

		# create a geometry collection
		results <- sf::st_sfc(polyList, crs="WGS84")

		if(sf.wrap.dateline) results <- sf::st_wrap_dateline(results)

	}

	return(results)

}

#' Identify whether a given set of points is in, on or out of a small circle
#'
#' This method calculates the plane of the small circle and uses the normal to evaluate whether points are in or outside the small circle.
#' 
#' @param x A set of points to be determined a matrix of longtidues and latitudes.
#' @param center The center of a small circle
#' @param r The surface radius of a small circle.
#' @param origin The origin of the sphere.
#' @return A numeric vector indicating whether the points are in (1), on (0) or outside of the small circle.
#' @examples
#' # set.seed(1)
#' # generate a random small circle
#' set.seed(1)
#' pol <- icosa::rpsphere(1000, output="polar")
#' # use the first as a small circle's center
#' center <- pol[1, , drop=FALSE]
#' radius <- 5000 # define radius in km
#' 
#' # draw a small circle
#' circle <- sc_shape(x=center, r=radius, breaks=200)
#' 
#' # visualize
#' plot(NULL, NULL, xlim=c(-180, 180), ylim=c(-90, 90))
#' points(pol, col="red", pch=3, cex=2)
#' arcs(circle, col="blue", lwd=2)
#' 
#' # find points in small circle (about 10x faster than distance evaluation)
#' inside <- sc_in(x=pol, center=center, r=radius)
#' points(pol[inside==1, ], col="green", pch=3, cex=2)
#' 
#' @export
sc_in <- function(x, center, r,  origin=c(0,0,0)){

	# if the center is given as longlat, change to cartesian
	if(ncol(x)==2){
		xCart <- icosa::PolToCar(x)
	}else{
		xCart <- x
	}
	# if the center is given as longlat, change to cartesian
	if(ncol(center)==2){
		centerCart <- icosa::PolToCar(center)
	}else{
		centerCart <- center
	}

	# create 3 points on the small circle to get the plane - could be made faster
	three <- sc_shape(x=centerCart, r=r, breaks=4, output="cartesian", origin=origin)[1:3, ]

	# direction of the center either positive or negative
	centerDir <- sign(centerCart[3]-planez(centerCart[1], centerCart[2], planePoints=three))

	# the points
	pointDir <- sign(xCart[,3] - planez(xCart[,1], xCart[,2], planePoints=three))

	# should give 0s when the small circle coords are fed to the function!
	res <- rep(NA, nrow(xCart))

	res[pointDir==centerDir] <- 1
	res[pointDir!=centerDir] <- -1
	res[pointDir==0] <- 0
	return(res)
}
