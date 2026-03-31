# getting a cross product
create3D <- function(x) head(c(x, rep(0, 3)), 3)
cross <- function(x, y, i=1:3) {
  x <- create3D(x)
  y <- create3D(y)
  j <- function(i) (i-1) %% 3+1
  return (x[j(i+1)]*y[j(i+2)] - x[j(i+2)]*y[j(i+1)])
}

#' Identify the center of a small circle based on three points of a sphere
#'
#' @param x A matrix of 3 points, either polar longitude and latitude coordinates, or XYZ Cartesian coordinates.
#' @param origin The origin of the sphere.
#' @param output If set to "polar" then the function will return the longitude and latitude of the small circle center.
#' @return A list with two elements, the center o the small circle and the surface radius of the small circle.
#' @examples
#' # generate 3 points on a sphere
#' ps <- icosa::rpsphere(3, output="polar")
#' small <- circumcenter(x=ps)
#' circle<- smallcircles(x=small$center, r=small$r, output="polar")
#'
#' plot(NULL, NULL, xlim=c(-180, 180), ylim=c(-90, 90))
#' points(ps, col="red", pch=3, cex=4)
#' points(small$center, col="blue", pch=16)
#' points(circle, col="green", pch=16)
#' @export
circumcenter <- function(x, origin=c(0,0,0), output="polar"){

#	origin <- c(0, 0,0)
	# translate the Cartesian
	if(ncol(x)==2){
		three <- icosa::PolToCar(x)
	}else{
		three <- x
	}


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

authRadius <- icosa:::authRadius

# create small circle
standardSmallCircle <- function(rad, n, radius=authRadius){
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
#' @param x Numeric matrix or character string. Coordinates of the centers of the small circles
#' @param n Integer or NULL. only used if "x= random"
#' @param breaks Integer. The number of points to create from the small circle.
#' @param r The surface radius of the small circle expressed as great circle distance of the small circle's points from the center of the small circles.
#' @param r.ex The extrinsic radius of the small circles, i.e. the radius of the circle in the plane of the small circle.
#' @param r.rad The angle of the small circle in radians.
#' @param r.deg The angle of the small circle in degrees.
#' @param radius The radius of the sphere, defaults to the authalic radius of Earth.
#' @param origin The center of the sphere (don't touch this, unless you know what you are doing!).
#' @param output Output structure for the function. The value output="polar" will return the polar (longitude-latitude) coordinates of the small circles in an array.
#' @param drop If there is a single small circle to generate, should its array wrapper be dropped?
#' The setting output="cartesian" will return the 3D cartesian coordinates of the small circles in an array. The option output="sf" will return an sfc geometry collection.
#' @return Either a numeric array or and sfc geometry collection.
#' @export
smallcircles <- function(x="random", r=NULL, r.ex=NULL, r.rad=NULL, r.deg=NULL, n=NULL, breaks=100,radius=authRadius, origin=c(0,0,0), output="polar", sf.type="polygon", sf.wrap.dateline=TRUE, drop=TRUE){

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

	# 1. create a standard small circle
	stand <- standardSmallCircle(r.rad,  n = breaks)

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
		if(is.null(n)) stop("You must provide the number of small circles to generate.")
		if(length(n)!=1 | !is.numeric(n)) stop("Provide a single integer 'n'.")
		if(length(n)<1 | n%%1!=0) stop("Provide a single integer 'n'.")


		# to hold the data
		results <- array(NA, dim=c(breaks, length(coordNames), n))
		dimnames(results)[[2]] <- coordNames

		# rotate the small circles
		for(i in 1:n){
			results[,,i]<- rotate(stand, angles="random", pivot=origin, output=gen)
		}

	}

	if(inherits(x, "numeric")) if(is.null(dim(x))) x <- matrix(x, nrow=1)
	if(inherits(x, "matrix")){
		# make sure it is longitude latitude
		if(ncol(x)==3) x <- CarToPol(x, norad=TRUE)

		# the number of points to generate
		n <- nrow(x)

		# to hold the data
		results <- array(NA, dim=c(breaks, length(coordNames), n))
		dimnames(results)[[2]] <- coordNames

		for(i in 1:n){
			results[,,i]<- rotate(stand, long=x[i,1], lat=x[i,2], reflong=0, pivot=origin, output=gen)
		}
		# drop array
		if(n==1 & drop){
			results <- results[,,1]
		}

	}

	# if the output is an sf object
	if(output=="sf"){
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
