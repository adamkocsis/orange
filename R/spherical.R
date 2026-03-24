
# Calculate new coordinate(s) from original point(s) azimuth(s) and distances(km)
# Assumes a spherical model!
#
# @param coords 2column matrix: longitudes and latitudes
# @azimuth azimuth in degrees
# @distance Distances(s) in kilometers
# @radius radius of sphere in kilometers
# PointFromVector(cbind(long = -82.21, lat = 33.3362), azimuth = 145.9, distance = 5)
CoordFromDistAzimuth <- function(coords, azimuth, distance, radius=6371){
	long <- coords[,1 ]
	lat <- coords[,2]

## lat <- 33.3362
## long <- -82.21
## azimuth <- 145.9
## distance <- 5

	# should give  33.2989617271686
	newlat <- asin(sin(lat*pi/180)*cos(distance/radius)+cos(lat*pi/180)*sin(distance/radius)*cos(azimuth*pi/180))*180/pi

	# wrong formula!!
	# shoudl give  -82.1798382123285
	## newlong <- (long*pi/180+atan2(
	## 	sin(azimuth*pi/180)*sin(distance/radius)*cos(newlat*pi/180),
	## 	cos(distance/radius)-sin(newlat*pi/180)*sin(newlat*pi/180)))*180/pi

	dlong <- atan2(sin(azimuth*pi/180)*sin(distance/radius)*cos(lat*pi/180),cos(distance/radius)-sin(lat*pi/180)*sin(newlat*pi/180))
	newlong <- ((long*pi/180-dlong + pi)%%(2*pi) - pi) /pi *180

	return(cbind(long=newlong, lat=newlat))
}

# Function to calcualte the Euler rotation of a point
#
# Interface improvements over fossil::euler
# @param x A matrix of longitude-latitude coordinates
# @param rotation  A single rotation value, or a vector of rotations values
# @param pole A matrix of longtidue-latitude coordinates
# @examples
## point <- rpsphere(1, output="polar")
## pole <- rpsphere(1, output="polar")
## temp <- eulerRotation(x=point, rotation=0:360, pole=pole)
## plot(hex)
## points(pole, col="red", pch=16, cex=3)
## points(point, col="blue", pch=16)
## points(temp)
eulerRotation <- function(x, rotation, pole){

	# lat2, long2
	# The rotated point
	pointLong <- x[,1]*pi/180
	pointLat <- x[,2]*pi/180

	# the rotation
	rotation <- rotation*pi/180

	# The euler pole
	poleLong <- pole[,1]*pi/180
	poleLat <- pole[,2]*pi/180

	# differences
	diffLong <- pointLong-poleLong
	diffLat <- pointLat-poleLat

	# latitude
	cosPoleLat <- cos(poleLat)
	cosPointLat <- cos(pointLat)
	a<-(sin(diffLat/2))^2+cosPoleLat*cosPointLat*(sin(diffLong/2))^2
	d<-2*atan2(sqrt(a),sqrt(1-a))

	# bearing
	sinPointLat <- sin(pointLat)
	sinPoleLat <- sin(poleLat)
	bear<-atan2(sin(diffLong)*cosPointLat,cosPoleLat*sinPointLat-sinPoleLat*cosPointLat*cos(diffLong))
	tc<-bear-rotation


	# new Coordinates
	sd<-sin(d)
	cd<-cos(d)
	newLat <-asin(sinPoleLat*cd+cosPoleLat*sd*cos(tc))
	ndlon<-atan2(sin(tc)*sd*cosPoleLat,cd-sinPoleLat*sin(newLat))
	newLong <-((poleLong+ndlon+pi)%%(2*pi))-pi

	# return as a matrix
	return(cbind(newLong/pi*180, newLat/pi*180))

}




# Calculate points along a great circle with given distance from a source point
#
# @param sourcePoint A matrix with a source point
# @param dists A vector of distances
# @param otherPoint Another point defining the great circle
# @example
## sourcePoint <- rpsphere(1, output="polar")
## otherPoint <- rpsphere(1, output="polar")
## dists <- seq(500, 8000, 500)
## plot(hex)
## points(sourcePoint, col="red", pch=16, cex=3)
## points(otherPoint, col="blue", pch=16)
## temp <- PushFromPoint(sourcePoint, dists, otherPoint)
## points(temp)
## icosa::arcDist(sourcePoint, temp[1, ])
pushFromPoint <- function(sourcePoint, dists, otherPoint, output="polar"){
		# transform both to Cartesian
		if(ncol(sourcePoint)==2){
			sourceCart <- icosa::PolToCar(sourcePoint)
		}else{
			sourceCart <- sourcePoint
		}
		if(ncol(otherPoint)==2){
			otherCart <- icosa::PolToCar(otherPoint)
		}else{
			otherCart <- otherPoint
		}

		# b and a should be the same!
		bLength <- sqrt(sum(sourceCart^2))
		cLength <- apply(otherCart, 1, function(x){ sqrt(sum(x^2))})
		aLength <- apply(otherCart, 1, function(x){ sqrt(sum((sourceCart-x)^2)) })
#		aLength <- sqrt(sum((sourceCart-otherCart)^2))

		# the angle between the points
		alpha <- acos((bLength^2+cLength^2-aLength^2)/(2*bLength*cLength))

		# translate distance to radian
		theta <- dists/bLength

		# the new coordinates
		newX <- (sin(alpha-theta)*sourceCart[,1]+sin(theta)*otherCart[,1])/sin(alpha)
		newY <- (sin(alpha-theta)*sourceCart[,2]+sin(theta)*otherCart[,2])/sin(alpha)
		newZ <- (sin(alpha-theta)*sourceCart[,3]+sin(theta)*otherCart[,3])/sin(alpha)

		# combine transform and return
		res <- cbind(newX, newY, newZ)
		if(output=="polar") res <- icosa::CarToPol(res)[, 1:2]

		return(res)

}


