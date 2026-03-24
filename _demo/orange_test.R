# Demo script to show the range estimator plans in orange
# Á.T. Kocsis 2025-03-05
# CC-BY

library(orange)

########################################
# Example 1. Data from OBIS - nice behavior!
########################################
# 1. Canvas
hex <- hexagrid(deg=5, sf=TRUE)
plot(hex, reset=FALSE, xlim=c(-15, 40), ylim=c(25, 63))

# 2. Records
data(pinna)

# just the coordinates
coordMat <- SimpleCoordinates(pinna, long="decimallongitude", lat="decimallatitude")
points(coordMat)

# 3. calculate and visualize
cent <- centroid_points(coordMat, plot=TRUE)

# secondary visualization
# points(cent[1], cent[2])

########################################
# grid cell occupancy
occ <- occupancy(coordMat, icosa=hex, plot=TRUE)

# plot(hex, occ$cells, add=TRUE, col="green")

########################################
# maximum great circle distance cell occupancy
mgcd <- mgcd(coordMat, plot=TRUE)
# single line visualization
# arcs(coordMat[mgcd$index, ])

########################################
# Centroid radius (withouth confidence)
centDist <- cenrad(coordMat, plot=TRUE, q=1)

# lines(rbind(centDist$centroid,
# 		coordMat[which(centDist$estimate==centDist$distances)[1],]),
# 		col="blue", lwd=3)

########################################
# the latitudinal range
latrange <- latrange(coordMat, plot=TRUE)

# abline(h=latrange$range, lty=2)

########################################
mst <- mstlength(coordMat, plot=TRUE)
# lines(mst$show)
# arcs(mst$show)

########################################
planarChull <- chullplane(coordMat, gcbound=FALSE, plot=TRUE)
# plot(planarChull$sf, col="#55000033", add=TRUE)

planarChullGC <- chullplane(coordMat, gcbound=TRUE, plot=TRUE)
# plot(planarChullGC$sf, col="#00550033", add=TRUE)

########################################
sphericalChull <- chullsphere(coordMat, plot=TRUE)
# plot(sphericalChull$sf, col="#55000033", add=TRUE)


########################################
# Example 2. Procuedural data highlighting some issues
# These examples rely on the 'biodome' extension
########################################
########################################
# coarser grid
hex <- hexagrid(deg=5, sf=TRUE)

# controls only the centroid, but not the actual distribution generation
set.seed(5)
coordMat <- biodome::presences(100, kappa=5)
colnames(coordMat) <- c("long", "lat")

plot(hex, border="gray")
points(coordMat, pch=3, col="blue")


# centroids
cent <- centroid_points(coordMat, plot=TRUE)


# secondary visualization
# points(cent[1], cent[2])

########################################
# grid cell occupancy
occ <- occupancy(coordMat, icosa=hex, plot=TRUE)

# plot(hex, occ$cells, add=TRUE, col="green")

########################################
# maximum great circle distance cell occupancy
# this can create interesting results
mgcd <- mgcd(coordMat, plot=TRUE)

# single line visualization
# arcs(coordMat[mgcd$index, ])

########################################
# Centroid distances
centroidRadius <- cenrad(coordMat, plot=TRUE, q=1)

# lines(rbind(centroidRadius$centroid,
# 		coordMat[which(centroidRadius$estimate==centroidRadius$distances)[1],]),
# 		col="blue", lwd=3)

#####################################
# the latitudinal range
latrange <- latrange(coordMat, plot=TRUE)

# abline(h=latrange$range, lty=2)

########################################
mst <- mstlength(coordMat, plot=TRUE)
arcs(mst$show)

########################################
# What does not work!
########################################
sf::sf_use_s2(FALSE)
planarChull <- chullplane(coordMat, gcbound=FALSE, plot=TRUE)
# plot(planarChull$sf, col="#55000033", add=TRUE)

# totally fucked up
planarChullGC <- chullplane(coordMat, gcbound=TRUE, plot=TRUE)
# plot(planarChullGC$sf, col="#00550033", add=TRUE)

########################################
sf::sf_use_s2(TRUE)
sphericalChull <- chullsphere(coordMat, plot=TRUE)

# plot(sphericalChull$sf, col="#55000033", add=TRUE)
