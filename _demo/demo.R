# Demo script to show the shape descriptors in orange
# Á.T. Kocsis 2026-07-29
# CC-BY

library(orange)

########################################
# Example 1. Data from OBIS - nice behavior!
########################################
# 1. Canvas
hex <- hexagrid(spacing=2, sf=TRUE)
plot(hex, reset=FALSE, xlim=c(-15, 40), ylim=c(25, 63), border="gray70")

# 2. Records
data(pinna)

# focus on one species
nobilis <- pinna[which(pinna$species=="Pinna nobilis"), ]

# just the coordinates, default column names
coords <- nobilis[, c("decimalLongitude", "decimalLatitude")]
colnames(coords) <- c("long", "lat")
points(coords, pch=3, col="darkgreen")

# 3. calculate and visualize
cent <- centroid(coords, plot=TRUE)

# secondary visualization
# points(cent[1], cent[2])

########################################
# grid cell occupancy
occ <- occupancy(coords, s=hex, plot=TRUE)
occ
occ <- occupancy(coords, s=hex, full=TRUE)

# secondary visualization
# plot(hex, occ$cells, add=TRUE, col="green")

########################################
# maximum great circle distance cell occupancy
mgcd <- mgcd(coords, plot=TRUE)
mgcd
mgcd <- mgcd(coords, full=TRUE)

# single line visualization
# arcs(coords[mgcd$index, ])

########################################
# Centroid radius (without confidence)
centDist <- radius(coords, plot=TRUE)
centDist

centDist <- radius(coords, full=TRUE)

# lines(rbind(centDist$centroid,
# 		coords[which(centDist$estimate==centDist$distances)[1],]),
# 		col="blue", lwd=3)

########################################
# the latitudinal range
latrange <- latrange(coords, plot=TRUE)
latrange 
latrange <- latrange(coords, full=TRUE)

# abline(h=latrange$range, lty=2)

########################################
mst <- mstlength(coords, plot=TRUE, plot.args=list(col="green", lwd=3))
mst
mst <- mstlength(coords, full=TRUE)
# lines(mst$show)
# arcs(mst$show)


########################################
# Example 2. Procuedural data highlighting some issues
########################################
########################################
# coarser grid
hex <- hexagrid(deg=5, sf=TRUE)

# controls only the centroid, but not the actual distribution generation
data(kentsamples)
coords <- kentsamples$dateline_l
colnames(coords) <- c("long", "lat")

plot(hex, border="gray")
points(coords, pch=3, col="blue")


# centroids
cent <- centroid(coords, plot=TRUE)


# secondary visualization
# points(cent[1], cent[2])

########################################
# grid cell occupancy
occ <- occupancy(coords, s=hex, plot=TRUE)

# plot(hex, occ$cells, add=TRUE, col="green")

########################################
# gappiness
gap <- gappiness(coords, s=hex, full=TRUE)

plot(hex, names(gap$holes), col="#FFFF0044", add=TRUE)

########################################
# maximum great circle distance cell occupancy
# this can create interesting results
mgcd <- mgcd(coords, plot=TRUE)

# single line visualization
# arcs(coords[mgcd$index, ])

# lines(rbind(centroidRadius$centroid,
# 		coords[which(centroidRadius$estimate==centroidRadius$distances)[1],]),
# 		col="blue", lwd=3)

#####################################
# the latitudinal range
latrange <- latrange(coords, plot=TRUE)

# abline(h=latrange$range, lty=2)

########################################
mst <- mstlength(coords, full=TRUE)
arcs(mst$show, col="green")

########################################
hu <- shull(coords, method="centroidcircle")
plot(hu, add=TRUE, col="orange")
shullarea(hu, metric="prop")
