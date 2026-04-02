library(tinytest)
library(orange)
suppressPackageStartupMessages(library(divDyn))

# the data
data(pinna)
data(corals)

# mke a hexagrid
suppressMessages(hex <- hexagrid(deg=2, sf=TRUE))
################################################################################
# Base:

# Matrix method
mat <- as.matrix(pinna[, c("decimallongitude", "decimallatitude")])
cellsManual <- unique(locate(hex, mat))
manual <- length(cellsManual)

# actual tests
expect_silent(occ <- occupancy(mat, hex))
expect_equal(occ, manual)

expect_silent(occFull <- occupancy(mat, hex, full=TRUE))
expect_true(inherits(occFull, "orange"))
expect_equal(occFull$cells, cellsManual)
expect_equal(occFull$estimate, occ)

# should be the same with duplicates
expect_silent(occDup <- occupancy(mat, hex, duplicates=TRUE))
expect_equal(occDup, occ)

# should give error if columns are partially given
expect_error(occupancy(mat, hex, lat="decimallatitude"))
expect_error(occupancy(mat, hex, long="decimallongitude"))
expect_error(occupancy(mat, hex, long="decimallongitude", lat="wrong"))
expect_silent(occNamed <- occupancy(mat, hex, long="decimallongitude", lat="decimallatitude"))
expect_equal(occNamed, occ)

# for now q is not yet implemented
expect_error(occupancy(mat, hex, q=0.7))


# error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE))


# normal plotting
plot(hex, border=NA, col=NA)
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE))

# Single occurrences
# actual tests
first <- mat[1,, drop=FALSE]
expect_silent(sing <- occupancy(first, hex, full=TRUE))
# with explicitly given names (dropping problem)
expect_silent(sing <- occupancy(first, icosa=hex, long="decimallongitude", lat="decimallatitude"))
# false names
empty <- matrix(c(NA,NA), nrow=1)
expect_error(none <- occupancy( empty, icosa=hex, long="decimallongitude", lat="decimallatitude"))
# no coordinates
expect_silent(none <- occupancy(empty , icosa=hex))
expect_equal(none, 0L)
# proper omission
expect_silent(singNone <- occupancy(rbind(empty, first), icosa=hex))
expect_equal(singNone, sing)



################################################################################
# Data.frame -  appropriate fallback to matrix method
expect_silent(occDF <- occupancy(pinna, icosa=hex, long="decimallongitude", lat="decimallatitude"))
expect_equal(occDF, occ)
expect_silent(occDFfull <- occupancy(pinna, icosa=hex, long="decimallongitude", lat="decimallatitude", full=TRUE))
expect_equal(occDFfull, occFull)

# appropriate defaults for coordinates - 2column df
df <- as.data.frame(mat)
colnames(df) <- c("long", "lat")
expect_silent(occDFdefault <- occupancy(df, icosa=hex))
expect_equal(occDFdefault, occ)

################################################################################
# Data.frame - taxon iteration
# the corals
corCell<- locate(hex, corals[, c("lng", "lat")])
corOccManual <- tapply(INDEX=corals$genus, X=corCell, function(x) length(levels(factor((x)))))

# somewhat slow!
expect_silent(corOcc <- occupancy(x=corals, icosa=hex, long="lng", lat="lat", tax="genus" ))
expect_equal(corOccManual, corOcc)

# The full output
expect_silent(corOccList <- occupancy(x=corals, icosa=hex, long="lng", lat="lat", tax="genus", full=TRUE ))
expect_true(is.list(corOccList))
