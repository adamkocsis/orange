library(tinytest)
library(orange)
suppressPackageStartupMessages(library(divDyn))

# the data
data(pinna)
data(corals)

################################################################################
# Base:
# Matrix method
################################################################################
mat <- as.matrix(pinna[, c("decimallongitude", "decimallatitude")])

# the focus
foc <-centroid(mat)

# calculate metric manually initials with q=1
dists <- icosa::arcdistmat(mat,  matrix(foc, ncol=2))

# the esimate
manual <- max(dists)

# the index
manIndex <- which(dists==manual)

# actual tests
expect_silent(rad <- radius(mat))
expect_equal(rad, manual)

# The Full output
expect_silent(radFull <- radius(mat, full=TRUE))
expect_true(inherits(radFull, "orange"))
expect_equal(radFull$estimate, manual)
expect_true(inherits(radFull$index, "integer"))

# test whether the index works properly
indices<- radFull$index
remanual <- arcdist(
	mat[indices[1], ],
	foc	
)
expect_equal(radFull$estimate, remanual)


# should be the same with duplicates, add a duplicate manually for testing
matDupl <- rbind(mat[1,], mat)
expect_silent(radDup <- radius(matDupl, duplicates=TRUE))
expect_equal(radDup, rad)
expect_silent(radDupFull <- radius(matDupl, duplicates=TRUE, full=TRUE))
expect_true(inherits(radDupFull$index, "integer"))

# should be the same despite the index shift!!!
indices<- radDupFull$index
remanual2 <- arcdist(
	matDupl[indices[1], ],
	foc	
)
expect_equal(radFull$estimate, remanual2)


# should give error if columns are partially given
expect_error(radius(mat, lat="decimallatitude"))
expect_error(radius(mat, long="decimallongitude"))
expect_error(radius(mat, long="decimallongitude", lat="wrong"))
expect_silent(radNamed <- radius(mat, long="decimallongitude", lat="decimallatitude"))
expect_equal(radNamed, rad)


################################################################################
# Matrix method with q!=1
################################################################################

# for now q is not yet implemented
expect_error(maxdist(mat, q=0.7))

################################################################################
# Plotting
################################################################################

# no error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(radPlot <- radius(mat, plot=TRUE))

# normal plotting
plot(mat, xlim=c(-20, 50), ylim=c(0, 70))
expect_silent(radPlot <- radius(mat, plot=TRUE))

# Single occurrences
# actual tests
first <- mat[1,, drop=FALSE]
expect_silent(sing <- radius(first, full=TRUE))
expect_equal(sing$estimate, 0)
# with explicitly given names (dropping problem)
expect_silent(sing <- radius(first, long="decimallongitude", lat="decimallatitude"))
# false names
empty <- matrix(c(NA,NA), nrow=1)
expect_error(none <- radius( empty, long="decimallongitude", lat="decimallatitude"))
# no coordinates
expect_silent(none <- radius(empty))
# proper omission
expect_silent(singNone <- radius(rbind(empty, first)))
expect_equal(singNone, sing)



################################################################################
# Data.frame -  appropriate fallback to matrix method
expect_silent(radDF <- radius(pinna, long="decimallongitude", lat="decimallatitude"))
expect_equal(radDF, rad)
expect_silent(radDFfull <- radius(pinna, long="decimallongitude", lat="decimallatitude", full=TRUE))
expect_equal(radDFfull, radFull)

# appropriate defaults for coordinates - 2column df
df <- as.data.frame(mat)
colnames(df) <- c("long", "lat")
expect_silent(radDFdefault <- radius(df))
expect_equal(radDFdefault, rad)

################################################################################
# Data.frame - taxon iteration
# the corals
gen <- levels(factor(corals$genus))
radTaxManual <- rep(NA, length(gen))
radTaxListManual <- list()
for(i in 1:length(gen)){
	# the slice
	thisSlice <- corals[which(corals$genus==gen[i]), ]
	expect_silent(one <- radius(thisSlice, long="lng", lat="lat", full=TRUE))
	if(!is.na(one$estimate)){
		oneManual <- arcdist(
			as.matrix(thisSlice[one$index[1], c("lng", "lat")]),
			matrix(centroid(thisSlice[, c("lng", "lat")], long="lng", lat="lat"), ncol=2)
		)
		expect_equal(one$estimate, oneManual)
	}
	radTaxManual[i] <- one$estimate
	radTaxListManual[[i]] <- one
	## cat(i, "\r")
	## flush.console()
}
names(radTaxManual) <- gen
dim(radTaxListManual) <- length(radTaxListManual) 
names(radTaxListManual) <- gen

# somewhat slow!
expect_silent(radTax <- radius(x=corals, long="lng", lat="lat", tax="genus" ))
expect_equal(radTax, radTaxManual)

## ################################################################################
## # Not yet implemented
## ################################################################################

# The full output - not yet
expect_error(radius(x=corals, long="lng", lat="lat", tax="genus", full=TRUE ))
# The indices of the output needs to be properly set
## expect_silent(mdTaxList <- maxdist(x=corals, long="lng", lat="lat", tax="genus", full=TRUE ))
## expect_true(is.list(mdTaxList))
## expect_equal(mdTaxList, mdTaxListManual)


# Full output with a given distance matrix
## corDist <- arcdistmat(corals[, c("lng", "lat")])
## expect_error(maxdist(x=corals, dm=corDist, long="lng", lat="lat", tax="genus", full=TRUE ))
