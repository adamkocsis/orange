library(tinytest)
library(orange)
suppressPackageStartupMessages(library(divDyn))

# the data
data(pinna)
data(corals)

################################################################################
# Base:

# Matrix method
mat <- as.matrix(pinna[, c("decimallongitude", "decimallatitude")])
manual <- max(icosa::arcdistmat(mat))


# actual tests
expect_silent(md <- maxdist(mat))
expect_equal(md, manual)

expect_silent(mdFull <- maxdist(mat, full=TRUE))
expect_true(inherits(mdFull, "orange"))
expect_equal(mdFull$estimate, manual)
expect_true(inherits(mdFull$index, "integer"))
indices<- as.numeric(mdFull$index)
remanual <- arcdist(
	mat[indices[1], ],
	mat[indices[2], ]
)
expect_equal(mdFull$estimate, remanual)


# should be the same with duplicates, add a duplicate manually for testing
matDupl <- rbind(mat[1,], mat)
expect_silent(mdDup <- maxdist(matDupl, duplicates=TRUE))
expect_equal(mdDup, md)
expect_silent(mdDupFull <- maxdist(matDupl, duplicates=TRUE, full=TRUE))
expect_true(inherits(mdDupFull$index, "integer"))
indices<- as.numeric(mdDupFull$index)
remanual2 <- arcdist(
	matDupl[indices[1], ],
	matDupl[indices[2], ]
)
expect_equal(mdFull$estimate, remanual2)



# should give error if columns are partially given
expect_error(maxdist(mat, lat="decimallatitude"))
expect_error(maxdist(mat, long="decimallongitude"))
expect_error(maxdist(mat, long="decimallongitude", lat="wrong"))
expect_silent(mdNamed <- maxdist(mat, long="decimallongitude", lat="decimallatitude"))
expect_equal(mdNamed, md)

# for now q is not yet implemented
expect_error(maxdist(mat, q=0.7))


# no error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(mdPlot <- maxdist(mat, plot=TRUE))


# normal plotting
plot(mat)
expect_silent(mdPlot <- maxdist(mat, plot=TRUE))

# Single occurrences
# actual tests
first <- mat[1,, drop=FALSE]
expect_silent(sing <- maxdist(first, full=TRUE))
expect_equal(sing$estimate, 0)
# with explicitly given names (dropping problem)
expect_silent(sing <- maxdist(first, long="decimallongitude", lat="decimallatitude"))
# false names
empty <- matrix(c(NA,NA), nrow=1)
expect_error(none <- maxdist( empty, long="decimallongitude", lat="decimallatitude"))
# no coordinates
expect_silent(none <- maxdist(empty))
# proper omission
expect_silent(singNone <- maxdist(rbind(empty, first)))
expect_equal(singNone, sing)


################################################################################
# Manually calculated distance matrix
dmManual <- arcdistmat(mat)
expect_silent(mdDM <- maxdist(x=mat, dm=dmManual))
expect_equal(mdDM, md)

#  



################################################################################
# Data.frame -  appropriate fallback to matrix method
expect_silent(mdDF <- maxdist(pinna, long="decimallongitude", lat="decimallatitude"))
expect_equal(mdDF, md)
expect_silent(mdDFfull <- maxdist(pinna, long="decimallongitude", lat="decimallatitude", full=TRUE))
expect_equal(mdDFfull, mdFull)

# appropriate defaults for coordinates - 2column df
df <- as.data.frame(mat)
colnames(df) <- c("long", "lat")
expect_silent(mdDFdefault <- maxdist(df))
expect_equal(mdDFdefault, md)

################################################################################
# Data.frame - taxon iteration
# the corals
gen <- levels(factor(corals$genus))
mdTaxManual <- rep(NA, length(gen))
mdTaxListManual <- list()
for(i in 1:length(gen)){
	# the slice
	thisSlice <- corals[which(corals$genus==gen[i]), ]
	expect_silent(one <- maxdist(thisSlice, long="lng", lat="lat", full=TRUE))
	if(!is.na(one$estimate)){
		oneManual <- arcdist(
			as.matrix(thisSlice[one$index[1], c("lng", "lat")]),
			as.matrix(thisSlice[one$index[2], c("lng", "lat")])
		)
		expect_equal(one$estimate, oneManual)
	}
	mdTaxManual[i] <- one$estimate
	mdTaxListManual[[i]] <- one
	## cat(i, "\r")
	## flush.console()
}
dim(mdTaxManual) <- length(mdTaxManual) 
names(mdTaxManual) <- gen
dim(mdTaxListManual) <- length(mdTaxListManual) 
names(mdTaxListManual) <- gen

# somewhat slow!
expect_silent(mdTax <- maxdist(x=corals, long="lng", lat="lat", tax="genus" ))
expect_equal(mdTax, mdTaxManual)

################################################################################
# Not yet implemented
################################################################################

# The full output - not yet
expect_error(maxdist(x=corals, long="lng", lat="lat", tax="genus", full=TRUE ))
# The indices of the output needs to be properly set
## expect_silent(mdTaxList <- maxdist(x=corals, long="lng", lat="lat", tax="genus", full=TRUE ))
## expect_true(is.list(mdTaxList))
## expect_equal(mdTaxList, mdTaxListManual)


# Full output with a given distance matrix
## corDist <- arcdistmat(corals[, c("lng", "lat")])
## expect_error(maxdist(x=corals, dm=corDist, long="lng", lat="lat", tax="genus", full=TRUE ))
