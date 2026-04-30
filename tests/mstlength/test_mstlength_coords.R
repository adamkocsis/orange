# in data frame
# mstlength: x: matrix/data.frame, s: missing
library(tinytest)
suppressPackageStartupMessages(library(divDyn))
suppressPackageStartupMessages(library(vegan))


################################################################################
# 1. Single-taxon dataset
################################################################################
data(pinna)

# the matrix
mat <- as.matrix(pinna[,c("decimallongitude", "decimallatitude")])

# calculate metric manually (no duplicates)
unmat <- unique(mat)
dm <- icosa::arcdistmat(unmat)
stre <- vegan::spantree(as.dist(dm))
manual <- sum(stre$dist)

# simple matrix solution
expect_silent(mstMat <- mstlength(mat))
expect_equal(manual, mstMat)

# with full outpu
expect_silent(mstMatFull <- mstlength(as.matrix(pinna[,c("decimallongitude", "decimallatitude")]), full=TRUE))
expect_equal(names(mstMatFull), c("estimate", "index", "show"))
expect_equal(mstMat, mstMatFull$estimate)

# duplicates make no difference
expect_silent(mstMatDupl <- mstlength(as.matrix(pinna[,c("decimallongitude", "decimallatitude")]), duplicates=TRUE))
expect_equal(mstMatDupl, mstMat)


# Single instances
# actual tests
first <- mat[1,, drop=FALSE]
expect_silent(sing <- mstlength(first, full=TRUE))
expect_equal(sing$estimate, 0)
# with explicitly given names (dropping problem)
expect_silent(sing <- mstlength(first, long="decimallongitude", lat="decimallatitude"))
expect_equal(sing, 0)

# false names
empty <- matrix(c(NA,NA), nrow=1)
expect_error(none <- mstlength( empty, long="decimallongitude", lat="decimallatitude"))
# no coordinates
expect_silent(none <- mstlength(empty))
# proper omission
expect_silent(singNone <- mstlength(rbind(empty, first)))
expect_equal(singNone, sing)


# missing values
# same as with missing values
matMissFront <- rbind(c(NA,NA), mat)
expect_silent(mstFront <- mstlength(matMissFront, long="decimallongitude", lat="decimallatitude"))
expect_equal(mstMat,mstFront)

matMissBack <- rbind(mat, c(NA,NA))
expect_silent(mstBack <- mstlength(matMissBack, long="decimallongitude", lat="decimallatitude"))
expect_equal(mstMat,mstBack)

matMissMid <- rbind(mat[1:3,], c(NA,NA), mat[4:nrow(mat),])
expect_silent(mstMid <- mstlength(matMissMid, long="decimallongitude", lat="decimallatitude"))
expect_equal(mstMat,mstMid)

# only missing, default names
expect_silent(mstMiss <- mstlength(matrix(NA, ncol=2, nrow=1)))
expect_equal(mstMiss,NA )


# Plotting
# expect_silent(namedlat <- mstlength(mat, long="decimallongitude", lat="decimallatitude", plot=TRUE))

################################################################################
# Data.frame method

# coordinate pairs
expect_silent(mstDF <- mstlength(pinna, long="decimallongitude", lat="decimallatitude"))
expect_equal(mstMat, mstDF)

# with full output
expect_silent(mstDFFull <- mstlength(pinna, long="decimallongitude", lat="decimallatitude", full=TRUE))
expect_equal(mstMat, mstDFFull$estimate)




################################################################################
# 2. Multi-taxon dataset
################################################################################

################################################################################
# Occupancy based on coordinates, multiple taxa
data(corals)

# occupancy in the corals dataset
expect_silent(msts <- mstlength(corals, tax="genus", long="lng", lat="lat"))
expect_equal(class(msts), "numeric")

# manually calculated
gen <- levels(factor(corals$genus))
manual <- rep(NA, length(gen))
names(manual) <- gen
for(i in 1:length(gen)){
	genus <-corals[which(corals$genus==gen[i]), ]
	# the latitudes
	if(sum(is.na(genus[, "lat"]))==nrow(genus)){
		thing <- NA
	}else{
		unmat <- unique(genus[, c("lng", "lat")])
		unmat <- na.omit(unmat)
		if(nrow(unmat)!=1){
			dm <- icosa::arcdistmat(unmat)
			stre <- vegan::spantree(as.dist(dm))
			thing <- sum(stre$dist)
		}else{
			thing <- 0
		}
	}
	manual[i] <- thing 

}

# only estimates: match?
expect_equal(manual, msts)



################################################################################
# Not implemented etc.

# check full output match
expect_error(mstsFull <- mstlegnth(corals, tax="genus", long="lng", lat="lat", full=TRUE))
expect_warning(mstPlot  <- mstlength(corals, tax="genus", long="lng", lat="lat", plot=TRUE))

expect_equal(mstPlot, msts)
expect_error(mstQ <- mstlength(mat, q=0.95))
