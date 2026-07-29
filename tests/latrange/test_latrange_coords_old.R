# in data frame
# occupancy: x: data.frame, sf: missing, icosa: missing
library(tinytest)
suppressPackageStartupMessages(library(divDyn))


################################################################################
# 1. Single-taxon dataset
################################################################################
data(pinna)

manual <- diff(range(pinna$decimallatitude))
mat <- as.matrix(pinna[,c("decimallongitude", "decimallatitude")])
# simple matrix solution
expect_silent(ranMat <- latrange(mat))
expect_equal(manual, ranMat)

# with full outpu
expect_silent(ranMatFull <- latrange(as.matrix(pinna[,c("decimallongitude", "decimallatitude")]), full=TRUE))
expect_equal(diff(ranMatFull$range), ranMatFull$estimate)
expect_equal(ranMat, ranMatFull$estimate)

# duplicates make no difference
expect_silent(ranMatDupl <- latrange(as.matrix(pinna[,c("decimallongitude", "decimallatitude")]), duplicates=TRUE))
expect_equal(ranMatDupl, ranMat)


# Single instances
# actual tests
first <- mat[1,, drop=FALSE]
expect_silent(sing <- latrange(first, full=TRUE))
expect_equal(sing$estimate, 0)
# with explicitly given names (dropping problem)
expect_silent(sing <- latrange(first, long="decimallongitude", lat="decimallatitude"))
expect_equal(sing, 0)

# false names
empty <- matrix(c(NA,NA), nrow=1)
expect_error(none <- latrange( empty, long="decimallongitude", lat="decimallatitude"))
# no coordinates
expect_silent(none <- latrange(empty))
# proper omission
expect_silent(singNone <- latrange(rbind(empty, first)))
expect_equal(singNone, sing)


# missing values
# same as with missing values
matMissFront <- rbind(c(NA,NA), mat)
expect_silent(latFront <- latrange(matMissFront, long="decimallongitude", lat="decimallatitude"))
expect_equal(ranMat,latFront)

matMissBack <- rbind(mat, c(NA,NA))
expect_silent(latBack <- latrange(matMissBack, long="decimallongitude", lat="decimallatitude"))
expect_equal(ranMat,latBack)

matMissMid <- rbind(mat[1:3,], c(NA,NA), mat[4:nrow(mat),])
expect_silent(latMid <- latrange(matMissMid, long="decimallongitude", lat="decimallatitude"))
expect_equal(ranMat,latMid)

# only missing, default names
expect_silent(latMiss <- latrange(matrix(NA, ncol=2, nrow=1)))
expect_equal(latMiss,as.numeric(NA) )


# Plotting
# expect_silent(namedlat <- latrange(mat, long="decimallongitude", lat="decimallatitude", plot=TRUE))

################################################################################
# Data.frame method

# coordinate pairs
expect_silent(ranDF <- latrange(pinna, long="decimallongitude", lat="decimallatitude"))
expect_equal(ranMat, ranDF)

# with full output
expect_silent(ranDFFull <- latrange(pinna, long="decimallongitude", lat="decimallatitude", full=TRUE))
expect_equal(diff(ranDFFull$range), ranDFFull$estimate)
expect_equal(ranMat, ranDFFull$estimate)




################################################################################
# 2. Multi-taxon dataset
################################################################################

################################################################################
# Occupancy based on coordinates, multiple taxa
data(corals)

# occupancy in the corals dataset
expect_silent(rans <- latrange(corals, tax="genus", long="lng", lat="lat"))
expect_equal(class(rans), "numeric")

# manually calculated
gen <- levels(factor(corals$genus))
manual <- rep(NA, length(gen))
ran_list <- list()
names(manual) <- gen
for(i in 1:length(gen)){
	genus <-corals[which(corals$genus==gen[i]), ]
	# the latitudes
	lats <- genus[, "lat"]
	if(sum(is.na(lats))==length(lats)){
		thing <- c(NA, NA)
	}else{
		thing <- range(lats, na.rm=TRUE)
	}
	manual[i] <- diff(thing)
	ran_list[[i]] <- list(estimate=as.numeric(manual[i]), range=thing)

}
names(ran_list) <- gen

# only estimates match?
expect_equal(manual, rans)

# check full output match
expect_silent(ransFull <- latrange(corals, tax="genus", long="lng", lat="lat", full=TRUE))
expect_equal(ransFull, ran_list)


################################################################################
# Not implemented etc.

expect_warning(ranPlot  <- latrange(corals, tax="genus", long="lng", lat="lat", plot=TRUE))
expect_equal(ranPlot, rans)
expect_error(ranQ <- latrange(mat, q=0.95))
