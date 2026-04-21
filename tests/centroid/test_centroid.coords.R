# in data frame
# centroid: x: data.frame, sf: missing, icosa: missing
library(tinytest)
library(divDyn)


################################################################################
# 1. Single-taxon dataset (matrix)
################################################################################
data(pinna)

#  matrix method
mat <- as.matrix(pinna[, c("decimallongitude", "decimallatitude")])
expect_silent(cent <- centroid(mat))

# named matrix
expect_error(cent <- centroid(mat, long="long", lat="lat"))
expect_silent(namedcent <- centroid(mat, long="decimallongitude", lat="decimallatitude"))
expect_equal(cent,namedcent)

# DEFAULT: duplicates=FALSE
expect_silent(centNoDupl <- centroid(mat,long="decimallongitude", lat="decimallatitude", duplicates=FALSE))
expect_equal(cent,centNoDupl)

# not the same as duplicates=TRUE
expect_silent(centDupl <- centroid(mat,long="decimallongitude", lat="decimallatitude", duplicates=TRUE))
expect_true(sum(abs(centDupl-centNoDupl))!=0)

# same as with missing values
matMissFront <- rbind(c(NA,NA), mat)
expect_silent(centFront <- centroid(matMissFront, long="decimallongitude", lat="decimallatitude"))
expect_equal(cent,centFront)

matMissBack <- rbind(mat, c(NA,NA))
expect_silent(centBack <- centroid(matMissBack, long="decimallongitude", lat="decimallatitude"))
expect_equal(cent,centBack)

matMissMid <- rbind(mat[1:3,], c(NA,NA), mat[4:nrow(mat),])
expect_silent(centMid <- centroid(matMissMid, long="decimallongitude", lat="decimallatitude"))
expect_equal(cent,centMid)

# only missing
expect_silent(centMiss <- centroid(matrix(NA, ncol=2, nrow=1)))
misser <- c(NA, NA)
names(misser) <- c("long", "lat")
expect_equal(centMiss,misser )


# plotting both on previous and on new
expect_silent(namedcent <- centroid(mat, long="decimallongitude", lat="decimallatitude", plot=TRUE))
expect_silent(namedcent2 <- centroid(mat, long="decimallongitude", lat="decimallatitude", plot=TRUE,
	duplicates=TRUE, plot.args=list(col="#0000BB", pch=3)))

################################################################################
# 1. Single-taxon dataset (data.frame)
################################################################################

# data.frame method
# coordinate pairs - fall back to matrix
expect_silent(centDF <- centroid(pinna, long="decimallongitude", lat="decimallatitude"))
expect_equal(cent, centDF)


################################################################################
# 2. Multi-taxon dataset
################################################################################

################################################################################
# Centroid for multiple taxa
data(corals)

gen <- levels(factor(corals$genus))

# occupancy in the corals dataset
expect_silent(cents <- centroid(corals, tax="genus", long="lng", lat="lat"))
expect_equal(class(cents), "numeric")
expect_equal(dim(cents))
# there is a value for every
expect_equal(gen, rownames(cents))

# do it manually
manual <- matrix(NA, ncol=2, nrow=length(gen))
rownames(manual) <- gen
colnames(manual) <- c("long", "lat")
for(i in 1:length(gen)){
	genus <-corals[which(corals$genus==gen[i]), ]
	manual[i, ]<- centroid(as.matrix(genus[, c("lng", "lat")]))
}
expect_equal(manual, cents)

# expect errors, when wrong long-lat is specified
expect_error(centroid(corals))

# the same with duplicates=TRUE
expect_silent(centsDupl <- centroid(corals, tax="genus", long="lng", lat="lat", duplicates=TRUE))

# should give a different result!!!
expect_true(sum(abs(cents-centsDupl), na.rm=TRUE)>0)

#################################################################################
# Known issues

# Multi-taxon plotting is not there yet
expect_warning(centsDupl <- centroid(corals, tax="genus", long="lng", lat="lat", duplicates=TRUE, plot=TRUE))
