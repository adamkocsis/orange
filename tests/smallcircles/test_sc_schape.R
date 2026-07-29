library(orange)
library(tinytest)


################################################################################
# I. Interface tests
################################################################################

# NO input
expect_error(sc_shape())
#-------------------------------------------------------------------------------
# 1. random circle generation
#-------------------------------------------------------------------------------

set.seed(0)
expect_silent(oneSmall <- sc_shape(r=100))

# One random circle
set.seed(0)
expect_silent(one <- sc_shape(n=1, r=100))

# output structure
expect_true(inherits(one, "matrix"))
expect_equal(dim(one), c(100, 2))

# reproducibility  
set.seed(0)
expect_silent(oneAgain <- sc_shape(n=1, r=100))
expect_equal(one, oneAgain)

# setting breaks
expect_silent(oneFine <- sc_shape(n=1, r=100, breaks=1000))
expect_true(inherits(oneFine, "matrix"))
expect_equal(dim(oneFine), c(1000, 2))

# coordinate type
expect_true(-180<=min(oneFine[,"long"]))
expect_true(180>=max(oneFine[,"long"]))
expect_true(-90<=min(oneFine[,"lat"]))
expect_true(90>=max(oneFine[,"lat"]))

# changing to cartesian
set.seed(0)
expect_silent(oneCart <- sc_shape(n=1, r=100, output="cartesian"))
expect_true(inherits(oneCart, "matrix"))
expect_equal(dim(oneCart), c(100, 3))
expect_equal(CarToPol(oneCart[,])[, 1:2], one)

# do not drop array wrapper!
set.seed(0)
expect_silent(oneArr <- sc_shape(n=1, r=100, drop=FALSE))
expect_equal(one, oneArr[,,1])


#-------------------------------------------------------------------------------
# 2. Multiple circle generation
#-------------------------------------------------------------------------------

# same radius with different origins
expect_silent(fifty200 <- sc_shape(n=50, r=200))
expect_true(inherits(fifty200, "array"))
expect_equal(dim(fifty200), c(100, 2,50))

# visual check!
## plot(NULL, xlim=c(-180, 180), ylim=c(-90,90), axes=FALSE, xlab="", ylab="")
## for(i in 1:(dim(fifty200)[3])){
## 	points(fifty200[,,i], col=i, pch=16)
## }


expect_silent(diff200 <- sc_shape(r=runif(200, 100, 10000)))
expect_true(inherits(diff200, "array"))
expect_equal(dim(diff200), c(100, 2,200))

# visual check!
## plot(NULL, xlim=c(-180, 180), ylim=c(-90,90), xlab="", ylab="", axes=F)
## for(i in 1:(dim(diff200)[3])){
## 	points(diff200[,,i], col=i,pch=16)
## }


set.seed(0)
cent <- rpsphere(1, output="polar")
expect_silent(outward <- sc_shape(cent, r=seq(100, 15000, length.out=20)))

# visual check!
## plot(NULL, xlim=c(-180, 180), ylim=c(-90,90), xlab="", ylab="", axes=F)
## for(i in 1:(dim(outward)[3])){
## 	points(outward[,,i], col=i,pch=16)
## }



################################################################################
# II. Stress test - are the points actually as far from the origin, as intended?
################################################################################

# accuracy is not excellent due to rounding of sine-cosine transformations..
epsilon <- 1e-3 ## (10 m level)

# generate random points
set.seed(0)
rand <- rpsphere(100, output="polar")
rads <- seq(0, 20000, 1000)

for(i in 1:nrow(rand)){
	for(j in 1:length(rads)){
		#this random
		one <-  rand[i, , drop=FALSE]
		expect_silent(oneShape <- sc_shape(one, r=rads[j]))

		# test whether the distance of these points is the same!
		allDist <- arcdistmat(one, oneShape)
		expect_true(all(abs(allDist-rads[j])<epsilon))
	}
}
