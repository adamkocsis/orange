library(tinytest)

# Testing the uniform distribution generator

p <- 0.3
s <- 30
k <- 150


# 1. Single probability, integer species and cells
expect_silent(one <- uniformOcc(p=p, s=s, k=k))
expect_equal(dim(one), c(k, s))
expect_equal(rownames(one), paste0("C", 1:k))
expect_equal(colnames(one), paste0("T", 1:s))
sums <- apply(one, 2, sum)

# wrong input
expect_error(one <- uniformOcc(p=NA, s=s, k=k))
expect_error(one <- uniformOcc(p=12, s=s, k=k))
expect_error(one <- uniformOcc(p=-3, s=s, k=k))
expect_error(one <- uniformOcc(p=0.2, s=-3, k=k))
expect_error(one <- uniformOcc(p=0.2, s=c(3,4), k=k))

# 2.  multiple probability values
p <- seq(0.1, 0.9, length.out=6)

# invalid input
expect_error(uniformOcc(p=p, s=s, k=k))
expect_silent(two <- uniformOcc(p=p,k=k))
expect_equal(dim(two), c(k, length(p)))
expect_equal(rownames(two), paste0("C", 1:k))
expect_equal(colnames(two), paste0("T", 1:length(p)))
sums <- apply(two, 2, sum)
expect_true(all(diff(sums)>0))


# 3. Named probabilitiy values
names(p) <- paste0("Tax", 1:length(p))
expect_silent(twoNamed <- uniformOcc(p=p, k=k))
expect_equal(colnames(twoNamed), names(p))



# I. Using an icosa grid
hex <- icosa::hexagrid(deg=10)

# II/1. Single probabilitiy
p <- 0.3
expect_silent(three <- uniformOcc(g=hex, p=p, s=s))
expect_equal(dim(three), c(as.numeric(length(hex)), s))
expect_equal(rownames(three), faces(hex))
expect_equal(colnames(three), paste0("T", 1:s))

# II/2. Muliple probabilitiy
p <- seq(0.1, 0.9, length.out=6)
expect_silent(four <- uniformOcc(g=hex, p=p))
expect_equal(dim(four), c(as.numeric(length(hex)), length(p)))
expect_equal(rownames(four), faces(hex))
expect_equal(colnames(four), paste0("T", 1:length(p)))
