library(tinytest)
library(biodome)

# create a grid
hex <- hexagrid(2, sf=TRUE)

# an example shape
shape <- paste0("F", c(4, 5, 11, 13, 15, 21, 24, 26, 32, 33, 34, 35, 36))


# character-trigrid method
# the gappiness
expect_silent(one <- ranges_gappiness(shape, hex))
expect_equal(class(one), "list")
expect_equal(length(one), 2L)
expect_equal(one$estimate, 5/18)

# no holes
shape <- paste0("F", c(4, 5, 11))
expect_silent(one <- ranges_gappiness(shape, hex))
expect_equal(0, one$estimate)
expect_null(one$holes)


# full control with the seed
set.seed(5)
expect_silent(kentOne <- kentPresence(10, centers=cbind(0,0), kappa=5))

# look up the cells
cells <- locate(hex, kentOne)

# in case needs to be seen
## plot(hex)
## points(kentOne)
## plot(hex, cells, col="#FF550033", add=TRUE)
expect_silent(two <- ranges_gappiness(kentOne, hex))
expect_equal(class(two), "list")
expect_equal(length(two), 2L)
expect_equal(two$estimate, 1/9)
hole <- 1L
names(hole) <- "F23"
expect_equal(two$holes, hole)
