library(tinytest)

# create a grid
hex <- hexagrid(2, sf=TRUE)

################################################################################
# Known output with basic and full

# an example shape
shape <- paste0("F", c(4, 5, 11, 13, 15, 21, 24, 26, 32, 33, 34, 35, 36))

# character-trigrid method
# the gappiness
expect_silent(one <- gappiness(shape, hex))

expect_true(inherits(one, "numeric"))
expect_true(length(one)==1)
expect_equal(one, 5/18)

# wih full output
expect_silent(oneFull <- gappiness(shape, hex, full=TRUE))
expect_equal(class(oneFull), "list")
expect_equal(oneFull$estimate, one)
expect_equal(length(oneFull), 2L)

# character-trigrid method
# the gappiness
expect_silent(oneEx <- gappiness(shape, hex, exclude="F4"))
expect_equal(oneEx, 4/17)

################################################################################
# Known output with no holes basic and full
# no holes
shape <- paste0("F", c(4, 5, 11))
expect_silent(one <- gappiness(shape, hex))
expect_equal(0, one)
# full output
expect_silent(oneFull <- gappiness(shape, hex, full=TRUE))
expect_equal(class(oneFull), "list")
expect_null(oneFull$holes)

#  single cell
shape <- paste0("F", c(4))
expect_silent(one <- gappiness(shape, hex))
expect_equal(0, one)
# full output
expect_silent(oneFull <- gappiness(shape, hex, full=TRUE))
expect_equal(class(oneFull), "list")
expect_null(oneFull$holes)


# character-trigrid method
# the gappiness
expect_silent(oneEx <- gappiness(shape, hex, exclude="F4"))
expect_equal(oneEx, NA)
expect_silent(oneExFull <- gappiness(shape, hex, exclude="F4", full=TRUE))
expect_true(inherits(oneExFull, "list") )











## # full control with the seed
## set.seed(5)
## expect_silent(kentOne <- kentPresence(10, centers=cbind(0,0), kappa=5))

## # look up the cells
## cells <- locate(hex, kentOne)

## # in case needs to be seen
## ## plot(hex)
## ## points(kentOne)
## ## plot(hex, cells, col="#FF550033", add=TRUE)
## expect_silent(two <- ranges_gappiness(kentOne, hex))
## expect_equal(class(two), "list")
## expect_equal(length(two), 2L)
## expect_equal(two$estimate, 1/9)
## hole <- 1L
## names(hole) <- "F23"
## expect_equal(two$holes, hole)
