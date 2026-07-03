# Occupancy Tests for for coordinate inputs (x=coordMat coordDF)
# and pre-specified localities in the dataset(s=character)
library(tinytest)
library(orange)

# setup
data(pinna)
nobilis <- pinna[pinna$species=="Pinna nobilis", ]

# use an icosahedral grid for localities (precursor to icosa-specific methods)
suppressMessages(hex<- hexagrid(spacing=2))
nobilis$cell <- locate(hex, nobilis[, c("decimalLongitude", "decimalLatitude")])



## # diagnostics
## hex<- hexagrid(spacing=2, sf=TRUE) 
## plot(hex, xlim=c(-15, 40), ylim=c(30, 50), border="gray80")
## plot(hex, names(table(nobilis$cell)), col="#BB000066", add=TRUE, border="white")
## points(nobilis[, c("decimalLongitude", "decimalLatitude")])

# missing locality!

################################################################################
# 1. Single data.frame methods
################################################################################

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occup <- occupancy(nobilis, s="cell"))
expect_equal(class(occup), "numeric")
occupManual <- length(levels(factor(nobilis$cell)))
expect_equal(occupManual, occup)

## With missing values
### First
dfFirstNA <- nobilis
dfFirstNA[1, ] <- NA
expect_silent(occupFirstNA <- occupancy(dfFirstNA, s="cell"))
expect_equal(length(levels(factor(dfFirstNA$cell))), occupFirstNA)

### Second
dfSecondNA <- nobilis
dfSecondNA[2, ] <- NA
expect_silent(occupSecondNA <- occupancy(dfSecondNA, s="cell"))
expect_equal(length(levels(factor(dfSecondNA$cell))), occupSecondNA)

### Last
dfLastNA <- nobilis
dfLastNA[nrow(nobilis), ] <- NA
expect_silent(occupLastNA <- occupancy(dfLastNA, s="cell"))
expect_equal(length(levels(factor(dfLastNA$cell))), occupLastNA)

## 0-length data 
dfZero <- nobilis[0,]
expect_silent(occupZero <- occupancy(dfZero, s="cell"))
expect_equal(occupZero, 0L)

## Singular data
dfSingle <- nobilis[1,]
expect_silent(occupSingle <- occupancy(dfSingle, s="cell"))
expect_equal(occupSingle, 1L)

dfMissS <- nobilis
dfMissS$cell[1] <- NA
expect_silent(occupMissS <- occupancy(dfMissS, s="cell"))
expect_equal(length(levels(factor(dfMissS$cell))), occupMissS)

 
#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occupFull <- occupancy(nobilis, s="cell", full=TRUE))
expect_true(inherits(occupFull, "list"))
expect_equal(names(occupFull), c("estimate", "occupied"))
expect_equal(occup, occupFull$estimate)
expect_equal(length(occupFull$occupied), occupFull$estimate)

## With missing values
### First
expect_silent(occupFirstNAfull <- occupancy(dfFirstNA, s="cell", full=TRUE))
expect_equal(occupFirstNAfull$estimate, occupFirstNA)
expect_equal(length(occupFirstNAfull$occupied), occupFirstNAfull$estimate)

### Second
expect_silent(occupSecondNAfull <- occupancy(dfSecondNA, s="cell", full=TRUE))
expect_equal(occupSecondNAfull$estimate, occupSecondNA)
expect_equal(length(occupSecondNAfull$occupied), occupSecondNAfull$estimate)

### Last
expect_silent(occupLastNAfull <- occupancy(dfLastNA, s="cell", full=TRUE))
expect_equal(occupLastNAfull$estimate, occupLastNA)
expect_equal(length(occupLastNAfull$occupied), occupLastNAfull$estimate)

## 0-length data 
expect_silent(occupZeroFull <- occupancy(dfZero, s="cell", full=TRUE))
expect_equal(occupZeroFull$estimate, 0L)
expect_equal(length(occupZeroFull$occupied), occupZeroFull$estimate)

## Singular data
expect_silent(occupSingleFull <- occupancy(dfSingle, s="cell", full=TRUE))
expect_equal(occupSingleFull$estimate, 1L)
expect_equal(length(occupSingleFull$occupied), occupSingleFull$estimate)


## Plotting
expect_error(occPlot <- occupancy(pinna, s="cell",  plot=TRUE))

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------

# not required longitude and latitude
expect_error(occupancy(nobilis, s="cell", long="decimalLongitude", lat="decimalLatitude"))
expect_error(occupancy(nobilis, s="WRONG"))

################################################################################
# 2. sfc - not yet! 
################################################################################

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------
