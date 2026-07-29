# Proportional occupancy testing
# For coordinate inputs (x=coordMat coordDF) -
library(tinytest)
library(orange)

# setup
data(pinna)
nobilis <- pinna[pinna$species=="Pinna nobilis", ]

################################################################################
# 1. Matrix methods
################################################################################
mat <- as.matrix(nobilis[,c("decimalLongitude", "decimalLatitude")])

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

## Complete dataset
expect_error(occMat <- occupancy(mat, prop="garbage"))
expect_error(occMat <- occupancy(mat, prop="global"))
expect_silent(occMat <- occupancy(mat, prop="relative"))
expect_equal(occMat, 1L)

## With missing values
### First
matFirstNA <- mat
matFirstNA[1, ] <- NA
expect_error(occupancy(matFirstNA, prop="global"))
expect_silent(occFirstNA <- occupancy(matFirstNA, prop="relative"))
expect_equal(occFirstNA, 1L)

### Second
matSecondNA <- mat
matSecondNA[2, ] <- NA
expect_error(occupancy(matSecondNA, prop="global"))
expect_silent(occSecondNA <- occupancy(matSecondNA, prop="relative"))
expect_equal(occSecondNA, 1L)

### Last 
matLastNA <- mat
matLastNA[nrow(mat), ] <- NA
expect_error(occupancy(occLastNA, prop="global"))
expect_silent(occLastNA <- occupancy(matLastNA, prop="relative"))
expect_equal(occLastNA, 1L)

## 0-length data 
matZero <- mat[0,]
expect_error(occupancy(matZero, prop="global"))
expect_silent(occZero <- occupancy(matZero, prop="relative"))
expect_equal(occZero, NA)

## Singular data
matSingle <- mat[1,, drop=FALSE]
expect_error(occupancy(matSingle, prop="global"))
expect_silent(occSingle <- occupancy(matSingle, prop="relative"))
expect_equal(occSingle, 1L)

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------
## Duplicates=TRUE
# Should not run
expect_error(occupancy(mat, duplicates=TRUE, prop="global"))
expect_error(occMatDupl <- occupancy(mat, duplicates=TRUE, prop="relative"))

 
# one coordinate missing
matOneMiss <- mat
matOneMiss[1,1] <- NA
expect_warning(occupancy(matOneMiss, prop="global"))
expect_warning(occOneMiss <- occupancy(matOneMiss, prop="relative"))
expect_equal(occOneMiss, occFirstNA)


## Plotting
expect_error(occPlot <- occupancy(mat, plot=TRUE, prop="relative"))


#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------

## Complete dataset
expect_error(occupancy(mat, full=TRUE, prop="global"))
expect_silent(occMatFull <- occupancy(mat, full=TRUE, prop="relative"))
expect_true(inherits(occMatFull, "list"))
expect_equal(names(occMatFull), c("estimate", "occupied"))
expect_equal(occMatFull$estimate, 1L)
# same as previous
expect_equal(occMatFull$estimate, occMat)


## With missing values
### First
expect_error(occupancy(matFirstNA, full=TRUE, prop="global"))
expect_silent(occFirstNAfull <- occupancy(matFirstNA, full=TRUE, prop="relative"))
expect_equal(occFirstNA, occFirstNAfull$estimate)

### Second
expect_error(occupancy(matSecondNA, full=TRUE, prop="global"))
expect_silent(occSecondNAfull <- occupancy(matSecondNA, full=TRUE, prop="relative"))
expect_equal(occSecondNA, occSecondNAfull$estimate)

### Last 
expect_error(occupancy(matLastNA, full=TRUE, prop="global"))
expect_silent(occLastNAfull <- occupancy(matLastNA, full=TRUE, prop="relative"))
expect_equal(occLastNA, occLastNAfull$estimate)

## 0-length data 
expect_error(occupancy(matZero, full=TRUE, prop="global"))
expect_silent(occZeroFull <- occupancy(matZero, full=TRUE, prop="relative"))
expect_equal(occZeroFull$estimate, NA)

## Singular data
expect_error(occupancy(matSingle, full=TRUE, prop="global"))
expect_silent(occSingleFull <- occupancy(matSingle, full=TRUE, prop="relative"))
expect_equal(occSingleFull$estimate, 1L)


#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------


################################################################################
# 2. Single data.frame methods
################################################################################

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

## Complete dataset
expect_error(occupancy(nobilis, long="decimalLongitude", lat="decimalLatitude", prop="global"))
expect_silent(occ <- occupancy(nobilis, long="decimalLongitude", lat="decimalLatitude", prop="relative"))
expect_equal(occ, occMat)

### First
dfFirstNA <- nobilis
dfFirstNA[1, ] <- NA
expect_error( occupancy(dfFirstNA, long="decimalLongitude", lat="decimalLatitude", prop="global"))
expect_silent(occFirstNAdf<- occupancy(dfFirstNA, long="decimalLongitude", lat="decimalLatitude", prop="relative"))
expect_equal(occFirstNA, occFirstNAdf)

### Second
dfSecondNA <- nobilis
dfSecondNA[2, ] <- NA
expect_error( occupancy(dfSecondNA, long="decimalLongitude", lat="decimalLatitude", prop="global"))
expect_silent(occSecondNAdf <- occupancy(dfSecondNA, long="decimalLongitude", lat="decimalLatitude", prop="relative"))
expect_equal(occSecondNA, occSecondNAdf)

### Last 
dfLastNA <- nobilis
dfLastNA[nrow(mat), ] <- NA
expect_error( occupancy(dfLastNA, long="decimalLongitude", lat="decimalLatitude", prop="global"))
expect_silent(occLastNAdf <- occupancy(dfLastNA, long="decimalLongitude", lat="decimalLatitude", prop="relative"))
expect_equal(occLastNA, occLastNAdf)

## 0-length data 
dfZero <- nobilis[0,]
expect_error( occupancy(dfZero, long="decimalLongitude", lat="decimalLatitude", prop="global"))
expect_silent(occZerodf <- occupancy(dfZero, long="decimalLongitude", lat="decimalLatitude", prop="relative"))
expect_equal(occZerodf, NA)

## Singular data
dfSingle <- nobilis[1,, drop=FALSE]
expect_error( occupancy(dfSingle, long="decimalLongitude", lat="decimalLatitude", prop="global"))
expect_silent(occSingledf <- occupancy(dfSingle, long="decimalLongitude", lat="decimalLatitude", prop="relative"))
expect_equal(occSingledf, 1L)

#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------

## With missing values
### First
expect_error( occupancy(dfFirstNA, long="decimalLongitude", lat="decimalLatitude", full=TRUE, prop="global"))
expect_silent(occFirstNAfullDF <- occupancy(dfFirstNA, long="decimalLongitude", lat="decimalLatitude", full=TRUE, prop="relative"))
expect_equal(occFirstNA, occFirstNAfullDF$estimate)

### Second
expect_error( occupancy(dfSecondNA, long="decimalLongitude", lat="decimalLatitude", full=TRUE, prop="global"))
expect_silent(occSecondNAfullDF <- occupancy(dfSecondNA, long="decimalLongitude", lat="decimalLatitude", full=TRUE, prop="relative"))
expect_equal(occSecondNA, occSecondNAfullDF$estimate)

### Last 
expect_error( occupancy(dfLastNA, long="decimalLongitude", lat="decimalLatitude", full=TRUE, prop="global"))
expect_silent(occLastNAfullDF <- occupancy(dfLastNA, long="decimalLongitude", lat="decimalLatitude", full=TRUE, prop="relative"))
expect_equal(occLastNA, occLastNAfullDF$estimate)

## 0-length data 
expect_error( occupancy(dfZero, long="decimalLongitude", lat="decimalLatitude", full=TRUE, prop="global"))
expect_silent(occZeroFullDF <- occupancy(dfZero, long="decimalLongitude", lat="decimalLatitude", full=TRUE, prop="relative"))
expect_equal(occZeroFullDF$estimate, NA)

## Singular data
expect_error( occupancy(dfSingle, long="decimalLongitude", lat="decimalLatitude", full=TRUE, prop="global"))
expect_silent(occSingleFullDF <- occupancy(dfSingle, long="decimalLongitude", lat="decimalLatitude", full=TRUE, prop="relative"))
expect_equal(occSingleFullDF$estimate, 1L)

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------
expect_error(occupancy(nobilis, prop="relative"))
expect_error(occupancy(nobilis, long="gibberish", lat="waste", prop="relative"))

## Plotting
expect_error(occPlot <- occupancy(nobilis, long="decimalLongitude", lat="decimalLatitude", plot=TRUE, prop="relative"))

################################################################################
# 3. sfc - not yet! 
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
