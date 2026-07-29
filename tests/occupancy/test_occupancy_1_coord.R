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

# manually
occRows <- nrow(unique(nobilis[, c(long="decimalLongitude", lat="decimalLatitude")]))
#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occMat <- occupancy(mat))
expect_equal(occMat, occRows)

## With missing values
### First
matFirstNA <- mat
matFirstNA[1, ] <- NA
expect_silent(occFirstNA <- occupancy(matFirstNA))
occFirstNAmanual <- nrow(unique(matFirstNA[!is.na(matFirstNA[,1]), ]))
expect_equal(occFirstNA, occFirstNAmanual)

### Second
matSecondNA <- mat
matSecondNA[2, ] <- NA
expect_silent(occSecondNA <- occupancy(matSecondNA))
occSecondNAmanual <- nrow(unique(matSecondNA[!is.na(matSecondNA[,1]), ]))
expect_equal(occSecondNA, occSecondNAmanual)

### Last 
matLastNA <- mat
matLastNA[nrow(mat), ] <- NA
expect_silent(occLastNA <- occupancy(matLastNA))
occLastNAmanual <- nrow(unique(matLastNA[!is.na(matLastNA[,1]), ]))
expect_equal(occLastNA, occLastNAmanual)

## 0-length data 
matZero <- mat[0,]
expect_silent(occZero <- occupancy(matZero))
expect_equal(occZero, 0L)

## Singular data
matSingle <- mat[1,, drop=FALSE]
expect_silent(occSingle <- occupancy(matSingle))
expect_equal(occSingle, 1L)

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------
## Duplicates=TRUE
# Should not run
expect_error(occMatDupl <- occupancy(mat, duplicates=TRUE))

 
# one coordinate missing
matOneMiss <- mat
matOneMiss[1,1] <- NA
expect_warning(occOneMiss <- occupancy(matOneMiss))
expect_equal(occOneMiss, occFirstNA)


## Plotting
expect_error(occPlot <- occupancy(mat, plot=TRUE))


#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occMatFull <- occupancy(mat, full=TRUE))
expect_true(inherits(occMatFull, "list"))
expect_equal(names(occMatFull), c("estimate", "occupied"))
expect_equal(occMatFull$estimate, nrow(occMatFull$occupied))
# same as previous
expect_equal(occMatFull$estimate, occMat)


## With missing values
### First
expect_silent(occFirstNAfull <- occupancy(matFirstNA, full=TRUE))
expect_equal(occFirstNA, occFirstNAfull$estimate)
expect_equal(occFirstNAfull$estimate, nrow(occFirstNAfull$occupied))

### Second
expect_silent(occSecondNAfull <- occupancy(matSecondNA, full=TRUE))
expect_equal(occSecondNA, occSecondNAfull$estimate)
expect_equal(occSecondNAfull$estimate, nrow(occSecondNAfull$occupied))

### Last 
expect_silent(occLastNAfull <- occupancy(matLastNA, full=TRUE))
expect_equal(occLastNA, occLastNAfull$estimate)
expect_equal(occLastNAfull$estimate, nrow(occLastNAfull$occupied))

## 0-length data 
expect_silent(occZeroFull <- occupancy(matZero, full=TRUE))
expect_equal(occZeroFull$estimate, 0L)
expect_equal(nrow(occZeroFull$occupied), 0L)

## Singular data
expect_silent(occSingleFull <- occupancy(matSingle, full=TRUE))
expect_equal(occSingleFull$estimate, 1L)
expect_equal(nrow(occSingleFull$occupied), 1L)


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
expect_silent(occ <- occupancy(nobilis, long="decimalLongitude", lat="decimalLatitude"))
expect_equal(occ, occMat)

### First
dfFirstNA <- nobilis
dfFirstNA[1, ] <- NA
expect_silent(occFirstNAdf<- occupancy(dfFirstNA, long="decimalLongitude", lat="decimalLatitude"))
expect_equal(occFirstNA, occFirstNAdf)

### Second
dfSecondNA <- nobilis
dfSecondNA[2, ] <- NA
expect_silent(occSecondNAdf <- occupancy(dfSecondNA, long="decimalLongitude", lat="decimalLatitude"))
expect_equal(occSecondNA, occSecondNAdf)

### Last 
dfLastNA <- nobilis
dfLastNA[nrow(mat), ] <- NA
expect_silent(occLastNAdf <- occupancy(dfLastNA, long="decimalLongitude", lat="decimalLatitude"))
expect_equal(occLastNA, occLastNAdf)

## 0-length data 
dfZero <- nobilis[0,]
expect_silent(occZerodf <- occupancy(dfZero, long="decimalLongitude", lat="decimalLatitude"))
expect_equal(occZerodf, 0L)

## Singular data
dfSingle <- nobilis[1,, drop=FALSE]
expect_silent(occSingledf <- occupancy(dfSingle, long="decimalLongitude", lat="decimalLatitude"))
expect_equal(occSingledf, 1L)

#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------

## With missing values
### First
expect_silent(occFirstNAfullDF <- occupancy(dfFirstNA, long="decimalLongitude", lat="decimalLatitude", full=TRUE))
expect_equal(occFirstNA, occFirstNAfullDF$estimate)
expect_equal(occFirstNAfullDF$estimate, nrow(occFirstNAfullDF$occupied))

### Second
expect_silent(occSecondNAfullDF <- occupancy(dfSecondNA, long="decimalLongitude", lat="decimalLatitude", full=TRUE))
expect_equal(occSecondNA, occSecondNAfullDF$estimate)
expect_equal(occSecondNAfullDF$estimate, nrow(occSecondNAfullDF$occupied))

### Last 
expect_silent(occLastNAfullDF <- occupancy(dfLastNA, long="decimalLongitude", lat="decimalLatitude", full=TRUE))
expect_equal(occLastNA, occLastNAfullDF$estimate)
expect_equal(occLastNAfullDF$estimate, nrow(occLastNAfullDF$occupied))

## 0-length data 
 occupancy(as.matrix(dfZero[, c("decimalLongitude","decimalLatitude")]), full=TRUE)
expect_silent(occZeroFullDF <- occupancy(dfZero, long="decimalLongitude", lat="decimalLatitude", full=TRUE))
expect_equal(occZeroFullDF$estimate, 0L)
expect_equal(nrow(occZeroFullDF$occupied), 0L)

## Singular data
expect_silent(occSingleFullDF <- occupancy(dfSingle, long="decimalLongitude", lat="decimalLatitude", full=TRUE))
expect_equal(occSingleFullDF$estimate, 1L)
expect_equal(nrow(occSingleFullDF$occupied), 1L)

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------
expect_error(occupancy(nobilis))
expect_error(occupancy(nobilis, long="gibberish", lat="waste"))

## Plotting
expect_error(occPlot <- occupancy(nobilis, long="decimalLongitude", lat="decimalLatitude", plot=TRUE))

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
