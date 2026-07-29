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


################################################################################
# 1. Single data.frame methods
################################################################################

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

## Complete dataset
expect_error(occupancy(nobilis, s="cell", prop="global"))
expect_silent(occup <- occupancy(nobilis, s="cell", prop="relative"))
expect_equal(class(occup), "numeric")
occupManual <- 1
expect_equal(occupManual, occup)

## With missing values
### First
dfFirstNA <- nobilis
dfFirstNA[1, ] <- NA
expect_error(occupFirstNA <- occupancy(dfFirstNA, s="cell", prop="global"))
expect_silent(occupFirstNA <- occupancy(dfFirstNA, s="cell", prop="relative"))
expect_equal(1, occupFirstNA)

### Second
dfSecondNA <- nobilis
dfSecondNA[2, ] <- NA
expect_error(occupSecondNA <- occupancy(dfSecondNA, s="cell", prop="global"))
expect_silent(occupSecondNA <- occupancy(dfSecondNA, s="cell", prop="relative"))
expect_equal(1, occupSecondNA)

### Last
dfLastNA <- nobilis
dfLastNA[nrow(nobilis), ] <- NA
expect_error(occupLastNA <- occupancy(dfLastNA, s="cell", prop="global"))
expect_silent(occupLastNA <- occupancy(dfLastNA, s="cell", prop="relative"))
expect_equal(1, occupLastNA)

## 0-length data 
dfZero <- nobilis[0,]
expect_error(occupZero <- occupancy(dfZero, s="cell", prop="global"))
expect_silent(occupZero <- occupancy(dfZero, s="cell", prop="relative"))
expect_equal(occupZero, NA)

## Singular data
dfSingle <- nobilis[1,]
expect_error(occupSingle <- occupancy(dfSingle, s="cell",prop="global"))
expect_silent(occupSingle <- occupancy(dfSingle, s="cell",prop="relative"))
expect_equal(occupSingle, 1L)

dfMissS <- nobilis
dfMissS$cell[1] <- NA
expect_error(occupMissS <- occupancy(dfMissS, s="cell", prop="global"))
expect_silent(occupMissS <- occupancy(dfMissS, s="cell", prop="relative"))
expect_equal(1, occupMissS)
 
#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------

occManual <- function(x){
	levels(factor(x$cell))
}

## Complete dataset
expect_error(occupFull <- occupancy(nobilis, s="cell", full=TRUE, prop="global"))
expect_silent(occupFull <- occupancy(nobilis, s="cell", full=TRUE, prop="relative"))
expect_true(inherits(occupFull, "orange"))
expect_equal(names(occupFull), c("estimate", "occupied"))
expect_equal(occup, occupFull$estimate)
expect_equal(occManual(nobilis), occupFull$occupied)

## With missing values
### First
expect_error(occupFirstNAfull <- occupancy(dfFirstNA, s="cell", full=TRUE, prop="global"))
expect_silent(occupFirstNAfull <- occupancy(dfFirstNA, s="cell", full=TRUE, prop="relative"))
expect_true(inherits(occupFirstNAfull, "orange"))
expect_equal(occupFirstNAfull$estimate, occupFirstNA)
expect_equal(occManual(dfFirstNA), occupFirstNAfull$occupied)

### Second
expect_error(occupSecondNAfull <- occupancy(dfSecondNA, s="cell", full=TRUE, prop="global"))
expect_silent(occupSecondNAfull <- occupancy(dfSecondNA, s="cell", full=TRUE, prop="relative"))
expect_true(inherits(occupSecondNAfull, "orange"))
expect_equal(occupSecondNAfull$estimate, occupSecondNA)
expect_equal(occManual(dfSecondNA), occupSecondNAfull$occupied)

### Last
expect_error(occupLastNAfull <- occupancy(dfLastNA, s="cell", full=TRUE, prop="global"))
expect_silent(occupLastNAfull <- occupancy(dfLastNA, s="cell", full=TRUE, prop="relative"))
expect_true(inherits(occupLastNAfull, "orange"))
expect_equal(occupLastNAfull$estimate, occupLastNA)
expect_equal(occManual(dfLastNA), occupLastNAfull$occupied)

## 0-length data 
expect_error(occupancy(dfZero, s="cell", full=TRUE, prop="global"))
expect_silent(occupZeroFull <- occupancy(dfZero, s="cell", full=TRUE, prop="relative"))
expect_true(inherits(occupZeroFull, "orange"))
expect_equal(occupZeroFull$estimate, NA)
expect_equal(length(occupZeroFull$occupied), 0)

## Singular data
expect_error(occupSingleFull <- occupancy(dfSingle, s="cell", full=TRUE, prop="global"))
expect_silent(occupSingleFull <- occupancy(dfSingle, s="cell", full=TRUE, prop="relative"))
expect_equal(occupSingleFull$estimate, 1L)
expect_true(inherits(occupSingleFull, "orange"))
expect_equal(occManual(dfSingle), occupSingleFull$occupied)


## Plotting
expect_error(occPlot <- occupancy(pinna, s="cell",  plot=TRUE, prop="global"))
expect_error(occPlot <- occupancy(pinna, s="cell",  plot=TRUE, prop="relative"))

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------

# not required longitude and latitude
expect_error(occupancy(nobilis, s="cell", long="decimalLongitude", lat="decimalLatitude", prop="global"))
expect_error(occupancy(nobilis, s="cell", long="decimalLongitude", lat="decimalLatitude", prop="relative"))
expect_error(occupancy(nobilis, s="WRONG", prop="relative"))
expect_error(occupancy(nobilis, s="WRONG", prop="global"))

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
