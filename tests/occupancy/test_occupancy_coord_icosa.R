# For coordinate inputs (x=coordMat coordDF) + s="icosa", not iteration

# Setup
library(tinytest)
library(orange)

# setup
data(pinna)
nobilis <- pinna[pinna$species=="Pinna nobilis", ]
mat <- as.matrix(nobilis[, c("decimalLongitude", "decimalLatitude")])

# mke a hexagrid
suppressMessages(hex <- hexagrid(deg=2, sf=TRUE))

################################################################################
# 1. Matrix methods
################################################################################

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occ <- occupancy(x=mat, s=hex))
manual <- function(x, s){
	# get rid of missing vals
	cellsManual <- levels(factor((locate(s,x ))))
	length(cellsManual)
}
expect_equal(manual(x=mat, s=hex), occ)

## With missing values
### First
matFirstNA <- mat
matFirstNA[1, ] <- NA
expect_silent(occFirstNA <- occupancy(matFirstNA, s=hex))
expect_equal(occFirstNA, manual(matFirstNA, s=hex))

### Second
matSecondNA <- mat
matSecondNA[2, ] <- NA
expect_silent(occSecondNA <- occupancy(matSecondNA, s=hex))
expect_equal(occSecondNA, manual(matSecondNA, s=hex))

### Last 
matLastNA <- mat
matLastNA[nrow(mat), ] <- NA
expect_silent(occLastNA <- occupancy(matLastNA, s=hex))
expect_equal(occLastNA, manual(matLastNA, s=hex))

## 0-length data 
matZero <- mat[0,]
expect_silent(occZero <- occupancy(matZero, s=hex))
expect_equal(occZero, 0L)

## Singular data
matSingle <- mat[1,,drop=FALSE]
expect_silent(occSingle <- occupancy(matSingle, s=hex))
expect_equal(occSingle, 1L)

## Duplicates=TRUE
expect_error(occupancy(matSingle, s=hex, duplicates=TRUE))
expect_error(occupancy(matSingle, s=hex, duplicates=FALSE))
 
## Plotting
# error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE))



#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occFull <- occupancy(x=mat, s=hex, full=TRUE))
expect_true(inherits(occFull, "orange"))
expect_equal(occ, occFull$estimate)
expect_equal(occFull$estimate, length(occFull$occupied))

### First
expect_silent(occFirstNAfull <- occupancy(matFirstNA, s=hex, full=TRUE))
expect_equal(occFirstNA, occFirstNAfull$estimate)
expect_equal(length(occFirstNAfull$occupied), occFirstNAfull$estimate)

### Second
expect_silent(occSecondNAfull <- occupancy(matSecondNA, s=hex, full=TRUE))
expect_equal(occSecondNA, occSecondNAfull$estimate)
expect_equal(length(occSecondNAfull$occupied), occSecondNAfull$estimate)

### Last 
expect_silent(occLastNAfull <- occupancy(matLastNA, s=hex, full=TRUE))
expect_equal(occLastNA, occLastNAfull$estimate)
expect_equal(length(occLastNAfull$occupied), occLastNAfull$estimate)

## 0-length data 
expect_silent(occZeroFull <- occupancy(matZero, s=hex, full=TRUE))
expect_equal(occZeroFull$estimate, 0L)
expect_equal(length(occZeroFull$occupied), occZeroFull$estimate)

## Singular data
expect_silent(occSingleFull <- occupancy(matSingle, s=hex, full=TRUE))
expect_equal(occSingleFull$estimate, 1L)
expect_equal(length(occSingleFull$occupied), occSingleFull$estimate)

## Duplicates=TRUE
expect_error(occupancy(matSingle, s=hex, duplicates=TRUE, full=TRUE))
expect_error(occupancy(matSingle, s=hex, duplicates=FALSE, full=TRUE))

## Plotting
# error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE, full=TRUE))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE, full=TRUE))

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------

# wrong column names
expect_silent(occupancy(x=mat, s=hex, long="decimalLongitude",lat="decimalLatitude"))
expect_error(occupancy(x=mat, s=hex, long="sdf"))
expect_error(occupancy(x=mat, s=hex, lat="sdf"))
expect_error(occupancy(x=mat, s=hex, long="sdf",lat="sdf"))

# not implemented yet
expect_error(occupancy(mat, hex, q=0.7))


################################################################################
# 2. Single data.frame methods
################################################################################

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occDF <- occupancy(x=nobilis, long="decimalLongitude",lat="decimalLatitude", s=hex))
manual <- function(x, s){
	# get rid of missing vals
	cellsManual <- levels(factor((locate(s,x[, c("decimalLongitude", "decimalLatitude")] ))))
	length(cellsManual)
}
expect_equal(manual(x=nobilis, s=hex), occDF)
expect_equal(occ, occDF)

## With missing values
### First
dfFirstNA <- nobilis
dfFirstNA[1, ] <- NA
expect_silent(occFirstNAdf <- occupancy(dfFirstNA, long="decimalLongitude",lat="decimalLatitude",s=hex))
expect_equal(occFirstNAdf, manual(dfFirstNA, s=hex))
expect_silent(occFirstNA,occFirstNAdf)

### Second
dfSecondNA <- nobilis
dfSecondNA[2, ] <- NA
expect_silent(occSecondNAdf <- occupancy(dfSecondNA, long="decimalLongitude",lat="decimalLatitude",s=hex))
expect_equal(occSecondNAdf, manual(dfSecondNA, s=hex))
expect_silent(occSecondNA,occSecondNAdf)

### Last 
dfLastNA <- nobilis
dfLastNA[nrow(nobilis), ] <- NA
expect_silent(occLastNAdf <- occupancy(dfLastNA, long="decimalLongitude",lat="decimalLatitude",s=hex))
expect_equal(occLastNAdf, manual(dfLastNA, s=hex))
expect_silent(occLastNA,occLastNAdf)

## 0-length data 
dfZero <- nobilis[0,]
expect_silent(occZerodf <- occupancy(dfZero, long="decimalLongitude",lat="decimalLatitude",s=hex))
expect_equal(occZerodf, 0L)

## Singular data
dfSingle <- nobilis[1,,drop=FALSE]
expect_silent(occSingledf <- occupancy(dfSingle, long="decimalLongitude",lat="decimalLatitude",s=hex))
expect_equal(occSingledf, 1L)

## Duplicates=TRUE
expect_error(occupancy(dfSingle, long="decimalLongitude",lat="decimalLatitude", s=hex, duplicates=TRUE))
expect_error(occupancy(dfSingle, long="decimalLongitude",lat="decimalLatitude", s=hex, duplicates=FALSE))
 
## Plotting
# error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(nobilis, long="decimalLongitude",lat="decimalLatitude", hex, plot=TRUE))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(nobilis, long="decimalLongitude",lat="decimalLatitude", hex, plot=TRUE))

# reasonable column defaults for the data.frame-method
df <- as.data.frame(mat)
colnames(df) <- c("long", "lat")
expect_silent(occDFdefault <- occupancy(df, s=hex))
expect_equal(occDFdefault, occ)

#-------------------------------------------------------------------------------
# Full output/traceability 
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occDFFull <- occupancy(x=nobilis, long="decimalLongitude",lat="decimalLatitude", s=hex, full=TRUE))
expect_true(inherits(occDFFull, "orange"))
expect_equal(occDF, occDFFull$estimate)
expect_equal(occDFFull$estimate, length(occDFFull$occupied))

## With missing values
### First
expect_silent(occFirstNAdfFull <- occupancy(dfFirstNA, long="decimalLongitude",lat="decimalLatitude",s=hex, full=TRUE))
expect_equal(occFirstNAdf, occFirstNAdfFull$estimate)
expect_equal(length(occFirstNAdfFull$occupied), occFirstNAdfFull$estimate)

### Second
expect_silent(occSecondNAdfFull <- occupancy(dfSecondNA, long="decimalLongitude",lat="decimalLatitude",s=hex, full=TRUE))
expect_equal(occSecondNAdf, occSecondNAdfFull$estimate)
expect_equal(length(occSecondNAdfFull$occupied), occSecondNAdfFull$estimate)

### Last 
expect_silent(occLastNAdfFull <- occupancy(dfLastNA, long="decimalLongitude",lat="decimalLatitude",s=hex, full=TRUE))
expect_equal(occLastNAdf, occLastNAdfFull$estimate)
expect_equal(length(occLastNAdfFull$occupied), occLastNAdfFull$estimate)

## 0-length data 
expect_silent(occZerodfFull <- occupancy(dfZero, long="decimalLongitude",lat="decimalLatitude",s=hex, full=TRUE))
expect_equal(occZerodfFull$estimate, 0L)
expect_equal(length(occZerodfFull$occupied), occZerodfFull$estimate)

## Singular data
expect_silent(occSingledfFull <- occupancy(dfSingle, long="decimalLongitude",lat="decimalLatitude",s=hex, full=TRUE))
expect_equal(occSingleFull$estimate, 1L)
expect_equal(length(occSingleFull$occupied), occSingleFull$estimate)

 
# Plotting
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(nobilis, long="decimalLongitude",lat="decimalLatitude", hex, plot=TRUE, full=TRUE))

# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(nobilis, long="decimalLongitude",lat="decimalLatitude", hex, plot=TRUE, full=TRUE))

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------

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
