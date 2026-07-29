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

# manual checking
manualGlobal <- function(x, s){
	# get rid of missing vals
	cellsManual <- levels(factor((locate(s,x ))))
	length(cellsManual)/as.numeric(length(s))
}

manualRelative <- function(x, s){
	# get rid of missing vals
	cellsManual <- levels(factor((locate(s,x ))))
	length(cellsManual)/length(cellsManual)
}

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occCount <- occupancy(x=mat, s=hex))
expect_silent(occGlobal <- occupancy(x=mat, s=hex, prop="global"))
expect_equal(manualGlobal(x=mat, s=hex), occGlobal)
expect_silent(occRelative <- occupancy(x=mat, s=hex, prop="relative"))
expect_equal(occRelative, 1)

## With missing values
### First
matFirstNA <- mat
matFirstNA[1, ] <- NA
expect_silent(occFirstNAglobal <- occupancy(matFirstNA, s=hex, prop="global"))
expect_equal(occFirstNAglobal, manualGlobal(matFirstNA, s=hex))
expect_silent(occFirstNArelative <- occupancy(matFirstNA, s=hex, prop="relative"))
expect_equal(occFirstNArelative, 1)

### Second
matSecondNA <- mat
matSecondNA[2, ] <- NA
expect_silent(occSecondNAglobal <- occupancy(matSecondNA, s=hex, prop="global"))
expect_equal(occSecondNAglobal, manualGlobal(matSecondNA, s=hex))
expect_silent(occSecondNArelative <- occupancy(matSecondNA, s=hex, prop="relative"))
expect_equal(occSecondNArelative, 1)

### Last 
matLastNA <- mat
matLastNA[nrow(mat), ] <- NA
expect_silent(occLastNAglobal <- occupancy(matLastNA, s=hex, prop="global"))
expect_equal(occLastNAglobal, manualGlobal(matLastNA, s=hex))
expect_silent(occLastNArelative <- occupancy(matLastNA, s=hex, prop="relative"))
expect_equal(occLastNArelative, 1)

## 0-length data 
matZero <- mat[0,]
expect_silent(occZeroGlobal <- occupancy(matZero, s=hex, prop="global"))
expect_equal(occZeroGlobal, 0)
expect_silent(occZeroRelative <- occupancy(matZero, s=hex, prop="relative"))
expect_equal(occZeroRelative, NA)

## Singular data
matSingle <- mat[1,,drop=FALSE]
expect_silent(occSingleGlobal <- occupancy(matSingle, s=hex, prop="global"))
expect_equal(occSingleGlobal, as.numeric((1/length(hex))))
expect_silent(occSingleRelative <- occupancy(matSingle, s=hex, prop="relative"))
expect_equal(occSingleRelative, 1)

## Duplicates=TRUE
expect_error(occupancy(matSingle, s=hex, duplicates=TRUE, prop="global"))
expect_error(occupancy(matSingle, s=hex, duplicates=TRUE, prop="relative"))
expect_error(occupancy(matSingle, s=hex, duplicates=FALSE, prop="global"))
expect_error(occupancy(matSingle, s=hex, duplicates=FALSE, prop="relative"))
 
## Plotting
# error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE, prop="global"))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE, prop="global"))

# error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE, prop="relative"))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE, prop="relative"))



#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occFullGlobal <- occupancy(x=mat, s=hex, full=TRUE, prop="global"))
expect_true(inherits(occFullGlobal, "orange"))
expect_equal(occGlobal, occFullGlobal$estimate)
expect_equal(occFullGlobal$estimate, length(occFullGlobal$occupied)/as.numeric(length(hex)))

expect_silent(occFullRelative <- occupancy(x=mat, s=hex, full=TRUE, prop="relative"))
expect_true(inherits(occFullRelative, "orange"))
expect_equal(occRelative, occFullRelative$estimate)
expect_equal(occFullRelative$estimate, length(occFullRelative$occupied)/as.numeric(length(occFullRelative$occupied)))

### First
expect_silent(occFirstNAfullGlobal <- occupancy(matFirstNA, s=hex, full=TRUE, prop="global"))
expect_equal(occFirstNAglobal, occFirstNAfullGlobal$estimate)
expect_equal(length(occFirstNAfullGlobal$occupied)/as.numeric(length(hex)), occFirstNAfullGlobal$estimate)

expect_silent(occFirstNAfullRelative <- occupancy(matFirstNA, s=hex, full=TRUE, prop="relative"))
expect_equal(occFirstNArelative, occFirstNAfullRelative$estimate)
expect_equal(1, occFirstNAfullRelative$estimate)

### Second
expect_silent(occSecondNAfullGlobal <- occupancy(matSecondNA, s=hex, full=TRUE, prop="global"))
expect_equal(occSecondNAglobal, occSecondNAfullGlobal$estimate)
expect_equal(length(occSecondNAfullGlobal$occupied)/as.numeric(length(hex)), occSecondNAfullGlobal$estimate)

expect_silent(occSecondNAfullRelative <- occupancy(matSecondNA, s=hex, full=TRUE, prop="relative"))
expect_equal(occSecondNArelative, occSecondNAfullRelative$estimate)
expect_equal(1, occSecondNAfullRelative$estimate)

### Last 
expect_silent(occLastNAfullGlobal <- occupancy(matLastNA, s=hex, full=TRUE, prop="global"))
expect_equal(occLastNAglobal, occLastNAfullGlobal$estimate)
expect_equal(length(occLastNAfullGlobal$occupied)/as.numeric(length(hex)), occLastNAfullGlobal$estimate)

expect_silent(occLastNAfullRelative <- occupancy(matLastNA, s=hex, full=TRUE, prop="relative"))
expect_equal(occLastNArelative, occLastNAfullRelative$estimate)
expect_equal(1, occLastNAfullRelative$estimate)

## 0-length data 
expect_silent(occZeroFullGlobal <- occupancy(matZero, s=hex, full=TRUE, prop="global"))
expect_equal(occZeroGlobal, occZeroFullGlobal$estimate)
expect_equal(0, occZeroFullGlobal$estimate)

expect_silent(occZeroFullRelative <- occupancy(matZero, s=hex, full=TRUE, prop="relative"))
expect_equal(occZeroRelative, occZeroFullRelative$estimate)
expect_equal(NA, occZeroFullRelative$estimate)

## Singular data
expect_silent(occSingleFullGlobal <- occupancy(matSingle, s=hex, full=TRUE, prop="global"))
expect_equal(occSingleGlobal, occSingleFullGlobal$estimate)
expect_equal(length(occSingleFullGlobal$occupied)/as.numeric(length(hex)), occSingleFullGlobal$estimate)

expect_silent(occSingleFullRelative <- occupancy(matSingle, s=hex, full=TRUE, prop="relative"))
expect_equal(occSingleRelative, occSingleFullRelative$estimate)
expect_equal(1, occSingleFullRelative$estimate)

## Duplicates=TRUE
expect_error(occupancy(mat, s=hex, duplicates=TRUE, full=TRUE, prop="global"))
expect_error(occupancy(mat, s=hex, duplicates=FALSE, full=TRUE, prop="global"))
expect_error(occupancy(mat, s=hex, duplicates=TRUE, full=TRUE, prop="relative"))
expect_error(occupancy(mat, s=hex, duplicates=FALSE, full=TRUE, prop="relative"))

## Plotting
# error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE, full=TRUE, prop="global"))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE, full=TRUE, prop="global"))

p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE, full=TRUE, prop="relative"))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(mat, hex, plot=TRUE, full=TRUE, prop="relative"))

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------

# wrong column names
expect_silent(occupancy(x=mat, s=hex, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_error(occupancy(x=mat, s=hex, long="sdf", prop="global"))
expect_error(occupancy(x=mat, s=hex, lat="sdf", prop="global"))
expect_error(occupancy(x=mat, s=hex, long="sdf",lat="sdf", prop="global"))
expect_silent(occupancy(x=mat, s=hex, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_error(occupancy(x=mat, s=hex, long="sdf", prop="relative"))
expect_error(occupancy(x=mat, s=hex, lat="sdf", prop="relative"))
expect_error(occupancy(x=mat, s=hex, long="sdf",lat="sdf", prop="relative"))

# not implemented yet
expect_error(occupancy(mat, hex, q=0.7, prop="global"))


################################################################################
# 2. Single data.frame methods
################################################################################


# manual checking
manualGlobal <- function(x, s){
	x <- x[, c("decimalLongitude", "decimalLatitude")]
	# get rid of missing vals
	cellsManual <- levels(factor((locate(s,x ))))
	length(cellsManual)/as.numeric(length(s))
}

manualRelative <- function(x, s){
	x <- x[, c("decimalLongitude", "decimalLatitude")]
	# get rid of missing vals
	cellsManual <- levels(factor((locate(s,x ))))
	length(cellsManual)/length(cellsManual)
}

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occCount <- occupancy(x=nobilis, s=hex, long="decimalLongitude",lat="decimalLatitude"))
expect_silent(occGlobal <- occupancy(x=nobilis, s=hex, long="decimalLongitude",lat="decimalLatitude",  prop="global"))
expect_equal(manualGlobal(x=nobilis, s=hex), occGlobal)
expect_silent(occRelative <- occupancy(x=nobilis, s=hex, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_equal(occRelative, 1)

## With missing values
### First
dfFirstNA <- nobilis
dfFirstNA[1, ] <- NA
expect_silent(occFirstNAglobal <- occupancy(dfFirstNA, s=hex, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_equal(occFirstNAglobal, manualGlobal(dfFirstNA, s=hex))
expect_silent(occFirstNArelative <- occupancy(dfFirstNA, s=hex, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_equal(occFirstNArelative, 1)

### Second
dfSecondNA <- nobilis
dfSecondNA[2, ] <- NA
expect_silent(occSecondNAglobal <- occupancy(dfSecondNA, long="decimalLongitude",lat="decimalLatitude", s=hex, prop="global"))
expect_equal(occSecondNAglobal, manualGlobal(dfSecondNA, s=hex))
expect_silent(occSecondNArelative <- occupancy(dfSecondNA, long="decimalLongitude",lat="decimalLatitude", s=hex, prop="relative"))
expect_equal(occSecondNArelative, 1)

### Last 
dfLastNA <- nobilis
dfLastNA[nrow(nobilis), ] <- NA
expect_silent(occLastNAglobal <- occupancy(dfLastNA, long="decimalLongitude",lat="decimalLatitude", s=hex, prop="global"))
expect_equal(occLastNAglobal, manualGlobal(dfLastNA, s=hex))
expect_silent(occLastNArelative <- occupancy(dfLastNA, long="decimalLongitude",lat="decimalLatitude", s=hex, prop="relative"))
expect_equal(occLastNArelative, 1)

## 0-length data 
dfZero <- nobilis[0,]
expect_silent(occZeroGlobal <- occupancy(dfZero, s=hex, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_equal(occZeroGlobal, 0)
expect_silent(occZeroRelative <- occupancy(dfZero, s=hex, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_equal(occZeroRelative, NA)

## Singular data
dfSingle <- nobilis[1,,drop=FALSE]
expect_silent(occSingleGlobal <- occupancy(dfSingle, s=hex, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_equal(occSingleGlobal, as.numeric((1/length(hex))))
expect_silent(occSingleRelative <- occupancy(dfSingle, s=hex, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_equal(occSingleRelative, 1)

## Duplicates=TRUE
expect_error(occupancy(nobilis, s=hex, duplicates=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_error(occupancy(nobilis, s=hex, duplicates=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_error(occupancy(nobilis, s=hex, duplicates=FALSE, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_error(occupancy(nobilis, s=hex, duplicates=FALSE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
 
## Plotting
# error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(nobilis, hex, plot=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(nobilis, hex, plot=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))

# error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(nobilis, hex, plot=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(nobilis, hex, plot=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))



#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occFullGlobal <- occupancy(x=nobilis, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_true(inherits(occFullGlobal, "orange"))
expect_equal(occGlobal, occFullGlobal$estimate)
expect_equal(occFullGlobal$estimate, length(occFullGlobal$occupied)/as.numeric(length(hex)))

expect_silent(occFullRelative <- occupancy(x=nobilis, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_true(inherits(occFullRelative, "orange"))
expect_equal(occRelative, occFullRelative$estimate)
expect_equal(occFullRelative$estimate, length(occFullRelative$occupied)/as.numeric(length(occFullRelative$occupied)))

### First
expect_silent(occFirstNAfullGlobal <- occupancy(dfFirstNA, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_equal(occFirstNAglobal, occFirstNAfullGlobal$estimate)
expect_equal(length(occFirstNAfullGlobal$occupied)/as.numeric(length(hex)), occFirstNAfullGlobal$estimate)

expect_silent(occFirstNAfullRelative <- occupancy(dfFirstNA, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_equal(occFirstNArelative, occFirstNAfullRelative$estimate)
expect_equal(1, occFirstNAfullRelative$estimate)

### Second
expect_silent(occSecondNAfullGlobal <- occupancy(dfSecondNA, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_equal(occSecondNAglobal, occSecondNAfullGlobal$estimate)
expect_equal(length(occSecondNAfullGlobal$occupied)/as.numeric(length(hex)), occSecondNAfullGlobal$estimate)

expect_silent(occSecondNAfullRelative <- occupancy(dfSecondNA, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_equal(occSecondNArelative, occSecondNAfullRelative$estimate)
expect_equal(1, occSecondNAfullRelative$estimate)

### Last 
expect_silent(occLastNAfullGlobal <- occupancy(dfLastNA, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_equal(occLastNAglobal, occLastNAfullGlobal$estimate)
expect_equal(length(occLastNAfullGlobal$occupied)/as.numeric(length(hex)), occLastNAfullGlobal$estimate)

expect_silent(occLastNAfullRelative <- occupancy(dfLastNA, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_equal(occLastNArelative, occLastNAfullRelative$estimate)
expect_equal(1, occLastNAfullRelative$estimate)

## 0-length data 
expect_silent(occZeroFullGlobal <- occupancy(dfZero, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_equal(occZeroGlobal, occZeroFullGlobal$estimate)
expect_equal(0, occZeroFullGlobal$estimate)

expect_silent(occZeroFullRelative <- occupancy(dfZero, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_equal(occZeroRelative, occZeroFullRelative$estimate)
expect_equal(NA, occZeroFullRelative$estimate)

## Singular data
expect_silent(occSingleFullGlobal <- occupancy(dfSingle, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_equal(occSingleGlobal, occSingleFullGlobal$estimate)
expect_equal(length(occSingleFullGlobal$occupied)/as.numeric(length(hex)), occSingleFullGlobal$estimate)

expect_silent(occSingleFullRelative <- occupancy(dfSingle, s=hex, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_equal(occSingleRelative, occSingleFullRelative$estimate)
expect_equal(1, occSingleFullRelative$estimate)

## Duplicates=TRUE
expect_error(occupancy(nobilis, s=hex, duplicates=TRUE, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_error(occupancy(nobilis, s=hex, duplicates=FALSE, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))
expect_error(occupancy(nobilis, s=hex, duplicates=TRUE, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))
expect_error(occupancy(nobilis, s=hex, duplicates=FALSE, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))

## Plotting
# error without something to draw on
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(nobilis, hex, plot=TRUE, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(nobilis, hex, plot=TRUE, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="global"))

p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(nobilis, hex, plot=TRUE, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(nobilis, hex, plot=TRUE, full=TRUE, long="decimalLongitude",lat="decimalLatitude", prop="relative"))

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
