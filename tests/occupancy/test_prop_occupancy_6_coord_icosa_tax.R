# Proportional Occupancy tests
# x: points, s: icosa, iterated with tax

# Setup
library(tinytest)
library(orange)

# setup
data(pinna)

# mke a hexagrid
suppressMessages(hex <- hexagrid(deg=2, sf=TRUE))

################################################################################
# 1A. Data.frame with tax  - Global
################################################################################

#-------------------------------------------------------------------------------
# With tax=NULL - wrong argumentation 
#-------------------------------------------------------------------------------
expect_silent(occAll <- occupancy(pinna, long="decimalLongitude", lat="decimalLatitude", s=hex, prop="global"))

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

# function to calculate - as different as possilbe..
manualEstimateGlobal <- function(x, s){
	manual_occup <- unique(x[, c("species", "decimalLongitude", "decimalLatitude")])
	manual_occup$cell <- locate(s, manual_occup[, c("decimalLongitude", "decimalLatitude")])

	# create contiainer
	ta <-levels(factor(manual_occup$species))
	res <- rep(NA, length(ta))
	names(res) <- ta

	# for loop to make sure!
	for(i in 1:length(ta)){
		thisTa <- manual_occup[which(manual_occup$species==ta[i]), "cell"]
		# omit missing values
		thisTa <- unique(thisTa[!is.na(thisTa)])
		res[i] <- length(thisTa)

	}
	res <- res/as.numeric(length(s))

	return(res)
}

## Complete dataset
expect_silent(occ <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude",tax="species", prop="global"))
expect_true(inherits(occ, "numeric"))
expect_true(is.null(dim(occ)))
expect_equal(occ, manualEstimateGlobal(x=pinna, s=hex))

## With missing values - first
naFirst <- pinna
naFirst[1:10, "cell"] <- NA
expect_silent(occFirst <- occupancy(naFirst, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="global"))
expect_equal(manualEstimateGlobal(x=naFirst, s=hex), occFirst)

## With missing values - second
naSecond <- pinna
naSecond[11:20, "cell"] <- NA
expect_silent(occSecond <- occupancy(naSecond, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="global"))
expect_equal(manualEstimateGlobal(naSecond, s=hex), occSecond)

## With missing values - last
naLast <- pinna
naLast[nrow(pinna), "cell"] <- NA
expect_silent(occLast <- occupancy(naLast, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="global"))
expect_equal(manualEstimateGlobal(naLast, s=hex), occLast)

## Singular data in some
naSingle <- pinna[1,]
expect_silent(occSingle <- occupancy(naSingle, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="global"))
expect_equal(manualEstimateGlobal(naSingle, s=hex), occSingle)
expect_equivalent(1L/as.numeric(length(hex)), occSingle)

# 0-length
naZero <- pinna[0,]
expect_silent(occZero <- occupancy(naZero, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="global"))
expect_equal(occZero, numeric())

## Duplicates=TRUE
expect_error(occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", duplicates=TRUE, prop="global"))
## Missing in taxon name 
missTax <- pinna
missTax$species[missTax$species=="Pinna rotundata"] <- NA
expect_silent(occMissTax <- occupancy(missTax, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="global"))
expect_equal(length(occMissTax)+1, length(occ))
expect_equal(occMissTax, occ[names(occ)!="Pinna rotundata"])

missTaxOne <- pinna
missTaxOne$species[1] <- NA
expect_silent(occMissTaxOne <- occupancy(missTaxOne, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="global"))
expect_equal(length(occMissTaxOne), length(occ))
 
## Plotting
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", plot=TRUE, prop="global"))
expect_equal(occPlot, occ)

# plot is still different only last one shown!


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", plot=TRUE, prop="global"))


#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------
auxocc <- function(b) length(b$occupied)/as.numeric(length(hex))

## Complete dataset
expect_silent(occFull <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", full=TRUE, prop="global"))
expect_true(inherits(occFull, "array"))
expect_equal(names(occFull), names(occ)) # data for same taxa
expect_equal(occ, sapply(occFull,auxocc )) # results can be recreated 

## With missing values - first
expect_silent(occFullFirst <- occupancy(naFirst, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", full=TRUE, prop="global"))
expect_equal(names(occFullFirst), names(occFirst)) # data for same taxa
expect_equal(occFirst, sapply(occFullFirst, auxocc)) # results can be recreated 

## With missing values - second
expect_silent(occFullSecond <- occupancy(naSecond, s=hex, long="decimalLongitude", lat="decimalLatitude",tax="species", full=TRUE,prop="global"))
expect_equal(names(occFullSecond), names(occSecond)) # data for same taxa
expect_equal(occSecond, sapply(occFullSecond, auxocc)) # results can be recreated 

## With missing values - last
expect_silent(occFullLast <- occupancy(naLast, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", full=TRUE, prop="global"))
expect_equal(names(occFullLast), names(occLast)) # data for same taxa
expect_equal(occLast, sapply(occFullLast, auxocc )) # results can be recreated 

## 0-length data 
expect_silent(occZero <- occupancy(naZero, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", full=TRUE, prop="global"))
expect_equal(length(occZero), 0L)

## Singular data
expect_silent(occSingle <- occupancy(naSingle, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", full=TRUE, prop="global"))
expect_equal(length(occSingle), 1L)
expect_equal(occSingle[[1]]$estimate, 1/as.numeric(length(hex)))
expect_equal(occSingle[[1]]$estimate, length(occSingle[[1]]$occupied)/as.numeric(length(hex)))

## Plotting
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlotFull <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", plot=TRUE, full=TRUE, prop="global"))
expect_equal(occPlotFull, occFull)


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", plot=TRUE, prop="global"))

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------

expect_error(occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude",  tax="wrong", full=TRUE, prop="global"))

################################################################################
# 1B. Data.frame with tax - Relative proportions
################################################################################

#-------------------------------------------------------------------------------
# With tax=NULL - wrong argumentation 
#-------------------------------------------------------------------------------
expect_silent(occAll <- occupancy(pinna, long="decimalLongitude", lat="decimalLatitude", s=hex, prop="relative"))

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

# function to calculate - as different as possilbe..
manualEstimateRelative <- function(x, s){
	manual_occup <- unique(x[, c("species", "decimalLongitude", "decimalLatitude")])
	manual_occup$cell <- locate(s, manual_occup[, c("decimalLongitude", "decimalLatitude")])

	# create contiainer
	ta <-levels(factor(manual_occup$species))
	res <- rep(NA, length(ta))
	names(res) <- ta

	# for loop to make sure!
	for(i in 1:length(ta)){
		thisTa <- manual_occup[which(manual_occup$species==ta[i]), "cell"]
		# omit missing values
		thisTa <- unique(thisTa[!is.na(thisTa)])
		res[i] <- length(thisTa)

	}
	res <- res
	allCells <-levels(factor(manual_occup$cell))
	res <- res/length(allCells)

	return(res)
}

## Complete dataset
expect_silent(occ <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude",tax="species", prop="relative"))
expect_true(inherits(occ, "numeric"))
expect_true(is.null(dim(occ)))
expect_equal(occ, manualEstimateRelative(x=pinna, s=hex))

## With missing values - first
naFirst <- pinna
naFirst[1:10, "cell"] <- NA
expect_silent(occFirst <- occupancy(naFirst, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="relative"))
expect_equal(manualEstimateRelative(x=naFirst, s=hex), occFirst)

## With missing values - second
naSecond <- pinna
naSecond[11:20, "cell"] <- NA
expect_silent(occSecond <- occupancy(naSecond, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="relative"))
expect_equal(manualEstimateRelative(naSecond, s=hex), occSecond)

## With missing values - last
naLast <- pinna
naLast[nrow(pinna), "cell"] <- NA
expect_silent(occLast <- occupancy(naLast, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="relative"))
expect_equal(manualEstimateRelative(naLast, s=hex), occLast)

## Singular data in some
naSingle <- pinna[1,]
expect_silent(occSingle <- occupancy(naSingle, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="relative"))
expect_equal(manualEstimateRelative(naSingle, s=hex), occSingle)
expect_equivalent(1L, occSingle)

# 0-length
naZero <- pinna[0,]
expect_silent(occZero <- occupancy(naZero, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="relative"))
expect_equal(occZero, numeric())

## Duplicates=TRUE
expect_error(occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", duplicates=TRUE, prop="relative"))
## Missing in taxon name 
missTax <- pinna
missTax$species[missTax$species=="Pinna rotundata"] <- NA
expect_silent(occMissTax <- occupancy(missTax, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="relative"))
expect_equal(length(occMissTax)+1, length(occ))
expect_equal(occMissTax, manualEstimateRelative(missTax, hex))

missTaxOne <- pinna
missTaxOne$species[1] <- NA
expect_silent(occMissTaxOne <- occupancy(missTaxOne, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", prop="relative"))
expect_equal(length(occMissTaxOne), length(occ))
expect_equal(occMissTaxOne, manualEstimateRelative(missTaxOne, hex))
 
## Plotting
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlot <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", plot=TRUE, prop="relative"))
expect_equal(occPlot, occ)

# plot is still different only last one shown!


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", plot=TRUE, prop="relative"))


#-------------------------------------------------------------------------------
# Full output/tracability 
#-------------------------------------------------------------------------------
auxocc <- function(b) b$estimate

## Complete dataset
expect_silent(occFull <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", full=TRUE, prop="relative"))
expect_true(inherits(occFull, "array"))
expect_equal(names(occFull), names(occ)) # data for same taxa
expect_equal(occ, sapply(occFull,auxocc )) # results can be recreated 

## With missing values - first
expect_silent(occFullFirst <- occupancy(naFirst, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", full=TRUE, prop="relative"))
expect_equal(names(occFullFirst), names(occFirst)) # data for same taxa
expect_equal(occFirst, sapply(occFullFirst, auxocc)) # results can be recreated 

## With missing values - second
expect_silent(occFullSecond <- occupancy(naSecond, s=hex, long="decimalLongitude", lat="decimalLatitude",tax="species", full=TRUE,prop="relative"))
expect_equal(names(occFullSecond), names(occSecond)) # data for same taxa
expect_equal(occSecond, sapply(occFullSecond, auxocc)) # results can be recreated 

## With missing values - last
expect_silent(occFullLast <- occupancy(naLast, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", full=TRUE, prop="relative"))
expect_equal(names(occFullLast), names(occLast)) # data for same taxa
expect_equal(occLast, sapply(occFullLast, auxocc )) # results can be recreated 

## 0-length data 
expect_silent(occZero <- occupancy(naZero, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", full=TRUE, prop="relative"))
expect_equal(length(occZero), 0L)

## Singular data
expect_silent(occSingle <- occupancy(naSingle, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", full=TRUE, prop="relative"))
expect_equal(length(occSingle), 1L)
expect_equal(occSingle[[1]]$estimate, 1L)

## Plotting
p <- dev.cur()
while(p!=1){
	dev.off()
	p <- dev.cur()
}
expect_silent(occPlotFull <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", plot=TRUE, full=TRUE, prop="relative"))
expect_equal(occPlotFull, occFull)


# normal plotting
plot(NULL, xlim=c(-180, 180), ylim=c(-90, 90))
expect_silent(occPlot <- occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude", tax="species", plot=TRUE, prop="relative"))

#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------

expect_error(occupancy(pinna, s=hex, long="decimalLongitude", lat="decimalLatitude",  tax="wrong", full=TRUE, prop="relative"))


################################################################################
# 2. sf data with tax 
################################################################################
