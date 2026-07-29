# Coordinate occupancy of localities, iterated with tax

# Setup
library(tinytest)
library(orange)

# setup
data(pinna)

# use an icosahedral grid for localities (precursor to icosa-specific methods)
suppressMessages(hex<- hexagrid(spacing=2))
pinna$cell <- locate(hex, pinna[, c("decimalLongitude", "decimalLatitude")])

################################################################################
# 1. Data.frame with tax 
################################################################################

#-------------------------------------------------------------------------------
# With tax=NULL - wrong argumentation 
#-------------------------------------------------------------------------------
expect_error(occAll <- occupancy(pinna, s="cell", prop="global"))
expect_silent(occAll <- occupancy(pinna, s="cell", prop="relative"))

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

# function to calculate
manualEstimate <- function(x){
	manual_occup <- unique(x[, c("species", "cell")])
	# make sure that missing values are omitted
	manual_occup <- manual_occup[!is.na(manual_occup$species) & !is.na(manual_occup$cell), ]
	manual_occup_tab <- table(manual_occup$species)
	manu <- as.numeric(manual_occup_tab)
	names(manu) <- names(manual_occup_tab)
	# the total coccupancy
	manu <-manu/length(levels(factor(manual_occup$cell)))
	return(manu)
}

## Complete dataset
expect_error( occupancy(pinna, s="cell", tax="species", prop="global"))
expect_silent(occ <- occupancy(pinna, s="cell", tax="species", prop="relative"))
expect_true(inherits(occ, "numeric"))
expect_true(is.null(dim(occ)))
expect_equal(manualEstimate(pinna), occ)

## With missing values - first
naFirst <- pinna
naFirst[1:10, "cell"] <- NA
expect_error(occupancy(naFirst, s="cell", tax="species", prop="global"))
expect_silent(occFirst <- occupancy(naFirst, s="cell", tax="species", prop="relative"))
expect_equal(manualEstimate(naFirst), occFirst)

## With missing values - second
naSecond <- pinna
naSecond[11:20, "cell"] <- NA
expect_error(occupancy(naSecond, s="cell", tax="species", prop="global"))
expect_silent(occSecond <- occupancy(naSecond, s="cell", tax="species", prop="relative"))
expect_equal(manualEstimate(naSecond), occSecond)

## With missing values - last
naLast <- pinna
naLast[nrow(pinna), "cell"] <- NA
expect_error(occupancy(naLast, s="cell", tax="species", prop="global"))
expect_silent(occLast <- occupancy(naLast, s="cell", tax="species", prop="relative"))
expect_equal(manualEstimate(naLast), occLast)

## Singular data in some
naSingle <- pinna[1,]
expect_error(occupancy(naSingle, s="cell", tax="species", prop="global"))
expect_silent(occSingle <- occupancy(naSingle, s="cell", tax="species", prop="relative"))
expect_equal(manualEstimate(naSingle), occSingle)
expect_equivalent(1L, occSingle)

# 0-length
naZero <- pinna[0,]
expect_error(occupancy(naZero, s="cell", tax="species", prop="global"))
expect_silent(occZero <- occupancy(naZero, s="cell", tax="species", prop="relative"))
expect_equal(occZero, numeric())

## Duplicates=TRUE
expect_error(occupancy(pinna, s="cell", tax="species", duplicates=TRUE, prop="global"))
expect_error(occupancy(pinna, s="cell", tax="species", duplicates=TRUE, prop="relative"))

## Missing in taxon name 
missTax <- pinna
missTax$species[missTax$species=="Pinna rotundata"] <- NA
expect_error(occMissTax <- occupancy(missTax, s="cell", tax="species", prop="global"))
expect_silent(occMissTax <- occupancy(missTax, s="cell", tax="species", prop="relative"))
expect_equal(length(occMissTax)+1, length(occ))
expect_equal(occMissTax, occ[names(occ)!="Pinna rotundata"])

missTaxOne <- pinna
missTaxOne$species[1] <- NA
expect_silent(occMissTaxOne <- occupancy(missTaxOne, s="cell", tax="species", prop="relative"))
expect_equal(length(occMissTaxOne), length(occ))
expect_equal(manualEstimate(missTaxOne), occMissTaxOne)

#-------------------------------------------------------------------------------
# Full output/tracability, listout=TRUE
#-------------------------------------------------------------------------------

## Complete dataset
expect_error(occupancy(pinna, s="cell", tax="species", full=TRUE, prop="global"))
expect_silent(occFull <- occupancy(pinna, s="cell", tax="species", full=TRUE, prop="relative"))
expect_silent(occFullNoProp <- occupancy(pinna, s="cell", tax="species", full=TRUE))
expect_true(inherits(occFull, "array"))
expect_true(inherits(occFull[[1]], "orange"))
expect_equal(names(occFull), names(occ)) # data for same taxa
expect_equal(occ, sapply(occFull, function(b) b$estimate)) # results are the same 
expect_equal(lapply(occFullNoProp, function(b) b$occupied), lapply(occFull, function(b) b$occupied)) # occupied stuff are the same 

## With missing values - first
expect_error(occupancy(naFirst, s="cell", tax="species", full=TRUE, prop="global"))
expect_silent(occFullFirst <- occupancy(naFirst, s="cell", tax="species", full=TRUE, prop="relative"))
expect_silent(occFullFirstNoProp <- occupancy(naFirst, s="cell", tax="species", full=TRUE))
expect_true(inherits(occFullFirst[[1]], "orange"))
expect_equal(names(occFullFirst), names(occFirst)) # data for same taxa
expect_equal(occFirst, sapply(occFullFirst, function(b) b$estimate)) # results are the same 
expect_equal(lapply(occFullFirstNoProp, function(b) b$occupied), lapply(occFullFirst, function(b) b$occupied)) # occupied stuff are the same 

## With missing values - second
expect_error(occupancy(naSecond, s="cell", tax="species", full=TRUE, prop="global"))
expect_silent(occFullSecond <- occupancy(naSecond, s="cell", tax="species", full=TRUE, prop="relative"))
expect_silent(occFullSecondNoProp <- occupancy(naSecond, s="cell", tax="species", full=TRUE))
expect_true(inherits(occFullSecond[[1]], "orange"))
expect_equal(names(occFullSecond), names(occSecond)) # data for same taxa
expect_equal(occSecond, sapply(occFullSecond, function(b) b$estimate)) # results are the same 
expect_equal(lapply(occFullSecondNoProp, function(b) b$occupied), lapply(occFullSecond, function(b) b$occupied)) # occupied stuff are the same 

## With missing values - last
expect_error(occupancy(naLast, s="cell", tax="species", full=TRUE, prop="global"))
expect_silent(occFullLast <- occupancy(naLast, s="cell", tax="species", full=TRUE, prop="relative"))
expect_silent(occFullLastNoProp <- occupancy(naLast, s="cell", tax="species", full=TRUE))
expect_true(inherits(occFullSecond[[1]], "orange"))
expect_equal(names(occFullLast), names(occLast)) # data for same taxa
expect_equal(occLast, sapply(occFullLast, function(b) b$estimate)) # results are the same 
expect_equal(lapply(occFullLastNoProp, function(b) b$occupied), lapply(occFullLast, function(b) b$occupied)) # occupied stuff are the same 

## 0-length data 
expect_error(occupancy(naZero, s="cell", tax="species", full=TRUE, prop="global"))
expect_silent(occZero <- occupancy(naZero, s="cell", tax="species", full=TRUE, prop="relative"))
expect_equal(length(occZero), 0L)

## Singular data
expect_error(occupancy(naSingle, s="cell", tax="species", full=TRUE, prop="global"))
expect_silent(occSingle <- occupancy(naSingle, s="cell", tax="species", full=TRUE, prop="relative"))
expect_true(inherits(occSingle[[1]], "orange"))
expect_equal(length(occSingle), 1L)
expect_equal(occSingle[[1]]$estimate, 1L)
expect_equal(occSingle[[1]]$estimate, length(occSingle[[1]]$occupied))


#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------

expect_error(occupancy(pinna, s="cell", tax="wrong", full=TRUE, prop="global"))
expect_error(occupancy(pinna, s="cell", tax="wrong", full=TRUE, prop="relative"))

## Plotting
expect_error(occPlot <- occupancy(pinna, s="cell", tax="species",  plot=TRUE))

################################################################################
# 2. sf data with tax 
################################################################################
