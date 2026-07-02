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
expect_silent(occAll <- occupancy(pinna, s="cell"))

#-------------------------------------------------------------------------------
# Estimate only
#-------------------------------------------------------------------------------

# function to calculate
manualEstimate <- function(x){
	manual_occup <- unique(x[, c("species", "cell")])
	manual_occup <- manual_occup[!is.na(manual_occup$species) & !is.na(manual_occup$cell), ]
	manual_occup <- table(manual_occup$species)
	manu <- as.numeric(manual_occup)
	names(manu) <- names(manual_occup)
	return(manu)
}

## Complete dataset
expect_silent(occ <- occupancy(pinna, s="cell", tax="species"))
expect_true(inherits(occ, "numeric"))
expect_true(is.null(dim(occ)))
expect_equal(manualEstimate(pinna), occ)

## With missing values - first
naFirst <- pinna
naFirst[1:10, "cell"] <- NA
expect_silent(occFirst <- occupancy(naFirst, s="cell", tax="species"))
expect_equal(manualEstimate(naFirst), occFirst)

## With missing values - second
naSecond <- pinna
naSecond[11:20, "cell"] <- NA
expect_silent(occSecond <- occupancy(naSecond, s="cell", tax="species"))
expect_equal(manualEstimate(naSecond), occSecond)

## With missing values - last
naLast <- pinna
naLast[nrow(pinna), "cell"] <- NA
expect_silent(occLast <- occupancy(naLast, s="cell", tax="species"))
expect_equal(manualEstimate(naLast), occLast)

## Singular data in some
naSingle <- pinna[1,]
expect_silent(occSingle <- occupancy(naSingle, s="cell", tax="species"))
expect_equal(manualEstimate(naSingle), occSingle)
expect_equivalent(1L, occSingle)

# 0-length
naZero <- pinna[0,]
expect_silent(occZero <- occupancy(naZero, s="cell", tax="species"))
expect_equal(occZero, numeric())

## Duplicates=TRUE
expect_error(occupancy(pinna, s="cell", tax="species", duplicates=TRUE))
## Missing in taxon name 
missTax <- pinna
missTax$species[missTax$species=="Pinna rotundata"] <- NA
expect_silent(occMissTax <- occupancy(missTax, s="cell", tax="species"))
expect_equal(length(occMissTax)+1, length(occ))
expect_equal(occMissTax, occ[names(occ)!="Pinna rotundata"])

missTaxOne <- pinna
missTaxOne$species[1] <- NA
expect_silent(occMissTaxOne <- occupancy(missTaxOne, s="cell", tax="species"))
expect_equal(length(occMissTaxOne), length(occ))

#-------------------------------------------------------------------------------
# Full output/tracability, listout=TRUE
#-------------------------------------------------------------------------------

## Complete dataset
expect_silent(occFull <- occupancy(pinna, s="cell", tax="species", full=TRUE))
expect_true(inherits(occFull, "array"))
expect_equal(names(occFull), names(occ)) # data for same taxa
expect_equal(occ, sapply(occFull, function(b) length(b$occupied))) # results can be recreated 

## With missing values - first
expect_silent(occFullFirst <- occupancy(naFirst, s="cell", tax="species", full=TRUE))
expect_equal(names(occFullFirst), names(occFirst)) # data for same taxa
expect_equal(occFirst, sapply(occFullFirst, function(b) length(b$occupied))) # results can be recreated 

## With missing values - second
expect_silent(occFullSecond <- occupancy(naSecond, s="cell", tax="species", full=TRUE))
expect_equal(names(occFullSecond), names(occSecond)) # data for same taxa
expect_equal(occSecond, sapply(occFullSecond, function(b) length(b$occupied))) # results can be recreated 

## With missing values - last
expect_silent(occFullLast <- occupancy(naLast, s="cell", tax="species", full=TRUE))
expect_equal(names(occFullLast), names(occLast)) # data for same taxa
expect_equal(occLast, sapply(occFullLast, function(b) length(b$occupied))) # results can be recreated 

## 0-length data 
expect_silent(occZero <- occupancy(naZero, s="cell", tax="species", full=TRUE))
expect_equal(length(occZero), 0L)

## Singular data
expect_silent(occSingle <- occupancy(naSingle, s="cell", tax="species", full=TRUE))
expect_equal(length(occSingle), 1L)
expect_equal(occSingle[[1]]$estimate, 1L)
expect_equal(occSingle[[1]]$estimate, length(occSingle[[1]]$occupied))


#-------------------------------------------------------------------------------
# Wrong argumnets, appropriate defaults, plotting
#-------------------------------------------------------------------------------

expect_error(occupancy(pinna, s="cell", tax="wrong", full=TRUE))

################################################################################
# 2. sf data with tax 
################################################################################
