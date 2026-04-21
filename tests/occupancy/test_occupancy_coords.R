# in data frame
# occupancy: x: data.frame, sf: missing, icosa: missing
library(tinytest)
suppressPackageStartupMessages(library(divDyn))


################################################################################
# 1. Single-taxon dataset
################################################################################
data(pinna)

# coordinate pairs
expect_silent(occ <- occupancy(pinna, long="decimallongitude", lat="decimallatitude"))

# manually
occRows <- nrow(unique(pinna[, c(long="decimallongitude", lat="decimallatitude")]))
expect_equal(occ, occRows)

# simple matrix solution
expect_silent(occMat <- occupancy(as.matrix(pinna[,c("decimallongitude", "decimallatitude")])))
expect_equal(occMat, occRows)


################################################################################
# 2. Multi-taxon dataset
################################################################################

################################################################################
# Occupancy based on coordinates, multiple taxa
data(corals)
 
# occupancy in the corals dataset
expect_silent(occup <- occupancy(corals, tax="genus", long="lng", lat="lat"))
expect_equal(class(occup), "numeric")

# manually calculated
manual_occup <- unique(corals[, c("genus", "lng", "lat")])
manual_occup <- table(manual_occup$genus)
manu <- as.numeric(manual_occup)
names(manu) <- names(manual_occup)
expect_equal(manu, occup)

################################################################################
# Occupancy based on collections
expect_error(occup_coll <- occupancy(corals, tax="genus", s="WRONG"))
expect_silent(occup_coll <- occupancy(corals, tax="genus", s="collection_no"))
expect_equal(class(occup_coll), "numeric")

# manually recreat
core<- corals[, c("genus", "collection_no")]
core <- unique(core[!is.na(core$genus) & !is.na(core$collection_no), ])
tab <- table(core$genus)
tabNum <- as.numeric(tab)
names(tabNum) <- names(tab)
expect_equal(occup_coll, tabNum)

# expect errors, when wrong long-lat is specified
expect_error(occup_coll <- occupancy(corals))
