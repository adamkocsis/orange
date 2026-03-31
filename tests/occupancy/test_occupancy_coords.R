# in data frame
# occupancy: x: data.frame, sf: missing, icosa: missing
library(tinytest)
library(divDyn)


data(corals)

################################################################################
# Occupancy based on coordinates
 
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
expect_silent(occup_coll <- occupancy(corals, tax="genus", loc="collection_no"))
expect_equal(class(occup_coll), "numeric")

# expect errors
expect_error(occup_coll <- occupancy(corals))

