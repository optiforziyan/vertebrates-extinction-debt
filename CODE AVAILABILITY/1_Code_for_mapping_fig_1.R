################################################################################ .
# The code is used for reproducing Figure 1 in the main body of the manuscript   #
# Input needed:                                                                  #
# 1_Vertebrates_species_richness_grouped_by_IUCN_category (FOLDER)               #
# 2_Ensemble_weighted_global_forest_cover (FOLDER)                               #
# 3_Global_list_of_forest_dependent_terrestrial_vertebrates.xlsx                 #
# The above data were publicly available at GitHub                               #
################################################################################ .
# ------------------------------------------------------------------------------ ----  
# (1)  Figure 1                                                                  ----

# Set the work directory
setwd("D:/OneDrive/SCI_Publications/Extinction_Debt_In_draft/DATA AVAILABILITY")                 

# Load the necessary packages
library(raster)

# Load the forest cover in 1500 CE (Ensemble approach)
forest_1500 <- raster("./2_Ensemble_weighted_global_forest_cover/Ensemble_M5_1500.tif")
forest_1992 <- raster("./2_Ensemble_weighted_global_forest_cover/Ensemble_M5_1992.tif")

# Load the species richness for all mammals, amphibians and reptiles 
# e.g. 0.5 degree spatial resolution
mammals_sr    <- raster("./1_Vertebrates_species_richness_grouped_by_IUCN_category/All_Forest_Mammals_0.5.tif")
amphibians_sr <- raster("./1_Vertebrates_species_richness_grouped_by_IUCN_category/All_Forest_Amphibians_0.5.tif")
reptiles_sr   <- raster("./1_Vertebrates_species_richness_grouped_by_IUCN_category/All_Forest_reptiles_0.5.tif")  

# set robinson projection
robinson <- CRS("+proj=robin +over")

# plot
par(mfcol=c(2,3))
plot(projectRaster(forest_1500, crs = as.character(robinson)), main = "(A) Estimated forest cover in 1500 CE", axes = FALSE)
plot(projectRaster(forest_1992, crs = as.character(robinson)), main = "(B) Estimated forest cover in 1900 CE", axes = FALSE)
plot(projectRaster(forest_1992, crs = as.character(robinson)) - projectRaster(forest_1500, crs = as.character(robinson)), main = "(C) Percent gross forest cover changes (1500 CE - 1992 CE)", axes = FALSE)
plot(projectRaster(mammals_sr, crs = as.character(robinson)), main = "(D) Contemporary mammal species richness", axes = FALSE)
plot(projectRaster(amphibians_sr, crs = as.character(robinson)), main = "(E) Contemporary amphibian species richness", axes = FALSE)
plot(projectRaster(reptiles_sr, crs = as.character(robinson)), main = "(F) Contemporary reptile species richness", axes = FALSE)

# ------------------------------------------------------------------------------ ---- 
# EOF