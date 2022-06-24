################################################################################ .
# The code is used for reproducing Figure 2 in the main body of the manuscript   #
# Input needed:                                                                  #
# 1_Vertebrates_species_richness_grouped_by_IUCN_category (FOLDER)               #
# 2_Ensemble_weighted_global_forest_cover (FOLDER)                               #
# 4_All_correlations.csv                                                         #
# 5_All_correlations_emissions.csv                                               #
# The above data were publicly available at GitHub                               #
################################################################################ .
# ------------------------------------------------------------------------------ ----  
# (1)  Function                                                                  ----

# Set the work directory
setwd("D:/OneDrive/SCI_Publications/Extinction_Debt_In_draft/DATA AVAILABILITY")                 

# Load the necessary packages
library(raster)

# For calculating Pearson's and the corresponding P values, we create a new function called "Cor_Analysis"

Cor_Analysis <- function(Taxon = c("Amphibians", "Mammals", "Reptiles"),
                         redlistCategory = c("All", "Critically_Endangered", "Endangered", 
                                        "Lease_Concern", "Near_Threatened", "Vulnerable"),
                         Resolution = c(0.5,1,1.5,2),
                         Methods = c("Forward", "Backward", "JP", "Ensemble"),
                         Year = seq(1500,1992,1)){
  
  # Load the species richness for giving Taxon, IUCN red list category, spatial resolution, 
  # year of the historical forest cover and the reconstruction method 

  Animal_SR <- raster(paste0("./1_Vertebrates_species_richness_grouped_by_IUCN_category/", 
                             redlistCategory,"_Forest_", Taxon,"_",Resolution,".tif"))
  
  # Read the forest cover raster
  if(Methods == "Forward"| Methods == "Backward"| Methods == "Ensemble"){
    Forest_Cover <- raster(list.files("./2_Ensemble_weighted_global_forest_cover/", 
                                      pattern = paste0(Methods,"_M5_", Year,".tif"), full.names = T))
    } else if (Methods == "JP"){
    Forest_Cover <- raster(list.files("./2_Ensemble_weighted_global_forest_cover/", 
                                      pattern = paste0("Global_forest_",Methods,"_", Year,".tif"), full.names = T))  
  }
  
  # Extract values from the above two layers
  coords <- coordinates(Animal_SR)
  x <- values(Animal_SR)
  y <- values(Forest_Cover)
  anaysis_table <- as.data.frame(cbind("Animal_SR" = x, "Forest_Cover" = y,"Coords" = coords))
  anaysis_table <- anaysis_table[complete.cases(anaysis_table),]
  
  # Calculate the Pearson's correlation coefficient
  w <- cor.test(anaysis_table[, 1], anaysis_table[,2])
  w
  
  # Store the result
  Res <- data.frame("Taxon" = Taxon, "redlistCategory" = redlistCategory, 'Resolution'= Resolution,"Methods" = Methods, 
                    'Year' = Year, "Pearson R" = round(w[["estimate"]][["cor"]],3),"P value" = w[["p.value"]])
  
  # Return result
  return(Res)
}

# Please note, we only provide the forest cover data reconstructed using ensemble weighted method
# Run the function and save the data as csv file
# e.g. 0.5 degree spatial resolution
Cor_Analysis(Taxon = "Amphibians", redlistCategory = "All", Resolution = 0.5, Methods = "Ensemble", Year = 1900)
Cor_Analysis(Taxon = "Mammals", redlistCategory = "Critically_Endangered", Resolution = 0.5, Methods = "Ensemble", Year = 1900)
Cor_Analysis(Taxon = "Reptiles", redlistCategory = "Endangered", Resolution = 0.5, Methods = "Ensemble", Year = 1500)

# ------------------------------------------------------------------------------ ----  
# (2)  Figure 2 A                                                                ----

# Load all the correlations produced by "Cor_Analysis" function
data1 <- read.csv("./4_All_correlations.csv", header = T)
names(data1)
summary(data1)
table(data1$Taxon)
table(data1$redlistCategory)

data1$redlistCategory <- factor(data1$redlistCategory, 
                                   levels = c("All species",
                                              "Least Concern",
                                              "Vulnerable",
                                              "Near Threatened",
                                              "Endangered",
                                              "Critically Endangered"))
data1$Taxon <- factor(data1$Taxon, 
                             levels = c("Mammals",
                                        "Amphibians",
                                        "Reptiles"))

data1$Methods <- factor(data1$Methods, 
                       levels = c("Ensemble",
                                  "Backward",
                                  "Forward",
                                  "JP"))
# install.packages("ggplot2")
library(ggplot2)

# install.packages("ggsci")
library(ggsci) # Scientific Journal and Sci-Fi Themed Color Palettes for 'ggplot2'

p1 <- ggplot(data1, aes(x = Year, y = Pearson.s.r, color = Taxon)) +
  labs(x = "Year", y = "Pearson's r") +
  facet_grid(Methods~redlistCategory) +
  geom_line() +
  scale_y_continuous(breaks=seq(0, 0.6, 0.1)) +
  scale_color_aaas() +
  theme_test()

p1 + geom_vline(aes(xintercept = 1784), colour = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = 1870), colour = "coral", linetype = "dashed") + 
  geom_vline(aes(xintercept = 1969), colour = "blue", linetype = "dashed") 

# ------------------------------------------------------------------------------ ----  
# (3)  Figure 2 B                                                                ----

# create all the combinations of variables
# Comb_var <- expand.grid(c("Amphibians", "Mammals", "Reptiles"),
#                         c("All", "Critically_Endanged", "Endanged", "Lease_Concern", "Near_Threatened", "Vulnerable"),
#                        c("Forward", "Backward", "JP"), 
#                        seq(1500,1992,1))
# names(Comb_var) <- c("Taxon", "redlistCategory", "Methods", "Years")

# a data frame for recording all results
Res <- data.frame("Taxon" = NULL, "redlistCategory" = NULL, 'Resolution'= NULL,
                  "Methods" = NULL, 'Year' = NULL, "Def.R" = NULL)
# loop
for (Taxon in c("Mammals", "Amphibians", "Reptiles")) {
  for (redlistCategory in c("All species", "Critically Endangered", "Endangered", 
                       "Least Concern", "Near Threatened", "Vulnerable")) {
    for (Methods in c("Ensemble", "Backward","Forward", "JP")) {
        data_sub <- data1[which(data1$Methods == Methods & data1$Taxon == Taxon & data1$redlistCategory == redlistCategory), ]
     
        for (Year in seq(1510,1990,10)) {
        Def.R <- data_sub[which(data_sub$Year == Year),"Pearson.s.r"]- data_sub[which(data_sub$Year == Year) - 10, "Pearson.s.r"]
        res <- data.frame("Taxon" = Taxon, "redlistCategory" = redlistCategory, 'Resolution'= 0.5,
                           "Methods" = Methods, 'Year' = Year, "Def.R" = round(Def.R,3))
        
        Res <- rbind(Res,res)
        cat("Taxon: ", Taxon, " redlistCategory: ",redlistCategory," Methods: ", Methods, " Year: ", Year, " Finished...", "\n")   
      }
    }
  }
}

# write.csv(Res,"Rep_all.csv")


Res$redlistCategory <- factor(Res$redlistCategory, levels = c("All species",
                                                              "Least Concern",
                                                              "Vulnerable",
                                                              "Near Threatened",
                                                              "Endangered",
                                             "Critically Endangered"))
Res$Taxon <- factor(Res$Taxon, levels = c("Mammals",
                                       "Amphibians",
                                       "Reptiles"))

Res$Methods <- factor(Res$Methods, 
                      levels = c("Ensemble",
                                 "Backward",
                                 "Forward",
                                 "JP"))


p2 <- ggplot(Res, aes(x = Year, y = Def.R, color = Taxon)) +
  labs(x = "Year", y = "Pearson's r") +
  facet_grid(Methods ~ redlistCategory) +
  geom_line() +
  scale_y_continuous(breaks=seq(-0.03, 0, 0.005)) +
  scale_color_aaas() +
  theme_test()

p2 + geom_vline(aes(xintercept = 1784), colour = "black", linetype = "dashed") +
  geom_vline(aes(xintercept = 1870), colour = "coral", linetype = "dashed") + 
  geom_vline(aes(xintercept = 1969), colour = "blue", linetype = "dashed")

# ------------------------------------------------------------------------------ ----  
# (3)  Figure 2 C                                                                ----

# Now the MTCO2e values were added into the correlation table
data2 <- read.csv("5_All_correlations_emissions.csv", header = T)
names(data2)
summary(data2)


data2$redlistCategory <- factor(data2$redlistCategory, levels = c("All species",
                                                                "Least Concern",
                                                                "Vulnerable",
                                                                "Near Threatened",
                                                                "Endangered",
                                                                "Critically Endangered"))
data2$Taxon <- factor(data2$Taxon, levels = c("Mammals",
                                            "Amphibians",
                                            "Reptiles"))

data2$Methods <- factor(data2$Methods, levels = c("Ensemble",
                                                "Backward",
                                                "Forward",
                                                "JP"))

library(ggpmisc)
formula <- y ~ x

p3 <- ggplot(data2, aes(x = MTCO2e, y = Pearson.s.r, color = Taxon, shape = Taxon )) +
  labs(x = expression("MTCO"["2"]~"e"), y = "Pearson's r") +
  facet_grid(Methods ~ redlistCategory) +
  # facet_wrap(~Conservation.status, scales = "free_y") +
  # geom_line() +
  geom_point(alpha = 0.1, size =2) +
  scale_y_continuous(breaks=seq(0, 0.7, 0.1)) +
  scale_x_continuous(breaks =seq(0,70000,10000), labels = c(0,1,2,3,4,5,6,7))+
  scale_color_aaas() +
  stat_smooth(method="lm", formula = formula, se = F)+
  stat_poly_eq(aes(label = paste(stat(p.value.label),..rr.label.., sep = "~~~")), 
              formula = formula, parse = TRUE, size = 3) +
  scale_shape_manual(values = c(6,1,2))+ 
  theme_test() 

p3


# ------------------------------------------------------------------------------ ---- 
# EOF