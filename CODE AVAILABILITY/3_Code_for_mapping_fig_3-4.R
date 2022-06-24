################################################################################ .
# The code is used for reproducing Figure 3 & 4 in the main body of the ms       #
# Input needed:                                                                  #
# 1_Vertebrates_species_richness_grouped_by_IUCN_category (FOLDER)               #
# 2_Ensemble_weighted_global_forest_cover (FOLDER)                               #
# 6_Other_necessary_input_data (FOLDER)                                          #
# Please note that according to the terms of use of WDPA and RESOLVE             #
# we cannot directly share global protected areas data                           #
# However, researchers interested in our analysis can replicate all the results  #
# themselves using the present code after downloading the date from WDPA         #                                                  #
################################################################################ .
# ------------------------------------------------------------------------------ ---- 
# (1)  WDPA Analysis                                                             ----

# load packages
library(sf)
library(raster)
library(rgeos)
library(rgdal)

# load the polygon shapefile of WDPA
part1 <- st_read("./Input/Protected_areas/WDPA_Jan2022_Public_shp/WDPA_Jan2022_Public_shp_0/WDPA_Jan2022_Public_shp-polygons.shp")
class(part1)
part1 <- as(part1, "Spatial")  
class(part1)

part2 <- st_read("./Input/Protected_areas/WDPA_Jan2022_Public_shp/WDPA_Jan2022_Public_shp_1/WDPA_Jan2022_Public_shp-polygons.shp")
class(part2)
part2 <- as(part2, "Spatial")  
class(part2)

part3 <- st_read("./Input/Protected_areas/WDPA_Jan2022_Public_shp/WDPA_Jan2022_Public_shp_2/WDPA_Jan2022_Public_shp-polygons.shp")
class(part3)
part3 <- as(part3, "Spatial")  
class(part3)

# add Coordinate Reference System (CRS) projection. 
myCRS1 <- CRS("+init=epsg:4326") # WGS 84
crs(part1) <- myCRS1
crs(part2) <- myCRS1
crs(part3) <- myCRS1

# Extract PAs before 1992 CE
period_2_p1 <- part1[which(part1$STATUS_YR <= 1992 & part1$STATUS_YR > 1754), ]
period_2_p2 <- part2[which(part2$STATUS_YR <= 1992 & part2$STATUS_YR > 1754), ]
period_2_p3 <- part3[which(part3$STATUS_YR <= 1992 & part3$STATUS_YR > 1754), ]
b <- bind(period_2_p1, period_2_p2, period_2_p3)
b <- b[which(b$MARINE == 0),]
head(b)
table(b$STATUS_YR)

# Select PAs within each BIOME
for (eco_id in 1: length(BIOME_names)){
  # for each ecoregions
  cat("Biome", eco_id, "started...", "\n")
  eco_shp <- shapefile(paste0("./Output/Biome_scale/Forest_cover_change_per_biome/", eco_id , "_range.shp"))
  eco_shp$area_sqkm <- area(eco_shp) / 1000000
  pAREA <- NULL
  nPAs <- NULL
  # Expands the given geometry to include the area within the specified width with specific styling options.
  eco_shp_buffer <- buffer(eco_shp, width=2)
  b_crop <- crop(b, eco_shp_buffer)
  shapefile(b_crop[,c("IUCN_CAT","STATUS_YR")], paste0("./Output/Biome_scale/PAs_within_biome/PAs_within_biome_",eco_id,".shp"), overwrite =T)
  
  if(identical(b_crop, NULL)){
    pAREA <- rep(0, length(seq(1800,1992,1)))
    nPAs <- rep(0, length(seq(1800,1992,1)))
  }else{
    plot(eco_shp_buffer, main= Ecoregions_forested$ECO_NAME[eco_id])
    plot(eco_shp, add=T , col ="blue")
    # plot(b_crop, add =T ,col="red")
    for ( year in seq(1800,1992,1)) {
      cat("year", year, "started...", "\n")
      eco_shp_sel <- b_crop[which(b_crop$STATUS_YR == year),c("IUCN_CAT")]
      if(length(eco_shp_sel) == 0){
        pAREA <- rbind(pAREA, 0) 
        nPAs <-  rbind(nPAs, 0)
      }else{
        eco_int <- intersect(eco_shp, eco_shp_sel)
        if(identical(eco_int, NULL)){
          pAREA <- rbind(pAREA, 0)
          nPAs <-  rbind(nPAs, 0)
        }else{
          #plot(eco_shp, main = paste0("Estiblished PAs before year ", year))
          #mtext(paste0("Numbers of PAs: ", length(eco_int)))
          #mtext(side = 1, Ecoregions_forested$ECO_NAME[eco_id])
          #plot(eco_int, add =T, col = "red", )
          # Calculate square km of its range map
          eco_int$area_sqkm <- area(eco_int) / 1000000
          pAREA <- rbind(pAREA, round((sum(eco_int$area_sqkm) / eco_shp$area_sqkm) * 100, 3))
          nPAs <-  rbind(nPAs, length(eco_shp_sel))
        }
      }
    }
  }
  res <- as.data.frame(cbind(pAREA, nPAs, seq(1800,1992,1)))
  names(res) <- c("Percent_areas", "number_of_PAs", "years")
  saveRDS(res, paste0("./Output/PAs_within_biome/Percent_areas_of_PAs_within_Biome", eco_id, ".rds"))
}

# Calculating basic statistics of PA data
PAs <- list()
for (i in 1:7){
  PAs[[i]] <- readRDS(paste0("./Output/PAs_within_biome/","Percent_areas_of_PAs_within_Biome", i,".rds"))
  PAs[[i]] <- as.data.frame(PAs[[i]])
  tail(PAs[[i]])
  PAs_res_1 <- NULL
  PAs_res_2 <- NULL
  for (j in 1:193) {
    PAs_res_1 <- rbind(PAs_res_1, sum(PAs[[i]][1:j,1]))
    PAs_res_2 <- rbind(PAs_res_2, sum(PAs[[i]][1:j,2]))
  }
  PA_res <- as.data.frame(cbind(PAs_res_1, PAs_res_2,PAs[[i]][,3]))
  names(PA_res) <- c("Percent_areas","number_of_PAs", "years")
  tail(PA_res)
  saveRDS(PA_res,paste0("./Output/PAs_within_biome/Percent_areas_of_PAs_within_BIOME_",i, ".rds"))
}

# ------------------------------------------------------------------------------ ----
# (2)  BIOME Analysis and help functions                                         ----

# Set the work directory
setwd("D:/OneDrive/SCI_Publications/Extinction_Debt_In_draft/DATA AVAILABILITY")                 

# Load the necessary packages
library(raster)
library(sf)

# There are 867 terrestrial ecoregions, classified into 14 different biomes (WWF)
# For the purpose of this study, we select only seven forested biomes.
# The type of Forest depends mainly on location - that is, distance from equator and altitude - and climate
# From high latitudes to the equator:
# Boreal Forests/Taiga
# Temperate Conifer Forests
# Temperate Broadleaf & Mixed Forests
# Mediterranean Forests, Woodlands & Scrub
# Tropical & Subtropical Coniferous Forests
# Tropical & Subtropical Dry Broadleaf Forests
# Tropical & Subtropical Moist Broadleaf Forests

Biomes <- st_read("./6_Other_necessary_input_data/Forested_Biome_dissolve.shp")  #14458
class(Biomes) 
Biomes <- as(Biomes, "Spatial")  
class(Biomes) 

# Get the biome name
BIOME_names <- as.character(unique(Biomes$BIOME_NAME))

# plot
spplot(Biomes, "OBJECTID")

# A 0.5 arc-degree Global Grid system
sp.r <- st_read("./6_Other_necessary_input_data/global_grids_0.5.shp")
sp.r <- as(sp.r, "Spatial")  

# Load vertebrates species richness raster, e.g. 0.5 degree
Animal_SR <- stack(list.files("./1_Vertebrates_species_richness_grouped_by_IUCN_category/", pattern = "0.5.tif$", full.names = T))

# Load the forest cover raster, here only Ensemble appoarch
Forest_Cover <- stack(list.files("./2_Ensemble_weighted_global_forest_cover/", pattern = "Ensemble_M5", full.names = T))

# Read in ArcGIS processed data, the purpose of this process is to do intersect analysis of PA RANGE MAP
# BIOME range map and global 0.5 degree grid system. The output of this process in stored in following folder 
# Because this process would take longer time in R...
pA_NPA <- read.csv("./6_Other_necessary_input_data/PA_NonPA_Cells.csv", header = T)

# ID: id of the 0.5 degree grid system
# Bx: Biome x
# PA: cells within each biome that have been covered by at least one protected area
# NPA: cells within each biome that never been covered by any protected area
summary(pA_NPA)

# New function for cropping the forest cover raster data
fc_ras_crop <- function(Biomes = Biomes, Forest_Cover = Forest_Cover, biome_name = 1, ForID = 2){
  
  t1 <- Sys.time() 
  yr <- seq(1500,1992,1)
  cat("##", "Biome", biome_name, BIOME_names[biome_name], "- forest cover for yr", yr[ForID],"processing...", "\n")
  
  if(!file.exists("./Output")){dir.create("./Output")}
  if(!file.exists("./Output/Biome_scale")){dir.create("./Output/Biome_scale")}
  if(!file.exists("./Output/Biome_scale/Forest_cover_clipped")){dir.create("./Output/Biome_scale/Forest_cover_clipped")}
  Biome_ind <- Biomes[which(Biomes$BIOME_NAME == BIOME_names[biome_name]),]
  
  # for each biome range cells, we use ArcGIS for selecting the correponding cells covered by any PAs
  biome_range_cells <- shapefile(paste0("./6_Other_necessary_input_data/",biome_name,"_range_cells.shp"))
  length(biome_range_cells) # 12242
  
  # selected PAs aND Non-PAs cells
  ee_PA <- pA_NPA[which(pA_NPA[,paste0("B",biome_name,"_PA")]==1),"ID"] 
  length(ee_PA)
  ee_NPA <- pA_NPA[which(pA_NPA[,paste0("B",biome_name,"_NPA")]==1),"ID"] 
  length(ee_NPA)
  if(length(ee_PA) + length(ee_NPA) == length(biome_range_cells)){print("The ArcGIS processed data are correct")} else{print("Warning!!! ArcGIS processed data")}
  
  #(1) 2942 cells
  ee_final_PA <- sp.r[ee_PA,]
  ee_final_PA <- rasterize(ee_final_PA, raster(extent(sp.r), res=0.5, crs = '+proj=longlat +datum=WGS84'))
  ee_final_PA <- ee_final_PA > 0
  plot(ee_final_PA, main = "Cells covered by any of the protected areas from 1800 CE to 1992 CE", col = "red")
  mtext(paste0(sum(complete.cases(values(ee_final_PA)))," cells within biome: ", BIOME_names[biome_name]))
  
  #(2) 9300 cells
  ee_final_NPA <- sp.r[ee_NPA,]
  ee_final_NPA <- rasterize(ee_final_NPA, raster(extent(sp.r), res=0.5, crs = '+proj=longlat +datum=WGS84'))
  ee_final_NPA <- ee_final_NPA > 0
  plot(ee_final_NPA, main = "Cells never covered by any of the protected areas from 1800 CE to 1992 CE", col = "blue")
  mtext(paste0(sum(complete.cases(values(ee_final_NPA)))," cells within biome: ", BIOME_names[biome_name]))
  
  #(3) 12242 cells
  ee_final_ALL <- sp.r[c(ee_PA,ee_NPA),]
  ee_final_ALL <- rasterize(ee_final_ALL, raster(extent(sp.r), res=0.5, crs = '+proj=longlat +datum=WGS84'))
  ee_final_ALL <- ee_final_ALL > 0
  plot(ee_final_ALL, main = "All Cells within the following biome", col = "yellow")
  mtext(mtext(paste0(sum(complete.cases(values(ee_final_ALL)))," cells within biome: ", BIOME_names[biome_name])))
  
  # Check the numbers, should be TRUE
  if(sum(complete.cases(values(ee_final_PA))) + sum(complete.cases(values(ee_final_NPA))) ==  
     sum(complete.cases(values(ee_final_ALL)))){print("The rasterize process is correct")} else{print("Warning!!! Rasterize process")}
  
  
  # Forest Cover
  # using crop() and mask() functions will take ca. 90 mins for single one biome
  # using the following equations will take less than 20 mins
  # raster::beginCluster(n=myCores) # speeds up raster::crop and other raster functions
  
  Forest_Cover_PA   <- ee_final_PA  * Forest_Cover[[ForID]]
  Forest_Cover_NPA  <- ee_final_NPA * Forest_Cover[[ForID]]
  Forest_Cover_ALL  <- ee_final_ALL * Forest_Cover[[ForID]]
  
  names(Forest_Cover_PA)   <- names(Forest_Cover[[ForID]])
  names(Forest_Cover_NPA ) <- names(Forest_Cover[[ForID]])
  names(Forest_Cover_ALL ) <- names(Forest_Cover[[ForID]])
  
  plot(Forest_Cover_PA, main = names(Forest_Cover_PA))
  mtext(paste0(sum(complete.cases(values(Forest_Cover_PA)))," cells within biome: ", BIOME_names[biome_name]))
  
  plot(Forest_Cover_NPA, main = names(Forest_Cover_NPA))
  mtext(paste0(sum(complete.cases(values(Forest_Cover_NPA)))," cells within biome: ", BIOME_names[biome_name]))
  
  plot(Forest_Cover_ALL, main = names(Forest_Cover_ALL))
  mtext(paste0(sum(complete.cases(values(Forest_Cover_ALL)))," cells within biome: ", BIOME_names[biome_name]))
  
  # Check the numbers, should be TRUE
  if(sum(complete.cases(values(Forest_Cover_PA))) + sum(complete.cases(values(Forest_Cover_NPA))) ==  
     sum(complete.cases(values(Forest_Cover_ALL)))){print("The raster calculation is correct")} else{print("Warning!!!")}
  
  
  writeRaster(Forest_Cover_PA, paste0("./Output/Biome_scale/Forest_cover_clipped/",names(Forest_Cover_PA),"_Biome",biome_name,"_PA_forest_cover.tif"), overwrite =T)
  writeRaster(Forest_Cover_NPA, paste0("./Output/Biome_scale/Forest_cover_clipped/",names(Forest_Cover_NPA),"_Biome",biome_name,"_NPA_forest_cover.tif"), overwrite =T)
  writeRaster(Forest_Cover_ALL, paste0("./Output/Biome_scale/Forest_cover_clipped/",names(Forest_Cover_ALL),"_Biome",biome_name,"_ALL_forest_cover.tif"), overwrite =T)
  Forest_Cover_PA_Value <- data.frame(values(Forest_Cover_PA))
  Forest_Cover_NPA_Value <- data.frame(values(Forest_Cover_NPA))
  Forest_Cover_ALL_Value <- data.frame(values(Forest_Cover_ALL))
  names(Forest_Cover_PA_Value) <- paste0(names(Forest_Cover_PA),"_Biome",biome_name,"_PA_forest_cover")
  names(Forest_Cover_NPA_Value) <- paste0(names(Forest_Cover_NPA),"_Biome",biome_name,"_NPA_forest_cover")
  names(Forest_Cover_ALL_Value) <- paste0(names(Forest_Cover_ALL),"_Biome",biome_name,"_ALL_forest_cover")
  head(Forest_Cover_PA_Value)
  head(Forest_Cover_NPA_Value)
  head(Forest_Cover_ALL_Value)
  saveRDS(Forest_Cover_PA_Value, paste0("./Output/Biome_scale/Forest_cover_clipped/",names(Forest_Cover_PA),"_Biome",biome_name,"_PA_forest_cover.rds"))
  saveRDS(Forest_Cover_NPA_Value, paste0("./Output/Biome_scale/Forest_cover_clipped/",names(Forest_Cover_NPA),"_Biome",biome_name,"_NPA_forest_cover.rds"))
  saveRDS(Forest_Cover_ALL_Value, paste0("./Output/Biome_scale/Forest_cover_clipped/",names(Forest_Cover_ALL),"_Biome",biome_name,"_ALL_forest_cover.rds"))
  
  t2 <- Sys.time()
  print(t2-t1)
  cat("##", "Biome", biome_name, BIOME_names[biome_name], "- forest cover for yr", yr[ForID],"done...", "\n")
}

# for example, biome 1 - Boreal Forests/Taiga, year 1992
fc_ras_crop(Biomes = Biomes, Forest_Cover = Forest_Cover, biome_name = 1, ForID = 493)

# start parallel computation 
library(doParallel)
cl <- makeCluster(96)
registerDoParallel(cl)

biome_name = 1:7; ForID = 1:493
var_comb <- expand.grid(biome_name, ForID)
names(var_comb) <- c("biome_name", "ForID")

foreach(i = 1:length(var_comb$biome_name), .packages=c("raster")) %dopar% {
  fc_ras_crop(Biome = Biome,  Forest_Cover = Forest_Cover, biome_name = var_comb$biome_name[i], ForID = var_comb$ForID[i])
}

# function for processing animal richness rasters
ar_ras_crop <- function(Biomes = Biomes, Animal_SR = Animal_SR, biome_name = 1, AnsrID = 1:18){
  
  cat("##", "Biome", biome_name, BIOME_names[biome_name], "processing...", "\n")
  
  t1 <- Sys.time()   
  if(!file.exists("./Output/Biome_scale/Animal_richness_clipped")) dir.create("./Output/Biome_scale/Animal_richness_clipped")
  
  Biome_ind <- Biomes[which(Biomes$BIOME_NAME == BIOME_names[biome_name]),]
  
  # for each biome range cells, we use ArcGIS for selecting the correponding cells covered by any PAs
  biome_range_cells <- shapefile(paste0("./6_Other_necessary_input_data/",biome_name,"_range_cells.shp"))
  length(biome_range_cells) # 12242
  
  # selected PAs aND Non-PAs cells
  ee_PA <- pA_NPA[which(pA_NPA[,paste0("B",biome_name,"_PA")]==1),"ID"] 
  length(ee_PA)
  ee_NPA <- pA_NPA[which(pA_NPA[,paste0("B",biome_name,"_NPA")]==1),"ID"] 
  length(ee_NPA)
  if(length(ee_PA) + length(ee_NPA) == length(biome_range_cells)){print("The ArcGIS processed data are correct")} else{print("Warning!!! ArcGIS processed data")}
  
  #(1) 2942 cells
  ee_final_PA <- sp.r[ee_PA,]
  ee_final_PA <- rasterize(ee_final_PA, raster(extent(sp.r), res=0.5, crs = '+proj=longlat +datum=WGS84'))
  ee_final_PA <- ee_final_PA > 0
  plot(ee_final_PA, main = "Cells covered by any of the protected areas from 1800 CE to 1992 CE", col = "red")
  mtext(paste0(sum(complete.cases(values(ee_final_PA)))," cells within biome: ", BIOME_names[biome_name]))
  
  #(2) 9300 cells
  ee_final_NPA <- sp.r[ee_NPA,]
  ee_final_NPA <- rasterize(ee_final_NPA, raster(extent(sp.r), res=0.5, crs = '+proj=longlat +datum=WGS84'))
  ee_final_NPA <- ee_final_NPA > 0
  plot(ee_final_NPA, main = "Cells never covered by any of the protected areas from 1800 CE to 1992 CE", col = "blue")
  mtext(paste0(sum(complete.cases(values(ee_final_NPA)))," cells within biome: ", BIOME_names[biome_name]))
  
  
  #(3) 12242 cells
  ee_final_ALL <- sp.r[c(ee_PA,ee_NPA),]
  ee_final_ALL <- rasterize(ee_final_ALL, raster(extent(sp.r), res=0.5, crs = '+proj=longlat +datum=WGS84'))
  ee_final_ALL <- ee_final_ALL > 0
  plot(ee_final_ALL, main = "All Cells within the following biome", col = "yellow")
  mtext(mtext(paste0(sum(complete.cases(values(ee_final_ALL)))," cells within biome: ", BIOME_names[biome_name])))
  
  # Check the numbers, should be TRUE
  if(sum(complete.cases(values(ee_final_PA))) + sum(complete.cases(values(ee_final_NPA))) ==  
     sum(complete.cases(values(ee_final_ALL)))){print("The rasterize process is correct")} else{print("Warning!!! Rasterize process")}
  
  # Vertebrates richness
  # using crop() and mask() functions will take ca. 90 mins for single one biome
  # using the following equations will take less than 20 mins
  # raster::beginCluster(n=myCores) # speeds up raster::crop and other raster functions
  
  Animal_SR_PA   <- ee_final_PA  * Animal_SR[[AnsrID]]
  Animal_SR_NPA  <- ee_final_NPA * Animal_SR[[AnsrID]]
  Animal_SR_ALL  <- ee_final_ALL * Animal_SR[[AnsrID]]
  
  names(Animal_SR_PA)   <- names(Animal_SR[[AnsrID]])
  names(Animal_SR_NPA)  <- names(Animal_SR[[AnsrID]])
  names(Animal_SR_ALL)  <- names(Animal_SR[[AnsrID]])
  
  plot(Animal_SR_PA, main = names(Animal_SR[[AnsrID]]))
  # mtext(paste0(sum(complete.cases(values(Animal_SR_PA)))," cells within biome: ", BIOME_names[biome_name]))
  
  plot(Animal_SR_NPA, main = names(Animal_SR[[AnsrID]]))
  # mtext(paste0(sum(complete.cases(values(Animal_SR_NPA)))," cells within biome: ", BIOME_names[biome_name]))
  
  plot(Animal_SR_ALL, main = names(Animal_SR[[AnsrID]]))
  # mtext(paste0(sum(complete.cases(values(Animal_SR_ALL)))," cells within biome: ", BIOME_names[biome_name]))
 
  # Check the numbers, should be TRUE
  if(sum(complete.cases(values(Animal_SR_PA))) + sum(complete.cases(values(Animal_SR_NPA))) ==  
     sum(complete.cases(values(Animal_SR_ALL)))){print("The raster calculation is correct")} else{print("Warning!!!")}
  
  for (i in  AnsrID ) {
    writeRaster(Animal_SR_PA[[i]], paste0("./Output/Biome_scale/Animal_richness_clipped/",names(Animal_SR_PA[[i]]),"_Biome",biome_name,"_PA_SR.tif"), overwrite =T)
    writeRaster(Animal_SR_NPA[[i]], paste0("./Output/Biome_scale/Animal_richness_clipped/",names(Animal_SR_NPA[[i]]),"_Biome",biome_name,"_NPA_SR.tif"), overwrite =T)
    writeRaster(Animal_SR_ALL[[i]], paste0("./Output/Biome_scale/Animal_richness_clipped/",names(Animal_SR_ALL[[i]]),"_Biome",biome_name,"_ALL_SR.tif"), overwrite =T)
    Animal_SR_PA_Value <- data.frame(values(Animal_SR_PA[[i]]))
    Animal_SR_NPA_Value <- data.frame(values(Animal_SR_NPA[[i]]))
    Animal_SR_ALL_Value <- data.frame(values(Animal_SR_ALL[[i]]))
    names(Animal_SR_PA_Value) <- paste0(names(Animal_SR_PA[[i]]),"_Biome",biome_name,"_PA_SR")
    names(Animal_SR_NPA_Value) <- paste0(names(Animal_SR_NPA[[i]]),"_Biome",biome_name,"_NPA_SR")
    names(Animal_SR_ALL_Value) <- paste0(names(Animal_SR_ALL[[i]]),"_Biome",biome_name,"_ALL_SR")
    head(Animal_SR_PA_Value)
    head(Animal_SR_NPA_Value)
    head(Animal_SR_ALL_Value)
    saveRDS(Animal_SR_PA_Value, paste0("./Output/Biome_scale/Animal_richness_clipped/",names(Animal_SR_PA[[i]]),"_Biome",biome_name,"_PA_SR.rds"))
    saveRDS(Animal_SR_NPA_Value, paste0("./Output/Biome_scale/Animal_richness_clipped/",names(Animal_SR_NPA[[i]]),"_Biome",biome_name,"_NPA_SR.rds"))
    saveRDS(Animal_SR_ALL_Value, paste0("./Output/Biome_scale/Animal_richness_clipped/",names(Animal_SR_ALL[[i]]),"_Biome",biome_name,"_ALL_SR.rds"))
  }
  
  t2 <- Sys.time()
  print(t2-t1)
  cat("##", "Biome", biome_name, BIOME_names[biome_name], "done...", "\n")
  cat("##", "All amphibian, reptile and mammal richness layers in this biome are completed", "\n")
}

# for example, biome 1 - Boreal Forests/Taiga, year 1992
ar_ras_crop(Biomes = Biomes, Animal_SR = Animal_SR, biome_name = 1, AnsrID = 1:18)

# loop
for (biome_name in 1:7) {ar_ras_crop(Biomes = Biomes, Animal_SR = Animal_SR, biome_name = biome_name, AnsrID = 1:18)}

# raw data function
Bio_stat <- function(biome_name = 1, Methods = "Ensemble"){
  
  t1 <- Sys.time()
  
  # Load the clipped animals richness values
  Animal_SR_PA_Value  <- list.files("./Output/Biome_scale/Animal_richness_clipped/",pattern = paste0("_Biome",biome_name,"_PA_SR.rds"),full.names =T)
  Animal_SR_NPA_Value  <- list.files("./Output/Biome_scale/Animal_richness_clipped/",pattern = paste0("_Biome",biome_name,"_NPA_SR.rds"),full.names =T)
  Animal_SR_ALL_Value  <- list.files("./Output/Biome_scale/Animal_richness_clipped/",pattern = paste0("_Biome",biome_name,"_ALL_SR.rds"),full.names =T)
  
  # Load the clipped forest cover values
  Forest_Cover_PA_Value <- list.files("./Output/Biome_scale/Forest_cover_clipped/",pattern = paste0("_Biome",biome_name,"_PA_forest_cover.rds"), full.names = T)
  Forest_Cover_NPA_Value <- list.files("./Output/Biome_scale/Forest_cover_clipped/",pattern = paste0("_Biome",biome_name,"_NPA_forest_cover.rds"), full.names = T)
  Forest_Cover_ALL_Value <- list.files("./Output/Biome_scale/Forest_cover_clipped/",pattern = paste0("_Biome",biome_name,"_ALL_forest_cover.rds"), full.names = T)
  
  res_all_PA <-cbind(Reduce(cbind,lapply(Animal_SR_PA_Value , readRDS)),
                     Reduce(cbind,lapply(Forest_Cover_PA_Value , readRDS)))
  
  res_all_NPA <-cbind(Reduce(cbind,lapply(Animal_SR_NPA_Value , readRDS)),
                      Reduce(cbind,lapply(Forest_Cover_NPA_Value , readRDS)))
  
  res_all_ALL <-cbind(Reduce(cbind,lapply(Animal_SR_ALL_Value , readRDS)),
                      Reduce(cbind,lapply(Forest_Cover_ALL_Value , readRDS)))
  
  # Extract values from the above two rasters
  x <- area(Animal_SR[[1]])
  cord <- coordinates(x)
  plot(x)
  x <- values(x)
  length(x)
  
  anaysis_table_PA  <- as.data.frame(cbind("id" = 1:length(x), "sqkm" = x, cord, res_all_PA))
  anaysis_table_NPA <- as.data.frame(cbind("id" = 1:length(x), "sqkm" = x, cord, res_all_NPA))
  anaysis_table_ALL <- as.data.frame(cbind("id" = 1:length(x), "sqkm" = x, cord, res_all_ALL))
  anaysis_table_PA  <- anaysis_table_PA[complete.cases(anaysis_table_PA),]
  anaysis_table_NPA <- anaysis_table_NPA[complete.cases(anaysis_table_NPA),]
  anaysis_table_ALL <- anaysis_table_ALL[complete.cases(anaysis_table_ALL),]
  

  save(anaysis_table_PA,anaysis_table_NPA, anaysis_table_ALL, file = paste0("./Output/Biome_scale/3Y_data/", biome_name, "_raw_data.RData"))
  
  # To load the data again
  # load(paste0("./Output/Biome_scale/3Y_data/", biome_name, "_raw_data.RData"))
  
  t2 <- Sys.time()
  print(t2-t1)
  
}

# loop
for (biome_name in 1:7){Bio_stat(biome_name, Methods = "Ensemble")}


# ------------------------------------------------------------------------------ ----  
# (3)  Figure 3                                                                  ----

# Figure 3B 
# Forest cover change per biome - PAs vs Non-PAs 

library(ggsci) 
library(ggplot2)
library(reshape2)

IUCN_level = "All"
Resolution = 0.5
Methods = "Ensemble"
BIOME = 1
Categroy = "Mammals"


data <- list()

for (BIOME in 1:7) {
  data[[BIOME]] <- readRDS(paste0("./6_Other_necessary_input_data/3Y_data/",IUCN_level,"_",Categroy,"_",Resolution,"_",
                                  Methods,"_BIOME_", BIOME ,".rds"))
  # data[[sel]]$biome <- sel
  data[[BIOME]]$BIOME <- BIOME 
}

data_new <- Reduce(rbind, data)
summary(data_new)

pl <- data_new[which(data_new$number_protected_area!=0),]

pldata<-pl[,c("years","forest_cover_pa","forest_cover_npa","BIOME")]
pll<-melt(pldata,id=c('years','BIOME'))
pll$variable <-as.factor(pll$variable)
pll$BIOME <- factor(pll$BIOME, levels = 1:7)
p <- ggplot(pll,aes(years, value, color = BIOME,linetype = variable))+
  labs(x = "Year", y = "Forest cover (%)") +
  geom_line() +
  scale_color_d3()+
  scale_y_continuous(limits = c(0,70,10), breaks = seq(0,70,10))+
  scale_x_continuous(breaks = seq(1800, 2000, 40)) +
  theme_test()+
  geom_vline(aes(xintercept = 1870), colour = "coral", linetype = "dashed") + 
  geom_vline(aes(xintercept = 1969), colour = "blue", linetype = "dashed") 
p

# Figure 3C-E Correlation coefficients - Mammals, Amphibians and Reptiles

IUCN_level = "All"
Resolution = 0.5
Methods = "Ensemble"
BIOME = 1:7
Categroy = c("Mammals","Amphibians","Reptiles")


for (Categroy in Categroy) {
  
  data <- list()
  
  for (sel in 1:7) {
    data[[sel]] <- readRDS(paste0("./6_Other_necessary_input_data/3Y_data/",IUCN_level,"_",Categroy,"_",Resolution,"_",
                                  Methods,"_BIOME_", BIOME[sel] ,".rds"))
    data[[sel]]$biome <- sel
    # data[[sel]]$Categroy <- Categroy[sel]
  }
  
  data_new <- Reduce(rbind, data)
  summary(data_new)
  
  pl <- data_new[which(data_new$number_protected_area>=0),]
  
  pldata<-pl[,c("years","biome_cor_all","biome")]
  pll<-melt(pldata,id=c('years','biome'))
  pll$variable <-as.factor(pll$variable)
  pll$biome <- factor(pll$biome, levels = 1:7)
  p <- ggplot(pll,aes(years, value, color = biome, linetype = variable))+
    labs(x = "Year", y = "Pearson's r") +
    geom_line() +
    scale_y_continuous(limits = c(-0.05,0.7,0.1), breaks = seq(0, 0.8, 0.1)) +
    scale_x_continuous(breaks = seq(1500, 2000, 50)) +
    scale_color_d3()+
    theme_test()+
    geom_vline(aes(xintercept = 1784), colour = "black", linetype = "dashed") +
    geom_vline(aes(xintercept = 1870), colour = "coral", linetype = "dashed") + 
    geom_vline(aes(xintercept = 1969), colour = "blue", linetype = "dashed") 
  
  p

}

# ------------------------------------------------------------------------------ ----  
# (4)  Figure 4                                                                  ----

# Correlation coefficients per biome - PAs vs Non-PAs 
library(ggsci) 
library(ggplot2)
library(reshape2)

IUCN_level = "All"
Resolution = 0.5
Methods = "Ensemble"
BIOME = 1
Categroy = c("Mammals","Amphibians","Reptiles")


for (BIOME in 1:7) {
  
  data <- list()
  
  for (sel in 1:3) {
    data[[sel]] <- readRDS(paste0("./6_Other_necessary_input_data/3Y_data/",IUCN_level,"_",Categroy[sel],"_",Resolution,"_",
                                  Methods,"_BIOME_", BIOME ,".rds"))
    data[[sel]]$Categroy <- Categroy[sel]
  }
  
  data_new <- Reduce(rbind, data)
  summary(data)
  
  pl <- data_new[which(data_new$number_protected_area!=0),]
  
  pldata<-pl[,c("years","biome_cor_pa","biome_cor_npa","Categroy")]
  pll<-melt(pldata,id=c('years','Categroy'))
  pll$variable <-as.factor(pll$variable)
  pll$Categroy <- factor(pll$Categroy, levels = c("Mammals",
                                                  "Amphibians",
                                                  "Reptiles"))
  p <- ggplot(pll,aes(years, value, color = Categroy,linetype = variable))+
    labs(x = "Year", y = "Pearson's r") +
    geom_line() +
    scale_y_continuous(breaks = seq(0, 0.8, 0.1)) +
    scale_x_continuous(limits = c(1800,2000,40),breaks = seq(1800, 2000, 40)) +
    scale_color_aaas()+
    theme_test()+
    geom_vline(aes(xintercept = 1870), colour = "coral", linetype = "dashed") + 
    geom_vline(aes(xintercept = 1969), colour = "blue", linetype = "dashed") 
  
  p
}

## forest cover change vs. protected areas

# Load the R package, where ggplot2 is used for plotting, 
# gtable is used for extracting image attributes, and grid is used for merging figures
library(ggplot2)
library(gtable)
library(grid)

# Define functions to combine ggplot2 plot results to build dual axes
# Citation: https://stackoverflow.com/questions/36754891/ggplot2-adding-secondary-y-axis-on-top-of-a-plot

y2_plot <- function(p1, p2) {
  p1 <- ggplotGrob(p1)
  p2 <- ggplotGrob(p2)
  
  # Get the location of the plot panel in p1.
  # These are used later when transformed elements of p2 are put back into p1
  pp <- c(subset(p1$layout, name == 'panel', se = t:r))
  
  # Overlap panel for second plot on that of the first plot
  p1 <- gtable_add_grob(p1, p2$grobs[[which(p2$layout$name == 'panel')]], pp$t, pp$l, pp$b, pp$l)
  
  # Then proceed as before:
  
  # ggplot contains many labels that are themselves complex grob; 
  # usually a text grob surrounded by margins.
  # When moving the grobs from, say, the left to the right of a plot,
  # Make sure the margins and the justifications are swapped around.
  # The function below does the swapping.
  # Taken from the cowplot package:
  # https://github.com/wilkelab/cowplot/blob/master/R/switch_axis.R 
  
  hinvert_title_grob <- function(grob){
    
    # Swap the widths
    widths <- grob$widths
    grob$widths[1] <- widths[3]
    grob$widths[3] <- widths[1]
    grob$vp[[1]]$layout$widths[1] <- widths[3]
    grob$vp[[1]]$layout$widths[3] <- widths[1]
    
    # Fix the justification
    grob$children[[1]]$hjust <- 1 - grob$children[[1]]$hjust 
    grob$children[[1]]$vjust <- 1 - grob$children[[1]]$vjust 
    grob$children[[1]]$x <- unit(1, 'npc') - grob$children[[1]]$x
    grob
  }
  
  # Get the y axis title from p2
  index <- which(p2$layout$name == 'ylab-l') # Which grob contains the y axis title?
  ylab <- p2$grobs[[index]]                # Extract that grob
  ylab <- hinvert_title_grob(ylab)         # Swap margins and fix justifications
  
  # Put the transformed label on the right side of p1
  p1 <- gtable_add_cols(p1, p2$widths[p2$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, ylab, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'ylab-r')
  
  # Get the y axis from p2 (axis line, tick marks, and tick mark labels)
  index <- which(p2$layout$name == 'axis-l')  # Which grob
  yaxis <- p2$grobs[[index]]                  # Extract the grob
  
  # yaxis is a complex of grobs containing the axis line, the tick marks, and the tick mark labels.
  # The relevant grobs are contained in axis$children:
  #   axis$children[[1]] contains the axis line;
  #   axis$children[[2]] contains the tick marks and tick mark labels.
  
  # First, move the axis line to the left
  yaxis$children[[1]]$x <- unit.c(unit(0, 'npc'), unit(0, 'npc'))
  
  # Second, swap tick marks and tick mark labels
  ticks <- yaxis$children[[2]]
  ticks$widths <- rev(ticks$widths)
  ticks$grobs <- rev(ticks$grobs)
  
  # Third, move the tick marks
  ticks$grobs[[1]]$x <- ticks$grobs[[1]]$x - unit(1, 'npc') + unit(3, 'pt')
  
  # Fourth, swap margins and fix justifications for the tick mark labels
  ticks$grobs[[2]] <- hinvert_title_grob(ticks$grobs[[2]])
  
  # Fifth, put ticks back into yaxis
  yaxis$children[[2]] <- ticks
  
  # Put the transformed yaxis on the right side of p1
  p1 <- gtable_add_cols(p1, p2$widths[p2$layout[index, ]$l], pp$r)
  p1 <- gtable_add_grob(p1, yaxis, pp$t, pp$r + 1, pp$b, pp$r + 1, clip = 'off', name = 'axis-r')
  grid.newpage()
  grid.draw(p1)
}


IUCN_level = "All"
Resolution = 0.5
Methods = "Ensemble"
BIOME = 1
Categroy = "Mammals"

for (BIOME in 1:7) {
  
  data <- list()
  
  for (sel in 1) {
    data[[sel]] <- readRDS(paste0("./6_Other_necessary_input_data/3Y_data/",IUCN_level,"_",Categroy[sel],"_",Resolution,"_",
                                  Methods,"_BIOME_", BIOME ,".rds"))
    # data[[sel]]$biome <- sel
    data[[sel]]$Categroy <- Categroy[sel]
  }
  
  data_new <- Reduce(rbind, data)
  summary(data_new)
  
  pl <- data_new[which(data_new$number_protected_area>=0),]
  
  pldata<-pl[,c("years","forest_cover_all","Categroy")]
  pll<-melt(pldata,id=c('years','Categroy'))
  pll$variable <-as.factor(pll$variable)
  pll$Categroy <- factor(pll$Categroy)
  
  psdata<-pl[,c("years","P_protected_area","Categroy")]
  psl<-melt(psdata,id=c('years','Categroy'))
  psl$variable <-as.factor(psl$variable)
  psl$Categroy <- factor(psl$Categroy)
  
  
  p1 <- ggplot(pll, aes(years, value, color = Categroy))+
    labs(x = "Year", y = "Forest cover (%)") +
    geom_line(color = "blue") +
    scale_x_continuous(breaks = c(1500, 1600, 1700, 1800, 1850, 1900,1950,2000)) +
    scale_y_continuous(limits = c(0,70,5), breaks = seq(0,70,5))+
    scale_color_d3()+
    theme_test()+
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, color = 'black'), 
          axis.text.y = element_text(color = 'blue'), axis.ticks.y = element_line(color = 'blue'), 
          axis.title.y = element_text(color = 'blue')) +
    geom_vline(aes(xintercept = 1784), colour = "black", linetype = "dashed") +
    geom_vline(aes(xintercept = 1870), colour = "coral", linetype = "dashed") + 
    geom_vline(aes(xintercept = 1969), colour = "blue", linetype = "dashed") 
  p1
  
  p2 <- ggplot(psl, aes(years, value, color = Categroy, linetype = variable)) +
    geom_line(color = 'red') +
    scale_x_continuous(breaks = c(1500, 1600, 1700, 1800, 1850, 1900,1950,2000)) +
    scale_y_continuous(limits  = c(0,10, 1))+
    theme(panel.grid = element_blank(), panel.background = element_rect(fill = NA, color = 'black'),
          axis.text.y = element_text(color = 'red'), axis.ticks.y = element_line(color = 'red'), 
          axis.title.y = element_text(color = 'red')) +
    labs(y = 'Protected areas cover (%)')
  
  p2

  y2_plot(p1, p2)
  
  
}


# ------------------------------------------------------------------------------ ---- 
# EOF