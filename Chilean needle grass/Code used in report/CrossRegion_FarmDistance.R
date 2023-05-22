##################################################
###Generates distance matrices between farms in different regions
###Loops through each pair of regions and stores results for each pair separately
###Calculates minimum Euclidean distance between farm boundaries
###so that farms with parcels in multiple regions may have close neighbours in different regions
##################################################

library(INA)
library(tidyverse)
library(dplyr)
library(sp)
library(actuar)
memory.limit(size=600000) 

ResultsDir= paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
for(region1 in 1:(length(Regions)-1))
{
Region1Farms = sf::st_read(paste0(ClimateScenarios[cs],"/Inputs/",Regions[region1],"SheepBeefAtRisk",ClimateScenarios[cs], ".wgs84.shp"))
Region1FarmsSpatial <-sf::as_Spatial(Region1Farms)
Region1FarmsSpatial <-sp::spTransform(Region1FarmsSpatial, CRS("+proj=utm +zone=60 +south +datum=WGS84 +units=m +no_defs"))
for(region2 in (region1+1):length(Regions))
{
#----load data-------------------------------------------------------------

Region2Farms = sf::st_read(paste0(ClimateScenarios[cs],"/Inputs/",Regions[region2],"SheepBeefAtRisk",ClimateScenarios[cs], ".wgs84.shp"))

Region2FarmsSpatial <-sf::as_Spatial(Region2Farms)
Region2FarmsSpatial <-sp::spTransform(Region2FarmsSpatial, CRS("+proj=utm +zone=60 +south +datum=WGS84 +units=m +no_defs"))

###Calculate cross regions farm distances
###With Region1 farms as columns
Region1_Region2FarmDist<- rgeos::gDistance(Region1FarmsSpatial,Region2FarmsSpatial, byid = TRUE)
colnames(Region1_Region2FarmDist) = Region1Farms$farm_id
row.names(Region1_Region2FarmDist) = Region2Farms$farm_id
dist_mat<-as.matrix(Region1_Region2FarmDist)
saveRDS(dist_mat, paste0(ResultsDir,Regions[region1],"_",Regions[region2],"_dist_mat_farm_LDD.rds"))
}
}


