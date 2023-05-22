##############################################
###Identifies climatically suitable farms for each region
###Converts ecoclimatic index values to establishment probability
###Calculates distance between farms within each region
###Prints maps ecoclimatic index values and establishment probability for each region
#############################################

val.path <- r"[N:\\Projects\BaseData\RPackages\R-4.1.1]"
if (dir.exists(val.path))
{
  .libPaths(val.path)
  .libPaths()
} else
{
  cat(val.path, "\r\n")
  stop( "Cannot find the lib path specified")
}

library(sf)
library(rgdal)
library(rgeos)
library(tidyverse)

############################################
###Load agribase data into R session
############################################
Agribase<-sf::st_read(AgribaseLocation)
AgribaseWS84 = st_transform(Agribase, "+proj=longlat +datum=WGS84")
st_crs(AgribaseWS84)
sf::st_bbox(AgribaseWS84)

############################################
###Identify regions then loop through each region to identify suitable farms 
###and generate distance matrices between each pair of suitable farms
############################################
Regions = sort(unique(AgribaseWS84$region))
Regions = Regions[Regions != "CHAT"]
Regions = Regions[is.na(Regions) == F]
for(region in 1:length(Regions))
{
RegionFarms.wgs84<-AgribaseWS84[AgribaseWS84$region==Regions[region],]
RegionFarms.wgs84<-sf::st_make_valid(RegionFarms.wgs84)
RegionFarms.wgs84<-sf::st_transform(RegionFarms.wgs84, CRS("+init=epsg:4326"))
sf::st_crs(RegionFarms.wgs84)
bb<-sf::st_bbox(RegionFarms.wgs84)
paste(bb)

############################################
###Group multiple polygons per farm into single multipart polygon per farm 
############################################
RegionFarms.wgs84v2<-RegionFarms.wgs84%>% 
  group_by(farm_id) %>% 
  summarize(size_ha= sum(size_ha), bef_nos=sum(bef_nos), shp_nos=sum(shp_nos), postal_twn=first(postal_twn), locality=first(locality), postal_cod=first(postal_cod), farm_type=first(farm_type))

############################################
###Select only sheep and beef farms
############################################
RegionFarms.wgs84SheepBeef<-RegionFarms.wgs84v2 %>% 
  mutate(sheep_or_beef=ifelse(farm_type=="BEF"
                              |farm_type=="SNB"
                              |farm_type=="SHP"
                               #|bef_nos>1000
                              , 1,0)) %>% 
					  filter(sheep_or_beef==1)

############################################
###Define centroids for each farm such that the centroid is forced to be in the farm even if the shape is irregular.
############################################
RegionFarms.wgs84_centroids<-sf::st_point_on_surface(RegionFarms.wgs84SheepBeef)
RegionFarms.wgs84_centroids2<- as.data.frame(sf::st_coordinates(RegionFarms.wgs84_centroids$geometry))
RegionFarms.wgs84SheepBeef<-bind_cols(RegionFarms.wgs84SheepBeef, RegionFarms.wgs84_centroids2) 
##Error in longitude for at least 1 farm in wellington region
RegionFarms.wgs84SheepBeef = RegionFarms.wgs84SheepBeef[RegionFarms.wgs84SheepBeef$X>0,]
bb<- sf::st_bbox(RegionFarms.wgs84SheepBeef)

###Save sheep and beef farms to shape file 
sf::st_write(RegionFarms.wgs84SheepBeef, paste0("Inputs/",Regions[region],"SheepBeef.wgs84.shp"),append=FALSE)

######################################################
###Get the ecoclimatic index per farm per CLIMEX model
######################################################
climex<-sf::st_read(ClimexFile)
climex<-sf::st_transform(climex, CRS("+init=epsg:4326"))
st_crs(climex)
st_bbox(climex)
climex_crop<-sf::st_crop(climex,bb)


st_bbox(climex_crop)

############################################
###Calculate mean EI for each farm
###Select only farms with mean EI >=6
###Separate code needed for each climate scenario (cs)
###As variable name differs for each
############################################
if(cs == 1)
{
Region_intersects_EI_6<-aggregate(climex_crop["EI"], RegionFarms.wgs84SheepBeef, mean, simplify=T, join = st_intersects, do_union=F)
Region_intersects_EI_6$EI<-ifelse(Region_intersects_EI_6$EI<6, NA,Region_intersects_EI_6$EI) 
RegionFarms.wgs84SheepBeef$Avg_EI<-as.numeric(Region_intersects_EI_6$EI)
RegionFarms.wgs84SheepBeef = RegionFarms.wgs84SheepBeef[is.na(RegionFarms.wgs84SheepBeef$Avg_EI) == F,]
}


if(cs == 2)
{
Region_intersects_EI_6<-aggregate(climex_crop["EI_Niwa204"], RegionFarms.wgs84SheepBeef, mean, simplify=T, join = st_intersects, do_union=F)
Region_intersects_EI_6$EI_Niwa204<-ifelse(Region_intersects_EI_6$EI_Niwa204<6, NA,Region_intersects_EI_6$EI_Niwa204) 
RegionFarms.wgs84SheepBeef$Avg_EI<-as.numeric(Region_intersects_EI_6$EI_Niwa204)
RegionFarms.wgs84SheepBeef = RegionFarms.wgs84SheepBeef[is.na(RegionFarms.wgs84SheepBeef$Avg_EI) == F,]
}

if(cs == 3)
{
Region_intersects_EI_6<-aggregate(climex_crop["EI_Niwa209"], RegionFarms.wgs84SheepBeef, mean, simplify=T, join = st_intersects, do_union=F)
Region_intersects_EI_6$EI_Niwa209<-ifelse(Region_intersects_EI_6$EI_Niwa209<6, NA,Region_intersects_EI_6$EI_Niwa209) 
RegionFarms.wgs84SheepBeef$Avg_EI<-as.numeric(Region_intersects_EI_6$EI_Niwa209)
RegionFarms.wgs84SheepBeef = RegionFarms.wgs84SheepBeef[is.na(RegionFarms.wgs84SheepBeef$Avg_EI) == F,]
}



if(nrow(RegionFarms.wgs84SheepBeef)>0)
{
################################################################################
###Allow different types of curve linking EI score to establishment
###Default is logit curve forced through x = 30, y = 0.8
###Other option provided is "Split linear" with different gradient above and below EI = 30
##################################################################################

############################################
###Different logit curve shapes forced through point EI = 30, prob = 0.8 are possible
###Use one which gives rapid increase in establighment at low EI values
###This is somewhat arbitrary
###Alternative is to use climate matching index which varies between 0 and 1 as establishment probability 
############################################
###Slow rate of increase at low values
#Intercept = -10
#Slope = 0.3795

###Moderate rate of increase at low values
#Intercept = -8
#Slope = 0.31287 

###High rate of increase at low values
Intercept = -5
Slope = 0.21287
Farm_infestation_status_INA<-st_set_geometry((RegionFarms.wgs84SheepBeef %>% 
                                                mutate(Probability_Estab=plogis(Intercept+Avg_EI*Slope))), NULL)
RegionFarms.wgs84SheepBeef$Probability_Estab = plogis(Intercept+RegionFarms.wgs84SheepBeef$Avg_EI*Slope)



#############################################
###Split linear option
#############################################
Prob= vector(length = length(RegionFarms.wgs84SheepBeef$Avg_EI))
Prob[RegionFarms.wgs84SheepBeef$Avg_EI<=30] = RegionFarms.wgs84SheepBeef$Avg_EI[RegionFarms.wgs84SheepBeef$Avg_EI<=30]*(80/30)/100 
Prob[RegionFarms.wgs84SheepBeef$Avg_EI>30] = 0.8+(RegionFarms.wgs84SheepBeef$Avg_EI[RegionFarms.wgs84SheepBeef$Avg_EI>30]-30)*(20/70)/100 
Farm_infestation_status_INA<-st_set_geometry((RegionFarms.wgs84SheepBeef %>% 
                                                mutate(Probability_Estab=Prob)), NULL)
###Add split linear option to dataframe
RegionFarms.wgs84SheepBeef$Probability_Estab_SplitLinear = Prob

#############################################
####Write out shape file of at risk farms simple points for INA and establishment probability
#############################################
write_csv(Farm_infestation_status_INA, paste0(ClimateScenarios[cs],"/Inputs/",Regions[region],"_Simple_points_for_INA.csv"))
sf::st_write(RegionFarms.wgs84SheepBeef, paste0(ClimateScenarios[cs],"/Inputs/",Regions[region],"SheepBeefAtRisk",ClimateScenarios[cs], ".wgs84.shp"),append=FALSE)
#hist(RegionFarms.wgs84SheepBeef$Avg_EI, main="Number of farms by ecoclimatic score")
#dev.off()

#########################################
######PLOT FARMS
######Scale fill to EI or establishment probability values
#########################################
p = ggplot()+
  geom_sf(data=RegionFarms.wgs84SheepBeef, aes(fill=Avg_EI))+ scale_fill_viridis_b(breaks=c(5,10,20,30,40,60,100),"CLIMEX \necoclimatic index")
ggsave(paste0(ClimateScenarios[cs],"/",Regions[region],"_SheepBeef_EI_",ClimateScenarios[cs],".jpg"))
p = ggplot()+
  geom_sf(data=RegionFarms.wgs84SheepBeef, aes(fill=Probability_Estab))+ scale_fill_viridis_b(breaks=c(0.05,0.1,0.3,0.5,0.7,0.9,1),"Establishment probability")
ggsave(paste0(ClimateScenarios[cs],"/",Regions[region],"_SheepBeef_ProbEstab_",ClimateScenarios[cs],".jpg"))


########################################################
###Generate distance matrix between each pair of at-risk farms
###Calculates minimum Euclidean distance between farms
###Using farm boundaries. 
###Farms with multiple far-fling parcels may have more close neighbours
########################################################

###gDistance seems to require this other spatial format for rgeos package
For_distance <-sf::as_Spatial(RegionFarms.wgs84SheepBeef)
For_distance <-sp::spTransform(For_distance, CRS("+proj=utm +zone=60 +south +datum=WGS84 +units=m +no_defs"))

###distance matrix all farms
Farm_Dist_Mat_Region<- rgeos::gDistance(For_distance, byid = TRUE)

###Use Farm ID as row and column names
farm_ids<-RegionFarms.wgs84SheepBeef$farm_id
rownames(Farm_Dist_Mat_Region) <- farm_ids
colnames(Farm_Dist_Mat_Region) <- farm_ids
dist_mat_farm<-as.matrix(Farm_Dist_Mat_Region)

###############################################################
###Store distance matrix as .RDS file
###############################################################
saveRDS(dist_mat_farm, paste0(ClimateScenarios[cs],"/Inputs/",Regions[region],"_dist_mat_farm.rds"))
#write.csv(Farm_Dist_Mat_Region, "Inputs/dist_mat_farm.csv")


}
}
###############################################################
###End of region loop
###############################################################

