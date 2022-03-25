#Get farms of interest, and associate CNG infestations with them

#Code for working with Agribase and Climex outputs.

library(sf)
library(rgdal)
library(rgeos)
library(tidyverse)
library(readxl)
#already subsetted Agribase for farms in Hawkes Bay but could use same code to obtain only Hawkes Bay region from the larger Agribase database
Agribase<-sf::st_read("Inputs/HB_wgs84.shp")


unique(Agribase$region)
HB.wgs84<-Agribase [Agribase$region=="HBAY",]

HB.wgs84<-sf::st_make_valid(HB.wgs84)
HB.wgs84<-sf::st_transform(HB.wgs84, CRS("+init=epsg:4326"))
#check class and projection, plus bounding box
class(HB.wgs84)
sf::st_crs(HB.wgs84)
bb<-sf::st_bbox(HB.wgs84)
paste(bb)

#run code like this if you are doing analysis region by region

#HB.wgs84<-sf::st_write(HB.wgs84, "Inputs/HB_wgs84.shp")

#group multiple polygons per farm into single multipart polygon per farm 
HB.wgs84v2<-HB.wgs84%>% 
  group_by(farm_id) %>% 
  summarize(size_ha= sum(size_ha), bef_nos=sum(bef_nos), shp_nos=sum(shp_nos), postal_twn=first(postal_twn), locality=first(locality), postal_cod=first(postal_cod), farm_type=first(farm_type))


#check crs and bounding again
bb<- sf::st_bbox(HB.wgs84v2)
bb


#centroids for each farm such that the centroid is forced to be in the farm even if the shape is irregular.

HB.wgs84_centroids<-sf::st_point_on_surface(HB.wgs84v2)

HB.wgs84_centroids2<- as.data.frame(sf::st_coordinates(HB.wgs84_centroids$geometry))
HB.wgs84v2<-bind_cols(HB.wgs84v2, HB.wgs84_centroids2) 
bb<- sf::st_bbox(HB.wgs84v2)
bb
#get the ecoenvironmental index per farm per CLIMEX model

climex<-sf::st_read("Inputs/NZ_fine_scale_fishnet_join.shp")

climex<-sf::st_transform(climex, CRS("+init=epsg:4326"))
st_crs(climex)

climex_crop<-sf::st_crop(climex,bb)
#st_bbox(climex)


#identify mean EI values for farms map, anything over a 5 is suitable

HB_intersects_EI_6<-aggregate(climex_crop["Avg_EI"], HB.wgs84v2, mean, simplify=T, join = st_intersects, do_union=F)
HB_intersects_EI_6$Avg_EI<-ifelse(HB_intersects_EI_6$Avg_EI<6, NA,HB_intersects_EI_6$Avg_EI) 
#add data back to main shape file
HB.wgs84v2$Avg_EI<-as.numeric(HB_intersects_EI_6$Avg_EI)

#load CNG data
crscheck<-sf::st_crs(HB.wgs84v2)
CNGpoints<-read_xlsx("Inputs/HBRC_CNG_points.xlsx")%>% 
  filter(Longitude>176)

CNGpoints2<-sf::st_as_sf(CNGpoints, coords = c("Longitude", "Latitude"), crs = crscheck)



#Identify farms with CNG points in them

HB_risk2<-sf::st_intersects(HB.wgs84v2,CNGpoints2, prepared = FALSE)
HB_risk2_values<-summary(HB_risk2)
check<-HB_risk2_values[,1]
HB.wgs84v2$CNGpresent<-as.numeric(check)
class(HB.wgs84v2)

#HB.wgs84v2$CNGpresent  #uncomment to check 

#some points aren't on farms, count records if points are within 20 m, for example along roads
HB_risk2<-sf::st_is_within_distance(HB.wgs84v2, CNGpoints2, 20)
HB_risk2_values<-summary(HB_risk2)
check<-HB_risk2_values[,1]
HB.wgs84v2$CNGclose<-as.numeric(check)
#class(HB.wgs84v2)
#HB.wgs84v2$CNGclose #uncomment to check

#remove low EI records, keep sheep and beef, and farms with CNG. This reduces the number of nodes in the network and makes the 

HB.wgs84v2<-HB.wgs84v2 %>% 
  mutate(sheep_or_beef_or_CNG=ifelse(farm_type=="BEF"
                              |farm_type=="SNB"
                              |farm_type=="SHP"
                              |CNGclose>0
                              #|bef_nos>1000
                              , 1,0)) %>% 
  filter(sheep_or_beef_or_CNG==1)
HB.wgs84v2<-HB.wgs84v2 %>% filter(Avg_EI!="NA")##determine higher risk farm types selecting sheep and/or beef

hist(HB.wgs84v2$Avg_EI, main="Number of farms by ecoclimatic score")
#write out the CSV for location without geometry
#Also fit a slope between 0 and 40 to set establishment probability which is based on a slope of one between the lowest and highest values of eh AVerage EI score.
Farm_infestation_status_INA<-st_set_geometry((HB.wgs84v2 %>% 
                                                mutate(Probability_Estab=Avg_EI*2.5/100)), NULL)

write_csv(Farm_infestation_status_INA, "Inputs/Simple_points_for_INA.csv")

#MAP OCCURRENCE DATA TO FARMS INCLUDE THOSE ADJACENT TO FARMS

#load CNG data
crscheck<-sf::st_crs(HB.wgs84v2)
CNGpoints<-read_xlsx("Inputs/HBRC_CNG_points.xlsx")%>% 
  filter(Longitude>176)

CNGpoints2<-sf::st_as_sf(CNGpoints, coords = c("Longitude", "Latitude"), crs = crscheck)








## The next step is to make a distance matrix, using polygons. Centroids would be simpler but using them to estimate farm to farm distances does not address farms that share a boundary. The goal is to use a dispersal kernel that takes into account farm distances a the boundary.

#gDistance seems to require this other spatial format for rgeos package

For_distance <-sf::as_Spatial(HB.wgs84v2)
For_distance <-sp::spTransform(For_distance, CRS("+proj=utm +zone=60 +south +datum=WGS84 +units=m +no_defs"))

#summary(For_distance)

# #distance matrix all farms
Farm_Dist_Mat_HB<- rgeos::gDistance(For_distance, byid = TRUE)
#row.names(Farm_Dist_Mat_HB)

farm_ids<-HB.wgs84v2$farm_id

rownames(Farm_Dist_Mat_HB) <- farm_ids
colnames(Farm_Dist_Mat_HB) <- farm_ids
dist_mat_farm<-as.matrix(Farm_Dist_Mat_HB)
saveRDS(dist_mat_farm, "Inputs/dist_mat_farm.rds")
#write.csv(Farm_Dist_Mat_HB, "Inputs/dist_mat_farm.csv")

## PLOT FARMS and OCCURRENCE INFORMATION

#CNG points with Avg EI, echo=FALSE, fig.cap= "Chilean needle grass records (red points) overlayed on farms with the corresponding ecoclimatic index."}

ggplot()+
  geom_sf(data=HB.wgs84v2, aes(fill=Avg_EI))+ scale_fill_viridis_b(breaks=c(5,10,20,30,40),"CLIMEX \necoclimatic index")+
  geom_point(data=CNGpoints, aes(x=Longitude, y=Latitude), shape=18, size=1, colour="red")
