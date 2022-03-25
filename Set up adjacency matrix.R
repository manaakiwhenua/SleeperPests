#This sets up an adjacency matrix for network analysis


library(tidyverse)
library(dplyr)
library(igraph)
library(ggraph)
library(actuar)


#First decide on the relationship between dispersal distance and probability of dispersal

k= 7 #shape parameter
b= 500 #scale parameter
thresh= 1/10 #value for spread between adjacent farms
d=seq(from=0, to=50000, by=1 )

#library(actuar)

prob2<- thresh*actuar::ppareto(d, shape=k,  scale=b, lower.tail=F,log.p = F)

for_plot<-tibble(distance=d, probability=prob2, shape_par=k, scale_par=b)

#check probability at a few distances
for_plot %>% 
  filter(distance==0|distance==1|distance==150|distance==500|distance==50000|probability<1/10000000)

for_plot %>% 
  filter(distance<500) %>% 
  ggplot(aes(x=distance, y=probability))+ geom_line()


#Get the distance matrix and change to an adjacency matrix showing the probabilities of dispersal between nodes..


#READ DISTANCE MATRIX
dist_mat_farm<-readRDS("Inputs/dist_mat_farm.rds")

#check the number of farms to farm distances that are in a certain range
length(dist_mat_farm[dist_mat_farm[ ] ==0 ])
length(dist_mat_farm[dist_mat_farm[ ]>=50000 ])

#can make the distances that are shown as a zero a very small number (zeros work with this pareto distribution but not for powerlaw or negative exponential distributions)

dist_mat_farm[dist_mat_farm==0]<-0.01
farm2farm_probs <- thresh*actuar::ppareto(dist_mat_farm, shape=k,  scale=b, lower.tail=F,log.p = F)


# Because of the long fat tail you need to decide if you want to remove some low probability links. One way is to choose a probability associated with a distance, or another is to remove based on a low arbitrary value of 1/10000 or similar. Here is a way to calculate the probability associated with a distance and remove values less than that. A key concern is the more edges in the network the slower the analysis. Also the network becomes difficult to visualize.

#determine probability for 20 km and remove all probabilities lower. Later I set up a threshold of 1/10000 as the lowest probability of dispersal that we are interested in.

p20<-thresh*actuar::ppareto(20000, shape=k,  scale=b, lower.tail=F,log.p = F)

#how many farm links are removed

length(farm2farm_probs[farm2farm_probs[ ]<=p20 ])
length(farm2farm_probs[farm2farm_probs[ ]==0 ])
#replace low probs with zero
farm2farm_probs[farm2farm_probs < p20] <- 0

length(farm2farm_probs[farm2farm_probs[ ]==0 ])
length(farm2farm_probs[farm2farm_probs[ ]>0 ])
#set diagonal to 1   
diag(farm2farm_probs)<-1

saveRDS(farm2farm_probs, "Inputs/farm2farm_probs.rds")
