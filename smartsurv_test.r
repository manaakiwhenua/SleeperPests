# set up for smart_surv output

library(tidyverse)
library(ggraph)
library(igraph)
library(INA)
library(nzcensr)
library(ggspatial)
library(sf)
library(maps)

data<-read_csv("Inputs/Simple_points_for_INA.csv")
adj<-readRDS("Inputs/farm2farm_probs.rds")


#truncate the kernel so the fat tail is smaller and fewer links are available
# length(adj[adj[]>0])
# adj[adj < 1/10000]<-0
# length(adj[adj[]>0])

#Surveillance graphs
dim(adj)

surv_out<-smartsurv(adjmat = adj, stoch = F, nrealz=1)

#surv_out_stoch<-smartsurv(adjmat = adj, stoch = T, nrealz=10) #stop function for low probability nodes causes a problem 
saveRDS(surv_out, "Inputs/surv_out.rds")



# Graph a network


adj<-readRDS("Inputs/farm2farm_probs.rds")
length(adj[adj[]>0])
adj[adj < 1/10000]<-0
length(adj[adj[]>0])
surv_out<-readRDS("Inputs/surv_out.rds")
net<-graph_from_adjacency_matrix(adj, mode="directed", weighted = TRUE, diag = F)
#too many edges for plotting remove 
length(E(net))
net<-net %>% delete_edges(sample(seq(1, length(E(net))), size=5000))
length(V(net))
length(E(net))
layout<-create_layout(net, layout = "kk")
length(data$X)
blah<-colMeans(surv_out$meanarr)
length(blah)
layout$x<-data$X
layout$y<-data$Y
layout$smart<-blah




hist(blah)

p<-ggraph(layout) + 
  geom_edge_link(aes(color=weight, width=weight), alpha=0.25, show.legend = F) + 
  geom_node_point(aes(color=smart, size=smart)  )+
  scale_color_viridis_c(breaks=c(200, 400, 600, 800, 1000, 1200))+
  scale_edge_width(range = c(0, 0.3), guide="none")+
  scale_size_continuous(range = c(0, 3), breaks=c(200, 400, 600, 800, 1000, 1200))+
  guides(size=guide_legend("Uninfested nodes \nat detection"), color=guide_legend("Uninfested nodes \nat detection"))

p

ggsave("Uninfested nodes at detection2.jpg", width=8, height=8)

#map https://remiller1450.github.io/s230s19/Intro_maps.html




