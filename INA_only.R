# library(devtools)
#devtools::install_github("GarrettLab/INA")
library(INA)
library(tidyverse)
library(dplyr)


## ----load data-------------------------------------------------------------
memory.limit(size=112000) 
data<-read_csv("Inputs/Simple_points_for_INA.csv")
adj<-readRDS("Inputs/farm2farm_probs.rds")

#dimnames(adj)<-NULL
class(data)
class(adj)

length(adj[adj[]>0])
#adj[adj < 1/10000]<-0

#Farm HB08197 was the original farm that was infested - it is now a pine block.
row.names(adj)[742]

#big matrix version
dimnames(adj)<-NULL #remove dimnames because they increase memory use
geocoords2<-matrix(c(data$X, data$Y), byrow=F, ncol=2)
prob_est2<-as.vector(data$Probability_Estab)
prob_est2<-1-(1-prob_est2)*20/100 #option for varying probability of establishment based on a 20% reduction in establishment +between the best and worst sites with an ecoclimatic value >5..
class(prob_est2)
initbio2<- rep(0, length(data$farm_id))
initinfo2<-rep(1, length(data$farm_id))
initbio2[742]<-1 # original infested farm in Waipawa
class(initbio2)
initbio2
## ----example Large Matrix, error=TRUE, echo=TRUE--------------------------------------
Sys.time()
CNGvLarge <-
  INAscene(
    nreals = 30,
    ntimesteps = 2022-1962,
    doplot = F,
    outputvol = "less",
    readgeocoords = T,
    geocoords = geocoords2,
    numnodes = NA,
    xrange = NA,
    yrange = NA,
    randgeo = F,
    readinitinfo = F,
    initinfo = NA,#all nodes have information
    initinfo.norp = 'num',
    initinfo.n = 1830,
    initinfo.p = NA,
    initinfo.dist = 'random',
    readinitbio = F,
    initbio = NA, #initbio2 if you want to start the invasion in particular places
    initbio.norp = "num",
    initbio.n = 2, #starts in 1 or 2 nodes
    initbio.p = NA,
    initbio.dist = "random",
    readseam = F,
    seam = NA,
    seamdist = 'random',
    seamrandp = 1,
    seampla = NA,
    seamplb = NA,
    readbpam = T,
    bpam = adj,
    bpamdist = F,
    bpamrandp = NA,
    bpampla = NA,
    bpamplb = NA,
    readprobadoptvec = F,
    probadoptvec = NA,
    probadoptmean = seq(from = 0, to = 1, by = 0.2),
    probadoptsd = 0.2,
    readprobestabvec = F,
    probestabvec = NA,
    probestabmean = 1,
    probestabsd = 0.05,
    maneffdir = 'decrease_estab',
    maneffmean = seq(from = 0, to = 1, by = 0.2),
    maneffsd = 0.2,
    usethreshman = F,
    maneffthresh = NA,
    sampeffort = NA
  )
Sys.time()

saveRDS(CNGvLarge, "Inputs/Large_model_1in10.rds")
 




