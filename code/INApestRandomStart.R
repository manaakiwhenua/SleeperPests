###########################################################################
###########################################################################
###Declares a function simulating pest spread between areas of suitable climate and land use
###Key inputs are: 
###1) matrix of dispersal probibilities between each pair of sites (i.e. nodes of the network). Matrix can be non-symmetrical (i.e. can have source and sink nodes)
###2) Envionmentally-determined establishment probability for each node. Set all values to 1 for no environmental limitation on establishment.
###3) Management parameters
### a) Annual detection probability
### b) Annual management adoption probability subsequent to detection
### c) Annual erradication probability
### d) Spread reduction when management adopted
###4) A vector of invasion risk probabilities
###In this version initial invasion points are assigned through randomisation weighted by invasion risk
###This may be based on proximity to areas of suitable habitat (e.g. farms) in neighbouring regions where the pest is present
###Or proximity to likey border incursion points like ports and airports
###Key outputs are:
###4-dimensional arrays of invasion and management status for each node in each year of each permuation of each initial invasion configuration
###2-dimensional arrays of invasion probability (i.e. proportion of permutations pest present) for each node in each year
###Line graphs summarising number of nodes infested and under management against time  
###########################################################################
###########################################################################

INApestRandomStart = function(
Nreals,                 #Number of random initial invasion configurations
Nperm,                  #Number of permutations per initial configuratrion
Nyears,                 #Simulation duration
DetectionProb,          #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
ManageProb,             #Annual Probability of node adopting management upon detection
AnnualErradicationProb, #Annual probability of erradication (must be between 0 and 1) when management adopted
SpreadReduction,        #Reduction in dispersal probability when management adopted. Must be betweenn 0 (no spread reduction) and 1 (complete prevention of spread)
LDDrate,                #Mean annual per-node long distance dispersal rate. Value reflects rate for fully invaded network and zero
                        #environmental limitation on establishment
InitBioP,		#Proportion of nodes infested at start of simulations
InvasionRisk,           #Vector of probabilities for weighting random assignment of initial invasion occurrences
EnvEstabProb,           #Environmentally-determined establishment probability
BPAM,                   #Biophysical adjaceny matrix - disperal probability between each pair of nodes
SEAM = NA, 		#Option to provide socioeconomic adjacency matrix for information spread
LDDmatrix=NA,           #Option to provide long distance dispersal matrix instead of distance-independent dispesal rate
			#e.g. could be weighted by law of human visitation or data on stock movements
geocoords,              #XY points for INAscene
OutputDir		#Directory for storing results	
)
{
###Declare array for tracking management adoption status 
###of individual nodes in each year of each realisation
###This is a measure of potential disruption to farm businesses
###or ongoing management burden (surveillance and removal)
###for publicly-owned lands
ManagingResults = array(dim = c(nrow(BPAM),Nyears,Nperm,Nreals))

###Declare array tracking detection status 
###of individual nodes in each year of each realisation
DetectedResults = ManagingResults

###Declare array tracking infestation status 
###of individual nodes in each year of each realisation
InvasionResults = ManagingResults


###Share long distance dispersal events evenly across source nodes
###If no long-distance dispersal matrix is provided
if(is.na(LDDmatrix) == T)
      {
      ###Convert annual LDDrate to probability of achieving >0 disperal events
      LDDprob = ppois(0,LDDrate/(nrow(BPAM)-1),lower.tail = F)
      
      ###Treat SDD and LDD as separate trials when combining to obtain a single
      ###dispersal probability matrix
      BPAM =1-(1-BPAM)*(1-LDDprob)  
      }
if(is.na(LDDmatrix) == F)
      {
      ###Convert annual LDDrate to probability of achieving >0 disperal events
      LDDprob = ppois(0,LDDmatrix,lower.tail = F)

      ###Treat SDD and LDD as separate trials when combining to obtain a single
      ###dispersal probability matrix
      BPAM =1-(1-BPAM)*(1-LDDmatrix)
      }

###Weight dispersal by environmental establishment probability
BPAM = sweep(BPAM,2,EnvEstabProb,`*`)

###Make sure biophysical adjaceny matrix diagonals = 1
diag(BPAM) = 1

###SEAM is needed to run INA but zero info transfer between nodes is assumed unless SEAM provided
if(is.na(SEAM) == T)
{
SEAM = BPAM
SEAM[,] = 0
}


for(real in 1:Nreals)
 {
 ###Declare initial infestation vector
 ###and assign infestations based on invasion risk
 ###If no invasion risk vector is given initial invasion is equiprobable for all nodes
 InitBio = rep(0,times = nrow(BPAM))
 if(length(InvasionRisk) == nrow(BPAM))
 	Infested = sample(1:nrow(BPAM),size = ceiling(nrow(BPAM)*InitBioP),prob = InvasionRisk)
 if(length(InvasionRisk) != nrow(BPAM))
 	Infested = sample(1:nrow(BPAM),size = ceiling(nrow(BPAM)*InitBioP))
 InitBio[Infested] = 1

###Probability of info at start of simulation depends on
###Presence of pest and detection probability
InitInfoProb = InitBio*DetectionProb

###Loop through each realisation
###Nodes with information (which have detected infestation)
###Vary across realisations 
for(perm in 1:Nperm)
{

###Randomly assign probability of mangement adoption upon detection of infestation
ManageProbSD = 0.01
NodeManageProb = rnorm(ManageProb,ManageProbSD,n = nrow(BPAM))

###Randomly assign spread reduction factor when management adopted
SpreadReductionSD = 0.00001
NodeSpreadReduction = rnorm(SpreadReduction,SpreadReductionSD,n = nrow(BPAM))

###Select nodes that have detected infestation 
InitInfo = vector(length = nrow(BPAM))
RandInitInfoProb = runif(n=nrow(BPAM),0,1)

InitInfo = ifelse(RandInitInfoProb>InitInfoProb,0,1) 
Invaded = InitBio 
HaveInfo = InitInfo
###Loop through years
###This allows adoption of management based on detection of infestation
###Non-infested nodes never manage, but nodes that have been infested may manage
###even post erradication
for(year in 1:Nyears)
  {
  cat("\r","Config ",real, " Realisation ", perm, " Year ", year, "...")
  RandManageProb = runif(0,1,n=nrow(BPAM))
  Managing = vector(length = nrow(BPAM))
 ###This assumes that individual farmers or managers may not adopt management every single year
 ###Allows for some inconsistency in implementation for whatever reason
 Managing = ifelse(RandManageProb>NodeManageProb,0,1)	
  ###Management is only applied subsequent to detection
  Managing = Managing*HaveInfo
  ###Use this vector to moderate spread probabilities in biophysical adjacency matrix
  ManagementEffect =  Managing*NodeSpreadReduction
  Detected = Invaded*HaveInfo
  CNGvLarge <-
  INAscene(
    nreals = 1,
    ntimesteps = 1,
    doplot = F,
    outputvol = "more",
    readgeocoords = T,
    geocoords = geocoords,
    numnodes = NA,
    xrange = NA,
    yrange = NA,
    randgeo = F,
    readinitinfo = T,
    initinfo = HaveInfo,##Input nodes with info
    initinfo.norp = NA,
    initinfo.n = NA,
    initinfo.p = NA,
    initinfo.dist = NA,
    readinitbio = T,
    initbio = Invaded, ##Input infested nodes
    initbio.norp = NA,
    initbio.n = NA,
    initbio.p = NA,
    initbio.dist = NA,
    readseam = T,
    seam = SEAM*Detected, ##Only allow information spread from detected extant infestations
    seamdist = NA,
    seamrandp = NA,
    seampla = NA,
    seamplb = NA,
    readbpam = T,
    bpam =  BPAM*(1-ManagementEffect), ##biophysical adjacency matrix moderated by spread reduction in managing nodes
    bpamdist = F,
    bpamrandp = NA,
    bpampla = NA,
    bpamplb = NA,
    readprobadoptvec = T,
    probadoptvec = Managing, ##Use actual adoption calculated externally instead of adoption probability
    probadoptmean = NA,
    probadoptsd = NA,
    readprobestabvec = F,
    probestabvec = NA,
    probestabmean = 1,
    probestabsd = 0.00000001,
    maneffdir = 'decrease_estab',
    maneffmean = AnnualErradicationProb, ##Set mean management efficacy - effectively annual erradication prob
    maneffsd = 0.00000001,
    usethreshman = F,
    maneffthresh = NA,
    sampeffort = NA
  )

  LargeOut = CNGvLarge$multdetails
  Estab = as.array(LargeOut[[1]][[1]][[1]]$estabvecL)
  ###Update infestation vector
  Invaded = ifelse(Estab[[1]]==FALSE,0,1)
  
  InfoOut = as.vector(LargeOut[[1]]$multout[[1]]$vect1cL[[2]])
  HaveInfo[HaveInfo == 0] = InfoOut[HaveInfo == 0]

  ###Record nodes adopting management
  ManagingResults[,year,perm,real] = Managing
  ###Record infested nodes
  InvasionResults[,year,perm,real] = Invaded
  
  ###Select new nodes where infestation detected
  NewInfoProb = Invaded*DetectionProbs[detprob]
  RandNewInfoProb = runif(n=nrow(BPAM),0,1)
  NewHaveInfo = vector(length= length(HaveInfo))
  NewHaveInfo =  ifelse(RandNewInfoProb>NewInfoProb,0,1) 

 ###Add newly detected infestation to info vector
 HaveInfo[HaveInfo==0] = NewHaveInfo[HaveInfo==0]  

 ###Record nodes where infestation detected
 DetectedResults[,year,perm,real] = HaveInfo*Invaded
 }
}
}
#############################################
###Save results for post-hoc stats
#############################################
###Detection probability, erradication probability and spread reduction recorded in filename
###Use standard filename format for ease of reading results to produce heat maps 
###and conduct post-hoc stats comparing managment settings/scenarios 
FileNameStem = paste0(OutputDir,"DetProb_",DetectionProb,"_ErradProb_",round(AnnualErradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")

saveRDS(ManagingResults, paste0(FileNameStem,"ManagingLargeOut.rds"))
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))

##########################################################
###Store annual farm-level invasion probs for heat maps
###and estimation of invasion threat to other regions
##########################################################

InvasionProb = matrix(ncol = Nyears, nrow = nrow(BPAM))
for(year in 1:Nyears)
{
YearData = InvasionResults[,year,,]
InvasionProb[,year] = rowSums(YearData)/(Nperm*Nreals)
}
saveRDS(InvasionProb, paste0(FileNameStem,"InvasionProb.rds"))



###########################################################
###Produce summary figs when processing completed
###########################################################

Title = paste0("Detection Prob. ",DetectionProb," Erradication Prob. ",round(AnnualErradicationProb,digits =3),
		";\nSpread Reduction ",SpreadReduction)

InvasionSummary = as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(InvasionSummary) = c("Config",  "Realisation",   "Year",  "NodesInfested")


for(real in 1:Nreals)
for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,perm,real]
NodesInfested = colSums(InvasionData)
Config = rep(real,times = Nyears)
Realisation = perm 
Year = 1:Nyears
Results = data.frame(Config,Realisation,Year,NodesInfested)
InvasionSummary = rbind(InvasionSummary,Results)
}

Filename = paste0(FileNameStem,"InvasionRaw.png")
png(Filename)
plot(InvasionSummary$Year,InvasionSummary$NodesInfested,ylim = c(0,max(InvasionSummary$NodesInfested)),pch = NA
, xlab = "Time since incursion detected (years)",
ylab = "Number of nodes infested",main = Title)
for(real in 1:Nreals)
for(perm in 1:Nperm)
{
Sub = InvasionSummary[InvasionSummary$Config == real,]
Sub = Sub[Sub$Realisation == perm,]
lines(Sub$Year,Sub$NodesInfested,col  = real)
}
dev.off()

Quantiles = as.data.frame(aggregate(InvasionSummary$NodesInfested, by = list(InvasionSummary$Year),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"InvasionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (years)",
ylab = "Number of nodes infested",main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()

ManagingSummary = as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(InvasionSummary) = c("Config",  "Realisation",   "Year",  "NodesManaging")
for(real in 1:Nreals)
for(perm in 1:Nperm)
{
ManagingData = ManagingResults[,,perm,real]
NodesManaging = colSums(ManagingData)
Config = rep(real,times = Nyears)
Realisation = perm 
Year = 1:Nyears
Results = data.frame(Config,Realisation,Year,NodesManaging)
ManagingSummary = rbind(ManagingSummary,Results)
}

Filename = paste0(FileNameStem,"ManagingRaw.png")
png(Filename)
plot(ManagingSummary$Year,ManagingSummary$NodesManaging,ylim = c(0,max(ManagingSummary$NodesManaging)),pch = NA
, xlab = "Time since incursion detected (years)",
ylab = "Number of nodes under management",main = Title)
for(real in 1:Nreals)
for(perm in 1:Nperm)
{
Sub = ManagingSummary[ManagingSummary$Config == real,]
Sub = Sub[Sub$Realisation == perm,]
lines(Sub$Year,Sub$NodesManaging,col  = real)
}
dev.off()

Filename = paste0(FileNameStem,"ManagingSummary.png")
png(Filename)
Quantiles = as.data.frame(aggregate(ManagingSummary$NodesManaging, by = list(ManagingSummary$Year),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (years)",
ylab = "Number of nodes under management",main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()

DetectedSummary = as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(DetectedSummary) = c("Config",  "Realisation",   "Year",  "NodesDetected")
for(real in 1:Nreals)
for(perm in 1:Nperm)
{
DetectedData = DetectedResults[,,perm,real]
NodesDetected = colSums(DetectedData)
Config = rep(real,times = Nyears)
Realisation = perm 
Year = 1:Nyears
Results = data.frame(Config,Realisation,Year,NodesDetected)
DetectedSummary = rbind(DetectedSummary,Results)
}

Filename = paste0(FileNameStem,"DetectedRaw.png")
png(Filename)
plot(DetectedSummary$Year,DetectedSummary$NodesDetected,ylim = c(0,max(DetectedSummary$NodesDetected)),pch = NA
, xlab = "Time since incursion detected (years)",
ylab = "Number of nodes pest detected",main = Title)
for(real in 1:Nreals)
for(perm in 1:Nperm)
{
Sub = DetectedSummary[DetectedSummary$Config == real,]
Sub = Sub[Sub$Realisation == perm,]
lines(Sub$Year,Sub$NodesDetected,col  = real)
}
dev.off()

Filename = paste0(FileNameStem,"DetectedSummary.png")
png(Filename)
Quantiles = as.data.frame(aggregate(DetectedSummary$NodesDetected, by = list(DetectedSummary$Year),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (years)",
ylab = "Number of nodes pest detected",main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


DetectedProportionSummary = as.data.frame(matrix(ncol = 4, nrow = 0))
colnames(DetectedProportionSummary) = c("Config",  "Realisation",   "Year",  "DetectedProportion")
for(real in 1:Nreals)
for(perm in 1:Nperm)
{
DetectedData = DetectedResults[,,perm,real]
InvasionData = InvasionResults[,,perm,real]
NodesDetected = colSums(DetectedData)
NodesInvaded = colSums(InvasionData)
DetectedProportion = NodesDetected/NodesInvaded
DetectedProportion[is.na(DetectedProportion)==T] = 1
Config = rep(real,times = Nyears)
Realisation = perm 
Year = 1:Nyears
Results = data.frame(Config,Realisation,Year,DetectedProportion)
DetectedProportionSummary = rbind(DetectedProportionSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedProportionRaw.png")
png(Filename)
plot(DetectedProportionSummary$Year,DetectedProportionSummary$DetectedProportion,ylim = c(0,max(DetectedProportionSummary$DetectedProportion)),pch = NA
, xlab = "Time since incursion detected (years)",
ylab = "Proportion of infested nodes detected",main = Title)
for(real in 1:Nreals)
for(perm in 1:Nperm)
{
Sub = DetectedProportionSummary[DetectedProportionSummary$Config == real,]
Sub = Sub[Sub$Realisation == perm,]
lines(Sub$Year,Sub$DetectedProportion,col  = real)
}
dev.off()

Filename = paste0(FileNameStem,"DetectedProportionSummary.png")
png(Filename)
Quantiles = as.data.frame(aggregate(DetectedProportionSummary$DetectedProportion, by = list(DetectedProportionSummary$Year),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (years)",
ylab = "Proportion of infested nodes detected",main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()

}


################################################################
################################################################
###End of function
################################################################
################################################################
