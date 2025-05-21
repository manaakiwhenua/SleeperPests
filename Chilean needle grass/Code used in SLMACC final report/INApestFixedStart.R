###########################################################################
###########################################################################
###Declares a function simulating pest spread between areas of suitable climate and land use
###Key inputs are: 
###1) matrix of dispersal probibilities between each pair of sites (i.e. nodes of the network). Matrix can be non-symmetrical (i.e. can have source and sink nodes) 
###2) Envionmentally-determined establishment probability for each node. Set all values to 1 for no environmental limitation on establishment.
###3) Annual per-farm long distance dispersal rate or matrix of annual LDD rates for each pair of farms 
###4) Management parameters
### a) Annual detection probability
### b) Annual management adoption probability subsequent to detection
### c) Annual Eradication probability
### d) Spread reduction when management adopted
###5) A binary vector of intially infested nodes
###In this version initial invasion points are fixed
###This will usually apply when some current or historic distribution data
###are available for the pest species in question within the area of interest
###Key outputs are:
###3-dimensional arrays of invasion and management status for each node in each year of each permuation
###2-dimensional arrays of invasion probability (i.e. proportion of permutations pest present) for each node in each year
###Line graphs summarising number of nodes infested, where infestations
###are detected, the proportion of infestations detected and nodes under management against  time 
###########################################################################
###########################################################################

####################################################
###Development notes - Combine fixed and random start in single function
###Re-reun all analyses with eradication spell error fixed
#####################################################

INApestFixedStart = function(
Nperm,                  #Number of permutations per paramteter combination
Nyears,                 #Simulation duration
DetectionProb,          #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
ManageProb,             #Annual Probability or vector of probabilities vector length nrow(BPAM)of node adopting management upon detection
AnnualEradicationProb, #Annual probability of Eradication (must be between 0 and 1) when management adopted. Can be single value or vector length nrow(BPAM)
SpreadReduction,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread)
LDDrate,                #Mean annual per-node distant-independent dispersal rate. Value reflects rate for fully invaded network and zero
                        #environmental limitation on establishment
InitialInvasion,        #Nodes infested at start of simulations
EnvEstabProb,           #Environmentally determined establishment probability
BPAM,                   #Biophysical adjaceny matrix - disperal probability between each pair of nodes
SEAM = NA,			#Option to provide socioeconomic adjacency matrix for information spread
LDDmatrix = NA,         #Option to provide long distance dispersal matrix instead of distance-independent dispesal rate
			      #e.g. could be weighted by law of human visitation or data on stock movements
geocoords,              #XY points for INAscene
OutputDir		      #Directory for storing results	
)
{
###Declare array tracking infestation status 
###of individual nodes in each year of each realisation
InvasionResults = array(dim = c(nrow(BPAM),Nyears,Nperm))

###Declare array tracking detection status 
###of individual nodes in each year of each realisation
DetectedResults = InvasionResults

###Declare array for tracking management adoption status 
###of individual nodes in each year of each realisation
###This is a measure of potential disruption to farm businesses
###or ongoing management burden (surveillance and removal)
###for publicly-owned lands
ManagingResults = InvasionResults



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

###Probability of info at start of simulation depends on
###Presence of pest and detection probability
InitInfoProb = InitialInvasion*DetectionProb

###Loop through each realisation
###Nodes with information (which have detected infestation)
###Vary across realisations 
for(perm in 1:Nperm)
{
###Randomly assign probability of mangement adoption upon detection of infestation
ManageProbSD = mean(ManageProb)/10
NodeManageProb = rnorm(ManageProb,ManageProbSD,n = nrow(BPAM))
NodeManageProb[NodeManageProb<0] = 0
NodeManageProb[NodeManageProb>1] = 1
###Randomly assign spread reduction factor when management adopted
SpreadReductionSD = (1-mean(SpreadReduction))/10
NodeSpreadReduction = rnorm(SpreadReduction,SpreadReductionSD,n = nrow(BPAM))
NodeSpreadReduction[NodeSpreadReduction<0] = 0
NodeSpreadReduction[NodeSpreadReduction>1] = 1


###Select nodes that have detected infestation 
InitInfo = vector(length = nrow(BPAM))
RandInitInfoProb = runif(n=nrow(BPAM),0,1)
InitInfo = ifelse(RandInitInfoProb>InitInfoProb,0,1) 
Invaded = InitialInvasion 
HaveInfo = InitInfo
###Loop through years
###This allows allocation of info to farmers or managers based on detection of infestation
###Non-infested nodes never manage, but nodes that have been infested may manage
###even post Eradication
for(year in 1:Nyears)
  {
  cat("\r", "Realisation ", perm, "Year ", year, "...")
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
    seam = SEAM*Detected, ##Only allow information spread from dtected extant infestations
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
    maneffmean = AnnualEradicationProb, ##Set mean management efficacy - effectively annual Eradication prob
    maneffsd = 0.00000001,
    usethreshman = F,
    maneffthresh = NA,
    sampeffort = NA
  )

  LargeOut = CNGvLarge$multdetails
  ###Update info vector for any info spread (if SEAM supplied)
  InfoOut = as.vector(LargeOut[[1]]$multout[[1]]$vect1cL[[2]])
  HaveInfo[HaveInfo == 0] = InfoOut[HaveInfo == 0]
  Estab = as.array(LargeOut[[1]][[1]][[1]]$estabvecL)
  ###Update infestation vector
  Invaded = ifelse(Estab[[1]]==FALSE,0,1)

  ###Record nodes adopting management
  ManagingResults[,year,perm] = Managing
  ###Record infested nodes
  InvasionResults[,year,perm] = Invaded
  ###Select new nodes where infestation detected
  NewInfoProb = Invaded*DetectionProb
  RandNewInfoProb = runif(n=nrow(BPAM),0,1)
  NewHaveInfo = vector(length= length(HaveInfo))
  NewHaveInfo =  ifelse(RandNewInfoProb>NewInfoProb,0,1) 
 ###Add newly detected infestation to info vector
 HaveInfo[HaveInfo==0] = NewHaveInfo[HaveInfo==0]  
 ###Record detection status
 DetectedResults[,year,perm] = HaveInfo*Invaded
 }
}
###########################################################
###Save results for post-hoc analyses
###########################################################
###Detection probability, Eradication probability and spread reduction recorded in filename
###Use standard format for ease of reading results to produce heat maps 
###and conduct post-hoc stats comparing managment settings/scenarios 

FileNameStem = paste0(OutputDir,"DetProb_",DetectionProb,"_erradProb_",round(AnnualEradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
saveRDS(ManagingResults, paste0(FileNameStem,"InfoLargeOut.rds"))
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))

##########################################################
###Store annual farm-level invasion probs for heat maps
###and estimation of invasion threat to other regions
##########################################################

InvasionProb = matrix(ncol = Nyears, nrow = nrow(BPAM))
for(year in 1:Nyears)
{
YearData = InvasionResults[,year,]
InvasionProb[,year] = rowSums(YearData)/Nperm
}
saveRDS(InvasionProb, paste0(FileNameStem,"InvasionProb.rds"))


###########################################################
###Produce summary figs when processing completed
###########################################################

Title = paste0("Detection Prob. ",DetectionProb," Eradication Prob. ",round(AnnualEradicationProb,digits =3),
		"\nSpread Reduction ",SpreadReduction)


InvasionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(InvasionSummary) = c("Realisation",   "Year",  "NodesInfested")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,perm]
dim(InvasionData)
NodesInfested = colSums(InvasionData)
Realisation = perm 
Year = 1:Nyears
Results = data.frame(Realisation,Year,NodesInfested)
InvasionSummary = rbind(InvasionSummary,Results)
}

Filename = paste0(FileNameStem,"InvasionRaw.png")
png(Filename)
plot(InvasionSummary$Year,InvasionSummary$NodesInfested,ylim = c(0,max(InvasionSummary$NodesInfested)),pch = NA
, xlab = "Time since incursion detected (years)",
ylab = "Number of nodes infested", main = Title)

for(perm in 1:Nperm)
{
Sub = InvasionSummary[InvasionSummary$Realisation == perm,]
lines(Sub$Year,Sub$NodesInfested,col  = perm)
}
if(DetectionProb == 0)
{
#points(13,56*0.67,col = 1,cex = 2,pch = 19) ###Historic N farms from Bell 2006 https://nzpps.org/_journal/index.php/nzpp/article/view/4417/4245
					    ###67% of infested paddocks under grazing
#points(18,(96-20)*0.67,col = 1,cex = 2,pch = 19) ###Subtract 20 here as 20 new incursions due to subdivisions 
}

dev.off()

Quantiles = as.data.frame(aggregate(InvasionSummary$NodesInfested, by = list(InvasionSummary$Year),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"InvasionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (years)",
ylab = "Number of nodes infested", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
if(DetectionProb == 0)
{
#points(13,56*0.67,col = 1,cex = 2,pch = 19) ###Historic N farms from Bell 2006 https://nzpps.org/_journal/index.php/nzpp/article/view/4417/4245
					    ###67% of infested paddocks under grazing
#points(18,(96-20)*0.67,col = 1,cex = 2,pch = 19) ###Subtract 20 here as 20 new incursions due to subdivisions 
}
dev.off()


ManagingSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(ManagingSummary) = c("Realisation",   "Year",  "NodesManaging")

for(perm in 1:Nperm)
{
ManagingData = ManagingResults[,,perm]
dim(ManagingData)
NodesManaging = colSums(ManagingData)
Realisation = perm 
Year = 1:Nyears
Results = data.frame(Realisation,Year,NodesManaging)
ManagingSummary = rbind(ManagingSummary,Results)
}


Filename = paste0(FileNameStem,"ManagingRaw.png")
png(Filename)
plot(ManagingSummary$Year,ManagingSummary$NodesManaging,ylim = c(0,max(ManagingSummary$NodesManaging)),pch = NA
, xlab = "Time since incursion detected (years)",
ylab = "Nodes under management", main = Title)

for(perm in 1:Nperm)
{
Sub = ManagingSummary[ManagingSummary$Realisation == perm,]
lines(Sub$Year,Sub$NodesManaging,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(ManagingSummary$NodesManaging, by = list(ManagingSummary$Year),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"ManagingSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (years)",
ylab = "Nodes under management", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


DetectedSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedSummary) = c("Realisation",   "Year",  "NodesDetected")

for(perm in 1:Nperm)
{
DetectedData = DetectedResults[,,perm]
dim(DetectedData)
NodesDetected = colSums(DetectedData)
Realisation = perm 
Year = 1:Nyears
Results = data.frame(Realisation,Year,NodesDetected)
DetectedSummary = rbind(DetectedSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedRaw.png")
png(Filename)
plot(DetectedSummary$Year,DetectedSummary$NodesDetected,ylim = c(0,max(DetectedSummary$NodesDetected)),pch = NA
, xlab = "Time since incursion detected (years)",
ylab = "Nodes pest detected", main = Title)

for(perm in 1:Nperm)
{
Sub = DetectedSummary[DetectedSummary$Realisation == perm,]
lines(Sub$Year,Sub$NodesDetected,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(DetectedSummary$NodesDetected, by = list(DetectedSummary$Year),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"DetectedSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (years)",
ylab = "Nodes pest detected", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


DetectedProportionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedProportionSummary) = c("Realisation",   "Year",  "DetectedProportion")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,perm]
DetectedData = DetectedResults[,,perm]
NodesDetected = colSums(DetectedData)
NodesInvaded = colSums(InvasionData)
DetectedProportion = NodesDetected/NodesInvaded
DetectedProportion[is.na(DetectedProportion)==T] = 1
Realisation = perm 
Year = 1:Nyears
Results = data.frame(Realisation,Year,DetectedProportion)
DetectedProportionSummary = rbind(DetectedProportionSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedProportionRaw.png")
png(Filename)
plot(DetectedProportionSummary$Year,DetectedProportionSummary$DetectedProportion,ylim = c(0,max(DetectedProportionSummary$DetectedProportion)),pch = NA
, xlab = "Time since incursion detected (years)",
ylab = "Proportion of infested nodes detected", main = Title)

for(perm in 1:Nperm)
{
Sub = DetectedProportionSummary[DetectedProportionSummary$Realisation == perm,]
lines(Sub$Year,Sub$DetectedProportion,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(DetectedProportionSummary$DetectedProportion, by = list(DetectedProportionSummary$Year),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"DetectedProportionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (years)",
ylab = "Proportion infested nodes detected", main = Title)
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

