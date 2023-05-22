###########################################################################
###########################################################################
###Declares a function simulating pest spread and management between areas of suitable climate and land use
###Key inputs are: 
###1) matrix of dispersal probibilities between each pair of sites (i.e. nodes of the network). Matrix can be non-symmetrical (i.e. can have source and sink nodes) 
###2) Envionmentally-determined establishment probability. Set all values to 1 for no environmental limitation on establishment.
###3) Annual per-farm long distance dispersal rate or matrix of annual LDD rates for each pair of farms 
###4) Management parameters
### a) Annual detection probability
### b) Annual management adoption probability subsequent to detection
### c) Annual eradication probability
### d) Spread reduction when management adopted
###5) A binary vector of intially infested nodes OR
###   A vector of invasion risk probabilities and/or proportion of nodes initially invaded
###Key outputs are:
###3-dimensional arrays of invasion, management and detection status for each node in each year of each permuation
###2-dimensional array of invasion probability (i.e. proportion of permutations pest present) for each node in each year
###Line graphs summarising number of nodes infested, where infestations
###are detected, the proportion of infestations detected and nodes under management against time 
###########################################################################
###########################################################################


INApest = function(
ModelName,              #Name for storing results to file 
Nperm,                  #Number of permutations per parameter combination
Nyears,                 #Simulation duration
DetectionProb,          #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb,             #Annual Probability or vector of probabilities vector length nrow(SDDprob) of node adopting management upon detection
ManageSD = NULL, #Option to provide standard deviation for management probability can be vector (nodes)
AnnualEradicationProb, #Annual probability of eradication (must be between 0 and 1) when management adopted. Only single value possibe
AnnualEradicationSD = NULL, #Option to provide standard deviation for annual probability of eradicatio. Only single value possibe
SpreadReduction,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
SpreadReductionSD = NULL,        #Option to provide standard deviation for spread reduction can be vector (nodes)
InitialInvasion = NA,        #Nodes infested at start of simulations
InitBioP = NA,		#Proportion of nodes infested at start of simulations
InvasionRisk = NA,           #Vector (nodes) or matrix (nodes x years) of invasion risk from external sources
EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x years)
LocalExtinction = 0,           #Annual local population extinction probability. Set to zero for no local extinction. Can be vector (nodes) or matrix (nodes x years)
SDDprob,                   #Biophysical adjaceny matrix - disperal probability between each pair of nodes
SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
LDDprob = 0,         #Option to provide long distance dispersal matrix instead of distance-independent dispesal rate
			      #e.g. could be weighted by law of human visitation or data on stock movements
geocoords,              #XY points for INAscene
OngoingExternal = F,   ##Option to include ongoing invasion from external sources
OutputDir = NA,		      #Directory for storing results
DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
{
###POTENTIAL ADDITIONS
###1) Allow provision of ManageProb, DetectionProb and SpreadReduction as matrices (nodes x years) to allow temporal variation in detection and management adoption 
###   (e.g. as awareness of the pest increases), and spread reduction (e.g. as hygeine practices improve or decline) 
###2) Allow provision of AnnualEradicationProb as a vector (years) to allow temporal varation in eradication 
###   (as new technologies emerge e.g. pesticides and biocontrol or the pest develops resistance) 
###Allow provision of management variables (and SD values), EnvEstabProb and LocalExtinction as vectors (length = Nyears), for simple temporal variation


###Declare array tracking infestation status 
###of individual nodes in each year of each realisation
InvasionResults = array(dim = c(nrow(SDDprob),Nyears,Nperm))

###Declare array tracking detection status 
###of individual nodes in each year of each realisation
DetectedResults = InvasionResults

###Declare array for tracking management adoption status 
###of individual nodes in each year of each realisation
###This is a measure of potential disruption to farm businesses
###or ongoing management burden (surveillance and removal)
###for publicly-owned lands
ManagingResults = InvasionResults


###If LDD matrix provided combine long and short distance dispersal probability 
###into a single dispersal probability matrix 
DispProb = SDDprob
if(is.matrix(LDDprob) == T)
  DispProb =1-(1-SDDprob)*(1-LDDprob)
    

###Weight dispersal by environmental establishment probability
if(is.matrix(EnvEstabProb) == F)
  BPAM = sweep(DispProb,2,EnvEstabProb,`*`)


###SEAM is needed to run INA but zero info transfer between nodes is assumed unless SEAM provided
if(is.matrix(SEAM) == F)
{
SEAM = BPAM
SEAM[,] = 0
}

###Assign standard deviation value to management in no value provided
if(is.null(ManageSD) == T)
	ManageSD = mean(ManageProb)/10
if(is.null(SpreadReductionSD) == T)
	SpreadReductionSD = (1-mean(SpreadReduction))/10
if(is.null(DetectionSD) == T)
	DetectionSD = mean(DetectionProb)/10
if(is.null(AnnualEradicationSD) == T)
	AnnualEradicationSD = AnnualEradicationProb/10

###Set annual local population survival rate
if(is.matrix(LocalExtinction) == F)
  Survival = 1-LocalExtinction

###########################################################
###Start of simulation loop
###########################################################

###Loop through each realisation
###Nodes with information (which have detected infestation) at start of simulations
###vary across realisations 
for(perm in 1:Nperm)
{
###Assign initial infestations according either to "InitialInvasion" binary vector OR
###"InvasionRisk" probabilities and/or initial proportion of nodes infested ("InitBioP") OR
###just "InitBioP" if neither "InitialInvasion" or "InvasionRisk" supplied by user
InitBio = rep(0,times = nrow(SDDprob))

###Binary vector of initial infestation status not provided?
if(length(InitialInvasion) != nrow(SDDprob))
{
###Vector of invasion risk probabilities used to select initial ifestations 
if(length(InvasionRisk) == nrow(SDDprob))
        {
 	###If initial infestation proportion provided, use wieghted randomisation to select the desired number of infested nodes
        if(is.na(InitBioP) == F)
	  Infested = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitBioP),prob = InvasionRisk)
	###If initial infestation proportion not provided, use random binomial process to definie initial infestation status
        if(is.na(InitBioP) == T)
          {
	  Infested = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
          Infested = which(Infested == 1)
	  } 
	}
###Invasion risk probabilities either not provided or provided as matrix (nodes x years)
 if(length(InvasionRisk) != nrow(SDDprob))
        {
 	###If initial infestation proportion provided with no invasion risk probability vector
        ###use unweighted randomisation to select infested nodes
        if(is.matrix(InvasionRisk) == F)
	 Infested = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitBioP))
       ###If invasion risk supplied as matrix, use first column to randomly select initial infestations via random binomial process
       if(is.matrix(InvasionRisk) == T)
          {
          Infested = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,1])
          Infested = which(Infested == 1)
	  }
	}
InitBio[Infested] = 1
}

###Assign initial infestations using binary vector
if(length(InitialInvasion) == nrow(SDDprob))
	InitBio = InitialInvasion


###Randomly assign annual detection probability, based on mean and sd
NodeDetectionProb = rnorm(DetectionProb,DetectionSD,n = nrow(SDDprob))
NodeDetectionProb[NodeDetectionProb<0] = 0
NodeDetectionProb[NodeDetectionProb>1] = 1


###Randomly assign probability of mangement adoption upon detection of infestation
NodeManageProb = rnorm(ManageProb,ManageSD,n = nrow(SDDprob))
NodeManageProb[NodeManageProb<0] = 0
NodeManageProb[NodeManageProb>1] = 1

###Randomly assign spread reduction factor when management adopted
NodeSpreadReduction = rnorm(SpreadReduction,SpreadReductionSD,n = nrow(SDDprob))
NodeSpreadReduction[NodeSpreadReduction<0] = 0
NodeSpreadReduction[NodeSpreadReduction>1] = 1

###Probability of info at start of simulation depends on
###Presence of pest and detection probability
InitInfoProb = InitBio*NodeDetectionProb

###Select nodes that have detected infestation 
InitInfo = vector(length = nrow(SDDprob))
RandInitInfoProb = runif(n=nrow(SDDprob),0,1)
InitInfo = ifelse(RandInitInfoProb>InitInfoProb,0,1) 

###Populate Invasion and info status vectors ahead of year loop
Invaded = InitBio 
HaveInfo = InitInfo

###Loop through years
###This allows allocation of info to farmers or managers based on detection of infestation
###Non-infested nodes don't manage unless they receive information from nodes with known extant infestations
###Nodes may continue to manage even post eradication
for(year in 1:Nyears)
  {
  ###Print progress
  cat("\r", "Realisation ", perm, "Year ", year, "...")
  
  ###Allow for variation in establishment through time
  ###e.g. annual climate change predictions
  ###Note: could be done outside loop, but would take heaps of memory to store 
  if(is.matrix(EnvEstabProb) == T)
    BPAM = sweep(DispProb,2,EnvEstabProb[,year],`*`)
  	  
    
  ###Set annual local population survival rate
  if(is.matrix(LocalExtinction) == T)
    Survival = 1-LocalExtinction[,year]
  SurvivalIsVector = NA
  ###Need to declare these variables to use survival appropriately in INAscene call
  if(length(Survival) == 1)
      {
      SurvivalIsVector = F
      SurvivalMean = Survival
      SurvivalSD = 0.00000001
      }
  if(length(Survival) == nrow(SDDprob))
      {
      SurvivalIsVector = T
      SurvivalVector = Survival
      SurvivalMean = NA
      SurvivalSD = NA
      }

  ###Assign management status to nodes   
  ###This assumes that individual farmers or managers may not adopt management every single year
  ###Allows for inconsistency in implementation 
  RandManageProb = runif(0,1,n=nrow(SDDprob))
  Managing = vector(length = nrow(SDDprob))
  Managing = ifelse(RandManageProb>NodeManageProb,0,1)	
  
  ###Management is only applied to nodes which have information
  ###i.e. where pest has been detected or following communication of information
  ###from neighbouring infested farms 
  Managing = Managing*HaveInfo
  
  ###Identify nodes with known extant infestations 
  Detected = Invaded*HaveInfo

  INAsceneLarge <-
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
    bpam =  BPAM*(1-Managing*NodeSpreadReduction), ##biophysical adjacency matrix moderated by spread reduction in managing nodes
    bpamdist = F,
    bpamrandp = NA,
    bpampla = NA,
    bpamplb = NA,
    readprobadoptvec = T,
    probadoptvec = Managing, ##Use actual adoption calculated externally instead of adoption probability
    probadoptmean = NA,
    probadoptsd = NA,
    readprobestabvec = SurvivalIsVector,
    probestabvec = SurvivalVector,
    probestabmean = SurvivalMean,###set this to <1 to allow pops to die out on their own.
    probestabsd = SurvivalSD,
    maneffdir = 'decrease_estab',
    maneffmean = AnnualEradicationProb, ##Set mean management efficacy - effectively annual eradication prob
    maneffsd = AnnualEradicationSD,
    usethreshman = F,
    maneffthresh = NA,
    sampeffort = NA
  )
 ###Extract INAscene results for each node  
 LargeOut = INAsceneLarge$multdetails
  
 ###Update info vector for any info spread (if SEAM supplied)
 ###Note once nodes obtain info they always have info (only zero values updated)
 InfoOut = as.vector(LargeOut[[1]]$multout[[1]]$vect1cL[[2]])
 HaveInfo[HaveInfo == 0] = InfoOut[HaveInfo == 0]
 
 ###Update infestation vector
 Estab = as.array(LargeOut[[1]][[1]][[1]]$estabvecL)
 Invaded = ifelse(Estab[[1]]==FALSE,0,1)
 
 ###Add invasion resulting from colonisation from external sources
 if(OngoingExternal == T)
  {
  if(is.matrix(InvasionRisk) == F)
    ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
  ###If invasion risk supplied as matrix, use year column to randomly assign new infestations via random binomial process
  if(is.matrix(InvasionRisk) == T)
    ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,year])
  Invaded[Invaded == 0] = ExternalInvasion[Invaded==0]
  }
 ###Record nodes adopting management
 ManagingResults[,year,perm] = Managing
  
 ###Record infested nodes
 InvasionResults[,year,perm] = Invaded
  
 ###Select new nodes where infestation detected
 NewInfoProb = Invaded*DetectionProb
 RandNewInfoProb = runif(n=nrow(SDDprob),0,1)
 NewHaveInfo = vector(length= length(HaveInfo))
 NewHaveInfo =  ifelse(RandNewInfoProb>NewInfoProb,0,1) 
 
 ###Add newly detected infestations to info vector
 ###Note once nodes obtain info they always have info (only zero values updated)
 HaveInfo[HaveInfo==0] = NewHaveInfo[HaveInfo==0]  
 
 ###Record detection status
 DetectedResults[,year,perm] = HaveInfo*Invaded
 }
}
###########################################################
###End of simulation loop
###########################################################

###########################################################
###Save results for post-hoc analyses
###########################################################
###ModelName used to generate filenames
###Use standard format for ease of reading results to produce heat maps 
###and conduct post-hoc stats comparing managment settings/scenarios 
if(is.na(OutputDir) == T)
	OutputDir = ""
FileNameStem = paste0(OutputDir,ModelName)

###These are 3D arrays with dimensions (Nodes,Years,Realisations)
saveRDS(ManagingResults, paste0(FileNameStem,"InfoLargeOut.rds"))
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))

##########################################################
###Store annual farm-level invasion probs for heat maps
###and estimation of invasion threat to other regions
##########################################################


InvasionProb = matrix(ncol = Nyears, nrow = nrow(SDDprob))
for(year in 1:Nyears)
{
YearData = InvasionResults[,year,]
InvasionProb[,year] = rowSums(YearData)/Nperm
}
saveRDS(InvasionProb, paste0(FileNameStem,"InvasionProb.rds"))

if(DoPlots == T)
{
###########################################################
###Produce summary figs when processing completed
###########################################################

Title = paste0("Detection Prob. ",round(mean(DetectionProb),digits = 3),
	" Manage Prob. ",round(mean(ManageProb),digits = 3),
       "\nEradication Prob. ",round(AnnualEradicationProb,digits =3),
		" Spread Reduction ",round(mean(SpreadReduction),digits=3))


###Change in number of nodes infested with time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

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

###Option to plot points for Marlborough historic data
if(mean(DetectionProb) == 0)
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

###Option to plot points for Marlborough historic data
if(mean(DetectionProb) == 0)
{
#points(13,56*0.67,col = 1,cex = 2,pch = 19) ###Historic N farms from Bell 2006 https://nzpps.org/_journal/index.php/nzpp/article/view/4417/4245
					    ###67% of infested paddocks under grazing
#points(18,(96-20)*0.67,col = 1,cex = 2,pch = 19) ###Subtract 20 here as 20 new incursions due to subdivisions 
}
dev.off()


###Change in number of farms managing through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided

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


###Change in number of known extant infestations through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

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


###Change in proportion of extant infestations detected through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided
 
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
}


################################################################
################################################################
###End of function
################################################################
################################################################

