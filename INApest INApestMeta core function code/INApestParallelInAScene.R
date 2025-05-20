###########################################################################
###########################################################################
###Declares a function simulating pest spread and management between areas of suitable climate and land use
###Key inputs are: 
###1) matrix of short-distance (self-mediated) dispersal probibilities between each pair of sites (i.e. nodes of the network). Matrix can be non-symmetrical (i.e. can have source and sink nodes) 
###2) Envionmentally-determined establishment probability. Set all values to 1 for no environmental limitation on establishment.
###3) matrix of long-distance (human-mediated) dispersal probibilities between each pair of sites (i.e. nodes of the network). Matrix can be non-symmetrical 
###4) Management parameters
### a)  detection probability
### b)  management adoption probability subsequent to detection
### c)  eradication probability
### d) Spread reduction when management adopted
###5) A binary vector of intially infested nodes OR
###   A vector of invasion risk probabilities and/or proportion of nodes initially invaded
###Key outputs are:
###3-dimensional arrays of invasion, management and detection status for each node in each timestep of each permuation
###2-dimensional array of invasion probability (i.e. proportion of permutations pest present) for each node in each timestep
###Line graphs summarising number of nodes infested, where infestations
###are detected, the proportion of infestations detected and nodes under management against time 
###########################################################################
###########################################################################

#######################################################################
###This implements parallel processing in foreach for the permutation loop.
###See file "ParallelSetup.r" for notes on steps for setting up parallel processing
#######################################################################

###Load the required packages 
library(abind)
library(INA)
library(doParallel)


INApestParallelINAscene = function(
ModelName,              #Name for storing results to disk 
Nperm,                  #Number of permutations per parameter combination
Ntimesteps,                 #Simulation duration
DetectionProb,          # detection probability (must be between 0 and 1). Can be single number, vector (nodes) or matrix (nodes x timesteps)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be single number or vector (nodes)
ManageProb,             # Probability adopting management upon detection. Can be single number, vector (nodes) or matrix (nodes x timesteps)
ManageSD = NULL, #Option to provide standard deviation for management probability can be single number or vector (nodes)
EradicationProb, # probability of eradication (must be between 0 and 1) when management adopted. Can be single number, vector (nodes) or matrix (nodes x timesteps)
EradicationSD = NULL, #Option to provide standard deviation for  probability of eradication. Can be single number or vector (nodes)
SpreadReduction,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector (nodes) or matrix (nodes x timesteps)
SpreadReductionSD = NULL,        #Option to provide standard deviation for spread reduction can be single number or vector (nodes)
InitialInvasion = NA,        #Nodes infested at start of simulations
InitBioP = NA,		#Proportion of nodes infested at start of simulations
InvasionRisk = NA,           #Vector (nodes) or matrix (nodes x timesteps) of invasion risk from external sources
InitialInfo = NA,        #Vector or of nodes with information at start of simulations
InitInfoP = NA,		#Proportion of nodes with information at start of simulations
ExternalInfoProb = NA,           #Vector of probabilities of communication from external sources
EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
SDDprob,                   #Short-distance (self-mediated) disperal probability between each pair of nodes
SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
LDDprob = 0,         #Option to provide long-distance (human-mediated) dispersal probability matrix
			      #e.g. could be weighted by law of human visitation or data on stock movements
geocoords,              #XY points for INAscene
OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
OngoingExternalInfo = F,   ##Option to include ongoing communication from external sources
OutputDir = NA,		      #Directory for storing results
DoPlots = TRUE	     #Option to omit printing of line graphs. Default is to print.
)
{
###POTENTIAL ADDITIONS
###1) Allow for increased/decreased management adoption if neighbours, or contacts in the social network are managing
###2) Might also be interesting to allow for other aspects of management response to vary depending on information held by neighbours
####e.g. detection probability might increase if neighbours area managing

###Changes for final report
###1) Change local extinction to survival
###2) Change timesteps to steps 

###Declare array tracking infestation status 
###of individual nodes in each timestep of each realisation
InvasionResults = array(dim = c(nrow(SDDprob),Ntimesteps,Nperm))

###Declare array tracking detection status 
###of individual nodes in each timestep of each realisation
DetectedResults = InvasionResults

###Declare array for tracking management adoption status 
###of individual nodes in each timestep of each realisation
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
if(is.null(EradicationSD) == T)
	EradicationSD = mean(EradicationProb)/10

###Set  local population survival rate
if(is.matrix(Survival) == F &&(length(Survival) == 1 ||length(Survival) == nrow(SDDprob)))
   NodeSurvival = Survival

###########################################################
###Start of simulation loop
###########################################################

#########################################################
###Use parallel processing in foreach loop to increase processing speed
###Combines results for invasion, managing and detection status in 4d array
#########################################################

totalCores = detectCores()
#Leave one core to avoid overload your computer
cluster <- makeCluster(totalCores[1]-1) 
registerDoParallel(cluster)
###Loop through each realisation
###Nodes with information (which have detected infestation) at start of simulations
###vary across realisations 

###Define a function for combining results of each permutation
acomb <- function(...) abind(..., along=4)

###Need to include required packages in the .packages arguement of the foreach call
PermOut <- foreach(1:Nperm, .combine = 'acomb',.packages=c("abind","INA")) %dopar% 
{
InvasionResultsLoop <- array(dim = c(nrow(SDDprob),Ntimesteps))
ManagingResultsLoop <- InvasionResultsLoop
DetectedResultsLoop <- InvasionResultsLoop

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
###Invasion risk probabilities either not provided or provided as matrix (nodes x timesteps)
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

###Select nodes with information at start of simulation  according either to "InitialInfo" binary vector OR
###"ExternalInfoProb" probabilities and/or initial proportion of nodes with information ("InitInfoP") OR
###just "InitInfoP" if neither "InitialInfo" or "ExternalInfoProb" supplied by user.
###If no initial info variables provided, no nodes have info at start of simulations
InitInfo = rep(0,times = nrow(SDDprob))
if(is.na(sum(InitialInfo))== F || is.na(InitInfoP) == F || is.na(sum(ExternalInfoProb)) == F )
{
  if(length(InitialInfo) != nrow(SDDprob))
  {
    if(length(ExternalInfoProb) == nrow(SDDprob))
    {
      if(is.na(InitInfoP) == F)
        Info = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitInfoP),prob = ExternalInfoProb)
      if(is.na(InitInfoP) == T)
      {
        Info = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb) 
        Info = which(Info == 1)
      } 
    }
    if(length(ExternalInfoProb) != nrow(SDDprob))
    {
      if(is.matrix(ExternalInfoProb) == F)
        Info = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitInfoP))
      if(is.matrix(ExternalInfoProb) == T)
      {
        Info = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb[,1])
        Info = which(Info == 1)
      }
    }
    InitInfo[Info] = 1
    
  }
  if(length(InitialInfo) == nrow(SDDprob))
    InitInfo = InitialInfo  
}


###Assign initial infestations using binary vector
if(length(InitialInvasion) == nrow(SDDprob))
	InitBio = InitialInvasion

###Randomly assign  detection probability, based on mean and sd
###If DetectionProb given as single value or vector (nodes)
if(is.matrix(DetectionProb)==FALSE &&(length(DetectionProb) == 1 ||length(DetectionProb) == nrow(SDDprob) ))
      {
	NodeDetectionProb = rnorm(DetectionProb,DetectionSD,n = nrow(SDDprob))
      NodeDetectionProb[NodeDetectionProb<0] = 0
      NodeDetectionProb[NodeDetectionProb>1] = 1
      }

###If DetectionProb given as matrix (nodes x timesteps) use values for first timestep to get initial detections
if(is.matrix(DetectionProb)==TRUE && nrow(DetectionProb) == nrow(SDDprob) && ncol(DetectionProb) == Ntimesteps)
      {
	NodeDetectionProb = rnorm(DetectionProb[,1],DetectionSD,n = nrow(SDDprob))
      NodeDetectionProb[NodeDetectionProb<0] = 0
      NodeDetectionProb[NodeDetectionProb>1] = 1
      }


###Randomly assign probability of mangement adoption upon detection of infestation
###If ManageProb given as single value or vector (nodes)
if(is.matrix(ManageProb)==FALSE &&(length(ManageProb) == 1 ||length(ManageProb) == nrow(SDDprob) ))
      {
	NodeManageProb = rnorm(ManageProb,ManageSD,n = nrow(SDDprob))
      NodeManageProb[NodeManageProb<0] = 0
      NodeManageProb[NodeManageProb>1] = 1
      }

###Randomly assign spread reduction factor when management adopted
###If SpreadReduction given as single value or vector (nodes)
if(is.matrix(SpreadReduction)==FALSE &&(length(SpreadReduction) == 1 ||length(SpreadReduction) == nrow(SDDprob) ))
      {
	NodeSpreadReduction = rnorm(SpreadReduction,ManageSD,n = nrow(SDDprob))
      NodeSpreadReduction[NodeSpreadReduction<0] = 0
      NodeSpreadReduction[NodeSpreadReduction>1] = 1
      }

###Randomly assign  eradication probability when management applied
###If EradicationProb given as single value or vector (nodes)
if(is.matrix(EradicationProb)==FALSE &&(length(EradicationProb) == 1 ||length(EradicationProb) == nrow(SDDprob) ))
      {
      NodeEradicationProb = rnorm(EradicationProb,EradicationSD,n = nrow(SDDprob))
      NodeEradicationProb[NodeEradicationProb<0] = 0
      NodeEradicationProb[NodeEradicationProb>1] = 1
      }



###Probability of info at start of simulation depends on
###Presence of pest and detection probability
###Select nodes that have detected infestation 
InitDetection = rbinom(1:nrow(SDDprob),size = 1,prob = InitBio*NodeDetectionProb)
###Add detections to nodes which already have info (e.g. pre-emptive control and hygiene measures)
InitInfo[InitInfo == 0] = InitDetection[InitInfo == 0]

###Populate Invasion and info status vectors ahead of timestep loop
Invaded = InitBio 
HaveInfo = InitInfo

###Probability of info at start of simulation depends on
###Presence of pest and detection probability
###Select nodes that have detected infestation 
InitDetection = rbinom(1:nrow(SDDprob),size = 1,prob = Invaded*NodeDetectionProb)
###Add detections to nodes which already have info (e.g. pre-emptive control and hygiene measures)
InitInfo[InitInfo == 0] = InitDetection[InitInfo == 0]

###Loop through timesteps
###This allows allocation of info to farmers or managers based on detection of infestation
###Non-infested nodes don't manage unless they receive information from nodes with known extant infestations
###Nodes may continue to manage even post eradication
for(timestep in 1:Ntimesteps)
  {
  ###Print progress
  #cat("\r", "Realisation ", perm, "Timestep ", timestep, "...")
  
  ###Allow for variation in establishment through time
  ###e.g.  climate change predictions
  ###Note: could be done outside loop, but would take heaps of memory to store 
  if(is.matrix(EnvEstabProb) == T)
    BPAM = sweep(DispProb,2,EnvEstabProb[,timestep],`*`)
  	  
    
  ###Set  local population survival rate
  if(is.matrix(Survival) == T && nrow(Survival) == nrow(SDDprob) && ncol(Survival) == Ntimesteps)
    NodeSurvival = Survival[,timestep]

     

  ###Randomly assign  detection probability, based on mean and sd
  ###If DetectionProb given as matrix (nodes x timesteps)
  if(is.matrix(DetectionProb)==TRUE && nrow(DetectionProb) == nrow(SDDprob) && ncol(DetectionProb) == Ntimesteps)
   {	
   NodeDetectionProb = rnorm(DetectionProb[,timestep],DetectionSD,n = nrow(SDDprob))
   NodeDetectionProb[NodeDetectionProb<0] = 0
   NodeDetectionProb[NodeDetectionProb>1] = 1
   }

  ###Randomly assign probability of mangement adoption upon detection of infestation
  ###If ManageProb given as matrix (nodes x timesteps)
  if(is.matrix(ManageProb)==TRUE && nrow(ManageProb) == nrow(SDDprob) && ncol(ManageProb) == Ntimesteps)
   {	
   NodeManageProb = rnorm(ManageProb[,timestep],ManageSD,n = nrow(SDDprob))
   NodeManageProb[NodeManageProb<0] = 0
   NodeManageProb[NodeManageProb>1] = 1
   }

  ###Randomly assign spread reduction factor when management adopted
  ###If SpreadReduction given as matrix (nodes x timesteps)
  if(is.matrix(SpreadReduction)==TRUE && nrow(SpreadReduction) == nrow(SDDprob) && ncol(SpreadReduction) == Ntimesteps)
   {	
   NodeSpreadReduction = rnorm(SpreadReduction[,timestep],ManageSD,n = nrow(SDDprob))
   NodeSpreadReduction[NodeSpreadReduction<0] = 0
   NodeSpreadReduction[NodeSpreadReduction>1] = 1
   }
  
  ###Randomly assign  eradication probability when management applied
  ###If EradicationProb given as matrix (nodes x timesteps)
  if(is.matrix(EradicationProb)==TRUE && nrow(EradicationProb) == nrow(SDDprob) && ncol(EradicationProb) == Ntimesteps)
      {
	NodeEradicationProb = rnorm(EradicationProb[,timestep],EradicationSD,n = nrow(SDDprob))
      NodeEradicationProb[NodeEradicationProb<0] = 0
      NodeEradicationProb[NodeEradicationProb>1] = 1
      }

 
  ###Management is only applied to nodes which have information
  ###i.e. where pest has been detected or following communication of information
  ###from neighbouring infested farms 
  Managing = rbinom(1:nrow(SDDprob),size = 1,prob = NodeManageProb*HaveInfo)
  
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
    readprobestabvec = T,
    probestabvec = NodeSurvival*(1-NodeEradicationProb*Managing),###set SurvivalVector values to <1 to allow pops to die out on their own.
                   ###Moderate survival by node-specific eradication  
    probestabmean = NA,
    probestabsd = NA,
    maneffdir = 'decrease_estab',
    maneffmean = 0, ##Set management efficacy to zero - eradication applied by moderating survival
    maneffsd = 0.00000000000001,
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
 if(OngoingExternalInvasion == T)
  {
  if(is.matrix(InvasionRisk) == F)
    ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
  ###If invasion risk supplied as matrix, use timestep column to randomly assign new infestations via random binomial process
  if(is.matrix(InvasionRisk) == T)
    ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,timestep])
  Invaded[Invaded == 0] = ExternalInvasion[Invaded==0]
  }
 
 ###Add nodes with information resulting from external sources
 if(OngoingExternalInfo == T)
   {
   if(is.matrix(ExternalInfoProb) == F)
     ExternalInfo = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb)
   if(is.matrix(ExternalInfoProb) == T)
     ExternalInfo = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb[,timestep])
   HaveInfo[HaveInfo == 0] = ExternalInfo[HaveInfo==0]
   }
 
 ###Record nodes adopting management
 ManagingResultsLoop[,timestep] = Managing
  
 ###Record infested nodes
 InvasionResultsLoop[,timestep] = Invaded
  
 ###Select new nodes where infestation detected
 NewHaveInfo =  rbinom(1:nrow(SDDprob),size = 1,prob = Invaded*NodeDetectionProb)
 
 ###Add newly detected infestations to info vector
 ###Note once nodes obtain info they always have info (only zero values updated)
 HaveInfo[HaveInfo==0] = NewHaveInfo[HaveInfo==0]  
 
 ###Record detection status
 DetectedResultsLoop[,timestep] = HaveInfo*Invaded
}
abind(InvasionResultsLoop,ManagingResultsLoop,DetectedResultsLoop,along = 3)
}
stopCluster(cluster)

###Extract results for invasion, managing and dectection status into separate 3d arrays
InvasionResults <- PermOut[,,1,]
ManagingResults <- PermOut[,,2,]
DetectedResults <- PermOut[,,3,]
dim(InvasionResults)
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

###These are 3D arrays with dimensions (Nodes,Timesteps,Realisations)
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(ManagingResults, paste0(FileNameStem,"ManagingLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))

##########################################################
###Store  farm-level invasion probs for heat maps
###and estimation of invasion threat to other regions
##########################################################


InvasionProb = matrix(ncol = Ntimesteps, nrow = nrow(SDDprob))
for(timestep in 1:Ntimesteps)
{
TimestepData = InvasionResults[,timestep,]
InvasionProb[,timestep] = rowSums(TimestepData)/Nperm
}
saveRDS(InvasionProb, paste0(FileNameStem,"InvasionProb.rds"))

if(DoPlots == T)
{
###########################################################
###Produce summary figs when processing completed
###########################################################

Title = paste0("Detection Prob. ",round(mean(DetectionProb),digits = 3),
	" Manage Prob. ",round(mean(ManageProb),digits = 3),
       "\nEradication Prob. ",round(mean(EradicationProb),digits =3),
		" Spread Reduction ",round(mean(SpreadReduction),digits=3))


###Change in number of nodes infested with time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

InvasionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(InvasionSummary) = c("Realisation",   "Timestep",  "NodesInfested")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,perm]
dim(InvasionData)
NodesInfested = colSums(InvasionData)
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesInfested)
InvasionSummary = rbind(InvasionSummary,Results)
}


Filename = paste0(FileNameStem,"InvasionRaw.png")
png(Filename)
plot(InvasionSummary$Timestep,InvasionSummary$NodesInfested,ylim = c(0,max(InvasionSummary$NodesInfested)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Number of nodes infested", main = Title)

for(perm in 1:Nperm)
{
Sub = InvasionSummary[InvasionSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesInfested,col  = perm)
}

###Option to plot points for Marlborough historic data
if(mean(DetectionProb) == 0)
{
#points(13,56*0.67,col = 1,cex = 2,pch = 19) ###Historic N farms from Bell 2006 https://nzpps.org/_journal/index.php/nzpp/article/view/4417/4245
					    ###67% of infested paddocks under grazing
#points(18,(96-20)*0.67,col = 1,cex = 2,pch = 19) ###Subtract 20 here as 20 new incursions due to subdivisions 
}

dev.off()

Quantiles = as.data.frame(aggregate(InvasionSummary$NodesInfested, by = list(InvasionSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"InvasionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
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
colnames(ManagingSummary) = c("Realisation",   "Timestep",  "NodesManaging")

for(perm in 1:Nperm)
{
ManagingData = ManagingResults[,,perm]
dim(ManagingData)
NodesManaging = colSums(ManagingData)
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesManaging)
ManagingSummary = rbind(ManagingSummary,Results)
}



 
Filename = paste0(FileNameStem,"ManagingRaw.png")
png(Filename)
plot(ManagingSummary$Timestep,ManagingSummary$NodesManaging,ylim = c(0,max(ManagingSummary$NodesManaging)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes under management", main = Title)

for(perm in 1:Nperm)
{
Sub = ManagingSummary[ManagingSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesManaging,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(ManagingSummary$NodesManaging, by = list(ManagingSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"ManagingSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes under management", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


###Change in number of known extant infestations through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

DetectedSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedSummary) = c("Realisation",   "Timestep",  "NodesDetected")

for(perm in 1:Nperm)
{
DetectedData = DetectedResults[,,perm]
dim(DetectedData)
NodesDetected = colSums(DetectedData)
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesDetected)
DetectedSummary = rbind(DetectedSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedRaw.png")
png(Filename)
plot(DetectedSummary$Timestep,DetectedSummary$NodesDetected,ylim = c(0,max(DetectedSummary$NodesDetected)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes pest detected", main = Title)

for(perm in 1:Nperm)
{
Sub = DetectedSummary[DetectedSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesDetected,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(DetectedSummary$NodesDetected, by = list(DetectedSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"DetectedSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes pest detected", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


###Change in proportion of extant infestations detected through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided
 
DetectedProportionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedProportionSummary) = c("Realisation",   "Timestep",  "DetectedProportion")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,perm]
DetectedData = DetectedResults[,,perm]
NodesDetected = colSums(DetectedData)
NodesInvaded = colSums(InvasionData)
DetectedProportion = NodesDetected/NodesInvaded
DetectedProportion[is.na(DetectedProportion)==T] = 1
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,DetectedProportion)
DetectedProportionSummary = rbind(DetectedProportionSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedProportionRaw.png")
png(Filename)
plot(DetectedProportionSummary$Timestep,DetectedProportionSummary$DetectedProportion,ylim = c(0,max(DetectedProportionSummary$DetectedProportion)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Proportion of infested nodes detected", main = Title)

for(perm in 1:Nperm)
{
Sub = DetectedProportionSummary[DetectedProportionSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$DetectedProportion,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(DetectedProportionSummary$DetectedProportion, by = list(DetectedProportionSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"DetectedProportionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
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

