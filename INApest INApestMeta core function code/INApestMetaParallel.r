###########################################################################
###########################################################################
###Declares a function overlaying management on a metapopulation spread model 
###Key inputs are: 
###1) matrix of natural dispersal probibilities between each pair of sites (i.e. nodes of the network). Matrix can be non-symmetrical (i.e. can have source and sink nodes) 
###2) Envionmentally-determined per-capita propagule production
###3) Envionmentally-determined carrying capacity (K)
###4) Management parameters
### a) Annual detection probability
### b) Annual management adoption probability subsequent to detection
### c) Annual mortality probability under management
### d) Spread reduction when management adopted
###5) A binary vector of intially infested nodes OR
###   A vector of invasion risk probabilities and proportion of nodes initially invaded
###Key outputs are:
###3-dimensional arrays of invasion, management and detection status for each node in each timestep of each permuation
###2-dimensional array of invasion probability (i.e. proportion of permutations pest present) for each node in each timestep
###Line graphs summarising number of total population (as a proportion of carrying capacity) nodes infested, where infestations
###are detected, the proportion of infestations detected and nodes under management against time 
###########################################################################
###########################################################################

#######################################################################
###This version implements parallel processing in foreach for the permutation loop.
###See file "ParallelSetup.r" for notes on steps for setting up parallel processing
#######################################################################

library(abind)
library(doParallel)


INApestMetaParallel = function(
ModelName, #Name for storing results to file 
Nperm,                  #Number of permutations per parameter combination
Ntimesteps,                 #Simulation duration timesteps can be any length of time
DetectionProb,          #Per-individual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be single number or vector (nodes)
ManageProb,             #Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
ManageSD = NULL, #Option to provide standard deviation for management probability. Can be single number or vector (nodes)
MortalityProb,           #Mortality probability under management
MortalitySD = NULL, #Option to provide standard deviation for mortality probability. Can be single number or vector (nodes)
SpreadReduction,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
SpreadReductionSD = NULL, #Option to provide standard deviation for spread reduction. Can be single number or vector (nodes)
InitialPopulation = NA,        #Vector of population sizes at start of simulations
InitBioP = NA,		#Proportion of nodes infested at start of simulations
InvasionRisk = NA,           #Vector or matrix (nodes x timesteps) of probabilities of invasion from external sources
InitialInfo = NA,        #Vector or of nodes with information at start of simulations
InitInfoP = NA,		#Proportion of nodes with information at start of simulations
ExternalInfoProb = NA,           #Vector of probabilities of communication from external sources
EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
K,		       #Population carrying capacity - vector (nodes)
PropaguleProduction, #Propagules produced per individual, can be single value, vector (nodes) or matrix (nodes x years)
PropaguleEstablishment, #Propagules establishment probability. The likelihood of a dispersing propagule encountering a single
                        #host plant or establishment site within a node. Can be a ratio of search radius or patch size to node area
IncursionStartPop=NA,      #option to set population size for new incursions
SDDprob,                   #Natural disperal probability between each pair of nodes
SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
LDDprob = NA,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
			      #e.g. could be weighted by law of human visitation or data on stock movements
LDDrate = 0,         #Proportion of available propagules entering LDD
OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
OngoingExternalInfo = F,   ##Option to include ongoing communication from external sources
OutputDir = NA,		      #Directory for storing results
DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
{
###POTENTIAL ADDITIONS
###1) Make detection prob a function of population size. Could be based on individual detection prob so that DetectionProb = 1-(1-DPindividual)^N)
###   DPindividual could vary between nodes


# pre-evaluate some variables for efficiency
if(is.matrix(K) == FALSE)
{
K_is_0 <- K<=0
inv_K <- 1 / sum(K)
NodeK = K
}

if(is.matrix(PropaguleProduction) == FALSE)
 NodePropaguleProduction = PropaguleProduction

if(is.matrix(PropaguleEstablishment) == FALSE)
  NodePropaguleEstablishment = PropaguleEstablishment

if(is.matrix(EnvEstabProb) == F)
  NodeEnvEstabProb <- EnvEstabProb

if(is.matrix(Survival) == F)
  NodeSurvival <- Survival


###Declare matrix for information spread simulations
if(is.matrix(SEAM) == T)
     {
     diag(SEAM) = 0
     RandSEAM <- matrix(NA,nrow = nrow(SDDprob),ncol=nrow(SDDprob))
     }

###Assign standard deviation value to management in no value provided
if(is.null(ManageSD) == T)
	ManageSD = mean(ManageProb)/10
if(is.null(SpreadReductionSD) == T)
	SpreadReductionSD = (1-mean(SpreadReduction))/10
if(is.null(DetectionSD) == T)
	DetectionSD = mean(DetectionProb)/10
if(is.null(MortalitySD) == T)
	MortalitySD = mean(MortalityProb)/10


###########################################################
###Start of simulation
###########################################################
    
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
PermOut <- foreach(1:Nperm, .combine = 'acomb',.packages=c("abind")) %dopar% 
  {
  ###Max integer for propagule dispersal using rmultinom
  MaxInteger <- .Machine$integer.max  
  
  InvasionResultsLoop <- array(dim = c(nrow(SDDprob),Ntimesteps))
  PopulationResultsLoop <- InvasionResultsLoop
  ManagingResultsLoop <- InvasionResultsLoop
  DetectedResultsLoop <- InvasionResultsLoop

  
  ###If carrying capacity provided as matrix assign values from first timestep for population initialisation
  if(is.matrix(K) == TRUE)
    {
    K_is_0 <- K[,1]<=0
    inv_K <- 1 / sum(K[,1])
    NodeK = K[,1] 
    }      
  
  
###Assign initial infestations according either to "InitialInvasion" binary vector OR
###"InvasionRisk" probabilities and/or initial proportion of nodes infested ("InitBioP") OR
###just "InitBioP" if neither "InitialInvasion" or "InvasionRisk" supplied by user
InitBio = rep(0,times = nrow(SDDprob))
if(length(InitialPopulation) != nrow(SDDprob))
{
 if(length(InvasionRisk) == nrow(SDDprob))
        {
 	if(is.na(InitBioP) == F)
	  Infested = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitBioP),prob = InvasionRisk)
	if(is.na(InitBioP) == T)
          {
	  Infested = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
          Infested = which(Infested == 1)
	  } 
	}
 if(length(InvasionRisk) != nrow(SDDprob))
        {
 	if(is.matrix(InvasionRisk) == F)
	 Infested = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitBioP))
        if(is.matrix(InvasionRisk) == T)
          {
          Infested = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,1])
          Infested = which(Infested == 1)
	  }
	}
if(is.na(IncursionStartPop) == T) 
	InitBio[Infested] = 1
if(is.na(IncursionStartPop) == F) 
	InitBio[Infested] = IncursionStartPop
}

if(length(InitialPopulation) == nrow(SDDprob))
  InitBio = InitialPopulation

###Ensure initial population not greater than carrying capacity
InitBio[InitBio > NodeK] = NodeK[InitBio > NodeK] 

# initialise the population
N <- InitBio
if(sum(N) == 0 && OngoingExternalInvasion == F)
  warning("No initial populations and no future external invasions")

###Select nodes with information at start of simulation  according either to "InitialInfo" binary vector OR
###"ExternalInfoProb" probabilities and/or initial proportion of nodes with information ("InitInfoP") OR
###just "InitInfoP" if neither "InitialInfo" or "ExternalInfoProb" supplied by user.
###If no initial info variables provided, no nodes have info at start of simulations
InitInfo = rep(0,times = nrow(SDDprob))
if(length(InitialInfo) == nrow(SDDprob) || (is.na(InitInfoP) == F && InitInfoP>0) || is.na(sum(ExternalInfoProb)) == F )
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

###Randomly assign annual detection probability, based on mean and sd
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

###Randomly assign mortality probability when management applied
###If MortalityProb given as single value or vector (nodes)
if(is.matrix(MortalityProb)==FALSE &&(length(MortalityProb) == 1 ||length(MortalityProb) == nrow(SDDprob) ))
      {
      NodeMortalityProb = rnorm(MortalityProb,MortalitySD,n = nrow(SDDprob))
      NodeMortalityProb[NodeMortalityProb<0] = 0
      NodeMortalityProb[NodeMortalityProb>1] = 1
      }

###Populate invasion status vector ahead of timestep loop
Invaded = ifelse(InitBio>0,1,0) 

###Probability of info at start of simulation depends on
###Presence of pest and detection probability
###Select nodes that have detected infestation 
InitDetection = rbinom(1:nrow(SDDprob),size = 1,prob = 1-(1-NodeDetectionProb)^InitBio)
InitInfo[InitInfo == 0] = InitDetection[InitInfo == 0]
###Populate information status vector ahead of timestep loop
HaveInfo = InitInfo

    
  # run simulation
for (timestep in 1:Ntimesteps) 
  { 
 
  ###Allow for variation in establishment through time
  ###e.g.  climate change predictions
  ###Note: could be done outside loop, but would take heaps of memory to store 
  if(is.matrix(EnvEstabProb) == T)
    NodeEnvEstabProb <- EnvEstabProb[,timestep]
   
  if(is.matrix(Survival) == T)
    NodeSurvival <- Survival[,timestep]
    
    
  ###If carrying capacity provided as matrix assign values for relevant timestep
  if(is.matrix(K) == TRUE)
    {
    K_is_0 <- K[,timestep]<=0
    inv_K <- 1 / sum(K[,timestep])
    NodeK = K[,timestep] 
    }  

  ###If propagule production provided as matrix assign values for relevant timestep
  if(is.matrix(PropaguleProduction) == TRUE)
    NodePropaguleProduction = PropaguleProduction[,timestep] 
  
  if(is.matrix(PropaguleEstablishment) == TRUE)
    NodePropaguleEstablishment = PropaguleEstablishment[,timestep]
      
  ###Randomly assign annual detection probability, based on mean and sd
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
  
  ###Randomly assign annual mortality probability when management applied
  ###If MortalityProb given as matrix (nodes x timesteps)
  if(is.matrix(MortalityProb)==TRUE && nrow(MortalityProb) == nrow(SDDprob) && ncol(MortalityProb) == Ntimesteps)
      {
      NodeMortalityProb = rnorm(MortalityProb[,timestep],MortalitySD,n = nrow(SDDprob))
      NodeMortalityProb[NodeMortalityProb<0] = 0
      NodeMortalityProb[NodeMortalityProb>1] = 1
      }


  ###Assign management status to nodes   
  ###Management is only applied to nodes which have information
  ###i.e. where pest has been detected or following communication of information
  ###from neighbouring infested farms 
  Managing = rbinom(1:nrow(SDDprob),size = 1,prob = NodeManageProb*HaveInfo)
  
  ###Identify nodes with known extant infestations 
  Detected = Invaded*HaveInfo
  
  ###Adjust starting population for natural and managed mortality
  N0 = rbinom(nrow(SDDprob),N,NodeSurvival*(1-NodeMortalityProb*Managing))
  if(sum(N0)<=0 )
    N = N0
  Pin <-0
  Qin <- 0  
    # natural dispersal 
  if(sum(N0)>0 ) 
  {
    
  Propagules <- rpois(nrow(SDDprob), NodePropaguleProduction * N0)# propagules are produced
  
  ###natural dispersal
  ###adjust propagule spread for environmentally determined establishment probability
  ###in receiving nodes
  Pout <- Propagules*(1-LDDrate)
  if(sum(Pout)>0 && sum(Pout)<= MaxInteger) 
    Pin <- t(rmultinom(1, size=sum(Pout*rowSums(SDDprob)), prob=Pout %*% SDDprob))  # propagules are dispersed
  if(sum(Pout) > MaxInteger) 
    Pin <- colSums(sweep(SDDprob,1,Pout,`*`))
  ###human-mediated spread
  ###adjust propagule spread for environmentally determined establishment probability
  ###in receiving nodes
  if (is.matrix(LDDprob)==T) 
    {
    Qout  = Propagules*LDDrate *(1-NodeSpreadReduction*Managing)       
    if(sum(Qout)>0 && sum(Qout) < MaxInteger)  
      Qin <- t(rmultinom(1, size=sum(Qout*rowSums(LDDprob)), prob=Qout %*% LDDprob))    # propagules are dispersed
    if(sum(Qout) > MaxInteger) 
      Qin <- colSums(sweep(LDDprob,1,Qout,`*`))
    }
  
  
  #Propagule establishment probability can be viewed as a spatial process indicating the likelihood
  #that a dispersing propagule will encounter and individual site or host.
  #The simplest approach is to express the search area as a proportion of the node area
  #For a weed this could be the area of individual sites (patches) relative to node area
  #Where K is the number of available patches within the node
  #This also incorporates density dependence since propagule success depends on availability of uninfested host plants 
  #or unoccupied patches
  ###Note can produce NAs if NodeK-N0 very large and using old version of R
  N <- ifelse(K_is_0, 0, N0 + rbinom(nrow(SDDprob), NodeK-N0, 1 - exp(-NodePropaguleEstablishment*NodeEnvEstabProb*(Pin+Qin))))

  } 
 ###Update info vector for any info spread (if SEAM supplied)
 ###Note once nodes obtain info they always have info (only zero values updated)
 if(is.matrix(SEAM) == T)
  {
  RandSEAM[] <- rbinom(n=nrow(SDDprob)^2, size=1, prob = SEAM*Detected)
  InfoTransferred = ifelse(colSums(RandSEAM)>0,1,0)
  HaveInfo[HaveInfo == 0] = InfoTransferred[HaveInfo == 0]
  }
 
 
 ###Add invasion resulting from colonisation from external sources
 if(OngoingExternalInvasion == T)
  {
  if(is.matrix(InvasionRisk) == F)
   ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
  if(is.matrix(InvasionRisk) == T)
   ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,timestep])
  Invaded[Invaded == 0] = ExternalInvasion[Invaded==0]
  if(is.na(IncursionStartPop) == T) 
	N = N+ExternalInvasion
  if(is.na(IncursionStartPop) == F) 
	N = N+ExternalInvasion*IncursionStartPop
  N[N > NodeK] = NodeK[N > NodeK] 
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
 ###Update infestation vector
 Invaded = ifelse(N>0,1,0)
 
 ###Record nodes adopting management
 ManagingResultsLoop[,timestep] = Managing
  
 ###Record infested nodes
 InvasionResultsLoop[,timestep] = Invaded

 ###Record populations
 PopulationResultsLoop[,timestep] = N

 ###Select new nodes where infestation detected
 NewHaveInfo =  rbinom(1:nrow(SDDprob),size = 1,prob = 1-(1-NodeDetectionProb)^N)
 
 ###Add newly detected infestations to info vector
 ###Note once nodes obtain info they always have info (only zero values updated)
 HaveInfo[HaveInfo==0] = NewHaveInfo[HaveInfo==0]  
 
 ###Record detection status
 DetectedResultsLoop[,timestep] = HaveInfo*Invaded 
 }
 abind(InvasionResultsLoop,PopulationResultsLoop,ManagingResultsLoop,DetectedResultsLoop,along = 3)
}
stopCluster(cluster)
###########################################################
###End of Simulation
###########################################################

###Extract results for invasion, managing and dectection status into separate 3d arrays
InvasionResults <- PermOut[,,1,]
PopulationResults <- PermOut[,,2,]
ManagingResults <- PermOut[,,3,]
DetectedResults <- PermOut[,,4,]

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
saveRDS(ManagingResults, paste0(FileNameStem,"InfoLargeOut.rds"))
saveRDS(PopulationResults, paste0(FileNameStem,"PopulationLargeOut.rds"))
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))

##########################################################
###Store annual node-level invasion probs for heat maps
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

Title = ModelName


###Change in total population with time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

PopulationSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(PopulationSummary) = c("Realisation",   "Timestep",  "NodesInfested")

if(is.matrix(K) == TRUE)
    inv_K <- 1 / colSums(K)

for(perm in 1:Nperm)
{
PopulationData = PopulationResults[,,perm]
dim(PopulationData)
NodesInfested = colSums(PopulationData)*inv_K
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesInfested)
PopulationSummary = rbind(PopulationSummary,Results)
}


Filename = paste0(FileNameStem,"PopulationRaw.png")
png(Filename)
plot(PopulationSummary$Timestep,PopulationSummary$NodesInfested,ylim = c(0,1),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Total population (proportion of K)", main = Title)

for(perm in 1:Nperm)
{
Sub = PopulationSummary[PopulationSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesInfested,col  = perm)
}
dev.off()

Quantiles = as.data.frame(aggregate(PopulationSummary$NodesInfested, by = list(PopulationSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"PopulationSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,1), xlab = "Time since incursion detected (timesteps)",
ylab = "Total population (proportion of K)", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()

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
