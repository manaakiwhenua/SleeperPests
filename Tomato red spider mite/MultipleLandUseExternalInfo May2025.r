#####################################################################################
###Code testing inclusion of pre-emptive management in nodes adjacent to detected infestations
###And ongoing risk of incursions at the border
#####################################################################################


#--------------------------------------------------------------------------------------------------------------------
# Tomato red spider mite metapopulation model
# by John Kean

# version 2023-05-5

 
#--------------------------------------------------------------------------------------------------------------------
# This model simulates tomato red spider mite (Tetranychus evansi) growth and spread on a landscape grid.
#
# Two dispersal functions are implemented: a distance-based exponential power kernel,
# (NB based on distance between cell centroids - I have explored the errors associate with this but ignore them here)
# and human population-based gravity model for human-vectored dispersal.
# All processes are integer-based, so stochastic.

# The input .csv file should have rows for each subpopulation and columns for:
#     Id: (integer) a unique identifier for each subpopulation
#      x: (real) the x-coordinate (in the same units as the y-coordinate)
#      y: (real) the y-coordinate (in the same units as the x-coordinate)
#  hosts: (integer) the estimated number of host plants available in each location
# humans: (real) the local human population size, used to drive a gravity dispersal component
#     N0: (integer) the starting population size
#     EI: (real) whether the population can overwinter: 0 if not, >0 if can (so can use CLIMEX EI)
#   gens: (real) generations per timestep at each site, estimated from a CLIMEX model


#--------------------------------------------------------------------------------------------------------------------
# preliminaries

# clear the memory
#rm(list = ls())

# load required libraries
#library(ggplot2)
library (gnorm)  # generalised normal distribution aka exponential power function



#--------------------------------------------------------------------------------------------------------------------
# define parameters and variables

# set csv file names
input_file   <- "Tomato red spider mite data all NZ.csv"
ResultsDir = "TestingMay2025/"
dir.create(ResultsDir)


# parameters that convert generations per timestep from CLIMEX into rate of increase
ROI_a <- -0.126
ROI_b <- 3.242
ROI_c <- -9.767

# natural propagule production rate
alpha <- 0.001

# natural dispersal kernel parameters
disp_scale <- 1.5   # distance scalar (km)
disp_shape <- 2     # <1: fat-tailed, 1: exponential, 2: Gaussian, 3: Ribbens, 10: ~uniform

# human-assisted dispersal parameters
human_prop <- 1e-6  # propagule pressure coefficient
human_acti <- 2.4   # human activity exponent
human_dist <- 2.5   # distance exponent - should be approx 2

# propagule establishment parameter
# this assumes that dispersing propagules can search an area 
# of 50m x 50m when seeking a host plant
estab <- 0.0001


#--------------------------------------------------------------------------------------------------------------------
# load variables from file and calculate various local parameters

# read csv file into data frame
d <- read.csv(input_file)

###Use part of dataset for faster testing
#d<-d[1:1000,]
#d$N0[1] = 10

if(DoClimateChange == TRUE)
{
# use climate change scenarios
d$EI <- d$EI_2090
d$gens <- d$gens_2090
}

# find dimensions of x and y
x_dist <- min(diff(sort(unique(d$x))))
y_dist <- min(diff(sort(unique(d$y))))

# convert generations per timestep into rate of increase
d$R0 <- exp(ROI_a*d$gens*d$gens + ROI_b*d$gens + ROI_c)

# generate carrying capacity depending on host availability and climate suitability
K <- ifelse(d$hosts<=0, 0, d$hosts)
K <- ifelse(d$EI<=0, 0, K)






#--------------------------------------------------------------------------------------------------------------------
# calculate distances between locations

# function distance_matrix
#   Returns a symmetrical matrix of distances, where the diagonals = 0
#   x and y (vectors) are the x and y coordinates of the grid location centres
#   x_dist and y_dist (real) are the size of location cells in the x and y directions
#   For this type of data, this is ~10x faster than R's dist() function

distance_matrix <- function (x, y, x_dist=0, y_dist=0) {
  
  # find dimensions of x and y if not supplied
  if (x_dist<=0) x_dist <- min(diff(sort(unique(x))))
  if (y_dist<=0) y_dist <- min(diff(sort(unique(y))))
  
  # find dimensions of x and y and invert for speed
  x_rate <- 1/x_dist
  y_rate <- 1/y_dist
  
  # find array size
  n <- length(x)
  
  # calculate result and return it
  Dist <- matrix(sapply(1:n, function(i) { sqrt(((x[i]-x)*x_rate)^2 + ((y[i]-y)*y_rate)^2) }), nrow=n, ncol=n)
  return (Dist)
  
}


# create a matrix of distances between locations
Dist <- distance_matrix(d$x, d$y, x_dist, y_dist)

#######################################################
###Estimate external invasion risk as product of 
###annual border incursion rate and proportion of total human population
#######################################################

PropTotalHumans = d$"humans"/sum(d$"humans")




########################################################
###Define a maximum threshold distance for informaton spread
###from known infestations
###and communication rate (e.g. by authorities and owners of infested farms) to farms within 
###the threshold distance
###Set socioeconomic adjacency matrix (SEAM) so that nodes within threshold distance
###are assigned a fixed probability of info transfer equal to the communication rate
###set all other SEAM values to zero
###For TRSM example scenario may reflect active surveilance in grid cells
###neighbouring known infestations 
#############################################################

InfoMaxDist = 1.5 #Info transfer threshold distance in grid square units
			#This transfers info to the 8 squares surrounding a detected infestation	
                  #Can be set to maximum annual (non-human aided) dispersal distance of pest
CommunicationRate = 0.9 #Proportion of nodes within threshold distance contacted each timestep
SEAM = Dist
SEAM = ifelse(Dist < InfoMaxDist,CommunicationRate,0)
diag(SEAM) = 0

#--------------------------------------------------------------------------------------------------------------------
# set up natural dispersal transition matrix

# function exp_power_dispersal_matrix
#   Calculates a dispersal transition matrix based on an exponential power function kernel
#   See e.g. Bullock et al 2017. J Ecol 105: 6-19. https://doi.org/10.1111/1365-2745.12666.
#   Returns a transition matrix
#   x and y (vectors) are the x and y coordinates of the grid location centres
#   dist (matrix) is the matrix of distances between locations
#   scale (real, >0) approximates the "average" dispersal distance
#   shape (real, >0) is the shape of the tail; <1: fat-tailed, 1: exponential, 2: Gaussian, 10: ~uniform
#   x_dist and y_dist (real) are the size of location cells in the x and y directions
#   remain (boolean) determines whether to allow dispersers to settle back in their source subpopulation
#   contain (boolean) determines whether to contain dispersers within the subpopulations; F allows loss into the sea etc
#   plot (boolean) is whether to plot the dispersal kernel

exp_power_dispersal_matrix <- function (x, y, Dist, scale=1, shape=2, remain=T, contain=F, plot=F) {
  
  # plot dispersal kernel
  if (plot) curve(exp(-(x/scale)^shape), 0, scale*3, xlab="Natural dispersal distance", ylab="Relative probability")
  
  # apply dispersal kernel
  over_scale <- 1/scale
  W <- exp(-(Dist*over_scale)^shape)
  
  # if dispersers must leave their source location
  if (!remain) {
    diag(W) <- 0
  }
  
  # normalise the matrix
  if (contain) {
    W <- W / rowSums(W)
  } else {
    W <- W / max(rowSums(W))  # NB assumes at least one location is well surrounded
  }
  W[is.nan(W)] <- 0         
  
  return (W)
}


# set up natural dispersal transition matrix
if (alpha > 0) {
  nd <- exp_power_dispersal_matrix(d$x, d$y, Dist, disp_scale, disp_shape, remain=T, contain=F, plot=T)
} else {
  nd <- NULL
}




#--------------------------------------------------------------------------------------------------------------------
# set up a human-vectored dispersal transition matrix

# function gravity_dispersal_matrix
#   Calculates a dispersal transition matrix based on a gravity model
#   See Maino et al. 2021. J Appl Ecol 58: 789-800. https://doi.org/10.1111/1365-2664.13812.
#   Returns a transition matrix
#   dist (matrix) is the matrix of distances between locations
#   pull (vector) is the gravity component e.g. human population
#   ge (real) is the gravity exponent
#   de (real) is the distance exponent

gravity_dispersal_matrix <- function (Dist, pull=0, ge=0, de=0) {
  
  # number of locations
  dim <- nrow(Dist)
  
  # calculate matrix
  G <- matrix(sapply(1:dim, function(i) { (pull[i] * pull)^ge }), nrow=dim, ncol=dim)
  G <- G / (Dist^de)
  diag(G) <- 0
  
  # normalise the matrix
  G <- G / rowSums(G)
  G[is.na(G)] <- 0       # occurs for locations where pull=0
  
  # fix diagonals
  diag(G) <- diag(G) + 1 - rowSums(G)
  
  return (G)                               
}
hd <- NULL
if(DoHumanSpread == T)
{
# set up human-vectored dispersal transition matrix
hd <- gravity_dispersal_matrix(Dist, d$humans, human_acti, human_dist)
}

# remove distance matrix to free up resources
rm(Dist)


############################################################################
###Provide K as matrix of nodes x land uses
###Just use artificial example to test new functionality
############################################################################
Planduses = c(0.5,0.5)
Klanduse = matrix(nrow = length(K),ncol = 2)
Klanduse[,1] = floor(K*Planduses[1])
Klanduse[,2] = ceiling(K*Planduses[2])
N0landuse  = cbind(floor(d$N0*Planduses[1]),ceiling(d$N0*Planduses[2]))

AnnualIncursionRate = 0.1
H_vectors <- pmin(1, human_prop * d$humans)
DetectionProbs = c(0,0.15,0.2,0.3,0.5,0.75,0.9)
detprob = length(DetectionProbs)


#ModelName = paste0("StandardCurrentClim_DetProb_Parallel",DetectionProbs[detprob])
#if(DoClimateChange == TRUE)
#  ModelName = paste0("StandardFutureClim_DetProb_Parallel",DetectionProbs[detprob])
Ntimesteps = 50
Nperm = 50
ModelName = "MultipleLandUseTest"
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallelMultipleLandUse(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  Nlanduses = 2,
  DetectionProb = c(0,0.5),  #Annual per-individual detection Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  ManageProb = c(0.5,0.99),             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = c(0.5,0.99),           #Annual mortality probability under management
  SpreadReduction = c(0,0.9),        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = N0landuse,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  #InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  #InitialInfo = sample(1:nrow(d),10),
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = Klanduse,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
ParallelTime <- End-Start


ModelName = "MultipleLandUseTestZeroMan"
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()

INApestMetaParallelMultipleLandUse(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  Nlanduses = 2,
  DetectionProb = c(0,0),  #Annual per-individual detection Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  ManageProb = c(0.5,0.99),             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = c(0.5,0.99),           #Annual mortality probability under management
  SpreadReduction = c(0,0.9),        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = N0landuse,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  #InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  #InitialInfo = sample(1:nrow(d),10),
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = Klanduse,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
ParallelTime <- End-Start


################################
###Detection prob 0.9 no external info
################################
ModelName = paste0("StandardCurrentClim_DetProb_Parallel","0.9")
if(DoClimateChange == TRUE)
  ModelName = paste0("StandardFutureClim_DetProb_Parallel","0.9")
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = 0.9,  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  #InitialInfo = sample(1:nrow(d),10),
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
ParallelTime <- End-Start

################################
###Initial info as binary vector
################################
InitialInfo = rep(0,times = nrow(d))
Sample <- sample(1:nrow(d),10)
InitialInfo[Sample] = 1
ModelName = paste0("StandardCurrentClim_DetProb_ParallelInitialInfo",DetectionProbs[detprob])
if(DoClimateChange == TRUE)
ModelName = paste0("StandardFutureClim_DetProb_ParallelInitialInfo",DetectionProbs[detprob])
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
ModelName = ModelName,
Nperm = Nperm,                  #Number of permutations per parameter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionProbs[detprob],  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
MortalityProb = 0.99,           #Annual mortality probability under management
SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
InitialPopulation = d$N0,        #Population size at start if simulations
InitBioP = NA,		#Proportion of nodes infested at start of simulations
InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
InitialInfo = InitialInfo, #
EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
K = K,		       #Population carrying capacity
PropaguleProduction = alpha*d$R0, #Propagules produced per individual
PropaguleEstablishment = estab, #Propagules establishment rate
IncursionStartPop=10,      #option to set population size for new incursions
SDDprob = nd,                   #Natural disperal probability between each pair of nodes
SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
			      #e.g. could be weighted by law of human visitation or data on stock movements
LDDrate = H_vectors,         #Proportion of available propagules entering LDD
OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
OngoingExternalInfo = F,   ##Option to include ongoing invasion from external sources
OutputDir = OutputDir,		      #Directory for storing results
DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
ParallelTime <- End-Start

################################
###Initial info as proportion of nodes
################################
InitInfoP = 10/nrow(d)
ModelName = paste0("StandardCurrentClim_DetProb_ParallelInitInfoP",DetectionProbs[detprob])
if(DoClimateChange == TRUE)
  ModelName = paste0("StandardFutureClim_DetProb_ParallelInitInfoP",DetectionProbs[detprob])
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProbs[detprob],  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  InitInfoP = InitInfoP, #
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OngoingExternalInfo = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
ParallelTime <- End-Start

################################
###Initial info as vector of probabilities
################################
ExternalInfoProb = rep(10/nrow(d),times = nrow(d))
ModelName = paste0("StandardCurrentClim_DetProb_ParallelExternalInfoProb",DetectionProbs[detprob])
if(DoClimateChange == TRUE)
  ModelName = paste0("StandardFutureClim_DetProb_ParallelExternalInfoProb",DetectionProbs[detprob])
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProbs[detprob],  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  ExternalInfoProb = ExternalInfoProb, #
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OngoingExternalInfo = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
ParallelTime <- End-Start


################################
###Initial info as vector of probabilities
###and proportion of nodes
################################
ModelName = paste0("StandardCurrentClim_DetProb_ParallelExternalInfoProbInitInfoP",DetectionProbs[detprob])
if(DoClimateChange == TRUE)
  ModelName = paste0("StandardFutureClim_DetProb_ParallelExternalInfoProbInitInfoP",DetectionProbs[detprob])
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProbs[detprob],  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  InitInfoP = 10/nrow(d), #
  ExternalInfoProb = ExternalInfoProb, #
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OngoingExternalInfo = F,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
ParallelTime <- End-Start

################################
###Initial info as vector of probabilities
###with ongoing external communication
################################
ModelName = paste0("StandardCurrentClim_DetProb_ParallelOngoingExternalInfo",DetectionProbs[detprob])
if(DoClimateChange == TRUE)
  ModelName = paste0("StandardFutureClim_DetProb_ParallelOngoingExternalInfo",DetectionProbs[detprob])
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProbs[detprob],  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  ExternalInfoProb = ExternalInfoProb, #
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OngoingExternalInfo = T,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
ParallelTime <- End-Start

################################
###Initial info as matrix of probabilities
###with ongoing external communication
################################
ExternalInfoProbMatrix = matrix(ncol = Ntimesteps, nrow  = nrow(d))
for(i in 1:20)
  ExternalInfoProbMatrix[,i] = ExternalInfoProb
for(i in 21:Ntimesteps)
  ExternalInfoProbMatrix[,i] = ExternalInfoProb*10
ModelName = paste0("StandardCurrentClim_DetProb_ParallelOngoingExternalInfoMatrix",DetectionProbs[detprob])
if(DoClimateChange == TRUE)
  ModelName = paste0("StandardFutureClim_DetProb_ParallelOngoingExternalInfoMatrix",DetectionProbs[detprob])
OutputDir = paste0(ResultsDir,ModelName,"/")
dir.create(OutputDir)
Start <- Sys.time()
INApestMetaParallel(
  ModelName = ModelName,
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProbs[detprob],  #Annual detection probability or vector of probabilties per node (e.g. farm) (must be between 0 and 1)
  ManageProb = 0.99,             #Annual Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
  MortalityProb = 0.99,           #Annual mortality probability under management
  SpreadReduction = 0.9,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
  InitialPopulation = d$N0,        #Population size at start if simulations
  InitBioP = NA,		#Proportion of nodes infested at start of simulations
  InvasionRisk = PropTotalHumans*AnnualIncursionRate,  #Vector of probabilities of external invasion risk
  ExternalInfoProb = ExternalInfoProbMatrix, #
  EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
  Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
  K = K,		       #Population carrying capacity
  PropaguleProduction = alpha*d$R0, #Propagules produced per individual
  PropaguleEstablishment = estab, #Propagules establishment rate
  IncursionStartPop=10,      #option to set population size for new incursions
  SDDprob = nd,                   #Natural disperal probability between each pair of nodes
  SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
  LDDprob = hd,         #Option to provide long distance (human-mediated) dispersal matrix instead of distance-independent dispesal rate
  #e.g. could be weighted by law of human visitation or data on stock movements
  LDDrate = H_vectors,         #Proportion of available propagules entering LDD
  OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
  OngoingExternalInfo = T,   ##Option to include ongoing invasion from external sources
  OutputDir = OutputDir,		      #Directory for storing results
  DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
End <- Sys.time()
ParallelTime <- End-Start