
####################################################
###Input paramaters for INA
####################################################
###Set number of realisations
Nperm = 100
###Set simulation duration
Ntimesteps = 80

###Probability of management adoption when infestation detected
ManageProb = 0.9

###Set spread reduction under management
###Assume 1000-fold reduction in spread prob under management
SpreadReduction = 0.999

###Set probability of detection
###50% farmers confident of CNG ID 
###Default assumes those with ID knowledge have 50% chance of detection annually
###Declare vector to explore other probs
DetectionProbs = c(0.01,0.025,0.05,0.075,0.1,0.125,0.15,0.2) 

###Declare initial infestation vector
###Start with infested farms from historic data
InitBio = rep(0,times = nrow(adj))
InitBio[FarmNames %in% HistoricInfestation$farm_id] = 1 

#source("Scripts/INApestParallel.R")





####################################################
###Generate different types of social nework
####################################################
###IMPORT OR GENERATE A SOCIAL NETWORK

#OPTION 1: import social network data

#Norman's code includes the following line:
#SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread

#convert the matrix to a graph object
mtx <- SEAM
SEAM1 <-graph_from_adjacency_matrix(mtx) #convert the matrix to an iGraph object

#OPTION 2: construct a random social network (Erdős–Rényi model) based on a given probability of land-owners being connected to each other.

n <- nrow(SEAM)#number of nodes in the network, i.e. how many land-owners (nodes) you want to include
SumComm <- sum(SEAM)-nrow(SEAM)
MeanComm = SumComm/(nrow(SEAM)^2-nrow(SEAM))
  #OR
  m <- #the number of links in the network
  # in the "type", choose type = c("gnp" or "gnm") based on the choice between p and m  
  
  #example:
  SEAM2 <- erdos.renyi.game( #igraph package
    n, #number of nodes in the network, i.e. how many land-owners (nodes) you want to include
    SumComm, #either the probability for drawing a link (connection) between two random nodes, or the number of links in the network
             ####Ensure same level of connectedness as standard
    type = c("gnm"), #see  line above
    directed = TRUE,
    loops = FALSE
  )
SEAM2mat <- as.matrix(SEAM2)
sum(SEAM2mat)


#OPTION 3: small-world model for a social network
dim <- #the dimension of the starting lattice, use 1
  size <- #number of land-owners in the network  
  nei <- #the neighborhood within which the links will be connected
  p <- #the rewiring probability
  
#example: 
SEAM3 <- sample_smallworld(1, nrow(SEAM), 2, 0.1)
SEAM3mat <- as.matrix(SEAM3)
###Ensure same level of connectedness as standard
SEAM3mat <- SEAM3mat*SumComm/sum(SEAM3mat)
#sum(SEAM3mat)
##################################################
###Call INApest function 
###Looping through different detection probabilities
##################################################
for(detprob in 1:length(DetectionProbs))
{
DetectionProb = DetectionProbs[detprob]
###Option 1: Standard settings - adjency-led communication
###Declare name for storing results
ModelName = paste0("StandardSEAM","DetProb_",DetectionProb)
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(ModelName = ModelName,      #Name for storing results to disk
  Nperm = Nperm,                  #Number of permutations per paramteter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
  DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
  ManageProb = ManageProb,         
  ManageSD = NULL,
  EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
  EradicationSD = NULL,
  SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
  SpreadReductionSD = NULL,
  InitialInvasion = InitBio,        #Nodes infested at start of simulations
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  OutputDir = OutputDir	#Directory for storing results to disk	
  )

###Option2: random SEAM
###Declare name for storing results
ModelName = paste0("RandomSEAM","DetProb_",DetectionProb)
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManageProb,         
ManageSD = NULL,
EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM2mat,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 
OutputDir = OutputDir	#Directory for storing results to disk	
)

####Option 3: small world
####Declare name for storing results
ModelName = paste0("SmallWorldSEAM","DetProb_",DetectionProb)
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
  ModelName = ModelName,      #Name for storing results to disk
  Nperm = Nperm,                  #Number of permutations per parameter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = DetectionProb,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
  DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
  ManageProb = ManageProb,         
  ManageSD = NULL,
  EradicationProb = EradicationProb, #Annual probability of eradication (must be between 0 and 1)
  EradicationSD = NULL,
  SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
  SpreadReductionSD = NULL,
  InitialInvasion = InitBio,        #Nodes infested at start of simulations
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = SEAM3mat,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  OutputDir = OutputDir	#Directory for storing results to disk	
)
}
