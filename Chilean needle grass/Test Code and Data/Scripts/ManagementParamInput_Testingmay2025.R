
####################################################
###Input paramaters for INA
####################################################
###Set number of realisations
Nperm = 30
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
DetectionProbs = c(0,0.05,0.2) 


###Declare initial infestation vector
###Start with infested farms from historic data
InitBio = rep(0,times = nrow(adj))
InitBio[FarmNames %in% HistoricInfestation$farm_id] = 1 

######################################################
###Generate hypothetical vector (nodes) for detection, eradication, management prob and spread reduction
###Weighted by initial risk 
######################################################
StandardDetection = 0.1
StandardManage = 0.1
StandardSpreadReduction = 0.1
StandardEradication = EradicationProb

DispProb =1-(1-adj)*(1-LDDprob)
BPAM = sweep(DispProb,2,prob_est,`*`)
diag(BPAM) = 1

Risk =  sweep(BPAM,1,InitBio,`*`)
FarmRisk = colSums(Risk)
FarmRisk = FarmRisk^(1/3.1)
FarmRisk = FarmRisk/sum(FarmRisk)*length(FarmRisk)
DetectionVector = StandardDetection*FarmRisk
hist(DetectionVector)
mean(DetectionVector)

FarmRisk = colSums(Risk)
FarmRisk = FarmRisk^(1/3.1)
FarmRisk = FarmRisk/sum(FarmRisk)*length(FarmRisk)
hist(FarmRisk)
ManageVector = StandardManage*FarmRisk
hist(ManageVector)
max(ManageVector)
mean(ManageVector)

FarmRisk = colSums(Risk)
FarmRisk = FarmRisk^(1/3.1)
FarmRisk = FarmRisk/sum(FarmRisk)*length(FarmRisk)
hist(FarmRisk)
SpreadReductionVector = StandardSpreadReduction*FarmRisk
hist(SpreadReductionVector)
max(SpreadReductionVector)
mean(SpreadReductionVector)

FarmRisk = colSums(Risk)
FarmRisk = FarmRisk^(1/2)
FarmRisk = FarmRisk/sum(FarmRisk)*length(FarmRisk)
hist(FarmRisk)
EradicationVector = StandardEradication*FarmRisk
hist(EradicationVector)
max(EradicationVector)
mean(EradicationVector)

#######################################################
###Generate hypothetical climate time series with increasing 
###Prob est. Try 1% annual increase
###Would be possible to explore extreme timesteps where 
###conditions particularly favourable for establishment
#######################################################
prob_est_matrix = matrix(nrow = nrow(adj), ncol = Ntimesteps)
AnnulIncrease = 0.01
prob_est_matrix[,1] = prob_est
for(timestep in 2:Ntimesteps)
{
prob_est_matrix[,timestep] = prob_est_matrix[,(timestep-1)]*(1+AnnulIncrease*(timestep-1))
prob_est_matrix[,timestep] = ifelse(prob_est_matrix[,timestep]<1,prob_est_matrix[,timestep],1)
}
#plot(1:Ntimesteps,colMeans(prob_est_matrix)) 

#######################################################
###Generate hypothetical detection and management probability time series
###Try 1% annual increase
###May reflect awareness raising scenario
#######################################################
DetectionMatrix = matrix(nrow = nrow(adj), ncol = Ntimesteps)
AnnualIncrease = 0.01
DetectionMatrix[,1] = rnorm(nrow(adj),0.025,0.025/10)
for(timestep in 2:Ntimesteps)
{
DetectionMatrix[,timestep] = DetectionMatrix[,(timestep-1)]+AnnualIncrease
}


ManagementMatrix = matrix(nrow = nrow(adj), ncol = Ntimesteps)
AnnualIncrease = 0.01
ManagementMatrix[,1] = rnorm(nrow(adj),0.1,0.1/10)
for(timestep in 2:Ntimesteps)
{
ManagementMatrix[,timestep] = ManagementMatrix[,(timestep-1)]+AnnualIncrease
}
#plot(1:Ntimesteps,colMeans(ManagementMatrix)) 

EradicationMatrix = matrix(nrow = nrow(adj), ncol = Ntimesteps)
AnnualIncrease = 0.01
EradicationMatrix[,1] = rnorm(nrow(adj),StandardEradication,StandardEradication/10)
for(timestep in 2:Ntimesteps)
{
EradicationMatrix[,timestep] = EradicationMatrix[,(timestep-1)]+AnnualIncrease
}
#plot(1:Ntimesteps,colMeans(EradicationMatrix))

#######################################################
###Generate hypothetical spread reduction timeseries
###Try logistic growth curve with R = 0.5
###May reflect rapid increase in hygeine of movement restricton
#######################################################

N0 = 10
Nt1 = N0
K = 100
R = 0.5
SR = Nt1/100
for(timesteps in 2:80)
{
Nt2 = Nt1+R*Nt1*((K-Nt1)/K)
SR = c(SR,Nt2/100)
Nt1 = Nt2
}
#plot(1:Ntimesteps,SR)
SpreadReductionMatrix = matrix(nrow = nrow(adj), ncol = Ntimesteps)
AnnualIncrease = 0.01
SpreadReductionMatrix[,1] = SR[1]
for(timestep in 2:Ntimesteps)
{
SpreadReductionMatrix[,timestep] = SR[timestep]
}
#plot(1:Ntimesteps,colMeans(SpreadReductionMatrix)) 


##################################################
###Call INApest function 
###Comparing different trials with a standad model
##################################################

###Standard settings
###Declare name for storing results
ModelName = paste0("Standard")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
EradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)

###matrix (nodes x timesteps) of establishment probs increasing with time
ModelName = paste0("EstabMatrixStandrardMan")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
  ModelName = ModelName,      #Name for storing results to disk
  Nperm = Nperm,                  #Number of permutations per paramteter combination
  Ntimesteps = Ntimesteps,                 #Simulation duration
  DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
  DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
  ManageProb = StandardManage,         
  ManageSD = NULL,
  EradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
  EradicationSD = NULL,
  SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
  SpreadReductionSD = NULL,
  InitialInvasion = InitBio,        #Nodes infested at start of simulations
  EnvEstabProb = prob_est,           #Environmentally determined establishment probability
  SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
  SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
  LDDprob = LDDprob,         #Long distance dispersal probability matrix 
  
  OutputDir = OutputDir	#Directory for storing results to disk	
)

###Vector (nodes) of annual eradication prob weighted by inital risk
###Declare name for storing results
ModelName = paste0("EradicationVector")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
EradicationProb = EradicationVector, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)




###Vector (nodes) of detection prob weighted by inital risk
ModelName = paste0("DetectionVector")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionVector,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
EradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)

###Vector (nodes) of manage prob weighted by inital risk
ModelName = paste0("ManageVector")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManageVector,         
ManageSD = NULL,
EradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)

###Vector (nodes) of spread reduction weighted by inital risk
ModelName = paste0("SpreadReductionVector")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
EradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = SpreadReductionVector,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)


###Vector (nodes) of detection, management, eradication and spread reduction weighted by inital risk
###Prioritises all spatially varying managemet variables by risk
ModelName = paste0("AllManagementPrioritised")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionVector,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManageVector,         
ManageSD = NULL,
EradicationProb = EradicationVector, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = SpreadReductionVector,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)


###Matrix (nodes x timesteps) of annual eradication prob randomly allocated to nodes, 
###increasing over timeme
###Declare name for storing results
ModelName = paste0("EradicationMatrix")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
EradicationProb = EradicationMatrix, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)



###Matrix (nodes x timesteps) of annually increasing detection probs 
ModelName = paste0("DetectionMatrix")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionMatrix,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
EradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)

###Matrix (nodes x timesteps) of annually increasing management probs
###Declare name for storing results
ModelName = paste0("ManagementMatrix")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManagementMatrix,         
ManageSD = NULL,
EradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = StandardSpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)

###Spread reduction matrix
###Declare name for storing results
ModelName = paste0("SpreadReductionMatrix")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = StandardDetection,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = StandardManage,         
ManageSD = NULL,
EradicationProb = StandardEradication, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = SpreadReductionMatrix,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)

###All management variables improving over time
###Declare name for storing results
ModelName = paste0("AllManagementIncreasing")
OutputDir = paste0(RegionResultsDir,ModelName,"/")
dir.create(OutputDir,showWarnings = F)
INApestParallel(
ModelName = ModelName,      #Name for storing results to disk
Nperm = Nperm,                  #Number of permutations per paramteter combination
Ntimesteps = Ntimesteps,                 #Simulation duration
DetectionProb = DetectionMatrix,   #Annual detection probability per node (e.g. farm) (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability can be vector (nodes)
ManageProb = ManagementMatrix,         
ManageSD = NULL,
EradicationProb = EradicationMatrix, #Annual probability of eradication (must be between 0 and 1)
EradicationSD = NULL,
SpreadReduction = SpreadReductionMatrix,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
SpreadReductionSD = NULL,
InitialInvasion = InitBio,        #Nodes infested at start of simulations
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
SDDprob = adj,                   #Biophysical adjaceny matrix - distance-based disperal probability
SEAM = SEAM,                  ##Socioeconomic adjaceny matrix - info transfer probability between farms
LDDprob = LDDprob,         #Long distance dispersal probability matrix 

OutputDir = OutputDir	#Directory for storing results to disk	
)


