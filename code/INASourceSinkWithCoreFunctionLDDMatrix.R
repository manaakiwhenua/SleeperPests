##################################################
###Simulate spread in sink region from starting invasions
###allocated by weighted randomisation.
###In this instance initial invasion probability
###is weighted by incursion risk from an infested region.
###Requires identification of source and sink regions (using Agribase region codes)
###and distance-based invasion probabilities between individual farms
###in source and sink regions. 
###These are obtained by running source files:
###"CrossRegion_FarmDistanceInvasionProb.R"
###"CrossRegion_FarmDistanceInvasionProbLDD.R"
###"CrossRegion_LDDWeights.R"
###Records farms infested, farms under management 
###and invasive threat to neighbouring regions
##################################################

####################################################
###Potential upgrades
####################################################

###1) could include ongoing invasion from source region in sink region simulations
###2) invasion threat from source region could be dynamic - e.g. read in simulation results from source region
###   to permit changes in number and identity of infested source farms through time.
###3) Explore matrix algebra options for deriving analytical solutions of management required for XX% erradication prob in XX years? 

library(INA)
library(tidyverse)
library(dplyr)
library(actuar)
memory.limit(size=600000) 

###########################################################################
###########################################################################
###Read in function from source file
###This simulates invasion from multiple randomised start points
###Start point probabilities can be weighted by a risk vector
###For instance, spatial or ownership links of farms in a "sink" region
###to farms in a "source region"
###########################################################################

source("INApestRandomStart.R")
source("CrossRegionInvasionThreat.R")


######################################################
######################################################
###Estimate annual erradication probability required to 
###achieve 0.95 prob of achieving erradication
###following existing model of seedbank decline
###with management (using glyphosate) to prevent seed production
###use this as management efficay in INAscene function call
######################################################
######################################################
###Seeds per metre square
Nseeds = 1500
###Seed loss from seed bank due to attrition
SeedBankLoss = 0.61
###Seed loss from seed bank due to germination
Germination = 0.064

###Model reduction in seed bank
Nyears = 50
Nseeds = vector(length = Nyears)
N0 = 1500
for(i in 1:Nyears)
 {
 Nseeds[i] = N0*(1-(SeedBankLoss+Germination))
 N0 = Nseeds[i] 
 }
###Find year where seed density <1 seed per hectare
###Consider this effective erradiation 
ErradicationYear = length(Nseeds[Nseeds>1/10^4])
 
###Annual erradication prob to achieve 0.95 errodication
###after 14 (year when Nseeds <1 per hectare) years
AnnualErradicationProb = 1-(1-0.95)^(1/ErradicationYear)

#############################################################
#############################################################

ResultsDir= paste0(main.dir,"/CrossRegionInvasionCoreFunctionLDDmatrix/",ClimateScenarios[cs],"/")
dir.create(ResultsDir,recursive = T)
Regions = c("AUCK", "EBOP", "CANT", "OTAG", "GISB", "WAIK", "HBAY","MNWG", "WELL", "MARL", "STHL", "NRTH", "TNKI",
"TASM", "NELS")
RegionResultsDir= paste0(ResultsDir,SinkRegion,"_From_",SourceRegion,"/")
dir.create(RegionResultsDir)
CrossRegionThreatDir = paste0(RegionResultsDir,"CrossRegionThreat/") 
dir.create(CrossRegionThreatDir)

########################################################
###Read in simple XY points for INA
###And Ecoclimatic index values for each farm
########################################################

data<-read_csv(paste0(ClimateScenarios[cs],"/Inputs/",SinkRegion,"_Simple_points_for_INA.csv"))

#Points for plotting in INA
geocoords<-matrix(c(data$X, data$Y), byrow=F, ncol=2)

########################################################
###Calculate distance-based dispersal probability
###and set long distance dispersal rate per farm
#############################################################

###Set mean annual number of outgoing long distance events per infested farm per year
LongDistProbPerFarm = 0.05

###Read in distance matrix
dist_mat = readRDS(paste0(ClimateScenarios[cs],"/Inputs/",SinkRegion,"_dist_mat_farm.rds"))

###Use custom distance kernel to obtain probabilities
k= 7 #shape parameter
b= 500 #scale parameter
thresh= 1/20 #1 in 20 year spread risk between adjacent farms
adj = thresh*actuar::ppareto(dist_mat, shape=k,  scale=b, lower.tail=F,log.p = F)

###Store names for nodes in network
FarmNames = row.names(adj)

######################################################
######################################################
###Read in LDD weights - based on NZ cattle movement data 
###Sourced from https://www.mpi.govt.nz/dmsdocument/49114-Analysis-of-New-Zealand-Stock-Movement-Data
###Figure 13a
###Accounts for probability of LD dispersal outside region
######################################################
######################################################

LDDweightDir = paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
LDDweight = readRDS(paste0(LDDweightDir,SinkRegion,"_",SinkRegion,"_InvasionProb_LDDweight.rds"))


##############################################################
###Assign climate based establishment prob based on ecoclimatic index
###Default function is logit
##############################################################
prob_est<-as.vector(data$Probability_Estab)
if(EI_Prob_CurveType == "SplitLinear")
	prob_est<-as.vector(data$Probability_Estab_SplitLinear)

##############################################
###Read in farm-level cross-region invasion risk
###Use to weight prob of initial infestation occurrence for each farm
###Also assign the proportion of farms invaded at the start of simulations
##############################################
CrossRegionSDDmatrix = readRDS(paste0(ClimateScenarios[cs],"/Inputs/CrossRegionDistance/",SourceRegion,"_",SinkRegion,"_InvasionProb.rds"))
CrossRegionLDDmatrix = readRDS(paste0(ClimateScenarios[cs],"/Inputs/CrossRegionDistance/",SourceRegion,"_",SinkRegion,"_InvasionProb_LDDweight.rds"))
FarmSDDRisk = rowSums(CrossRegionSDDmatrix)
FarmLDDRisk = rowSums(CrossRegionLDDmatrix)*LongDistProbPerFarm

###Farm risk weighted by cross-region incursion prob
###and climatic suitability
FarmRisk = FarmSDDRisk+FarmLDDRisk
FarmRisk = FarmRisk*prob_est
FarmRiskWeight = FarmRisk/sum(FarmRisk)

###Assign proportion of farms invaded at the start of simulations
InitBioP = 0.005

####################################################
###Input paramaters for INA
####################################################
###Set number of initial configurations
Nconfig = 10
###Set number of realisations
Nperm = 10
###Set simulation duration
Nyears = 80

###Probability of management adoption when infestation detected
ManageProb = 0.9

###Set spread reduction under management
###Assume 1000-fold reduction in spread prob under management
SpreadReduction = 0.999

###Management efficacy defined as annual erradication probability
ManEfficacy = AnnualErradicationProb

###Set probability of detection
###50% farmers confident of CNG ID 
###Default assumes those with ID knowledge have 50% chance of detection annually
###i.e. detection prob = 0.25
###Declare vector to explore other detection probabilities
DetectionProbs = c(0,0.05,0.1,0.15,0.2,0.25,0.3) 



##################################################
###Call INApest function 
###Looping through different detection probabilities
###Function incorporates management variables in filenames of outputs
##################################################

for(detprob in 1:length(DetectionProbs))
{
DetectionProb = DetectionProbs[detprob]

INApestRandomStart(
Nreals = Nconfig,       #Number of random initial invasion configurations
Nperm = Nperm,                  #Number of permutations per paramteter combination
Nyears = Nyears,                 #Simulation duration
DetectionProb = DetectionProb,          #Annual detection probability per node (e.g. farm) (must be betwee 0 and 1)
ManageProb = ManageProb,             #Annual Probability of node adopting management upon detection
AnnualErradicationProb = AnnualErradicationProb, #Annual probability of erradication (must be betwee 0 and 1)
SpreadReduction = SpreadReduction,        #Reduction in dispersal probability from nodes adopting management (must be between 0 and 1)
LDDrate = LongDistProbPerFarm,                #Mean annual per-farm distant-independent dispersal rate. Value reflects rate for fully invaded network and zero
                        #environmental limitation on establishment
InitBioP = InitBioP,			#Proportion of nodes infested at start of simulations
InvasionRisk = FarmRiskWeight,           #Vector of probabilities for weighting random assignment of initial invasion occurrences
EnvEstabProb = prob_est,           #Environmentally determined establishment probability
BPAM = adj,                   #Biophysical sdjaceny matrix - distance-based disperal probability
LDDmatrix = LDDweight*LongDistProbPerFarm,         #Option to provide long distance dispersal matrix instead of distance-independent dispesal rate
				#e.g. could be weighted by human visitation law or data on stock movements
geocoords = geocoords,              #XY points for INAscene
OutputDir = RegionResultsDir			#Directory for storing results	
)

##################################################
###Use INApest outputs to calculate 
###Invasion threat to other regions
##################################################

ThreatSource = SinkRegion
for(sink in 1:length(Regions))
if(Regions[sink] != ThreatSource)
{
ThreatSink = Regions[sink]
InputDir = paste0(main.dir,"/",ClimateScenarios[cs],"/Inputs/")
OutputDir = CrossRegionThreatDir
dir.create(OutputDir)
CrossRegionInvasionDir = paste0(main.dir,paste0("/",ClimateScenarios[cs],"/Inputs/CrossRegionDistance/"))
InvasionFileNameStem = paste0(RegionResultsDir,"DetProb_",DetectionProb,"_ErradProb_",round(AnnualErradicationProb,digits =2),
		"_SpreadReduction_",SpreadReduction,"_")
OutputFileNameStem = paste0(OutputDir,"DP_",DetectionProb,"_EP_",round(AnnualErradicationProb,digits =2),
		"_SR_",SpreadReduction,"_")
SourceInvasionProb = readRDS(paste0(InvasionFileNameStem,"InvasionProb.rds"))
CrossRegionInvasionThreat(
ThreatSource = ThreatSource, ###Source region
ThreatSink = ThreatSink, ###Sink region
LDDrate = LongDistProbPerFarm, ###annual long distance dispersal events per source farm
SourceInvasionProb = SourceInvasionProb, ###Annual farm-level invasion prob from simulations
CrossRegionInvasionDir=CrossRegionInvasionDir , ###Directory where cross-region invasion probs stored
EIDir = InputDir, ##Directory where Ecoclimatic Index data stored
OutputFileNameStem = OutputFileNameStem
)
}

##################################################
##################################################

}

##################################################
###End of detection prob loop
##################################################