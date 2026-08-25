# Regression checks for pathogen-detection-triggered information
source(file.path("..","src","INApestPathogen.R"))
source(file.path("..","src","INApest.R"))

# Binary: certain pathogen detection should create information; disabled trigger should not.
A <- diag(2)
p_on <- INApestPathogen(Model="Binary", InitialPresent=c(1,0), DetectionProb=1, DetectionTriggersInfo=TRUE)
p_off <- INApestPathogen(Model="Binary", InitialPresent=c(1,0), DetectionProb=1, DetectionTriggersInfo=FALSE)
common <- list(ModelName="x",Nperm=1,Ntimesteps=1,DetectionProb=0,DetectionSD=0,ManageProb=1,ManageSD=0,
               EradicationProb=0,EradicationSD=0,SpreadReduction=0,SpreadReductionSD=0,
               InitialInvasion=c(1,1),InitialInfo=c(0,0),InfoRetentionProb=1,EnvEstabProb=1,Survival=1,
               SDDprob=A,LDDprob=0,SEAM=0,SaveResults=FALSE,DoPlots=FALSE,Seed=1)
z_on <- do.call(INApest,c(common,list(Pathogen=p_on)))
z_off <- do.call(INApest,c(common,list(Pathogen=p_off)))
stopifnot(z_on$PathogenDetectedResults[1,1,1] == 1)
stopifnot(z_off$PathogenDetectedResults[1,1,1] == 1)
# Trigger affects information/management timing, while pathogen detection itself is recorded in both cases.
stopifnot(z_on$DetectedResults[1,1,1] == 1)
stopifnot(z_off$DetectedResults[1,1,1] == 0)

# Detection is impossible without infectious/present pathogen state.
stopifnot(z_on$PathogenDetectedResults[2,1,1] == 0)

cat("PASS: pathogen detection information-trigger checks\n")
