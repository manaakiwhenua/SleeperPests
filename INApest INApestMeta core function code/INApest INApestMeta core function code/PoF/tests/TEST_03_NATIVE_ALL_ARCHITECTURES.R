options(stringsAsFactors = FALSE)

script_path <- function() {
  z <- grep('^--file=', commandArgs(FALSE), value = TRUE)
  if (length(z)) return(normalizePath(sub('^--file=', '', z[1L]), winslash='/', mustWork=TRUE))
  if (!is.null(sys.frame(1)$ofile)) return(normalizePath(sys.frame(1)$ofile, winslash='/', mustWork=TRUE))
  stop('Could not determine test-script path.')
}
test_dir <- dirname(script_path())
root <- normalizePath(file.path(test_dir, '..'), winslash='/', mustWork=TRUE)
engine_dir <- file.path(root, 'engine')
source(file.path(root, 'src', 'INApestProofOfFreedom.R'))

pass <- function(label) cat('PASS -', label, '\n')
check <- function(ok, label, msg='condition was not true') {
  if (!isTRUE(ok)) stop('FAIL - ', label, ': ', msg, call.=FALSE)
  pass(label)
}
close <- function(x,y,tol=1e-12) isTRUE(all.equal(as.numeric(x),as.numeric(y),tolerance=tol,check.attributes=FALSE))

load_engine <- function(files, point=FALSE) {
  e <- new.env(parent=.GlobalEnv)
  sys.source(file.path(engine_dir, 'INApestBiocontrol.R'), envir=e)
  if (point) sys.source(file.path(engine_dir, 'INApestBiocontrolPointAdapter.R'), envir=e)
  for (f in files) sys.source(file.path(engine_dir, f), envir=e)
  e
}

check_pof_output <- function(result, architecture, expected_nperm, expected_nt=2L,
                             expected_no_detection=NULL) {
  fs <- INApestBiocontrolFreedomState(result, 'Q')
  check(fs$Nparticles == expected_nperm, paste(architecture, 'PoF particle count'))
  check(fs$Ntimesteps == expected_nt, paste(architecture, 'PoF timestep count'))
  check(all(!fs$Freedom), paste(architecture, 'Q present therefore not free'))
  surv <- INApestBiocontrolSurveillance(DetectionProb=0.5, DetectStages='adult', DetectAgents='Q')
  L <- INApestBiocontrolObservationLikelihood(result, surv, 1L, 0L, 'Q')
  if (!is.null(expected_no_detection))
    check(close(L, rep(expected_no_detection, expected_nperm)), paste(architecture, 'Q no-detection likelihood'))
  fit <- INApestBiocontrolPoF(result, ObservationHistory=c(0,NA), Surveillance=surv, BiocontrolAgents='Q')
  check(close(fit$Summary$PosteriorPoF_Biocontrol[1L], 0), paste(architecture, 'posterior Q PoF'))
  invisible(TRUE)
}

make_node_bc <- function(e, target_stage=NULL, Ntimesteps=2L) {
  release <- array(0L, dim=c(2L,2L,Ntimesteps), dimnames=list(NULL,c('juvenile','adult'),NULL))
  release[1L,'adult',1L] <- 2L
  movement <- matrix(c(0,1,0,1), nrow=2L, byrow=TRUE)
  agent <- e$INApestBiocontrolAgent(
    Name='Q', Stages=c('juvenile','adult'), InitialState=0,
    Release=release, ReleaseStage='adult', Transition=diag(2),
    Movement=movement, MovementStages=c('juvenile','adult'),
    AttackStage='adult', TargetStage=target_stage, AttackRate=1e6,
    RecruitStage='juvenile', OffspringPerAttack=1L)
  e$INApestBiocontrol(agent)
}

make_point_bc <- function(e, target_stage=NULL) {
  init <- matrix(c(0L,1L), nrow=1L, dimnames=list(NULL,c('juvenile','adult')))
  agent <- e$INApestBiocontrolAgent(
    Name='Q', Stages=c('juvenile','adult'), InitialState=init,
    Release=NULL, Transition=diag(2), Movement=NULL,
    AttackStage='adult', TargetStage=target_stage, AttackRate=0,
    RecruitStage='juvenile', OffspringPerAttack=1L)
  e$INApestBiocontrol(agent)
}

cat('INApest biocontrol-target PoF native architecture validation\n')
cat('R:', R.version.string, '\n')
cat('Platform:', R.version$platform, '\n\n')

SDD0 <- matrix(0,2,2)
A3 <- diag(3)
initial_tm <- matrix(c(0L,0L,0L, 2L,5L,1L), 2L,3L,byrow=TRUE)

# 1. Binary ------------------------------------------------------------------
e <- load_engine('INApest.R')
bc <- make_node_bc(e)
r <- e$INApest(ModelName='POF_BC_BINARY',Nperm=2L,Ntimesteps=2L,
  DetectionProb=0,ManageProb=0,EradicationProb=0,SpreadReduction=0,
  InitialInvasion=c(0L,1L),InitialInfo=c(0L,0L),ExternalInfoProb=0,
  EnvEstabProb=1,Survival=1,SDDprob=SDD0,LDDprob=0,
  OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,Pathogen=NULL,
  Biocontrol=bc,SaveResults=FALSE,DoPlots=FALSE,Seed=20260924L)
check_pof_output(r,'Binary',2L,expected_no_detection=0.25)

# 2. Meta --------------------------------------------------------------------
e <- load_engine('INApestMeta.r')
bc <- make_node_bc(e)
set.seed(20260924L)
r <- e$INApestMeta(ModelName='POF_BC_META',Nperm=2L,Ntimesteps=2L,
  DetectionProb=0,ManageProb=0,MortalityProb=0,FecundityReduction=0,SpreadReduction=0,
  InitialPopulation=c(0L,5L),InitialInfo=c(0L,0L),ExternalInfoProb=0,
  EnvEstabProb=1,Survival=1,K=c(10L,10L),PropaguleProduction=0,
  PropaguleEstablishment=0,SDDprob=SDD0,OngoingExternalInvasion=FALSE,
  OngoingExternalInfo=FALSE,Pathogen=NULL,Biocontrol=bc,OutputDir=tempdir(),
  DoPlots=FALSE,ReturnResults=TRUE)
check_pof_output(r,'Meta',2L,expected_no_detection=0.25)

# 3. Multiple Land Use -------------------------------------------------------
e <- load_engine('INApestMetaMultipleLandUse.r')
bc <- make_node_bc(e)
initial_mlu <- matrix(c(0L,0L,2L,3L),2L,2L,byrow=TRUE)
set.seed(20260924L)
r <- e$INApestMetaMultipleLandUse(ModelName='POF_BC_MLU',Nperm=2L,Ntimesteps=2L,Nlanduses=2L,
  DetectionProb=0,ManageProb=0,MortalityProb=0,FecundityReduction=0,SpreadReduction=0,
  InitialPopulation=initial_mlu,InitialInfo=c(0L,0L),ExternalInfoProb=0,EnvEstabProb=1,
  Survival=1,K=matrix(10L,2L,2L),PropaguleProduction=0,PropaguleEstablishment=0,
  SDDprob=SDD0,OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,Pathogen=NULL,
  Biocontrol=bc,OutputDir=tempdir(),DoPlots=FALSE,ReturnResults=TRUE)
check_pof_output(r,'Multiple Land Use',2L,expected_no_detection=0.25)

# 4. Transition Matrix -------------------------------------------------------
e <- load_engine('INApestMetaTransitionMatrix.r')
bc <- make_node_bc(e,target_stage=2L)
set.seed(20260924L)
r <- e$INApestMetaTransitionMatrix(ModelName='POF_BC_TM',Nperm=2L,Ntimesteps=2L,
  Nstages=3L,Weights=c(1,1,1),Transition=A3,DetectionProb=0,ManageProb=0,
  MortalityProb=0,FecundityReduction=0,SpreadReduction=0,InitialPopulation=initial_tm,
  InitialInfo=c(0L,0L),ExternalInfoProb=0,EnvEstabProb=1,K=c(20L,20L),
  SeedbankK=c(20L,20L),PropaguleEstablishment=0,SDDprob=SDD0,
  OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,Pathogen=NULL,Biocontrol=bc,
  OutputDir=tempdir(),DoPlots=FALSE,ReturnResults=TRUE)
check_pof_output(r,'Transition Matrix',2L,expected_no_detection=0.25)

# 5. Vertebrate Node ---------------------------------------------------------
e <- load_engine('INApestVertebrateNode.R')
bc <- make_node_bc(e,target_stage=2L)
node_interaction <- function(population,timestep,...) { population[1L,1L] <- population[1L,1L] + 1L; population }
node_vertebrate <- list(Interaction=list(Update=node_interaction))
set.seed(20260924L)
r <- e$INApestVertebrateNode(ModelName='POF_BC_VERTEBRATE_NODE',Nperm=2L,Ntimesteps=2L,
  Nstages=3L,Weights=c(1,1,1),Transition=A3,DetectionProb=0,ManageProb=0,
  MortalityProb=0,FecundityReduction=0,SpreadReduction=0,InitialPopulation=initial_tm,
  InitialInfo=c(0L,0L),ExternalInfoProb=0,EnvEstabProb=1,K=c(20L,20L),
  SeedbankK=c(20L,20L),PropaguleEstablishment=0,SDDprob=SDD0,
  OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,Vertebrate=node_vertebrate,
  Biocontrol=bc,OutputDir=tempdir(),DoPlots=FALSE,SaveResults=FALSE,DoProgress=FALSE)
check_pof_output(r,'Vertebrate Node',2L,expected_no_detection=0.25)

# Shared point settings ------------------------------------------------------
zero_kernel <- function(n,parents=NULL,timestep=NULL,perm=NULL) data.frame(dx=rep(0,n),dy=rep(0,n))
initial_points_meta <- data.frame(x=0.5,y=0.5)
initial_points_tm <- data.frame(x=0.5,y=0.5,stage=1L)

# 6. MetaPoint ---------------------------------------------------------------
e <- load_engine('INApestMetaPoint_Biocontrol_v0.1.R',point=TRUE)
support <- e$INApestBiocontrolPointSupport(xmin=0,xmax=1,ymin=0,ymax=1,nrow=1L,ncol=1L)
bc <- make_point_bc(e)
r <- e$INApestMetaPoint(ModelName='POF_BC_METAPOINT',Nperm=2L,Ntimesteps=2L,
  InitialPoints=initial_points_meta,Pathogen=NULL,Survival=1,PropaguleProduction=0,
  PropaguleEstablishment=1,EnvEstabProb=1,SDDkernel=zero_kernel,LDDrate=0,
  LocalK=Inf,DetectionProb=0,ManageProb=0,MortalityProb=0,FecundityReduction=0,
  SpreadReduction=0,InitialInfo=NULL,InfoRadius=0,InfoTransferProb=0,
  ExternalInfoProb=0,OngoingExternalInfo=FALSE,OngoingExternalInvasion=FALSE,
  SaveResults=FALSE,DoProgress=FALSE,Seed=20260924L,Biocontrol=bc,
  BiocontrolPointSupport=support)
check_pof_output(r,'MetaPoint',2L,expected_no_detection=0.5)

# 7. Point Transition Matrix -------------------------------------------------
e <- load_engine('INApestPointTransitionMatrix_Biocontrol_v0.1.R',point=TRUE)
support <- e$INApestBiocontrolPointSupport(xmin=0,xmax=1,ymin=0,ymax=1,nrow=1L,ncol=1L)
bc <- make_point_bc(e,target_stage=1L)
r <- e$INApestPointTransitionMatrix(ModelName='POF_BC_POINT_TM',Nperm=2L,Ntimesteps=2L,
  Nstages=2L,Weights=c(1,1),Transition=diag(2),InitialPoints=initial_points_tm,
  SDDkernel=zero_kernel,LDDkernel=NULL,LDDrate=0,PropaguleEstablishment=1,EnvEstabProb=1,
  DetectionProb=0,ManageProb=0,MortalityProb=0,FecundityReduction=0,SpreadReduction=0,
  InfoRadius=0,InfoTransferProb=0,ExternalInfoProb=0,OngoingExternalInfo=FALSE,
  OngoingExternalInvasion=FALSE,OutputDir=tempdir(),SaveResults=FALSE,DoProgress=FALSE,
  Seed=20260924L,Biocontrol=bc,BiocontrolPointSupport=support,Pathogen=NULL)
check_pof_output(r,'Point Transition Matrix',2L,expected_no_detection=0.5)

# 8. Vertebrate Point --------------------------------------------------------
e <- load_engine('INApestPointTransitionMatrix_Biocontrol_v0.1.R',point=TRUE)
support <- e$INApestBiocontrolPointSupport(xmin=-1,xmax=3,ymin=-1,ymax=1,nrow=1L,ncol=1L)
bc <- make_point_bc(e,target_stage=1L)
point_interaction <- function(points,timestep,...) { points$x <- points$x + 0.1; points }
point_vertebrate <- list(Interaction=list(Update=point_interaction))
r <- e$INApestVertebratePoint(ModelName='POF_BC_VERTEBRATE_POINT',Nperm=2L,Ntimesteps=2L,
  Nstages=2L,Weights=c(1,1),Transition=diag(2),InitialPoints=initial_points_tm,
  SDDkernel=zero_kernel,LDDkernel=NULL,LDDrate=0,PropaguleEstablishment=1,EnvEstabProb=1,
  DetectionProb=0,ManageProb=0,MortalityProb=0,FecundityReduction=0,SpreadReduction=0,
  InfoRadius=0,InfoTransferProb=0,ExternalInfoProb=0,OngoingExternalInfo=FALSE,
  OngoingExternalInvasion=FALSE,Vertebrate=point_vertebrate,OutputDir=tempdir(),
  SaveResults=FALSE,DoProgress=FALSE,Seed=20260924L,Biocontrol=bc,
  BiocontrolPointSupport=support)
check_pof_output(r,'Vertebrate Point',2L,expected_no_detection=0.5)
check(all(abs(r$FinalPoints$x - 0.7) < 1e-12), 'Vertebrate Point specialist interaction executed with active Q')

cat('\nArchitectures executed: 8\n')
cat('OVERALL: PASS\n')
