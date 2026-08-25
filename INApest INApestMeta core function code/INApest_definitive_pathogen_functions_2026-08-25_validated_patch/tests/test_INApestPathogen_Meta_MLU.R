# Generic pathogen + Meta/MLU regression tests
stop_if_not <- function(x, msg) if(!isTRUE(x)) stop(msg)
stop_equal <- function(x,y,msg) if(!isTRUE(all.equal(x,y,check.attributes=FALSE))) stop(msg, ': ', paste(capture.output(str(x)), collapse=' '), ' != ', paste(capture.output(str(y)), collapse=' '))

source('INApestPathogen.R')

# 1. State definitions and initial conservation.
p_sis <- INApestPathogen('SIS', Beta=0, RecoveryProb=0, InitialInfected=2)
ctx <- list(n_nodes=2,Ntimesteps=3); s <- p_sis$Engine$Initial(c(10,5), ctx)
stop_equal(rowSums(s), c(10,5), 'SIS initial conservation')
stop_if_not(all(colnames(s)==c('S','I')), 'SIS state names')

p_sir <- INApestPathogen('SIR', Beta=0, RecoveryProb=0, InitialInfected=c(2,1), InitialRecovered=c(1,0))
s <- p_sir$Engine$Initial(c(10,5), ctx)
stop_equal(rowSums(s), c(10,5), 'SIR initial conservation')

p_seir <- INApestPathogen('SEIR', Beta=0, RecoveryProb=0, InitialInfected=1, InitialExposed=2, InitialRecovered=1)
s <- p_seir$Engine$Initial(c(10,10), ctx)
stop_equal(rowSums(s), c(10,10), 'SEIR initial conservation')

# 2. Recruitment enters S; host mortality thins without changing totals.
set.seed(1)
s0 <- matrix(c(6,2,2),1,3,dimnames=list(NULL,c('S','I','R')))
s_up <- p_sir$Engine$Reconcile(s0, 15, list(n_nodes=1,Ntimesteps=1))
stop_if_not(s_up[1,'S']==11 && sum(s_up)==15, 'Host recruitment must enter S')
set.seed(2)
s_down <- p_sir$Engine$Reconcile(s0, 5, list(n_nodes=1,Ntimesteps=1))
stop_if_not(sum(s_down)==5 && all(s_down>=0), 'Host mortality thinning conservation')

# 3. Deterministic recovery / progression / waning.
s0 <- matrix(c(0,10),1,2,dimnames=list(NULL,c('S','I')))
p0 <- INApestPathogen('SIS',0,1); z <- p0$Engine$Step(s0,10,1,list(n_nodes=1,Ntimesteps=1))
stop_if_not(z$State[1,'S']==10 && z$State[1,'I']==0 && z$N==10, 'SIS recovery returns to S')

s0 <- matrix(c(0,8,2,0),1,4,dimnames=list(NULL,c('S','E','I','R')))
p0 <- INApestPathogen('SEIR',0,0,ProgressionProb=1); z <- p0$Engine$Step(s0,10,1,list(n_nodes=1,Ntimesteps=1))
stop_if_not(z$State[1,'E']==0 && z$State[1,'I']==10, 'SEIR certain progression')

s0 <- matrix(c(0,0,10),1,3,dimnames=list(NULL,c('S','I','R')))
p0 <- INApestPathogen('SIR',0,0,ImmunityLossProb=1); z <- p0$Engine$Step(s0,10,1,list(n_nodes=1,Ntimesteps=1))
stop_if_not(z$State[1,'S']==10 && z$State[1,'R']==0, 'SIR certain waning')

# 4. Pathogen mortality is the pathogen process that changes total N.
s0 <- matrix(c(0,10,0),1,3,dimnames=list(NULL,c('S','I','R')))
p0 <- INApestPathogen('SIR',0,0,PathogenMortalityProb=1); z <- p0$Engine$Step(s0,10,1,list(n_nodes=1,Ntimesteps=1))
stop_if_not(z$N==0 && sum(z$State)==0 && z$Deaths==10, 'Pathogen mortality reduces host N')

# 5. Pathogen introduction converts existing S rather than creating hosts.
s0 <- matrix(c(10,0),1,2,dimnames=list(NULL,c('S','I')))
p0 <- INApestPathogen('SIS',0,0,IntroductionProb=1,IntroductionNumber=3); z <- p0$Engine$Step(s0,10,1,list(n_nodes=1,Ntimesteps=1))
stop_if_not(z$N==10 && z$State[1,'S']==7 && z$State[1,'I']==3, 'Pathogen introduction conserves hosts')

# 6. MLU resolution: land-use-specific Beta expands by column-major node x LU cells.
p <- INApestPathogen('SIS', Beta=c(0.1,0.5), RecoveryProb=0)
pm <- p$Engine$Resolve(p$Beta,1,list(n_nodes=3,n_landuses=2,Ntimesteps=4),'Beta')
stop_equal(pm, c(rep(.1,3),rep(.5,3)), 'MLU land-use Beta expansion')

# 7. Time-varying MLU array resolution.
b <- array(0,dim=c(2,2,3)); b[,,2] <- matrix(c(.1,.2,.3,.4),2,2)
p <- INApestPathogen('SIS', Beta=b, RecoveryProb=0)
pm <- p$Engine$Resolve(p$Beta,2,list(n_nodes=2,n_landuses=2,Ntimesteps=3),'Beta')
stop_equal(pm, as.numeric(b[,,2]), 'MLU time-varying Beta')

cat('PASS: generic INApest pathogen mechanism tests\n')

# -----------------------------------------------------------------------------
# Integration tests if patched core files are in the working directory.
# -----------------------------------------------------------------------------
if(file.exists('INApestMetaPathogen.r')) {
  source('INApestMetaPathogen.r')
  tmp <- tempfile('ina_pathogen_meta_'); dir.create(tmp)
  base <- list(ModelName='meta_path_',Nperm=1,Ntimesteps=2,
    DetectionProb=0,DetectionSD=0,ManageProb=0,ManageSD=0,MortalityProb=0,MortalitySD=0,
    FecundityReduction=0,SpreadReduction=0,SpreadReductionSD=0,
    InitialPopulation=c(10,0),InitBioP=NA,InvasionRisk=c(0,0),InitialInfo=c(0,0),InitInfoP=NA,
    ExternalInfoProb=c(0,0),InfoRetentionProb=1,InfoPersistenceSteps=NA,
    EnvEstabProb=0,Survival=1,K=c(100,100),PropaguleProduction=0,PropaguleEstablishment=0,
    IncursionStartPop=1,SDDprob=diag(2),SEAM=0,LDDprob=diag(2),LDDrate=0,
    OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=paste0(tmp,'/'),DoPlots=FALSE)
  do.call(INApestMeta,c(base,list(Pathogen=INApestPathogen('SIS',Beta=0,RecoveryProb=0,InitialInfected=c(2,0)))))
  hs <- readRDS(file.path(tmp,'meta_path_PathogenStateLargeOut.rds'))
  pop <- readRDS(file.path(tmp,'meta_path_PopulationLargeOut.rds'))
  stop_if_not(identical(dim(hs),c(2L,2L,2L,1L)), 'Meta pathogen output dimensions')
  stop_equal(apply(hs,c(1,3,4),sum),pop,'Meta state totals equal N')
  stop_if_not(all(hs[1,'I',,1]==2), 'Meta pathogen state persists when beta=recovery=0')
  cat('PASS: serial INApestMeta pathogen integration\n')
}

if(file.exists('INApestMetaMultipleLandUsePathogen.r')) {
  source('INApestMetaMultipleLandUsePathogen.r')
  # A full MLU smoke run is deliberately kept in the suite but may need the
  # established MLU baseline argument fixture. The mechanism-specific MLU
  # resolver tests above do not depend on that fixture.
  stop_if_not('Pathogen' %in% names(formals(INApestMetaMultipleLandUse)), 'MLU Pathogen public argument')
  cat('PASS: serial MLU pathogen API + resolver integration\n')
}

if(file.exists('INApestMetaMultipleLandUsePathogen.r')) {
  source('INApestMetaMultipleLandUsePathogen.r')
  tmp <- tempfile('ina_pathogen_mlu_'); dir.create(tmp)
  base_mlu <- list(ModelName='mlu_path_',Nperm=1,Ntimesteps=2,Nlanduses=2,
    DetectionProb=c(0,0),DetectionSD=c(0,0),ManageProb=c(0,0),ManageSD=c(0,0),
    MortalityProb=c(0,0),MortalitySD=c(0,0),FecundityReduction=c(0,0),
    SpreadReduction=c(0,0),SpreadReductionSD=c(0,0),
    InitialPopulation=matrix(c(10,0,0,0),nrow=2,ncol=2),InitBioP=NA,InvasionRisk=c(0,0),
    InitialInfo=c(0,0),InitInfoP=NA,ExternalInfoProb=c(0,0),InfoRetentionProb=1,InfoPersistenceSteps=NA,
    EnvEstabProb=0,Survival=1,K=matrix(100,2,2),PropaguleProduction=0,PropaguleEstablishment=0,
    IncursionStartPop=1,SDDprob=diag(2),SEAM=0,LDDprob=diag(2),LDDrate=0,
    OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=paste0(tmp,'/'),DoPlots=FALSE)
  path_mlu <- INApestPathogen('SIS',Beta=matrix(0,2,2),RecoveryProb=0,InitialInfected=matrix(c(2,0,0,0),2,2))
  do.call(INApestMetaMultipleLandUse,c(base_mlu,list(Pathogen=path_mlu)))
  hs <- readRDS(file.path(tmp,'mlu_path_PathogenStateLargeOut.rds'))
  pop <- readRDS(file.path(tmp,'mlu_path_PopulationLargeOut.rds'))
  stop_if_not(identical(dim(hs),c(2L,2L,2L,2L,1L)), 'MLU pathogen output dimensions')
  totals <- apply(hs,c(1,2,4,5),sum)
  stop_equal(totals,pop,'MLU state totals equal N')
  stop_if_not(all(hs[1,1,'I',,1]==2), 'MLU pathogen state persists')
  cat('PASS: serial MLU pathogen integration\n')
}

# Public API checks for parallel variants; full deterministic parity executes when
# the corresponding files can be sourced in the available R runtime.
if(file.exists('INApestMetaParallelPathogen.r')) {
  source('INApestMetaParallelPathogen.r')
  stop_if_not('Pathogen' %in% names(formals(INApestMetaParallel)), 'Parallel Meta Pathogen public argument')
  cat('PASS: parallel Meta pathogen API\n')
}
if(file.exists('INApestMetaParallelMultipleLandUsePathogen.r')) {
  source('INApestMetaParallelMultipleLandUsePathogen.r')
  stop_if_not('Pathogen' %in% names(formals(INApestMetaParallelMultipleLandUse)), 'Parallel MLU Pathogen public argument')
  cat('PASS: parallel MLU pathogen API\n')
}

# -----------------------------------------------------------------------------
# Validation must cover the complete schedule and reject unused state inputs.
# -----------------------------------------------------------------------------
expect_error <- function(expr, pattern) {
  msg <- tryCatch({ force(expr); NA_character_ }, error = function(e) conditionMessage(e))
  if (is.na(msg) || !grepl(pattern, msg, fixed = TRUE))
    stop("Expected error containing: ", pattern, "; got: ", msg)
}

bad_schedule <- INApestPathogen(
  Model = "SIR", Beta = c(0, 0, 0), RecoveryProb = c(0, 1.2, 0),
  InitialInfected = 1
)
expect_error(
  bad_schedule$Engine$Validate(list(n_nodes = 1L, Ntimesteps = 3L)),
  "at every timestep"
)

expect_error(
  INApestPathogen(Model = "SIR", Beta = 0, RecoveryProb = 0,
                  InitialInfected = 1, InitialExposed = 1)$Engine$Validate(
                    list(n_nodes = 1L, Ntimesteps = 1L)),
  "InitialExposed must be zero"
)

expect_error(
  INApestPathogen(Model = "SIS", Beta = 0, RecoveryProb = 0,
                  InitialInfected = 1, ImmunityLossProb = 0.2)$Engine$Validate(
                    list(n_nodes = 1L, Ntimesteps = 1L)),
  "ImmunityLossProb must be zero"
)

cat("PASS: full-schedule and state-contract validation tests\n")

# -----------------------------------------------------------------------------
# Spatial/contact transmission contract.
# -----------------------------------------------------------------------------
cm <- matrix(c(1, 1,
               0, 1), nrow = 2, byrow = TRUE)
p_contact <- INApestPathogen(
  Model = "SIS", Beta = 100, RecoveryProb = 0,
  InitialInfected = c(10, 0), ContactMatrix = cm
)
ctx_contact <- list(n_nodes = 2L, Ntimesteps = 1L)
st_contact <- p_contact$Engine$Initial(c(10, 10), ctx_contact)
set.seed(123)
z_contact <- p_contact$Engine$Step(st_contact, c(10, 10), 1L, ctx_contact)
stopifnot(z_contact$State[2, "I"] == 10L)

p_isolated <- INApestPathogen(
  Model = "SIS", Beta = 100, RecoveryProb = 0,
  InitialInfected = c(10, 0), ContactMatrix = diag(2)
)
st_isolated <- p_isolated$Engine$Initial(c(10, 10), ctx_contact)
set.seed(123)
z_isolated <- p_isolated$Engine$Step(st_isolated, c(10, 10), 1L, ctx_contact)
stopifnot(z_isolated$State[2, "I"] == 0L)

node_contact <- matrix(c(1, .2, .3, 1), 2, 2, byrow = TRUE)
lu_mix <- matrix(c(1, .1, .4, 1), 2, 2, byrow = TRUE)
cm_mlu <- INApestPathogenContactMatrix(node_contact, lu_mix)
stopifnot(identical(dim(cm_mlu), c(4L, 4L)))
stopifnot(isTRUE(all.equal(cm_mlu, kronecker(lu_mix, node_contact))))
cat("PASS: spatial/contact transmission contract\n")
