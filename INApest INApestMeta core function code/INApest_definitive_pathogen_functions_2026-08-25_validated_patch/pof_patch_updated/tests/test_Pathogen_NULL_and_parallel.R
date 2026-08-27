stop_equal <- function(x,y,msg) if(!isTRUE(all.equal(x,y,check.attributes=TRUE))) stop(msg)

run_meta <- function(fun, model, outdir, pathogen_missing=TRUE, nperm=2) {
  args <- list(ModelName=model,Nperm=nperm,Ntimesteps=3,
    DetectionProb=.05,DetectionSD=0,ManageProb=.2,ManageSD=0,MortalityProb=.1,MortalitySD=0,
    FecundityReduction=0,SpreadReduction=.1,SpreadReductionSD=0,
    InitialPopulation=c(10,0,0),InitBioP=NA,InvasionRisk=c(0,0,0),InitialInfo=c(0,0,0),InitInfoP=NA,
    ExternalInfoProb=c(0,0,0),InfoRetentionProb=1,InfoPersistenceSteps=NA,
    EnvEstabProb=1,Survival=.95,K=c(100,100,100),PropaguleProduction=.5,PropaguleEstablishment=.1,
    IncursionStartPop=1,SDDprob=diag(3),SEAM=0,LDDprob=diag(3),LDDrate=0,
    OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=paste0(outdir,'/'),DoPlots=FALSE)
  if(!pathogen_missing) args$Pathogen <- NULL
  set.seed(1234); do.call(fun,args)
  list(pop=readRDS(file.path(outdir,paste0(model,'PopulationLargeOut.rds'))),
       inv=readRDS(file.path(outdir,paste0(model,'InvasionLargeOut.rds'))),
       info=readRDS(file.path(outdir,paste0(model,'InfoLargeOut.rds'))),
       det=readRDS(file.path(outdir,paste0(model,'DetectedLargeOut.rds'))))
}

source('baseline/INApestMeta.r'); old_meta <- INApestMeta
source('INApestMetaPathogen.r'); new_meta <- INApestMeta
od <- tempfile('oldmeta'); nd <- tempfile('newmeta'); dir.create(od); dir.create(nd)
a <- run_meta(old_meta,'m_',od,TRUE); b <- run_meta(new_meta,'m_',nd,FALSE)
stop_equal(a,b,'Pathogen=NULL must preserve serial Meta outputs under fixed seed')
cat('PASS: Meta Pathogen=NULL backward compatibility\n')

# Deterministic serial/parallel comparison exercises the same Pathogen object.
serial_dir <- tempfile('serialp'); parallel_dir <- tempfile('parallelp'); dir.create(serial_dir); dir.create(parallel_dir)
path <- INApestPathogen('SIR',Beta=0,RecoveryProb=0,InitialInfected=c(2,0,0))
args_path <- list(ModelName='p_',Nperm=2,Ntimesteps=2,Pathogen=path,
    DetectionProb=0,DetectionSD=0,ManageProb=0,ManageSD=0,MortalityProb=0,MortalitySD=0,
    FecundityReduction=0,SpreadReduction=0,SpreadReductionSD=0,
    InitialPopulation=c(10,0,0),InitBioP=NA,InvasionRisk=c(0,0,0),InitialInfo=c(0,0,0),InitInfoP=NA,
    ExternalInfoProb=c(0,0,0),InfoRetentionProb=1,InfoPersistenceSteps=NA,
    EnvEstabProb=0,Survival=1,K=c(100,100,100),PropaguleProduction=0,PropaguleEstablishment=0,
    IncursionStartPop=1,SDDprob=diag(3),SEAM=0,LDDprob=diag(3),LDDrate=0,
    OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,DoPlots=FALSE)
args_s <- args_path; args_s$OutputDir <- paste0(serial_dir,'/')
do.call(new_meta,args_s)
source('INApestMetaParallelPathogen.r'); par_meta <- INApestMetaParallel
args_p <- args_path; args_p$OutputDir <- paste0(parallel_dir,'/')
do.call(par_meta,args_p)
for(suffix in c('PopulationLargeOut.rds','PathogenStateLargeOut.rds','InvasionLargeOut.rds'))
  stop_equal(readRDS(file.path(serial_dir,paste0('p_',suffix))),readRDS(file.path(parallel_dir,paste0('p_',suffix))),paste('Meta serial/parallel',suffix))
cat('PASS: Meta serial/parallel deterministic pathogen parity\n')
