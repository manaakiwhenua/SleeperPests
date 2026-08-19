# Native-R validation for the standardised PSOCK/parLapply implementation.
# Run from the directory containing the updated INApest source files.
# This complements validation_fecundity_reduction.R, which was executed in WebR.

cores <- parallel::detectCores()
if (is.na(cores) || cores < 3L) {
  stop("This validation needs at least 3 detected cores so the INApest parallel functions choose >1 PSOCK worker.")
}

out_dir <- file.path(tempdir(), "inapest_native_psock")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- paste0(out_dir, .Platform$file.sep)

run_silent <- function(fun, args, seed) {
  set.seed(seed)
  invisible(capture.output(do.call(fun, args)))
}

summarise_final <- function(x, time_dim, perm_dim, timestep) {
  perm_n <- dim(x)[perm_dim]
  vals <- vapply(seq_len(perm_n), function(i) {
    idx <- rep(list(TRUE), length(dim(x)))
    idx[[time_dim]] <- timestep
    idx[[perm_dim]] <- i
    sum(do.call(`[`, c(list(x), idx, list(drop = TRUE))))
  }, numeric(1))
  c(mean = mean(vals), se = stats::sd(vals) / sqrt(length(vals)))
}

compare_mc <- function(serial_summary, parallel_summary) {
  denom <- sqrt(serial_summary["se"]^2 + parallel_summary["se"]^2)
  z <- if (denom == 0) ifelse(serial_summary["mean"] == parallel_summary["mean"], 0, Inf) else
    abs(serial_summary["mean"] - parallel_summary["mean"]) / denom
  c(serial_mean = serial_summary["mean"], parallel_mean = parallel_summary["mean"], z = z, pass = z <= 4)
}

source("INApest.R"); ina_serial <- INApest
source("INApestParallel.R"); ina_parallel <- INApestParallel
source("INApestMetaMultipleLandUse.r"); mlu_serial <- INApestMetaMultipleLandUse
source("INApestMetaParallelMultipleLandUse.r"); mlu_parallel <- INApestMetaParallelMultipleLandUse
source("INApestMetaTransitionMatrix.r"); tm_serial <- INApestMetaTransitionMatrix
source("INApestMetaTransitionMatrixParallel.r"); tm_parallel <- INApestMetaTransitionMatrixParallel

results <- data.frame(model=character(), serial_mean=numeric(), parallel_mean=numeric(), z=numeric(), pass=logical())

# Binary INApest -------------------------------------------------------
n <- 5L; nt <- 5L; np <- 120L
S <- matrix(c(.7,.2,0,0,0, .1,.6,.2,0,0, 0,.1,.6,.2,0, 0,0,.1,.6,.2, 0,0,0,.1,.7), n, n, byrow=TRUE)
args <- list(ModelName='native_ina_s_',Nperm=np,Ntimesteps=nt,DetectionProb=0,DetectionSD=0,ManageProb=1,ManageSD=0,
  EradicationProb=.25,EradicationSD=0,SpreadReduction=.2,SpreadReductionSD=0,InitialInvasion=c(1,0,0,0,0),InitBioP=NA,
  InvasionRisk=rep(0,n),InitialInfo=rep(1,n),InitInfoP=NA,ExternalInfoProb=rep(0,n),InfoRetentionProb=1,InfoPersistenceSteps=NA,
  EnvEstabProb=1,Survival=1,SDDprob=S,SEAM=0,LDDprob=0,OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,
  OutputDir=out_dir,DoPlots=FALSE)
run_silent(ina_serial,args,101L); xs <- readRDS(paste0(out_dir,'native_ina_s_InvasionLargeOut.rds'))
args$ModelName <- 'native_ina_p_'; run_silent(ina_parallel,args,102L); xp <- readRDS(paste0(out_dir,'native_ina_p_InvasionLargeOut.rds'))
z <- compare_mc(summarise_final(xs,2,3,nt),summarise_final(xp,2,3,nt))
results <- rbind(results,data.frame(model='INApest',serial_mean=z[1],parallel_mean=z[2],z=z[3],pass=as.logical(z[4])))

# Multiple-land-use Meta ---------------------------------------------
n <- 4L; lu <- 3L; nt <- 5L; np <- 120L; S <- diag(n)
Kmat <- matrix(c(300,400,500),nrow=n,ncol=lu,byrow=TRUE)
init <- matrix(c(20,12,8,16,10,6,12,8,4,8,5,3),nrow=n,byrow=TRUE)
args <- list(ModelName='native_mlu_s_',Nperm=np,Ntimesteps=nt,Nlanduses=lu,DetectionProb=rep(0,lu),DetectionSD=rep(0,lu),
  ManageProb=rep(1,lu),ManageSD=rep(0,lu),MortalityProb=rep(.2,lu),MortalitySD=rep(0,lu),FecundityReduction=c(.1,.5,.9),
  SpreadReduction=rep(0,lu),SpreadReductionSD=rep(0,lu),InitialPopulation=init,InitBioP=NA,InvasionRisk=rep(0,n),InitialInfo=rep(1,n),
  InitInfoP=0,ExternalInfoProb=rep(0,n),InfoRetentionProb=1,InfoPersistenceSteps=NA,EnvEstabProb=1,Survival=1,K=Kmat,
  PropaguleProduction=2,PropaguleEstablishment=.02,IncursionStartPop=NA,SDDprob=S,SEAM=0,LDDprob=NA,LDDrate=0,
  OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=out_dir,DoPlots=FALSE)
run_silent(mlu_serial,args,201L); xs <- readRDS(paste0(out_dir,'native_mlu_s_PopulationLargeOut.rds'))
args$ModelName <- 'native_mlu_p_'; run_silent(mlu_parallel,args,202L); xp <- readRDS(paste0(out_dir,'native_mlu_p_PopulationLargeOut.rds'))
z <- compare_mc(summarise_final(xs,3,4,nt),summarise_final(xp,3,4,nt))
results <- rbind(results,data.frame(model='INApestMetaMultipleLandUse',serial_mean=z[1],parallel_mean=z[2],z=z[3],pass=as.logical(z[4])))

# Transition-matrix Meta ---------------------------------------------
A <- matrix(c(.2,3,4,.5,.6,0,0,.2,.8),3,3,byrow=TRUE)
init_tm <- cbind(rep(0,n),c(12,10,8,6),c(4,4,3,2))
args <- list(ModelName='native_tm_s_',Nperm=np,Ntimesteps=nt,Nstages=3L,Weights=c(1,1,1),Transition=A,
  DetectionProb=c(0,0,0),DetectionSD=0,ManageProb=1,ManageSD=0,MortalityProb=c(.2,.2,.2),MortalitySD=0,FecundityReduction=.4,
  SpreadReduction=0,SpreadReductionSD=0,InitialPopulation=init_tm,InitBioP=NA,InvasionRisk=rep(0,n),InitialInfo=rep(1,n),InitInfoP=NA,
  ExternalInfoProb=rep(0,n),InfoRetentionProb=1,InfoPersistenceSteps=NA,EnvEstabProb=1,K=rep(1000,n),SeedbankK=rep(1000,n),
  PropaguleEstablishment=.02,IncursionStartPop=NA,SDDprob=S,SEAM=0,LDDprob=NA,LDDrate=0,DispersalDensityFactor=0,
  BlockedTransitionMortality=0,OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=out_dir,DoPlots=FALSE)
run_silent(tm_serial,args,301L); xs <- readRDS(paste0(out_dir,'native_tm_s_PopulationStageLargeOut.rds'))
args$ModelName <- 'native_tm_p_'; run_silent(tm_parallel,args,302L); xp <- readRDS(paste0(out_dir,'native_tm_p_PopulationStageLargeOut.rds'))
z <- compare_mc(summarise_final(xs,3,4,nt),summarise_final(xp,3,4,nt))
results <- rbind(results,data.frame(model='INApestMetaTransitionMatrix',serial_mean=z[1],parallel_mean=z[2],z=z[3],pass=as.logical(z[4])))

print(results)
if (!all(results$pass)) stop("At least one serial/parallel Monte Carlo comparison exceeded 4 combined standard errors.")
cat("PASS: true PSOCK execution completed and all tested parallel distributions were close to serial references.\n")
