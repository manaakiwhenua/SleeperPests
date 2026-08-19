# Reproduce the scientific FecundityReduction smoke tests and demonstration
# tables using serial current-source functions only. This path is independent
# of the number of CPU cores on the host and avoids parallel RNG differences.

bundle_root <- normalizePath(Sys.getenv("INAPEST_BUNDLE_ROOT", unset = "."), mustWork = TRUE)
out_dir <- Sys.getenv("INAPEST_TEST_OUT", unset = file.path(bundle_root, "results", "native_serial_reproduction"))
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
model_out <- tempdir(); if (!grepl("/$", model_out)) model_out <- paste0(model_out, "/")

run_silent <- function(fun, args, seed = 20260818L) {
  set.seed(seed)
  invisible(capture.output(do.call(fun, args)))
}
final_mean_3d <- function(x, timestep) {
  vals <- vapply(seq_len(dim(x)[3]), function(i) sum(x[, timestep, i]), numeric(1))
  c(mean = mean(vals), sd = sd(vals), q05 = unname(quantile(vals, .05)), q95 = unname(quantile(vals, .95)))
}
final_mean_4d <- function(x, timestep) {
  vals <- vapply(seq_len(dim(x)[4]), function(i) sum(x[, , timestep, i]), numeric(1))
  c(mean = mean(vals), sd = sd(vals), q05 = unname(quantile(vals, .05)), q95 = unname(quantile(vals, .95)))
}
final_mean_tm <- function(x, timestep) {
  vals <- vapply(seq_len(dim(x)[4]), function(i) sum(x[, , timestep, i]), numeric(1))
  c(mean = mean(vals), sd = sd(vals), q05 = unname(quantile(vals, .05)), q95 = unname(quantile(vals, .95)))
}

source_dir <- file.path(bundle_root, "source", "latest")
source(file.path(source_dir, "INApestMeta.r")); meta_serial <- INApestMeta
source(file.path(source_dir, "INApestMetaMultipleLandUse.r")); mlu_serial <- INApestMetaMultipleLandUse
source(file.path(source_dir, "INApestMetaTransitionMatrix.r")); tm_serial <- INApestMetaTransitionMatrix

# ---------------------------------------------------------------------
# 2. Parameterisation smoke tests
# ---------------------------------------------------------------------
smoke <- data.frame(model=character(),parameterisation=character(),final_mean=numeric(),final_sd=numeric(),status=character())

# Ordinary Meta: scalar, node vector, node x time matrix
n <- 4L; nt <- 3L; np <- 5L; S <- diag(n)
base_meta <- list(ModelName='',Nperm=np,Ntimesteps=nt,DetectionProb=0,DetectionSD=0,ManageProb=1,ManageSD=0,MortalityProb=.2,MortalitySD=0,
  SpreadReduction=0,SpreadReductionSD=0,InitialPopulation=c(20,16,12,8),InitBioP=NA,InvasionRisk=rep(0,n),InitialInfo=rep(1,n),InitInfoP=NA,
  ExternalInfoProb=rep(0,n),InfoRetentionProb=1,InfoPersistenceSteps=NA,EnvEstabProb=1,Survival=1,K=rep(500,n),PropaguleProduction=.6,
  PropaguleEstablishment=.0015,IncursionStartPop=NA,SDDprob=S,SEAM=0,LDDprob=NA,LDDrate=0,OngoingExternalInvasion=FALSE,
  OngoingExternalInfo=FALSE,OutputDir=model_out,DoPlots=FALSE)
meta_shapes <- list(
  scalar=.4,
  node_vector=c(0,.2,.4,.6),
  node_by_timestep=matrix(c(0,.2,.4,.6, .2,.4,.6,.8, .4,.6,.8,1),nrow=n,ncol=nt)
)
for(i in seq_along(meta_shapes)) {
  nm <- names(meta_shapes)[i]; a <- base_meta; a$ModelName <- paste0('shape_meta_',i,'_'); a$FecundityReduction <- meta_shapes[[i]]
  run_silent(meta_serial,a,1000L+i); p <- readRDS(paste0(model_out,a$ModelName,'PopulationLargeOut.rds')); z <- final_mean_3d(p,nt)
  smoke <- rbind(smoke,data.frame(model='INApestMeta',parameterisation=nm,final_mean=z['mean'],final_sd=z['sd'],status='PASS'))
}

# Multiple-land-use: scalar, land-use vector, node x land-use, node x land-use x time
n <- 4L; lu <- 3L; nt <- 3L; np <- 5L; S <- diag(n); Kmat <- matrix(300,n,lu); init <- matrix(c(15,10,5),n,lu,byrow=TRUE)
base_mlu <- list(ModelName='',Nperm=np,Ntimesteps=nt,Nlanduses=lu,DetectionProb=rep(0,lu),DetectionSD=rep(0,lu),ManageProb=rep(1,lu),ManageSD=rep(0,lu),
  MortalityProb=rep(.2,lu),MortalitySD=rep(0,lu),SpreadReduction=rep(0,lu),SpreadReductionSD=rep(0,lu),InitialPopulation=init,InitBioP=NA,
  InvasionRisk=rep(0,n),InitialInfo=rep(1,n),InitInfoP=0,ExternalInfoProb=rep(0,n),InfoRetentionProb=1,InfoPersistenceSteps=NA,
  EnvEstabProb=1,Survival=1,K=Kmat,PropaguleProduction=.6,PropaguleEstablishment=.0015,IncursionStartPop=NA,SDDprob=S,SEAM=0,LDDprob=NA,
  LDDrate=0,OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=model_out,DoPlots=FALSE)
mlu_matrix <- outer(seq(0,.6,length.out=n),c(0,.2,.4),'+'); mlu_matrix[mlu_matrix>1] <- 1
mlu_array <- array(0,dim=c(n,lu,nt)); for(t in 1:nt) mlu_array[,,t] <- pmin(1,mlu_matrix + (t-1)*.1)
mlu_shapes <- list(scalar=.4,landuse_vector=c(0,.4,.8),node_by_landuse=mlu_matrix,node_by_landuse_by_timestep=mlu_array)
for(i in seq_along(mlu_shapes)) {
  nm <- names(mlu_shapes)[i]; a <- base_mlu; a$ModelName <- paste0('shape_mlu_',i,'_'); a$FecundityReduction <- mlu_shapes[[i]]
  run_silent(mlu_serial,a,2000L+i); p <- readRDS(paste0(model_out,a$ModelName,'PopulationLargeOut.rds')); z <- final_mean_4d(p,nt)
  smoke <- rbind(smoke,data.frame(model='INApestMetaMultipleLandUse',parameterisation=nm,final_mean=z['mean'],final_sd=z['sd'],status='PASS'))
}

# Transition matrix: scalar, node vector, stage vector, node x time, node x stage x time
n <- 4L; nt <- 3L; np <- 5L; S <- diag(n); A <- matrix(c(.4,1.5,3,.4,.5,0,0,.3,.8),3,3,byrow=TRUE); init <- cbind(rep(0,n),rep(12,n),rep(4,n))
base_tm <- list(ModelName='',Nperm=np,Ntimesteps=nt,Nstages=3L,Weights=c(1,1,1),Transition=A,DetectionProb=c(0,0,0),DetectionSD=0,
  ManageProb=1,ManageSD=0,MortalityProb=c(.2,.2,.2),MortalitySD=0,SpreadReduction=0,SpreadReductionSD=0,InitialPopulation=init,InitBioP=NA,
  InvasionRisk=rep(0,n),InitialInfo=rep(1,n),InitInfoP=NA,ExternalInfoProb=rep(0,n),InfoRetentionProb=1,InfoPersistenceSteps=NA,
  EnvEstabProb=1,K=rep(1000,n),SeedbankK=rep(1000,n),PropaguleEstablishment=.01,IncursionStartPop=NA,SDDprob=S,SEAM=0,LDDprob=NA,
  LDDrate=0,DispersalDensityFactor=0,BlockedTransitionMortality=0,OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=model_out,DoPlots=FALSE)
tm_node_time <- matrix(c(0,.2,.4,.6, .2,.4,.6,.8, .4,.6,.8,1),nrow=n,ncol=nt)
tm_arr <- array(.2,dim=c(n,3,nt)); tm_arr[,3,] <- .7; tm_arr[,,3] <- pmin(1,tm_arr[,,3]+.2)
tm_shapes <- list(scalar=.4,node_vector=c(0,.2,.4,.6),stage_vector=c(0,.2,.8),node_by_timestep=tm_node_time,node_by_stage_by_timestep=tm_arr)
for(i in seq_along(tm_shapes)) {
  nm <- names(tm_shapes)[i]; a <- base_tm; a$ModelName <- paste0('shape_tm_',i,'_'); a$FecundityReduction <- tm_shapes[[i]]
  run_silent(tm_serial,a,3000L+i); p <- readRDS(paste0(model_out,a$ModelName,'PopulationStageLargeOut.rds')); z <- final_mean_tm(p,nt)
  smoke <- rbind(smoke,data.frame(model='INApestMetaTransitionMatrix',parameterisation=nm,final_mean=z['mean'],final_sd=z['sd'],status='PASS'))
}

# Stage 1 is intentionally unused for fecundity: changing only its reduction must not alter output.
a <- base_tm; a$ModelName <- 'stage1_a_'; a$FecundityReduction <- c(0,.3,.8); run_silent(tm_serial,a,4444L)
p1 <- readRDS(paste0(model_out,'stage1_a_PopulationStageLargeOut.rds'))
a$ModelName <- 'stage1_b_'; a$FecundityReduction <- c(.95,.3,.8); run_silent(tm_serial,a,4444L)
p2 <- readRDS(paste0(model_out,'stage1_b_PopulationStageLargeOut.rds'))
smoke <- rbind(smoke,data.frame(model='INApestMetaTransitionMatrix',parameterisation='stage1_value_unused_check',final_mean=max(abs(p1-p2)),final_sd=0,status=if(max(abs(p1-p2))==0)'PASS' else 'FAIL'))

# Ambiguous vector check when nodes == stages.
amb_args <- base_tm; amb_args$Nperm <- 1L; amb_args$Ntimesteps <- 1L; amb_args$Nstages <- 3L; amb_args$SDDprob <- diag(3); amb_args$InitialPopulation <- matrix(c(0,12,4),3,3,byrow=TRUE)
amb_args$InvasionRisk <- rep(0,3); amb_args$InitialInfo <- rep(1,3); amb_args$ExternalInfoProb <- rep(0,3); amb_args$K <- rep(1000,3); amb_args$SeedbankK <- rep(1000,3); amb_args$FecundityReduction <- c(.1,.2,.3); amb_args$ModelName <- 'amb_'
amb_msg <- tryCatch({run_silent(tm_serial,amb_args,5555L); ''},error=function(e) conditionMessage(e))
smoke <- rbind(smoke,data.frame(model='INApestMetaTransitionMatrix',parameterisation='ambiguous_node_stage_vector_error',final_mean=NA,final_sd=NA,status=if(grepl('ambiguous',amb_msg,ignore.case=TRUE))'PASS' else paste('FAIL',amb_msg)))
write.csv(smoke,file.path(out_dir,'parameterisation_smoke_results.csv'),row.names=FALSE)

# ---------------------------------------------------------------------
# 3. Trial simulations: scalar fecundity reduction x mortality
# ---------------------------------------------------------------------
n <- 4L; nt <- 6L; np <- 60L; S <- diag(n)
base_effect <- list(ModelName='',Nperm=np,Ntimesteps=nt,DetectionProb=0,DetectionSD=0,ManageProb=1,ManageSD=0,
  MortalitySD=0,SpreadReduction=0,SpreadReductionSD=0,InitialPopulation=rep(25,n),InitBioP=NA,InvasionRisk=rep(0,n),InitialInfo=rep(1,n),
  InitInfoP=NA,ExternalInfoProb=rep(0,n),InfoRetentionProb=1,InfoPersistenceSteps=NA,EnvEstabProb=1,Survival=1,K=rep(500,n),
  PropaguleProduction=.6,PropaguleEstablishment=.0015,IncursionStartPop=NA,SDDprob=S,SEAM=0,LDDprob=NA,LDDrate=0,
  OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=model_out,DoPlots=FALSE)
effect_meta <- data.frame()
k <- 0L
for(mort in c(0,.3,.6)) for(fr in c(0,.25,.5,.75,1)) {
  k <- k+1L; a <- base_effect; a$MortalityProb <- mort; a$FecundityReduction <- fr; a$ModelName <- paste0('effect_meta_',k,'_')
  run_silent(meta_serial,a,6000L+k); p <- readRDS(paste0(model_out,a$ModelName,'PopulationLargeOut.rds')); z <- final_mean_3d(p,nt)
  effect_meta <- rbind(effect_meta,data.frame(mortality=mort,fecundity_reduction=fr,final_population_mean=z['mean'],final_population_sd=z['sd'],q05=z['q05'],q95=z['q95']))
}
effect_meta$ratio_to_no_fecundity_reduction <- ave(effect_meta$final_population_mean,effect_meta$mortality,FUN=function(x)x/x[1])
write.csv(effect_meta,file.path(out_dir,'mortality_by_fecundity_results.csv'),row.names=FALSE)

# ---------------------------------------------------------------------
# 4. Trial simulations: land-use-specific reduction x mortality
# ---------------------------------------------------------------------
n <- 4L; lu <- 3L; nt <- 6L; np <- 60L; S <- diag(n)
Klu <- rbind(c(500,0,0),c(0,500,0),c(0,0,500),c(200,200,200))
initlu <- rbind(c(25,0,0),c(0,25,0),c(0,0,25),c(10,10,10))
base_lu_effect <- list(ModelName='',Nperm=np,Ntimesteps=nt,Nlanduses=lu,DetectionProb=rep(0,lu),DetectionSD=rep(0,lu),ManageProb=rep(1,lu),ManageSD=rep(0,lu),
  MortalitySD=rep(0,lu),SpreadReduction=rep(0,lu),SpreadReductionSD=rep(0,lu),InitialPopulation=initlu,InitBioP=NA,InvasionRisk=rep(0,n),InitialInfo=rep(1,n),
  InitInfoP=0,ExternalInfoProb=rep(0,n),InfoRetentionProb=1,InfoPersistenceSteps=NA,EnvEstabProb=1,Survival=1,K=Klu,PropaguleProduction=.6,
  PropaguleEstablishment=.0015,IncursionStartPop=NA,SDDprob=S,SEAM=0,LDDprob=NA,LDDrate=0,OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,
  OutputDir=model_out,DoPlots=FALSE)
patterns <- list(none=c(0,0,0),uniform_half=c(.5,.5,.5),landuse_gradient=c(0,.5,.9),uniform_strong=c(.9,.9,.9))
effect_lu <- data.frame(); k <- 0L
for(mort in c(0,.3,.6)) for(pat in names(patterns)) {
  k <- k+1L; a <- base_lu_effect; a$MortalityProb <- rep(mort,lu); a$FecundityReduction <- patterns[[pat]]; a$ModelName <- paste0('effect_lu_',k,'_')
  run_silent(mlu_serial,a,7000L+k); p <- readRDS(paste0(model_out,a$ModelName,'PopulationLargeOut.rds'))
  finals <- sapply(seq_len(np),function(i) rowSums(p[,,nt,i]))
  effect_lu <- rbind(effect_lu,data.frame(mortality=mort,pattern=pat,FR_LU1=patterns[[pat]][1],FR_LU2=patterns[[pat]][2],FR_LU3=patterns[[pat]][3],
    final_total_mean=mean(colSums(finals)),node1_LU1_mean=mean(finals[1,]),node2_LU2_mean=mean(finals[2,]),node3_LU3_mean=mean(finals[3,]),node4_mixed_mean=mean(finals[4,])))
}
write.csv(effect_lu,file.path(out_dir,'landuse_mortality_fecundity_results.csv'),row.names=FALSE)


# Compare the three scientific evidence tables against the archived reference.
compare_files <- c("parameterisation_smoke_results.csv", "mortality_by_fecundity_results.csv", "landuse_mortality_fecundity_results.csv")
cmp <- do.call(rbind, lapply(compare_files, function(f) {
  a <- read.csv(file.path(bundle_root, "results", "archived", f), check.names=FALSE)
  b <- read.csv(file.path(out_dir, f), check.names=FALSE)
  data.frame(file=f, all_equal=isTRUE(all.equal(a,b,check.attributes=FALSE)), stringsAsFactors=FALSE)
}))
write.csv(cmp, file.path(out_dir, "archived_scientific_table_comparison.csv"), row.names=FALSE)
print(cmp)
cat("Serial scientific reproduction completed.\n")
