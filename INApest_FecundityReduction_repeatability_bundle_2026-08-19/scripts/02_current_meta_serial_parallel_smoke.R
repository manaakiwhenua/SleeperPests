# Supplemental checks introduced after the original fecundity validation snapshot.
# These use the current source set, which now includes a serial INApestMeta and
# restored user-supplied LocalDynamics hooks.

bundle_root <- normalizePath(Sys.getenv("INAPEST_BUNDLE_ROOT", unset = "."), mustWork = TRUE)
out_dir <- Sys.getenv("INAPEST_TEST_OUT", unset = file.path(bundle_root, "results", "rerun"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
source_dir <- file.path(bundle_root, "source", "latest")
source(file.path(source_dir, "INApestMeta.r")); meta_serial <- INApestMeta; default_meta_ld <- local.dynamics
source(file.path(source_dir, "INApestMetaParallel.r")); meta_parallel <- INApestMetaParallel

model_out <- tempdir(); if (!grepl("/$", model_out)) model_out <- paste0(model_out, "/")
run_silent <- function(fun, args, seed) {
  set.seed(seed)
  invisible(capture.output(do.call(fun, args)))
}

n <- 4L; nt <- 4L; np <- 8L; S <- diag(n)
base <- list(
  ModelName='', Nperm=np, Ntimesteps=nt,
  DetectionProb=0, DetectionSD=0, ManageProb=1, ManageSD=0,
  MortalityProb=.2, MortalitySD=0,
  SpreadReduction=0, SpreadReductionSD=0,
  InitialPopulation=c(20,16,12,8), InitBioP=NA,
  InvasionRisk=rep(0,n), InitialInfo=rep(1,n), InitInfoP=NA,
  ExternalInfoProb=rep(0,n), InfoRetentionProb=1, InfoPersistenceSteps=NA,
  EnvEstabProb=1, Survival=1, K=rep(500,n), PropaguleProduction=.6,
  PropaguleEstablishment=.0015, IncursionStartPop=NA,
  SDDprob=S, SEAM=0, LDDprob=NA, LDDrate=0,
  OngoingExternalInvasion=FALSE, OngoingExternalInfo=FALSE,
  OutputDir=model_out, DoPlots=FALSE
)
shapes <- list(
  scalar=.4,
  node_vector=c(0,.2,.4,.6),
  node_by_timestep=matrix(c(0,.2,.4,.6, .2,.4,.6,.8, .4,.6,.8,1, .1,.3,.5,.7), nrow=n, ncol=nt)
)
res <- data.frame(parameterisation=character(), all_equal=logical(), max_abs_difference=numeric())
for (i in seq_along(shapes)) {
  a <- base; a$FecundityReduction <- shapes[[i]]
  a$ModelName <- paste0('cur_meta_s_', i, '_'); run_silent(meta_serial, a, 8100L+i)
  xs <- readRDS(paste0(model_out, a$ModelName, 'PopulationLargeOut.rds'))
  a$ModelName <- paste0('cur_meta_p_', i, '_'); run_silent(meta_parallel, a, 8100L+i)
  xp <- readRDS(paste0(model_out, a$ModelName, 'PopulationLargeOut.rds'))
  res <- rbind(res, data.frame(parameterisation=names(shapes)[i], all_equal=isTRUE(all.equal(xs,xp)), max_abs_difference=max(abs(xs-xp))))
}
write.csv(res, file.path(out_dir, 'current_meta_serial_parallel.csv'), row.names=FALSE)

# Custom LocalDynamics compatibility contract for fecundity management.
# Old-style custom dynamics without nodefecundityreduction remain usable when
# FecundityReduction is zero, but active reduction must not be silently ignored.
old_style_ld <- function(sddprob, nodepropaguleproduction, nodeenvestabprob, n,
                         lddprob, lddrate, k_is_0, nodeK,
                         nodepropaguleestablishment, nodespreadreduction,
                         managing, maxinteger) {
  default_meta_ld(
    sddprob=sddprob, nodepropaguleproduction=nodepropaguleproduction,
    nodeenvestabprob=nodeenvestabprob, n=n, lddprob=lddprob, lddrate=lddrate,
    k_is_0=k_is_0, nodeK=nodeK,
    nodepropaguleestablishment=nodepropaguleestablishment,
    nodespreadreduction=nodespreadreduction, nodefecundityreduction=0,
    managing=managing, maxinteger=maxinteger
  )
}
compat <- data.frame(check=character(), pass=logical(), message=character())
a <- base; a$Nperm <- 2L; a$Ntimesteps <- 2L; a$ModelName <- 'old_ld_zero_'; a$FecundityReduction <- 0; a$LocalDynamics <- old_style_ld
msg0 <- tryCatch({run_silent(meta_serial, a, 9001L); ''}, error=function(e) conditionMessage(e))
compat <- rbind(compat, data.frame(check='old_custom_function_FR_zero', pass=identical(msg0,''), message=msg0))
a$ModelName <- 'old_ld_active_'; a$FecundityReduction <- .4
msg1 <- tryCatch({run_silent(meta_serial, a, 9002L); ''}, error=function(e) conditionMessage(e))
compat <- rbind(compat, data.frame(check='old_custom_function_FR_active_errors', pass=nchar(msg1)>0 && grepl('nodefecundityreduction|FecundityReduction|LocalDynamics', msg1, ignore.case=TRUE), message=msg1))
write.csv(compat, file.path(out_dir, 'current_local_dynamics_contract.csv'), row.names=FALSE)

if (!all(res$all_equal) || !all(compat$pass)) stop('At least one supplemental current-source check failed.')
cat('PASS: current ordinary Meta serial/parallel and LocalDynamics fecundity contract checks completed.\n')
