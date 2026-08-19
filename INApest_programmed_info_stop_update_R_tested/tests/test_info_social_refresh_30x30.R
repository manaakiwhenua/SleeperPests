# Model-level tests for programmed information persistence / management stopping.
# 30 x 30 landscape (900 nodes), 4 timesteps, 2 stochastic realisations.
# Parallel functions use test-only sequential shims because WebR lacks the
# production parallel/INA package stack. Production source is not rewritten.

options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
bundle_root <- if (length(args) >= 1) normalizePath(args[1], mustWork = TRUE) else normalizePath('.', mustWork = TRUE)
baseline_dir <- file.path(bundle_root, 'baseline')
modified_dir <- file.path(bundle_root, 'modified')
out_dir <- file.path(bundle_root, 'test_output_info_persistence')
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_dir <- paste0(normalizePath(out_dir, mustWork = TRUE), .Platform$file.sep)

nr <- 30L; nc <- 30L; n_nodes <- nr * nc
n_timesteps <- 4L; n_perm <- 2L; n_landuses <- 2L
zero_route <- matrix(0, nrow = n_nodes, ncol = n_nodes)
zero_external <- rep(0, n_nodes)
all_info <- rep(1, n_nodes)
all_binary <- rep(1L, n_nodes)
all_population <- rep(20, n_nodes)
all_population_lu <- matrix(20, nrow = n_nodes, ncol = n_landuses)
all_population_stage <- matrix(0, nrow = n_nodes, ncol = 2L)
all_population_stage[, 2L] <- 20
geocoords <- cbind(x = rep(seq_len(nc), times = nr), y = rep(seq_len(nr), each = nc))
transition_matrix <- matrix(c(0,0,0,1), nrow = 2L, byrow = TRUE)

result_path <- function(model_name, suffix) file.path(out_dir, paste0(model_name, suffix))
run_and_read <- function(fun, call_args, model_name, suffixes, seed = 20260813L) {
  for (s in suffixes) unlink(result_path(model_name, s), force = TRUE)
  call_args$ModelName <- model_name; call_args$OutputDir <- out_dir
  set.seed(seed); invisible(do.call(fun, call_args))
  out <- lapply(suffixes, function(s) readRDS(result_path(model_name, s)))
  names(out) <- suffixes; out
}
run_and_read_warnings <- function(fun, call_args, model_name, suffixes, seed = 20260813L) {
  warnings <- character()
  out <- withCallingHandlers(
    run_and_read(fun, call_args, model_name, suffixes, seed),
    warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart('muffleWarning') }
  )
  list(out = out, warnings = warnings)
}
assert_identical_outputs <- function(a,b,label) {
  if (!identical(names(a), names(b))) stop(label, ': output names differ')
  for (nm in names(a)) if (!identical(a[[nm]], b[[nm]])) stop(label, ': output differed for ', nm)
  TRUE
}
management_profile <- function(x, mlu = FALSE) {
  if (!mlu) return(vapply(seq_len(n_timesteps), function(t) sum(x[,t,,drop=FALSE]) / n_perm, numeric(1)))
  vapply(seq_len(n_timesteps), function(t) sum(x[,,t,,drop=FALSE]) / n_perm, numeric(1))
}
assert_profile <- function(actual, expected, label, tol=1e-12) {
  if (length(actual) != length(expected) || any(abs(actual-expected) > tol))
    stop(label, ': expected ', paste(expected, collapse=', '), ' but got ', paste(round(actual,3), collapse=', '))
  TRUE
}
expect_error <- function(expr, pattern, label) {
  msg <- tryCatch({ force(expr); NA_character_ }, error=function(e) conditionMessage(e))
  if (is.na(msg)) stop(label, ': expected an error but none was raised')
  if (!grepl(pattern, msg, fixed=TRUE)) stop(label, ': wrong error: ', msg)
  TRUE
}

# Test-only shims for unavailable packages/PSOCK execution in WebR.
abind_test <- function(..., along) {
  xs <- list(...); d <- dim(xs[[1]])
  if (is.null(d) || along != length(d)+1L) stop('test abind shim only supports a new final dimension')
  out <- array(NA, dim=c(d,length(xs))); out[] <- unlist(xs,use.names=FALSE); out
}
foreach_test <- function(..., .combine, .packages=character()) {
  dots <- as.list(substitute(list(...)))[-1L]
  list(iter_expr=dots[[1L]], combine=.combine, env=parent.frame())
}
dopar_test <- function(obj, expr) {
  block <- substitute(expr); vals <- eval(obj$iter_expr,obj$env)
  ans <- lapply(vals,function(v) eval(block,envir=new.env(parent=obj$env)))
  comb <- get(obj$combine,envir=obj$env,inherits=TRUE); do.call(comb,ans)
}
ina_test <- function(..., initinfo, initbio, bpam, probestabvec, probadoptvec) {
  incoming <- as.numeric(as.vector(initbio) %*% bpam)
  estab <- as.logical((as.vector(initbio)>0 & as.vector(probestabvec)>0) | incoming>0)
  step <- list(vect1cL=list(NULL,as.vector(initinfo)),estabvecL=list(estab))
  list(multdetails=list(list(multout=list(step))))
}
load_model <- function(code_dir, filename, fun_name, parallel=FALSE, ina=FALSE, transition_parallel=FALSE) {
  env <- new.env(parent=.GlobalEnv)
  if (parallel) {
    env$abind <- abind_test; env$foreach <- foreach_test; env$`%dopar%` <- dopar_test
    env$detectCores <- function(...) 2L; env$makeCluster <- function(...) structure(list(),class='testCluster')
    env$registerDoParallel <- function(...) invisible(NULL); env$stopCluster <- function(...) invisible(NULL)
    env$clusterExport <- function(...) invisible(NULL); env$clusterEvalQ <- function(...) invisible(NULL)
    env$parLapply <- function(cl,X,fun,...) lapply(X,fun,...)
    if (ina) env$INAscene <- ina_test
    lines <- readLines(file.path(code_dir,filename),warn=FALSE)
    lines <- lines[!grepl('^[[:space:]]*library\\((abind|doParallel|INA)\\)[[:space:]]*$',lines)]
    if (transition_parallel) lines <- gsub('parallel::stopCluster(cl)','stopCluster(cl)',lines,fixed=TRUE)
    eval(parse(text=lines,keep.source=FALSE),envir=env)
  } else sys.source(file.path(code_dir,filename),envir=env)
  get(fun_name,envir=env,inherits=FALSE)
}

basic_args <- function(eradication=1, detection=0) list(
  Nperm=n_perm,Ntimesteps=n_timesteps,DetectionProb=detection,DetectionSD=0,
  ManageProb=1,ManageSD=0,EradicationProb=eradication,EradicationSD=0,
  SpreadReduction=0,SpreadReductionSD=0,InitialInvasion=all_binary,InitBioP=NA,
  InvasionRisk=zero_external,InitialInfo=all_info,InitInfoP=NA,ExternalInfoProb=zero_external,
  EnvEstabProb=1,Survival=1,SDDprob=zero_route,SEAM=0,LDDprob=zero_route,
  OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,DoPlots=FALSE)

meta_args <- function(mortality=1, detection=0) list(
  Nperm=n_perm,Ntimesteps=n_timesteps,DetectionProb=detection,DetectionSD=0,
  ManageProb=1,ManageSD=0,MortalityProb=mortality,MortalitySD=0,
  SpreadReduction=0,SpreadReductionSD=0,InitialPopulation=all_population,InitBioP=NA,
  InvasionRisk=zero_external,InitialInfo=all_info,InitInfoP=NA,ExternalInfoProb=zero_external,
  EnvEstabProb=1,Survival=1,K=rep(100,n_nodes),PropaguleProduction=0,PropaguleEstablishment=1,
  IncursionStartPop=NA,SDDprob=zero_route,SEAM=0,LDDprob=zero_route,LDDrate=0,
  OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,DoPlots=FALSE)

mlu_args <- function(mortality=c(1,1), detection=c(0,0), mortality_sd=c(0,0), detection_sd=c(0,0)) list(
  Nperm=n_perm,Ntimesteps=n_timesteps,Nlanduses=n_landuses,DetectionProb=detection,DetectionSD=detection_sd,
  ManageProb=c(1,1),ManageSD=c(0,0),MortalityProb=mortality,MortalitySD=mortality_sd,
  SpreadReduction=c(0,0),SpreadReductionSD=c(0,0),InitialPopulation=all_population_lu,InitBioP=NA,
  InvasionRisk=zero_external,InitialInfo=all_info,InitInfoP=NA,ExternalInfoProb=zero_external,
  EnvEstabProb=1,Survival=c(1,1),K=matrix(100,nrow=n_nodes,ncol=n_landuses),PropaguleProduction=0,
  PropaguleEstablishment=1,IncursionStartPop=NA,SDDprob=zero_route,SEAM=0,LDDprob=zero_route,LDDrate=0,
  OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,DoPlots=FALSE)

transition_args <- function(mortality=c(1,1), detection=c(0,0), mortality_sd=0, detection_sd=0) list(
  Nperm=n_perm,Ntimesteps=n_timesteps,Nstages=2L,Weights=c(1,1),Transition=transition_matrix,
  DetectionProb=detection,DetectionSD=detection_sd,ManageProb=1,ManageSD=0,MortalityProb=mortality,
  MortalitySD=mortality_sd,SpreadReduction=0,SpreadReductionSD=0,InitialPopulation=all_population_stage,
  InitBioP=NA,InvasionRisk=zero_external,InitialInfo=all_info,InitInfoP=NA,ExternalInfoProb=zero_external,
  EnvEstabProb=1,K=rep(100,n_nodes),SeedbankK=rep(100,n_nodes),PropaguleEstablishment=1,
  IncursionStartPop=NA,SDDprob=zero_route,SEAM=0,LDDprob=zero_route,LDDrate=0,
  DispersalDensityFactor=0,BlockedTransitionMortality=0,OngoingExternalInvasion=FALSE,
  OngoingExternalInfo=FALSE,DoPlots=FALSE)

model_specs <- list(
  list(label='INApest',file='INApest.R',fun='INApest',builder=basic_args,suffixes=c('InfoLargeOut.rds','InvasionLargeOut.rds','DetectedLargeOut.rds'),manage_suffix='InfoLargeOut.rds',mlu=FALSE,type='basic'),
  list(label='INApestParallel',file='INApestParallel.R',fun='INApestParallel',builder=basic_args,suffixes=c('ManagingLargeOut.rds','InvasionLargeOut.rds','DetectedLargeOut.rds'),manage_suffix='ManagingLargeOut.rds',mlu=FALSE,type='basic',parallel=TRUE),
  list(label='INApestParallelINAscene',file='INApestParallelInAScene.R',fun='INApestParallelINAscene',builder=basic_args,suffixes=c('ManagingLargeOut.rds','InvasionLargeOut.rds','DetectedLargeOut.rds'),manage_suffix='ManagingLargeOut.rds',mlu=FALSE,type='basic',parallel=TRUE,ina=TRUE),
  list(label='INApestMetaParallel',file='INApestMetaParallel.r',fun='INApestMetaParallel',builder=meta_args,suffixes=c('InfoLargeOut.rds','PopulationLargeOut.rds','InvasionLargeOut.rds','DetectedLargeOut.rds'),manage_suffix='InfoLargeOut.rds',mlu=FALSE,type='meta',parallel=TRUE),
  list(label='INApestMetaMultipleLandUse',file='INApestMetaMultipleLandUse.r',fun='INApestMetaMultipleLandUse',builder=mlu_args,suffixes=c('InfoLargeOut.rds','PopulationLargeOut.rds','InvasionLargeOut.rds','DetectedLargeOut.rds'),manage_suffix='InfoLargeOut.rds',mlu=TRUE,type='mlu'),
  list(label='INApestMetaParallelMultipleLandUse',file='INApestMetaParallelMultipleLandUse.r',fun='INApestMetaParallelMultipleLandUse',builder=mlu_args,suffixes=c('InfoLargeOut.rds','PopulationLargeOut.rds','InvasionLargeOut.rds','DetectedLargeOut.rds'),manage_suffix='InfoLargeOut.rds',mlu=TRUE,type='mlu',parallel=TRUE),
  list(label='INApestMetaTransitionMatrix',file='INApestMetaTransitionMatrix.r',fun='INApestMetaTransitionMatrix',builder=transition_args,suffixes=c('InfoLargeOut.rds','PopulationLargeOut.rds','PopulationStageLargeOut.rds','InvasionLargeOut.rds','DetectedLargeOut.rds'),manage_suffix='InfoLargeOut.rds',mlu=FALSE,type='transition'),
  list(label='INApestMetaTransitionMatrixParallel',file='INApestMetaTransitionMatrixParallel.r',fun='INApestMetaTransitionMatrixParallel',builder=transition_args,suffixes=c('InfoLargeOut.rds','PopulationLargeOut.rds','PopulationStageLargeOut.rds','InvasionLargeOut.rds','DetectedLargeOut.rds'),manage_suffix='InfoLargeOut.rds',mlu=FALSE,type='transition',parallel=TRUE,transition_parallel=TRUE)
)

# Time-varying evidence-reset arguments.
make_reset_args <- function(spec) {
  if (spec$type == 'basic') {
    erad <- matrix(0,nrow=n_nodes,ncol=n_timesteps); erad[,2] <- 1
    return(spec$builder(eradication=erad,detection=0))
  }
  if (spec$type == 'meta') {
    mort <- matrix(0,nrow=n_nodes,ncol=n_timesteps); mort[,2] <- 1
    det <- matrix(0,nrow=n_nodes,ncol=n_timesteps); det[,1] <- 1
    return(spec$builder(mortality=mort,detection=det))
  }
  if (spec$type == 'mlu') {
    mort <- array(0,dim=c(n_nodes,n_landuses,n_timesteps)); mort[,,2] <- 1
    det <- array(0,dim=c(n_nodes,n_landuses,n_timesteps)); det[,,1] <- 1
    return(spec$builder(mortality=mort,detection=det,mortality_sd=array(0,dim=dim(mort)),detection_sd=array(0,dim=dim(det))))
  }
  mort <- array(0,dim=c(n_nodes,2L,n_timesteps)); mort[,,2] <- 1
  det <- array(0,dim=c(n_nodes,2L,n_timesteps)); det[,,1] <- 1
  spec$builder(mortality=mort,detection=det,mortality_sd=array(0,dim=dim(mort)),detection_sd=array(0,dim=dim(det)))
}

make_no_kill_args <- function(spec) {
  if (spec$type == 'meta') {
    det <- matrix(0,nrow=n_nodes,ncol=n_timesteps); det[,1] <- 1
    return(spec$builder(mortality=0,detection=det))
  }
  if (spec$type == 'mlu') {
    det <- array(0,dim=c(n_nodes,n_landuses,n_timesteps)); det[,,1] <- 1
    return(spec$builder(mortality=c(0,0),detection=det,detection_sd=array(0,dim=dim(det))))
  }
  if (spec$type == 'transition') {
    det <- array(0,dim=c(n_nodes,2L,n_timesteps)); det[,,1] <- 1
    return(spec$builder(mortality=c(0,0),detection=det,detection_sd=array(0,dim=dim(det))))
  }
  spec$builder(eradication=0,detection=0)
}


# Dedicated test: information from another node refreshes a locally expired
# programmed-stop state. Node 1 is donor; node 2 is recipient.
modified_dir <- file.path(bundle_root, 'modified')
seam <- matrix(0, nrow=n_nodes, ncol=n_nodes)
seam[1,2] <- 1

recipient_profile <- function(x, mlu=FALSE, node=2L) {
  if (!mlu) return(vapply(seq_len(n_timesteps), function(t) sum(x[node,t,,drop=FALSE]) / n_perm, numeric(1)))
  vapply(seq_len(n_timesteps), function(t) sum(x[node,,t,,drop=FALSE]) / n_perm, numeric(1))
}

make_social_args <- function(spec, with_seam=TRUE) {
  if (spec$type == 'basic') {
    a <- spec$builder(eradication=0,detection=1)
    a$InitialInvasion <- rep(0L,n_nodes); a$InitialInvasion[1] <- 1L
  } else if (spec$type == 'meta') {
    a <- spec$builder(mortality=0,detection=1)
    a$InitialPopulation <- rep(0,n_nodes); a$InitialPopulation[1] <- 20
  } else if (spec$type == 'mlu') {
    a <- spec$builder(mortality=c(0,0),detection=c(1,1))
    a$InitialPopulation <- matrix(0,nrow=n_nodes,ncol=n_landuses); a$InitialPopulation[1,] <- 20
  } else {
    a <- spec$builder(mortality=c(0,0),detection=c(1,1))
    a$InitialPopulation <- matrix(0,nrow=n_nodes,ncol=2L); a$InitialPopulation[1,2] <- 20
  }
  a$InitialInfo <- rep(0,n_nodes); a$InitialInfo[c(1,2)] <- 1
  a$InfoPersistenceSteps <- 0
  a$InfoRetentionProb <- 1
  a$SEAM <- if (with_seam) seam else 0
  if (isTRUE(spec$ina)) a$geocoords <- geocoords
  a
}

rows <- list(); progress_file <- file.path(bundle_root,'social_refresh_progress.txt'); unlink(progress_file)
mark <- function(x) cat(x,'\n',file=progress_file,append=TRUE)
for (spec in model_specs) {
  label <- spec$label; cat('\n===',label,'===\n'); mark(paste('START',label))
  fun <- load_model(modified_dir,spec$file,spec$fun,parallel=isTRUE(spec$parallel),ina=isTRUE(spec$ina),transition_parallel=isTRUE(spec$transition_parallel))
  expected_on <- rep(if(spec$mlu) n_landuses else 1, n_timesteps)
  expected_off <- c(if(spec$mlu) n_landuses else 1, rep(0,n_timesteps-1L))

  on <- run_and_read(fun,make_social_args(spec,TRUE),paste0(label,'_social_on_'),spec$suffixes,71717L)
  off <- run_and_read(fun,make_social_args(spec,FALSE),paste0(label,'_social_off_'),spec$suffixes,71717L)
  p_on <- recipient_profile(on[[spec$manage_suffix]],spec$mlu)
  p_off <- recipient_profile(off[[spec$manage_suffix]],spec$mlu)
  assert_profile(p_on,expected_on,paste0(label,' cross-node refresh'))
  assert_profile(p_off,expected_off,paste0(label,' no cross-node refresh'))
  rows[[length(rows)+1L]] <- data.frame(function_name=label,
    seam_refresh=paste(p_on,collapse=' -> '), no_seam=paste(p_off,collapse=' -> '),
    stringsAsFactors=FALSE)
  mark(paste('PASS',label)); cat('SEAM refresh:',paste(p_on,collapse=' -> '),' | no SEAM:',paste(p_off,collapse=' -> '),'\n')
}
res <- do.call(rbind,rows)
write.csv(res,file.path(bundle_root,'social_refresh_results.csv'),row.names=FALSE)
cat('\nALL SOCIAL REFRESH TESTS PASSED\n')
print(res)
