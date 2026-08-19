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

model_specs <- model_specs[5:8]
summary_rows <- list(); progress_file <- file.path(bundle_root,'info_persistence_progress_g2.txt'); unlink(progress_file)
mark <- function(x) cat(x,'\n',file=progress_file,append=TRUE)

for (spec in model_specs) {
  label <- spec$label; cat('\n===',label,'===\n'); mark(paste('START',label))
  is_parallel <- isTRUE(spec$parallel)
  old_fun <- load_model(baseline_dir,spec$file,spec$fun,parallel=is_parallel,ina=isTRUE(spec$ina),transition_parallel=isTRUE(spec$transition_parallel))
  new_fun <- load_model(modified_dir,spec$file,spec$fun,parallel=is_parallel,ina=isTRUE(spec$ina),transition_parallel=isTRUE(spec$transition_parallel))
  prep <- function(a) { if (isTRUE(spec$ina)) a$geocoords <- geocoords; a }
  max_units <- n_nodes * if(spec$mlu) n_landuses else 1L; half_units <- max_units/2

  # 1) Exact backward compatibility when the new rule is omitted (NA default).
  base_args <- prep(spec$builder())
  old_out <- run_and_read(old_fun,base_args,paste0(label,'_baseline_'),spec$suffixes,10101L)
  new_out <- run_and_read(new_fun,base_args,paste0(label,'_default_'),spec$suffixes,10101L)
  assert_identical_outputs(old_out,new_out,paste0(label,' default compatibility'))
  p_default <- management_profile(new_out[[spec$manage_suffix]],spec$mlu)
  mark(paste('DEFAULT',label))

  # 2) Scalar programmed stop = 2 steps after the t1 known presence/kill.
  a <- prep(spec$builder()); a$InfoPersistenceSteps <- 2
  o <- run_and_read(new_fun,a,paste0(label,'_steps2_'),spec$suffixes,20202L)
  p_steps2 <- management_profile(o[[spec$manage_suffix]],spec$mlu)
  assert_profile(p_steps2,c(max_units,max_units,max_units,0),paste0(label,' scalar persistence'))
  mark(paste('STEPS2',label))

  # 3) Zero extra steps: management ends after the timestep with last known presence.
  a <- prep(spec$builder()); a$InfoPersistenceSteps <- 0
  o <- run_and_read(new_fun,a,paste0(label,'_steps0_'),spec$suffixes,30303L)
  p_steps0 <- management_profile(o[[spec$manage_suffix]],spec$mlu)
  assert_profile(p_steps0,c(max_units,0,0,0),paste0(label,' zero persistence'))
  mark(paste('STEPS0',label))

  # 4) Node vector: half persist 2 steps, half stop immediately.
  node_steps <- c(rep(2,n_nodes/2L),rep(0,n_nodes/2L))
  a <- prep(spec$builder()); a$InfoPersistenceSteps <- node_steps
  o <- run_and_read(new_fun,a,paste0(label,'_node_'),spec$suffixes,40404L)
  p_node <- management_profile(o[[spec$manage_suffix]],spec$mlu)
  assert_profile(p_node,c(max_units,half_units,half_units,0),paste0(label,' node persistence'))
  mark(paste('NODE',label))

  # 5) Node x timestep policy: at t2 first half changes to immediate stop.
  time_steps <- matrix(2,nrow=n_nodes,ncol=n_timesteps); time_steps[seq_len(n_nodes/2L),2] <- 0
  a <- prep(spec$builder()); a$InfoPersistenceSteps <- time_steps
  o <- run_and_read(new_fun,a,paste0(label,'_nodetime_'),spec$suffixes,50505L)
  p_time <- management_profile(o[[spec$manage_suffix]],spec$mlu)
  assert_profile(p_time,c(max_units,max_units,half_units,0),paste0(label,' node x timestep persistence'))
  mark(paste('NODETIME',label))

  # 6) Programmed stop takes priority over stochastic decay, with one warning.
  a <- prep(spec$builder()); a$InfoPersistenceSteps <- 2; a$InfoRetentionProb <- 0
  wr <- run_and_read_warnings(new_fun,a,paste0(label,'_priority_'),spec$suffixes,60606L)
  p_priority <- management_profile(wr$out[[spec$manage_suffix]],spec$mlu)
  assert_profile(p_priority,c(max_units,max_units,max_units,0),paste0(label,' priority over decay'))
  priority_warnings <- wr$warnings[grepl('Programmed stopping takes priority',wr$warnings,fixed=TRUE)]
  if (length(priority_warnings) != 1L) stop(label, ': expected exactly one priority warning; got ',length(priority_warnings))
  mark(paste('PRIORITY',label))

  # 7) Mixed mode: programmed stop on first half, stochastic decay=0 on second half.
  mixed_steps <- c(rep(2,n_nodes/2L),rep(NA,n_nodes/2L))
  a <- prep(spec$builder()); a$InfoPersistenceSteps <- mixed_steps; a$InfoRetentionProb <- 0
  o <- suppressWarnings(run_and_read(new_fun,a,paste0(label,'_mixed_'),spec$suffixes,70707L))
  p_mixed <- management_profile(o[[spec$manage_suffix]],spec$mlu)
  assert_profile(p_mixed,c(max_units,half_units,half_units,0),paste0(label,' mixed programmed/decay'))
  mark(paste('MIXED',label))

  # 8) Evidence resets the clock. Occupancy models use known extant presence;
  # Meta models specifically require a management kill at t2 to extend through t3.
  reset_args <- prep(make_reset_args(spec)); reset_args$InfoPersistenceSteps <- 1
  o <- tryCatch(run_and_read(new_fun,reset_args,paste0(label,'_reset_'),spec$suffixes,80808L), error=function(e) { mark(paste('RESET_ERROR',label,conditionMessage(e))); stop(e) })
  p_reset <- management_profile(o[[spec$manage_suffix]],spec$mlu)
  assert_profile(p_reset,c(max_units,max_units,max_units,0),paste0(label,' evidence reset'))
  mark(paste('RESET',label))

  if (spec$type != 'basic') {
    no_kill_args <- prep(make_no_kill_args(spec)); no_kill_args$InfoPersistenceSteps <- 1
    o0 <- run_and_read(new_fun,no_kill_args,paste0(label,'_nokill_'),spec$suffixes,80808L)
    p_no_kill <- management_profile(o0[[spec$manage_suffix]],spec$mlu)
    assert_profile(p_no_kill,c(max_units,max_units,0,0),paste0(label,' no management-kill reset'))
  } else {
    persist_args <- prep(make_no_kill_args(spec)); persist_args$InfoPersistenceSteps <- 1
    op <- run_and_read(new_fun,persist_args,paste0(label,'_presence_'),spec$suffixes,90909L)
    p_no_kill <- management_profile(op[[spec$manage_suffix]],spec$mlu)
    assert_profile(p_no_kill,c(max_units,max_units,max_units,max_units),paste0(label,' extant presence reset'))
  }
  mark(paste('COMPARE',label))

  # 9) Validation.
  bad <- prep(spec$builder()); bad$InfoPersistenceSteps <- c(1,1); bad$ModelName <- paste0(label,'_badlen_'); bad$OutputDir <- out_dir
  expect_error(do.call(new_fun,bad),'InfoPersistenceSteps must be a single value, vector of length nodes, or matrix nodes x Ntimesteps',paste0(label,' invalid persistence length'))
  bad <- prep(spec$builder()); bad$InfoPersistenceSteps <- 1.5; bad$ModelName <- paste0(label,'_badwhole_'); bad$OutputDir <- out_dir
  expect_error(do.call(new_fun,bad),'InfoPersistenceSteps values must be non-negative whole numbers or NA',paste0(label,' non-whole persistence'))
  mark(paste('VALIDATION',label))

  mark(paste('PASS',label)); cat('PASS',label,'\n')
  cat('  steps=2: ',paste(p_steps2,collapse=' -> '),'\n')
  cat('  steps=0: ',paste(p_steps0,collapse=' -> '),'\n')
  cat('  node:    ',paste(p_node,collapse=' -> '),'\n')
  cat('  node-time:',paste(p_time,collapse=' -> '),'\n')
  cat('  priority:',paste(p_priority,collapse=' -> '),' warning=',length(priority_warnings),'\n')
  cat('  mixed:   ',paste(p_mixed,collapse=' -> '),'\n')
  cat('  reset:   ',paste(p_reset,collapse=' -> '),'\n')
  cat('  compare: ',paste(p_no_kill,collapse=' -> '),'\n')

  summary_rows[[length(summary_rows)+1L]] <- data.frame(
    Function=label,DefaultCompatible=TRUE,
    Steps2_t1=p_steps2[1],Steps2_t2=p_steps2[2],Steps2_t3=p_steps2[3],Steps2_t4=p_steps2[4],
    Steps0_t1=p_steps0[1],Steps0_t2=p_steps0[2],Steps0_t3=p_steps0[3],Steps0_t4=p_steps0[4],
    Node_t1=p_node[1],Node_t2=p_node[2],Node_t3=p_node[3],Node_t4=p_node[4],
    NodeTime_t1=p_time[1],NodeTime_t2=p_time[2],NodeTime_t3=p_time[3],NodeTime_t4=p_time[4],
    Priority_t1=p_priority[1],Priority_t2=p_priority[2],Priority_t3=p_priority[3],Priority_t4=p_priority[4],PriorityWarningCount=length(priority_warnings),
    Mixed_t1=p_mixed[1],Mixed_t2=p_mixed[2],Mixed_t3=p_mixed[3],Mixed_t4=p_mixed[4],
    Reset_t1=p_reset[1],Reset_t2=p_reset[2],Reset_t3=p_reset[3],Reset_t4=p_reset[4],
    Compare_t1=p_no_kill[1],Compare_t2=p_no_kill[2],Compare_t3=p_no_kill[3],Compare_t4=p_no_kill[4],
    InvalidLengthRejected=TRUE,NonWholeRejected=TRUE,stringsAsFactors=FALSE)
}

summary_df <- do.call(rbind,summary_rows)
write.csv(summary_df,file.path(bundle_root,'info_persistence_results_g2.csv'),row.names=FALSE)
writeLines(capture.output(print(summary_df,row.names=FALSE)),file.path(bundle_root,'info_persistence_results_g2.txt'))
cat('\nALL INFORMATION-PERSISTENCE MODEL TESTS PASSED\n')
