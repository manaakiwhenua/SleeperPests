###############################################################################
### INApest PoA Layer 2 -- exact internal behavioural validation
###
### Purpose
### -------
### Validate the PoA inference mathematics against independently specified
### closed-form / finite-state results.  These tests deliberately use synthetic
### model outputs so ecological-engine complexity cannot hide inference errors.
###
### Sequence
###  L2-01 one surveillance round
###  L2-02 repeated zero-detection rounds
###  L2-03 reinvasion between rounds
###  L2-04 extinction + reinvasion between rounds
###  L2-05 Background + InfoTriggered surveillance
###  L2-06 observation-compatible history / management feedback
###  L2-07 high-PoA discretisation convergence
###  L2-08 unsupported posterior-class propagation safeguard
###  L2-09 simulation-mode Monte Carlo benchmark
###  L2-10 positive-detection binary likelihood
###############################################################################

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1L], mustWork=TRUE) else normalizePath(getwd(), mustWork=TRUE)
source_file <- file.path(root, "source", "INApestPoA.R")
if(!file.exists(source_file)) source_file <- file.path(root, "..", "..", "definitive_source", "INApestPoA.R")
if(!file.exists(source_file)) stop("PoA source not found for validation root: ", root)
source(source_file, local=.GlobalEnv)

out_dir <- file.path(root, "validation_output")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

results <- list(); metrics <- list(); convergence <- list(); kk <- 0L; mm <- 0L; cc <- 0L

assert <- function(x, msg="assertion failed") if(!isTRUE(x)) stop(msg, call.=FALSE)
assert_close <- function(x, y, tol=1e-12, msg=NULL) {
  if(length(x)!=1L || length(y)!=1L || !is.finite(x) || !is.finite(y) || abs(x-y)>tol)
    stop(if(is.null(msg)) sprintf("expected %.15g; got %.15g; |error|=%.6g", y, x, abs(x-y)) else msg, call.=FALSE)
}
expect_error <- function(expr, pattern=NULL) {
  got <- NULL
  tryCatch(force(expr), error=function(e) got <<- conditionMessage(e))
  if(is.null(got)) stop("expected error but expression succeeded", call.=FALSE)
  if(!is.null(pattern) && !grepl(pattern, got, fixed=TRUE))
    stop("wrong error: ", got, call.=FALSE)
  got
}
record_metric <- function(test_id, quantity, observed, expected, tolerance=1e-12, note="") {
  mm <<- mm + 1L
  metrics[[mm]] <<- data.frame(
    test_id=test_id, quantity=quantity, observed=as.numeric(observed), expected=as.numeric(expected),
    abs_error=abs(as.numeric(observed)-as.numeric(expected)), tolerance=tolerance,
    pass=isTRUE(abs(as.numeric(observed)-as.numeric(expected)) <= tolerance), note=note,
    stringsAsFactors=FALSE
  )
  assert_close(as.numeric(observed), as.numeric(expected), tolerance,
               paste0(test_id, " ", quantity, " mismatch"))
  invisible(observed)
}
run_test <- function(id, description, fun) {
  kk <<- kk + 1L
  t0 <- proc.time()[[3L]]; status <- "PASS"; detail <- ""
  tryCatch({ detail <- as.character(fun()) }, error=function(e) {
    status <<- "FAIL"; detail <<- conditionMessage(e)
  })
  results[[kk]] <<- data.frame(test_id=id, description=description, status=status,
    detail=detail, elapsed_seconds=proc.time()[[3L]]-t0, stringsAsFactors=FALSE)
  cat(sprintf("%-7s %-5s %s\n", id, status, description))
}

as_node_array <- function(x) {
  x <- as.matrix(x)
  np <- nrow(x); nt <- ncol(x)
  a <- array(0, dim=c(1L, nt, np))
  a[1,,] <- t(x)
  a
}
resolve_particle_time <- function(x, np, nt, name) {
  if(length(x)==1L) return(matrix(as.numeric(x), nrow=np, ncol=nt))
  x <- as.matrix(x)
  if(!identical(dim(x), c(np,nt))) stop(name, " must be scalar or permutations x timesteps")
  x
}
make_binary_trajectories <- function(present, bg_events=0, info_events=0, info_before=0,
                                     p_bg=0, p_info=0, model_name="layer2_binary") {
  present <- as.matrix(present); np <- nrow(present); nt <- ncol(present)
  bg_events <- resolve_particle_time(bg_events,np,nt,"bg_events")
  info_events <- resolve_particle_time(info_events,np,nt,"info_events")
  info_before <- resolve_particle_time(info_before,np,nt,"info_before")
  p_bg <- resolve_particle_time(p_bg,np,nt,"p_bg")
  p_info <- resolve_particle_time(p_info,np,nt,"p_info")
  list(
    ModelName=model_name,
    InvasionResults=as_node_array(present),
    BackgroundDetectedResults=as_node_array(bg_events),
    InfoTriggeredDetectedResults=as_node_array(info_events),
    DetectedResults=as_node_array(bg_events),
    InformationStateBeforeSurveillanceResults=as_node_array(info_before),
    BackgroundDetectionProbabilityResults=as_node_array(p_bg),
    InfoTriggeredDetectionProbabilityResults=as_node_array(p_info)
  )
}
make_population_one_round <- function(abundance, bg_events=0, p_bg=.2, model_name="layer2_population") {
  abundance <- as.numeric(abundance); np <- length(abundance)
  z <- matrix(abundance,nrow=np,ncol=1L)
  make <- function(v) as_node_array(matrix(v,nrow=np,ncol=1L))
  list(
    ModelName=model_name,
    PopulationResults=make(abundance),
    BackgroundDetectedResults=make(rep(bg_events,length.out=np)),
    InfoTriggeredDetectedResults=make(rep(0,np)),
    DetectedResults=make(rep(bg_events,length.out=np)),
    InformationStateBeforeSurveillanceResults=make(rep(0,np)),
    BackgroundDetectionProbabilityResults=make(rep(p_bg,length.out=np)),
    InfoTriggeredDetectionProbabilityResults=make(rep(0,np))
  )
}
obs_bg_zero <- function(nt) data.frame(Timestep=seq_len(nt), BackgroundDetections=0)
round_row <- function(z, r) z$PoASummary[z$PoASummary$Round==r,,drop=FALSE]

# -----------------------------------------------------------------------------
# L2-01. Single surveillance round: closed-form two-state Bayes
# -----------------------------------------------------------------------------
run_test("L2-01", "single zero-detection round matches closed-form Bayes", function() {
  x <- make_binary_trajectories(matrix(c(0,1),ncol=1), p_bg=.8)
  z <- INApestPoACore(x,.INApestPoAAdapterBinaryNode,PriorPoA=.8,
    ObservationHistory=data.frame(Timestep=1,BackgroundDetections=0),
    EvidenceSources="Background",PoAMethod="likelihood")
  expected <- .8/(.8+.2*.2)
  rr <- round_row(z,1)
  record_metric("L2-01","posterior_poa",rr$PosteriorPoA,expected)
  record_metric("L2-01","background_sse",rr$BackgroundSSe,.8)
  "posterior 0.952380952381; SSe 0.8"
})

# -----------------------------------------------------------------------------
# L2-02. Repeated zero-detection rounds
# -----------------------------------------------------------------------------
run_test("L2-02", "three repeated zero-detection rounds reproduce recursive Bayes", function() {
  x <- make_binary_trajectories(rbind(c(0,0,0),c(1,1,1)), p_bg=.8)
  z <- INApestPoACore(x,.INApestPoAAdapterBinaryNode,PriorPoA=.8,
    ObservationHistory=obs_bg_zero(3),EvidenceSources="Background",PoAMethod="likelihood")
  for(r in 1:3) {
    expected <- .8/(.8+.2*(.2^r))
    record_metric("L2-02",paste0("posterior_round_",r),round_row(z,r)$PosteriorPoA,expected)
  }
  "recursive posterior sequence reproduced"
})

# -----------------------------------------------------------------------------
# L2-03. Reinvasion between rounds
# -----------------------------------------------------------------------------
run_test("L2-03", "reinvasion between rounds matches exact finite-state update", function() {
  p1 <- c(rep(0,10),rep(1,10))
  p2 <- c(rep(0,9),1,rep(1,10))  # 10% of the t1-absent class reinvades
  x <- make_binary_trajectories(cbind(p1,p2), p_bg=.8)
  z <- INApestPoACore(x,.INApestPoAAdapterBinaryNode,PriorPoA=.8,
    ObservationHistory=obs_bg_zero(2),EvidenceSources="Background",PoAMethod="likelihood")
  post1 <- .8/(.8+.2*.2)
  prior2 <- post1*(1-.1)
  post2 <- prior2/(prior2+(1-prior2)*.2)
  record_metric("L2-03","posterior_round_1",round_row(z,1)$PosteriorPoA,post1)
  record_metric("L2-03","prior_round_2",round_row(z,2)$PriorPoA,prior2)
  record_metric("L2-03","posterior_round_2",round_row(z,2)$PosteriorPoA,post2)
  "10% reinvasion exact"
})

# -----------------------------------------------------------------------------
# L2-04. Extinction + reinvasion between rounds
# -----------------------------------------------------------------------------
run_test("L2-04", "extinction plus reinvasion matches exact finite-state transition", function() {
  p1 <- c(rep(0,10),rep(1,10))
  p2 <- c(rep(0,8),rep(1,2), rep(0,3),rep(1,7)) # r=0.2 from absent; e=0.3 from present
  x <- make_binary_trajectories(cbind(p1,p2), p_bg=.6)
  z <- INApestPoACore(x,.INApestPoAAdapterBinaryNode,PriorPoA=.7,
    ObservationHistory=obs_bg_zero(2),EvidenceSources="Background",PoAMethod="likelihood")
  post1 <- .7/(.7+.3*.4)
  prior2 <- post1*(1-.2)+(1-post1)*.3
  post2 <- prior2/(prior2+(1-prior2)*.4)
  record_metric("L2-04","posterior_round_1",round_row(z,1)$PosteriorPoA,post1)
  record_metric("L2-04","prior_round_2",round_row(z,2)$PriorPoA,prior2)
  record_metric("L2-04","posterior_round_2",round_row(z,2)$PosteriorPoA,post2)
  "reinvasion 0.2 + extinction 0.3 exact"
})

# -----------------------------------------------------------------------------
# L2-05. Background + information-triggered surveillance
# -----------------------------------------------------------------------------
run_test("L2-05", "two surveillance pathways retain separate and combined exact SSe", function() {
  present <- matrix(c(0,0,1,1),ncol=1)
  info <- matrix(c(0,0,0,1),ncol=1) # only one of two present trajectories was already informed
  x <- make_binary_trajectories(present,info_before=info,p_bg=.6,p_info=.5)
  z <- INApestPoACore(x,.INApestPoAAdapterBinaryNode,PriorPoA=.5,
    ObservationHistory=data.frame(Timestep=1,BackgroundDetections=0,InfoTriggeredDetections=0),
    EvidenceSources=c("Background","InfoTriggered"),PoAMethod="likelihood")
  rr <- round_row(z,1)
  qcomb <- mean(c(.4,.4*.5))
  expected_post <- .5/(.5+.5*qcomb)
  record_metric("L2-05","background_sse",rr$BackgroundSSe,.6)
  record_metric("L2-05","info_triggered_sse",rr$InfoTriggeredSSe,.25)
  record_metric("L2-05","combined_sse",rr$CombinedSSe,.7)
  record_metric("L2-05","posterior_poa",rr$PosteriorPoA,expected_post)
  "Background SSe 0.6; targeted SSe 0.25; combined SSe 0.7"
})

# -----------------------------------------------------------------------------
# L2-06. Detection changes later management: compatible histories are essential
# -----------------------------------------------------------------------------
run_test("L2-06", "compatible-history propagation matches exact management-feedback benchmark", function() {
  # Particle 3 is present at t1, is actually detected, and is absent at t2
  # (standing in for successful detection-triggered management). Particle 4 is
  # present but not detected and remains present. Observed evidence is zero,
  # so particle 3 must NOT be allowed to contribute to t2 biology.
  present <- rbind(c(0,0),c(0,0),c(1,0),c(1,1))
  bg <- rbind(c(0,0),c(0,0),c(1,0),c(0,0))
  x <- make_binary_trajectories(present,bg_events=bg,p_bg=.5)
  z <- INApestPoACore(x,.INApestPoAAdapterBinaryNode,PriorPoA=.5,
    ObservationHistory=obs_bg_zero(2),EvidenceSources="Background",PoAMethod="likelihood")
  post1 <- 2/3
  exact_prior2 <- 2/3
  exact_post2 <- .8
  naive_prior2 <- 5/6
  naive_post2 <- 10/11
  record_metric("L2-06","posterior_round_1",round_row(z,1)$PosteriorPoA,post1)
  record_metric("L2-06","prior_round_2_compatible_history",round_row(z,2)$PriorPoA,exact_prior2)
  record_metric("L2-06","posterior_round_2",round_row(z,2)$PosteriorPoA,exact_post2)
  assert(abs(round_row(z,2)$PriorPoA-naive_prior2)>.1,"result is suspiciously close to naive propagation")
  assert(abs(round_row(z,2)$PosteriorPoA-naive_post2)>.05,"posterior is suspiciously close to naive propagation")
  paste0("exact posterior=0.8; naive all-history posterior=",format(naive_post2,digits=10)," correctly rejected")
})

# -----------------------------------------------------------------------------
# L2-07. High-PoA discretisation convergence and non-100% behaviour
# -----------------------------------------------------------------------------
run_test("L2-07", "high-PoA finite-particle approximation converges without false certainty", function() {
  alpha <- sqrt(2)/2
  q_exact <- alpha*.8 + (1-alpha)*(.8^20)
  exact <- .999/(.999+.001*q_exact)
  Ns <- c(1000L,10000L,50000L,200000L)
  last_error <- NA_real_
  for(N in Ns) {
    npresent <- max(2L,as.integer(floor(sqrt(N))))
    nlow <- as.integer(round(npresent*alpha))
    abundance <- c(rep(0,N-npresent),rep(1,nlow),rep(20,npresent-nlow))
    x <- make_population_one_round(abundance,p_bg=.2,model_name=paste0("high_poa_",N))
    z <- INApestPoACore(x,.INApestPoAAdapterPopulationNode,PriorPoA=.999,
      ObservationHistory=data.frame(Timestep=1,BackgroundDetections=0),
      EvidenceSources="Background",PoAMethod="likelihood")
    obs <- round_row(z,1)$PosteriorPoA
    err <- abs(obs-exact)
    cc <<- cc + 1L
    convergence[[cc]] <<- data.frame(Nparticles=N,Npresent=npresent,
      conditional_fraction_N1=nlow/npresent,posterior_poa=obs,exact_posterior_poa=exact,
      abs_error=err,stringsAsFactors=FALSE)
    assert(obs < 1,"high-PoA benchmark spuriously reported 100%")
    last_error <- err
  }
  assert(last_error < 5e-6,"200,000-particle discretisation did not converge closely enough")
  record_metric("L2-07","posterior_N200000",tail(convergence,1)[[1]]$posterior_poa,exact,5e-6,
                "high-PoA approximation")
  paste0("exact=",format(exact,digits=12),"; 200k error=",format(last_error,digits=6))
})

# -----------------------------------------------------------------------------
# L2-08. No compatible support for a positive posterior class must stop
# -----------------------------------------------------------------------------
run_test("L2-08", "unsupported posterior class stops instead of becoming certainty", function() {
  present <- rbind(c(0,0),c(0,0),c(1,1),c(1,1))
  bg <- rbind(c(0,0),c(0,0),c(1,0),c(1,0)) # every present t1 history conflicts with observed zero
  x <- make_binary_trajectories(present,bg_events=bg,p_bg=.5)
  msg <- expect_error(INApestPoACore(x,.INApestPoAAdapterBinaryNode,PriorPoA=.5,
    ObservationHistory=obs_bg_zero(2),EvidenceSources="Background",PoAMethod="likelihood"),
    "positive posterior mass for the present class")
  paste0("guard triggered: ",msg)
})

# -----------------------------------------------------------------------------
# L2-09. Simulation-mode conditioning converges to the same exact Bayes result
# -----------------------------------------------------------------------------
run_test("L2-09", "simulation-mode zero-detection conditioning agrees with exact Bayes", function() {
  N <- 100000L; nhalf <- N%/%2L
  present <- matrix(c(rep(0,nhalf),rep(1,N-nhalf)),ncol=1)
  set.seed(20260831)
  bg <- matrix(0,nrow=N,ncol=1)
  bg[(nhalf+1L):N,1] <- rbinom(N-nhalf,1,.8)
  x <- make_binary_trajectories(present,bg_events=bg,p_bg=.8)
  z <- INApestPoACore(x,.INApestPoAAdapterBinaryNode,PriorPoA=.8,
    ObservationHistory=data.frame(Timestep=1,BackgroundDetections=0),
    EvidenceSources="Background",PoAMethod="simulation")
  exact <- .8/(.8+.2*.2)
  observed <- round_row(z,1)$PosteriorPoA
  record_metric("L2-09","posterior_poa",observed,exact,.003,"100,000-particle simulation-mode tolerance")
  assert(round_row(z,1)$PosteriorPoA < 1,"simulation mode spuriously returned certainty")
  paste0("observed=",format(observed,digits=10),"; exact=",format(exact,digits=10))
})

# -----------------------------------------------------------------------------
# L2-10. Positive detection in binary likelihood mode makes absence impossible
# -----------------------------------------------------------------------------
run_test("L2-10", "positive binary detection drives exact absence posterior to zero", function() {
  x <- make_binary_trajectories(matrix(c(0,1),ncol=1),bg_events=matrix(c(0,1),ncol=1),p_bg=.8)
  z <- INApestPoACore(x,.INApestPoAAdapterBinaryNode,PriorPoA=.8,
    ObservationHistory=data.frame(Timestep=1,BackgroundDetections=1),
    EvidenceSources="Background",ObservationMatch="binary",PoAMethod="likelihood")
  rr <- round_row(z,1)
  record_metric("L2-10","posterior_poa",rr$PosteriorPoA,0)
  record_metric("L2-10","observation_probability",rr$ObservationProbability,.16)
  "positive detection posterior PoA = 0"
})

res <- do.call(rbind,results)
met <- if(length(metrics)) do.call(rbind,metrics) else data.frame()
conv <- if(length(convergence)) do.call(rbind,convergence) else data.frame()
utils::write.csv(res,file.path(out_dir,"layer2_test_results.csv"),row.names=FALSE)
utils::write.csv(met,file.path(out_dir,"layer2_numeric_comparisons.csv"),row.names=FALSE)
utils::write.csv(conv,file.path(out_dir,"layer2_high_poa_convergence.csv"),row.names=FALSE)

summary <- c(
  "INApest PoA Layer 2 exact internal behavioural validation",
  paste0("Date/time: ",format(Sys.time(),tz="",usetz=TRUE)),
  paste0("R: ",R.version.string),
  paste0("Platform: ",R.version$platform),
  paste0("Tests PASS: ",sum(res$status=="PASS")," / ",nrow(res)),
  paste0("Tests FAIL: ",sum(res$status=="FAIL")),
  if(nrow(met)) paste0("Maximum asserted numeric error: ",format(max(met$abs_error),digits=8)) else "",
  "",
  paste(res$test_id,res$status,res$description,res$detail,sep=" | ")
)
writeLines(summary,file.path(out_dir,"layer2_summary.txt"))
capture.output(sessionInfo(),file=file.path(out_dir,"sessionInfo.txt"))
cat("\n",paste(summary,collapse="\n"),"\n",sep="")
if(any(res$status=="FAIL")) quit(status=1L)
