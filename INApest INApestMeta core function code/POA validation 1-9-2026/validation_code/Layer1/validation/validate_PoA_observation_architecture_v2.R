###############################################################################
### INApest PoA observation architecture v2 validation -- base R only
###
### Tests the common observation/inference contract independently of the full
### biological engines.  This includes the canonical Background/InfoTriggered
### split, pre-surveillance information gating, realised detection-probability
### likelihoods, compatibility aliases and adapter coverage.
###############################################################################

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1L], mustWork=TRUE) else normalizePath(getwd(), mustWork=TRUE)
source_file <- file.path(root, "source", "INApestPoA.R")
if(!file.exists(source_file)) source_file <- file.path(root, "..", "..", "definitive_source", "INApestPoA.R")
if(!file.exists(source_file)) source_file <- file.path(root, "INApestPoA.R")
if(!file.exists(source_file)) stop("PoA source not found for validation root: ", root)
source(source_file, local=.GlobalEnv)

out_dir <- file.path(root, "validation_output", "observation_architecture")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

results <- list(); k <- 0L
run_test <- function(id, group, description, fun) {
  k <<- k + 1L; t0 <- proc.time()[[3L]]; status <- "PASS"; detail <- ""
  tryCatch({ detail <- as.character(fun()) }, error=function(e) {
    status <<- "FAIL"; detail <<- conditionMessage(e)
  })
  results[[k]] <<- data.frame(test_id=id,group=group,description=description,
    status=status,detail=detail,elapsed_seconds=proc.time()[[3L]]-t0,
    stringsAsFactors=FALSE)
  cat(sprintf("%-8s %-5s %s\n",id,status,description))
}
assert <- function(x,msg="assertion failed") if(!isTRUE(x)) stop(msg,call.=FALSE)
assert_close <- function(x,y,tol=1e-12,msg=NULL) {
  if(length(x)!=length(y) || any(!is.finite(c(x,y))) || max(abs(x-y))>tol)
    stop(if(is.null(msg)) paste0("expected ",paste(y,collapse=","),"; got ",paste(x,collapse=",")) else msg,call.=FALSE)
}
expect_error <- function(expr, pattern=NULL) {
  got <- NULL; tryCatch(force(expr),error=function(e) got <<- conditionMessage(e))
  if(is.null(got)) stop("expected error but expression succeeded",call.=FALSE)
  if(!is.null(pattern)&&!grepl(pattern,got,fixed=TRUE)) stop("wrong error: ",got,call.=FALSE)
  invisible(got)
}
expect_warning <- function(expr, pattern=NULL) {
  got <- character(0); value <- withCallingHandlers(expr,warning=function(w){
    got <<- c(got,conditionMessage(w)); invokeRestart("muffleWarning")
  })
  if(!length(got)) stop("expected warning but none was emitted",call.=FALSE)
  if(!is.null(pattern)&&!any(grepl(pattern,got,fixed=TRUE))) stop("wrong warning(s): ",paste(got,collapse=" | "),call.=FALSE)
  value
}

make_binary <- function(explicit=TRUE,info_event=TRUE,include_info_before=TRUE,
                        legacy_only=FALSE,with_prob=TRUE) {
  inv <- array(0,dim=c(1,2,4))
  inv[1,,1] <- c(0,0); inv[1,,2] <- c(1,0); inv[1,,3] <- c(1,1); inv[1,,4] <- c(1,1)
  out <- list(ModelName="synthetic_binary",InvasionResults=inv)
  if(explicit && !legacy_only) {
    bg <- array(0,dim=dim(inv)); bg[1,2,3] <- 1
    out$BackgroundDetectedResults <- bg
    if(info_event) { it <- array(0,dim=dim(inv)); it[1,2,4] <- 1; out$InfoTriggeredDetectedResults <- it }
    out$DetectedResults <- bg
  } else out$DetectedResults <- inv
  if(include_info_before) {
    h <- array(0,dim=dim(inv)); h[1,2,4] <- 1
    out$InformationStateBeforeSurveillanceResults <- h
  }
  if(with_prob) {
    out$BackgroundDetectionProbabilityResults <- array(.8,dim=dim(inv))
    out$InfoTriggeredDetectionProbabilityResults <- array(.5,dim=dim(inv))
  }
  out
}
make_meta <- function() {
  b <- make_binary(); b$PopulationResults <- b$InvasionResults*3; b$InvasionResults <- NULL; b$ModelName <- "synthetic_meta"; b
}
make_mlu <- function() {
  b <- make_binary(); x <- array(0,dim=c(1,2,2,4))
  x[1,1,,] <- b$InvasionResults[1,,]; x[1,2,,] <- b$InvasionResults[1,,]
  bg <- array(0,dim=c(1,2,4)); bg[1,2,3] <- 1
  it <- array(0,dim=c(1,2,4)); it[1,2,4] <- 1
  h <- array(0,dim=c(1,2,4)); h[1,2,4] <- 1
  pb <- array(.2,dim=dim(x)); pi <- array(.3,dim=dim(x))
  list(ModelName="synthetic_mlu",PopulationResults=x,
       BackgroundDetectedResults=bg,InfoTriggeredDetectedResults=it,DetectedResults=bg,
       InformationStateBeforeSurveillanceResults=h,
       BackgroundDetectionProbabilityResults=pb,InfoTriggeredDetectionProbabilityResults=pi)
}
make_tm <- function() {
  b <- make_binary(); x <- array(0,dim=c(1,2,2,4))
  x[1,1,,] <- b$InvasionResults[1,,]; x[1,2,,] <- b$InvasionResults[1,,]
  bg <- array(0,dim=c(1,2,4)); bg[1,2,3] <- 1
  it <- array(0,dim=c(1,2,4)); it[1,2,4] <- 1
  h <- array(0,dim=c(1,2,4)); h[1,2,4] <- 1
  pb <- array(.2,dim=dim(x)); pi <- array(.3,dim=dim(x))
  list(ModelName="synthetic_tm",PopulationStageResults=x,
       BackgroundDetectedResults=bg,InfoTriggeredDetectedResults=it,DetectedResults=bg,
       InformationStateBeforeSurveillanceResults=h,
       BackgroundDetectionProbabilityResults=pb,InfoTriggeredDetectionProbabilityResults=pi)
}
make_point <- function(with_info_event=TRUE,with_prob=TRUE) {
  s <- expand.grid(perm=1:4,timestep=1:2,KEEP.OUT.ATTRS=FALSE)
  s <- s[order(s$perm,s$timestep),]
  s$n_end <- c(0,0, 1,0, 1,1, 1,1)
  s$n_new_detections <- 0; s$n_new_background_detections <- 0
  s$n_new_background_detections[s$perm==3&s$timestep==2] <- 1
  if(with_info_event) {
    s$n_new_info_triggered_detections <- 0
    s$n_new_info_triggered_detections[s$perm==4&s$timestep==2] <- 1
  }
  s$n_info_before_surveillance <- 0
  s$n_info_before_surveillance[s$perm==4&s$timestep==2] <- 1
  ph <- data.frame()
  if(with_prob) {
    rows <- which(s$n_end>0)
    ph <- data.frame(perm=s$perm[rows],timestep=s$timestep[rows],id=seq_along(rows),x=0,y=0,
      have_info_before_surveillance=(s$perm[rows]==4&s$timestep[rows]==2),
      background_detection_prob=.2,info_triggered_detection_prob=.5,
      stringsAsFactors=FALSE)
  }
  list(ModelName="synthetic_point",Summary=s,PointHistory=ph)
}
make_ptm <- function() { z<-make_point(); z$ModelName<-"synthetic_ptm"; z$Summary$stage1_end<-z$Summary$n_end; z$Summary$stage2_end<-0; z }

obs_bg0_t2 <- data.frame(Timestep=2L,BackgroundDetections=0,InfoTriggeredDetections=NA_real_)
obs_both0_t2 <- data.frame(Timestep=2L,BackgroundDetections=0,InfoTriggeredDetections=0)

run_test("A01","contract","HaveInfo is state, not observational evidence",function(){
  x<-make_binary(info_event=FALSE); a<-.INApestPoAAdapterBinaryNode(x)
  assert(all(a$InfoTriggeredDetections==0)); assert(a$InfoTriggeredEventContract=="none"); "state/evidence separation retained"
})
run_test("A02","contract","Background and InfoTriggered event streams remain separate",function(){
  a<-.INApestPoAAdapterBinaryNode(make_binary()); assert(a$BackgroundDetections[3,2]==1);assert(a$InfoTriggeredDetections[4,2]==1); "event streams separated"
})
run_test("A03","contract","Legacy persistent DetectedResults is labelled legacy state",function(){
  a<-suppressWarnings(.INApestPoAAdapterBinaryNode(make_binary(explicit=FALSE,legacy_only=TRUE)))
  assert(a$BackgroundEventContract=="legacy_state"); "legacy state not promoted to explicit event"
})
run_test("A04","contract","Same-round background detection cannot activate InfoTriggered search",function(){
  x<-make_binary(); mod<-INApestPoAIndependentDetectionModel(.8,.5); a<-.INApestPoAAdapterBinaryNode(x);q<-mod(a,x,2L,"InfoTriggered")
  assert_close(q[3],1);assert_close(q[4],.5); "pre-round gate applied"
})

run_test("B01","likelihood","Background-only zero detection matches exact Bayes",function(){
  x<-make_binary();mod<-INApestPoAIndependentDetectionModel(.8,0)
  z<-INApestPoACore(x,.INApestPoAAdapterBinaryNode,PoAStartTimestep=2,ObservationHistory=obs_bg0_t2,EvidenceSources="Background",ObservationModel=mod,PoAMethod="likelihood")
  assert_close(tail(z$PoASummary$PosteriorPoA,1),.5/(.5+.5*.2)); "exact Bayes reproduced"
})
run_test("B02","likelihood","InfoTriggered zero detection is gated by prior information",function(){
  x<-make_binary();mod<-INApestPoAIndependentDetectionModel(0,.5)
  z<-INApestPoACore(x,.INApestPoAAdapterBinaryNode,PoAStartTimestep=2,ObservationHistory=data.frame(Timestep=2,InfoTriggeredDetections=0),EvidenceSources="InfoTriggered",ObservationModel=mod,PoAMethod="likelihood")
  assert_close(tail(z$PoASummary$PosteriorPoA,1),.5/(.5+.5*.75)); "targeted Bayes reproduced"
})
run_test("B03","likelihood","Combined q multiplies Background and InfoTriggered pathways",function(){
  x<-make_binary();mod<-INApestPoAIndependentDetectionModel(.8,.5)
  z<-INApestPoACore(x,.INApestPoAAdapterBinaryNode,PoAStartTimestep=2,ObservationHistory=obs_both0_t2,EvidenceSources=c("Background","InfoTriggered"),ObservationModel=mod,PoAMethod="likelihood")
  assert_close(tail(z$PoASummary$PosteriorPoA,1),.5/(.5+.5*.15)); "combined q reproduced"
})
run_test("B04","likelihood","Abundance likelihood uses per-individual detection",function(){
  x<-make_meta();a<-.INApestPoAAdapterPopulationNode(x);mod<-INApestPoAIndependentDetectionModel(.2,0);q<-mod(a,x,2,"Background");assert_close(q[3:4],c(.8^3,.8^3));"q=(1-p)^N"
})
run_test("B05","likelihood","MLU partition preserves likelihood",function(){
  x<-make_mlu();a<-.INApestPoAAdapterMultipleLandUse(x);mod<-INApestPoAIndependentDetectionModel(.2,0);q<-mod(a,x,2,"Background");assert_close(q[3:4],c(.8^2,.8^2));"MLU partition consistent"
})
run_test("B06","likelihood","Stage partition preserves likelihood",function(){
  x<-make_tm();a<-.INApestPoAAdapterTransitionMatrix(x);mod<-INApestPoAIndependentDetectionModel(.2,0);q<-mod(a,x,2,"Background");assert_close(q[3:4],c(.8^2,.8^2));"stage partition consistent"
})
run_test("B07","likelihood","Point independent model uses information strata",function(){
  x<-make_point();a<-.INApestPoAAdapterPoint(x);mod<-INApestPoAIndependentDetectionModel(.2,.5);q<-mod(a,x,2,"InfoTriggered");assert_close(q[3],1);assert_close(q[4],.5);"point pre-info strata respected"
})

run_test("R01","recorded","Likelihood auto-wires recorded Background probabilities",function(){
  x<-make_binary();z<-INApestPoACore(x,.INApestPoAAdapterBinaryNode,PoAStartTimestep=2,ObservationHistory=obs_bg0_t2,EvidenceSources="Background",PoAMethod="likelihood")
  assert(z$Settings$ObservationModelMode=="recorded_engine_detection_probabilities");assert_close(tail(z$PoASummary$PosteriorPoA,1),.5/(.5+.5*.2));"auto-wired recorded p"
})
run_test("R02","recorded","Recorded particle-specific detectability is retained",function(){
  x<-make_binary(); x$BackgroundDetectionProbabilityResults[1,2,3]<-.2; x$BackgroundDetectionProbabilityResults[1,2,4]<-.8
  a<-.INApestPoAAdapterBinaryNode(x);q<-INApestPoARecordedDetectionModel()(a,x,2,"Background")
  assert_close(q[3:4],c(.8,.2));"particle-level DetectionSD outcomes preserved"
})
run_test("R03","recorded","Recorded InfoTriggered p remains pre-information gated",function(){
  x<-make_binary();a<-.INApestPoAAdapterBinaryNode(x);q<-INApestPoARecordedDetectionModel()(a,x,2,"InfoTriggered")
  assert_close(q[3],1);assert_close(q[4],.5);"recorded target p gated correctly"
})
run_test("R04","recorded","Recorded MLU probabilities combine across land uses",function(){
  x<-make_mlu();a<-.INApestPoAAdapterMultipleLandUse(x);q<-INApestPoARecordedDetectionModel()(a,x,2,"Background")
  assert_close(q[3:4],c(.8^2,.8^2));"MLU recorded p exact"
})
run_test("R05","recorded","Recorded stage probabilities combine across stages",function(){
  x<-make_tm();a<-.INApestPoAAdapterTransitionMatrix(x);q<-INApestPoARecordedDetectionModel()(a,x,2,"Background")
  assert_close(q[3:4],c(.8^2,.8^2));"TM recorded p exact"
})
run_test("R06","recorded","Recorded point probabilities are evaluated per point",function(){
  x<-make_point();a<-.INApestPoAAdapterPoint(x);mod<-INApestPoARecordedDetectionModel();qbg<-mod(a,x,2,"Background");qi<-mod(a,x,2,"InfoTriggered")
  assert_close(qbg[3:4],c(.8,.8));assert_close(qi[3:4],c(1,.5));"point-specific p exact"
})
run_test("R07","recorded","Combined auto likelihood uses both recorded streams",function(){
  x<-make_binary();z<-INApestPoACore(x,.INApestPoAAdapterBinaryNode,PoAStartTimestep=2,ObservationHistory=obs_both0_t2,EvidenceSources=c("Background","InfoTriggered"),PoAMethod="likelihood")
  assert_close(tail(z$PoASummary$PosteriorPoA,1),.5/(.5+.5*.15));"combined recorded likelihood exact"
})

run_test("C01","sequential","Sequential likelihood rejects legacy persistent detection state",function(){
  x<-make_binary(explicit=FALSE,legacy_only=TRUE);mod<-INApestPoAIndependentDetectionModel(.5,0);obs<-data.frame(Timestep=1:2,BackgroundDetections=0)
  expect_error(suppressWarnings(INApestPoACore(x,.INApestPoAAdapterBinaryNode,ObservationHistory=obs,EvidenceSources="Background",ObservationModel=mod,PoAMethod="likelihood")),"Sequential PoA requires explicit BackgroundDetectedResults");"unsafe propagation blocked"
})
run_test("C02","sequential","Sequential targeted evidence requires explicit targeted events",function(){
  x<-make_binary(info_event=FALSE);mod<-INApestPoAIndependentDetectionModel(.5,.5);obs<-data.frame(Timestep=1:2,BackgroundDetections=0,InfoTriggeredDetections=0)
  expect_error(INApestPoACore(x,.INApestPoAAdapterBinaryNode,ObservationHistory=obs,EvidenceSources=c("Background","InfoTriggered"),ObservationModel=mod,PoAMethod="likelihood"),"Sequential PoA with InfoTriggered evidence requires explicit");"target history guard enforced"
})
run_test("C03","compatibility","Deprecated Surveillance alias reproduces Background result",function(){
  x<-make_binary();mod<-INApestPoAIndependentDetectionModel(.8,0)
  z1<-INApestPoACore(x,.INApestPoAAdapterBinaryNode,PoAStartTimestep=2,ObservationHistory=obs_bg0_t2,EvidenceSources="Background",ObservationModel=mod,PoAMethod="likelihood")
  z2<-expect_warning(INApestPoACore(x,.INApestPoAAdapterBinaryNode,PoAStartTimestep=2,ObservationHistory=data.frame(Timestep=2,SurveillanceDetections=0),EvidenceSources="Surveillance",ObservationModel=mod,PoAMethod="likelihood"),"deprecated")
  assert_close(tail(z1$PoASummary$PosteriorPoA,1),tail(z2$PoASummary$PosteriorPoA,1));"legacy alias retained"
})
run_test("C04","compatibility","HaveInfo alone cannot substitute for pre-surveillance information",function(){
  x<-make_binary(include_info_before=FALSE);x$HaveInfoResults<-array(1,dim=dim(x$InvasionResults));mod<-INApestPoAIndependentDetectionModel(0,.5)
  expect_error(INApestPoACore(x,.INApestPoAAdapterBinaryNode,PoAStartTimestep=2,ObservationHistory=data.frame(Timestep=2,InfoTriggeredDetections=0),EvidenceSources="InfoTriggered",ObservationModel=mod,PoAMethod="likelihood"),"InfoTriggered surveillance requires InformationStateBeforeSurveillanceResults");"HaveInfo not silently substituted"
})
run_test("C05","compatibility","Serial Meta PoA wrapper is available",function(){
  assert(exists("INApestMetaPoA",mode="function"));z<-INApestMetaPoA(ModelResults=make_meta(),PoAStartTimestep=2,ObservationHistory=obs_bg0_t2,EvidenceSources="Background",PoAMethod="likelihood");assert(inherits(z,"INApestPoA"));"serial Meta wrapper restored"
})

adapter_cases <- list(
  INApest=list(x=make_binary(),f=.INApestPoAAdapterBinaryNode),
  INApestMeta=list(x=make_meta(),f=.INApestPoAAdapterPopulationNode),
  INApestMetaMultipleLandUse=list(x=make_mlu(),f=.INApestPoAAdapterMultipleLandUse),
  INApestMetaTransitionMatrix=list(x=make_tm(),f=.INApestPoAAdapterTransitionMatrix),
  INApestMetaPoint=list(x=make_point(),f=.INApestPoAAdapterPoint),
  INApestPointTransitionMatrix=list(x=make_ptm(),f=.INApestPoAAdapterPointTransitionMatrix)
)
for(nm in names(adapter_cases)) run_test(paste0("D",sprintf("%02d",match(nm,names(adapter_cases)))),"coverage",paste0(nm," exposes the v2 contract"),function(){
  a<-adapter_cases[[nm]]$f(adapter_cases[[nm]]$x)
  assert(all(c("BackgroundDetections","InfoTriggeredDetections","ObservationAbundance","BackgroundEventContract","InfoTriggeredEventContract")%in%names(a)))
  assert(identical(dim(a$ObservationAbundance)[2:3],c(a$Ntimesteps,a$Nperm)))
  "adapter contract present"
})

res<-do.call(rbind,results);utils::write.csv(res,file.path(out_dir,"poa_observation_architecture_v2_results.csv"),row.names=FALSE)
summary<-c("INApest PoA observation architecture v2 validation",paste0("Date/time: ",format(Sys.time(),tz="",usetz=TRUE)),paste0("R: ",R.version.string),paste0("Platform: ",R.version$platform),paste0("PASS: ",sum(res$status=="PASS")," / ",nrow(res)),paste0("FAIL: ",sum(res$status=="FAIL")),"",paste(res$test_id,res$status,res$description,sep=" | "))
writeLines(summary,file.path(out_dir,"poa_observation_architecture_v2_summary.txt"));capture.output(sessionInfo(),file=file.path(out_dir,"sessionInfo.txt"));cat("\n",paste(summary,collapse="\n"),"\n",sep="")
if(any(res$status=="FAIL")) quit(status=1L)
