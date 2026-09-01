###############################################################################
### INApest PoA engine observation patch -- architecture-wide runtime tests
### Base R only; no external packages required except parallel (base/recommended).
###
### Runs the two-stream observation contract through all six biological model
### architectures and their parallel implementations.  These are small,
### deterministic tests: background p=0, targeted p=1, DetectionSD=0.
###############################################################################
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1L],mustWork=TRUE) else normalizePath(getwd(),mustWork=TRUE)
src <- file.path(root,"source"); if(!dir.exists(src)) src <- file.path(root,"..","..","definitive_source"); if(!dir.exists(src)) src <- root
out_dir <- file.path(root,"validation_output","engine_observation_patch");dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

results<-list();k<-0L
run_test<-function(id,arch,description,fun){k<<-k+1L;t0<-proc.time()[[3L]];status<-"PASS";detail<-"";tryCatch({detail<-as.character(fun())},error=function(e){status<<-"FAIL";detail<<-conditionMessage(e)});results[[k]]<<-data.frame(test_id=id,architecture=arch,description=description,status=status,detail=detail,elapsed_seconds=proc.time()[[3L]]-t0,stringsAsFactors=FALSE);cat(sprintf("%-9s %-5s %-38s %s\n",id,status,arch,description))}
assert<-function(x,msg="assertion failed")if(!isTRUE(x))stop(msg,call.=FALSE)
assert_close<-function(x,y,tol=1e-12,msg="values differ")if(length(x)!=length(y)||max(abs(as.numeric(x)-as.numeric(y)))>tol)stop(msg,call.=FALSE)
new_env<-function(files){e<-new.env(parent=.GlobalEnv);for(f in files){p<-file.path(src,f);if(!file.exists(p))stop("Missing source: ",p);sys.source(p,envir=e)};e}

check_node_target<-function(z){
  req<-c("BackgroundDetectedResults","InfoTriggeredDetectedResults","BackgroundDetectionProbabilityResults","InfoTriggeredDetectionProbabilityResults","InformationStateBeforeSurveillanceResults","HaveInfoResults")
  assert(all(req%in%names(z)),"missing observation output(s)")
  assert(sum(z$BackgroundDetectedResults)==0,"background p=0 produced detection")
  assert(sum(z$InfoTriggeredDetectedResults)==1,"target p=1 with pre-info did not detect exactly once")
  assert(all(z$BackgroundDetectionProbabilityResults==0),"recorded background p is wrong")
  assert(all(z$InfoTriggeredDetectionProbabilityResults==1),"recorded targeted p is wrong")
  assert(all(z$InformationStateBeforeSurveillanceResults==1),"pre-surveillance info not retained")
  TRUE
}

binary_args<-function(info=1,Ntimesteps=1,DetectionProb=0){
 list(ModelName="poa_bin_",Nperm=1,Ntimesteps=Ntimesteps,DetectionProb=DetectionProb,DetectionSD=0,
  ManageProb=0,ManageSD=0,EradicationProb=0,EradicationSD=0,SpreadReduction=0,SpreadReductionSD=0,
  InitialInvasion=1,InitialInfo=info,InfoRetentionProb=1,EnvEstabProb=1,Survival=1,
  SDDprob=matrix(0,1,1),LDDprob=0,SEAM=0,OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,
  SaveResults=FALSE,DoPlots=FALSE,Seed=11,InformationAcquisition="host",InfoTriggeredDetectionProb=1,InfoTriggeredDetectionSD=0)
}
meta_args<-function(info=1,parallel=FALSE){
 td<-file.path(tempdir(),paste0("poa_meta_",if(parallel)"p"else"s","_",sample.int(1e7,1)));dir.create(td)
 list(ModelName="m_",Nperm=1,Ntimesteps=1,DetectionProb=0,DetectionSD=0,ManageProb=0,ManageSD=0,
  MortalityProb=0,MortalitySD=0,FecundityReduction=0,SpreadReduction=0,SpreadReductionSD=0,
  InitialPopulation=1,InitBioP=NA,InvasionRisk=0,InitialInfo=info,InitInfoP=NA,ExternalInfoProb=0,
  InfoRetentionProb=1,InfoPersistenceSteps=NA,EnvEstabProb=1,Survival=1,K=100,PropaguleProduction=0,
  PropaguleEstablishment=0,IncursionStartPop=1,SDDprob=matrix(0,1,1),SEAM=0,LDDprob=NA,LDDrate=0,
  OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=paste0(td,"/"),DoPlots=FALSE,
  InformationAcquisition="host",InfoTriggeredDetectionProb=1,InfoTriggeredDetectionSD=0,ReturnResults=TRUE,
  Cores=if(parallel)1 else NULL,Seed=if(parallel)11 else NULL)
}
mlu_args<-function(info=1,parallel=FALSE){
 td<-file.path(tempdir(),paste0("poa_mlu_",if(parallel)"p"else"s","_",sample.int(1e7,1)));dir.create(td)
 list(ModelName="mlu_",Nperm=1,Ntimesteps=1,Nlanduses=2,DetectionProb=c(0,0),DetectionSD=c(0,0),
  ManageProb=c(0,0),ManageSD=c(0,0),MortalityProb=c(0,0),MortalitySD=c(0,0),FecundityReduction=c(0,0),
  SpreadReduction=c(0,0),SpreadReductionSD=c(0,0),InitialPopulation=matrix(c(1,0),1,2),InitBioP=NA,
  InvasionRisk=0,InitialInfo=info,InitInfoP=NA,ExternalInfoProb=0,InfoRetentionProb=1,InfoPersistenceSteps=NA,
  EnvEstabProb=1,Survival=1,K=matrix(100,1,2),PropaguleProduction=0,PropaguleEstablishment=0,
  IncursionStartPop=1,SDDprob=matrix(0,1,1),SEAM=0,LDDprob=NA,LDDrate=0,OngoingExternalInvasion=FALSE,
  OngoingExternalInfo=FALSE,OutputDir=paste0(td,"/"),DoPlots=FALSE,InformationAcquisition="host",
  InfoTriggeredDetectionProb=c(1,1),InfoTriggeredDetectionSD=c(0,0),ReturnResults=TRUE,
  Cores=if(parallel)1 else NULL,Seed=if(parallel)11 else NULL)
}
tm_args<-function(info=1,parallel=FALSE){
 td<-file.path(tempdir(),paste0("poa_tm_",if(parallel)"p"else"s","_",sample.int(1e7,1)));dir.create(td)
 list(ModelName="tm_",Nperm=1,Ntimesteps=1,Nstages=2,Weights=c(1,1),Transition=diag(2),
  DetectionProb=c(0,0),DetectionSD=matrix(0,1,2),ManageProb=0,ManageSD=0,MortalityProb=c(0,0),MortalitySD=matrix(0,1,2),
  FecundityReduction=0,SpreadReduction=0,SpreadReductionSD=0,InitialPopulation=matrix(c(0,1),1,2),InitBioP=NA,
  InvasionRisk=0,InitialInfo=info,InitInfoP=NA,ExternalInfoProb=0,InfoRetentionProb=1,InfoPersistenceSteps=NA,
  EnvEstabProb=1,K=100,SeedbankK=100,PropaguleEstablishment=0,IncursionStartPop=1,SDDprob=matrix(0,1,1),
  SEAM=0,LDDprob=NA,LDDrate=0,OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=paste0(td,"/"),
  DoPlots=FALSE,InfoTriggeredDetectionProb=c(1,1),InfoTriggeredDetectionSD=0,ReturnResults=TRUE,
  Cores=if(parallel)1 else NULL,Seed=if(parallel)11 else NULL)
}
point_args<-function(info=TRUE,Ntimesteps=1,DetectionProb=0){
 list(ModelName="pt_",Nperm=1,Ntimesteps=Ntimesteps,InitialPoints=data.frame(x=0,y=0),Survival=1,
  PropaguleProduction=0,PropaguleEstablishment=0,SDDkernel=NULL,DetectionProb=DetectionProb,DetectionSD=0,
  ManageProb=0,ManageSD=0,MortalityProb=0,MortalitySD=0,InitialInfo=info,InfoRetentionProb=1,
  SaveResults=FALSE,DoProgress=FALSE,Seed=11,InfoTriggeredDetectionProb=1,InfoTriggeredDetectionSD=0)
}
ptm_args<-function(info=TRUE){
 list(ModelName="ptm_",Nperm=1,Ntimesteps=1,Nstages=2,Weights=c(1,1),Transition=diag(2),InitialPoints=data.frame(x=0,y=0,stage=2L),
  SDDkernel=NULL,PropaguleEstablishment=0,DetectionProb=0,DetectionSD=0,ManageProb=0,ManageSD=0,MortalityProb=0,MortalitySD=0,
  InitialInfo=info,InfoRetentionProb=1,SaveResults=FALSE,DoProgress=FALSE,Seed=11,
  InfoTriggeredDetectionProb=1,InfoTriggeredDetectionSD=0)
}

### Binary serial + parallel
e<-new_env(c("INApest.R","INApestParallel.R"))
run_test("E01","INApest","serial target stream + probability bookkeeping",function(){z<-do.call(e$INApest,binary_args(1));check_node_target(z);"serial contract passed"})
run_test("E02","INApest","target stream is inactive without pre-existing information",function(){z<-do.call(e$INApest,binary_args(0));assert(sum(z$InfoTriggeredDetectedResults)==0);assert(sum(z$InformationStateBeforeSurveillanceResults)==0);"pre-info gate passed"})
run_test("E03","INApest","same-round background detection does not trigger targeted search",function(){a<-binary_args(0,2,matrix(c(0,1),nrow=1));z<-do.call(e$INApest,a);assert(z$BackgroundDetectedResults[1,2,1]==1);assert(z$InfoTriggeredDetectedResults[1,2,1]==0);assert(z$InformationStateBeforeSurveillanceResults[1,2,1]==0);"same-round rule passed"})
run_test("E04","INApestParallel","parallel wrapper combines new outputs",function(){a<-binary_args(1);a$Cores<-1;a$Backend<-"psock";z<-do.call(e$INApestParallel,a);check_node_target(z);"parallel contract passed"})

### Meta serial + parallel
e<-new_env(c("INApestMeta.r","INApestMetaParallel.r"))
run_test("E05","INApestMeta","serial target stream + probability bookkeeping",function(){a<-meta_args(1,FALSE);a$Cores<-NULL;a$Seed<-NULL;z<-do.call(e$INApestMeta,a);check_node_target(z);"serial contract passed"})
run_test("E06","INApestMeta","target stream inactive without pre-info",function(){a<-meta_args(0,FALSE);a$Cores<-NULL;a$Seed<-NULL;z<-do.call(e$INApestMeta,a);assert(sum(z$InfoTriggeredDetectedResults)==0);"pre-info gate passed"})
run_test("E07","INApestMetaParallel","parallel target stream + probability bookkeeping",function(){z<-do.call(e$INApestMetaParallel,meta_args(1,TRUE));check_node_target(z);"parallel contract passed"})
run_test("E07b","INApestMeta","targeted-off probability bookkeeping is explicit zero",function(){a<-meta_args(1,FALSE);a$Cores<-NULL;a$Seed<-NULL;a$InfoTriggeredDetectionProb<-0;a$InfoTriggeredDetectionSD<-0;z<-do.call(e$INApestMeta,a);assert(all(z$InfoTriggeredDetectedResults==0));assert(all(z$InfoTriggeredDetectionProbabilityResults==0));"targeted-off zero contract passed"})

### MLU serial + parallel
e<-new_env(c("INApestMetaMultipleLandUse.r","INApestMetaParallelMultipleLandUse.r"))
run_test("E08","INApestMetaMultipleLandUse","serial node-event / land-use probability contract",function(){a<-mlu_args(1,FALSE);a$Cores<-NULL;a$Seed<-NULL;z<-do.call(e$INApestMetaMultipleLandUse,a);check_node_target(z);assert(identical(dim(z$BackgroundDetectionProbabilityResults),c(1L,2L,1L,1L)));"serial MLU contract passed"})
run_test("E09","INApestMetaParallelMultipleLandUse","parallel node-event / land-use probability contract",function(){z<-do.call(e$INApestMetaParallelMultipleLandUse,mlu_args(1,TRUE));check_node_target(z);"parallel MLU contract passed"})

### Transition matrix serial + parallel
e<-new_env(c("INApestMetaTransitionMatrix.r","INApestMetaTransitionMatrixParallel.r"))
run_test("E10","INApestMetaTransitionMatrix","serial node-event / stage-probability contract",function(){a<-tm_args(1,FALSE);a$Cores<-NULL;a$Seed<-NULL;z<-do.call(e$INApestMetaTransitionMatrix,a);check_node_target(z);assert(identical(dim(z$BackgroundDetectionProbabilityResults),c(1L,2L,1L,1L)));"serial TM contract passed"})
run_test("E11","INApestMetaTransitionMatrixParallel","parallel node-event / stage-probability contract",function(){z<-do.call(e$INApestMetaTransitionMatrixParallel,tm_args(1,TRUE));check_node_target(z);"parallel TM contract passed"})

### MetaPoint serial + parallel
e<-new_env(c("INApestMetaPoint.R","INApestMetaPointParallel.R")); fixed<-e$INApestPointKernelFixed(0)
run_test("E12","INApestMetaPoint","serial point observation fields",function(){a<-point_args(TRUE);a$SDDkernel<-fixed;z<-do.call(e$INApestMetaPoint,a);r<-z$Summary[1,];assert(r$n_new_background_detections==0);assert(r$n_new_info_triggered_detections==1);assert(r$n_info_before_surveillance==1);assert(all(z$PointHistory$background_detection_prob==0));assert(all(z$PointHistory$info_triggered_detection_prob==1));"serial point contract passed"})
run_test("E13","INApestMetaPoint","target stream inactive without pre-info",function(){a<-point_args(FALSE);a$SDDkernel<-fixed;z<-do.call(e$INApestMetaPoint,a);assert(z$Summary$n_new_info_triggered_detections[1]==0);"pre-info gate passed"})
run_test("E14","INApestMetaPointParallel","parallel point observation fields",function(){a<-point_args(TRUE);a$SDDkernel<-fixed;a$Cores<-1;a$Backend<-"psock";z<-do.call(e$INApestMetaPointParallel,a);assert(z$Summary$n_new_info_triggered_detections[1]==1);assert(all(c("background_detection_prob","info_triggered_detection_prob","have_info_before_surveillance")%in%names(z$PointHistory)));"parallel point contract passed"})

### Point transition serial + parallel
e<-new_env(c("INApestPointTransitionMatrix.R","INApestPointTransitionMatrixParallel.R"));fixed<-e$INApestPointKernelFixed(0)
run_test("E15","INApestPointTransitionMatrix","serial point-transition observation fields",function(){a<-ptm_args(TRUE);a$SDDkernel<-fixed;z<-do.call(e$INApestPointTransitionMatrix,a);assert(z$Summary$n_new_background_detections[1]==0);assert(z$Summary$n_new_info_triggered_detections[1]==1);assert(z$Summary$n_info_before_surveillance[1]==1);assert(all(z$PointHistory$info_triggered_detection_prob==1));"serial point-TM contract passed"})
run_test("E16","INApestPointTransitionMatrixParallel","parallel point-transition observation fields",function(){a<-ptm_args(TRUE);a$SDDkernel<-fixed;a$Cores<-1;a$Backend<-"psock";z<-do.call(e$INApestPointTransitionMatrixParallel,a);assert(z$Summary$n_new_info_triggered_detections[1]==1);assert(all(c("background_detection_prob","info_triggered_detection_prob","have_info_before_surveillance")%in%names(z$PointHistory)));"parallel point-TM contract passed"})

res<-do.call(rbind,results);utils::write.csv(res,file.path(out_dir,"poa_engine_observation_patch_results.csv"),row.names=FALSE)
summary<-c("INApest PoA engine observation patch runtime validation",paste0("Date/time: ",format(Sys.time(),tz="",usetz=TRUE)),paste0("R: ",R.version.string),paste0("Platform: ",R.version$platform),paste0("PASS: ",sum(res$status=="PASS")," / ",nrow(res)),paste0("FAIL: ",sum(res$status=="FAIL")),"",paste(res$test_id,res$status,res$architecture,res$description,sep=" | "))
writeLines(summary,file.path(out_dir,"poa_engine_observation_patch_summary.txt"));capture.output(sessionInfo(),file=file.path(out_dir,"sessionInfo.txt"));cat("\n",paste(summary,collapse="\n"),"\n",sep="")
if(any(res$status=="FAIL")) quit(status=1L)
