###############################################################################
### INApest PoA Layer 3 -- external published validation
###
### Staged ladder
###   3A  Published single-round Bayes / stopping calculations (Ramsey et al. 2023)
###   3B  Published repeated-survey Argentine-ant results (Ward et al. 2016)
###   3C  Published spatial detection kernels (Ward et al. 2016)
###   3D  Published broadscale / reintroduction framework (Anderson et al. 2017)
###   3E  Realistic multi-zone nutria application (Anderson et al. 2022)
###
### Status meanings
###   PASS      exact or published-rounding reproduction
###   SUPPORTED simplified public-data reconstruction behaves consistently with
###             a more detailed published model, but is not an exact replication
###   REVIEW    source-data boundary prevents full independent replication
###   FAIL      a reproducible target is missed
###############################################################################

args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) normalizePath(args[1L], mustWork=TRUE) else normalizePath(getwd(), mustWork=TRUE)
source_file <- file.path(root,"source","INApestPoA.R")
if(!file.exists(source_file)) source_file <- file.path(root,"..","..","definitive_source","INApestPoA.R")
if(!file.exists(source_file)) stop("PoA source not found for validation root: ",root)
source(source_file,local=.GlobalEnv)
out_dir <- file.path(root,"validation_output"); dir.create(out_dir,showWarnings=FALSE,recursive=TRUE)

results <- list(); metrics <- list(); kk <- 0L; mm <- 0L
assert <- function(x,msg="assertion failed") if(!isTRUE(x)) stop(msg,call.=FALSE)

record <- function(id,stage,paper,description,status,detail="") {
  kk <<- kk+1L
  results[[kk]] <<- data.frame(test_id=id,stage=stage,paper=paper,description=description,
    status=status,detail=detail,stringsAsFactors=FALSE)
  cat(sprintf("%-8s %-9s %-3s %s\n",id,status,stage,description))
}
metric <- function(id,quantity,observed,published,tolerance,evidence_level,note="") {
  mm <<- mm+1L
  metrics[[mm]] <<- data.frame(test_id=id,quantity=quantity,observed=as.numeric(observed),
    published_target=as.numeric(published),abs_difference=abs(as.numeric(observed)-as.numeric(published)),
    tolerance=as.numeric(tolerance),evidence_level=evidence_level,note=note,stringsAsFactors=FALSE)
}
run <- function(id,stage,paper,description,fun) {
  ans <- tryCatch(fun(),error=function(e) list(status="FAIL",detail=conditionMessage(e)))
  if(is.character(ans)) ans <- list(status="PASS",detail=ans)
  if(is.null(ans$status)) ans$status <- "PASS"
  if(is.null(ans$detail)) ans$detail <- ""
  record(id,stage,paper,description,ans$status,ans$detail)
}
close_or_fail <- function(id,quantity,observed,published,tol,evidence,note="") {
  metric(id,quantity,observed,published,tol,evidence,note)
  if(!is.finite(observed) || abs(observed-published)>tol)
    stop(sprintf("%s: published %.12g, observed %.12g, |difference| %.6g > tolerance %.6g",
      quantity,published,observed,abs(observed-published),tol),call.=FALSE)
}
close_supported <- function(id,quantity,observed,published,tol,evidence,note="") {
  metric(id,quantity,observed,published,tol,evidence,note)
  if(!is.finite(observed) || abs(observed-published)>tol)
    stop(sprintf("supporting approximation outside tolerance: published %.12g observed %.12g",published,observed),call.=FALSE)
  list(status="SUPPORTED",detail=sprintf("published %.5f; simplified public-data approximation %.5f",published,observed))
}

# Synthetic two-state carrier for exercising the actual INApest PoA core.
as_node_array <- function(x) {
  x <- as.matrix(x); np <- nrow(x); nt <- ncol(x)
  a <- array(0,dim=c(1L,nt,np)); a[1,,] <- t(x); a
}
make_binary <- function(sse,nt=1L) {
  present <- rbind(rep(0,nt),rep(1,nt))
  p <- rbind(rep(sse,nt),rep(sse,nt))
  zero <- matrix(0,nrow=2,ncol=nt)
  list(ModelName="layer3_external",InvasionResults=as_node_array(present),
       BackgroundDetectedResults=as_node_array(zero),InfoTriggeredDetectedResults=as_node_array(zero),
       DetectedResults=as_node_array(zero),InformationStateBeforeSurveillanceResults=as_node_array(zero),
       BackgroundDetectionProbabilityResults=as_node_array(p),InfoTriggeredDetectionProbabilityResults=as_node_array(zero))
}
core_one_round <- function(prior,sse) {
  x <- make_binary(sse)
  z <- INApestPoACore(x,.INApestPoAAdapterBinaryNode,PriorPoA=prior,
      ObservationHistory=data.frame(Timestep=1,BackgroundDetections=0),
      EvidenceSources="Background",PoAMethod="likelihood")
  rr <- z$PoASummary[z$PoASummary$Round == 1L,,drop=FALSE]
  if(nrow(rr) != 1L) stop("Layer-3 harness expected exactly one Round-1 posterior row")
  rr
}
bayes_zero <- function(prior,sse) prior/(1-sse*(1-prior))
required_sse <- function(prior,target) (target-prior)/(target*(1-prior))

###############################################################################
### Stage 3A -- Ramsey et al. 2023: simple published Bayes examples
###############################################################################
run("L3A-01","3A","Ramsey et al. 2023","Prior 0.5 + SSe 0.9 gives PoA > 0.90",function(){
  rr <- core_one_round(.5,.9); pub <- 10/11
  close_or_fail("L3A-01","posterior_poa",rr$PosteriorPoA,pub,1e-12,"exact published Bayes example")
  sprintf("PoA %.12f",rr$PosteriorPoA)
})
run("L3A-02","3A","Ramsey et al. 2023","Prior 0.8 + SSe 0.6 gives the same PoA > 0.90",function(){
  rr <- core_one_round(.8,.6); pub <- 10/11
  close_or_fail("L3A-02","posterior_poa",rr$PosteriorPoA,pub,1e-12,"exact published Bayes example")
  sprintf("PoA %.12f",rr$PosteriorPoA)
})
run("L3A-03","3A","Ramsey et al. 2023","Prior 0.9 requires SSe about 0.53 for target PoA 0.95",function(){
  s <- required_sse(.9,.95); rr <- core_one_round(.9,s)
  close_or_fail("L3A-03","required_sse",s,.53,.005,"published rounded stopping example")
  close_or_fail("L3A-03","posterior_at_exact_required_sse",rr$PosteriorPoA,.95,1e-12,"equation cross-check")
  sprintf("exact SSe %.9f; PoA %.9f",s,rr$PosteriorPoA)
})
run("L3A-04","3A","Ramsey et al. 2023","Prior 0.9 requires SSe about 0.91 for target PoA 0.99",function(){
  s <- required_sse(.9,.99); rr <- core_one_round(.9,s)
  close_or_fail("L3A-04","required_sse",s,.91,.005,"published rounded stopping example")
  close_or_fail("L3A-04","posterior_at_exact_required_sse",rr$PosteriorPoA,.99,1e-12,"equation cross-check")
  sprintf("exact SSe %.9f; PoA %.9f",s,rr$PosteriorPoA)
})

###############################################################################
### Stage 3B -- Ward et al. 2016: four sequential Argentine-ant surveys
###############################################################################
ward_sse <- c(.149,.713,.734,.736)
ward_post <- c(.312,.611,.855,.957)
ward_date <- c("March 2013","October 2013","November 2013","February 2014")
ward_prior_implied <- ward_post*(1-ward_sse)/(1-ward_post*ward_sse)
for(i in seq_along(ward_sse)) {
  id <- sprintf("L3B-%02d",i)
  run(id,"3B","Ward et al. 2016",paste0(ward_date[i]," published median SSe/PoE row"),local({ii<-i; iid<-id; function(){
    rr <- core_one_round(ward_prior_implied[ii],ward_sse[ii])
    close_or_fail(iid,"posterior_poe",rr$PosteriorPoA,ward_post[ii],5e-4,"published table/equation reproduction",
      sprintf("effective prior implied by rounded published row %.12f",ward_prior_implied[ii]))
    sprintf("SSe %.3f; effective prior %.6f; reproduced PoE %.3f",ward_sse[ii],ward_prior_implied[ii],rr$PosteriorPoA)
  }}))
}
run("L3B-05","3B","Ward et al. 2016","published PERT prior plus first-survey median SSe approximates first PoE independently",function(){
  # Standard beta-PERT parameterisation, lambda=4, min=0, mode=.25, max=.75.
  a<-0;m<-.25;b<-.75;lambda<-4
  alpha<-1+lambda*(m-a)/(b-a); beta<-1+lambda*(b-m)/(b-a)
  prior_med <- a+(b-a)*qbeta(.5,alpha,beta)
  rr <- core_one_round(prior_med,ward_sse[1])
  close_or_fail("L3B-05","first_survey_poe_from_pert_median",rr$PosteriorPoA,.312,.005,
    "independent published-input reconstruction","uses standard beta-PERT median and published median SSe")
  sprintf("PERT median prior %.6f -> PoE %.6f vs published 0.312",prior_med,rr$PosteriorPoA)
})
run("L3B-06","3B","Ward et al. 2016","between-survey effective priors show only small downward reintroduction discounts",function(){
  discounts <- 1-ward_prior_implied[2:4]/ward_post[1:3]
  assert(all(discounts>0 & discounts<.01),"effective inter-survey discounts are not all between 0 and 1%")
  list(status="PASS",detail=paste("effective discounts",paste(format(discounts,digits=6),collapse=", "),
    "consistent with the paper's very-low reintroduction assumption; not treated as an exact hazard reconstruction"))
})

###############################################################################
### Stage 3C -- Ward et al. 2016: spatially heterogeneous surveillance kernels
###############################################################################
halfnormal <- function(d,g0,sigma) g0*exp(-(d^2)/(2*sigma^2))
line_detection <- function(g0,sigma,spacing,n_each_side=100L) {
  d <- abs((-n_each_side:n_each_side)*spacing)
  1-prod(1-halfnormal(d,g0,sigma))
}
ward_kernel <- list(
  bait=list(g0=.548,sigma=1.331,spacing=2,published=.70,tol=.005),
  visual=list(g0=.733,sigma=.4,spacing=1,published=.75,tol=.005),
  dog=list(g0=.750,sigma=1.65,spacing=2,published=.90,tol=.01)
)
for(spec in list(c("L3C-01","bait"),c("L3C-02","visual"),c("L3C-03","dog"))) {
  id<-spec[1]; nm<-spec[2]; w<-ward_kernel[[nm]]
  run(id,"3C","Ward et al. 2016",paste0(nm," overlapping half-normal detection kernel reproduces published same-cell probability"),local({iid<-id;nn<-nm;ww<-w;function(){
    p <- line_detection(ww$g0,ww$sigma,ww$spacing)
    close_or_fail(iid,paste0(nn,"_same_line_detection"),p,ww$published,ww$tol,"spatial kernel reproduction")
    sprintf("calculated %.8f; published approximately %.2f",p,ww$published)
  }}))
}
run("L3C-04","3C","Ward et al. 2016","spatial kernel sensitivity can be passed unchanged to the PoA likelihood layer",function(){
  p <- line_detection(.733,.4,1)
  rr <- core_one_round(.5,p)
  close_or_fail("L3C-04","background_sse",rr$BackgroundSSe,p,1e-12,"INApest mapping of independently reconstructed spatial SSe")
  expected <- bayes_zero(.5,p)
  close_or_fail("L3C-04","posterior_poa",rr$PosteriorPoA,expected,1e-12,"Bayesian mapping")
  sprintf("spatial SSe %.8f -> PoA %.8f",p,rr$PosteriorPoA)
})

###############################################################################
### Stage 3D -- Anderson et al. 2017: broadscale and reintroduction framework
###############################################################################
run("L3D-01","3D","Anderson et al. 2017","published one-year Stage-I calibration reaches approximately 0.95 freedom",function(){
  pd<-.90; prp<-.98; pu<-1; prior<-.70
  se <- 1-(1-pd*prp)^pu
  rr <- core_one_round(prior,se)
  close_or_fail("L3D-01","stage1_posterior_freedom",rr$PosteriorPoA,.95,.005,"published worked calibration")
  metric("L3D-01","stage1_Se",se,.882,1e-12,"published Eq. 1 calculation")
  sprintf("Se %.3f; posterior %.8f",se,rr$PosteriorPoA)
})
run("L3D-02","3D","Anderson et al. 2017","ten independent 0.95-free management zones imply about 0.40 residual risk",function(){
  residual <- 1-.95^10
  close_or_fail("L3D-02","at_least_one_zone_not_free",residual,.40,.005,"published worked broadscale example")
  sprintf("exact residual risk %.9f; whole-area freedom %.9f",residual,1-residual)
})
run("L3D-03","3D","Anderson et al. 2017","published recursive Eq. 3/4 is reproduced over repeated negative surveys",function(){
  prior<-.25; se<-.35
  for(i in 1:4) {
    rr <- core_one_round(prior,se); expected<-bayes_zero(prior,se)
    close_or_fail("L3D-03",paste0("round_",i,"_posterior"),rr$PosteriorPoA,expected,1e-12,"published recurrence equation")
    prior <- rr$PosteriorPoA
  }
  sprintf("four-period recurrence ends at %.9f",prior)
})
run("L3D-04","3D","Anderson et al. 2017","reintroduction discount can be inserted between published Bayesian updates",function(){
  prior1<-.70; se1<-.60; intro<-.01; se2<-.60
  post1<-core_one_round(prior1,se1)$PosteriorPoA
  prior2<-post1*(1-intro)
  post2<-core_one_round(prior2,se2)$PosteriorPoA
  expected2<-bayes_zero(bayes_zero(prior1,se1)*(1-intro),se2)
  close_or_fail("L3D-04","posterior_after_discounted_second_round",post2,expected2,1e-12,"published reintroduction-discount structure")
  sprintf("round1 %.8f; 1%% discounted prior %.8f; round2 %.8f",post1,prior2,post2)
})

###############################################################################
### Stage 3E -- Anderson et al. 2022: realistic Delmarva nutria application
###############################################################################
starts <- c(Blackwater=1,Virginia=4,Maryland=9,Delaware=11)
pub2022 <- c(Blackwater=8,Virginia=11,Maryland=16,Delaware=18)
run("L3E-01","3E","Anderson et al. 2022","published occupied-cell growth rule reproduces all four 2022 zone values",function(){
  calc<-starts+7
  for(z in names(calc)) close_or_fail("L3E-01",paste0(z,"_occupied_cells_2022"),calc[z],pub2022[z],0,"published growth-rule reproduction")
  paste(names(calc),calc,collapse="; ")
})
public_p <- .03
run("L3E-02","3E","Anderson et al. 2022","Virginia public-only sensitivity is close to published full-model zone sensitivity",function(){
  p<-1-(1-public_p)^11
  close_supported("L3E-02","Virginia_2022_zone_sensitivity",p,.33,.06,"supporting consistency",
    "full published zone model includes uncertainty and all available surveillance")
})
run("L3E-03","3E","Anderson et al. 2022","Maryland public-only sensitivity is close to published full-model zone sensitivity",function(){
  p<-1-(1-public_p)^16
  close_supported("L3E-03","Maryland_2022_zone_sensitivity",p,.40,.03,"supporting consistency",
    "full published zone model includes uncertainty and all available surveillance")
})
run("L3E-04","3E","Anderson et al. 2022","Delaware public-only sensitivity is close to published full-model zone sensitivity",function(){
  p<-1-(1-public_p)^18
  close_supported("L3E-04","Delaware_2022_zone_sensitivity",p,.43,.02,"supporting consistency",
    "full published zone model includes uncertainty and all available surveillance")
})
run("L3E-05","3E","Anderson et al. 2022","Blackwater correctly requires active surveillance beyond the public-only component",function(){
  p_public<-1-(1-public_p)^8; published<-.90
  metric("L3E-05","Blackwater_public_only_vs_full_zone_sensitivity",p_public,published,Inf,"directional published comparison",
    "paper states active surveillance was concentrated in Blackwater")
  assert(p_public < published-.5,"public-only approximation did not remain well below published Blackwater total sensitivity")
  sprintf("public-only %.6f vs published full-model 0.90; large gap is expected from concentrated active surveillance",p_public)
})
run("L3E-06","3E","Anderson et al. 2022","published 0.01 prior to 0.75 PoA is exactly representable by an equivalent cumulative SSe",function(){
  prior<-.01; post<-.75
  q <- prior*(1-post)/(post*(1-prior)); sse<-1-q
  rr<-core_one_round(prior,sse)
  close_or_fail("L3E-06","equivalent_cumulative_sse",sse,.9966329966329966,1e-12,"published summary compatibility")
  close_or_fail("L3E-06","posterior_poa",rr$PosteriorPoA,.75,1e-12,"published summary compatibility")
  sprintf("equivalent cumulative SSe %.9f reproduces published mean PoA 0.75",sse)
})
run("L3E-07","3E","Anderson et al. 2022","full raw-data replication boundary is explicit",function(){
  list(status="REVIEW",detail=paste("The paper publishes code, parameters and summary outcomes but states that",
    "the underlying USDA surveillance data are not openly available. Therefore the exact 2015-2020",
    "spatial PoA trajectory cannot be independently rerun from public data; no synthetic substitute is treated as validation."))
})

res <- do.call(rbind,results)
met <- if(length(metrics)) do.call(rbind,metrics) else data.frame()
utils::write.csv(res,file.path(out_dir,"layer3_test_results.csv"),row.names=FALSE)
utils::write.csv(met,file.path(out_dir,"layer3_numeric_comparisons.csv"),row.names=FALSE)

n_pass<-sum(res$status=="PASS"); n_supported<-sum(res$status=="SUPPORTED"); n_review<-sum(res$status=="REVIEW"); n_fail<-sum(res$status=="FAIL")
summary <- c(
  "INApest PoA Layer 3 external published validation",
  paste0("Date/time: ",format(Sys.time(),tz="",usetz=TRUE)),
  paste0("R: ",R.version.string),
  paste0("Platform: ",R.version$platform),
  paste0("PASS: ",n_pass),
  paste0("SUPPORTED: ",n_supported),
  paste0("REVIEW / public-data boundary: ",n_review),
  paste0("FAIL: ",n_fail),
  "",
  paste(res$test_id,res$status,res$stage,res$description,res$detail,sep=" | ")
)
writeLines(summary,file.path(out_dir,"layer3_summary.txt"))
capture.output(sessionInfo(),file=file.path(out_dir,"sessionInfo.txt"))
cat("\n",paste(summary,collapse="\n"),"\n",sep="")
if(n_fail>0L) quit(status=1L)
