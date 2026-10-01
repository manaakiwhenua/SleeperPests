###############################################################################
# INApestTMRK4LocalDynamics.R -- host-only continuous local biology + discrete TM boundary
# v0.7 -- 30 September 2026
###############################################################################
INApestTMRK4LocalDynamics <- function(HostFluxFunction,TimestepLength=1,RKMaxStep=0.025,Parameters=NULL,StartTime=0,TransitionDynamics=NULL,RunWhenEmpty=FALSE,CapacityTolerance=1e-9){
  if(!is.function(HostFluxFunction))stop("HostFluxFunction must be a function")
  if(is.null(TransitionDynamics)){if(!exists("local.dynamics.transition.matrix",mode="function"))stop("Source current TM engine first");TransitionDynamics<-get("local.dynamics.transition.matrix",mode="function")}
  stochastic<-INApestStochasticLocalDynamics(HostFluxFunction,TimestepLength,RKMaxStep,Parameters,StartTime);disc<-TransitionDynamics;tol<-as.numeric(CapacityTolerance)[1L]
  f<-function(nodetransition,weights,sddprob,nodeenvestabprob,n0,lddprob=NA,lddrate=0,nodeK,node.seedbankK,nodepropaguleestablishment,nodespreadreduction,managing,MaxInteger=.Machine$integer.max,nodefecundityreduction=0,transition_sddprob=NULL,transition_lddprob=NULL,transition_lddrate=0,BlockedTransitionMortality=0,DispersalDensityFactor=0,timestep=NULL,Ntimesteps=NULL){
    N0<-as.matrix(n0);N1<-stochastic(n0=N0,timestep=timestep,Ntimesteps=Ntimesteps,nodeK=nodeK,node.seedbankK=node.seedbankK,weights=weights,nodetransition=nodetransition,managing=managing)
    if(!is.matrix(N1)||!identical(dim(N1),dim(N0)))stop("TM RK bridge did not preserve node x stage shape");.inapest_tm_capacity_check(N1,nodeK,node.seedbankK,weights,tol,"Host-only RK state")
    out<-disc(nodetransition=nodetransition,weights=weights,sddprob=sddprob,nodeenvestabprob=nodeenvestabprob,n0=N1,lddprob=lddprob,lddrate=lddrate,nodeK=nodeK,node.seedbankK=node.seedbankK,nodepropaguleestablishment=nodepropaguleestablishment,nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,managing=managing,MaxInteger=MaxInteger,transition_sddprob=transition_sddprob,transition_lddprob=transition_lddprob,transition_lddrate=transition_lddrate,BlockedTransitionMortality=BlockedTransitionMortality,DispersalDensityFactor=DispersalDensityFactor)
    attr(out,"INApestTMRK4")<-list(PreTransitionHost=N1,Timestep=timestep,Ntimesteps=Ntimesteps,Flux=attr(N1,"INApestStochasticFlux"));out
  }
  attr(f,"INApestTMRK4LocalDynamics")<-list(Version="0.7",StateMode="integer",Boundary="continuous biology then discrete TM transition")
  attr(f,"INApestRunWhenEmpty")<-isTRUE(RunWhenEmpty);class(f)<-c("INApestTMRK4LocalDynamics","function");f
}
