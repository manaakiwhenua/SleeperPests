###############################################################################
# INApestVNRK4LocalDynamics.R -- host-only continuous local biology + Vertebrate Node boundary
# v0.8 -- 1 October 2026
###############################################################################
INApestVNRK4LocalDynamics <- function(HostFluxFunction,TimestepLength=1,RKMaxStep=0.025,Parameters=NULL,StartTime=0,TransitionDynamics=NULL,RunWhenEmpty=FALSE,CapacityTolerance=1e-9){
  if(!is.function(HostFluxFunction))stop("HostFluxFunction must be a function")
  if(is.null(TransitionDynamics)){if(!exists(".iv_local_dynamics_transition_matrix",mode="function"))stop("Source current Vertebrate Node engine first");TransitionDynamics<-get(".iv_local_dynamics_transition_matrix",mode="function")}
  stochastic<-INApestStochasticLocalDynamics(HostFluxFunction,TimestepLength,RKMaxStep,Parameters,StartTime);disc<-TransitionDynamics;tol<-as.numeric(CapacityTolerance)[1L]
  f<-function(nodetransition,weights,sddprob,nodeenvestabprob,n0,lddprob=NA,lddrate=0,nodeK,node.seedbankK,nodepropaguleestablishment,nodespreadreduction,managing,MaxInteger=.Machine$integer.max,nodefecundityreduction=0,nodecontrolfecundityreduction=0,nodebirthmean=NULL,nodebirthmothers=NULL,transition_sddprob=NULL,transition_lddprob=NULL,transition_lddrate=0,BlockedTransitionMortality=0,DispersalDensityFactor=0,ApplyFootprintToTransitions=FALSE,timestep=NULL,Ntimesteps=NULL){
    N0<-as.matrix(n0)
    N1<-stochastic(n0=N0,timestep=timestep,Ntimesteps=Ntimesteps,nodeK=nodeK,node.seedbankK=node.seedbankK,weights=weights,nodetransition=nodetransition,managing=managing)
    if(!is.matrix(N1)||!identical(dim(N1),dim(N0)))stop("Vertebrate RK bridge did not preserve node x stage shape")
    .inapest_tm_capacity_check(N1,nodeK,node.seedbankK,weights,tol,"Vertebrate host-only RK state")
    out<-disc(nodetransition=nodetransition,weights=weights,sddprob=sddprob,nodeenvestabprob=nodeenvestabprob,n0=N1,lddprob=lddprob,lddrate=lddrate,nodeK=nodeK,node.seedbankK=node.seedbankK,nodepropaguleestablishment=nodepropaguleestablishment,nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,nodecontrolfecundityreduction=nodecontrolfecundityreduction,nodebirthmean=nodebirthmean,nodebirthmothers=nodebirthmothers,managing=managing,MaxInteger=MaxInteger,transition_sddprob=transition_sddprob,transition_lddprob=transition_lddprob,transition_lddrate=transition_lddrate,ApplyFootprintToTransitions=ApplyFootprintToTransitions,BlockedTransitionMortality=BlockedTransitionMortality,DispersalDensityFactor=DispersalDensityFactor)
    attr(out,"INApestVNRK4")<-list(PreTransitionHost=N1,Timestep=timestep,Ntimesteps=Ntimesteps,Flux=attr(N1,"INApestStochasticFlux"));out
  }
  attr(f,"INApestVNRK4LocalDynamics")<-list(Version="0.8",StateMode="integer",Boundary="continuous biology then discrete Vertebrate Node transition")
  attr(f,"INApestRunWhenEmpty")<-isTRUE(RunWhenEmpty);class(f)<-c("INApestVNRK4LocalDynamics","function");f
}
