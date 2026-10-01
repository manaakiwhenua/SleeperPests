###############################################################################
# INApestTMRK4HostPathogenLocalDynamics.R -- H+P continuous biology + discrete TM transport
# v0.7 -- 30 September 2026
###############################################################################
INApestTMRK4HostPathogenLocalDynamics <- function(HostFluxFunction,PathogenRates=INApestContinuousPathogenRates(),StageMixing=NULL,TimestepLength=1,RKMaxStep=0.025,HostParameters=NULL,StartTime=0,PathogenTransitionDynamics=NULL,RunWhenEmpty=FALSE,CapacityTolerance=1e-9){
  if(!is.function(HostFluxFunction)||!inherits(PathogenRates,"INApestContinuousPathogenRates"))stop("Invalid TM H+P specification")
  if(is.null(PathogenTransitionDynamics))PathogenTransitionDynamics<-get("local.dynamics.transition.matrix.pathogen",mode="function")
  hf<-HostFluxFunction;pr<-PathogenRates;pars<-HostParameters;sm<-StageMixing;dt<-as.numeric(TimestepLength)[1L];hmax<-as.numeric(RKMaxStep)[1L];t0<-as.numeric(StartTime)[1L];disc<-PathogenTransitionDynamics;tol<-as.numeric(CapacityTolerance)[1L]
  f<-function(nodetransition,weights,sddprob,nodeenvestabprob,n0,lddprob=NA,lddrate=0,nodeK,node.seedbankK,nodepropaguleestablishment,nodespreadreduction,managing,MaxInteger=.Machine$integer.max,nodefecundityreduction=0,transition_sddprob=NULL,transition_lddprob=NULL,transition_lddrate=0,BlockedTransitionMortality=0,DispersalDensityFactor=0,pathogen_state,Pathogen,timestep,Ntimesteps){
    N0<-as.matrix(n0);nn<-nrow(N0);D<-ncol(N0);si<-as.integer(timestep)[1L];nts<-as.integer(Ntimesteps)[1L]
    if(!inherits(Pathogen,"INApestPathogen")||Pathogen$Model=="Binary")stop("Continuous TM H+P requires non-binary INApestPathogen")
    if(any(.iptm_resolve(Pathogen$IntroductionProb,si,nn,D,nts,"IntroductionProb",TRUE)!=0))stop("Continuous TM H+P v0.7 requires Pathogen$IntroductionProb = 0")
    H0<-.inapest_tm_pathogen_flat(pathogen_state);if(any(rowSums(H0)!=.inapest_tm_flat(N0)))stop("TM pathogen_state must sum stage-wise to host abundance")
    runtime<-list(nodeK=nodeK,node.seedbankK=node.seedbankK,weights=weights,nodetransition=nodetransition,managing=managing,timestep=si,Ntimesteps=nts)
    rf<-function(t,state,...).inapest_tm_rate_model(t,state,hf,pars,pr,Pathogen,si,nts,nn,D,sm,NULL,runtime)
    z<-INApestStochasticCompartmentIntegrate(H0,rf,Duration=dt,MaxStep=hmax,StartTime=t0+(si-1L)*dt)
    H1<-z$State;N1<-.inapest_tm_unflat(rowSums(H1),nn,D);.inapest_tm_capacity_check(N1,nodeK,node.seedbankK,weights,tol,"H+P RK state")
    st1<-.inapest_tm_pathogen_array(H1,nn,D);Ptransport<-.inapest_tm_transport_pathogen(Pathogen)
    out<-disc(nodetransition=nodetransition,weights=weights,sddprob=sddprob,nodeenvestabprob=nodeenvestabprob,n0=N1,lddprob=lddprob,lddrate=lddrate,nodeK=nodeK,node.seedbankK=node.seedbankK,nodepropaguleestablishment=nodepropaguleestablishment,nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,managing=managing,MaxInteger=MaxInteger,transition_sddprob=transition_sddprob,transition_lddprob=transition_lddprob,transition_lddrate=transition_lddrate,BlockedTransitionMortality=BlockedTransitionMortality,DispersalDensityFactor=DispersalDensityFactor,pathogen_state=st1,Pathogen=Ptransport,timestep=si,Ntimesteps=nts,StageMixing=sm)
    pd<-if("pathogen"%in%colnames(z$Exits)).inapest_tm_unflat(z$Exits[,"pathogen"],nn,D)else matrix(0L,nn,D)
    out$PathogenDeaths<-matrix(as.integer(pd),nn,D);out$NewInfections<-matrix(NA_integer_,nn,D)
    out$CoupledRK<-list(PreTransitionHost=N1,PreTransitionPathogen=st1,Flux=z,NewInfectionsDiagnostic="not event-counted in continuous bridge")
    out
  }
  attr(f,"INApestTMRK4HostPathogenLocalDynamics")<-list(Version="0.7",StateMode="integer-compartment",Boundary="continuous H+P then discrete TM pathogen transport")
  attr(f,"INApestCouplesPathogen")<-TRUE;attr(f,"INApestRunWhenEmpty")<-isTRUE(RunWhenEmpty);class(f)<-c("INApestTMRK4HostPathogenLocalDynamics","function");f
}
