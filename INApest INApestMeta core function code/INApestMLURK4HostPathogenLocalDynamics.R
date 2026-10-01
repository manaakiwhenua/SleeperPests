###############################################################################
# INApestMLURK4HostPathogenLocalDynamics.R
# Coupled stochastic continuous-time host + pathogen LocalDynamics for MLU
# v0.6 -- 30 September 2026
###############################################################################
if (!exists("INApestStochasticCompartmentIntegrate", mode="function"))
  stop("Source INApestRK4.R and INApestStochasticCompartmentBridge.R first")

if(!exists("INApestContinuousPathogenRates",mode="function")) {
  INApestContinuousPathogenRates <- function(RecoveryRate=0,ProgressionRate=0,PathogenMortalityRate=0,ImmunityLossRate=0) {
    structure(list(RecoveryRate=RecoveryRate,ProgressionRate=ProgressionRate,PathogenMortalityRate=PathogenMortalityRate,ImmunityLossRate=ImmunityLossRate),class=c("INApestContinuousPathogenRates","list"))
  }
}

.inapest_mlu_expand <- function(x, template, label) {
  d <- dim(template); n <- length(template)
  if(is.function(x)) stop(label," resolver functions must be evaluated by HostFluxFunction")
  if(is.null(dim(x))) {
    z <- as.numeric(x)
    if(length(z)==1L) return(matrix(z,d[1],d[2]))
    if(length(z)==n) return(matrix(z,d[1],d[2]))
  } else if(identical(dim(x),d)) return(matrix(as.numeric(x),d[1],d[2]))
  stop(label," must be scalar, length node x land-use cells, or a node x land-use matrix")
}
.inapest_mlu_host_flux <- function(fun,t,Nmat,Parameters=NULL,runtime=list()) {
  ans <- .inapest_rk4_call_supported(fun,c(list(t=t,time=t,state=Nmat,State=Nmat,pars=Parameters,Parameters=Parameters),runtime))
  if(!is.list(ans)) stop("HostFluxFunction must return GainRate and Hazards")
  gain <- ans$GainRate; if(is.null(gain)) gain <- ans$Gains
  if(is.null(gain)) stop("HostFluxFunction must return GainRate")
  gain <- .inapest_mlu_expand(gain,Nmat,"Host GainRate")
  if(any(!is.finite(gain))||any(gain<0)) stop("Host GainRate must be finite and non-negative")
  H <- ans$Hazards
  if(is.null(H)||(is.list(H)&&length(H)==0L)) return(list(GainRate=c(gain),Hazards=matrix(numeric(),length(Nmat),0L)))
  if(!is.list(H)||is.null(names(H))||any(!nzchar(names(H)))||anyDuplicated(names(H))) stop("Host Hazards must be a uniquely named list")
  mat <- vapply(H,function(z)c(.inapest_mlu_expand(z,Nmat,"Host hazard")),numeric(length(Nmat)))
  if(is.null(dim(mat))) mat <- matrix(mat,ncol=1L,dimnames=list(NULL,names(H)))
  colnames(mat) <- names(H)
  if(any(!is.finite(mat))||any(mat<0)) stop("Host Hazards must be finite and non-negative")
  list(GainRate=c(gain),Hazards=mat)
}
.inapest_mlu_hp_rate_model <- function(t,State,HostFluxFunction,HostParameters,PathogenRates,Pathogen,PathogenEngine,PathogenContext,timestep,Ntimesteps,runtime=list()) {
  comps<-colnames(State); nu<-nrow(State); k<-ncol(State); nn<-PathogenContext$n_nodes; nl<-PathogenContext$n_landuses
  if(nu!=nn*nl||!all(c("S","I")%in%comps)) stop("Invalid MLU coupled pathogen state")
  N<-rowSums(State); Nmat<-matrix(N,nn,nl)
  hf<-.inapest_mlu_host_flux(HostFluxFunction,t,Nmat,HostParameters,runtime)
  gain<-matrix(0,nu,k,dimnames=dimnames(State)); gain[,"S"]<-hf$GainRate
  trans<-array(0,c(nu,k,k),dimnames=list(rownames(State),comps,comps))
  rr<-function(x,name){z<-as.numeric(PathogenEngine$Resolve(x,timestep,PathogenContext,name));if(length(z)!=nu||any(!is.finite(z))||any(z<0))stop(name," must resolve to finite non-negative rates per MLU cell");z}
  rec<-rr(PathogenRates$RecoveryRate,"RecoveryRate"); prog<-rr(PathogenRates$ProgressionRate,"ProgressionRate"); pmort<-rr(PathogenRates$PathogenMortalityRate,"PathogenMortalityRate"); wan<-rr(PathogenRates$ImmunityLossRate,"ImmunityLossRate")
  beta<-rr(Pathogen$Beta,"Beta"); ds<-rr(Pathogen$DensityScale,"DensityScale"); if(any(ds<=0))stop("DensityScale must be positive")
  C<-PathogenEngine$ContactMatrix(timestep,PathogenContext); I<-State[,"I"]; ip<-as.numeric(crossprod(I,C))
  if(Pathogen$Transmission=="frequency") {cp<-as.numeric(crossprod(N,C));lambda<-beta*ifelse(cp>0,ip/cp,0)} else lambda<-beta*ip/ds
  lambda<-pmax(0,lambda)
  if("E"%in%comps){trans[,"S","E"]<-lambda;trans[,"E","I"]<-prog}else{if(any(prog>0))stop("ProgressionRate must be zero for SIS/SIR");trans[,"S","I"]<-lambda}
  if(Pathogen$Model=="SIS"){trans[,"I","S"]<-rec;if(any(wan>0))stop("ImmunityLossRate must be zero for SIS")}else{trans[,"I","R"]<-rec;trans[,"R","S"]<-wan}
  hc<-colnames(hf$Hazards); causes<-c(if(length(hc))paste0("host:",hc) else character(),if(any(pmort>0))"pathogen" else character())
  exits<-array(0,c(nu,k,length(causes)),dimnames=list(rownames(State),comps,causes))
  if(length(hc)) for(j in seq_along(hc)) exits[,,paste0("host:",hc[j])]<-matrix(rep(hf$Hazards[,j],k),nu,k)
  if(any(pmort>0)) exits[,"I","pathogen"]<-pmort
  list(GainRate=gain,TransitionHazards=trans,ExitHazards=exits)
}

INApestMLURK4HostPathogenLocalDynamics <- function(HostFluxFunction,PathogenRates=INApestContinuousPathogenRates(),TimestepLength=1,RKMaxStep=0.025,HostParameters=NULL,StartTime=0,DispersalDynamics=NULL,RunWhenEmpty=FALSE,CapacityTolerance=1e-9) {
  if(!is.function(HostFluxFunction)||!inherits(PathogenRates,"INApestContinuousPathogenRates")) stop("Invalid host/pathogen RK specification")
  if(is.null(DispersalDynamics)){if(!exists("local.dynamicsLU",mode="function"))stop("Source current MLU engine first");DispersalDynamics<-get("local.dynamicsLU",mode="function")}
  host_fun<-HostFluxFunction; rates<-PathogenRates; pars<-HostParameters; dt<-as.numeric(TimestepLength)[1L]; hmax<-as.numeric(RKMaxStep)[1L]; t0<-as.numeric(StartTime)[1L]; disperse<-DispersalDynamics; cap_tol<-as.numeric(CapacityTolerance)[1L]
  if(!is.finite(dt)||dt<=0||!is.finite(hmax)||hmax<=0||!is.finite(t0)||!is.finite(cap_tol)||cap_tol<0)stop("Invalid RK settings")
  f<-function(sddprob,nodepropaguleproduction,nodeenvestabprob,n,lddprob,lddrate,k_is_0,nodeK,nodepropaguleestablishment,nodespreadreduction,nodefecundityreduction=0,managing,pathogen_state,pathogen,pathogen_engine,pathogen_context,timestep=NULL,Ntimesteps=NULL){
    if(is.null(pathogen)||!inherits(pathogen,"INApestPathogen")||pathogen$Model=="Binary")stop("Continuous MLU H+P requires non-binary Pathogen")
    nn<-nrow(n);nl<-ncol(n);if(!identical(dim(nodeK),dim(n)))stop("nodeK shape mismatch")
    if(!is.matrix(pathogen_state)||nrow(pathogen_state)!=length(n)||any(rowSums(pathogen_state)!=as.integer(c(n))))stop("pathogen_state must sum to flattened MLU host abundance")
    si<-if(is.null(timestep))1L else as.integer(timestep)[1L]; nts<-if(is.null(Ntimesteps))1L else as.integer(Ntimesteps)[1L]
    if(any(pathogen_engine$Resolve(pathogen$IntroductionProb,si,pathogen_context,"IntroductionProb")!=0))stop("Continuous MLU H+P v0.6 requires Pathogen$IntroductionProb = 0")
    runtime<-list(nodeK=nodeK,k_is_0=k_is_0,nodeenvestabprob=nodeenvestabprob,nodepropaguleproduction=nodepropaguleproduction,nodepropaguleestablishment=nodepropaguleestablishment,nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,managing=managing,sddprob=sddprob,lddprob=lddprob,lddrate=lddrate,timestep=si,Ntimesteps=nts)
    rf<-function(t,state,...) .inapest_mlu_hp_rate_model(t,state,host_fun,pars,rates,pathogen,pathogen_engine,pathogen_context,si,nts,runtime)
    z<-INApestStochasticCompartmentIntegrate(pathogen_state,rf,Duration=dt,MaxStep=hmax,StartTime=t0+(si-1L)*dt)
    p_after<-z$State;n_after<-matrix(as.integer(rowSums(p_after)),nn,nl)
    over<-which(n_after>nodeK+cap_tol*pmax(1,abs(nodeK)),arr.ind=TRUE);if(nrow(over))stop("Coupled RK host state exceeded MLU nodeK before dispersal")
    n_out<-disperse(sddprob=sddprob,nodepropaguleproduction=nodepropaguleproduction,nodeenvestabprob=nodeenvestabprob,n=n_after,lddprob=lddprob,lddrate=lddrate,k_is_0=k_is_0,nodeK=nodeK,nodepropaguleestablishment=nodepropaguleestablishment,nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,managing=managing)
    p_out<-pathogen_engine$Reconcile(p_after,n_out,pathogen_context);if(any(rowSums(p_out)!=as.integer(c(n_out))))stop("MLU H+P dispersal reconciliation failed")
    list(N=matrix(as.integer(n_out),nn,nl),PathogenState=p_out,CoupledRK=list(PreDispersalHost=n_after,PreDispersalPathogen=p_after,Flux=z,Timestep=si,Ntimesteps=nts))
  }
  attr(f,"INApestMLURK4HostPathogenLocalDynamics")<-list(Version="0.6",StateMode="integer-compartment",TimestepLength=dt,RKMaxStep=hmax,RunWhenEmpty=RunWhenEmpty)
  attr(f,"INApestCouplesPathogen")<-TRUE;attr(f,"INApestRunWhenEmpty")<-RunWhenEmpty
  class(f)<-c("INApestMLURK4HostPathogenLocalDynamics","function");f
}
