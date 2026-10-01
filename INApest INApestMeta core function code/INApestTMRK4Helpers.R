###############################################################################
# INApestTMRK4Helpers.R
# Shared helpers for stage-structured stochastic RK local biology
# v0.7 -- 30 September 2026
###############################################################################
if(!exists("INApestRK4Step",mode="function") || !exists("INApestStochasticCompartmentStep",mode="function"))
  stop("Source INApestRK4.R and INApestStochasticCompartmentBridge.R first")
if(!exists(".iptm_resolve",mode="function") || !exists("local.dynamics.transition.matrix.pathogen",mode="function"))
  stop("Source INApestPathogenTransitionMatrix.R first")

if(!exists("INApestContinuousPathogenRates",mode="function")) {
  INApestContinuousPathogenRates <- function(RecoveryRate=0,ProgressionRate=0,PathogenMortalityRate=0,ImmunityLossRate=0) {
    structure(list(RecoveryRate=RecoveryRate,ProgressionRate=ProgressionRate,
      PathogenMortalityRate=PathogenMortalityRate,ImmunityLossRate=ImmunityLossRate),
      class=c("INApestContinuousPathogenRates","list"))
  }
}
if(!exists("INApestContinuousBiocontrolAgentRates",mode="function")) {
  INApestContinuousBiocontrolAgentRates <- function(TransitionRates=NULL,ExitRates=0)
    structure(list(TransitionRates=TransitionRates,ExitRates=ExitRates),class=c("INApestContinuousBiocontrolAgentRates","list"))
  INApestContinuousBiocontrolRates <- function(AgentRates) {
    if(!is.list(AgentRates)||!length(AgentRates)||is.null(names(AgentRates))||any(!nzchar(names(AgentRates)))||anyDuplicated(names(AgentRates)))
      stop("AgentRates must be a uniquely named non-empty list")
    if(!all(vapply(AgentRates,inherits,logical(1),"INApestContinuousBiocontrolAgentRates")))
      stop("Every AgentRates entry must come from INApestContinuousBiocontrolAgentRates()")
    structure(AgentRates,class=c("INApestContinuousBiocontrolRates","list"))
  }
}

.inapest_tm_flat <- function(x) as.numeric(t(as.matrix(x)))
.inapest_tm_unflat <- function(x,n_nodes,n_stages) t(matrix(as.numeric(x),nrow=n_stages,ncol=n_nodes))

.inapest_tm_expand <- function(x,template,label) {
  nr<-nrow(template); ns<-ncol(template); n<-nr*ns
  d<-dim(x)
  if(is.null(d)) {
    z<-as.numeric(x)
    if(length(z)==1L) return(matrix(z,nr,ns))
    if(length(z)==ns && ns!=nr) return(matrix(rep(z,each=nr),nr,ns))
    if(length(z)==nr && nr!=ns) return(matrix(rep(z,ns),nr,ns))
    if(length(z)==n) return(matrix(z,nr,ns))
  } else if(length(d)==2L && identical(d,c(nr,ns))) return(as.matrix(x))
  stop(label," must be scalar, unambiguous node/stage vector, length nodes*stages, or nodes x stages matrix")
}

.inapest_tm_host_flux <- function(fun,t,Nmat,Parameters=NULL,runtime=list()) {
  ans<-.inapest_rk4_call_supported(fun,c(list(t=t,time=t,state=Nmat,State=Nmat,pars=Parameters,Parameters=Parameters),runtime))
  if(!is.list(ans)) stop("HostFluxFunction must return GainRate and Hazards")
  gain<-ans$GainRate;if(is.null(gain))gain<-ans$Gains;if(is.null(gain))stop("HostFluxFunction must return GainRate")
  gain<-.inapest_tm_expand(gain,Nmat,"Host GainRate");if(any(!is.finite(gain))||any(gain<0))stop("Host GainRate must be finite and non-negative")
  H<-ans$Hazards
  if(is.null(H)||(is.list(H)&&length(H)==0L))return(list(GainRate=.inapest_tm_flat(gain),Hazards=matrix(numeric(),length(Nmat),0L)))
  if(!is.list(H)||is.null(names(H))||any(!nzchar(names(H)))||anyDuplicated(names(H)))stop("Host Hazards must be a uniquely named list")
  mat<-vapply(H,function(z).inapest_tm_flat(.inapest_tm_expand(z,Nmat,"Host hazard")),numeric(length(Nmat)))
  if(is.null(dim(mat)))mat<-matrix(mat,ncol=1L,dimnames=list(NULL,names(H)));colnames(mat)<-names(H)
  if(any(!is.finite(mat))||any(mat<0))stop("Host Hazards must be finite and non-negative")
  list(GainRate=.inapest_tm_flat(gain),Hazards=mat)
}

.inapest_tm_stage_mixing <- function(StageMixing,D) {
  if(is.null(StageMixing))StageMixing<-matrix(1,D,D)
  if(!is.matrix(StageMixing)||!identical(dim(StageMixing),c(D,D))||any(!is.finite(StageMixing))||any(StageMixing<0))
    stop("StageMixing must be a finite non-negative Nstages x Nstages matrix")
  StageMixing
}
.inapest_tm_resolve_rate <- function(x,timestep,n_nodes,D,Ntimesteps,name) {
  z<-.iptm_resolve(x,timestep,n_nodes,D,Ntimesteps,name,FALSE)
  z<-.inapest_tm_flat(z)
  if(any(!is.finite(z))||any(z<0))stop(name," must resolve to finite non-negative rates by node x stage")
  z
}
.inapest_tm_pathogen_contact <- function(Pathogen,timestep,n_nodes,D,Ntimesteps,StageMixing) {
  kronecker(.iptm_node_contact(Pathogen,timestep,n_nodes,Ntimesteps),.inapest_tm_stage_mixing(StageMixing,D))
}
.inapest_tm_pathogen_flat <- function(state3) {
  d<-dim(state3);if(length(d)!=3L)stop("pathogen_state must be node x stage x pathogen-state")
  out<-matrix(0L,d[1]*d[2],d[3],dimnames=list(NULL,dimnames(state3)[[3]]))
  for(i in seq_len(d[1]))for(s in seq_len(d[2]))out[(i-1L)*d[2]+s,]<-state3[i,s,]
  storage.mode(out)<-"integer";out
}
.inapest_tm_pathogen_array <- function(flat,n_nodes,D) {
  k<-ncol(flat);out<-array(0L,c(n_nodes,D,k),dimnames=list(NULL,seq_len(D),colnames(flat)))
  for(i in seq_len(n_nodes))for(s in seq_len(D))out[i,s,]<-flat[(i-1L)*D+s,]
  out
}
.inapest_tm_transport_pathogen <- function(Pathogen) {
  p<-Pathogen
  p$Beta<-0;p$RecoveryProb<-0;p$ProgressionProb<-0;p$PathogenMortalityProb<-0;p$ImmunityLossProb<-0;p$IntroductionProb<-0;p$IntroductionNumber<-0
  p
}
.inapest_tm_capacity_check <- function(N,nodeK,node.seedbankK,weights,tol=1e-9,label="RK host state") {
  N<-as.matrix(N);nn<-nrow(N);D<-ncol(N);K<-rep_len(as.numeric(nodeK),nn);Ks<-rep_len(as.numeric(node.seedbankK),nn)
  W<-if(is.matrix(weights))weights else matrix(rep(as.numeric(weights),each=nn),nn,D)
  if(!identical(dim(W),c(nn,D))||any(!is.finite(W))||any(W<=0))stop("weights must resolve to positive nodes x stages values")
  bad1<-which(N[,1]>Ks+tol*pmax(1,abs(Ks)))
  if(length(bad1))stop(label," exceeded node.seedbankK at node(s): ",paste(bad1,collapse=", "))
  if(D>1L){pop<-rowSums(N[,2:D,drop=FALSE]*W[,2:D,drop=FALSE]);bad<-which(pop>K+tol*pmax(1,abs(K)));if(length(bad))stop(label," exceeded weighted nodeK at node(s): ",paste(bad,collapse=", "))}
  invisible(TRUE)
}

.inapest_tm_rate_model <- function(t,State,HostFluxFunction,HostParameters,PathogenRates=NULL,Pathogen=NULL,
                                   timestep,Ntimesteps,n_nodes,D,StageMixing=NULL,AttackHazardsUnit=NULL,runtime=list()) {
  comps<-colnames(State);nu<-nrow(State);k<-ncol(State);if(nu!=n_nodes*D)stop("TM unit count mismatch")
  N<-rowSums(State);Nmat<-.inapest_tm_unflat(N,n_nodes,D);hf<-.inapest_tm_host_flux(HostFluxFunction,t,Nmat,HostParameters,runtime)
  gain<-matrix(0,nu,k,dimnames=dimnames(State));gain[,if("S"%in%comps)"S"else comps[1L]]<-hf$GainRate
  trans<-array(0,c(nu,k,k),dimnames=list(rownames(State),comps,comps));pmort<-rep(0,nu)
  if(!is.null(PathogenRates)) {
    if(!all(c("S","I")%in%comps))stop("Coupled pathogen state must contain S and I")
    rec<-.inapest_tm_resolve_rate(PathogenRates$RecoveryRate,timestep,n_nodes,D,Ntimesteps,"RecoveryRate")
    prog<-.inapest_tm_resolve_rate(PathogenRates$ProgressionRate,timestep,n_nodes,D,Ntimesteps,"ProgressionRate")
    pmort<-.inapest_tm_resolve_rate(PathogenRates$PathogenMortalityRate,timestep,n_nodes,D,Ntimesteps,"PathogenMortalityRate")
    wan<-.inapest_tm_resolve_rate(PathogenRates$ImmunityLossRate,timestep,n_nodes,D,Ntimesteps,"ImmunityLossRate")
    beta<-.inapest_tm_resolve_rate(Pathogen$Beta,timestep,n_nodes,D,Ntimesteps,"Beta")
    ds<-.inapest_tm_resolve_rate(Pathogen$DensityScale,timestep,n_nodes,D,Ntimesteps,"DensityScale");if(any(ds<=0))stop("DensityScale must be positive")
    C<-.inapest_tm_pathogen_contact(Pathogen,timestep,n_nodes,D,Ntimesteps,StageMixing);I<-State[,"I"];ip<-as.numeric(crossprod(I,C))
    if(Pathogen$Transmission=="frequency"){den<-as.numeric(crossprod(N,C));lambda<-beta*ifelse(den>0,ip/den,0)}else lambda<-beta*ip/ds
    lambda<-pmax(0,lambda)
    if("E"%in%comps){trans[,"S","E"]<-lambda;trans[,"E","I"]<-prog}else{if(any(prog>0))stop("ProgressionRate must be zero for SIS/SIR");trans[,"S","I"]<-lambda}
    if(Pathogen$Model=="SIS"){trans[,"I","S"]<-rec;if(any(wan>0))stop("ImmunityLossRate must be zero for SIS")}else{trans[,"I","R"]<-rec;trans[,"R","S"]<-wan}
  }
  hc<-colnames(hf$Hazards);bio<-if(is.null(AttackHazardsUnit))character()else colnames(AttackHazardsUnit)
  causes<-c(if(length(hc))paste0("host:",hc)else character(),if(!is.null(PathogenRates))"pathogen"else character(),paste0("bio:",bio))
  exits<-array(0,c(nu,k,length(causes)),dimnames=list(rownames(State),comps,causes))
  if(length(hc))for(j in seq_along(hc))exits[,,paste0("host:",hc[j])]<-matrix(rep(hf$Hazards[,j],k),nu,k)
  if(!is.null(PathogenRates))exits[,"I","pathogen"]<-pmort
  if(length(bio))for(j in seq_along(bio))exits[,,paste0("bio:",bio[j])]<-matrix(rep(AttackHazardsUnit[,j],k),nu,k)
  list(GainRate=gain,TransitionHazards=trans,ExitHazards=exits)
}

# Biocontrol local-rate helpers, inherited in semantics from frozen Meta v0.5.
.inapest_tm_cbr_call <- function(x,t,timestep,context,agent,State=NULL){if(!is.function(x))return(x);fm<-names(formals(x));a<-list(t=t,time=t,timestep=timestep,context=context,agent=agent,state=State,State=State);if(!is.null(fm)&&!("..."%in%fm))a<-a[intersect(names(a),fm)];do.call(x,a)}
.inapest_tm_cbr_transition_rates <- function(x,t,timestep,context,agent,State){s<-length(agent$Stages);if(is.null(x))return(matrix(0,s,s,dimnames=list(agent$Stages,agent$Stages)));x<-.inapest_tm_cbr_call(x,t,timestep,context,agent,State);if(!is.matrix(x)||!identical(dim(x),c(s,s))||any(!is.finite(x))||any(x<0))stop("Continuous agent TransitionRates must be a finite non-negative stages x stages matrix");diag(x)<-0;dimnames(x)<-list(agent$Stages,agent$Stages);x}
.inapest_tm_cbr_exit_rates <- function(x,t,timestep,context,agent,State){n<-context$n_nodes;s<-length(agent$Stages);x<-.inapest_tm_cbr_call(x,t,timestep,context,agent,State);d<-dim(x);if(is.null(d)){z<-as.numeric(x);if(length(z)==1L)out<-matrix(z,n,s)else if(length(z)==s)out<-matrix(rep(z,each=n),n,s)else if(length(z)==n&&n!=s)out<-matrix(rep(z,s),n,s)else stop("Continuous agent ExitRates shape invalid")}else if(length(d)==2L&&identical(d,c(n,s)))out<-as.matrix(x)else stop("Continuous agent ExitRates shape invalid");if(any(!is.finite(out))||any(out<0))stop("Continuous agent ExitRates must be finite and non-negative");dimnames(out)<-list(rownames(State),agent$Stages);out}
.inapest_tm_cbr_agent_rate_function <- function(spec,agent,timestep,context){force(spec);force(agent);force(timestep);force(context);function(t,State,...){State<-as.matrix(State);n<-nrow(State);s<-ncol(State);tr<-.inapest_tm_cbr_transition_rates(spec$TransitionRates,t,timestep,context,agent,State);T<-array(0,c(n,s,s),dimnames=list(rownames(State),agent$Stages,agent$Stages));for(src in seq_len(s))for(dst in seq_len(s))if(src!=dst)T[,src,dst]<-tr[dst,src];ex<-.inapest_tm_cbr_exit_rates(spec$ExitRates,t,timestep,context,agent,State);E<-array(ex,c(n,s,1L),dimnames=list(rownames(State),agent$Stages,"agent_mortality"));list(GainRate=matrix(0,n,s,dimnames=dimnames(State)),TransitionHazards=T,ExitHazards=E)}}
.inapest_tm_cbr_attacker_exposure <- function(State,Time,Step,RateFunction,AttackStage){State<-as.matrix(State);n<-nrow(State);s<-ncol(State);z0<-c(as.numeric(State),rep(0,n));deriv<-function(t,state,...){X<-matrix(state[seq_len(n*s)],n,s,dimnames=dimnames(State));dX<-INApestCompartmentMeanDerivative(t,X,RateFunction);c(as.numeric(dX),as.numeric(X[,as.integer(AttackStage)]))};zend<-INApestRK4Step(z0,Time=Time,Step=Step,RateFunction=deriv,NonNegative="allow");e<-as.numeric(zend[n*s+seq_len(n)]);if(any(!is.finite(e))||any(e< -1e-8))stop("Continuous biocontrol attacker exposure left valid domain");pmax(0,e)}
.inapest_tm_target_stages <- function(agent,context){ts<-agent$TargetStage;if(is.character(ts))ts<-match(ts,context$host_stages);ts<-as.integer(ts);if(!length(ts)||anyNA(ts)||any(ts<1L|ts>length(context$host_stages)))stop("Invalid TargetStage for continuous TM biocontrol agent ",agent$Name);unique(ts)}

# Stable environment marker used by the parallel engine to export shared TM RK helpers.
INApestTMRK4HelpersMarker <- function() invisible(TRUE)
