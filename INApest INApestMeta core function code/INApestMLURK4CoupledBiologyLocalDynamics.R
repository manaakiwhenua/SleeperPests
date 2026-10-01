###############################################################################
# INApestMLURK4CoupledBiologyLocalDynamics.R
# MLU stochastic continuous-time host + biocontrol (+ optional pathogen)
# v0.6 -- 30 September 2026
###############################################################################
if(!exists("INApestStochasticCompartmentStep",mode="function")||!exists("INApestCompartmentMeanDerivative",mode="function")||!exists("INApestRK4Step",mode="function"))
  stop("Source INApestRK4.R and INApestStochasticCompartmentBridge.R first")
if(!exists("INApestBiocontrol",mode="function")) stop("Source INApestBiocontrol.R first")
if(!exists("INApestContinuousPathogenRates",mode="function")) {
  INApestContinuousPathogenRates<-function(RecoveryRate=0,ProgressionRate=0,PathogenMortalityRate=0,ImmunityLossRate=0)
    structure(list(RecoveryRate=RecoveryRate,ProgressionRate=ProgressionRate,PathogenMortalityRate=PathogenMortalityRate,ImmunityLossRate=ImmunityLossRate),class=c("INApestContinuousPathogenRates","list"))
}
if(!exists(".inapest_mlu_expand",mode="function")) {
  .inapest_mlu_expand<-function(x,template,label){d<-dim(template);n<-length(template);if(is.null(dim(x))){z<-as.numeric(x);if(length(z)==1L)return(matrix(z,d[1],d[2]));if(length(z)==n)return(matrix(z,d[1],d[2]))}else if(identical(dim(x),d))return(matrix(as.numeric(x),d[1],d[2]));stop(label," must be scalar, length cells, or node x land-use matrix")}
}
if(!exists(".inapest_mlu_host_flux",mode="function")) {
  .inapest_mlu_host_flux<-function(fun,t,Nmat,Parameters=NULL,runtime=list()){ans<-.inapest_rk4_call_supported(fun,c(list(t=t,time=t,state=Nmat,State=Nmat,pars=Parameters,Parameters=Parameters),runtime));if(!is.list(ans))stop("HostFluxFunction must return GainRate and Hazards");g<-ans$GainRate;if(is.null(g))g<-ans$Gains;if(is.null(g))stop("HostFluxFunction must return GainRate");g<-.inapest_mlu_expand(g,Nmat,"Host GainRate");if(any(!is.finite(g))||any(g<0))stop("Invalid Host GainRate");H<-ans$Hazards;if(is.null(H)||(is.list(H)&&length(H)==0L))return(list(GainRate=c(g),Hazards=matrix(numeric(),length(Nmat),0L)));if(!is.list(H)||is.null(names(H))||any(!nzchar(names(H)))||anyDuplicated(names(H)))stop("Host Hazards must be uniquely named list");m<-vapply(H,function(z)c(.inapest_mlu_expand(z,Nmat,"Host hazard")),numeric(length(Nmat)));if(is.null(dim(m)))m<-matrix(m,ncol=1L,dimnames=list(NULL,names(H)));colnames(m)<-names(H);if(any(!is.finite(m))||any(m<0))stop("Invalid Host Hazards");list(GainRate=c(g),Hazards=m)}
}

INApestContinuousBiocontrolAgentRates <- function(TransitionRates=NULL,ExitRates=0)
  structure(list(TransitionRates=TransitionRates,ExitRates=ExitRates),class=c("INApestContinuousBiocontrolAgentRates","list"))
INApestContinuousBiocontrolRates <- function(AgentRates){
  if(!is.list(AgentRates)||!length(AgentRates)||is.null(names(AgentRates))||any(!nzchar(names(AgentRates)))||anyDuplicated(names(AgentRates)))stop("AgentRates must be uniquely named non-empty list")
  if(!all(vapply(AgentRates,inherits,logical(1),"INApestContinuousBiocontrolAgentRates")))stop("Every AgentRates entry must come from INApestContinuousBiocontrolAgentRates()")
  structure(AgentRates,class=c("INApestContinuousBiocontrolRates","list"))
}
.inapest_cbr_call<-function(x,t,timestep,context,agent,State=NULL){if(!is.function(x))return(x);fm<-names(formals(x));a<-list(t=t,time=t,timestep=timestep,context=context,agent=agent,state=State,State=State);if(!is.null(fm)&&!("..."%in%fm))a<-a[intersect(names(a),fm)];do.call(x,a)}
.inapest_cbr_transition_rates<-function(x,t,timestep,context,agent,State){s<-length(agent$Stages);if(is.null(x))return(matrix(0,s,s,dimnames=list(agent$Stages,agent$Stages)));x<-.inapest_cbr_call(x,t,timestep,context,agent,State);if(length(dim(x))!=2L||!identical(dim(x),c(s,s)))stop("Continuous agent TransitionRates must be stages x stages");out<-as.matrix(x);if(any(!is.finite(out))||any(out<0))stop("Invalid continuous TransitionRates");diag(out)<-0;dimnames(out)<-list(agent$Stages,agent$Stages);out}
.inapest_cbr_exit_rates<-function(x,t,timestep,context,agent,State){n<-context$n_nodes;s<-length(agent$Stages);x<-.inapest_cbr_call(x,t,timestep,context,agent,State);d<-dim(x);if(is.null(d)){z<-as.numeric(x);if(length(z)==1L)out<-matrix(z,n,s)else if(length(z)==s)out<-matrix(rep(z,each=n),n,s)else if(length(z)==n&&n!=s)out<-matrix(rep(z,s),n,s)else stop("Continuous agent ExitRates shape unsupported")}else if(length(d)==2L&&identical(d,c(n,s)))out<-as.matrix(x)else stop("Continuous agent ExitRates shape unsupported");if(any(!is.finite(out))||any(out<0))stop("Invalid continuous ExitRates");dimnames(out)<-list(rownames(State),agent$Stages);out}
.inapest_cbr_agent_rate_function<-function(spec,agent,timestep,context){force(spec);force(agent);force(timestep);force(context);function(t,State,...){State<-as.matrix(State);n<-nrow(State);s<-ncol(State);tr<-.inapest_cbr_transition_rates(spec$TransitionRates,t,timestep,context,agent,State);T<-array(0,c(n,s,s),dimnames=list(rownames(State),agent$Stages,agent$Stages));for(src in seq_len(s))for(dst in seq_len(s))if(src!=dst)T[,src,dst]<-tr[dst,src];ex<-.inapest_cbr_exit_rates(spec$ExitRates,t,timestep,context,agent,State);E<-array(ex,c(n,s,1L),dimnames=list(rownames(State),agent$Stages,"agent_mortality"));list(GainRate=matrix(0,n,s,dimnames=dimnames(State)),TransitionHazards=T,ExitHazards=E)}}
.inapest_cbr_attacker_exposure<-function(State,Time,Step,RateFunction,AttackStage){State<-as.matrix(State);n<-nrow(State);s<-ncol(State);z0<-c(as.numeric(State),rep(0,n));deriv<-function(t,state,...){X<-matrix(state[seq_len(n*s)],n,s,dimnames=dimnames(State));dX<-INApestCompartmentMeanDerivative(t,X,RateFunction);c(as.numeric(dX),as.numeric(X[,as.integer(AttackStage)]))};zend<-INApestRK4Step(z0,Time=Time,Step=Step,RateFunction=deriv,NonNegative="allow");e<-as.numeric(zend[n*s+seq_len(n)]);if(any(!is.finite(e))||any(e< -1e-8))stop("Invalid attacker exposure");pmax(0,e)}

.inapest_mlu_coupled_host_rate<-function(t,State,HostFluxFunction,HostParameters,AttackHazardsNode,PathogenRates=NULL,Pathogen=NULL,PathogenEngine=NULL,PathogenContext=NULL,timestep,Ntimesteps,runtime=list()){
  comps<-colnames(State);nu<-nrow(State);k<-ncol(State);nn<-if(!is.null(PathogenContext))PathogenContext$n_nodes else nrow(runtime$nodeK);nl<-if(!is.null(PathogenContext))PathogenContext$n_landuses else ncol(runtime$nodeK)
  N<-rowSums(State);Nmat<-matrix(N,nn,nl);hf<-.inapest_mlu_host_flux(HostFluxFunction,t,Nmat,HostParameters,runtime)
  gain<-matrix(0,nu,k,dimnames=dimnames(State));gain[,if("S"%in%comps)"S" else comps[1L]]<-hf$GainRate
  trans<-array(0,c(nu,k,k),dimnames=list(rownames(State),comps,comps));pmort<-rep(0,nu)
  if(!is.null(PathogenRates)){
    rr<-function(x,name){z<-as.numeric(PathogenEngine$Resolve(x,timestep,PathogenContext,name));if(length(z)!=nu||any(!is.finite(z))||any(z<0))stop(name," must resolve per MLU cell");z}
    rec<-rr(PathogenRates$RecoveryRate,"RecoveryRate");prog<-rr(PathogenRates$ProgressionRate,"ProgressionRate");pmort<-rr(PathogenRates$PathogenMortalityRate,"PathogenMortalityRate");wan<-rr(PathogenRates$ImmunityLossRate,"ImmunityLossRate");beta<-rr(Pathogen$Beta,"Beta");ds<-rr(Pathogen$DensityScale,"DensityScale");if(any(ds<=0))stop("DensityScale must be positive")
    C<-PathogenEngine$ContactMatrix(timestep,PathogenContext);I<-State[,"I"];ip<-as.numeric(crossprod(I,C));if(Pathogen$Transmission=="frequency"){cp<-as.numeric(crossprod(N,C));lambda<-beta*ifelse(cp>0,ip/cp,0)}else lambda<-beta*ip/ds;lambda<-pmax(0,lambda)
    if("E"%in%comps){trans[,"S","E"]<-lambda;trans[,"E","I"]<-prog}else{if(any(prog>0))stop("ProgressionRate must be zero for SIS/SIR");trans[,"S","I"]<-lambda}
    if(Pathogen$Model=="SIS"){trans[,"I","S"]<-rec;if(any(wan>0))stop("ImmunityLossRate must be zero for SIS")}else{trans[,"I","R"]<-rec;trans[,"R","S"]<-wan}
  }
  bio_names<-colnames(AttackHazardsNode);attack_cell<-vapply(seq_along(bio_names),function(j)rep(AttackHazardsNode[,j],times=nl),numeric(nu));if(is.null(dim(attack_cell)))attack_cell<-matrix(attack_cell,ncol=1L);colnames(attack_cell)<-bio_names
  hc<-colnames(hf$Hazards);causes<-c(if(length(hc))paste0("host:",hc)else character(),if(!is.null(PathogenRates))"pathogen"else character(),paste0("bio:",bio_names))
  exits<-array(0,c(nu,k,length(causes)),dimnames=list(rownames(State),comps,causes));if(length(hc))for(j in seq_along(hc))exits[,,paste0("host:",hc[j])]<-matrix(rep(hf$Hazards[,j],k),nu,k);if(!is.null(PathogenRates))exits[,"I","pathogen"]<-pmort;for(j in seq_along(bio_names))exits[,,paste0("bio:",bio_names[j])]<-matrix(rep(attack_cell[,j],k),nu,k)
  list(GainRate=gain,TransitionHazards=trans,ExitHazards=exits)
}

INApestMLURK4CoupledBiologyLocalDynamics<-function(HostFluxFunction,BiocontrolRates,PathogenRates=NULL,TimestepLength=1,RKMaxStep=0.025,HostParameters=NULL,StartTime=0,DispersalDynamics=NULL,CapacityTolerance=1e-9){
  if(!is.function(HostFluxFunction)||!inherits(BiocontrolRates,"INApestContinuousBiocontrolRates"))stop("Invalid coupled-biology specification")
  if(!is.null(PathogenRates)&&!inherits(PathogenRates,"INApestContinuousPathogenRates"))stop("PathogenRates must be NULL or continuous pathogen rates")
  if(is.null(DispersalDynamics)){if(!exists("local.dynamicsLU",mode="function"))stop("Source current MLU engine first");DispersalDynamics<-get("local.dynamicsLU",mode="function")}
  host_fun<-HostFluxFunction;brates<-BiocontrolRates;prates<-PathogenRates;pars<-HostParameters;dt<-as.numeric(TimestepLength)[1L];hmax<-as.numeric(RKMaxStep)[1L];t0<-as.numeric(StartTime)[1L];disperse<-DispersalDynamics;cap_tol<-as.numeric(CapacityTolerance)[1L];couples_pathogen<-!is.null(prates)
  if(!is.finite(dt)||dt<=0||!is.finite(hmax)||hmax<=0||!is.finite(t0)||!is.finite(cap_tol)||cap_tol<0)stop("Invalid RK settings")
  f<-function(sddprob,nodepropaguleproduction,nodeenvestabprob,n,lddprob,lddrate,k_is_0,nodeK,nodepropaguleestablishment,nodespreadreduction,nodefecundityreduction=0,managing,biocontrol_state,biocontrol,biocontrol_engine,biocontrol_context,pathogen_state=NULL,pathogen=NULL,pathogen_engine=NULL,pathogen_context=NULL,timestep=NULL,Ntimesteps=NULL){
    if(is.null(biocontrol)||!inherits(biocontrol,"INApestBiocontrol"))stop("Coupled MLU biology requires Biocontrol")
    nn<-nrow(n);nl<-ncol(n);nu<-length(n);agents<-biocontrol$Agents;if(!identical(names(brates),names(agents)))stop("BiocontrolRates names must match agents");if(!is.list(biocontrol_state)||!identical(names(biocontrol_state),names(agents)))stop("biocontrol_state mismatch")
    si<-if(is.null(timestep))1L else as.integer(timestep)[1L];nts<-if(is.null(Ntimesteps))1L else as.integer(Ntimesteps)[1L]
    if(couples_pathogen){if(is.null(pathogen)||!inherits(pathogen,"INApestPathogen")||pathogen$Model=="Binary")stop("H+P+B requires non-binary Pathogen");if(!is.matrix(pathogen_state)||nrow(pathogen_state)!=nu||any(rowSums(pathogen_state)!=as.integer(c(n))))stop("pathogen_state must sum to MLU host abundance");if(any(pathogen_engine$Resolve(pathogen$IntroductionProb,si,pathogen_context,"IntroductionProb")!=0))stop("Continuous MLU H+P+B v0.6 requires Pathogen$IntroductionProb = 0");Hstate<-pathogen_state}else Hstate<-matrix(as.integer(c(n)),ncol=1L,dimnames=list(NULL,"H"))
    Astate<-biocontrol_state
    for(nm in names(agents)){a<-agents[[nm]];z<-as.matrix(Astate[[nm]]);rel<-.inabc_resolve_state_matrix(a$Release,si,biocontrol_context,a,"Release");z<-z+matrix(as.integer(round(rel)),nrow(z),ncol(z),dimnames=dimnames(z));storage.mode(z)<-"integer";Astate[[nm]]<-z}
    attacks_cell<-lapply(agents,function(a)integer(nu));names(attacks_cell)<-names(agents);nsub<-max(1L,as.integer(ceiling(dt/hmax-1e-14)));h<-dt/nsub;tt<-t0+(si-1L)*dt
    runtime<-list(nodeK=nodeK,k_is_0=k_is_0,nodeenvestabprob=nodeenvestabprob,nodepropaguleproduction=nodepropaguleproduction,nodepropaguleestablishment=nodepropaguleestablishment,nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,managing=managing,sddprob=sddprob,lddprob=lddprob,lddrate=lddrate,timestep=si,Ntimesteps=nts)
    for(ss in seq_len(nsub)){
      ah<-matrix(0,nn,length(agents),dimnames=list(NULL,names(agents)));next_agents<-Astate
      for(nm in names(agents)){a<-agents[[nm]];rf<-.inapest_cbr_agent_rate_function(brates[[nm]],a,si,biocontrol_context);expo<-.inapest_cbr_attacker_exposure(Astate[[nm]],tt,h,rf,a$AttackStage);ar<-.inabc_resolve_node(a$AttackRate,si,biocontrol_context,"AttackRate",a);ah[,nm]<-ar*(expo/h);az<-INApestStochasticCompartmentStep(Astate[[nm]],tt,h,rf);next_agents[[nm]]<-az$State}
      hfun<-function(t,state,...) .inapest_mlu_coupled_host_rate(t,state,host_fun,pars,ah,if(couples_pathogen)prates else NULL,if(couples_pathogen)pathogen else NULL,if(couples_pathogen)pathogen_engine else NULL,if(couples_pathogen)pathogen_context else NULL,si,nts,runtime)
      hz<-INApestStochasticCompartmentStep(Hstate,tt,h,hfun);Hstate<-hz$State
      for(nm in names(agents)){cause<-paste0("bio:",nm);killed<-if(cause%in%colnames(hz$Exits))as.integer(hz$Exits[,cause])else integer(nu);attacks_cell[[nm]]<-attacks_cell[[nm]]+killed;killmat<-matrix(killed,nn,nl);recnode<-rowSums(killmat)*as.integer(agents[[nm]]$OffspringPerAttack);next_agents[[nm]][,agents[[nm]]$RecruitStage]<-next_agents[[nm]][,agents[[nm]]$RecruitStage]+recnode}
      Astate<-next_agents;tt<-tt+h
    }
    for(nm in names(agents)){a<-agents[[nm]];M<-.inabc_resolve_movement(a$Movement,si,biocontrol_context,a);Astate[[nm]]<-.inabc_move_agent(Astate[[nm]],M,a$MovementStages);storage.mode(Astate[[nm]])<-"integer"}
    n_after<-matrix(as.integer(rowSums(Hstate)),nn,nl);over<-which(n_after>nodeK+cap_tol*pmax(1,abs(nodeK)),arr.ind=TRUE);if(nrow(over))stop("Coupled RK host state exceeded MLU nodeK before dispersal")
    n_out<-disperse(sddprob=sddprob,nodepropaguleproduction=nodepropaguleproduction,nodeenvestabprob=nodeenvestabprob,n=n_after,lddprob=lddprob,lddrate=lddrate,k_is_0=k_is_0,nodeK=nodeK,nodepropaguleestablishment=nodepropaguleestablishment,nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,managing=managing)
    Hout<-if(couples_pathogen)pathogen_engine$Reconcile(Hstate,n_out,pathogen_context)else NULL
    impact_att<-lapply(names(agents),function(nm)matrix(attacks_cell[[nm]],nn,nl,dimnames=list(NULL,paste0("landuse_",seq_len(nl)))));names(impact_att)<-names(agents)
    rec_node<-lapply(names(agents),function(nm)as.integer(rowSums(matrix(attacks_cell[[nm]],nn,nl))*agents[[nm]]$OffspringPerAttack));names(rec_node)<-names(agents);rec_total<-vapply(rec_node,sum,integer(1))
    impact<-list(AttacksByAgent=impact_att,RecruitsByAgent=rec_total,RecruitsByAgentNode=rec_node)
    out<-list(N=matrix(as.integer(n_out),nn,nl),BiocontrolState=Astate,BiocontrolImpact=impact,CoupledRK=list(PreDispersalHost=n_after,Attacks=attacks_cell,Timestep=si,Ntimesteps=nts,NSubsteps=nsub,InternalStep=h));if(couples_pathogen)out$PathogenState<-Hout;out
  }
  attr(f,"INApestMLURK4CoupledBiologyLocalDynamics")<-list(Version="0.6",StateMode="integer-coupled",TimestepLength=dt,RKMaxStep=hmax,CouplesPathogen=couples_pathogen)
  attr(f,"INApestCouplesBiocontrol")<-TRUE;attr(f,"INApestCouplesPathogen")<-couples_pathogen;attr(f,"INApestRunWhenEmpty")<-TRUE
  class(f)<-c("INApestMLURK4CoupledBiologyLocalDynamics","function");f
}
