###############################################################################
# INApestMetaRK4CoupledBiologyLocalDynamics.R
# Meta stochastic continuous-time coupled host + pathogen + biocontrol biology
# v0.5 -- 30 September 2026
#
# Scope:
#   * H + B or H + P + B local biology in one stochastic RK framework.
#   * Host demography is supplied as gross gains + per-capita loss hazards.
#   * Pathogen SIS/SIR/SEIR uses linked compartment transfers from v0.4.
#   * Biocontrol attack is an explicit competing host-exit cause. Each realised
#     attack therefore removes exactly one host (and, when a pathogen is active,
#     the correct S/E/I/R individual) and creates the corresponding integer
#     biocontrol recruits.
#   * Biocontrol stage transitions/mortality are continuous-rate stochastic
#     compartment dynamics. Scheduled releases and spatial movement retain the
#     established INApest companion ordering at parent-timestep boundaries.
#
# This file does not alter the definitive discrete INApestBiocontrol.R companion.
###############################################################################

if (!exists("INApestStochasticCompartmentStep", mode="function") ||
    !exists("INApestCompartmentMeanDerivative", mode="function") ||
    !exists("INApestRK4Step", mode="function"))
  stop("Source INApestRK4.R and INApestStochasticCompartmentBridge.R first")

INApestContinuousBiocontrolAgentRates <- function(
    TransitionRates = NULL,
    ExitRates = 0) {
  structure(list(TransitionRates=TransitionRates, ExitRates=ExitRates),
            class=c("INApestContinuousBiocontrolAgentRates","list"))
}

INApestContinuousBiocontrolRates <- function(AgentRates) {
  if (inherits(AgentRates,"INApestContinuousBiocontrolAgentRates"))
    stop("AgentRates must be a named list, one entry per biocontrol agent")
  if (!is.list(AgentRates) || !length(AgentRates) || is.null(names(AgentRates)) ||
      any(!nzchar(names(AgentRates))) || anyDuplicated(names(AgentRates)))
    stop("AgentRates must be a uniquely named non-empty list")
  if (!all(vapply(AgentRates,inherits,logical(1),"INApestContinuousBiocontrolAgentRates")))
    stop("Every AgentRates entry must come from INApestContinuousBiocontrolAgentRates()")
  structure(AgentRates,class=c("INApestContinuousBiocontrolRates","list"))
}

.inapest_cbr_call <- function(x, t, timestep, context, agent, State=NULL) {
  if (!is.function(x)) return(x)
  fm <- names(formals(x))
  a <- list(t=t,time=t,timestep=timestep,context=context,agent=agent,
            state=State,State=State)
  if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a),fm)]
  do.call(x,a)
}

# Continuous stage-transition rates use the same visual orientation as the
# definitive discrete agent Transition matrix: rows=destination, columns=source.
.inapest_cbr_transition_rates <- function(x, t, timestep, context, agent, State) {
  s <- length(agent$Stages)
  if (is.null(x)) return(matrix(0,s,s,dimnames=list(agent$Stages,agent$Stages)))
  x <- .inapest_cbr_call(x,t,timestep,context,agent,State)
  d <- dim(x)
  if (length(d)==2L && identical(d,c(s,s))) out <- as.matrix(x)
  else stop("Continuous agent TransitionRates must be a stages x stages matrix or resolver function")
  if (any(!is.finite(out)) || any(out<0)) stop("Continuous agent TransitionRates must be finite and non-negative")
  diag(out) <- 0
  dimnames(out) <- list(agent$Stages,agent$Stages)
  out
}

.inapest_cbr_exit_rates <- function(x, t, timestep, context, agent, State) {
  n <- context$n_nodes; s <- length(agent$Stages)
  x <- .inapest_cbr_call(x,t,timestep,context,agent,State)
  d <- dim(x)
  if (is.null(d)) {
    z <- as.numeric(x)
    if (length(z)==1L) out <- matrix(z,n,s)
    else if (length(z)==s) out <- matrix(rep(z,each=n),n,s)
    else if (length(z)==n && n!=s) out <- matrix(rep(z,s),n,s)
    else stop("Continuous agent ExitRates must be scalar, length stages, nodes x stages, or resolver function")
  } else if (length(d)==2L && identical(d,c(n,s))) out <- as.matrix(x)
  else stop("Continuous agent ExitRates must be scalar, length stages, nodes x stages, or resolver function")
  if (any(!is.finite(out)) || any(out<0)) stop("Continuous agent ExitRates must be finite and non-negative")
  dimnames(out) <- list(rownames(State),agent$Stages)
  out
}

.inapest_cbr_agent_rate_function <- function(spec, agent, timestep, context) {
  force(spec); force(agent); force(timestep); force(context)
  function(t, State, ...) {
    State <- as.matrix(State); n <- nrow(State); s <- ncol(State)
    tr <- .inapest_cbr_transition_rates(spec$TransitionRates,t,timestep,context,agent,State)
    T <- array(0,c(n,s,s),dimnames=list(rownames(State),agent$Stages,agent$Stages))
    # stored matrix is destination x source; bridge wants source x destination.
    for (src in seq_len(s)) for (dst in seq_len(s)) if (src!=dst)
      T[,src,dst] <- tr[dst,src]
    ex <- .inapest_cbr_exit_rates(spec$ExitRates,t,timestep,context,agent,State)
    E <- array(ex,c(n,s,1L),dimnames=list(rownames(State),agent$Stages,"agent_mortality"))
    list(GainRate=matrix(0,n,s,dimnames=dimnames(State)),
         TransitionHazards=T, ExitHazards=E)
  }
}

# Deterministic RK4 integral of the attacking-stage abundance over one internal
# interval. The resulting count*time exposure drives the continuous attack
# hazard while the realised agent endpoint remains an integer stochastic draw.
.inapest_cbr_attacker_exposure <- function(State, Time, Step, RateFunction, AttackStage) {
  State <- as.matrix(State); n <- nrow(State); s <- ncol(State)
  AttackStage <- as.integer(AttackStage)
  z0 <- c(as.numeric(State),rep(0,n))
  deriv <- function(t,state,...) {
    z <- state
    X <- matrix(z[seq_len(n*s)],n,s,dimnames=dimnames(State))
    dX <- INApestCompartmentMeanDerivative(t,X,RateFunction)
    c(as.numeric(dX),as.numeric(X[,AttackStage]))
  }
  zend <- INApestRK4Step(z0,Time=Time,Step=Step,RateFunction=deriv,NonNegative="allow")
  exposure <- as.numeric(zend[n*s + seq_len(n)])
  if (any(!is.finite(exposure)) || any(exposure < -1e-8))
    stop("Continuous biocontrol attacker exposure left the valid biological domain")
  pmax(0,exposure)
}

.inapest_cbr_host_flux <- function(fun,t,N,Parameters=NULL,runtime=list()) {
  ans <- .inapest_rk4_call_supported(
    fun,c(list(t=t,time=t,state=N,State=N,pars=Parameters,Parameters=Parameters),runtime)
  )
  if (!is.list(ans)) stop("HostFluxFunction must return GainRate and Hazards")
  gain <- ans$GainRate; if (is.null(gain)) gain <- ans$Gains
  if (is.null(gain)) stop("HostFluxFunction must return GainRate")
  gain <- as.numeric(gain); if(length(gain)==1L) gain <- rep(gain,length(N))
  if(length(gain)!=length(N)||any(!is.finite(gain))||any(gain<0))
    stop("Host GainRate must resolve to finite non-negative values per node")
  H <- ans$Hazards
  if (is.null(H) || (is.list(H)&&length(H)==0L))
    return(list(GainRate=gain,Hazards=matrix(numeric(),length(N),0L)))
  if (is.list(H)) {
    if(is.null(names(H))||any(!nzchar(names(H)))||anyDuplicated(names(H)))
      stop("Host hazard list must have unique non-empty names")
    mat <- vapply(H,function(z){z<-as.numeric(z);if(length(z)==1L)z<-rep(z,length(N));if(length(z)!=length(N))stop("Each host hazard must be scalar or length nodes");z},numeric(length(N)))
    if(is.null(dim(mat))) mat<-matrix(mat,ncol=1L,dimnames=list(NULL,names(H)))
    colnames(mat)<-names(H)
  } else {
    mat<-as.matrix(H);if(nrow(mat)==1L&&length(N)>1L)mat<-mat[rep(1L,length(N)),,drop=FALSE]
    if(nrow(mat)!=length(N))stop("Host Hazards must have one row per node")
    if(is.null(colnames(mat)))colnames(mat)<-paste0("loss",seq_len(ncol(mat)))
  }
  if(any(!is.finite(mat))||any(mat<0))stop("Host Hazards must be finite and non-negative")
  list(GainRate=gain,Hazards=mat)
}

.inapest_cbr_host_rate_model <- function(t, State, HostFluxFunction, HostParameters,
                                         AttackHazards,
                                         PathogenRates=NULL, Pathogen=NULL,
                                         PathogenEngine=NULL, PathogenContext=NULL,
                                         timestep, Ntimesteps, runtime=list()) {
  comps <- colnames(State); n <- nrow(State); k <- ncol(State); N <- rowSums(State)
  hf <- .inapest_cbr_host_flux(HostFluxFunction,t,N,HostParameters,runtime)
  gain <- matrix(0,n,k,dimnames=dimnames(State))
  gain[,if("S" %in% comps) "S" else comps[1L]] <- hf$GainRate
  trans <- array(0,c(n,k,k),dimnames=list(rownames(State),comps,comps))
  pmort <- rep(0,n)

  if (!is.null(PathogenRates)) {
    if (!all(c("S","I") %in% comps)) stop("Coupled pathogen state must contain S and I")
    resolve_rate <- function(x,name) {
      z <- as.numeric(PathogenEngine$Resolve(x,timestep,PathogenContext,name))
      if(length(z)!=n||any(!is.finite(z))||any(z<0)) stop(name," must resolve to finite non-negative rates per node")
      z
    }
    rec <- resolve_rate(PathogenRates$RecoveryRate,"RecoveryRate")
    prog <- resolve_rate(PathogenRates$ProgressionRate,"ProgressionRate")
    pmort <- resolve_rate(PathogenRates$PathogenMortalityRate,"PathogenMortalityRate")
    wan <- resolve_rate(PathogenRates$ImmunityLossRate,"ImmunityLossRate")
    beta <- as.numeric(PathogenEngine$Resolve(Pathogen$Beta,timestep,PathogenContext,"Beta"))
    ds <- as.numeric(PathogenEngine$Resolve(Pathogen$DensityScale,timestep,PathogenContext,"DensityScale"))
    C <- PathogenEngine$ContactMatrix(timestep,PathogenContext)
    I <- State[,"I"]; ip <- as.numeric(crossprod(I,C))
    if(Pathogen$Transmission=="frequency") {
      cp <- as.numeric(crossprod(N,C)); pressure <- ifelse(cp>0,ip/cp,0); lambda <- beta*pressure
    } else lambda <- beta*ip/ds
    lambda <- pmax(0,lambda)
    if("E" %in% comps) {trans[,"S","E"]<-lambda;trans[,"E","I"]<-prog}
    else {if(any(prog>0))stop("ProgressionRate must be zero for SIS/SIR");trans[,"S","I"]<-lambda}
    if(Pathogen$Model=="SIS") {trans[,"I","S"]<-rec;if(any(wan>0))stop("ImmunityLossRate must be zero for SIS")}
    else {trans[,"I","R"]<-rec;trans[,"R","S"]<-wan}
  }

  host_causes <- colnames(hf$Hazards)
  bio_names <- colnames(AttackHazards)
  cause_names <- c(if(length(host_causes))paste0("host:",host_causes) else character(),
                   if(!is.null(PathogenRates))"pathogen" else character(),
                   paste0("bio:",bio_names))
  exits <- array(0,c(n,k,length(cause_names)),dimnames=list(rownames(State),comps,cause_names))
  if(length(host_causes)) for(j in seq_along(host_causes))
    exits[,,paste0("host:",host_causes[j])] <- matrix(rep(hf$Hazards[,j],k),n,k)
  if(!is.null(PathogenRates)) exits[,"I","pathogen"] <- pmort
  for(j in seq_along(bio_names))
    exits[,,paste0("bio:",bio_names[j])] <- matrix(rep(AttackHazards[,j],k),n,k)
  list(GainRate=gain,TransitionHazards=trans,ExitHazards=exits)
}

INApestMetaRK4CoupledBiologyLocalDynamics <- function(
    HostFluxFunction,
    BiocontrolRates,
    PathogenRates = NULL,
    TimestepLength = 1,
    RKMaxStep = 0.025,
    HostParameters = NULL,
    StartTime = 0,
    DispersalDynamics = NULL,
    CapacityTolerance = 1e-9) {
  if(!is.function(HostFluxFunction)) stop("HostFluxFunction must be a function")
  if(!inherits(BiocontrolRates,"INApestContinuousBiocontrolRates"))
    stop("BiocontrolRates must come from INApestContinuousBiocontrolRates()")
  if(!is.null(PathogenRates) && !inherits(PathogenRates,"INApestContinuousPathogenRates"))
    stop("PathogenRates must be NULL or INApestContinuousPathogenRates()")
  TimestepLength<-as.numeric(TimestepLength)[1L];RKMaxStep<-as.numeric(RKMaxStep)[1L];StartTime<-as.numeric(StartTime)[1L]
  if(!is.finite(TimestepLength)||TimestepLength<=0||!is.finite(RKMaxStep)||RKMaxStep<=0||!is.finite(StartTime))stop("Invalid timestep/RK settings")
  if(is.null(DispersalDynamics)){if(!exists("local.dynamics",mode="function"))stop("Source current INApestMeta.r first");DispersalDynamics<-get("local.dynamics",mode="function")}
  host_fun<-HostFluxFunction; brates<-BiocontrolRates; prates<-PathogenRates; pars<-HostParameters
  dt<-TimestepLength;hmax<-RKMaxStep;t0<-StartTime;disperse<-DispersalDynamics;cap_tol<-CapacityTolerance
  couples_pathogen <- !is.null(prates)

  f <- function(sddprob,nodepropaguleproduction,nodeenvestabprob,n,lddprob,lddrate,
                k_is_0,nodeK,nodepropaguleestablishment,nodespreadreduction,
                nodefecundityreduction=0,managing,maxinteger,
                biocontrol_state,biocontrol,biocontrol_engine,biocontrol_context,
                pathogen_state=NULL,pathogen=NULL,pathogen_engine=NULL,pathogen_context=NULL,
                timestep=NULL,Ntimesteps=NULL) {
    if(is.null(biocontrol)||!inherits(biocontrol,"INApestBiocontrol"))stop("Coupled biology LocalDynamics requires Biocontrol")
    agents <- biocontrol$Agents
    if(!identical(names(brates),names(agents)))stop("BiocontrolRates names must exactly match Biocontrol agent names")
    if(!is.list(biocontrol_state)||!identical(names(biocontrol_state),names(agents)))stop("biocontrol_state does not match configured agents")
    step_index<-if(is.null(timestep))1L else as.integer(timestep)[1L];nts<-if(is.null(Ntimesteps))1L else as.integer(Ntimesteps)[1L]
    if(couples_pathogen) {
      if(is.null(pathogen)||!inherits(pathogen,"INApestPathogen")||pathogen$Model=="Binary")stop("H+P+B coupling requires non-binary INApestPathogen")
      if(!is.matrix(pathogen_state)||any(rowSums(pathogen_state)!=as.integer(n)))stop("pathogen_state must sum to host abundance")
      intro<-pathogen_engine$Resolve(pathogen$IntroductionProb,step_index,pathogen_context,"IntroductionProb")
      if(any(intro!=0))stop("Continuous Meta H+P+B coupling v0.5 requires Pathogen$IntroductionProb = 0")
      Hstate<-pathogen_state
    } else Hstate<-matrix(as.integer(n),ncol=1L,dimnames=list(NULL,"H"))

    Astate <- biocontrol_state
    # Scheduled releases occur at the parent boundary before attack, matching
    # the definitive companion. Continuous within-step demography follows.
    for(nm in names(agents)) {
      a<-agents[[nm]]; z<-as.matrix(Astate[[nm]])
      if(!identical(dim(z),c(biocontrol_context$n_nodes,length(a$Stages))))stop("Biocontrol state shape mismatch for ",nm)
      rel<-.inabc_resolve_state_matrix(a$Release,step_index,biocontrol_context,a,"Release")
      if(any(!is.finite(rel))||any(rel<0))stop("Biocontrol Release must be finite and non-negative")
      z<-z+matrix(as.integer(round(rel)),nrow(z),ncol(z),dimnames=dimnames(z));storage.mode(z)<-"integer";Astate[[nm]]<-z
    }

    attacks_node<-lapply(agents,function(a)integer(biocontrol_context$n_nodes));names(attacks_node)<-names(agents)
    nsub<-max(1L,as.integer(ceiling(dt/hmax-1e-14)));h<-dt/nsub;tt<-t0+(step_index-1L)*dt
    runtime<-list(nodeK=nodeK,k_is_0=k_is_0,nodeenvestabprob=nodeenvestabprob,
      nodepropaguleproduction=nodepropaguleproduction,nodepropaguleestablishment=nodepropaguleestablishment,
      nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,
      managing=managing,sddprob=sddprob,lddprob=lddprob,lddrate=lddrate,timestep=step_index,Ntimesteps=nts)

    for(ss in seq_len(nsub)) {
      attack_h<-matrix(0,biocontrol_context$n_nodes,length(agents),dimnames=list(NULL,names(agents)))
      next_agents<-Astate
      for(nm in names(agents)) {
        a<-agents[[nm]]; rf<-.inapest_cbr_agent_rate_function(brates[[nm]],a,step_index,biocontrol_context)
        exposure<-.inapest_cbr_attacker_exposure(Astate[[nm]],tt,h,rf,a$AttackStage)
        ar<-.inabc_resolve_node(a$AttackRate,step_index,biocontrol_context,"AttackRate",a)
        if(any(!is.finite(ar))||any(ar<0))stop("AttackRate must resolve finite non-negative values")
        attack_h[,nm]<-ar*(exposure/h)
        az<-INApestStochasticCompartmentStep(Astate[[nm]],tt,h,rf)
        next_agents[[nm]]<-az$State
      }
      hfun<-function(t,state,...) .inapest_cbr_host_rate_model(t,state,host_fun,pars,attack_h,
        if(couples_pathogen)prates else NULL,if(couples_pathogen)pathogen else NULL,
        if(couples_pathogen)pathogen_engine else NULL,if(couples_pathogen)pathogen_context else NULL,
        step_index,nts,runtime)
      hz<-INApestStochasticCompartmentStep(Hstate,tt,h,hfun)
      Hstate<-hz$State
      for(nm in names(agents)) {
        cause<-paste0("bio:",nm); killed<-if(cause%in%colnames(hz$Exits))as.integer(hz$Exits[,cause]) else integer(biocontrol_context$n_nodes)
        attacks_node[[nm]]<-attacks_node[[nm]]+killed
        rec<-killed*as.integer(agents[[nm]]$OffspringPerAttack)
        next_agents[[nm]][,agents[[nm]]$RecruitStage]<-next_agents[[nm]][,agents[[nm]]$RecruitStage]+rec
      }
      Astate<-next_agents;tt<-tt+h
    }

    # Established companion movement remains the final biocontrol event.
    for(nm in names(agents)) {
      a<-agents[[nm]]; M<-.inabc_resolve_movement(a$Movement,step_index,biocontrol_context,a)
      Astate[[nm]]<-.inabc_move_agent(Astate[[nm]],M,a$MovementStages);storage.mode(Astate[[nm]])<-"integer"
    }
    n_after<-as.integer(rowSums(Hstate))
    Kvec<-as.numeric(nodeK);if(length(Kvec)==1L)Kvec<-rep(Kvec,length(n_after))
    over<-which(n_after>Kvec+cap_tol*pmax(1,abs(Kvec)))
    if(length(over))stop("Coupled RK host state exceeded Meta nodeK before dispersal at node(s): ",paste(over,collapse=", "))
    n_out<-disperse(sddprob=sddprob,nodepropaguleproduction=nodepropaguleproduction,nodeenvestabprob=nodeenvestabprob,
      n=n_after,lddprob=lddprob,lddrate=lddrate,k_is_0=k_is_0,nodeK=nodeK,nodepropaguleestablishment=nodepropaguleestablishment,
      nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,managing=managing,maxinteger=maxinteger)
    if(couples_pathogen) Hout<-pathogen_engine$Reconcile(Hstate,n_out,pathogen_context) else Hout<-NULL
    impact_att<-lapply(names(agents),function(nm)matrix(attacks_node[[nm]],ncol=1L,dimnames=list(NULL,"all")));names(impact_att)<-names(agents)
    rec_node<-lapply(names(agents),function(nm)as.integer(attacks_node[[nm]]*agents[[nm]]$OffspringPerAttack));names(rec_node)<-names(agents)
    rec_total<-vapply(rec_node,sum,integer(1))
    impact<-list(AttacksByAgent=impact_att,RecruitsByAgent=rec_total,RecruitsByAgentNode=rec_node)
    out<-list(N=as.integer(n_out),BiocontrolState=Astate,BiocontrolImpact=impact,
      CoupledRK=list(PreDispersalHost=n_after,Attacks=attacks_node,Timestep=step_index,Ntimesteps=nts,NSubsteps=nsub,InternalStep=h))
    if(couples_pathogen) out$PathogenState<-Hout
    out
  }
  attr(f,"INApestMetaRK4CoupledBiologyLocalDynamics")<-list(Version="0.5",StateMode="integer-coupled",
    TimestepLength=dt,RKMaxStep=hmax,CouplesPathogen=couples_pathogen,
    Ordering=c("parent survival/management mortality","scheduled biocontrol release","joint stochastic RK local biology",
      "agent movement","Meta host dispersal/establishment","external invasion","surveillance"))
  attr(f,"INApestCouplesBiocontrol")<-TRUE
  attr(f,"INApestCouplesPathogen")<-couples_pathogen
  # Required so biocontrol stages continue to age/move even if hosts reach zero.
  attr(f,"INApestRunWhenEmpty")<-TRUE
  class(f)<-c("INApestMetaRK4CoupledBiologyLocalDynamics","function")
  f
}
