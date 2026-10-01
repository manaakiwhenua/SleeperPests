###############################################################################
# INApestVertebratePointContinuousBiology.R
# Continuous pathogen / biocontrol coupling for explicit Vertebrate Point hosts
# v0.9 -- 1 October 2026
#
# Point identity, birth, demographic stage progression, spatial movement,
# vertebrate control and generic Interaction remain parent-engine processes.
# This companion integrates only within-timestep pathogen state and node-based
# biocontrol demography/attack, using the stochastic RK compartment bridge.
###############################################################################

.ivpcrk_stop <- function(...) stop(..., call.=FALSE)

.ivpcrk_resolve_point <- function(x, points, timestep, perm, name, nonnegative=TRUE) {
  n <- nrow(points)
  if (!n) return(numeric())
  if (is.function(x)) {
    fm <- names(formals(x))
    a <- list(points=points, timestep=timestep, perm=perm, n_points=n)
    if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a),fm)]
    x <- do.call(x,a)
  }
  x <- as.numeric(x)
  if (length(x)==1L) x <- rep(x,n)
  if (length(x)!=n || any(!is.finite(x))) .ivpcrk_stop(name," must resolve to one finite value per point")
  if (nonnegative && any(x<0)) .ivpcrk_stop(name," must resolve to non-negative values")
  x
}

.ivpcrk_states <- function(spec) spec$Pathogen$States

INApestContinuousPointPathogen <- function(
    Pathogen,
    PathogenRates=INApestContinuousPathogenRates(),
    ContactRadius=Inf,
    ContactKernel=NULL,
    StateField="pathogen_state",
    DefaultState="S") {
  if (!inherits(Pathogen,"INApestPathogen") || identical(Pathogen$Model,"Binary"))
    .ivpcrk_stop("Continuous point pathogen requires a non-binary INApestPathogen object")
  if (!inherits(PathogenRates,"INApestContinuousPathogenRates"))
    .ivpcrk_stop("PathogenRates must come from INApestContinuousPathogenRates()")
  if (!is.numeric(ContactRadius) || length(ContactRadius)!=1L || is.na(ContactRadius) || ContactRadius<0)
    .ivpcrk_stop("ContactRadius must be one non-negative value or Inf")
  if (!is.null(ContactKernel) && !is.function(ContactKernel))
    .ivpcrk_stop("ContactKernel must be NULL or a function")
  if (!is.character(StateField) || length(StateField)!=1L || !nzchar(StateField))
    .ivpcrk_stop("StateField must be one non-empty name")
  if (!(DefaultState %in% Pathogen$States)) .ivpcrk_stop("DefaultState must be a state in Pathogen$States")
  structure(list(Pathogen=Pathogen,PathogenRates=PathogenRates,ContactRadius=as.numeric(ContactRadius),
                 ContactKernel=ContactKernel,StateField=StateField,DefaultState=DefaultState),
            class=c("INApestContinuousPointPathogen","list"))
}

INApestVertebratePointContinuousBiology <- function(
    Pathogen=NULL,
    BiocontrolRates=NULL,
    TimestepLength=1,
    RKMaxStep=0.025,
    StartTime=0) {
  if (!is.null(Pathogen) && !inherits(Pathogen,"INApestContinuousPointPathogen"))
    .ivpcrk_stop("Pathogen must be NULL or INApestContinuousPointPathogen()")
  if (!is.null(BiocontrolRates) && !inherits(BiocontrolRates,"INApestContinuousBiocontrolRates"))
    .ivpcrk_stop("BiocontrolRates must be NULL or INApestContinuousBiocontrolRates()")
  dt <- as.numeric(TimestepLength)[1L]; h <- as.numeric(RKMaxStep)[1L]; t0 <- as.numeric(StartTime)[1L]
  if (!is.finite(dt) || dt<=0 || !is.finite(h) || h<=0 || !is.finite(t0))
    .ivpcrk_stop("TimestepLength/RKMaxStep must be positive and StartTime finite")
  structure(list(Pathogen=Pathogen,BiocontrolRates=BiocontrolRates,TimestepLength=dt,RKMaxStep=h,StartTime=t0),
            class=c("INApestVertebratePointContinuousBiology","list"))
}

.ivpcrk_ensure_pathogen_state <- function(points, spec, perm=1L, seed_initial=FALSE) {
  if (is.null(spec) || !nrow(points)) return(points)
  fld <- spec$StateField; states <- .ivpcrk_states(spec)
  explicit <- fld %in% names(points) && any(!is.na(points[[fld]]) & nzchar(as.character(points[[fld]])))
  if (!(fld %in% names(points))) points[[fld]] <- spec$DefaultState
  z <- as.character(points[[fld]])
  z[is.na(z)|!nzchar(z)] <- spec$DefaultState
  bad <- setdiff(unique(z),states); if(length(bad)) .ivpcrk_stop("Unknown point pathogen state(s): ",paste(bad,collapse=", "))
  points[[fld]] <- z
  if (!seed_initial || explicit) return(points)
  p <- spec$Pathogen
  get_count <- function(x,name){
    if(is.function(x)||length(x)!=1L||!is.finite(x)||x<0||x!=floor(x))
      .ivpcrk_stop(name," must be one non-negative whole-number count for automatic point seeding")
    as.integer(x)
  }
  ni<-get_count(p$InitialInfected,"InitialInfected")
  ne<-if("E"%in%states)get_count(p$InitialExposed,"InitialExposed")else 0L
  nr<-if("R"%in%states)get_count(p$InitialRecovered,"InitialRecovered")else 0L
  if(ni+ne+nr>nrow(points)) .ivpcrk_stop("Initial pathogen counts exceed InitialPoints")
  points[[fld]] <- spec$DefaultState
  available <- seq_len(nrow(points))
  assign_state <- function(k,s){if(!k)return();pick<-sample(available,k,FALSE);points[[fld]][pick]<<-s;available<<-setdiff(available,pick)}
  assign_state(ne,"E");assign_state(ni,"I");assign_state(nr,"R")
  points
}

INApestVertebratePointContinuousInitialize <- function(points, ContinuousBiology, perm=1L) {
  if (!inherits(ContinuousBiology,"INApestVertebratePointContinuousBiology"))
    .ivpcrk_stop("ContinuousBiology must be INApestVertebratePointContinuousBiology()")
  .ivpcrk_ensure_pathogen_state(points,ContinuousBiology$Pathogen,perm,seed_initial=TRUE)
}

.ivpcrk_contact_matrix <- function(points, pspec, timestep, perm) {
  n<-nrow(points); W<-matrix(0,n,n)
  if(n<2L)return(W)
  dx<-outer(points$x,points$x,"-");dy<-outer(points$y,points$y,"-");d<-sqrt(dx^2+dy^2)
  W[d<=pspec$ContactRadius] <- 1; diag(W)<-0
  if(!is.null(pspec$ContactKernel)) {
    ij<-which(W>0,arr.ind=TRUE)
    if(nrow(ij)) {
      dist<-d[ij]
      src<-points[ij[,1],,drop=FALSE];tgt<-points[ij[,2],,drop=FALSE]
      fm<-names(formals(pspec$ContactKernel));a<-list(distance=dist,source=src,target=tgt,timestep=timestep,perm=perm)
      if(!is.null(fm)&&!("..."%in%fm))a<-a[intersect(names(a),fm)]
      km<-do.call(pspec$ContactKernel,a);if(length(km)==1L)km<-rep(km,nrow(ij))
      if(length(km)!=nrow(ij)||any(!is.finite(km))||any(km<0)) .ivpcrk_stop("ContactKernel must return finite non-negative multipliers")
      W[]<-0;W[ij]<-km
    }
  }
  W
}

.ivpcrk_pathogen_rate_arrays <- function(points, pspec, timestep, perm) {
  pr<-pspec$PathogenRates;p<-pspec$Pathogen
  list(
    beta=.ivpcrk_resolve_point(p$Beta,points,timestep,perm,"Beta"),
    density=.ivpcrk_resolve_point(p$DensityScale,points,timestep,perm,"DensityScale"),
    rec=.ivpcrk_resolve_point(pr$RecoveryRate,points,timestep,perm,"RecoveryRate"),
    prog=.ivpcrk_resolve_point(pr$ProgressionRate,points,timestep,perm,"ProgressionRate"),
    mort=.ivpcrk_resolve_point(pr$PathogenMortalityRate,points,timestep,perm,"PathogenMortalityRate"),
    wan=.ivpcrk_resolve_point(pr$ImmunityLossRate,points,timestep,perm,"ImmunityLossRate")
  )
}

.ivpcrk_host_rate_function <- function(points, pspec, rates, W, attack_hazards) {
  cp<-!is.null(pspec); causes<-c(if(cp)"pathogen" else character(),if(!is.null(attack_hazards))paste0("bio:",colnames(attack_hazards))else character())
  function(t,State,...) {
    State<-as.matrix(State);n<-nrow(State);comps<-colnames(State);k<-ncol(State)
    gain<-matrix(0,n,k,dimnames=dimnames(State));trans<-array(0,c(n,k,k),dimnames=list(rownames(State),comps,comps))
    exits<-array(0,c(n,k,length(causes)),dimnames=list(rownames(State),comps,causes))
    if(cp) {
      p<-pspec$Pathogen
      if(!all(c("S","I")%in%comps)) .ivpcrk_stop("Continuous point pathogen state must contain S and I")
      I<-State[,"I"]; pressure<-as.numeric(crossprod(I,W))
      lambda<-rates$beta*pressure
      if(p$Transmission=="density") lambda<-lambda/rates$density
      lambda<-pmax(0,lambda)
      if("E"%in%comps){trans[,"S","E"]<-lambda;trans[,"E","I"]<-rates$prog}else{if(any(rates$prog>0)) .ivpcrk_stop("ProgressionRate must be zero without E");trans[,"S","I"]<-lambda}
      if(p$Model=="SIS"){trans[,"I","S"]<-rates$rec;if(any(rates$wan>0)) .ivpcrk_stop("ImmunityLossRate must be zero for SIS")}else{trans[,"I","R"]<-rates$rec;trans[,"R","S"]<-rates$wan}
      exits[,"I","pathogen"]<-rates$mort
    }
    if(!is.null(attack_hazards)&&nrow(attack_hazards))for(j in seq_len(ncol(attack_hazards)))
      exits[,,paste0("bio:",colnames(attack_hazards)[j])]<-matrix(rep(attack_hazards[,j],k),n,k)
    list(GainRate=gain,TransitionHazards=trans,ExitHazards=exits)
  }
}

.ivpcrk_agent_exposure <- function(State,Time,Step,RateFunction,AttackStages) {
  AttackStages<-as.integer(AttackStages);out<-rep(0,nrow(State))
  for(s in AttackStages)out<-out+.inapest_tm_cbr_attacker_exposure(State,Time,Step,RateFunction,s)
  out
}

.ivpcrk_add_recruits <- function(state,node,killed,offspring,recruit_stage) {
  if(!length(killed)||!any(killed>0)||offspring<=0L)return(state)
  rec<-integer(nrow(state));idx<-which(killed>0 & !is.na(node))
  if(length(idx))for(i in idx)rec[node[i]]<-rec[node[i]]+as.integer(killed[i])*offspring
  state[,recruit_stage]<-state[,recruit_stage]+rec;state
}

.ivpcrk_detection <- function(points, pspec, timestep, perm) {
  if(is.null(pspec)||!nrow(points))return(logical(nrow(points)))
  p<-pspec$Pathogen;fld<-pspec$StateField
  pd<-.ivpcrk_resolve_point(p$DetectionProb,points,timestep,perm,"DetectionProb",nonnegative=FALSE)
  if(any(pd<0|pd>1)) .ivpcrk_stop("DetectionProb must resolve to [0,1]")
  points[[fld]]=="I" & (rbinom(nrow(points),1L,pd)==1L)
}

INApestVertebratePointContinuousBiologyStep <- function(
    points, PointBiocontrolState=NULL, Biocontrol=NULL, PointSupport=NULL,
    ContinuousBiology, timestep, perm=1L, HostStages=NULL) {
  if(!inherits(ContinuousBiology,"INApestVertebratePointContinuousBiology")) .ivpcrk_stop("Invalid ContinuousBiology")
  if(!is.data.frame(points)||!all(c("id","x","y","stage")%in%names(points))) .ivpcrk_stop("Continuous point biology requires id/x/y/stage")
  pspec<-ContinuousBiology$Pathogen;br<-ContinuousBiology$BiocontrolRates;cp<-!is.null(pspec);cb<-!is.null(br)
  dt<-ContinuousBiology$TimestepLength;hmax<-ContinuousBiology$RKMaxStep;tt<-ContinuousBiology$StartTime+(as.integer(timestep)-1L)*dt
  if(cp) {
    if(any(.ivpcrk_resolve_point(pspec$Pathogen$IntroductionProb,points,timestep,perm,"IntroductionProb",FALSE)!=0))
      .ivpcrk_stop("Continuous Vertebrate Point v0.9 requires Pathogen$IntroductionProb = 0")
    points<-.ivpcrk_ensure_pathogen_state(points,pspec,perm,FALSE)
    # New offspring enter the susceptible/default state; external arrivals retain any explicit supplied state.
    if(all(c("birth_timestep","parent_id")%in%names(points))){nb<-points$birth_timestep==timestep & !is.na(points$parent_id);nb[is.na(nb)]<-FALSE;if(any(nb))points[[pspec$StateField]][nb]<-pspec$DefaultState}
  }
  agents<-NULL;Astate<-NULL;context<-NULL
  if(cb) {
    if(is.null(Biocontrol)||!inherits(Biocontrol,"INApestBiocontrol")) .ivpcrk_stop("Biocontrol object required when BiocontrolRates are active")
    if(is.null(PointBiocontrolState)||!inherits(PointBiocontrolState,"INApestPointBiocontrolState")) .ivpcrk_stop("PointBiocontrolState is required")
    if(is.null(PointSupport)||!inherits(PointSupport,"INApestBiocontrolPointSupport")) .ivpcrk_stop("PointSupport is required")
    if(abs(.ibp_dt(Biocontrol)-dt)>1e-12) .ivpcrk_stop("Biocontrol TimestepLength must equal ContinuousBiology TimestepLength")
    agents<-.ibp_agents(Biocontrol);if(!identical(names(br),names(PointBiocontrolState$State))) .ivpcrk_stop("BiocontrolRates names must match point biocontrol state")
    Astate<-PointBiocontrolState$State;context<-PointBiocontrolState$Context;context$timestep<-timestep;context$perm<-perm
    for(a in seq_along(agents)){rel<-.ibp_resolve_release(agents[[a]],context,timestep,perm);Astate[[a]]<-Astate[[a]]+rel}
  }
  n<-nrow(points);node<-if(cb)INApestBiocontrolPointNode(points$x,points$y,PointSupport)else rep(NA_integer_,n)
  if(cp){states<-.ivpcrk_states(pspec);Hstate<-matrix(0L,n,length(states),dimnames=list(as.character(points$id),states));Hstate[cbind(seq_len(n),match(points[[pspec$StateField]],states))]<-1L;W<-.ivpcrk_contact_matrix(points,pspec,timestep,perm);rates<-.ivpcrk_pathogen_rate_arrays(points,pspec,timestep,perm)}else{Hstate<-matrix(1L,n,1L,dimnames=list(as.character(points$id),"H"));W<-matrix(0,n,n);rates<-NULL}
  exits_total<-matrix(0L,n,(if(cp)1L else 0L)+(if(cb)length(agents)else 0L));cn<-c(if(cp)"pathogen" else character(),if(cb)paste0("bio:",names(Astate))else character());colnames(exits_total)<-cn
  nsub<-max(1L,as.integer(ceiling(dt/hmax-1e-14)));h<-dt/nsub
  for(ss in seq_len(nsub)) {
    ah<-if(cb)matrix(0,n,length(agents),dimnames=list(NULL,names(Astate)))else NULL;nextA<-Astate
    if(cb)for(a in seq_along(agents)){
      nm<-names(Astate)[a];ag<-agents[[a]];rfA<-.inapest_tm_cbr_agent_rate_function(br[[nm]],ag,timestep,context)
      expo<-.ivpcrk_agent_exposure(Astate[[nm]],tt,h,rfA,.ibp_attack_stage_index(ag));ar<-.ibp_resolve_rate(ag,timestep,perm,context)
      eligible<-.ibp_target_mask(points,ag,"PointTransitionMatrix",HostStages)&!is.na(node)&rowSums(Hstate)>0
      ii<-which(eligible);if(length(ii))ah[ii,nm]<-ar[node[ii]]*(expo[node[ii]]/h)
      nextA[[nm]]<-INApestStochasticCompartmentStep(Astate[[nm]],tt,h,rfA)$State
    }
    if(n>0L){rfH<-.ivpcrk_host_rate_function(points,pspec,rates,W,ah);hz<-INApestStochasticCompartmentStep(Hstate,tt,h,rfH);Hstate<-hz$State;if(length(cn))exits_total<-exits_total+hz$Exits
      if(cb)for(nm in names(Astate)){cause<-paste0("bio:",nm);killed<-if(cause%in%colnames(hz$Exits))hz$Exits[,cause]else integer(n);ag<-agents[[match(nm,names(Astate))]];nextA[[nm]]<-.ivpcrk_add_recruits(nextA[[nm]],node,killed,.ibp_offspring(ag),.ibp_recruit_stage(ag))}}
    Astate<-nextA;tt<-tt+h
  }
  # Agent movement is the final biocontrol operation, after all continuous local attacks/recruitment.
  if(cb)for(a in seq_along(agents)){nm<-names(Astate)[a];M<-.ibp_resolve_movement(agents[[a]],timestep,perm,context);Astate[[nm]]<-.ibp_move_state(Astate[[nm]],M)}
  dead<-rowSums(Hstate)==0L;pathdeath<-if(cp&&"pathogen"%in%colnames(exits_total))exits_total[,"pathogen"]>0 else rep(FALSE,n)
  attack_events<-data.frame();if(cb&&n){parts<-list();kk<-0L;for(nm in names(Astate)){cause<-paste0("bio:",nm);idx<-which(exits_total[,cause]>0);if(length(idx)){kk<-kk+1L;ae<-points[idx,intersect(c("id","parent_id","x","y","stage"),names(points)),drop=FALSE];names(ae)[names(ae)=="id"]<-"point_id";if(!"parent_id"%in%names(ae))ae$parent_id<-NA_integer_;ae$perm<-perm;ae$timestep<-timestep;ae$node<-node[idx];ae$agent<-nm;ae$detail<-"continuous_biocontrol_mortality";parts[[kk]]<-ae[,c("perm","timestep","point_id","parent_id","x","y","stage","node","agent","detail"),drop=FALSE]}};if(length(parts))attack_events<-do.call(rbind,parts)}
  pathogen_events<-data.frame();if(cp&&n){fld<-pspec$StateField;from<-as.character(points[[fld]]);to<-rep(NA_character_,n);surv<-which(!dead);if(length(surv))to[surv]<-colnames(Hstate)[max.col(Hstate[surv,,drop=FALSE],ties.method="first")];pe_idx<-which(dead|from!=to);if(length(pe_idx))pathogen_events<-data.frame(perm=perm,timestep=timestep,point_id=points$id[pe_idx],from_state=from[pe_idx],to_state=to[pe_idx],pathogen_death=pathdeath[pe_idx],diagnostic="endpoint_change_not_exact_transition_count",stringsAsFactors=FALSE);if(length(surv))points[[fld]][surv]<-to[surv]}
  if(any(dead))points<-points[!dead,,drop=FALSE]
  detected<-if(cp).ivpcrk_detection(points,pspec,timestep,perm)else logical(nrow(points));if(cp&&any(detected)){p<-pspec$Pathogen;if(isTRUE(p$DetectionTriggersInfo)){if("have_info"%in%names(points))points$have_info[detected]<-TRUE;if("last_known_timestep"%in%names(points))points$last_known_timestep[detected]<-timestep};de<-data.frame(perm=perm,timestep=timestep,point_id=points$id[detected],from_state=points[[pspec$StateField]][detected],to_state=points[[pspec$StateField]][detected],pathogen_death=FALSE,diagnostic="pathogen_detection",stringsAsFactors=FALSE);pathogen_events<-if(nrow(pathogen_events))rbind(pathogen_events,de)else de}
  if(cb){PointBiocontrolState$State<-Astate;hist<-INApestPointBiocontrolStateHistory(PointBiocontrolState,perm,timestep)}else hist<-data.frame()
  list(points=points,state=PointBiocontrolState,attacks=attack_events,history=hist,pathogen_events=pathogen_events,
       n_biocontrol_deaths=if(nrow(attack_events))nrow(attack_events)else 0L,n_pathogen_deaths=sum(pathdeath),
       diagnostics=list(NSubsteps=nsub,InternalStep=h,PointIdentity="surviving IDs preserved",HostDynamics="parent point engine",AgentMovement=if(cb)"last" else "inactive"))
}
