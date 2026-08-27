###############################################################################
### INApest point-pathogen analytical methods
###############################################################################

.ina_point_resolve <- function(x,points,timestep=1L,perm=1L,name="parameter"){
  n<-nrow(points);if(is.function(x)){
    fm<-names(formals(x));a<-list(points=points,timestep=timestep,perm=perm)
    if(!is.null(fm)&&!"..."%in%fm)a<-a[intersect(names(a),fm)]
    z<-do.call(x,a)
  } else z<-if(length(x)==1L)rep(as.numeric(x),n)else as.numeric(x)
  if(length(z)==1L)z<-rep(z,n)
  if(length(z)!=n||any(!is.finite(z)))stop(name," must resolve to one finite value per point")
  z
}

INApestPointPathogenEdgeMatrix <- function(
    Points,Pathogen,ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,
    timestep=1L,perm=1L,IncludeDiagonal=FALSE){
  if(!is.data.frame(Points)||!all(c("x","y")%in%names(Points)))stop("Points must contain x and y")
  if(inherits(Pathogen,"INApestPointPathogenInteraction"))Pathogen<-Pathogen$Pathogen
  if(!inherits(Pathogen,"INApestPathogen"))stop("Pathogen must be created by INApestPathogen() or supplied through INApestPointPathogenInteraction()")
  n<-nrow(Points);Q<-matrix(0,n,n)
  if(n<1L)return(Q)
  if(n<2L&&!isTRUE(IncludeDiagonal))return(Q)
  beta<-.ina_point_resolve(Pathogen$Beta,Points,timestep,perm,"Beta")
  if(any(beta<0))stop("Beta must be non-negative")
  ds<-.ina_point_resolve(Pathogen$DensityScale,Points,timestep,perm,"DensityScale")
  if(any(ds<=0))stop("DensityScale must be positive")
  pairs<-expand.grid(i=seq_len(n),j=seq_len(n),KEEP.OUT.ATTRS=FALSE)
  if(!isTRUE(IncludeDiagonal))pairs<-pairs[pairs$i!=pairs$j,,drop=FALSE]
  dx<-Points$x[pairs$i]-Points$x[pairs$j];dy<-Points$y[pairs$i]-Points$y[pairs$j];d<-sqrt(dx^2+dy^2)
  keep<-d<=ContactRadius;pairs<-pairs[keep,,drop=FALSE];d<-d[keep]
  if(!nrow(pairs))return(Q)
  if(is.function(ContactProb)){
    fm<-names(formals(ContactProb));a<-list(distance=d,source=Points[pairs$i,,drop=FALSE],target=Points[pairs$j,,drop=FALSE],timestep=timestep,perm=perm)
    if(!is.null(fm)&&!"..."%in%fm)a<-a[intersect(names(a),fm)];cp<-do.call(ContactProb,a)
  } else cp<-ContactProb
  if(length(cp)==1L)cp<-rep(cp,nrow(pairs));if(length(cp)!=nrow(pairs)||any(!is.finite(cp)))stop("ContactProb must resolve per candidate pair")
  if(!is.null(ContactKernel)){
    if(!is.function(ContactKernel))stop("ContactKernel must be NULL or a function")
    fm<-names(formals(ContactKernel));a<-list(distance=d,source=Points[pairs$i,,drop=FALSE],target=Points[pairs$j,,drop=FALSE],timestep=timestep,perm=perm)
    if(!is.null(fm)&&!"..."%in%fm)a<-a[intersect(names(a),fm)];km<-do.call(ContactKernel,a)
    if(length(km)==1L)km<-rep(km,nrow(pairs));if(length(km)!=nrow(pairs)||any(!is.finite(km))||any(km<0))stop("ContactKernel must return finite non-negative multipliers")
    cp<-cp*km
  }
  cp<-pmin(1,pmax(0,cp))
  if(Pathogen$Transmission=="frequency")edge<-cp*pmin(1,pmax(0,beta[pairs$j]))
  else edge<-cp*(-expm1(-beta[pairs$j]/ds[pairs$j]))
  Q[cbind(pairs$j,pairs$i)]<-pmin(1,pmax(0,edge)) # row target, col source
  dimnames(Q)<-list(target=if("id"%in%names(Points))Points$id else seq_len(n),source=if("id"%in%names(Points))Points$id else seq_len(n))
  Q
}

INApestPointPathogenGrowthOperator <- function(
    Model=c("SIS","SIR","SEIR"),Points,Pathogen,
    ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,timestep=1L,perm=1L){
  Model<-match.arg(Model);if(Model!=Pathogen$Model)warning("Model differs from Pathogen$Model; using requested analytical compartment structure")
  n<-nrow(Points);Q<-INApestPointPathogenEdgeMatrix(Points,Pathogen,ContactRadius,ContactKernel,ContactProb,timestep,perm)
  rec<-.ina_point_resolve(Pathogen$RecoveryProb,Points,timestep,perm,"RecoveryProb");mort<-.ina_point_resolve(Pathogen$PathogenMortalityProb,Points,timestep,perm,"PathogenMortalityProb")
  prog<-.ina_point_resolve(Pathogen$ProgressionProb,Points,timestep,perm,"ProgressionProb")
  if(any(rec<0|rec>1)||any(mort<0|mort>1)||any(rec+mort>1+1e-12)||any(prog<0|prog>1))stop("Invalid pathogen probabilities")
  stayI<-diag(1-rec-mort,n)
  if(Model%in%c("SIS","SIR"))G<-stayI+Q else G<-rbind(cbind(diag(1-prog,n),Q),cbind(diag(prog,n),stayI))
  lam<-max(Mod(eigen(G,only.values=TRUE)$values))
  list(Model="INApestMetaPoint",PathogenModel=Model,Operator=G,TransmissionEdgeMatrix=Q,
       IntrinsicRarePathogenMultiplier=lam,Classification=if(lam>1+1e-12)"growing"else if(lam<1-1e-12)"declining"else"threshold",
       Diagnostics=c(
         "The realised-contact layer is integrated out exactly: each source-target pair contributes an independent one-source transmission probability.",
         "Newly infected points do not recover/progress in the same timestep, matching the point pathogen engine.",
         if(Model=="SEIR")"SEIR preserves the one-timestep latent-state delay."else"SIS/SIR active-I linearisation."))
}

INApestPointPathogenNextGeneration <- function(GrowthObject){
  x<-GrowthObject;G<-x$Operator;Q<-x$TransmissionEdgeMatrix;n<-nrow(Q);model<-x$PathogenModel;Z<-matrix(0,n,n)
  if(model%in%c("SIS","SIR")){T<-G-Q;F<-Q}else{F<-rbind(cbind(Z,Q),cbind(Z,Z));T<-G-F}
  rhoT<-max(Mod(eigen(T,only.values=TRUE)$values))
  if(rhoT>=1-1e-12){K<-matrix(NA,nrow(T),ncol(T));R0<-Inf}else{K<-F%*%solve(diag(nrow(T))-T);R0<-max(Mod(eigen(K,only.values=TRUE)$values))}
  list(TransitionWithoutNewInfection=T,NewInfectionOperator=F,NextGenerationOperator=K,R0=R0,
       PerTimestepLambda=x$IntrinsicRarePathogenMultiplier,
       ThresholdAgreement=if(is.finite(R0))sign(R0-1)==sign(x$IntrinsicRarePathogenMultiplier-1)else x$IntrinsicRarePathogenMultiplier>=1,
       Diagnostics="Point-model R0 is the lifetime secondary-infection spectral radius after integrating over repeated infectious persistence at the same point.")
}

.ina_point_state_table <- function(Model,n){
  vals<-switch(Model,SIS=c("S","I"),SIR=c("S","I","R"),SEIR=c("S","E","I","R"));m<-length(vals)
  g<-expand.grid(rep(list(seq_len(m)),n),KEEP.OUT.ATTRS=FALSE);z<-matrix(vals[as.matrix(g)],nrow=nrow(g),ncol=n);colnames(z)<-paste0("p",seq_len(n));z
}

.ina_point_state_key <- function(z)paste(z,collapse="|")

.ina_point_exact_operator <- function(Model,Points,Pathogen,Q,MaxStates=100000L,timestep=1L,perm=1L){
  n<-nrow(Points);st<-.ina_point_state_table(Model,n);if(nrow(st)>MaxStates)stop("Exact point pathogen state space exceeds MaxStates")
  keys<-apply(st,1,.ina_point_state_key);idx<-setNames(seq_len(nrow(st)),keys);T<-matrix(0,nrow(st),nrow(st),dimnames=list(keys,keys))
  rec<-.ina_point_resolve(Pathogen$RecoveryProb,Points,timestep,perm,"RecoveryProb");mort<-.ina_point_resolve(Pathogen$PathogenMortalityProb,Points,timestep,perm,"PathogenMortalityProb")
  if(any(mort!=0))stop("Exact fixed-point operator requires PathogenMortalityProb = 0 because pathogen death removes host points")
  prog<-.ina_point_resolve(Pathogen$ProgressionProb,Points,timestep,perm,"ProgressionProb");wan<-.ina_point_resolve(Pathogen$ImmunityLossProb,Points,timestep,perm,"ImmunityLossProb");ip<-.ina_point_resolve(Pathogen$IntroductionProb,Points,timestep,perm,"IntroductionProb")
  if(any(rec<0|rec>1)||any(prog<0|prog>1)||any(wan<0|wan>1)||any(ip<0|ip>1))stop("Invalid pathogen probabilities")
  # Enumerate independent per-point end-state distributions conditional on the start state.
  for(r in seq_len(nrow(st))){s0<-st[r,];inf<-which(s0=="I");outs<-vector("list",n)
    for(j in seq_len(n)){
      if(s0[j]=="S"){
        pinf<-if(length(inf))1-prod(1-Q[j,inf])else 0
        # Introduction applies only if not newly infected.
        if(Model=="SEIR"){
          pE<-pinf+(1-pinf)*ip[j];outs[[j]]<-c(S=1-pE,E=pE)
        } else {pI<-pinf+(1-pinf)*ip[j];outs[[j]]<-c(S=1-pI,I=pI)}
      } else if(s0[j]=="I"){
        if(Model=="SIS")outs[[j]]<-c(S=rec[j],I=1-rec[j])
        else outs[[j]]<-c(I=1-rec[j],R=rec[j])
      } else if(s0[j]=="E")outs[[j]]<-c(E=1-prog[j],I=prog[j])
      else if(s0[j]=="R")outs[[j]]<-c(S=wan[j],R=1-wan[j])
    }
    # Cartesian product of point outcomes.
    namesv<-lapply(outs,names);grid<-expand.grid(lapply(namesv,seq_along),KEEP.OUT.ATTRS=FALSE)
    for(k in seq_len(nrow(grid))){dest<-character(n);pr<-1
      for(j in seq_len(n)){ii<-grid[k,j];dest[j]<-namesv[[j]][ii];pr<-pr*outs[[j]][ii]}
      if(pr>0)T[r,idx[[.ina_point_state_key(dest)]]]<-T[r,idx[[.ina_point_state_key(dest)]]]+pr
    }
  }
  if(max(abs(rowSums(T)-1))>1e-11)stop("Internal exact point operator error")
  attr(T,"States")<-st;T
}

.ina_point_hitting_probability <- function(T,target,avoid=rep(FALSE,nrow(T)),tol=1e-13,maxiter=100000L){
  target<-as.logical(target);avoid<-as.logical(avoid);if(length(target)!=nrow(T)||length(avoid)!=nrow(T))stop("hitting-set length mismatch")
  unknown<-!(target|avoid);u<-as.numeric(target)
  if(any(unknown)){
    Q<-T[unknown,unknown,drop=FALSE];r<-rowSums(T[unknown,target,drop=FALSE]);z<-rep(0,sum(unknown))
    for(k in seq_len(maxiter)){zn<-as.numeric(r+Q%*%z);if(max(abs(zn-z))<tol){z<-zn;break};z<-zn}
    u[unknown]<-pmin(1,pmax(0,z))
  }
  u
}

.ina_point_static_absorption <- function(T,States,model,outside,initial_index){
  active_state<-apply(States,1,function(z)any(z=="I"|(model=="SEIR"&z=="E")))
  extinct<-!active_state
  escaped<-if(any(outside))apply(States[,outside,drop=FALSE],1,function(z)any(z=="I"|(model=="SEIR"&z=="E")))else rep(FALSE,nrow(States))
  pe<-if(any(escaped)).ina_point_hitting_probability(T,escaped) else rep(0,nrow(T))
  # With no exogenous reintroduction, extinction is an absorbing pathogen event;
  # treating escape as the competing absorbing event gives the exact ordering.
  pext_first<-.ina_point_hitting_probability(T,extinct,escaped)
  pesc_first<-if(any(escaped)).ina_point_hitting_probability(T,escaped,extinct) else rep(0,nrow(T))
  transient<-!(extinct|escaped);etime<-NA_real_
  if(transient[initial_index]){
    Q<-T[transient,transient,drop=FALSE];rho<-if(nrow(Q))max(Mod(eigen(Q,only.values=TRUE)$values))else 0
    if(rho<1-1e-12){tt<-solve(diag(nrow(Q))-Q,rep(1,nrow(Q)));etime<-tt[which(which(transient)==initial_index)]}
  } else etime<-0
  list(EventualEscapeProbability=pe[initial_index],EscapeBeforeExtinctionProbability=pesc_first[initial_index],
       ExtinctionBeforeEscapeProbability=pext_first[initial_index],ExpectedTimestepsToEscapeOrExtinction=etime,
       PerStateEventualEscape=pe)
}

INApestPointPathogenExactFixedGeometry <- function(
    Ntimesteps,Points,Pathogen,ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,
    InitialState=NULL,Outside=NULL,MaxStates=100000L,ReturnOperator=FALSE){
  cr_missing<-missing(ContactRadius);ck_missing<-missing(ContactKernel);cp_missing<-missing(ContactProb)
  interaction<-inherits(Pathogen,"INApestPointPathogenInteraction")
  if(interaction&&cr_missing)ContactRadius<-.ina_point_interaction_contact(Pathogen,"ContactRadius",Inf)
  if(interaction&&ck_missing)ContactKernel<-.ina_point_interaction_contact(Pathogen,"ContactKernel",NULL)
  if(interaction&&cp_missing)ContactProb<-.ina_point_interaction_contact(Pathogen,"ContactProb",1)
  Pathogen<-.ina_point_pathogen_spec(Pathogen)
  model<-Pathogen$Model;if(!model%in%c("SIS","SIR","SEIR"))stop("Exact point pathogen analysis supports SIS/SIR/SEIR")
  Ntimesteps<-as.integer(Ntimesteps);if(length(Ntimesteps)!=1L||is.na(Ntimesteps)||Ntimesteps<1L)stop("Ntimesteps must be a positive integer")
  n<-nrow(Points);st<-.ina_point_state_table(model,n);if(nrow(st)>MaxStates)stop("Exact point pathogen state space exceeds MaxStates")
  keys<-apply(st,1,.ina_point_state_key)
  if(is.null(InitialState)){
    z<-rep("S",n);ni<-as.integer(Pathogen$InitialInfected);ne<-if(model=="SEIR")as.integer(Pathogen$InitialExposed)else 0L;nr<-if(model!="SIS")as.integer(Pathogen$InitialRecovered)else 0L
    if(ni+ne+nr>n)stop("Initial pathogen counts exceed points")
    pos<-seq_len(n);if(ne){z[pos[seq_len(ne)]]<-"E";pos<-pos[-seq_len(ne)]};if(ni){z[pos[seq_len(ni)]]<-"I";pos<-pos[-seq_len(ni)]};if(nr)z[pos[seq_len(nr)]]<-"R";InitialState<-z
  }
  InitialState<-as.character(InitialState);if(length(InitialState)!=n)stop("InitialState must have one pathogen state per point")
  dist<-numeric(nrow(st));ii<-match(.ina_point_state_key(InitialState),keys);if(is.na(ii))stop("InitialState contains invalid state");dist[ii]<-1
  full<-matrix(0,Ntimesteps+1,nrow(st),dimnames=list(timestep=0:Ntimesteps,state=keys));full[1,]<-dist
  active<-function(row)row=="I"|(row=="E"&model=="SEIR")
  EI<-matrix(0,Ntimesteps+1,n);EI[1,]<-InitialState=="I";EE<-if(model=="SEIR")matrix(0,Ntimesteps+1,n)else NULL;if(model=="SEIR")EE[1,]<-InitialState=="E"
  ext<-numeric(Ntimesteps+1);ext[1]<-as.numeric(!any(active(InitialState)))
  outside<-if(is.null(Outside))rep(FALSE,n)else as.logical(Outside);if(length(outside)!=n)stop("Outside must be logical length nrow(Points)")
  escstates<-if(any(outside))apply(st[,outside,drop=FALSE],1,function(z)any(z=="I"|(model=="SEIR"&z=="E")))else rep(FALSE,nrow(st))
  escape<-numeric(Ntimesteps+1);safe<-dist
  if(any(escstates&dist>0)){escape[1]<-sum(dist[escstates]);safe[escstates]<-0}else safe[escstates]<-0
  Tlist<-vector("list",Ntimesteps);growths<-vector("list",Ntimesteps);ngs<-vector("list",Ntimesteps)
  for(tt in seq_len(Ntimesteps)){
    Qtt<-INApestPointPathogenEdgeMatrix(Points,Pathogen,ContactRadius,ContactKernel,ContactProb,tt,1L)
    Ttt<-.ina_point_exact_operator(model,Points,Pathogen,Qtt,MaxStates,tt,1L)
    if(!identical(rownames(Ttt),keys))stop("Internal point exact state ordering changed across timesteps")
    Tlist[[tt]]<-Ttt
    growths[[tt]]<-INApestPointPathogenGrowthOperator(model,Points,Pathogen,ContactRadius,ContactKernel,ContactProb,tt,1L)
    ngs[[tt]]<-INApestPointPathogenNextGeneration(growths[[tt]])
    dist<-as.numeric(dist%*%Ttt);full[tt+1,]<-dist
    EI[tt+1,]<-as.numeric(dist%*%(st=="I"));if(model=="SEIR")EE[tt+1,]<-as.numeric(dist%*%(st=="E"))
    ext[tt+1]<-sum(dist[apply(st,1,function(z)!any(z=="I"|(model=="SEIR"&z=="E")))])
    if(any(outside)){
      nxt<-as.numeric(safe%*%Ttt);newesc<-sum(nxt[escstates]);escape[tt+1]<-escape[tt]+newesc;nxt[escstates]<-0;safe<-nxt
    }
  }
  ord<-INApestTransitionPathogenPeriodicGrowth(lapply(growths,`[[`,"Operator"))
  static_inputs<-!any(vapply(list(Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$IntroductionProb,Pathogen$DensityScale,ContactProb,ContactKernel),is.function,logical(1)))
  static_operator<-static_inputs&&all(vapply(Tlist[-1L],function(M)isTRUE(all.equal(M,Tlist[[1L]],tolerance=0)),logical(1)))
  if(length(Tlist)==1L)static_operator<-static_inputs
  intro_zero<-all(.ina_point_resolve(Pathogen$IntroductionProb,Points,1L,1L,"IntroductionProb")==0)
  absorption<-if(static_operator&&intro_zero).ina_point_static_absorption(Tlist[[1L]],st,model,outside,ii)else NULL
  out<-list(Model="INApestMetaPoint",PathogenModel=model,Exact=TRUE,Points=Points,
    StepGrowth=growths,StepNextGeneration=ngs,OrderedGrowth=ord,StaticAbsorption=absorption,
    ExpectedInfected=EI,ExpectedExposed=EE,ExtinctionProbability=ext,
    Escape=if(any(outside))list(Outside=outside,ProbabilityByTimestep=escape,ProbabilityByHorizon=tail(escape,1))else NULL,
    StateTable=st,StateDistribution=full,
    Diagnostics=c(
      "Exact finite-state pathogen solution conditional on a fixed point set and geometry.",
      "The exact operator integrates over realised contacts analytically and preserves the engine's one-state-transition-per-timestep timing.",
      "Timestep-specific pathogen/contact parameters are handled by an ordered sequence of exact transition operators rather than by averaging them.",
      if(!is.null(absorption)) "For static no-reintroduction systems, eventual escape, extinction-before-escape and mean absorption time are solved from the exact finite-state chain." else "Eventual static absorption metrics are omitted when parameters are time-varying or pathogen reintroduction is active.",
      "Pathogen mortality is excluded here because it removes host points and therefore changes geometry."))
  if(ReturnOperator)out$Operators<-Tlist
  class(out)<-c("INApestPointPathogenExactFixedGeometry","list");out
}

INApestPointTransitionPathogenGrowthOperator <- function(
    Model=c("SIS","SIR","SEIR"),Representatives,HostTypePopulationDiseaseStep,
    ParentCarrierOperator,Pathogen,ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,
    timestep=1L,perm=1L){
  Model<-match.arg(Model);R<-as.data.frame(Representatives);nt<-nrow(R);N<-as.numeric(HostTypePopulationDiseaseStep);H<-as.matrix(ParentCarrierOperator)
  if(length(N)!=nt||any(!is.finite(N))||any(N<0))stop("HostTypePopulationDiseaseStep must be finite non-negative length nrow(Representatives)")
  if(!all(dim(H)==c(nt,nt))||any(H<0)||any(!is.finite(H))||any(colSums(H)>1+1e-10))stop("ParentCarrierOperator must be a substochastic target x source matrix")
  Q<-INApestPointPathogenEdgeMatrix(R,Pathogen,ContactRadius,ContactKernel,ContactProb,timestep,perm,IncludeDiagonal=TRUE)
  # One-carrier offspring mean.  A carrier in source type i sees N_j
  # susceptible hosts in every other target type j, but only N_i - 1
  # susceptible hosts in its own type because the carrier itself occupies one host.
  SusceptibleTargets<-matrix(N,nrow=nt,ncol=nt)
  diag(SusceptibleTargets)<-pmax(0,N-1)
  B<-Q*SusceptibleTargets
  rec<-.ina_point_resolve(Pathogen$RecoveryProb,R,timestep,perm,"RecoveryProb");mort<-.ina_point_resolve(Pathogen$PathogenMortalityProb,R,timestep,perm,"PathogenMortalityProb");prog<-.ina_point_resolve(Pathogen$ProgressionProb,R,timestep,perm,"ProgressionProb")
  stayI<-diag(1-rec-mort,nt)
  if(Model%in%c("SIS","SIR")){D<-stayI+B;G<-D%*%H}
  else {D<-rbind(cbind(diag(1-prog,nt),B),cbind(diag(prog,nt),stayI));Z<-matrix(0,nt,nt);H2<-rbind(cbind(H,Z),cbind(Z,H));G<-D%*%H2}
  lam<-max(Mod(eigen(G,only.values=TRUE)$values))
  out<-list(Model="INApestPointTransitionMatrix",PathogenModel=Model,Operator=G,ParentCarrierOperator=H,
            DiseaseOperator=D,TransmissionEdgeMatrix=Q,TransmissionBlock=B,SusceptibleTargetsByCarrier=SusceptibleTargets,HostTypePopulationDiseaseStep=N,
            Representatives=R,Pathogen=Pathogen,
            ResolvedDiseaseParameters=list(RecoveryProb=rec,PathogenMortalityProb=mort,ProgressionProb=prog),
            IntrinsicRarePathogenMultiplier=lam,
            Classification=if(lam>1+1e-12)"growing"else if(lam<1-1e-12)"declining"else"threshold",
            Diagnostics=c(
              "Continuous point-stage movement is represented through the host analytical Parent carrier matrix; reproductive recruits are excluded from pathogen carriage because offspring are susceptible.",
              "Transmission is applied after stage movement/reproduction using the post-host-step type abundance and the point contact geometry.",
              "Within-type transmission is retained for contracted types containing multiple hosts; the infectious carrier itself is removed from its own susceptible target count.",
              "When Representatives are analysis-grid cell centres, within-cell contact is evaluated at zero representative distance and is therefore a grid-resolution approximation to the continuous-space contact process."))
  class(out)<-c("INApestPointTransitionPathogenGrowthOperator","list");out
}

INApestPointTransitionPathogenNextGeneration <- function(GrowthObject){
  x<-GrowthObject;H<-x$ParentCarrierOperator;B<-x$TransmissionBlock;D<-x$DiseaseOperator;nt<-nrow(H);Z<-matrix(0,nt,nt)
  if(x$PathogenModel%in%c("SIS","SIR")){T<-(D-B)%*%H;F<-B%*%H}
  else {Fdis<-rbind(cbind(Z,B),cbind(Z,Z));Tdis<-D-Fdis;H2<-rbind(cbind(H,Z),cbind(Z,H));T<-Tdis%*%H2;F<-Fdis%*%H2}
  rt<-max(Mod(eigen(T,only.values=TRUE)$values));if(rt>=1-1e-12){K<-matrix(NA,nrow(T),ncol(T));R0<-Inf}else{K<-F%*%solve(diag(nrow(T))-T);R0<-max(Mod(eigen(K,only.values=TRUE)$values))}
  list(R0=R0,NextGenerationOperator=K,TransitionWithoutNewInfection=T,NewInfectionOperator=F,
       PerTimestepLambda=x$IntrinsicRarePathogenMultiplier,
       ThresholdAgreement=if(is.finite(R0))sign(R0-1)==sign(x$IntrinsicRarePathogenMultiplier-1)else x$IntrinsicRarePathogenMultiplier>=1)
}

INApestPointTransitionPathogenBranching <- function(GrowthObject,timesteps,OutsideTypes=integer(0)){
  x<-GrowthObject;model<-x$PathogenModel;H<-x$ParentCarrierOperator;Q<-x$TransmissionEdgeMatrix;N<-x$HostTypePopulationDiseaseStep;nt<-nrow(H);Tn<-as.integer(timesteps)
  rec<-x$ResolvedDiseaseParameters$RecoveryProb;mort<-x$ResolvedDiseaseParameters$PathogenMortalityProb;prog<-x$ResolvedDiseaseParameters$ProgressionProb
  outside<-unique(as.integer(OutsideTypes));if(any(!outside%in%seq_len(nt)))stop("OutsideTypes must index analytical host types")
  integerN<-all(abs(N-round(N))<1e-10)
  trans_factor<-function(v,zE){
    Ssus<-pmax(0,N-as.numeric(seq_len(nt)==v));q<-Q[,v]
    if(integerN)prod((1-q+q*zE)^as.integer(round(Ssus))) else exp(sum(Ssus*q*(zE-1)))
  }
  if(model%in%c("SIS","SIR")){
    q<-rep(0,nt);h<-rep(1,nt);if(length(outside))h[outside]<-0;qh<-hh<-matrix(NA,nt,Tn)
    stepfun<-function(z){zn<-numeric(nt);for(u in seq_len(nt)){
      val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0)){stay<-1-rec[v]-mort[v];val<-val+H[v,u]*((1-stay)+stay*z[v])*trans_factor(v,z)};zn[u]<-val};pmin(1,pmax(0,zn))}
    for(tt in seq_len(Tn)){q<-stepfun(q);h<-stepfun(h);if(length(outside))h[outside]<-0;qh[,tt]<-q;hh[,tt]<-h}
  }else{
    q<-rep(0,2*nt);h<-rep(1,2*nt);if(length(outside))h[c(outside,nt+outside)]<-0;qh<-hh<-matrix(NA,2*nt,Tn)
    stepfun<-function(z){zn<-numeric(2*nt)
      for(u in seq_len(nt)){val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0))val<-val+H[v,u]*((1-prog[v])*z[v]+prog[v]*z[nt+v]);zn[u]<-val}
      for(u in seq_len(nt)){val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0)){stay<-1-rec[v]-mort[v];val<-val+H[v,u]*((1-stay)+stay*z[nt+v])*trans_factor(v,z[seq_len(nt)])};zn[nt+u]<-val};pmin(1,pmax(0,zn))}
    for(tt in seq_len(Tn)){q<-stepfun(q);h<-stepfun(h);if(length(outside))h[c(outside,nt+outside)]<-0;qh[,tt]<-q;hh[,tt]<-h}
  }
  inside<-setdiff(seq_len(nt),outside);itypes<-if(model=="SEIR")c(inside,nt+inside)else inside
  esc<-1-h
  list(Model="INApestPointTransitionMatrix",PathogenModel=model,ExtinctionByHorizon=q,EscapeByHorizon=esc,
       ExtinctionHistory=qh,EscapeHistory=1-hh,OutsideTypes=outside,InsideStartingTypes=itypes,
       MaxEscapeByHorizonFromInside=if(length(itypes))max(esc[itypes])else NA_real_,
       Approximation=if(integerN)"finite-count one-lineage branching on contracted point types"else"mean-matched Poisson branching on contracted point types",
       Diagnostics=c("Branching recursion composes stage/space carrier movement with point-contact transmission after movement.",
         "When analytical type abundances are non-integer expectations, transmission offspring use a mean-matched Poisson PGF; integer type abundances use the finite-count binomial PGF.",
         "For containment, include explicit outside analysis-grid cells in OutsideTypes; movement that leaves the analytical grid entirely is treated as loss rather than silently counted as pathogen escape."))
}

INApestPointTransitionPathogenAnalytical <- function(
    Ntimesteps=10,Nstages,Transition,InitialPoints,Pathogen,
    InformationMode=c("none","all_informed"),
    SDDkernel,LDDkernel=NULL,LDDrate=0,PropaguleEstablishment=1,EnvEstabProb=1,
    TransitionKernels=NULL,TransitionHabitatSearch=FALSE,ApplyHabitatToTransitions=FALSE,
    TransitionEstablishment=1,BlockedTransitionMortality=0,
    HabitatSuitability=NULL,HabitatSearchRadius=0,HabitatSearchCandidates=128,LocalK=Inf,KRadius=0,
    DetectionProb=0,DetectionSpatial=NULL,ManageProb=0,ManageSpatial=NULL,
    MortalityProb=0,MortalitySpatial=NULL,FecundityReduction=0,FecundityReductionSpatial=NULL,
    SpreadReduction=0,SpreadReductionSpatial=NULL,SpreadReductionAppliesTo=c("LDD","all"),
    InfoRetentionProb=1,InfoRadius=0,InfoTransferProb=0,InfoKernel=NULL,
    PointAnalysisGrid=NULL,MaxAnalysisCells=400L,KernelSamples=1000L,PointSeed=1L,
    OutsideEstablishmentProb=1,
    ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,
    OutsideCells=integer(0),ReturnOperators=FALSE){
  cr_missing<-missing(ContactRadius);ck_missing<-missing(ContactKernel);cp_missing<-missing(ContactProb)
  interaction<-inherits(Pathogen,"INApestPointPathogenInteraction")
  Pspec<-.ina_point_pathogen_spec(Pathogen)
  if(interaction&&cr_missing)ContactRadius<-.ina_point_interaction_contact(Pathogen,"ContactRadius",Inf)
  if(interaction&&ck_missing)ContactKernel<-.ina_point_interaction_contact(Pathogen,"ContactKernel",NULL)
  if(interaction&&cp_missing)ContactProb<-.ina_point_interaction_contact(Pathogen,"ContactProb",1)
  Pathogen<-Pspec
  InformationMode<-match.arg(InformationMode);SpreadReductionAppliesTo<-match.arg(SpreadReductionAppliesTo)
  Ntimesteps<-as.integer(Ntimesteps);Nstages<-as.integer(Nstages)
  if(length(Ntimesteps)!=1L||is.na(Ntimesteps)||Ntimesteps<1L)stop("Ntimesteps must be a positive integer")
  if(length(Nstages)!=1L||is.na(Nstages)||Nstages<1L)stop("Nstages must be a positive integer")
  model<-Pathogen$Model;if(!model%in%c("SIS","SIR","SEIR"))stop("Point-transition pathogen analytical methods support SIS/SIR/SEIR")
  .ina_pt_require_point_helpers(TRUE)
  analysis_grid<-.ina_pt_analysis_grid(PointAnalysisGrid,HabitatSuitability,DetectionSpatial,ManageSpatial,MortalitySpatial,FecundityReductionSpatial,SpreadReductionSpatial)
  if(is.null(analysis_grid)&&(!is.null(TransitionKernels)||is.function(ContactKernel)||is.finite(ContactRadius)))
    stop("Spatial transition/contact analysis requires PointAnalysisGrid so movement and contact geometry can be represented")
  if(!is.null(analysis_grid)&&is.null(PointAnalysisGrid)&&analysis_grid$nrow*analysis_grid$ncol>MaxAnalysisCells)stop("Inherited analysis grid exceeds MaxAnalysisCells; supply a coarser PointAnalysisGrid")
  reps<-.ina_pt_transition_representatives(InitialPoints,analysis_grid,Nstages)
  x<-.ina_pt_transition_initial(InitialPoints,analysis_grid,Nstages)
  mode<-if(InformationMode=="none")"none"else"all_informed"
  steps<-vector("list",Ntimesteps);growths<-vector("list",Ntimesteps);ngs<-vector("list",Ntimesteps);branches<-vector("list",Ntimesteps)
  hosttraj<-matrix(0,length(x),Ntimesteps+1L);hosttraj[,1L]<-x
  outside_types<-integer(0)
  if(length(OutsideCells)){
    if(is.null(analysis_grid))stop("OutsideCells requires PointAnalysisGrid")
    nc<-analysis_grid$nrow*analysis_grid$ncol;if(any(!OutsideCells%in%seq_len(nc)))stop("OutsideCells contains invalid analysis cells")
    outside_types<-unlist(lapply(seq_len(Nstages),function(s)(s-1L)*nc+OutsideCells),use.names=FALSE)
  }
  for(tt in seq_len(Ntimesteps)){
    st<-.ina_pt_transition_operator_step(tt,Ntimesteps,reps,analysis_grid,Transition,Nstages,
      SDDkernel,LDDkernel,LDDrate,PropaguleEstablishment,EnvEstabProb,
      TransitionKernels,TransitionHabitatSearch,ApplyHabitatToTransitions,TransitionEstablishment,
      BlockedTransitionMortality,HabitatSuitability,HabitatSearchRadius,HabitatSearchCandidates,
      DetectionProb,DetectionSpatial,ManageProb,ManageSpatial,MortalityProb,MortalitySpatial,
      FecundityReduction,FecundityReductionSpatial,SpreadReduction,SpreadReductionSpatial,
      SpreadReductionAppliesTo,InfoRetentionProb,InfoRadius,InfoTransferProb,InfoKernel,
      KernelSamples,PointSeed,mode,OutsideEstablishmentProb)
    Ghost<-if(mode=="none")st$G0 else st$GH;H<-if(mode=="none")st$Parent0 else st$ParentH
    x1<-as.numeric(Ghost%*%x);hosttraj[,tt+1L]<-x1
    pg<-INApestPointTransitionPathogenGrowthOperator(model,reps,x1,H,Pathogen,ContactRadius,ContactKernel,ContactProb,tt,1L)
    growths[[tt]]<-pg;ngs[[tt]]<-INApestPointTransitionPathogenNextGeneration(pg)
    branches[[tt]]<-INApestPointTransitionPathogenBranching(pg,1L,outside_types)
    steps[[tt]]<-st;x<-x1
  }
  ops<-lapply(growths,`[[`,"Operator");cyc<-INApestTransitionPathogenPeriodicGrowth(ops)
  brseq<-INApestPointPathogenBranchingSequence(growths,outside_types)
  result<-list(Model="INApestPointTransitionMatrix",PathogenModel=model,InformationMode=InformationMode,
    Representatives=reps,HostTypeTrajectory=hosttraj,StepGrowth=growths,StepNextGeneration=ngs,
    CycleGrowth=cyc,Branching=brseq,StepBranching=branches,OutsideTypes=outside_types,
    Diagnostics=c(
      "Host stage/space movement and reproduction are contracted with the existing INApest point-transition analytical kernel machinery; only the persisting/progressing Parent operator carries pathogen state.",
      "Reproductive offspring contribute to the susceptible host background before pathogen contact but never inherit pathogen state.",
      "Pathogen contact between contracted stage x grid types is evaluated at representative cell-centre geometry; refine PointAnalysisGrid where contact kernels vary sharply with distance.",
      "Dynamic detection-driven information is intentionally not collapsed into an average management probability here; use InformationMode='none' or 'all_informed' for analytical screening and the stochastic point model for endogenous information feedback.",
      if(!is.infinite(LocalK)||KRadius>0) "Finite LocalK/KRadius is a nonlinear local-crowding process and is excluded from the rare-pathogen carrier/background operator; use stochastic simulation once crowding is material." else "Local crowding is inactive in this analytical screening configuration."))
  if(ReturnOperators)result$HostStepData<-steps
  class(result)<-c("INApestPointTransitionPathogenAnalytical","list");result
}

print.INApestPointTransitionPathogenAnalytical <- function(x,...){
  cat("INApest point-transition pathogen analytical result\n")
  cat("  Pathogen model:",x$PathogenModel,"\n")
  cat("  Equivalent per-timestep multiplier:",format(x$CycleGrowth$PerTimestepMultiplier,digits=7),"\n")
  invisible(x)
}

###############################################################################
### Time-inhomogeneous branching and high-level INApestMetaPoint pathogen path
###############################################################################

.ina_point_pathogen_spec <- function(Pathogen){
  if(inherits(Pathogen,"INApestPointPathogenInteraction")){
    if(is.null(Pathogen$Pathogen)||!inherits(Pathogen$Pathogen,"INApestPathogen"))
      stop("Point pathogen interaction does not contain a valid INApestPathogen specification")
    return(Pathogen$Pathogen)
  }
  if(inherits(Pathogen,"INApestPathogen"))return(Pathogen)
  stop("Pathogen must be an INApestPathogen or INApestPointPathogenInteraction")
}

.ina_point_interaction_contact <- function(Pathogen,name,default){
  if(!inherits(Pathogen,"INApestPointPathogenInteraction"))return(default)
  e<-environment(Pathogen$Contact)
  if(!is.null(e)&&exists(name,envir=e,inherits=TRUE))get(name,envir=e,inherits=TRUE) else default
}

.ina_point_branch_step <- function(GrowthObject,z){
  x<-GrowthObject;model<-x$PathogenModel;H<-x$ParentCarrierOperator;Q<-x$TransmissionEdgeMatrix
  N<-x$HostTypePopulationDiseaseStep;nt<-nrow(H)
  rec<-x$ResolvedDiseaseParameters$RecoveryProb;mort<-x$ResolvedDiseaseParameters$PathogenMortalityProb
  prog<-x$ResolvedDiseaseParameters$ProgressionProb
  integerN<-all(abs(N-round(N))<1e-10)
  trans_factor<-function(v,zE){
    Ssus<-pmax(0,N-as.numeric(seq_len(nt)==v));q<-Q[,v]
    if(integerN)prod((1-q+q*zE)^as.integer(round(Ssus))) else exp(sum(Ssus*q*(zE-1)))
  }
  if(model%in%c("SIS","SIR")){
    if(length(z)!=nt)stop("Branching state length mismatch")
    zn<-numeric(nt)
    for(u in seq_len(nt)){
      val<-max(0,1-sum(H[,u]))
      for(v in which(H[,u]>0)){
        stay<-1-rec[v]-mort[v]
        val<-val+H[v,u]*((1-stay)+stay*z[v])*trans_factor(v,z)
      }
      zn[u]<-val
    }
  } else {
    if(length(z)!=2L*nt)stop("SEIR branching state length mismatch")
    zn<-numeric(2L*nt)
    for(u in seq_len(nt)){
      val<-max(0,1-sum(H[,u]))
      for(v in which(H[,u]>0))val<-val+H[v,u]*((1-prog[v])*z[v]+prog[v]*z[nt+v])
      zn[u]<-val
    }
    for(u in seq_len(nt)){
      val<-max(0,1-sum(H[,u]))
      for(v in which(H[,u]>0)){
        stay<-1-rec[v]-mort[v]
        val<-val+H[v,u]*((1-stay)+stay*z[nt+v])*trans_factor(v,z[seq_len(nt)])
      }
      zn[nt+u]<-val
    }
  }
  pmin(1,pmax(0,zn))
}

INApestPointPathogenBranchingSequence <- function(GrowthObjects,OutsideTypes=integer(0)){
  if(!is.list(GrowthObjects)||!length(GrowthObjects))stop("GrowthObjects must be a non-empty list")
  model<-GrowthObjects[[1L]]$PathogenModel;nt<-nrow(GrowthObjects[[1L]]$ParentCarrierOperator)
  if(any(vapply(GrowthObjects,function(x)x$PathogenModel!=model||nrow(x$ParentCarrierOperator)!=nt,logical(1))))
    stop("All GrowthObjects must use the same pathogen model and analytical type set")
  d<-if(model=="SEIR")2L*nt else nt
  outside<-unique(as.integer(OutsideTypes));if(any(!outside%in%seq_len(nt)))stop("OutsideTypes must index analytical host types")
  outside_active<-if(model=="SEIR")c(outside,nt+outside)else outside
  Tn<-length(GrowthObjects);qh<-hh<-matrix(NA_real_,d,Tn)
  for(hor in seq_len(Tn)){
    q<-rep(0,d)
    for(tt in rev(seq_len(hor)))q<-.ina_point_branch_step(GrowthObjects[[tt]],q)
    qh[,hor]<-q
    z<-rep(1,d);if(length(outside_active))z[outside_active]<-0
    for(tt in rev(seq_len(hor))){z<-.ina_point_branch_step(GrowthObjects[[tt]],z);if(length(outside_active))z[outside_active]<-0}
    hh[,hor]<-z
  }
  inside<-setdiff(seq_len(nt),outside);inside_active<-if(model=="SEIR")c(inside,nt+inside)else inside
  list(Model=GrowthObjects[[1L]]$Model,PathogenModel=model,
       ExtinctionHistory=qh,EscapeHistory=1-hh,
       ExtinctionByHorizon=qh[,Tn],EscapeByHorizon=1-hh[,Tn],
       OutsideTypes=outside,InsideStartingTypes=inside_active,
       MaxEscapeByHorizonFromInside=if(length(inside_active))max((1-hh[,Tn])[inside_active])else NA_real_,
       Diagnostics=c(
         "Time-varying branching PGFs are composed in chronological order: F1(F2(...FT(z)...)), rather than by averaging timestep-specific extinction or escape probabilities.",
         "Integer contracted host counts use finite-count binomial transmission offspring; fractional analytical host backgrounds use a mean-matched Poisson PGF.",
         "Escape is first passage to an active pathogen state in an outside analytical type at a timestep boundary."))
}

INApestMetaPointPathogenAnalytical <- function(
    Ntimesteps=10,InitialPoints,Pathogen,
    InformationMode=c("none","all_informed"),
    Survival=1,PropaguleProduction,PropaguleEstablishment=1,EnvEstabProb=1,
    SDDkernel,LDDkernel=NULL,LDDrate=0,
    HabitatSuitability=NULL,HabitatSearchRadius=0,HabitatSearchCandidates=128,LocalK=Inf,KRadius=0,
    DetectionProb=0,DetectionSpatial=NULL,ManageProb=0,ManageSpatial=NULL,
    MortalityProb=0,MortalitySpatial=NULL,FecundityReduction=0,FecundityReductionSpatial=NULL,
    SpreadReduction=0,SpreadReductionSpatial=NULL,SpreadReductionAppliesTo=c("LDD","all"),
    InfoRetentionProb=1,InfoRadius=0,InfoTransferProb=0,InfoKernel=NULL,
    PointAnalysisGrid=NULL,MaxAnalysisCells=400L,KernelSamples=1000L,PointSeed=1L,
    OutsideEstablishmentProb=1,
    ContactRadius=Inf,ContactKernel=NULL,ContactProb=1,
    OutsideCells=integer(0),ReturnOperators=FALSE){
  cr_missing<-missing(ContactRadius);ck_missing<-missing(ContactKernel);cp_missing<-missing(ContactProb)
  interaction<-inherits(Pathogen,"INApestPointPathogenInteraction")
  Pspec<-.ina_point_pathogen_spec(Pathogen)
  if(interaction&&cr_missing)ContactRadius<-.ina_point_interaction_contact(Pathogen,"ContactRadius",Inf)
  if(interaction&&ck_missing)ContactKernel<-.ina_point_interaction_contact(Pathogen,"ContactKernel",NULL)
  if(interaction&&cp_missing)ContactProb<-.ina_point_interaction_contact(Pathogen,"ContactProb",1)
  Pathogen<-Pspec
  InformationMode<-match.arg(InformationMode);SpreadReductionAppliesTo<-match.arg(SpreadReductionAppliesTo)
  Ntimesteps<-as.integer(Ntimesteps);if(length(Ntimesteps)!=1L||is.na(Ntimesteps)||Ntimesteps<1L)stop("Ntimesteps must be a positive integer")
  .ina_pt_require_point_helpers(FALSE)
  analysis_grid<-.ina_pt_analysis_grid(PointAnalysisGrid,HabitatSuitability,DetectionSpatial,ManageSpatial,MortalitySpatial,FecundityReductionSpatial,SpreadReductionSpatial)
  if(is.null(analysis_grid)&&(is.function(ContactKernel)||is.finite(ContactRadius)))
    stop("Finite-radius or geometry-dependent pathogen contact requires PointAnalysisGrid for the scalable moving-point analytical contraction")
  if(!is.null(analysis_grid)&&is.null(PointAnalysisGrid)&&analysis_grid$nrow*analysis_grid$ncol>MaxAnalysisCells)
    stop("Inherited analysis grid exceeds MaxAnalysisCells; supply a coarser PointAnalysisGrid")
  reps<-.ina_pt_metapoint_representatives(InitialPoints,analysis_grid)
  x<-.ina_pt_metapoint_initial(InitialPoints,analysis_grid)
  mode<-if(InformationMode=="none")"none"else"all_informed"
  steps<-growths<-ngs<-vector("list",Ntimesteps)
  hosttraj<-matrix(0,length(x),Ntimesteps+1L);hosttraj[,1L]<-x
  outside_types<-integer(0)
  if(length(OutsideCells)){
    if(is.null(analysis_grid))stop("OutsideCells requires PointAnalysisGrid")
    nc<-analysis_grid$nrow*analysis_grid$ncol;if(any(!OutsideCells%in%seq_len(nc)))stop("OutsideCells contains invalid analysis cells")
    outside_types<-unique(as.integer(OutsideCells))
  }
  for(tt in seq_len(Ntimesteps)){
    st<-.ina_pt_metapoint_operator_step(tt,Ntimesteps,reps,analysis_grid,
      Survival,PropaguleProduction,PropaguleEstablishment,EnvEstabProb,
      SDDkernel,LDDkernel,LDDrate,HabitatSuitability,HabitatSearchRadius,HabitatSearchCandidates,
      DetectionProb,DetectionSpatial,ManageProb,ManageSpatial,MortalityProb,MortalitySpatial,
      FecundityReduction,FecundityReductionSpatial,SpreadReduction,SpreadReductionSpatial,
      SpreadReductionAppliesTo,InfoRetentionProb,InfoRadius,InfoTransferProb,InfoKernel,
      KernelSamples,PointSeed,mode,OutsideEstablishmentProb)
    Ghost<-if(mode=="none")st$G0 else st$GH;H<-if(mode=="none")st$Parent0 else st$ParentH
    x1<-as.numeric(Ghost%*%x);hosttraj[,tt+1L]<-x1
    pg<-INApestPointTransitionPathogenGrowthOperator(Pathogen$Model,reps,x1,H,Pathogen,ContactRadius,ContactKernel,ContactProb,tt,1L)
    pg$Model<-"INApestMetaPoint";growths[[tt]]<-pg;ngs[[tt]]<-INApestPointTransitionPathogenNextGeneration(pg)
    steps[[tt]]<-st;x<-x1
  }
  cyc<-INApestTransitionPathogenPeriodicGrowth(lapply(growths,`[[`,"Operator"))
  br<-INApestPointPathogenBranchingSequence(growths,outside_types)
  out<-list(Model="INApestMetaPoint",PathogenModel=Pathogen$Model,InformationMode=InformationMode,
    Representatives=reps,HostTypeTrajectory=hosttraj,StepGrowth=growths,StepNextGeneration=ngs,
    OrderedGrowth=cyc,Branching=br,OutsideTypes=outside_types,
    Diagnostics=c(
      "Existing infected points are carried only by the host Parent survival/management operator; Poisson reproductive offspring enter the susceptible host background and never inherit pathogen state.",
      "Pathogen contact is applied after survival, management and recruitment, matching INApestMetaPoint event order.",
      "The scalable moving-point solution supports none/all-informed management backgrounds; endogenous detection-information feedback remains a correlated point process and is intentionally left to stochastic simulation.",
      if(!is.infinite(LocalK)||KRadius>0) "Finite LocalK/KRadius is a nonlinear local-crowding process and is excluded from the rare-pathogen background operator; use stochastic simulation once crowding is material." else "Local crowding is inactive in this analytical screening configuration.",
      if(is.null(analysis_grid)) "Homogeneous contact contraction: all hosts occupy one analytical type." else "Spatial contact and host dispersal are contracted to PointAnalysisGrid cell-centre types; refine the grid where contact kernels or habitat vary sharply."))
  if(ReturnOperators)out$HostStepData<-steps
  class(out)<-c("INApestMetaPointPathogenAnalytical","list");out
}

print.INApestMetaPointPathogenAnalytical <- function(x,...){
  cat("INApest MetaPoint pathogen analytical result\n")
  cat("  Pathogen model:",x$PathogenModel,"\n")
  cat("  Equivalent per-timestep multiplier:",format(x$OrderedGrowth$PerTimestepMultiplier,digits=7),"\n")
  if(length(x$OutsideTypes))cat("  Branching escape by horizon (max inside starting type):",format(x$Branching$MaxEscapeByHorizonFromInside,digits=7),"\n")
  invisible(x)
}
