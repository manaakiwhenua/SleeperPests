###############################################################################
### INApestMetaTransitionMatrix pathogen analytical methods
### Core mathematics: demographic carrier transport -> pathogen state process
###############################################################################

.ina_tm_ns_matrix <- function(x,n,S,name,prob=FALSE){
  if(is.matrix(x)){
    if(!all(dim(x)==c(n,S))) stop(name," matrix must be nodes x stages")
    z<-x
  } else {
    x<-as.numeric(x)
    if(length(x)==1L) z<-matrix(x,n,S)
    else if(length(x)==S && S!=n) z<-matrix(rep(x,each=n),n,S)
    else if(length(x)==n && n!=S) z<-matrix(rep(x,S),n,S)
    else if(length(x)==n*S) z<-matrix(x,n,S)
    else stop(name," must be scalar, stage vector, node vector, or nodes x stages")
  }
  if(any(!is.finite(z))) stop(name," must be finite")
  if(prob && any(z<0|z>1)) stop(name," must lie in [0,1]")
  z
}

.ina_tm_transition_movement <- function(x,S,n,name){
  if(S<=1L) return(list())
  if(is.null(x)) return(rep(list(NULL),S-1L))
  z<-if(is.list(x))x else rep(list(x),S-1L)
  if(length(z)!=S-1L) stop(name," must be a matrix or list of length Nstages-1")
  for(k in seq_along(z)) if(!is.null(z[[k]])){
    z[[k]]<-as.matrix(z[[k]])
    if(!all(dim(z[[k]])==c(n,n))) stop(name," matrices must be nodes x nodes")
    if(any(!is.finite(z[[k]]))||any(z[[k]]<0)||any(rowSums(z[[k]])>1+1e-12))
      stop(name," rows must contain finite non-negative probabilities summing to <=1")
  }
  z
}

INApestTransitionPathogenCarrierOperator <- function(
    Transition,Nstages,n_nodes,
    ManageProb=0,MortalityProb=0,
    TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0){
  S<-as.integer(Nstages);n<-as.integer(n_nodes)
  if(S<1L||n<1L)stop("Nstages and n_nodes must be positive")
  A<-.transition_list(Transition,n,S)
  a<-.ina_recycle(ManageProb,n,"ManageProb")
  M<-.ina_tm_ns_matrix(MortalityProb,n,S,"MortalityProb",TRUE)
  if(any(a<0|a>1))stop("ManageProb must lie in [0,1]")
  Ps<-.ina_tm_transition_movement(TransitionSDDprob,S,n,"TransitionSDDprob")
  Pl<-.ina_tm_transition_movement(TransitionLDDprob,S,n,"TransitionLDDprob")
  rr<-rep_len(as.numeric(TransitionLDDrate),max(1L,S-1L))
  if(any(!is.finite(rr))||any(rr<0|rr>1))stop("TransitionLDDrate must lie in [0,1]")
  U<-n*S;idx<-function(i,s)(i-1L)*S+s
  H<-matrix(0,U,U)
  # Management occurs before demographic progression. Integrating over the
  # Bernoulli adoption event gives the one-carrier survival factor below.
  q<-1-M*sweep(matrix(1,n,S),1,a,`*`)
  for(i in seq_len(n))for(k in seq_len(S)){
    src<-idx(i,k);Ai<-A[[i]];qq<-q[i,k]
    if(k<S){
      H[idx(i,k),src]<-H[idx(i,k),src]+qq*Ai[k,k]
      tr<-Ai[k+1L,k]
      if(tr>0){
        if(is.null(Ps[[k]])&&is.null(Pl[[k]])) P<-diag(n)
        else if(!is.null(Ps[[k]])&&!is.null(Pl[[k]])) P<-(1-rr[k])*Ps[[k]]+rr[k]*Pl[[k]]
        else if(!is.null(Ps[[k]])) P<-Ps[[k]] else P<-Pl[[k]]
        for(j in seq_len(n))if(P[i,j]>0)
          H[idx(j,k+1L),src]<-H[idx(j,k+1L),src]+qq*tr*P[i,j]
      }
    } else H[idx(i,S),src]<-H[idx(i,S),src]+qq*Ai[S,S]
  }
  map<-expand.grid(stage=seq_len(S),node=seq_len(n),KEEP.OUT.ATTRS=FALSE)
  # expand.grid stage-fast order already matches (node-1)*S+stage.
  labels<-paste0("n",map$node,"_s",map$stage)
  dimnames(H)<-list(target=labels,source=labels)
  list(Operator=H,UnitMap=map,ColumnSurvival=rowSums(t(H)),
       Diagnostics=c(
         "Carrier operator tracks only existing hosts: demographic fecundity is excluded because all offspring enter pathogen state S.",
         "Management mortality is applied before demographic stage survival/progression, matching INApestMetaTransitionMatrix.",
         "Stage-transition movement is included exactly in the uncrowded limit; residual row probability is export/loss.",
         "Capacity blocking and BlockedTransitionMortality are intentionally excluded from this rare-carrier operator and require a crowded-host correction."))
}

INApestTransitionPathogenDiseaseLinearOperator <- function(
    Model=c("SIS","SIR","SEIR"),HostPopulationDiseaseStep,
    Beta,RecoveryProb,ProgressionProb=1,PathogenMortalityProb=0,
    ContactMatrix=NULL,StageMixing=NULL,Transmission=c("frequency","density"),DensityScale=1){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission)
  Nmat<-as.matrix(HostPopulationDiseaseStep);n<-nrow(Nmat);S<-ncol(Nmat);U<-n*S
  if(any(!is.finite(Nmat))||any(Nmat<0))stop("HostPopulationDiseaseStep must be finite non-negative nodes x stages")
  # Stochastic engine is node-major, stage-fast.
  N<-as.numeric(t(Nmat))
  beta<-as.numeric(t(.ina_tm_ns_matrix(Beta,n,S,"Beta")))
  rec<-as.numeric(t(.ina_tm_ns_matrix(RecoveryProb,n,S,"RecoveryProb",TRUE)))
  prog<-as.numeric(t(.ina_tm_ns_matrix(ProgressionProb,n,S,"ProgressionProb",TRUE)))
  mort<-as.numeric(t(.ina_tm_ns_matrix(PathogenMortalityProb,n,S,"PathogenMortalityProb",TRUE)))
  ds<-as.numeric(t(.ina_tm_ns_matrix(DensityScale,n,S,"DensityScale")))
  if(any(beta<0)||any(ds<=0)||any(rec+mort>1+1e-12))stop("Invalid pathogen parameters")
  Cn<-if(is.null(ContactMatrix))diag(n)else as.matrix(ContactMatrix)
  if(!all(dim(Cn)==c(n,n))||any(!is.finite(Cn))||any(Cn<0))stop("ContactMatrix must be finite non-negative nodes x nodes")
  Sm<-if(is.null(StageMixing))matrix(1,S,S)else as.matrix(StageMixing)
  if(!all(dim(Sm)==c(S,S))||any(!is.finite(Sm))||any(Sm<0))stop("StageMixing must be finite non-negative stages x stages")
  C<-kronecker(Cn,Sm) # source x target, node-major stage-fast
  B<-matrix(0,U,U)
  if(Transmission=="frequency"){
    den<-as.numeric(crossprod(N,C))
    for(w in seq_len(U))if(N[w]>0&&den[w]>0&&beta[w]>0)B[w,]<-N[w]*beta[w]*C[,w]/den[w]
  }else for(w in seq_len(U))if(N[w]>0&&beta[w]>0)B[w,]<-N[w]*beta[w]*C[,w]/ds[w]
  stayI<-diag(1-rec-mort,U)
  if(Model%in%c("SIS","SIR")){
    G<-stayI+B
    note<-if(Model=="SIR")"SIR has the same active-I disease-step linearisation as SIS at the pathogen-free state; R is inactive to first order."else"SIS disease-step linearisation."
  } else {
    P<-diag(prog,U);stayE<-diag(1-prog,U)
    G<-rbind(cbind(stayE,B),cbind(P,stayI))
    note<-"SEIR preserves the discrete latent delay: infections generated this step enter E and cannot progress until a later pathogen step."
  }
  list(Operator=G,TransmissionBlock=B,CombinedContact=C,HostPopulationDiseaseStep=Nmat,
       Model=Model,Transmission=Transmission,Diagnostics=c(
         "Rows are recipient active pathogen states; columns are source active pathogen states.",
         "Frequency-transmission denominators use the post-demographic disease-step host distribution, matching the stochastic engine.",note))
}

INApestTransitionPathogenGrowthOperator <- function(
    Model=c("SIS","SIR","SEIR"),Transition,Nstages,HostPopulationDiseaseStep,
    Beta,RecoveryProb,ProgressionProb=1,PathogenMortalityProb=0,
    ContactMatrix=NULL,StageMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
    ManageProb=0,MortalityProb=0,
    TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission)
  Nmat<-as.matrix(HostPopulationDiseaseStep);n<-nrow(Nmat);S<-ncol(Nmat)
  if(S!=Nstages)stop("HostPopulationDiseaseStep columns must equal Nstages")
  car<-INApestTransitionPathogenCarrierOperator(Transition,S,n,ManageProb,MortalityProb,
                                                 TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  dis<-INApestTransitionPathogenDiseaseLinearOperator(Model,Nmat,Beta,RecoveryProb,ProgressionProb,
                                                       PathogenMortalityProb,ContactMatrix,StageMixing,
                                                       Transmission,DensityScale)
  H<-car$Operator;U<-n*S
  if(Model%in%c("SIS","SIR"))G<-dis$Operator%*%H
  else {
    H2<-rbind(cbind(H,matrix(0,U,U)),cbind(matrix(0,U,U),H))
    G<-dis$Operator%*%H2
  }
  ev<-eigen(G,only.values=TRUE)$values;lambda<-max(Mod(ev))
  list(Model="INApestMetaTransitionMatrix",PathogenModel=Model,Operator=G,
       CarrierOperator=H,DiseaseOperator=dis$Operator,TransmissionBlock=dis$TransmissionBlock,
       IntrinsicRarePathogenMultiplier=lambda,
       Classification=if(lambda>1+1e-12)"growing"else if(lambda<1-1e-12)"declining"else"threshold",
       HostPopulationDiseaseStep=Nmat,UnitMap=car$UnitMap,
       Diagnostics=c(car$Diagnostics,dis$Diagnostics,
         "Full-timestep active-pathogen operator is disease-step operator composed after demographic carrier transport.",
         "This is a rare-pathogen/uncrowded carrier result. Demographic births affect the host background and transmission denominators but do not directly create infected offspring."))
}

INApestTransitionPathogenPeriodicGrowth <- function(Operators){
  if(!is.list(Operators)||!length(Operators))stop("Operators must be a non-empty list of same-sized square matrices")
  mats<-lapply(Operators,function(x)if(is.list(x)&&!is.null(x$Operator))x$Operator else as.matrix(x))
  d<-dim(mats[[1L]]);if(d[1]!=d[2]||any(vapply(mats,function(M)!all(dim(M)==d),logical(1))))stop("All operators must have the same square dimension")
  P<-diag(d[1]);for(G in mats)P<-G%*%P
  rho<-max(Mod(eigen(P,only.values=TRUE)$values));g<-rho^(1/length(mats))
  list(CycleOperator=P,CycleMultiplier=rho,PerTimestepMultiplier=g,
       Classification=if(g>1+1e-12)"growing"else if(g<1-1e-12)"declining"else"threshold",
       Diagnostics="For a periodic environment the correct threshold is the spectral radius of the ordered product, not the mean of timestep-specific eigenvalues.")
}

.ina_tm_single_carrier_pinf <- function(v,N,beta,C,Transmission,ds){
  U<-length(N);out<-numeric(U)
  if(Transmission=="frequency"){
    den<-as.numeric(crossprod(N,C))
    for(w in seq_len(U))if(den[w]>0&&beta[w]>0&&C[v,w]>0)out[w]<--expm1(-beta[w]*C[v,w]/den[w])
  } else for(w in seq_len(U))if(beta[w]>0&&C[v,w]>0)out[w]<--expm1(-beta[w]*C[v,w]/ds[w])
  pmin(1,pmax(0,out))
}

INApestTransitionPathogenBranching <- function(
    Model=c("SIS","SIR","SEIR"),timesteps,Transition,Nstages,HostPopulationDiseaseStep,
    Beta,RecoveryProb,ProgressionProb=1,PathogenMortalityProb=0,
    ContactMatrix=NULL,StageMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
    ManageProb=0,MortalityProb=0,TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0,
    OutsideNodes=integer(0),OutsideUnits=integer(0)){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission);T<-as.integer(timesteps)
  Nmat<-as.matrix(HostPopulationDiseaseStep);n<-nrow(Nmat);S<-ncol(Nmat);U<-n*S;N<-as.numeric(t(Nmat))
  car<-INApestTransitionPathogenCarrierOperator(Transition,S,n,ManageProb,MortalityProb,
                                                 TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  H<-car$Operator
  beta<-as.numeric(t(.ina_tm_ns_matrix(Beta,n,S,"Beta")));rec<-as.numeric(t(.ina_tm_ns_matrix(RecoveryProb,n,S,"RecoveryProb",TRUE)))
  prog<-as.numeric(t(.ina_tm_ns_matrix(ProgressionProb,n,S,"ProgressionProb",TRUE)));mort<-as.numeric(t(.ina_tm_ns_matrix(PathogenMortalityProb,n,S,"PathogenMortalityProb",TRUE)))
  ds<-as.numeric(t(.ina_tm_ns_matrix(DensityScale,n,S,"DensityScale")))
  Cn<-if(is.null(ContactMatrix))diag(n)else as.matrix(ContactMatrix);Sm<-if(is.null(StageMixing))matrix(1,S,S)else as.matrix(StageMixing);C<-kronecker(Cn,Sm)
  outside<-unique(as.integer(OutsideUnits));if(length(OutsideNodes))for(i in OutsideNodes)outside<-c(outside,(i-1L)*S+seq_len(S));outside<-unique(outside)
  if(any(!outside%in%seq_len(U)))stop("Outside node/unit index invalid")
  integerN<-all(abs(N-round(N))<1e-10)
  transmission_pgf<-function(pinf,Ssus,z){
    if(integerN)prod((1-pinf+pinf*z)^as.integer(round(Ssus)))
    else exp(sum(Ssus*pinf*(z-1)))
  }
  # q = active-lineage extinction-by-horizon; h = no first-passage outside active state by horizon.
  if(Model%in%c("SIS","SIR")){q<-rep(0,U);h<-rep(1,U);if(length(outside))h[outside]<-0;qh<-hh<-matrix(NA,U,T)
    stepfun<-function(z,noescape=FALSE){zn<-numeric(U);for(u in seq_len(U)){
      val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0)){
        pinf<-.ina_tm_single_carrier_pinf(v,N,beta,C,Transmission,ds);Ssus<-pmax(0,N-as.numeric(seq_len(U)==v));trans<-transmission_pgf(pinf,Ssus,z);stay<-1-rec[v]-mort[v];carrier<-(1-stay)+stay*z[v];val<-val+H[v,u]*carrier*trans}
      zn[u]<-val};pmin(1,pmax(0,zn))}
    for(tt in seq_len(T)){q<-stepfun(q);h<-stepfun(h,TRUE);if(length(outside))h[outside]<-0;qh[,tt]<-q;hh[,tt]<-h}
    meanop<-matrix(0,U,U);for(u in seq_len(U))for(v in which(H[,u]>0)){
      pinf<-.ina_tm_single_carrier_pinf(v,N,beta,C,Transmission,ds);Ssus<-pmax(0,N-as.numeric(seq_len(U)==v));meanop[,u]<-meanop[,u]+H[v,u]*(Ssus*pinf);meanop[v,u]<-meanop[v,u]+H[v,u]*(1-rec[v]-mort[v])}
  } else {
    q<-rep(0,2*U);h<-rep(1,2*U);if(length(outside))h[c(outside,U+outside)]<-0;qh<-hh<-matrix(NA,2*U,T)
    stepfun<-function(z){zn<-numeric(2*U)
      # E source: demographic carrier then E stay/progress; no transmission.
      for(u in seq_len(U)){val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0))val<-val+H[v,u]*((1-prog[v])*z[v]+prog[v]*z[U+v]);zn[u]<-val}
      # I source: demographic carrier, transmission creates E, carrier may stay I.
      for(u in seq_len(U)){val<-max(0,1-sum(H[,u]));for(v in which(H[,u]>0)){
        pinf<-.ina_tm_single_carrier_pinf(v,N,beta,C,Transmission,ds);Ssus<-pmax(0,N-as.numeric(seq_len(U)==v));trans<-transmission_pgf(pinf,Ssus,z[seq_len(U)]);stay<-1-rec[v]-mort[v];val<-val+H[v,u]*((1-stay)+stay*z[U+v])*trans};zn[U+u]<-val}
      pmin(1,pmax(0,zn))}
    for(tt in seq_len(T)){q<-stepfun(q);h<-stepfun(h);if(length(outside))h[c(outside,U+outside)]<-0;qh[,tt]<-q;hh[,tt]<-h}
    meanop<-matrix(0,2*U,2*U)
    for(u in seq_len(U))for(v in which(H[,u]>0)){
      meanop[v,u]<-meanop[v,u]+H[v,u]*(1-prog[v]);meanop[U+v,u]<-meanop[U+v,u]+H[v,u]*prog[v]
      pinf<-.ina_tm_single_carrier_pinf(v,N,beta,C,Transmission,ds);Ssus<-pmax(0,N-as.numeric(seq_len(U)==v));meanop[seq_len(U),U+u]<-meanop[seq_len(U),U+u]+H[v,u]*(Ssus*pinf);meanop[U+v,U+u]<-meanop[U+v,U+u]+H[v,u]*(1-rec[v]-mort[v])
    }
  }
  rho<-max(Mod(eigen(meanop,only.values=TRUE)$values))
  inside<-setdiff(seq_len(U),outside)
  if(Model=="SEIR") inside_types<-c(inside,U+inside) else inside_types<-inside
  esc<-1-h
  list(Model="INApestMetaTransitionMatrix",PathogenModel=Model,Exact=FALSE,
       OneCarrierMeanOperator=meanop,OneCarrierBranchingMultiplier=rho,
       ExtinctionByHorizon=q,EscapeByHorizon=esc,
       InsideStartingTypes=inside_types,
       MaxEscapeByHorizonFromInside=if(length(inside_types))max(esc[inside_types])else NA_real_,
       ExtinctionHistory=qh,EscapeHistory=1-hh,OutsideUnits=outside,UnitMap=car$UnitMap,
       Diagnostics=c(
         "Multitype branching approximation is exact for descendants of a single lineage until different pathogen lineages compete for the same susceptible hosts.",
         "The one-carrier offspring distribution uses the stochastic engine's exact exponential infection probability, not only its infinitesimal Jacobian.",
         if(integerN) "Integer host backgrounds use the exact finite-count binomial offspring PGF." else "Fractional analytical host backgrounds use a mean-matched Poisson offspring PGF; this avoids treating a fractional host count as a binomial exponent.",
         "Escape is defined at timestep boundaries by an active pathogen state outside containment; transient carrier movement followed by recovery/death within the same timestep is not counted unless it leaves an active descendant outside."))
}

INApestTransitionHostMeanOperator <- function(
    Transition,Nstages,SDDprob,LDDprob=0,LDDrate=0,
    EnvEstabProb=1,PropaguleEstablishment=1,
    ManageProb=0,MortalityProb=0,SpreadReduction=0,FecundityReduction=0,
    DispersalDensityFactor=0,K=1,SeedbankK=1,
    TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0){
  SDD<-as.matrix(SDDprob);n<-nrow(SDD);S<-as.integer(Nstages)
  z<-transition_components(Transition,S,SDD,LDDprob,LDDrate,EnvEstabProb,
                           PropaguleEstablishment,ManageProb,MortalityProb,
                           SpreadReduction,DispersalDensityFactor,K,SeedbankK,
                           FecundityReduction)
  car<-INApestTransitionPathogenCarrierOperator(Transition,S,n,ManageProb,MortalityProb,
                                                 TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  H<-car$Operator;R<-matrix(0,n*S,n*S);idx<-function(i,s)(i-1L)*S+s
  F<-if(!is.null(z$fecundity_reduction))z$fecundity_reduction else matrix(0,n,S)
  for(i in seq_len(n)){
    Ai<-z$A[[i]];a<-z$adoption[i];g<-z$spread_reduction[i]
    for(k in 2:S){
      fec<-Ai[1,k];if(fec<=0)next
      m<-z$mortality[i,k];fr<-F[i,k]
      for(j in seq_len(n)){
        nat<-(1-z$LDDrate)*z$sdd_enabled[i]*z$SDD[i,j]
        hum0<-z$LDDrate*z$LDD[i,j]
        hum1<-z$LDDrate*(1-g)*z$LDD[i,j]
        w<-fec*z$recruit_success[j]*((1-a)*(nat+hum0)+a*(1-m)*(1-fr)*(nat+hum1))
        R[idx(j,1L),idx(i,k)]<-R[idx(j,1L),idx(i,k)]+w
      }
    }
  }
  G<-H+R
  dimnames(R)<-dimnames(H);dimnames(G)<-dimnames(H)
  list(Operator=G,CarrierOperator=H,RecruitmentOperator=R,UnitMap=car$UnitMap,
       Diagnostics=c(car$Diagnostics,
         "Host mean operator is the uncrowded first-moment operator: existing hosts are carried by H and reproduction contributes susceptible stage-1 recruits through R.",
         "Management mortality, fecundity reduction and LDD spread reduction are integrated over the node-level management adoption event in the reproductive term.",
         "Because recruitment and crowding are nonlinear in finite populations, this operator is intended for disease-free mean-background construction in the low-density/large-host regime."))
}

INApestTransitionPathogenGrowthFromHostStart <- function(
    Model=c("SIS","SIR","SEIR"),Transition,Nstages,InitialHostPopulation,
    SDDprob,LDDprob=0,LDDrate=0,EnvEstabProb=1,PropaguleEstablishment=1,
    ManageProb=0,MortalityProb=0,SpreadReduction=0,FecundityReduction=0,
    DispersalDensityFactor=0,K=1,SeedbankK=1,
    TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0,
    Beta,RecoveryProb,ProgressionProb=1,PathogenMortalityProb=0,
    ContactMatrix=NULL,StageMixing=NULL,Transmission=c("frequency","density"),DensityScale=1){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission)
  X<-as.matrix(InitialHostPopulation);n<-nrow(X);S<-ncol(X);if(S!=Nstages)stop("InitialHostPopulation columns must equal Nstages")
  hm<-INApestTransitionHostMeanOperator(Transition,S,SDDprob,LDDprob,LDDrate,EnvEstabProb,
    PropaguleEstablishment,ManageProb,MortalityProb,SpreadReduction,FecundityReduction,
    DispersalDensityFactor,K,SeedbankK,TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  x0<-as.numeric(t(X));x1<-as.numeric(hm$Operator%*%x0);N1<-matrix(x1,nrow=n,ncol=S,byrow=TRUE)
  pg<-INApestTransitionPathogenGrowthOperator(Model,Transition,S,N1,Beta,RecoveryProb,ProgressionProb,
    PathogenMortalityProb,ContactMatrix,StageMixing,Transmission,DensityScale,ManageProb,MortalityProb,
    TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  pg$InitialHostPopulation<-X;pg$DerivedHostPopulationDiseaseStep<-N1;pg$HostMeanOperator<-hm$Operator
  pg$Diagnostics<-c(pg$Diagnostics,
    "HostPopulationDiseaseStep was derived from the same-timestep uncrowded disease-free host mean operator. For strongly crowded or highly stochastic host populations, provide/estimate the disease-step host background directly instead.")
  pg
}

INApestTransitionPathogenExactDiseaseStepFixedN <- function(
    Model=c("SIS","SIR","SEIR"),Ntimesteps,HostPopulationDiseaseStep,
    InitialInfected,InitialExposed=0,InitialRecovered=0,
    Beta,RecoveryProb,ProgressionProb=1,ImmunityLossProb=0,
    IntroductionProb=0,IntroductionNumber=1,
    ContactMatrix=NULL,StageMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
    OutsideNodes=integer(0),OutsideUnits=integer(0),MaxStates=50000L,ReturnOperator=FALSE){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission)
  Nmat<-as.matrix(HostPopulationDiseaseStep);n<-nrow(Nmat);S<-ncol(Nmat);U<-n*S
  if(any(Nmat<0|Nmat!=floor(Nmat)))stop("HostPopulationDiseaseStep must contain non-negative whole numbers for exact finite-state analysis")
  N<-as.integer(as.numeric(t(Nmat)))
  rv<-function(x,name,integer=FALSE){z<-.ina_tm_ns_matrix(x,n,S,name);z<-as.numeric(t(z));if(integer)z<-as.integer(z);z}
  I0<-rv(InitialInfected,"InitialInfected",TRUE);E0<-rv(InitialExposed,"InitialExposed",TRUE);R0<-rv(InitialRecovered,"InitialRecovered",TRUE)
  if(any(I0+E0+R0>N))stop("Initial pathogen-state counts exceed host counts")
  beta<-rv(Beta,"Beta");rec<-rv(RecoveryProb,"RecoveryProb");prog<-rv(ProgressionProb,"ProgressionProb");wan<-rv(ImmunityLossProb,"ImmunityLossProb")
  ip<-rv(IntroductionProb,"IntroductionProb");inum<-rv(IntroductionNumber,"IntroductionNumber",TRUE);ds<-rv(DensityScale,"DensityScale")
  Cn<-if(is.null(ContactMatrix))diag(n)else as.matrix(ContactMatrix);Sm<-if(is.null(StageMixing))matrix(1,S,S)else as.matrix(StageMixing);C<-kronecker(Cn,Sm)
  outside<-unique(as.integer(OutsideUnits));if(length(OutsideNodes))for(i in OutsideNodes)outside<-c(outside,(i-1L)*S+seq_len(S));outside<-unique(outside)
  if(Model=="SIS")ans<-INApestMetaSISExactNetworkFixedN(Ntimesteps,N,I0,beta,rec,C,ip,inum,Transmission,ds,OutsideNodes=outside,ReturnOperator=ReturnOperator)
  else ans<-INApestMetaCompartmentExactNetworkFixedN(Model,Ntimesteps,N,I0,E0,R0,beta,rec,prog,wan,ip,inum,C,Transmission,ds,OutsideNodes=outside,MaxStates=MaxStates,ReturnOperator=ReturnOperator)
  ans$Model<-"INApestMetaTransitionMatrix";ans$HostArchitecture<-"conditional fixed node x demographic-stage abundance";ans$PathogenModel<-Model
  ans$HostPopulationDiseaseStep<-Nmat;ans$OutsideUnits<-outside
  ans$Diagnostics<-c(ans$Diagnostics,
    "This exact finite-state branch conditions on a fixed node x stage host background and validates the pathogen step itself.",
    "It deliberately does not apply demographic stage transitions between pathogen steps; full Transition-Matrix invasion growth/escape uses the carrier-composed operators.")
  class(ans)<-c("INApestTransitionPathogenExactDiseaseStepFixedN","list");ans
}

INApestTransitionPathogenNextGeneration <- function(GrowthObject){
  x<-GrowthObject
  if(!is.list(x)||is.null(x$CarrierOperator)||is.null(x$TransmissionBlock)||is.null(x$PathogenModel))stop("Supply an INApestTransitionPathogenGrowthOperator result")
  H<-x$CarrierOperator;B<-x$TransmissionBlock;U<-nrow(H);model<-x$PathogenModel
  # Recover disease-only persistence/progression from the already composed
  # disease operator. This avoids re-resolving parameters and preserves exact
  # orientation.
  D<-x$DiseaseOperator
  if(model%in%c("SIS","SIR")){
    # D = D_I + B
    DI<-D-B;T<-DI%*%H;F<-B%*%H
  } else {
    # disease block rows/cols = E,I. New-infection contribution is only I -> E.
    Z<-matrix(0,U,U);Fdis<-rbind(cbind(Z,B),cbind(Z,Z));Tdis<-D-Fdis
    H2<-rbind(cbind(H,Z),cbind(Z,H));T<-Tdis%*%H2;F<-Fdis%*%H2
  }
  rhoT<-max(Mod(eigen(T,only.values=TRUE)$values))
  if(rhoT>=1-1e-12){
    K<-matrix(NA_real_,nrow(T),ncol(T));R0<-Inf
    note<-"The no-new-infection active-state process is non-transient (spectral radius >= 1), so the expected lifetime next-generation sum diverges."
  } else {
    K<-F%*%solve(diag(nrow(T))-T);R0<-max(Mod(eigen(K,only.values=TRUE)$values))
    note<-"K = F (I-T)^(-1) sums new active infections produced over the full remaining demographic/pathogen lifetime of an initial active-state cohort."
  }
  list(Model="INApestMetaTransitionMatrix",PathogenModel=model,TransitionWithoutNewInfection=T,
       NewInfectionOperator=F,NextGenerationOperator=K,R0=R0,
       PerTimestepLambda=x$IntrinsicRarePathogenMultiplier,
       ThresholdAgreement=if(is.finite(R0))sign(R0-1)==sign(x$IntrinsicRarePathogenMultiplier-1)else x$IntrinsicRarePathogenMultiplier>=1,
       Diagnostics=c(note,
         "R0 and lambda answer different questions: R0 is lifetime secondary active infection production, whereas lambda is the asymptotic per-timestep multiplier.",
         "Under the usual non-negative transient-state conditions their invasion thresholds agree even though their numerical values differ."))
}

INApestTransitionPathogenAnalytical <- function(
    Ntimesteps=10,Transition,Nstages,Pathogen,
    InformationMode=c("none","all_informed"),
    HostPopulationDiseaseStep=NULL,InitialHostPopulation=NULL,
    SDDprob=NULL,LDDprob=0,LDDrate=0,EnvEstabProb=1,PropaguleEstablishment=1,
    K=1,SeedbankK=1,ManageProb=0,MortalityProb=0,SpreadReduction=0,FecundityReduction=0,
    DispersalDensityFactor=0,TransitionSDDprob=NULL,TransitionLDDprob=NULL,TransitionLDDrate=0,
    StageMixing=NULL,InitialInfected=0,InitialExposed=0,InitialRecovered=0,
    OutsideNodes=integer(0),OutsideUnits=integer(0),Exact=TRUE,ExactMaxStates=50000L,
    ReturnOperators=FALSE){
  if(!inherits(Pathogen,"INApestPathogen"))stop("Pathogen must be created by INApestPathogen()")
  InformationMode<-match.arg(InformationMode)
  model<-as.character(Pathogen$Model)[1L];if(!model%in%c("SIS","SIR","SEIR"))stop("Transition-Matrix analytical pathogen methods currently support SIS, SIR and SEIR")
  ManageProbEff<-if(InformationMode=="all_informed")ManageProb else 0
  if(is.null(HostPopulationDiseaseStep)){
    if(is.null(InitialHostPopulation)||is.null(SDDprob))stop("Supply HostPopulationDiseaseStep, or InitialHostPopulation plus SDDprob to derive the disease-step host mean background")
    bg<-INApestTransitionPathogenGrowthFromHostStart(model,Transition,Nstages,InitialHostPopulation,SDDprob,LDDprob,LDDrate,
      EnvEstabProb,PropaguleEstablishment,ManageProbEff,MortalityProb,SpreadReduction,FecundityReduction,
      DispersalDensityFactor,K,SeedbankK,TransitionSDDprob,TransitionLDDprob,TransitionLDDrate,
      Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$PathogenMortalityProb,
      Pathogen$ContactMatrix,StageMixing,Pathogen$Transmission,Pathogen$DensityScale)
    HostPopulationDiseaseStep<-bg$DerivedHostPopulationDiseaseStep;growth<-bg
  } else {
    HostPopulationDiseaseStep<-as.matrix(HostPopulationDiseaseStep)
    growth<-INApestTransitionPathogenGrowthOperator(model,Transition,Nstages,HostPopulationDiseaseStep,
      Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$PathogenMortalityProb,
      Pathogen$ContactMatrix,StageMixing,Pathogen$Transmission,Pathogen$DensityScale,
      ManageProbEff,MortalityProb,TransitionSDDprob,TransitionLDDprob,TransitionLDDrate)
  }
  ng<-INApestTransitionPathogenNextGeneration(growth)
  br<-INApestTransitionPathogenBranching(model,Ntimesteps,Transition,Nstages,HostPopulationDiseaseStep,
    Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$PathogenMortalityProb,
    Pathogen$ContactMatrix,StageMixing,Pathogen$Transmission,Pathogen$DensityScale,
    ManageProbEff,MortalityProb,TransitionSDDprob,TransitionLDDprob,TransitionLDDrate,OutsideNodes,OutsideUnits)
  exact<-NULL
  exact_reason<-NULL
  if(isTRUE(Exact)){
    if(any(as.numeric(Pathogen$PathogenMortalityProb)!=0,na.rm=TRUE)) exact_reason<-"Conditional fixed-host exact disease-step analysis requires PathogenMortalityProb = 0 because pathogen deaths change host abundance."
    else if(any(HostPopulationDiseaseStep!=floor(HostPopulationDiseaseStep))) exact_reason<-"Conditional exact disease-step analysis requires whole-number host counts."
    else {
      exact<-tryCatch(INApestTransitionPathogenExactDiseaseStepFixedN(model,Ntimesteps,HostPopulationDiseaseStep,
        InitialInfected,InitialExposed,InitialRecovered,Pathogen$Beta,Pathogen$RecoveryProb,
        Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,
        Pathogen$ContactMatrix,StageMixing,Pathogen$Transmission,Pathogen$DensityScale,OutsideNodes,OutsideUnits,
        ExactMaxStates,ReturnOperators),error=function(e){exact_reason<<-conditionMessage(e);NULL})
    }
  }
  out<-list(Model="INApestMetaTransitionMatrix",PathogenModel=model,InformationMode=InformationMode,
    HostPopulationDiseaseStep=HostPopulationDiseaseStep,Growth=growth,NextGeneration=ng,Branching=br,
    ExactConditionalDiseaseStep=exact,ExactConditionalUnavailableReason=exact_reason,
    Diagnostics=c(
      "Transition-Matrix pathogen analysis separates demographic carrier transport from the pathogen state process.",
      if(InformationMode=="none") "Management-dependent vital rates are evaluated on an uninformed background (no management adoption)." else "Management-dependent vital rates are evaluated on an all-informed background; ManageProb remains the conditional adoption probability.",
      "Endogenous host/pathogen detection, information transfer and programmed-stop feedback are shared node-level states and are not closed by an independent-lineage branching process; use none/all-informed analytical envelopes and stochastic simulation for that feedback.",
      "The scalable headline invasion metrics are the full-timestep rare-pathogen lambda, lifetime next-generation R0, and multitype branching extinction/escape.",
      "The finite-state exact branch conditions on fixed node x stage host abundance; it validates disease mechanics but does not replace the carrier-composed growth/branching solution for evolving demographic stages."))
  class(out)<-c("INApestTransitionPathogenAnalytical","list");out
}

print.INApestTransitionPathogenAnalytical <- function(x,...){
  cat("INApest Transition-Matrix pathogen analytical result\n")
  cat("  Pathogen model:",x$PathogenModel,"\n")
  cat("  Per-timestep lambda:",format(x$Growth$IntrinsicRarePathogenMultiplier,digits=7),"\n")
  cat("  Lifetime R0:",format(x$NextGeneration$R0,digits=7),"\n")
  if(length(x$Branching$OutsideUnits))cat("  Branching escape by horizon (max inside starting type):",format(x$Branching$MaxEscapeByHorizonFromInside,digits=7),"\n")
  if(!is.null(x$ExactConditionalDiseaseStep))cat("  Conditional finite-state disease-step solution: available\n")
  else if(!is.null(x$ExactConditionalUnavailableReason))cat("  Conditional finite-state disease-step solution:",x$ExactConditionalUnavailableReason,"\n")
  invisible(x)
}
