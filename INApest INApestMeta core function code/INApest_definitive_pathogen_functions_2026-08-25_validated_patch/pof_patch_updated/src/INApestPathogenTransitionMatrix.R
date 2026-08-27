###############################################################################
### Generic pathogen x demographic-stage support for INApestMetaTransitionMatrix
###############################################################################

.iptm_take <- function(counts, take) {
  counts <- as.integer(counts); take <- as.integer(take); total <- sum(counts)
  if (take <= 0L) return(integer(length(counts)))
  if (take >= total) return(counts)
  out <- integer(length(counts)); rem_take <- take; rem_total <- total
  if (length(counts) > 1L) for (i in seq_len(length(counts)-1L)) {
    out[i] <- rhyper(1L, counts[i], rem_total-counts[i], rem_take)
    rem_take <- rem_take-out[i]; rem_total <- rem_total-counts[i]
  }
  out[length(counts)] <- rem_take
  out
}

.iptm_flow <- function(source_counts, P) {
  n <- length(source_counts); flows <- matrix(0L,n,n); exported <- integer(n)
  for(i in which(source_counts>0)) {
    pr <- c(P[i,], max(0,1-sum(P[i,])))
    z <- as.integer(rmultinom(1L,source_counts[i],pr)); flows[i,] <- z[seq_len(n)]; exported[i] <- z[n+1L]
  }
  list(flows=flows, exported=exported)
}

.iptm_reconcile <- function(state, N) {
  if (length(dim(state)) != 3L) stop("PathogenStageState must be nodes x demographic stages x pathogen states")
  if (!all(dim(state)[1:2] == dim(N))) stop("PathogenStageState and N dimensions differ")
  P <- dim(state)[3]
  for (i in seq_len(nrow(N))) for (s in seq_len(ncol(N))) {
    z <- as.integer(state[i,s,]); target <- as.integer(N[i,s]); total <- sum(z)
    if (target < total) state[i,s,] <- .iptm_take(z, target)
    else if (target > total) state[i,s,1L] <- state[i,s,1L] + target-total
  }
  storage.mode(state) <- "integer"; state
}

INApestPathogenStageState <- function(N, Pathogen, InitialState = NULL, Ntimesteps = 1L) {
  if (!inherits(Pathogen, "INApestPathogen")) stop("Pathogen must be created by INApestPathogen()")
  N <- as.matrix(N); D <- ncol(N); P <- length(Pathogen$States); n <- nrow(N)
  if (!is.null(InitialState)) {
    if (!identical(dim(InitialState), c(n,D,P))) stop("InitialState must be nodes x demographic stages x pathogen states")
    if (any(InitialState < 0) || any(InitialState != floor(InitialState))) stop("InitialState must contain non-negative whole numbers")
    if (any(apply(InitialState,c(1,2),sum) != N)) stop("InitialState must sum to InitialPopulation within every node x demographic stage")
    dimnames(InitialState)[[3]] <- Pathogen$States
    return(InitialState)
  }
  out <- array(0L, c(n,D,P), dimnames=list(NULL,NULL,Pathogen$States)); out[, , "S"] <- N
  # Convenience seeding: node-level initial counts are assigned randomly across
  # demographic stages in proportion to the hosts actually present.
  ctx <- list(n_nodes=n, Ntimesteps=as.integer(Ntimesteps))
  tmp <- Pathogen$Engine$Initial(rowSums(N), ctx)
  for (i in seq_len(n)) {
    remaining <- as.integer(N[i,])
    for (st in setdiff(Pathogen$States, "S")) {
      k <- as.integer(tmp[i,st]); if (k <= 0L) next
      alloc <- .iptm_take(remaining, k); out[i,,st] <- alloc; out[i,,"S"] <- out[i,,"S"]-alloc; remaining <- remaining-alloc
    }
  }
  out
}

.iptm_resolve <- function(x, timestep, n_nodes, D, Ntimesteps, name, probability=FALSE) {
  if (is.function(x)) {
    fm <- names(formals(x)); a <- list(timestep=timestep,n_nodes=n_nodes,Nstages=D,Ntimesteps=Ntimesteps)
    if (!is.null(fm) && !"..." %in% fm) a <- a[intersect(names(a),fm)]
    x <- do.call(x,a)
  }
  d <- dim(x)
  if (!is.null(d)) {
    if (length(d)==2L && all(d==c(n_nodes,D))) z <- as.numeric(x)
    else if (length(d)==3L && all(d==c(n_nodes,D,Ntimesteps))) z <- as.numeric(x[,,timestep])
    else stop(name," must be scalar, stage vector, node vector, nodes x stages, or nodes x stages x timesteps")
  } else {
    x <- as.numeric(x)
    if (length(x)==1L) z <- rep(x,n_nodes*D)
    else if (length(x)==D && D!=n_nodes) z <- rep(x, each=n_nodes)
    else if (length(x)==n_nodes && D!=n_nodes) z <- rep(x,D)
    else if (length(x)==n_nodes*D) z <- x
    else if (length(x)==Ntimesteps && Ntimesteps!=D && Ntimesteps!=n_nodes) z <- rep(x[timestep],n_nodes*D)
    else stop(name," has an ambiguous or unsupported shape; use nodes x stages explicitly")
  }
  z <- matrix(z,nrow=n_nodes,ncol=D)
  if (probability && (any(!is.finite(z)) || any(z<0|z>1))) stop(name," must resolve to [0,1]")
  z
}

.iptm_node_contact <- function(Pathogen, timestep, n_nodes, Ntimesteps) {
  x <- Pathogen$ContactMatrix
  if (is.null(x)) return(diag(n_nodes))
  if (is.function(x)) {
    fm <- names(formals(x)); a <- list(timestep=timestep,n_nodes=n_nodes,Ntimesteps=Ntimesteps)
    if (!is.null(fm) && !"..." %in% fm) a <- a[intersect(names(a),fm)]
    x <- do.call(x,a)
  }
  d <- dim(x)
  if (length(d)==2L && all(d==c(n_nodes,n_nodes))) out <- x
  else if (length(d)==3L && all(d==c(n_nodes,n_nodes,Ntimesteps))) out <- x[,,timestep]
  else stop("For transition-matrix pathogen models ContactMatrix must be nodes x nodes or nodes x nodes x Ntimesteps")
  if (any(!is.finite(out)) || any(out<0)) stop("ContactMatrix entries must be finite and non-negative")
  as.matrix(out)
}

.iptm_pathogen_step <- function(state, Pathogen, timestep, Ntimesteps, StageMixing=NULL) {
  n_nodes <- dim(state)[1]; D <- dim(state)[2]; states <- Pathogen$States; P <- length(states)
  if (is.null(StageMixing)) StageMixing <- matrix(1,D,D)
  if (!is.matrix(StageMixing) || !all(dim(StageMixing)==c(D,D)) || any(StageMixing<0) || any(!is.finite(StageMixing))) stop("StageMixing must be a finite non-negative Nstages x Nstages matrix")
  # node-major, stage-fast unit matrix
  M <- matrix(0L,n_nodes*D,P,dimnames=list(NULL,states))
  for (i in seq_len(n_nodes)) for (s in seq_len(D)) M[(i-1L)*D+s,] <- state[i,s,]
  getp <- function(nm,prob=FALSE) as.numeric(t(.iptm_resolve(Pathogen[[nm]],timestep,n_nodes,D,Ntimesteps,nm,prob)))
  beta <- getp("Beta"); rec <- getp("RecoveryProb",TRUE); mort <- getp("PathogenMortalityProb",TRUE)
  prog <- getp("ProgressionProb",TRUE); wan <- getp("ImmunityLossProb",TRUE); ip <- getp("IntroductionProb",TRUE)
  inum <- as.integer(getp("IntroductionNumber")); ds <- getp("DensityScale")
  if (any(rec+mort>1+1e-12)) stop("RecoveryProb + PathogenMortalityProb must not exceed 1")
  if (any(beta<0|!is.finite(beta)) || any(ds<=0|!is.finite(ds))) stop("Beta must be non-negative and DensityScale positive")
  C <- kronecker(.iptm_node_contact(Pathogen,timestep,n_nodes,Ntimesteps), StageMixing)
  S0 <- M[,"S"]; I0 <- M[,"I"]; E0 <- if("E"%in%states) M[,"E"] else integer(nrow(M)); R0 <- if("R"%in%states) M[,"R"] else integer(nrow(M)); live <- rowSums(M)
  infp <- as.numeric(crossprod(I0,C))
  if (Pathogen$Transmission=="frequency") { den <- as.numeric(crossprod(live,C)); force <- beta*ifelse(den>0,infp/den,0) } else force <- beta*infp/ds
  pinf <- pmin(1,pmax(0,-expm1(-force))); newinf <- rbinom(length(S0),S0,pinf)
  introduced <- integer(length(S0)); ev <- rbinom(length(S0),1L,ip); ii <- which(ev>0 & S0-newinf>0); if(length(ii)) introduced[ii] <- pmin(inum[ii],S0[ii]-newinf[ii])
  progress <- if("E"%in%states) rbinom(length(E0),E0,prog) else integer(length(E0)); stay <- recover <- deaths <- integer(length(I0))
  for (i in which(I0>0)) { z<-as.numeric(rmultinom(1,I0[i],c(1-rec[i]-mort[i],rec[i],mort[i]))); stay[i]<-z[1]; recover[i]<-z[2]; deaths[i]<-z[3] }
  lose <- if("R"%in%states) rbinom(length(R0),R0,wan) else integer(length(R0))
  M[,"S"] <- S0-newinf-introduced+lose
  if("E"%in%states) { M[,"E"]<-E0-progress+newinf+introduced; M[,"I"]<-stay+progress } else M[,"I"]<-stay+newinf+introduced
  if("R"%in%states) M[,"R"]<-R0-lose+recover else M[,"S"]<-M[,"S"]+recover
  for (i in seq_len(n_nodes)) for(s in seq_len(D)) state[i,s,] <- M[(i-1L)*D+s,]
  list(State=state, Deaths=matrix(deaths,nrow=n_nodes,ncol=D), NewInfections=matrix(newinf,nrow=n_nodes,ncol=D), Introduced=matrix(introduced,nrow=n_nodes,ncol=D))
}

local.dynamics.transition.matrix.pathogen <- function(
    nodetransition, weights, sddprob, nodeenvestabprob, n0,
    lddprob=NA, lddrate=0, nodeK, node.seedbankK,
    nodepropaguleestablishment, nodespreadreduction, managing,
    MaxInteger=.Machine$integer.max, nodefecundityreduction=0,
    transition_sddprob=NULL, transition_lddprob=NULL, transition_lddrate=0,
    BlockedTransitionMortality=0, DispersalDensityFactor=0,
    pathogen_state, Pathogen, timestep, Ntimesteps, StageMixing=NULL, ...)
{
  if (!inherits(Pathogen,"INApestPathogen")) stop("Pathogen must be created by INApestPathogen()")
  N <- as.matrix(n0); n_nodes<-nrow(N); D<-ncol(N); states<-Pathogen$States; P<-length(states)
  state <- .iptm_reconcile(pathogen_state,N)
  A_list <- if(is.list(nodetransition)) nodetransition else replicate(n_nodes,nodetransition,simplify=FALSE)
  if(length(A_list)!=n_nodes || any(vapply(A_list,function(A)!all(dim(A)==c(D,D)),logical(1)))) stop("Demographic Transition must be Nstages x Nstages (or one such matrix per node)")
  w <- if(is.matrix(weights)) weights else matrix(rep(weights,each=n_nodes),nrow=n_nodes,ncol=D)
  if(!all(dim(w)==c(n_nodes,D))) stop("weights must describe demographic stages")
  nodeK<-rep_len(nodeK,n_nodes); node.seedbankK<-rep_len(node.seedbankK,n_nodes); env<-rep_len(nodeenvestabprob,n_nodes); foot<-rep_len(nodepropaguleestablishment,n_nodes)
  fr <- .iptm_resolve(nodefecundityreduction,1L,n_nodes,D,1L,"nodefecundityreduction",TRUE); efffr <- 1-fr*managing
  bm <- if(length(BlockedTransitionMortality)==1L) matrix(BlockedTransitionMortality,n_nodes,max(1,D-1L)) else matrix(rep_len(BlockedTransitionMortality,n_nodes*max(1,D-1L)),n_nodes)
  norm_tm <- function(x,name) {
    if(is.null(x)) return(rep(list(NULL),max(0,D-1L)))
    z <- if(is.list(x)) x else rep(list(x),D-1L)
    if(length(z)!=D-1L) stop(name," must be a matrix or list length Nstages-1")
    for(k in seq_along(z)) if(!is.null(z[[k]])) { if(!is.matrix(z[[k]]) || !all(dim(z[[k]])==c(n_nodes,n_nodes))) stop(name," matrices must be nodes x nodes"); if(any(z[[k]]<0)||any(rowSums(z[[k]])>1+1e-10)) stop(name," rows must contain non-negative probabilities summing to <=1") }
    z
  }
  tsdd <- norm_tm(transition_sddprob,"transition_sddprob"); tldd <- norm_tm(transition_lddprob,"transition_lddprob")
  trrate <- rep_len(transition_lddrate,max(1,D-1L)); if(any(trrate<0|trrate>1)) stop("transition_lddrate must be in [0,1]")
  # Fecundity from all pathogen states; offspring are susceptible.
  fecmeans <- numeric(n_nodes); mothers <- numeric(n_nodes)
  for(i in seq_len(n_nodes)) { A<-A_list[[i]]; f<-A[1,-1]; fecmeans[i]<-sum(f*N[i,-1]*efffr[i,-1]); mothers[i]<-sum(N[i,which(c(FALSE,f>0))]) }
  propagules <- rpois(n_nodes,pmax(0,fecmeans))
  if(!is.na(DispersalDensityFactor) && DispersalDensityFactor!=0 && any(propagules>0)) {
    if(!is.finite(DispersalDensityFactor) || DispersalDensityFactor<=0) stop("DispersalDensityFactor must be positive, or 0/NA")
    prepop <- rowSums(N[,-1,drop=FALSE]*w[,-1,drop=FALSE]); rel <- rep(1,n_nodes); ok<-nodeK>0; rel[ok]<-pmin(1,pmax(0,prepop[ok]/nodeK[ok])); attr<-pmax(0,1-rel)^DispersalDensityFactor
    export <- pmax(0,1-rowSums(sddprob)); choice <- sddprob*rep(attr,each=n_nodes); den<-rowSums(choice)+export; mult<-ifelse(den>0,1/den,0); sddprob<-choice*mult
  }
  # Terminal survival by pathogen state.
  for(i in seq_len(n_nodes)) { p<-A_list[[i]][D,D]; for(q in seq_len(P)) state[i,D,q]<-rbinom(1,state[i,D,q],p) }
  # Descending adjacent demographic progression with shared capacity.
  capacity_above <- numeric(n_nodes)
  if(D>1L) for(s in D:2L) {
    # Draw survival and progression candidates by node x pathogen state.
    staymat <- candmat <- matrix(0L,n_nodes,P)
    for(i in seq_len(n_nodes)) {
      A<-A_list[[i]]; tr<-A[s,s-1L]; st<-A[s-1L,s-1L]; surv<-min(1,tr+st); cond<-if(surv>0) tr/surv else 0
      src<-as.integer(state[i,s-1L,]); ns<-vapply(src,function(z) rbinom(1,z,surv),integer(1)); ca<-vapply(ns,function(z) rbinom(1,z,cond),integer(1)); staymat[i,]<-ns-ca; candmat[i,]<-ca
    }
    stage_weight<-w[,s]; target_now<-apply(state[,s,,drop=FALSE],1,sum); total_pop<-capacity_above+target_now*stage_weight; slots<-pmax(0,floor((nodeK-total_pop)/stage_weight)); maxslots<-nodeK/stage_weight; cp<-ifelse(maxslots>0,pmin(1,pmax(0,slots/maxslots)),0)
    moving <- !is.null(tsdd[[s-1L]]) || !is.null(tldd[[s-1L]])
    additions <- matrix(0L,n_nodes,P); blocked <- matrix(0L,n_nodes,P)
    if(!moving) {
      for(i in seq_len(n_nodes)) { accn<-min(slots[i],rbinom(1,sum(candmat[i,]),cp[i])); additions[i,]<-.iptm_take(candmat[i,],accn); blocked[i,]<-candmat[i,]-additions[i,] }
    } else {
      flows <- array(0L,c(n_nodes,n_nodes,P)); exported <- matrix(0L,n_nodes,P)
      for(q in seq_len(P)) {
        cvec<-candmat[,q]; Ps<-tsdd[[s-1L]]; Pl<-tldd[[s-1L]]
        if(!is.null(Ps) && !is.null(Pl)) { nl<-if(trrate[s-1L]<=0) integer(n_nodes) else if(trrate[s-1L]>=1) cvec else rbinom(n_nodes,cvec,trrate[s-1L]); ns<-cvec-nl } else if(!is.null(Ps)) { ns<-cvec; nl<-integer(n_nodes) } else { ns<-integer(n_nodes); nl<-cvec }
        fs<-if(!is.null(Ps)) .iptm_flow(ns,Ps) else list(flows=matrix(0L,n_nodes,n_nodes),exported=integer(n_nodes)); fl<-if(!is.null(Pl)) .iptm_flow(nl,Pl) else list(flows=matrix(0L,n_nodes,n_nodes),exported=integer(n_nodes))
        flows[,,q]<-fs$flows+fl$flows; exported[,q]<-fs$exported+fl$exported
      }
      accepted_source <- matrix(0L,n_nodes,P)
      for(j in seq_len(n_nodes)) {
        arrq<-vapply(seq_len(P),function(q) sum(flows[,j,q]),integer(1)); nt<-sum(arrq); if(nt<=0) next
        na<-min(slots[j],rbinom(1,nt,cp[j])); aq<-.iptm_take(arrq,na); additions[j,]<-aq
        for(q in which(aq>0)) accepted_source[,q] <- accepted_source[,q] + .iptm_take(flows[,j,q],aq[q])
      }
      internal_source <- matrix(0L,n_nodes,P); for(q in seq_len(P)) internal_source[,q]<-rowSums(flows[,,q]); blocked<-pmax(internal_source-accepted_source,0L); storage.mode(blocked)<-"integer"
    }
    for(i in seq_len(n_nodes)) {
      bprob<-bm[i,s-1L]; bsurv<-vapply(blocked[i,],function(z) if(bprob<=0) z else if(bprob>=1) 0L else rbinom(1,z,1-bprob),integer(1)); state[i,s-1L,]<-staymat[i,]+bsurv; state[i,s,]<-state[i,s,]+additions[i,]; capacity_above[i]<-total_pop[i]+sum(additions[i,])*stage_weight[i]
    }
  }
  N <- apply(state,c(1,2),sum)
  # Reproductive dispersal and susceptible recruitment, following the ordinary
  # SDD/LDD semantics but without density-dependent reweighting in this first version.
  Pin<-Qin<-numeric(n_nodes)
  if(any(propagules>0)) {
    ss<-propagules*(1-lddrate); ww<-as.numeric(ss %*% sddprob); if(sum(ww)>0) Pin<-as.numeric(rmultinom(1,floor(sum(ww)),ww))
    if(is.matrix(lddprob)) { red<-pmin(1,pmax(0,rep_len(nodespreadreduction,n_nodes)*managing)); ls<-propagules*lddrate*(1-red); ww<-as.numeric(ls %*% lddprob); if(sum(ww)>0) Qin<-as.numeric(rmultinom(1,floor(sum(ww)),ww)) }
  }
  slots<-pmax(0,floor(node.seedbankK-N[,1])); Pin<-pmax(0,floor(Pin)); Qin<-pmax(0,floor(Qin)); env<-pmin(1,pmax(0,env))
  unrestricted<-all(is.finite(foot)&foot>=1)
  if(unrestricted) accessible<-as.numeric(Pin>0) else { area<-foot*nodeK; press<-numeric(n_nodes); if(any(Pin>0)&&any(mothers>0)) { num<-as.numeric((mothers*area)%*%sddprob); ok<-nodeK>0; press[ok]<-num[ok]/nodeK[ok] }; accessible<-pmin(1,pmax(0,-expm1(-press))) }
  nslots<-floor(slots*accessible); nslots<-ifelse(Pin>0&accessible>0,pmax(1,nslots),0); nslots<-pmin(nslots,slots); oslots<-pmax(0,slots-nslots)
  lp<-ifelse(nslots>0,env*Pin/pmax(nslots,1),0); lq<-ifelse(slots>0,env*Qin/pmax(slots,1),0); recs<-rbinom(n_nodes,nslots,-expm1(-(lp+lq)))+rbinom(n_nodes,oslots,-expm1(-lq)); recs<-pmin(recs,Pin+Qin)
  state[,1,"S"] <- state[,1,"S"]+recs
  # Pathogen update occurs after births, so new hosts enter S and can be exposed in this timestep.
  ps <- .iptm_pathogen_step(state,Pathogen,timestep,Ntimesteps,StageMixing); state<-ps$State; N<-apply(state,c(1,2),sum)
  list(N=N, PathogenState=state, PathogenDeaths=ps$Deaths, NewInfections=ps$NewInfections, PathogenIntroduced=ps$Introduced)
}
