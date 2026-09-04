###############################################################################
### INApestVertebrateNodePathogenSupport -- pathogen transport for Vertebrate Node
###
### Supplies the pathogen-aware local-dynamics adapter used by the frozen
### Vertebrate Node pathogen engine. It carries node x stage x pathogen-state
### counts through host survival, progression, movement, reproduction and
### external host arrivals, then applies pathogen dynamics and detection.
###
### This file is support code, not a separate user-facing simulation architecture.
###############################################################################

# Carry pathogen-state counts through Vertebrate Node host dynamics.
.iv_node_pathogen_transport <- function(
    nodetransition, weights, sddprob, nodeenvestabprob, n0,
    lddprob=NA, lddrate=0, nodeK, node.seedbankK,
    nodepropaguleestablishment, nodespreadreduction, managing,
    MaxInteger=.Machine$integer.max, nodefecundityreduction=0,
    nodecontrolfecundityreduction=0, nodebirthmean=NULL, nodebirthmothers=NULL,
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
  if(any(!is.finite(w)) || any(w <= 0)) stop("All stage weights must be finite and greater than zero")
  nodeK<-rep_len(nodeK,n_nodes); node.seedbankK<-rep_len(node.seedbankK,n_nodes); env<-rep_len(nodeenvestabprob,n_nodes); foot<-rep_len(nodepropaguleestablishment,n_nodes)
  fr <- .iptm_resolve(nodefecundityreduction,1L,n_nodes,D,1L,"nodefecundityreduction",TRUE)
  cfr <- .iptm_resolve(nodecontrolfecundityreduction,1L,n_nodes,D,1L,"nodecontrolfecundityreduction",TRUE)
  efffr <- (1-fr*managing)*(1-cfr)
  if(D > 1L) {
    if(is.matrix(BlockedTransitionMortality)) {
      if(!identical(dim(BlockedTransitionMortality),c(n_nodes,D-1L)))
        stop("BlockedTransitionMortality must be scalar, length Nstages - 1, or nodes x (Nstages - 1)")
      bm <- BlockedTransitionMortality
    } else if(length(BlockedTransitionMortality)==1L) {
      bm <- matrix(BlockedTransitionMortality,n_nodes,D-1L)
    } else if(length(BlockedTransitionMortality)==D-1L) {
      bm <- matrix(rep(BlockedTransitionMortality,each=n_nodes),n_nodes,D-1L)
    } else stop("BlockedTransitionMortality must be scalar, length Nstages - 1, or nodes x (Nstages - 1)")
    if(any(!is.finite(bm)) || any(bm < 0 | bm > 1)) stop("BlockedTransitionMortality probabilities must be between 0 and 1")
  } else bm <- matrix(numeric(0),n_nodes,0L)
  norm_tm <- function(x,name) {
    if(is.null(x)) return(rep(list(NULL),max(0,D-1L)))
    z <- if(is.list(x)) x else rep(list(x),D-1L)
    if(length(z)!=D-1L) stop(name," must be a matrix or list length Nstages-1")
    for(k in seq_along(z)) if(!is.null(z[[k]])) { if(!is.matrix(z[[k]]) || !all(dim(z[[k]])==c(n_nodes,n_nodes))) stop(name," matrices must be nodes x nodes"); if(any(z[[k]]<0)||any(rowSums(z[[k]])>1+1e-10)) stop(name," rows must contain non-negative probabilities summing to <=1") }
    z
  }
  tsdd <- norm_tm(transition_sddprob,"transition_sddprob"); tldd <- norm_tm(transition_lddprob,"transition_lddprob")
  if(D > 1L) {
    if(!(length(transition_lddrate) %in% c(1L,D-1L))) stop("transition_lddrate must be scalar or length Nstages - 1")
    trrate <- rep_len(transition_lddrate,D-1L)
    if(any(!is.finite(trrate)) || any(trrate<0|trrate>1)) stop("transition_lddrate must be in [0,1]")
  } else trrate <- numeric(0)
  if(length(DispersalDensityFactor)!=1L) stop("DispersalDensityFactor must be a single value")
  # Fecundity from all pathogen states; offspring are susceptible.
  fecmeans <- numeric(n_nodes); mothers <- numeric(n_nodes)
  for(i in seq_len(n_nodes)) { A<-A_list[[i]]; f<-A[1,-1]; fecmeans[i]<-sum(f*N[i,-1]*efffr[i,-1]); mothers[i]<-sum(N[i,which(c(FALSE,f>0))]) }
  # Match the specialist vertebrate Birth hook contract from the frozen host
  # engine. A custom birth mean replaces transition-matrix fecundity, while
  # mother counts retain the host engine's footprint semantics.
  if(!is.null(nodebirthmean)) {
    nodebirthmean <- as.numeric(nodebirthmean)
    if(length(nodebirthmean)==1L) nodebirthmean <- rep(nodebirthmean,n_nodes)
    if(length(nodebirthmean)!=n_nodes || any(!is.finite(nodebirthmean)) || any(nodebirthmean<0))
      stop("nodebirthmean must contain one finite non-negative expected birth count per node")
    fecmeans <- nodebirthmean
  }
  if(!is.null(nodebirthmothers)) {
    nodebirthmothers <- as.numeric(nodebirthmothers)
    if(length(nodebirthmothers)==1L) nodebirthmothers <- rep(nodebirthmothers,n_nodes)
    if(length(nodebirthmothers)!=n_nodes || any(!is.finite(nodebirthmothers)) || any(nodebirthmothers<0))
      stop("nodebirthmothers must contain one finite non-negative value per node")
    mothers <- nodebirthmothers
  }
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
  # Reproductive dispersal and susceptible recruitment. Route realised integer
  # propagules using the same pest-validated TM semantics: first branch each
  # realised propagule to SDD/LDD, then route source-by-source to internal
  # destinations or the residual outside-landscape category.
  Pin <- numeric(n_nodes)
  Qin <- numeric(n_nodes)
  # Route reproductive propagules while keeping new offspring susceptible.
  .route_integer_propagules_pathogen_tm <- function(source_counts, Pmat, label) {
    arrivals <- numeric(n_nodes)
    if(!is.matrix(Pmat)) return(arrivals)
    if(!identical(dim(Pmat),c(n_nodes,n_nodes))) stop(label," must be an n_nodes x n_nodes matrix")
    active_sources <- which(source_counts > 0)
    if(!length(active_sources)) return(arrivals)
    for(ii in active_sources) {
      n_i <- floor(source_counts[ii])
      if(n_i <= 0) next
      p_internal <- pmax(0,as.numeric(Pmat[ii,]))
      rs <- sum(p_internal)
      if(!is.finite(rs) || rs > 1 + 1e-10) stop(label," source-row probabilities must be finite, non-negative, and sum to at most 1")
      if(rs > 1) p_internal <- p_internal / rs
      probs <- c(p_internal,max(0,1-sum(p_internal)))
      if(sum(probs) <= 0) next
      if(n_i < MaxInteger) {
        allocation <- as.numeric(rmultinom(1,size=n_i,prob=probs))
      } else {
        allocation <- floor(n_i * probs)
        allocation[length(allocation)] <- allocation[length(allocation)] + (n_i-sum(allocation))
      }
      arrivals <- arrivals + allocation[seq_len(n_nodes)]
    }
    arrivals
  }
  if(any(propagules > 0)) {
    if(lddrate <= 0) {
      ldd_sources <- numeric(n_nodes); sdd_sources <- propagules
    } else if(lddrate >= 1) {
      ldd_sources <- propagules; sdd_sources <- numeric(n_nodes)
    } else {
      ldd_sources <- rbinom(n_nodes,size=propagules,prob=lddrate)
      sdd_sources <- propagules-ldd_sources
    }
    Pin <- .route_integer_propagules_pathogen_tm(sdd_sources,sddprob,"sddprob")
    if(is.matrix(lddprob) && any(ldd_sources > 0)) {
      spread_reduction <- pmin(1,pmax(0,rep_len(nodespreadreduction,n_nodes)*rep_len(managing,n_nodes)))
      keep_prob <- 1-spread_reduction
      keep_sources <- ldd_sources
      stochastic_keep <- ldd_sources > 0 & keep_prob > 0 & keep_prob < 1
      keep_sources[keep_prob <= 0] <- 0
      if(any(stochastic_keep)) keep_sources[stochastic_keep] <- rbinom(sum(stochastic_keep),size=ldd_sources[stochastic_keep],prob=keep_prob[stochastic_keep])
      Qin <- .route_integer_propagules_pathogen_tm(keep_sources,lddprob,"lddprob")
    }
  }
  slots<-pmax(0,floor(node.seedbankK-N[,1])); Pin<-pmax(0,floor(Pin)); Qin<-pmax(0,floor(Qin)); env<-pmin(1,pmax(0,env))
  unrestricted<-all(is.finite(foot)&foot>=1)
  if(unrestricted) accessible<-as.numeric(Pin>0) else { area<-foot*nodeK; press<-numeric(n_nodes); if(any(Pin>0)&&any(mothers>0)) { num<-as.numeric((mothers*area)%*%sddprob); ok<-nodeK>0; press[ok]<-num[ok]/nodeK[ok] }; accessible<-pmin(1,pmax(0,-expm1(-press))) }
  nslots<-floor(slots*accessible); nslots<-ifelse(Pin>0&accessible>0,pmax(1,nslots),0); nslots<-pmin(nslots,slots); oslots<-pmax(0,slots-nslots)
  lp<-ifelse(nslots>0,env*Pin/pmax(nslots,1),0); lq<-ifelse(slots>0,env*Qin/pmax(slots,1),0); recs<-rbinom(n_nodes,nslots,-expm1(-(lp+lq)))+rbinom(n_nodes,oslots,-expm1(-lq)); recs<-pmin(recs,Pin+Qin)
  state[,1,"S"] <- state[,1,"S"]+recs
  # This helper performs host demography/transport only. Disease transitions
  # are deliberately applied later in the parent node timestep, after external
  # incursion and the specialist aggregate Interaction hook.
  N <- apply(state,c(1,2),sum)
  list(N=N, PathogenState=state)
}

# Assign pathogen states to newly arriving external vertebrate hosts.
.iv_node_external_pathogen <- function(state, N_before, N_after, ExternalPathogenStateProb,
                                       Pathogen, timestep, Ntimesteps) {
  state <- .iptm_reconcile(state, as.matrix(N_before))
  N_before <- as.matrix(N_before); N_after <- as.matrix(N_after)
  if (!all(dim(N_before) == dim(N_after))) stop("N_before and N_after dimensions differ")
  add <- pmax(0L, as.integer(N_after[,1] - N_before[,1]))
  state <- .iptm_reconcile(state, N_after)
  if (is.null(ExternalPathogenStateProb) || !any(add > 0L)) return(list(State=state, External=NULL))

  n_nodes <- nrow(N_after); states <- Pathogen$States; P <- length(states)
  x <- ExternalPathogenStateProb
  if (is.function(x)) {
    fm <- names(formals(x)); a <- list(timestep=timestep,n_nodes=n_nodes,states=states,Ntimesteps=Ntimesteps)
    if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a),fm)]
    x <- do.call(x,a)
  }
  if (is.null(dim(x))) {
    if (!is.numeric(x) || length(x) != P) stop("ExternalPathogenStateProb must have one probability per pathogen state")
    if (!is.null(names(x))) {
      if (!all(states %in% names(x))) stop("Named ExternalPathogenStateProb must contain every pathogen state")
      x <- x[states]
    }
    probs <- matrix(rep(as.numeric(x), each=n_nodes), nrow=n_nodes, ncol=P,
                    dimnames=list(NULL,states))
  } else {
    probs <- as.matrix(x)
    if (!identical(dim(probs), c(n_nodes,P))) stop("ExternalPathogenStateProb matrix must be nodes x pathogen states")
    if (!is.null(colnames(probs))) {
      if (!all(states %in% colnames(probs))) stop("ExternalPathogenStateProb columns must contain every pathogen state")
      probs <- probs[,states,drop=FALSE]
    }
  }
  if (any(!is.finite(probs)) || any(probs < 0) || any(rowSums(probs) <= 0))
    stop("ExternalPathogenStateProb must contain finite non-negative probabilities with positive row sums")
  probs <- probs / rowSums(probs)

  ext <- array(0L,c(n_nodes,ncol(N_after),P),dimnames=list(NULL,NULL,states))
  for (i in which(add > 0L)) {
    # Remove the default susceptible additions created by reconcile, then put
    # the same realised external hosts into the requested pathogen states.
    state[i,1,"S"] <- state[i,1,"S"] - add[i]
    z <- as.integer(rmultinom(1L,add[i],probs[i,]))
    state[i,1,] <- state[i,1,] + z
    ext[i,1,] <- z
  }
  list(State=state, External=ext)
}

# Draw pathogen detections from infectious hosts by node and stage.
.iv_node_pathogen_detect <- function(state, Pathogen, timestep, Ntimesteps) {
  n_nodes <- dim(state)[1]; D <- dim(state)[2]
  pdet <- .iptm_resolve(Pathogen$DetectionProb,timestep,n_nodes,D,Ntimesteps,"DetectionProb",TRUE)
  I <- state[, , "I", drop=FALSE][,,1]
  # Independent per-infectious-host detection across demographic stages.
  pnode <- 1 - apply((1-pdet)^I,1,prod)
  as.integer(rbinom(n_nodes,1L,pmin(1,pmax(0,pnode))))
}
