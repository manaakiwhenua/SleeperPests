###############################################################################
### INApest multiple-land-use pathogen analytical methods
###
### The MLU pathogen engine treats each node x land-use cell as a pathogen unit.
### R's column-major matrix flattening is preserved exactly: unit index
###   u = node + (land_use - 1) * n_nodes.
### This module therefore reuses the validated Meta exact network operators on
### the flattened pathogen-unit graph, then reshapes results back to node x LU.
###############################################################################

.ina_mlu_check_matrix <- function(x, n_nodes = NULL, n_landuses = NULL, name = deparse(substitute(x)), integer = FALSE) {
  if (!is.matrix(x)) stop(name, " must be a nodes x land-uses matrix")
  if (!is.null(n_nodes) && nrow(x) != n_nodes) stop(name, " has the wrong number of nodes")
  if (!is.null(n_landuses) && ncol(x) != n_landuses) stop(name, " has the wrong number of land uses")
  if (any(!is.finite(x))) stop(name, " must contain finite values")
  if (integer && any(x < 0 | x != floor(x))) stop(name, " must contain non-negative whole numbers")
  x
}

.ina_mlu_unit_index <- function(n_nodes, n_landuses) {
  z <- expand.grid(node = seq_len(n_nodes), land_use = seq_len(n_landuses), KEEP.OUT.ATTRS = FALSE)
  z$unit <- z$node + (z$land_use - 1L) * n_nodes
  z <- z[order(z$unit), , drop = FALSE]
  rownames(z) <- NULL
  z
}

.ina_mlu_resolve_unit <- function(x, n_nodes, n_landuses, name, integer = FALSE) {
  n_units <- n_nodes * n_landuses
  if (is.matrix(x)) {
    if (!all(dim(x) == c(n_nodes, n_landuses))) stop(name, " matrix must be nodes x land uses")
    z <- as.numeric(x)
  } else {
    z <- as.numeric(x)
    if (length(z) == 1L) z <- rep(z, n_units)
    else if (length(z) == n_units) z <- z
    else if (length(z) == n_landuses && n_landuses != n_nodes)
      z <- as.numeric(matrix(rep(z, each = n_nodes), nrow = n_nodes))
    else if (length(z) == n_nodes && n_nodes != n_landuses)
      z <- as.numeric(matrix(rep(z, n_landuses), nrow = n_nodes))
    else stop(name, " must be scalar, nodes x land uses, or an unambiguous node/LU/unit vector")
  }
  if (any(!is.finite(z))) stop(name, " must be finite")
  if (integer && any(z < 0 | z != floor(z))) stop(name, " must contain non-negative whole numbers")
  z
}

INApestMLUPathogenContactMatrix <- function(NodeContact, LandUseMixing = NULL) {
  NodeContact <- as.matrix(NodeContact)
  if (nrow(NodeContact) != ncol(NodeContact) || any(!is.finite(NodeContact)) || any(NodeContact < 0))
    stop("NodeContact must be a finite non-negative square source x target matrix")
  if (is.null(LandUseMixing)) return(NodeContact)
  LandUseMixing <- as.matrix(LandUseMixing)
  if (nrow(LandUseMixing) != ncol(LandUseMixing) || any(!is.finite(LandUseMixing)) || any(LandUseMixing < 0))
    stop("LandUseMixing must be a finite non-negative square source x target matrix")
  kronecker(LandUseMixing, NodeContact)
}

.ina_mlu_contact <- function(n_nodes, n_landuses, ContactMatrix = NULL,
                             NodeContact = NULL, LandUseMixing = NULL) {
  n_units <- n_nodes * n_landuses
  if (!is.null(ContactMatrix)) {
    C <- as.matrix(ContactMatrix)
    if (!all(dim(C) == c(n_units, n_units)) || any(!is.finite(C)) || any(C < 0))
      stop("ContactMatrix must be a finite non-negative pathogen-units x pathogen-units matrix")
    return(C)
  }
  if (is.null(NodeContact)) NodeContact <- diag(n_nodes)
  if (is.null(LandUseMixing)) LandUseMixing <- diag(n_landuses)
  C <- INApestMLUPathogenContactMatrix(NodeContact, LandUseMixing)
  if (!all(dim(C) == c(n_units, n_units))) stop("Combined contact matrix dimensions do not match nodes x land uses")
  C
}

.ina_mlu_reshape_time_unit <- function(x, n_nodes, n_landuses) {
  # x is timesteps x units. Return nodes x LU x timesteps.
  tt <- nrow(x)
  out <- array(0, dim = c(n_nodes, n_landuses, tt),
               dimnames = list(node = seq_len(n_nodes), land_use = seq_len(n_landuses), timestep = 0:(tt - 1L)))
  for (t in seq_len(tt)) out[, , t] <- matrix(x[t, ], nrow = n_nodes, ncol = n_landuses)
  out
}

INApestMLUPathogenGrowthOperator <- function(
    Model = c("SIS", "SIR", "SEIR"),
    HostPopulation,
    Beta,
    RecoveryProb,
    ProgressionProb = 1,
    PathogenMortalityProb = 0,
    ContactMatrix = NULL,
    NodeContact = NULL,
    LandUseMixing = NULL,
    Transmission = c("frequency", "density"),
    DensityScale = 1) {
  Model <- match.arg(Model); Transmission <- match.arg(Transmission)
  HostPopulation <- .ina_mlu_check_matrix(HostPopulation, name = "HostPopulation", integer = FALSE)
  n_nodes <- nrow(HostPopulation); n_landuses <- ncol(HostPopulation); n_units <- n_nodes * n_landuses
  N <- as.numeric(HostPopulation)
  beta <- .ina_mlu_resolve_unit(Beta, n_nodes, n_landuses, "Beta")
  rec <- .ina_mlu_resolve_unit(RecoveryProb, n_nodes, n_landuses, "RecoveryProb")
  prog <- .ina_mlu_resolve_unit(ProgressionProb, n_nodes, n_landuses, "ProgressionProb")
  mort <- .ina_mlu_resolve_unit(PathogenMortalityProb, n_nodes, n_landuses, "PathogenMortalityProb")
  ds <- .ina_mlu_resolve_unit(DensityScale, n_nodes, n_landuses, "DensityScale")
  C <- .ina_mlu_contact(n_nodes, n_landuses, ContactMatrix, NodeContact, LandUseMixing)
  g <- INApestMetaPathogenGrowthOperator(Model, N, beta, rec, prog, mort, C, Transmission, ds)
  map <- .ina_mlu_unit_index(n_nodes, n_landuses)
  active_states <- if (Model == "SEIR") 2L else 1L
  labels <- if (Model == "SEIR") {
    c(paste0("E_n", map$node, "_lu", map$land_use), paste0("I_n", map$node, "_lu", map$land_use))
  } else paste0("I_n", map$node, "_lu", map$land_use)
  dimnames(g$Operator) <- list(labels, labels)
  if (!is.null(g$TransmissionBlock))
    dimnames(g$TransmissionBlock) <- list(paste0("recipient_n", map$node, "_lu", map$land_use), paste0("source_n", map$node, "_lu", map$land_use))
  g$Model <- "INApestMetaMultipleLandUse"
  g$HostArchitecture <- "node x land-use pathogen units"
  g$UnitMap <- map
  g$HostPopulation <- HostPopulation
  g$Diagnostics <- c(g$Diagnostics,
                     "MLU flattening matches the stochastic engine: all nodes in land use 1, then all nodes in land use 2, and so on.",
                     "Land-use mixing and node contact can be supplied separately; the combined source x target contact matrix is kronecker(LandUseMixing, NodeContact).")
  class(g) <- c("INApestMLUPathogenGrowthOperator", "list")
  g
}

INApestMLUPathogenExactFixedN <- function(
    Model = c("SIS", "SIR", "SEIR"),
    Ntimesteps,
    HostPopulation,
    InitialInfected,
    InitialExposed = 0,
    InitialRecovered = 0,
    Beta,
    RecoveryProb,
    ProgressionProb = 1,
    ImmunityLossProb = 0,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    ContactMatrix = NULL,
    NodeContact = NULL,
    LandUseMixing = NULL,
    Transmission = c("frequency", "density"),
    DensityScale = 1,
    OutsideNodes = integer(0),
    OutsideUnits = integer(0),
    MaxStates = 50000L,
    ReturnOperator = FALSE) {
  Model <- match.arg(Model); Transmission <- match.arg(Transmission)
  Nmat <- .ina_mlu_check_matrix(HostPopulation, name = "HostPopulation", integer = TRUE)
  n_nodes <- nrow(Nmat); n_landuses <- ncol(Nmat); n_units <- n_nodes * n_landuses
  N <- as.integer(as.numeric(Nmat))
  I0 <- as.integer(.ina_mlu_resolve_unit(InitialInfected, n_nodes, n_landuses, "InitialInfected", TRUE))
  E0 <- as.integer(.ina_mlu_resolve_unit(InitialExposed, n_nodes, n_landuses, "InitialExposed", TRUE))
  R0 <- as.integer(.ina_mlu_resolve_unit(InitialRecovered, n_nodes, n_landuses, "InitialRecovered", TRUE))
  if (any(I0 + E0 + R0 > N)) stop("Initial pathogen-state counts exceed HostPopulation")
  beta <- .ina_mlu_resolve_unit(Beta, n_nodes, n_landuses, "Beta")
  rec <- .ina_mlu_resolve_unit(RecoveryProb, n_nodes, n_landuses, "RecoveryProb")
  prog <- .ina_mlu_resolve_unit(ProgressionProb, n_nodes, n_landuses, "ProgressionProb")
  wan <- .ina_mlu_resolve_unit(ImmunityLossProb, n_nodes, n_landuses, "ImmunityLossProb")
  ip <- .ina_mlu_resolve_unit(IntroductionProb, n_nodes, n_landuses, "IntroductionProb")
  inum <- as.integer(.ina_mlu_resolve_unit(IntroductionNumber, n_nodes, n_landuses, "IntroductionNumber", TRUE))
  ds <- .ina_mlu_resolve_unit(DensityScale, n_nodes, n_landuses, "DensityScale")
  C <- .ina_mlu_contact(n_nodes, n_landuses, ContactMatrix, NodeContact, LandUseMixing)
  outside <- unique(as.integer(OutsideUnits))
  if (length(OutsideNodes)) {
    if (any(!OutsideNodes %in% seq_len(n_nodes))) stop("OutsideNodes contains invalid node indices")
    map <- .ina_mlu_unit_index(n_nodes, n_landuses)
    outside <- unique(c(outside, map$unit[map$node %in% OutsideNodes]))
  }
  if (any(!outside %in% seq_len(n_units))) stop("OutsideUnits contains invalid pathogen-unit indices")
  if (Model == "SIS") {
    ans <- INApestMetaSISExactNetworkFixedN(Ntimesteps, N, I0, beta, rec, C, ip, inum,
                                            Transmission, ds, OutsideNodes = outside,
                                            ReturnOperator = ReturnOperator)
    EI <- ans$ExpectedInfectedByNode
    EE <- NULL
  } else {
    ans <- INApestMetaCompartmentExactNetworkFixedN(Model, Ntimesteps, N, I0, E0, R0,
                                                     beta, rec, prog, wan, ip, inum, C,
                                                     Transmission, ds, OutsideNodes = outside,
                                                     MaxStates = MaxStates,
                                                     ReturnOperator = ReturnOperator)
    EI <- ans$ExpectedInfectedByNode
    EE <- ans$ExpectedExposedByNode
  }
  map <- .ina_mlu_unit_index(n_nodes, n_landuses)
  out <- ans
  out$Model <- "INApestMetaMultipleLandUse"
  out$HostArchitecture <- "fixed abundance by node x land use"
  out$PathogenModel <- Model
  out$UnitMap <- map
  out$HostPopulation <- Nmat
  out$ExpectedInfectedByUnit <- EI
  out$ExpectedInfected <- .ina_mlu_reshape_time_unit(EI, n_nodes, n_landuses)
  if (!is.null(EE)) {
    out$ExpectedExposedByUnit <- EE
    out$ExpectedExposed <- .ina_mlu_reshape_time_unit(EE, n_nodes, n_landuses)
  }
  if (!is.null(out$Escape)) {
    out$Escape$OutsideUnits <- outside
    out$Escape$OutsideNodes <- as.integer(OutsideNodes)
  }
  out$Growth <- INApestMLUPathogenGrowthOperator(Model, Nmat, beta, rec, prog,
                                                  PathogenMortalityProb = if (Model == "SIS") 0 else 0,
                                                  ContactMatrix = C, Transmission = Transmission,
                                                  DensityScale = ds)
  out$Diagnostics <- c(out$Diagnostics,
                       "Exact MLU solution is obtained on the flattened node x land-use pathogen-unit graph and then reshaped back to nodes x land uses.",
                       "This fixed-host exact branch is intended as a small-system validator; host turnover and management are added in separate MLU analytical branches.")
  class(out) <- c("INApestMLUPathogenExactFixedN", "list")
  out
}

print.INApestMLUPathogenExactFixedN <- function(x, ...) {
  cat("INApest MLU pathogen exact fixed-host solution\n")
  cat("Pathogen model:", x$PathogenModel, "\n")
  cat("Nodes x land uses:", nrow(x$HostPopulation), "x", ncol(x$HostPopulation), "\n")
  cat("Rare-pathogen multiplier:", format(x$Growth$IntrinsicRarePathogenMultiplier, digits = 7), "\n")
  if (!is.null(x$Escape)) cat("Escape probability by horizon:", format(x$Escape$ProbabilityByHorizon, digits = 7), "\n")
  invisible(x)
}

###############################################################################
### Exact one-node MLU SIS with explicit gross host turnover
###############################################################################

.ina_mlu_sis_local_ni_states <- function(K) {
  rows <- list(); k <- 0L
  for (N in 0:as.integer(K)) for (I in 0:N) { k <- k + 1L; rows[[k]] <- c(N=N,I=I) }
  z <- as.data.frame(do.call(rbind, rows)); z$N <- as.integer(z$N); z$I <- as.integer(z$I); z
}

.ina_mlu_sis_turnover_state_table <- function(K) {
  K <- as.integer(K); L <- length(K)
  local <- lapply(K, .ina_mlu_sis_local_ni_states)
  dims <- vapply(local, nrow, integer(1))
  grid <- as.matrix(expand.grid(lapply(dims, seq_len), KEEP.OUT.ATTRS=FALSE)); storage.mode(grid) <- "integer"
  out <- matrix(0L, nrow(grid), 2L*L)
  cn <- character(2L*L)
  for (lu in seq_len(L)) {
    out[,2L*lu-1L] <- local[[lu]]$N[grid[,lu]]
    out[,2L*lu] <- local[[lu]]$I[grid[,lu]]
    cn[2L*lu-1L] <- paste0("N_lu",lu); cn[2L*lu] <- paste0("I_lu",lu)
  }
  colnames(out) <- cn
  attr(out,"LocalStates") <- local
  out
}

.ina_mlu_sis_turnover_operator_one_node <- function(
    K,
    Survival = 1,
    RecruitToCapacity = FALSE,
    Beta = 0,
    RecoveryProb = 0,
    PathogenMortalityProb = 0,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    LandUseMixing = NULL,
    Transmission = c("frequency","density"),
    DensityScale = 1) {
  Transmission <- match.arg(Transmission); K <- as.integer(K); L <- length(K)
  if(any(K<0))stop("K must be non-negative whole numbers")
  res <- function(x,name,int=FALSE){z<-.ina_mlu_resolve_unit(x,1L,L,name,int);z}
  surv<-res(Survival,"Survival");beta<-res(Beta,"Beta");rec<-res(RecoveryProb,"RecoveryProb");pm<-res(PathogenMortalityProb,"PathogenMortalityProb");ip<-res(IntroductionProb,"IntroductionProb");inum<-as.integer(res(IntroductionNumber,"IntroductionNumber",TRUE));ds<-res(DensityScale,"DensityScale")
  if(any(surv<0|surv>1)||any(beta<0)||any(rec<0|rec>1)||any(pm<0|pm>1)||any(rec+pm>1+1e-12)||any(ip<0|ip>1)||any(ds<=0))stop("Invalid turnover/pathogen parameters")
  C <- .ina_mlu_contact(1L,L,NodeContact=matrix(1,1,1),LandUseMixing=if(is.null(LandUseMixing))diag(L)else LandUseMixing)
  st <- .ina_mlu_sis_turnover_state_table(K); keys <- apply(st,1,paste,collapse=":"); idx <- setNames(seq_len(nrow(st)),keys); T <- matrix(0,nrow(st),nrow(st),dimnames=list(keys,keys))
  for(row in seq_len(nrow(st))) {
    N0 <- as.integer(st[row,seq(1,2L*L,by=2)]); I0 <- as.integer(st[row,seq(2,2L*L,by=2)]); S0 <- N0-I0
    host_lists <- vector("list",L)
    for(lu in seq_len(L)) {
      z <- list(); kk <- 0L
      for(ss in 0:S0[lu]) for(ii in 0:I0[lu]) {
        p <- dbinom(ss,S0[lu],surv[lu])*dbinom(ii,I0[lu],surv[lu]); if(p==0)next
        Npre <- ss+ii
        if(isTRUE(RecruitToCapacity)) Npre <- K[lu]
        kk<-kk+1L;z[[kk]]<-c(N=Npre,I=ii,prob=p)
      }
      host_lists[[lu]] <- as.data.frame(do.call(rbind,z))
    }
    hg <- as.matrix(expand.grid(lapply(host_lists,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(hg)<-"integer"
    for(hb in seq_len(nrow(hg))) {
      Npre<-Ipre<-numeric(L);ph<-1
      for(lu in seq_len(L)){z<-host_lists[[lu]][hg[hb,lu],];Npre[lu]<-z$N;Ipre[lu]<-z$I;ph<-ph*z$prob}
      if(ph==0)next
      pressure<-as.numeric(crossprod(Ipre,C));den<-as.numeric(crossprod(Npre,C));foi<-if(Transmission=="frequency")beta*ifelse(den>0,pressure/den,0)else beta*pressure/ds;pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)))
      dlist<-vector("list",L)
      for(lu in seq_len(L)) dlist[[lu]] <- as.data.frame(.ina_meta_sis_disease_outcomes_pinf(as.integer(Npre[lu]),as.integer(Ipre[lu]),pinf[lu],rec[lu],pm[lu],ip[lu],inum[lu]))
      dg<-as.matrix(expand.grid(lapply(dlist,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(dg)<-"integer"
      for(db in seq_len(nrow(dg))) {
        dest<-integer(2L*L);pd<-ph
        for(lu in seq_len(L)){z<-dlist[[lu]][dg[db,lu],];dest[2L*lu-1L]<-as.integer(z$N);dest[2L*lu]<-as.integer(z$I);pd<-pd*z$prob}
        ky<-paste(dest,collapse=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+pd
      }
    }
  }
  er<-max(abs(rowSums(T)-1));if(er>1e-10)stop("MLU turnover operator rows do not sum to 1; max error=",er)
  attr(T,"States")<-st;T
}

INApestMLUSISExactTurnoverOneNode <- function(
    Ntimesteps,
    K,
    InitialHostPopulation,
    InitialInfected,
    Survival = 1,
    RecruitToCapacity = FALSE,
    Beta = 0,
    RecoveryProb = 0,
    PathogenMortalityProb = 0,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    LandUseMixing = NULL,
    Transmission = c("frequency","density"),
    DensityScale = 1,
    ReturnOperator = FALSE) {
  Transmission<-match.arg(Transmission);K<-as.integer(K);L<-length(K)
  N0<-as.integer(.ina_mlu_resolve_unit(InitialHostPopulation,1L,L,"InitialHostPopulation",TRUE));I0<-as.integer(.ina_mlu_resolve_unit(InitialInfected,1L,L,"InitialInfected",TRUE));if(any(N0>K)||any(I0>N0))stop("Invalid initial MLU turnover state")
  T<-.ina_mlu_sis_turnover_operator_one_node(K,Survival,RecruitToCapacity,Beta,RecoveryProb,PathogenMortalityProb,IntroductionProb,IntroductionNumber,LandUseMixing,Transmission,DensityScale)
  st<-attr(T,"States");keys<-rownames(T);init<-as.vector(rbind(N0,I0));dist<-numeric(nrow(st));dist[match(paste(init,collapse=":"),keys)]<-1
  EN<-EI<-matrix(0,Ntimesteps+1L,L);EN[1,]<-N0;EI[1,]<-I0;presence<-numeric(Ntimesteps+1L);presence[1]<-as.numeric(any(I0>0));full<-matrix(0,Ntimesteps+1L,nrow(st));full[1,]<-dist
  Ncols<-seq(1,2L*L,by=2);Icols<-seq(2,2L*L,by=2)
  if(Ntimesteps>0)for(tt in seq_len(Ntimesteps)){dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;EN[tt+1,]<-as.numeric(dist%*%st[,Ncols,drop=FALSE]);EI[tt+1,]<-as.numeric(dist%*%st[,Icols,drop=FALSE]);presence[tt+1]<-sum(dist[rowSums(st[,Icols,drop=FALSE])>0])}
  out<-list(Model="INApestMetaMultipleLandUse",PathogenModel="SIS",HostAssumption=if(RecruitToCapacity)"gross survival followed by susceptible replacement to K in each land use"else"gross survival with no recruitment",Exact=TRUE,ExpectedHostByLandUse=EN,ExpectedInfectedByLandUse=EI,PathogenPresenceProbability=presence,PathogenFreedomProbability=1-presence,StateTable=st,StateDistribution=full,Diagnostics=c("Exact one-node MLU SIS turnover validator.","Host survival is applied to actual S/I hosts before any replacement recruitment; replacement hosts enter S.","Pathogen transmission occurs only after host turnover and uses the post-turnover land-use abundances."));if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMLUSISExactTurnoverOneNode","list");out
}

###############################################################################
### Shared node-level information / management across multiple land uses
###
### The stochastic MLU engine stores host/pathogen states by node x land use,
### but HaveInfo and LastKnownPresence are node-level states. These exact
### validators preserve that architecture rather than treating land uses as
### independent information units.
###############################################################################

.ina_mlu_model_states <- function(Model) {
  Model <- match.arg(Model, c("SIS","SIR","SEIR"))
  if (Model == "SIS") c("S","I") else if (Model == "SIR") c("S","I","R") else c("S","E","I","R")
}

.ina_mlu_local_compartment_states <- function(Model, K) {
  Model <- match.arg(Model, c("SIS","SIR","SEIR")); K <- as.integer(K)
  if (K < 0L) stop("K must be non-negative")
  if (Model == "SIS") {
    out <- list(); kk <- 0L
    for (N in 0:K) for (I in 0:N) {
      kk <- kk + 1L; out[[kk]] <- c(S=N-I, I=I)
    }
    z <- as.data.frame(do.call(rbind,out)); z$S <- as.integer(z$S); z$I <- as.integer(z$I); z$N <- z$S+z$I
    return(z)
  }
  .ina_meta_compartment_states_uptoK(Model, K)
}

.ina_mlu_info_age_rows <- function(InfoPersistenceSteps=NA) {
  if (is.na(InfoPersistenceSteps)) return(data.frame(H=c(0L,1L), A=c(-2L,-2L)))
  Kinfo <- as.integer(InfoPersistenceSteps)
  if (Kinfo < 0L || Kinfo != InfoPersistenceSteps) stop("InfoPersistenceSteps must be NA or a non-negative whole number")
  ages <- unique(c(-1L, 0:max(0L,Kinfo-1L)))
  rbind(data.frame(H=0L,A=-2L), data.frame(H=1L,A=ages))
}

.ina_mlu_one_node_info_states <- function(Model, K, InfoPersistenceSteps=NA, MaxStates=50000L) {
  Model <- match.arg(Model,c("SIS","SIR","SEIR")); K <- as.integer(K); L <- length(K)
  if (!L || any(K < 0L)) stop("K must contain non-negative land-use capacities")
  local <- lapply(K, function(k) .ina_mlu_local_compartment_states(Model,k))
  dims <- vapply(local,nrow,integer(1)); bio_n <- prod(dims)
  ages <- .ina_mlu_info_age_rows(InfoPersistenceSteps)
  nst <- bio_n*nrow(ages)
  if (nst > MaxStates) stop("Exact one-node MLU information state space has ",nst," states; increase MaxStates deliberately or use stochastic MLU")
  g <- as.matrix(expand.grid(lapply(dims,seq_len),KEEP.OUT.ATTRS=FALSE)); storage.mode(g)<-"integer"
  states <- .ina_mlu_model_states(Model)
  bio <- matrix(0L,nrow(g),length(states)*L)
  cn <- character(ncol(bio))
  for (lu in seq_len(L)) for (j in seq_along(states)) {
    col <- (lu-1L)*length(states)+j
    bio[,col] <- local[[lu]][[states[j]]][g[,lu]]
    cn[col] <- paste0(states[j],"_lu",lu)
  }
  colnames(bio) <- cn
  out <- bio[rep(seq_len(nrow(bio)), each=nrow(ages)),,drop=FALSE]
  aa <- ages[rep(seq_len(nrow(ages)), times=nrow(bio)),,drop=FALSE]
  out <- cbind(out,H=as.integer(aa$H),A=as.integer(aa$A))
  storage.mode(out)<-"integer"
  attr(out,"ModelStates")<-states; attr(out,"K")<-K; out
}

.ina_mlu_compartment_disease_outcomes_pinf <- function(state, Model, p_inf,
                                                        RecoveryProb=0, ProgressionProb=1,
                                                        ImmunityLossProb=0, PathogenMortalityProb=0,
                                                        IntroductionProb=0, IntroductionNumber=1) {
  Model <- match.arg(Model,c("SIR","SEIR")); states <- .ina_mlu_model_states(Model)
  S0 <- as.integer(state[["S"]]); I0 <- as.integer(state[["I"]]); R0 <- as.integer(state[["R"]]); E0 <- if(Model=="SEIR")as.integer(state[["E"]])else 0L
  rec <- as.numeric(RecoveryProb); prog <- as.numeric(ProgressionProb); wan <- as.numeric(ImmunityLossProb); mort <- as.numeric(PathogenMortalityProb)
  ip <- as.numeric(IntroductionProb); inum <- as.integer(IntroductionNumber); p_inf <- .ina_meta_sis_clip01(p_inf)
  if(any(!is.finite(c(rec,prog,wan,mort,ip)))||any(c(rec,prog,wan,mort,ip)<0)||any(c(rec,prog,wan,mort,ip)>1)||rec+mort>1+1e-12) stop("Invalid disease probabilities")
  out<-list();kk<-0L
  for(x in 0:S0){px<-dbinom(x,S0,p_inf);remS<-S0-x;mx<-min(inum,remS);ivals<-if(ip>0&&mx>0)c(0L,mx)else 0L;ips<-if(ip>0&&mx>0)c(1-ip,ip)else 1
    for(iv in seq_along(ivals)){intro<-ivals[iv];pi<-ips[iv]
      for(pg in if(Model=="SEIR")0:E0 else 0L){pp<-if(Model=="SEIR")dbinom(pg,E0,prog)else 1
        for(lose in 0:R0){pl<-dbinom(lose,R0,wan)
          for(stay in 0:I0)for(recover in 0:(I0-stay)){die<-I0-stay-recover;pr<-dmultinom(c(stay,recover,die),prob=c(1-rec-mort,rec,mort));if(pr==0)next
            if(Model=="SIR") z<-c(S=S0-x-intro+lose,I=stay+x+intro,R=R0-lose+recover)
            else z<-c(S=S0-x-intro+lose,E=E0-pg+x+intro,I=stay+pg,R=R0-lose+recover)
            kk<-kk+1L;out[[kk]]<-c(z,prob=unname(px*pi*pp*pl*pr))
  }}}}}
  z<-as.data.frame(do.call(rbind,out));for(nm in states)z[[nm]]<-as.integer(z[[nm]]);z$prob<-as.numeric(z$prob);z<-z[z$prob>0,,drop=FALSE]
  key<-apply(z[,states,drop=FALSE],1,paste,collapse=":");ps<-tapply(z$prob,key,sum);first<-match(names(ps),key);ans<-z[first,states,drop=FALSE];ans$prob<-as.numeric(ps);rownames(ans)<-NULL;ans
}

.ina_mlu_local_management_mortality_branches <- function(counts, states, Survival, ManageProb, MortalityProb, H) {
  out<-list();kk<-0L;mprob<-ManageProb*H
  for(managing in 0:1){pm<-if(managing==1)mprob else 1-mprob;if(pm<=0)next;mm<-MortalityProb*managing
    ls<-lapply(states,function(nm).ina_meta_class_mortality_outcomes(as.integer(counts[[nm]]),Survival,mm))
    g<-as.matrix(expand.grid(lapply(ls,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(g)<-"integer"
    for(rr in seq_len(nrow(g))){surv<-setNames(integer(length(states)),states);md<-0L;p<-pm
      for(j in seq_along(states)){z<-ls[[j]][g[rr,j],];surv[j]<-as.integer(z["survive"]);md<-md+as.integer(z["manage_death"]);p<-p*unname(z["prob"])}
      if(p>0){kk<-kk+1L;out[[kk]]<-c(surv,manage_deaths=md,managing=managing,prob=p)}
    }
  }
  z<-as.data.frame(do.call(rbind,out));for(nm in states)z[[nm]]<-as.integer(z[[nm]]);z$manage_deaths<-as.integer(z$manage_deaths);z$managing<-as.integer(z$managing);z$prob<-as.numeric(z$prob);z
}

.ina_mlu_one_node_info_prebranches <- function(state_row, Model, L, Survival, ManageProb, MortalityProb,
                                                host_acq, InfoPersistenceSteps, InfoRetentionProb) {
  states <- .ina_mlu_model_states(Model); getv<-function(nm)as.integer(if(is.matrix(state_row)||is.data.frame(state_row))state_row[1,nm]else state_row[[nm]]); H0<-getv("H");A0<-getv("A");finite<-!is.na(InfoPersistenceSteps);Kinfo<-if(finite)as.integer(InfoPersistenceSteps)else NA_integer_
  local<-vector("list",L)
  for(lu in seq_len(L)){
    cc<-setNames(vapply(states,function(nm)getv(paste0(nm,"_lu",lu)),integer(1)),states)
    local[[lu]]<-.ina_mlu_local_management_mortality_branches(cc,states,Survival[lu],ManageProb[lu],MortalityProb[lu],H0)
  }
  g<-as.matrix(expand.grid(lapply(local,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(g)<-"integer";out<-list();kk<-0L
  for(rr in seq_len(nrow(g))){vals<-integer(length(states)*L);names(vals)<-as.vector(vapply(seq_len(L),function(lu)paste0(states,"_lu",lu),character(length(states))));p<-1;md<-0L;mg<-integer(L)
    for(lu in seq_len(L)){z<-local[[lu]][g[rr,lu],];for(nm in states)vals[paste0(nm,"_lu",lu)]<-as.integer(z[[nm]]);md<-md+as.integer(z$manage_deaths);mg[lu]<-as.integer(z$managing);p<-p*z$prob}
    if(p==0)next;evidence<-host_acq&&md>0
    info<-list()
    if(finite){
      if(H0==0)info[[1]]<-c(H=0,A=-2,p=1)
      else if(evidence){if(Kinfo<=0)info[[1]]<-c(H=0,A=-2,p=1)else info[[1]]<-c(H=1,A=0,p=1)}
      else if(A0<0)info[[1]]<-c(H=0,A=-2,p=1)
      else {ag<-A0+1L;if(ag>=Kinfo)info[[1]]<-c(H=0,A=-2,p=1)else info[[1]]<-c(H=1,A=ag,p=1)}
    } else {
      if(H0==1L&&InfoRetentionProb<1){info[[1]]<-c(H=1,A=-2,p=InfoRetentionProb);info[[2]]<-c(H=0,A=-2,p=1-InfoRetentionProb)} else info[[1]]<-c(H=H0,A=-2,p=1)
    }
    for(q in info){kk<-kk+1L;out[[kk]]<-c(vals,H=as.integer(q["H"]),A=as.integer(q["A"]),manage_pattern=sum(mg*2^(seq_len(L)-1L)),prob=unname(p*q["p"]))}
  }
  z<-as.data.frame(do.call(rbind,out));for(nm in c(names(vals),"H","A","manage_pattern"))z[[nm]]<-as.integer(z[[nm]]);z$prob<-as.numeric(z$prob);z
}

.ina_mlu_info_operator_one_node <- function(Model,K,Survival=1,ManageProb=0,MortalityProb=0,
                                             HostDetectionProb=0,PathogenDetectionProb=0,
                                             InformationAcquisition=c("host","pathogen","both"),
                                             InfoPersistenceSteps=NA,InfoRetentionProb=1,
                                             Beta=0,RecoveryProb=0,ProgressionProb=1,ImmunityLossProb=0,
                                             PathogenMortalityProb=0,IntroductionProb=0,IntroductionNumber=1,
                                             LandUseMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
                                             MaxStates=50000L) {
  Model<-match.arg(Model,c("SIS","SIR","SEIR"));InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission);K<-as.integer(K);L<-length(K);states<-.ina_mlu_model_states(Model)
  res<-function(x,nm,int=FALSE).ina_mlu_resolve_unit(x,1L,L,nm,int)
  surv<-res(Survival,"Survival");mng<-res(ManageProb,"ManageProb");mort<-res(MortalityProb,"MortalityProb");hd<-res(HostDetectionProb,"HostDetectionProb");pd<-res(PathogenDetectionProb,"PathogenDetectionProb");beta<-res(Beta,"Beta");rec<-res(RecoveryProb,"RecoveryProb");prog<-res(ProgressionProb,"ProgressionProb");wan<-res(ImmunityLossProb,"ImmunityLossProb");pmort<-res(PathogenMortalityProb,"PathogenMortalityProb");ip<-res(IntroductionProb,"IntroductionProb");inum<-as.integer(res(IntroductionNumber,"IntroductionNumber",TRUE));ds<-res(DensityScale,"DensityScale")
  probs<-c(surv,mng,mort,hd,pd,rec,prog,wan,pmort,ip,InfoRetentionProb);if(any(!is.finite(probs))||any(probs<0)||any(probs>1)||any(rec+pmort>1+1e-12)||any(beta<0)||any(ds<=0))stop("Invalid MLU information/pathogen parameters")
  C<-.ina_mlu_contact(1L,L,NodeContact=matrix(1,1,1),LandUseMixing=if(is.null(LandUseMixing))diag(L)else LandUseMixing)
  st<-.ina_mlu_one_node_info_states(Model,K,InfoPersistenceSteps,MaxStates);bio_cols<-setdiff(colnames(st),c("H","A"));keys<-apply(st[,c(bio_cols,"H","A"),drop=FALSE],1,paste,collapse=":");idx<-setNames(seq_len(nrow(st)),keys);T<-matrix(0,nrow(st),nrow(st),dimnames=list(keys,keys));host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both");finite<-!is.na(InfoPersistenceSteps)
  for(row in seq_len(nrow(st))){pre<-.ina_mlu_one_node_info_prebranches(st[row,],Model,L,surv,mng,mort,host_acq,InfoPersistenceSteps,InfoRetentionProb)
    for(pb in seq_len(nrow(pre))){Ipre<-vapply(seq_len(L),function(lu)as.numeric(pre[pb,paste0("I_lu",lu)]),numeric(1));Npre<-vapply(seq_len(L),function(lu)sum(as.numeric(pre[pb,paste0(states,"_lu",lu)])),numeric(1));pressure<-as.numeric(crossprod(Ipre,C));den<-as.numeric(crossprod(Npre,C));foi<-if(Transmission=="frequency")beta*ifelse(den>0,pressure/den,0)else beta*pressure/ds;pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)))
      dl<-vector("list",L)
      for(lu in seq_len(L)){
        if(Model=="SIS")dl[[lu]]<-as.data.frame(.ina_meta_sis_disease_outcomes_pinf(as.integer(Npre[lu]),as.integer(Ipre[lu]),pinf[lu],rec[lu],pmort[lu],ip[lu],inum[lu]))
        else {cc<-setNames(vapply(states,function(nm)as.integer(pre[pb,paste0(nm,"_lu",lu)]),integer(1)),states);dl[[lu]]<-.ina_mlu_compartment_disease_outcomes_pinf(cc,Model,pinf[lu],rec[lu],prog[lu],wan[lu],pmort[lu],ip[lu],inum[lu])}
      }
      dg<-as.matrix(expand.grid(lapply(dl,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(dg)<-"integer"
      for(dd in seq_len(nrow(dg))){dest<-integer(length(bio_cols));names(dest)<-bio_cols;pbase<-pre$prob[pb]
        N1<-I1<-numeric(L)
        for(lu in seq_len(L)){z<-dl[[lu]][dg[dd,lu],];pbase<-pbase*z$prob
          if(Model=="SIS"){dest[paste0("S_lu",lu)]<-as.integer(z$N-z$I);dest[paste0("I_lu",lu)]<-as.integer(z$I)}else for(nm in states)dest[paste0(nm,"_lu",lu)]<-as.integer(z[[nm]])
          N1[lu]<-sum(dest[paste0(states,"_lu",lu)]);I1[lu]<-dest[paste0("I_lu",lu)]
        }
        if(pbase==0)next;qh<-if(host_acq)1-prod((1-hd)^N1)else 0;qp<-if(path_acq)1-prod((1-pd)^I1)else 0;qev<-1-(1-qh)*(1-qp);Hm<-as.integer(pre$H[pb]);Am<-as.integer(pre$A[pb])
        if(qev>0){ky<-paste(c(dest,1L,if(finite)0L else -2L),collapse=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+pbase*qev}
        if(qev<1){ky<-paste(c(dest,Hm,Am),collapse=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+pbase*(1-qev)}
      }
    }
  }
  er<-max(abs(rowSums(T)-1));if(er>1e-10)stop("MLU one-node information operator row error ",er);attr(T,"States")<-st;attr(T,"LandUseCount")<-L;T
}

INApestMLUPathogenExactInformationOneNode <- function(Model=c("SIS","SIR","SEIR"),Ntimesteps,K,
    InitialHostPopulation,InitialInfected=0,InitialExposed=0,InitialRecovered=0,InitialInfo=0,
    Survival=1,ManageProb=0,MortalityProb=0,HostDetectionProb=0,PathogenDetectionProb=0,
    InformationAcquisition=c("host","pathogen","both"),InfoPersistenceSteps=NA,InfoRetentionProb=1,
    Beta=0,RecoveryProb=0,ProgressionProb=1,ImmunityLossProb=0,PathogenMortalityProb=0,
    IntroductionProb=0,IntroductionNumber=1,LandUseMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
    MaxStates=50000L,ReturnOperator=FALSE) {
  Model<-match.arg(Model);InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission);K<-as.integer(K);L<-length(K);states<-.ina_mlu_model_states(Model)
  resi<-function(x,nm).ina_mlu_resolve_unit(x,1L,L,nm,TRUE);N0<-as.integer(resi(InitialHostPopulation,"InitialHostPopulation"));I0<-as.integer(resi(InitialInfected,"InitialInfected"));E0<-if(Model=="SEIR")as.integer(resi(InitialExposed,"InitialExposed"))else integer(L);R0<-if(Model%in%c("SIR","SEIR"))as.integer(resi(InitialRecovered,"InitialRecovered"))else integer(L)
  if(any(N0>K)||any(I0+E0+R0>N0)||!InitialInfo%in%c(0,1))stop("Invalid initial MLU information state")
  T<-.ina_mlu_info_operator_one_node(Model,K,Survival,ManageProb,MortalityProb,HostDetectionProb,PathogenDetectionProb,InformationAcquisition,InfoPersistenceSteps,InfoRetentionProb,Beta,RecoveryProb,ProgressionProb,ImmunityLossProb,PathogenMortalityProb,IntroductionProb,IntroductionNumber,LandUseMixing,Transmission,DensityScale,MaxStates)
  st<-attr(T,"States");bio_cols<-setdiff(colnames(st),c("H","A"));keys<-rownames(T);finite<-!is.na(InfoPersistenceSteps);hd<-.ina_mlu_resolve_unit(HostDetectionProb,1L,L,"HostDetectionProb");pd<-.ina_mlu_resolve_unit(PathogenDetectionProb,1L,L,"PathogenDetectionProb");host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both")
  init<-integer(length(bio_cols));names(init)<-bio_cols
  for(lu in seq_len(L)){init[paste0("S_lu",lu)]<-N0[lu]-I0[lu]-E0[lu]-R0[lu];init[paste0("I_lu",lu)]<-I0[lu];if(Model=="SEIR")init[paste0("E_lu",lu)]<-E0[lu];if(Model%in%c("SIR","SEIR"))init[paste0("R_lu",lu)]<-R0[lu]}
  qh<-if(host_acq)1-prod((1-hd)^N0)else 0;qp<-if(path_acq)1-prod((1-pd)^I0)else 0;qev<-1-(1-qh)*(1-qp);dist<-numeric(nrow(st));
  if(qev>0){ky<-paste(c(init,1L,if(finite)0L else -2L),collapse=":");dist[match(ky,keys)]<-qev}
  if(qev<1){H<-as.integer(InitialInfo);A<-if(H==1L&&finite)-1L else -2L;ky<-paste(c(init,H,A),collapse=":");ii<-match(ky,keys);dist[ii]<-dist[ii]+1-qev}
  full<-matrix(0,Ntimesteps+1L,nrow(st));full[1,]<-dist;PH<-numeric(Ntimesteps+1L);PH[1]<-sum(dist*st[,"H"]);mng<-.ina_mlu_resolve_unit(ManageProb,1L,L,"ManageProb");PM<-matrix(0,Ntimesteps,L);means<-lapply(states,function(nm)matrix(0,Ntimesteps+1L,L));names(means)<-states
  for(nm in states)for(lu in seq_len(L))means[[nm]][1,lu]<-sum(dist*st[,paste0(nm,"_lu",lu)])
  for(tt in seq_len(Ntimesteps)){PM[tt,]<-mng*sum(dist*st[,"H"]);dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;PH[tt+1]<-sum(dist*st[,"H"]);for(nm in states)for(lu in seq_len(L))means[[nm]][tt+1,lu]<-sum(dist*st[,paste0(nm,"_lu",lu)])}
  EN<-matrix(0,Ntimesteps+1L,L);for(nm in states)EN<-EN+means[[nm]]
  out<-list(Model="INApestMetaMultipleLandUse",PathogenModel=Model,Exact=TRUE,InformationAcquisition=InformationAcquisition,ExpectedHostByLandUse=EN,ExpectedStatesByLandUse=means,ExpectedInfectedByLandUse=means$I,InformationProbability=PH,ManagingProbabilityByLandUse=PM,StateTable=st,StateDistribution=full,Diagnostics=c("Exact one-node MLU joint pathogen/information validator with a single node-level HaveInfo state shared across land uses.","Management adoption is drawn independently by land use conditional on shared HaveInfo; actual management deaths are host evidence only when host acquisition is enabled.",if(finite)"Programmed stop is represented by explicit node-level time-since-evidence states."else"Information loss follows the memoryless retention pathway."))
  if(Model=="SEIR")out$ExpectedExposedByLandUse<-means$E;if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMLUPathogenExactInformationOneNode","list");out
}

###############################################################################
### Exact one-node MLU SIR / SEIR gross-turnover extension
###############################################################################

.ina_mlu_compartment_turnover_operator_one_node <- function(Model,K,Survival=1,RecruitToCapacity=FALSE,
    Beta=0,RecoveryProb=0,ProgressionProb=1,ImmunityLossProb=0,PathogenMortalityProb=0,
    IntroductionProb=0,IntroductionNumber=1,LandUseMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,MaxStates=50000L) {
  Model<-match.arg(Model,c("SIR","SEIR"));Transmission<-match.arg(Transmission);K<-as.integer(K);L<-length(K);states<-.ina_mlu_model_states(Model)
  local<-lapply(K,function(k).ina_meta_compartment_states_uptoK(Model,k));dims<-vapply(local,nrow,integer(1));nst<-prod(dims);if(nst>MaxStates)stop("Exact MLU turnover state space has ",nst," states")
  g<-as.matrix(expand.grid(lapply(dims,seq_len),KEEP.OUT.ATTRS=FALSE));storage.mode(g)<-"integer";st<-matrix(0L,nst,length(states)*L);cn<-character(ncol(st));for(lu in seq_len(L))for(j in seq_along(states)){cc<-(lu-1L)*length(states)+j;st[,cc]<-local[[lu]][[states[j]]][g[,lu]];cn[cc]<-paste0(states[j],"_lu",lu)};colnames(st)<-cn;keys<-apply(st,1,paste,collapse=":");idx<-setNames(seq_len(nst),keys);T<-matrix(0,nst,nst,dimnames=list(keys,keys))
  res<-function(x,nm,int=FALSE).ina_mlu_resolve_unit(x,1L,L,nm,int);surv<-res(Survival,"Survival");beta<-res(Beta,"Beta");rec<-res(RecoveryProb,"RecoveryProb");prog<-res(ProgressionProb,"ProgressionProb");wan<-res(ImmunityLossProb,"ImmunityLossProb");pmort<-res(PathogenMortalityProb,"PathogenMortalityProb");ip<-res(IntroductionProb,"IntroductionProb");inum<-as.integer(res(IntroductionNumber,"IntroductionNumber",TRUE));ds<-res(DensityScale,"DensityScale");C<-.ina_mlu_contact(1L,L,NodeContact=matrix(1,1,1),LandUseMixing=if(is.null(LandUseMixing))diag(L)else LandUseMixing)
  for(row in seq_len(nst)){hl<-vector("list",L);for(lu in seq_len(L)){cc<-setNames(vapply(states,function(nm)as.integer(st[row,paste0(nm,"_lu",lu)]),integer(1)),states);hl[[lu]]<-.ina_meta_compartment_host_outcomes(cc,Model,K[lu],surv[lu],RecruitToCapacity,0,0,NULL)};hg<-as.matrix(expand.grid(lapply(hl,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(hg)<-"integer"
    for(hb in seq_len(nrow(hg))){pre<-vector("list",L);ph<-1;Ipre<-Npre<-numeric(L);for(lu in seq_len(L)){z<-hl[[lu]][hg[hb,lu],];pre[[lu]]<-z;ph<-ph*z$prob;Ipre[lu]<-z$I;Npre[lu]<-sum(as.numeric(z[states]))};pressure<-as.numeric(crossprod(Ipre,C));den<-as.numeric(crossprod(Npre,C));foi<-if(Transmission=="frequency")beta*ifelse(den>0,pressure/den,0)else beta*pressure/ds;pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)));dl<-vector("list",L);for(lu in seq_len(L))dl[[lu]]<-.ina_mlu_compartment_disease_outcomes_pinf(pre[[lu]],Model,pinf[lu],rec[lu],prog[lu],wan[lu],pmort[lu],ip[lu],inum[lu]);dg<-as.matrix(expand.grid(lapply(dl,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(dg)<-"integer"
      for(dd in seq_len(nrow(dg))){dest<-integer(length(cn));names(dest)<-cn;p<-ph;for(lu in seq_len(L)){z<-dl[[lu]][dg[dd,lu],];p<-p*z$prob;for(nm in states)dest[paste0(nm,"_lu",lu)]<-as.integer(z[[nm]])};ky<-paste(dest,collapse=":");T[row,idx[[ky]]]<-T[row,idx[[ky]]]+p}
    }
  }
  er<-max(abs(rowSums(T)-1));if(er>1e-10)stop(Model," MLU turnover operator row error ",er);attr(T,"States")<-st;T
}

INApestMLUCompartmentExactTurnoverOneNode <- function(Model=c("SIR","SEIR"),Ntimesteps,K,InitialHostPopulation,InitialInfected=0,InitialExposed=0,InitialRecovered=0,Survival=1,RecruitToCapacity=FALSE,Beta=0,RecoveryProb=0,ProgressionProb=1,ImmunityLossProb=0,PathogenMortalityProb=0,IntroductionProb=0,IntroductionNumber=1,LandUseMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,MaxStates=50000L,ReturnOperator=FALSE){
  Model<-match.arg(Model);Transmission<-match.arg(Transmission);K<-as.integer(K);L<-length(K);states<-.ina_mlu_model_states(Model);resi<-function(x,nm).ina_mlu_resolve_unit(x,1L,L,nm,TRUE);N0<-as.integer(resi(InitialHostPopulation,"InitialHostPopulation"));I0<-as.integer(resi(InitialInfected,"InitialInfected"));E0<-if(Model=="SEIR")as.integer(resi(InitialExposed,"InitialExposed"))else integer(L);R0<-as.integer(resi(InitialRecovered,"InitialRecovered"));if(any(N0>K)||any(I0+E0+R0>N0))stop("Invalid initial counts")
  T<-.ina_mlu_compartment_turnover_operator_one_node(Model,K,Survival,RecruitToCapacity,Beta,RecoveryProb,ProgressionProb,ImmunityLossProb,PathogenMortalityProb,IntroductionProb,IntroductionNumber,LandUseMixing,Transmission,DensityScale,MaxStates);st<-attr(T,"States");keys<-rownames(T);init<-integer(length(states)*L);names(init)<-colnames(st);for(lu in seq_len(L)){init[paste0("S_lu",lu)]<-N0[lu]-I0[lu]-E0[lu]-R0[lu];init[paste0("I_lu",lu)]<-I0[lu];if(Model=="SEIR")init[paste0("E_lu",lu)]<-E0[lu];init[paste0("R_lu",lu)]<-R0[lu]};dist<-numeric(nrow(st));dist[match(paste(init,collapse=":"),keys)]<-1;full<-matrix(0,Ntimesteps+1L,nrow(st));full[1,]<-dist;means<-lapply(states,function(nm)matrix(0,Ntimesteps+1L,L));names(means)<-states;for(nm in states)for(lu in seq_len(L))means[[nm]][1,lu]<-init[paste0(nm,"_lu",lu)];active<-numeric(Ntimesteps+1L);active[1]<-as.numeric(any(I0+E0>0))
  for(tt in seq_len(Ntimesteps)){dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;for(nm in states)for(lu in seq_len(L))means[[nm]][tt+1,lu]<-sum(dist*st[,paste0(nm,"_lu",lu)]);active[tt+1]<-if(Model=="SEIR")sum(dist[rowSums(st[,grep("^(E|I)_lu",colnames(st)),drop=FALSE])>0])else sum(dist[rowSums(st[,grep("^I_lu",colnames(st)),drop=FALSE])>0])};EN<-matrix(0,Ntimesteps+1L,L);for(nm in states)EN<-EN+means[[nm]];out<-list(Model="INApestMetaMultipleLandUse",PathogenModel=Model,HostAssumption="explicit one-node MLU turnover",Exact=TRUE,ExpectedHostByLandUse=EN,ExpectedStatesByLandUse=means,ExpectedInfectedByLandUse=means$I,ActivePathogenProbability=active,PathogenFreedomProbability=1-active,StateTable=st,StateDistribution=full,Diagnostics=c("Exact one-node MLU compartment turnover validator.","Gross host survival thins pathogen states before optional susceptible replacement recruitment.","Pathogen transmission then occurs across land uses using the supplied land-use mixing matrix."));if(Model=="SEIR")out$ExpectedExposedByLandUse<-means$E;if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMLUCompartmentExactTurnoverOneNode","list");out
}

###############################################################################
### MLU pathogen analytical router and unified dispatcher extension
###############################################################################

.ina_mlu_initial_matrix <- function(x,n_nodes,n_landuses,name){
  if(is.matrix(x)){if(!all(dim(x)==c(n_nodes,n_landuses)))stop(name," must be nodes x land uses");z<-x}else{v<-as.numeric(x);if(length(v)==1L)z<-matrix(v,n_nodes,n_landuses)else if(length(v)==n_nodes*n_landuses)z<-matrix(v,n_nodes,n_landuses)else if(length(v)==n_landuses)z<-matrix(rep(v,each=n_nodes),n_nodes,n_landuses)else stop(name," has incompatible dimensions")};if(any(!is.finite(z))||any(z<0)||any(z!=floor(z)))stop(name," must contain non-negative whole numbers");storage.mode(z)<-"integer";z
}

INApestMLUPathogenAnalytical <- function(Ntimesteps,HostPopulation,Pathogen,
    HostAssumption=c("fixed","event_contract_turnover"),K=HostPopulation,Survival=1,RecruitToCapacity=FALSE,
    HostDetectionProb=0,InitialInfo=0,ManageProb=0,MortalityProb=0,InformationAcquisition=NULL,
    InfoPersistenceSteps=NA,InfoRetentionProb=1,NodeContact=NULL,LandUseMixing=NULL,
    OutsideNodes=integer(0),OutsideUnits=integer(0),Exact=TRUE,ExactMaxStates=50000L,ReturnOperators=FALSE){
  if(!inherits(Pathogen,"INApestPathogen"))stop("Pathogen must be returned by INApestPathogen()");Model<-as.character(Pathogen$Model)[1L];if(!Model%in%c("SIS","SIR","SEIR"))stop("MLU analytical router supports SIS, SIR and SEIR");HostAssumption<-match.arg(HostAssumption);N<-.ina_mlu_check_matrix(HostPopulation,name="HostPopulation",integer=TRUE);n_nodes<-nrow(N);L<-ncol(N);Kmat<-if(is.matrix(K))K else matrix(.ina_mlu_resolve_unit(K,n_nodes,L,"K",TRUE),n_nodes,L);if(!all(dim(Kmat)==dim(N))||any(Kmat<N))stop("K must match HostPopulation dimensions and be >= HostPopulation")
  I0<-.ina_mlu_initial_matrix(Pathogen$InitialInfected,n_nodes,L,"InitialInfected");E0<-if(Model=="SEIR").ina_mlu_initial_matrix(Pathogen$InitialExposed,n_nodes,L,"InitialExposed")else matrix(0L,n_nodes,L);R0<-if(Model%in%c("SIR","SEIR")).ina_mlu_initial_matrix(Pathogen$InitialRecovered,n_nodes,L,"InitialRecovered")else matrix(0L,n_nodes,L);if(any(I0+E0+R0>N))stop("Initial pathogen states exceed HostPopulation")
  acq<-InformationAcquisition;if(is.null(acq))acq<-if(isTRUE(Pathogen$DetectionTriggersInfo))"both"else"host";acq<-match.arg(as.character(acq),c("host","pathogen","both"));info_requested<-any(as.numeric(ManageProb)!=0,na.rm=TRUE)||any(as.numeric(InitialInfo)!=0,na.rm=TRUE)||any(as.numeric(HostDetectionProb)!=0,na.rm=TRUE)||any(as.numeric(Pathogen$DetectionProb)!=0,na.rm=TRUE)||any(!is.na(InfoPersistenceSteps))||any(as.numeric(InfoRetentionProb)!=1,na.rm=TRUE)
  C<-.ina_mlu_contact(n_nodes,L,ContactMatrix=Pathogen$ContactMatrix,NodeContact=NodeContact,LandUseMixing=LandUseMixing)
  if(!isTRUE(Exact)){g<-INApestMLUPathogenGrowthOperator(Model,N,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$PathogenMortalityProb,ContactMatrix=C,Transmission=Pathogen$Transmission,DensityScale=Pathogen$DensityScale);return(structure(list(Model="INApestMetaMultipleLandUse",PathogenModel=Model,Exact=FALSE,Growth=g,Diagnostics=c("Scalable MLU result is the rare-pathogen growth operator on the flattened node x land-use graph.","Use stochastic MLU for finite-prevalence managed landscapes once exact state expansion is impractical.")),class=c("INApestMLUPathogenAnalytical","list")))}
  if(HostAssumption=="event_contract_turnover"){
    if(n_nodes!=1L)stop("Exact MLU event-contract turnover is currently a one-node validation solution")
    if(info_requested)stop("MLU turnover and shared-information validators are exact but exposed separately; combined turnover+information is not routed generically")
    if(Model=="SIS")return(INApestMLUSISExactTurnoverOneNode(Ntimesteps,as.integer(Kmat[1,]),as.integer(N[1,]),as.integer(I0[1,]),Survival,RecruitToCapacity,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,if(is.null(LandUseMixing))diag(L)else LandUseMixing,Pathogen$Transmission,Pathogen$DensityScale,ReturnOperators))
    return(INApestMLUCompartmentExactTurnoverOneNode(Model,Ntimesteps,as.integer(Kmat[1,]),as.integer(N[1,]),as.integer(I0[1,]),as.integer(E0[1,]),as.integer(R0[1,]),Survival,RecruitToCapacity,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,if(is.null(LandUseMixing))diag(L)else LandUseMixing,Pathogen$Transmission,Pathogen$DensityScale,ExactMaxStates,ReturnOperators))
  }
  if(info_requested){
    if(n_nodes!=1L)stop("Exact shared-information MLU expansion is currently exposed for one-node validation systems; use stochastic MLU for larger managed networks")
    return(INApestMLUPathogenExactInformationOneNode(Model,Ntimesteps,as.integer(Kmat[1,]),as.integer(N[1,]),as.integer(I0[1,]),as.integer(E0[1,]),as.integer(R0[1,]),as.integer(InitialInfo)[1],Survival,ManageProb,MortalityProb,HostDetectionProb,Pathogen$DetectionProb,acq,InfoPersistenceSteps,InfoRetentionProb,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,if(is.null(LandUseMixing))diag(L)else LandUseMixing,Pathogen$Transmission,Pathogen$DensityScale,ExactMaxStates,ReturnOperators))
  }
  if(any(as.numeric(Pathogen$PathogenMortalityProb)!=0))stop("Fixed-host exact MLU branch requires PathogenMortalityProb=0; use HostAssumption='event_contract_turnover'")
  INApestMLUPathogenExactFixedN(Model,Ntimesteps,N,I0,E0,R0,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$ProgressionProb,Pathogen$ImmunityLossProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,ContactMatrix=C,Transmission=Pathogen$Transmission,DensityScale=Pathogen$DensityScale,OutsideNodes=OutsideNodes,OutsideUnits=OutsideUnits,MaxStates=ExactMaxStates,ReturnOperator=ReturnOperators)
}

print.INApestMLUPathogenAnalytical <- function(x,...){cat("INApest MLU pathogen analytical result\n");cat("  Pathogen model:",x$PathogenModel,"\n");cat("  Exact:",isTRUE(x$Exact),"\n");if(!is.null(x$Growth$IntrinsicRarePathogenMultiplier))cat("  Rare-pathogen multiplier:",round(x$Growth$IntrinsicRarePathogenMultiplier,6),"\n");invisible(x)}

INApestAnalytical_pre_mlu_pathogen <- INApestAnalytical
INApestAnalytical <- function(...) {
  args<-list(...);Model<-if(!is.null(args$Model))as.character(args$Model)[1L]else"INApest";Pathogen<-args$Pathogen
  if(!identical(Model,"INApestMetaMultipleLandUse")||is.null(Pathogen))return(do.call(INApestAnalytical_pre_mlu_pathogen,args))
  if(is.null(args$HostPopulation)){if(is.null(args$InitialState))stop("MLU pathogen analytical calls require HostPopulation or InitialState");args$HostPopulation<-args$InitialState}
  keep<-c("Ntimesteps","HostPopulation","Pathogen","HostAssumption","K","Survival","RecruitToCapacity","HostDetectionProb","InitialInfo","ManageProb","MortalityProb","InformationAcquisition","InfoPersistenceSteps","InfoRetentionProb","NodeContact","LandUseMixing","OutsideNodes","OutsideUnits","Exact","ExactMaxStates","ReturnOperators")
  do.call(INApestMLUPathogenAnalytical,args[intersect(names(args),keep)])
}

###############################################################################
### Exact tiny-network MLU SIS + shared node information + managed escape
###############################################################################

.ina_mlu_node_vector <- function(x,n_nodes,name,integer=FALSE){z<-as.numeric(x);if(length(z)==1L)z<-rep(z,n_nodes);if(length(z)!=n_nodes||any(!is.finite(z)))stop(name," must be scalar or length nodes");if(integer&&any(z!=floor(z)))stop(name," must contain whole numbers");z}

INApestMLUSISExactNetworkInformation <- function(Ntimesteps,K,InitialHostPopulation,InitialInfected,InitialInfo=0,
    Survival=1,ManageProb=0,MortalityProb=0,HostDetectionProb=0,PathogenDetectionProb=0,
    InformationAcquisition=c("host","pathogen","both"),InfoPersistenceSteps=NA,InfoRetentionProb=1,
    Beta=0,RecoveryProb=0,PathogenMortalityProb=0,IntroductionProb=0,IntroductionNumber=1,
    ContactMatrix=NULL,NodeContact=NULL,LandUseMixing=NULL,Transmission=c("frequency","density"),DensityScale=1,
    OutsideNodes=integer(0),MaxStates=5000L,ReturnOperator=FALSE){
  InformationAcquisition<-match.arg(InformationAcquisition);Transmission<-match.arg(Transmission)
  K<-.ina_mlu_check_matrix(K,name="K",integer=TRUE);N0<-.ina_mlu_check_matrix(InitialHostPopulation,nrow(K),ncol(K),"InitialHostPopulation",TRUE);I0<-.ina_mlu_check_matrix(InitialInfected,nrow(K),ncol(K),"InitialInfected",TRUE);if(any(N0>K)||any(I0>N0))stop("Invalid initial network state")
  n<-nrow(K);L<-ncol(K);n_units<-n*L;H0<-as.integer(.ina_mlu_node_vector(InitialInfo,n,"InitialInfo",TRUE));if(any(!H0%in%c(0L,1L)))stop("InitialInfo must be 0/1")
  ipers<-as.numeric(InfoPersistenceSteps);if(length(ipers)==1L)ipers<-rep(ipers,n);if(length(ipers)!=n||any(!is.na(ipers)&(ipers<0|ipers!=floor(ipers))))stop("InfoPersistenceSteps must be scalar or length nodes")
  iret<-.ina_mlu_node_vector(InfoRetentionProb,n,"InfoRetentionProb");if(any(iret<0|iret>1))stop("InfoRetentionProb must be in [0,1]")
  to_mat<-function(x,nm,int=FALSE){matrix(.ina_mlu_resolve_unit(x,n,L,nm,int),n,L)}
  surv<-to_mat(Survival,"Survival");mng<-to_mat(ManageProb,"ManageProb");mort<-to_mat(MortalityProb,"MortalityProb");hd<-to_mat(HostDetectionProb,"HostDetectionProb");pd<-to_mat(PathogenDetectionProb,"PathogenDetectionProb");beta<-as.numeric(.ina_mlu_resolve_unit(Beta,n,L,"Beta"));rec<-as.numeric(.ina_mlu_resolve_unit(RecoveryProb,n,L,"RecoveryProb"));pmort<-as.numeric(.ina_mlu_resolve_unit(PathogenMortalityProb,n,L,"PathogenMortalityProb"));intro<-as.numeric(.ina_mlu_resolve_unit(IntroductionProb,n,L,"IntroductionProb"));inum<-as.integer(.ina_mlu_resolve_unit(IntroductionNumber,n,L,"IntroductionNumber",TRUE));ds<-as.numeric(.ina_mlu_resolve_unit(DensityScale,n,L,"DensityScale"))
  if(any(c(surv,mng,mort,hd,pd,rec,pmort,intro)<0)||any(c(surv,mng,mort,hd,pd,rec,pmort,intro)>1)||any(rec+pmort>1+1e-12)||any(beta<0)||any(ds<=0))stop("Invalid network parameters")
  C<-.ina_mlu_contact(n,L,ContactMatrix,NodeContact,LandUseMixing);host_acq<-InformationAcquisition%in%c("host","both");path_acq<-InformationAcquisition%in%c("pathogen","both")
  local<-lapply(seq_len(n),function(j).ina_mlu_one_node_info_states("SIS",as.integer(K[j,]),ipers[j],MaxStates));dims<-vapply(local,nrow,integer(1));nst<-prod(dims);if(nst>MaxStates)stop("Exact MLU network information state space has ",nst," states; reduce capacities/nodes or increase MaxStates deliberately")
  grids<-as.matrix(expand.grid(lapply(dims,seq_len),KEEP.OUT.ATTRS=FALSE));storage.mode(grids)<-"integer";gkey<-apply(grids,1,paste,collapse=":");gidx<-setNames(seq_len(nst),gkey);local_keys<-lapply(local,function(z)apply(z,1,paste,collapse=":"));local_idx<-lapply(local,function(z)setNames(seq_len(nrow(z)),apply(z,1,paste,collapse=":")))
  T<-matrix(0,nst,nst)
  for(row in seq_len(nst)){
    preL<-vector("list",n)
    for(j in seq_len(n))preL[[j]]<-.ina_mlu_one_node_info_prebranches(local[[j]][grids[row,j],,drop=FALSE],"SIS",L,surv[j,],mng[j,],mort[j,],host_acq,ipers[j],iret[j])
    pg<-as.matrix(expand.grid(lapply(preL,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(pg)<-"integer"
    for(pb in seq_len(nrow(pg))){Npre<-Ipre<-matrix(0,n,L);Hmid<-Amid<-integer(n);pp<-1
      for(j in seq_len(n)){z<-preL[[j]][pg[pb,j],];for(lu in seq_len(L)){Ipre[j,lu]<-z[[paste0("I_lu",lu)]];Npre[j,lu]<-Ipre[j,lu]+z[[paste0("S_lu",lu)]]};Hmid[j]<-z$H;Amid[j]<-z$A;pp<-pp*z$prob};if(pp==0)next
      Iv<-as.numeric(Ipre);Nv<-as.numeric(Npre);pressure<-as.numeric(crossprod(Iv,C));den<-as.numeric(crossprod(Nv,C));foi<-if(Transmission=="frequency")beta*ifelse(den>0,pressure/den,0)else beta*pressure/ds;pinf<-.ina_meta_sis_clip01(-expm1(-pmax(0,foi)))
      dl<-vector("list",n_units);for(u in seq_len(n_units))dl[[u]]<-as.data.frame(.ina_meta_sis_disease_outcomes_pinf(as.integer(Nv[u]),as.integer(Iv[u]),pinf[u],rec[u],pmort[u],intro[u],inum[u]));dg<-as.matrix(expand.grid(lapply(dl,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(dg)<-"integer"
      for(dd in seq_len(nrow(dg))){N1<-I1<-numeric(n_units);pdis<-pp;for(u in seq_len(n_units)){z<-dl[[u]][dg[dd,u],];N1[u]<-z$N;I1[u]<-z$I;pdis<-pdis*z$prob};if(pdis==0)next;N1m<-matrix(N1,n,L);I1m<-matrix(I1,n,L)
        db<-vector("list",n);for(j in seq_len(n)){qh<-if(host_acq)1-prod((1-hd[j,])^N1m[j,])else 0;qp<-if(path_acq)1-prod((1-pd[j,])^I1m[j,])else 0;qev<-1-(1-qh)*(1-qp);if(qev>0&&qev<1)db[[j]]<-rbind(c(H=1,A=if(is.na(ipers[j]))-2 else 0,p=qev),c(H=Hmid[j],A=Amid[j],p=1-qev))else if(qev>=1)db[[j]]<-matrix(c(1,if(is.na(ipers[j]))-2 else 0,1),1,dimnames=list(NULL,c("H","A","p")))else db[[j]]<-matrix(c(Hmid[j],Amid[j],1),1,dimnames=list(NULL,c("H","A","p")))}
        bg<-as.matrix(expand.grid(lapply(db,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(bg)<-"integer"
        for(bb in seq_len(nrow(bg))){li<-integer(n);pfin<-pdis;for(j in seq_len(n)){q<-db[[j]][bg[bb,j],];bio<-integer(2L*L);for(lu in seq_len(L)){u<-j+(lu-1L)*n;bio[2L*lu-1L]<-as.integer(N1[u]-I1[u]);bio[2L*lu]<-as.integer(I1[u])};ky<-paste(c(bio,as.integer(q["H"]),as.integer(q["A"])),collapse=":");li[j]<-local_idx[[j]][[ky]];pfin<-pfin*unname(q["p"])};gi<-gidx[[paste(li,collapse=":")]];T[row,gi]<-T[row,gi]+pfin}
      }
    }
  }
  er<-max(abs(rowSums(T)-1));if(er>1e-10)stop("MLU network information operator row error ",er)
  # Initial distribution includes direct local evidence before timestep 1.
  local_init<-vector("list",n)
  for(j in seq_len(n)){qh<-if(host_acq)1-prod((1-hd[j,])^N0[j,])else 0;qp<-if(path_acq)1-prod((1-pd[j,])^I0[j,])else 0;qev<-1-(1-qh)*(1-qp);bi<-integer(2L*L);for(lu in seq_len(L)){bi[2L*lu-1L]<-N0[j,lu]-I0[j,lu];bi[2L*lu]<-I0[j,lu]};lst<-list();kk<-0L;if(qev>0){kk<-kk+1L;lst[[kk]]<-c(idx=local_idx[[j]][[paste(c(bi,1,if(is.na(ipers[j]))-2 else 0),collapse=":")]],p=qev)};if(qev<1){A<-if(H0[j]==1&&!is.na(ipers[j]))-1 else -2;kk<-kk+1L;lst[[kk]]<-c(idx=local_idx[[j]][[paste(c(bi,H0[j],A),collapse=":")]],p=1-qev)};local_init[[j]]<-do.call(rbind,lst)}
  ig<-as.matrix(expand.grid(lapply(local_init,function(z)seq_len(nrow(z))),KEEP.OUT.ATTRS=FALSE));storage.mode(ig)<-"integer";dist<-numeric(nst);for(rr in seq_len(nrow(ig))){li<-integer(n);p<-1;for(j in seq_len(n)){z<-local_init[[j]][ig[rr,j],];li[j]<-as.integer(z["idx"]);p<-p*z["p"]};dist[gidx[[paste(li,collapse=":")]]]<-dist[gidx[[paste(li,collapse=":")]]]+p}
  EI<-array(0,dim=c(n,L,Ntimesteps+1L));EN<-array(0,dim=c(n,L,Ntimesteps+1L));PH<-matrix(0,Ntimesteps+1L,n);PM<-array(0,dim=c(n,L,Ntimesteps));full<-matrix(0,Ntimesteps+1L,nst);full[1,]<-dist
  summarize<-function(d,tt){for(j in seq_len(n)){h<-numeric(nst);for(r in seq_len(nst)){lr<-local[[j]][grids[r,j],];h[r]<-lr["H"];for(lu in seq_len(L)){EI[j,lu,tt]<<-EI[j,lu,tt]+d[r]*lr[paste0("I_lu",lu)];EN[j,lu,tt]<<-EN[j,lu,tt]+d[r]*(lr[paste0("S_lu",lu)]+lr[paste0("I_lu",lu)])}};PH[tt,j]<<-sum(d*h)}}
  summarize(dist,1L);for(tt in seq_len(Ntimesteps)){for(j in seq_len(n))for(lu in seq_len(L))PM[j,lu,tt]<-mng[j,lu]*PH[tt,j];dist<-as.numeric(dist%*%T);full[tt+1,]<-dist;summarize(dist,tt+1L)}
  outside<-unique(as.integer(OutsideNodes));if(any(!outside%in%seq_len(n)))stop("OutsideNodes contains invalid node indices");escape<-NULL
  if(length(outside)){
    safe<-logical(nst);for(r in seq_len(nst)){ok<-TRUE;for(j in outside){lr<-local[[j]][grids[r,j],];if(sum(lr[grep("^I_lu",names(lr))])>0){ok<-FALSE;break}};safe[r]<-ok};sd<-full[1,safe];Q<-T[safe,safe,drop=FALSE];esc<-numeric(Ntimesteps+1L);esc[1]<-1-sum(sd);if(Ntimesteps>0)for(tt in seq_len(Ntimesteps)){sd<-as.numeric(sd%*%Q);esc[tt+1]<-1-sum(sd)};escape<-list(OutsideNodes=outside,ProbabilityByTimestep=esc,ProbabilityByHorizon=tail(esc,1))
  }
  out<-list(Model="INApestMetaMultipleLandUse",PathogenModel="SIS",Exact=TRUE,InformationAcquisition=InformationAcquisition,ExpectedHost=EN,ExpectedInfected=EI,InformationProbabilityByNode=PH,ManagingProbability=PM,Escape=escape,StateCount=nst,StateDistribution=full,Diagnostics=c("Exact tiny-network MLU SIS operator with land-use host/pathogen states and node-level shared information.","First-passage escape is computed with a safe-state sub-operator, so later clearance outside containment does not erase an earlier escape.","This branch is intentionally state-space limited and serves as validation truth; use stochastic MLU for realistic managed landscapes."));if(ReturnOperator)out$Operator<-T;class(out)<-c("INApestMLUSISExactNetworkInformation","list");out
}

# Replace only the MLU router with the final network-information-aware version.
INApestMLUPathogenAnalytical_v1 <- INApestMLUPathogenAnalytical
INApestMLUPathogenAnalytical <- function(Ntimesteps,HostPopulation,Pathogen,
    HostAssumption=c("fixed","event_contract_turnover"),K=HostPopulation,Survival=1,RecruitToCapacity=FALSE,
    HostDetectionProb=0,InitialInfo=0,ManageProb=0,MortalityProb=0,InformationAcquisition=NULL,
    InfoPersistenceSteps=NA,InfoRetentionProb=1,NodeContact=NULL,LandUseMixing=NULL,
    OutsideNodes=integer(0),OutsideUnits=integer(0),Exact=TRUE,ExactMaxStates=50000L,ReturnOperators=FALSE){
  N<-.ina_mlu_check_matrix(HostPopulation,name="HostPopulation",integer=TRUE);n<-nrow(N);L<-ncol(N);Model<-as.character(Pathogen$Model)[1L];acq<-InformationAcquisition;if(is.null(acq))acq<-if(isTRUE(Pathogen$DetectionTriggersInfo))"both"else"host";info_requested<-any(as.numeric(ManageProb)!=0,na.rm=TRUE)||any(as.numeric(InitialInfo)!=0,na.rm=TRUE)||any(as.numeric(HostDetectionProb)!=0,na.rm=TRUE)||any(as.numeric(Pathogen$DetectionProb)!=0,na.rm=TRUE)||any(!is.na(InfoPersistenceSteps))||any(as.numeric(InfoRetentionProb)!=1,na.rm=TRUE)
  if(isTRUE(Exact)&&match.arg(HostAssumption)=="fixed"&&info_requested&&Model=="SIS"&&n>1L){Kmat<-if(is.matrix(K))K else matrix(.ina_mlu_resolve_unit(K,n,L,"K",TRUE),n,L);I0<-.ina_mlu_initial_matrix(Pathogen$InitialInfected,n,L,"InitialInfected");C<-.ina_mlu_contact(n,L,ContactMatrix=Pathogen$ContactMatrix,NodeContact=NodeContact,LandUseMixing=LandUseMixing);return(INApestMLUSISExactNetworkInformation(Ntimesteps,Kmat,N,I0,InitialInfo,Survival,ManageProb,MortalityProb,HostDetectionProb,Pathogen$DetectionProb,acq,InfoPersistenceSteps,InfoRetentionProb,Pathogen$Beta,Pathogen$RecoveryProb,Pathogen$PathogenMortalityProb,Pathogen$IntroductionProb,Pathogen$IntroductionNumber,ContactMatrix=C,Transmission=Pathogen$Transmission,DensityScale=Pathogen$DensityScale,OutsideNodes=OutsideNodes,MaxStates=ExactMaxStates,ReturnOperator=ReturnOperators))}
  INApestMLUPathogenAnalytical_v1(Ntimesteps,HostPopulation,Pathogen,HostAssumption,K,Survival,RecruitToCapacity,HostDetectionProb,InitialInfo,ManageProb,MortalityProb,InformationAcquisition,InfoPersistenceSteps,InfoRetentionProb,NodeContact,LandUseMixing,OutsideNodes,OutsideUnits,Exact,ExactMaxStates,ReturnOperators)
}
