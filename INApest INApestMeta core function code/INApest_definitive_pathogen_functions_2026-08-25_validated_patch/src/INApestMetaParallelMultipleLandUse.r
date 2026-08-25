###########################################################################
###########################################################################
###Declares a function overlaying management on a metapopulation spread model 
###Key inputs are: 
###1) matrix of natural dispersal probibilities between each pair of sites (i.e. nodes of the network). Matrix can be non-symmetrical (i.e. can have source and sink nodes) 
###2) Envionmentally-determined per-capita propagule production
###3) Envionmentally-determined carrying capacity (K)
###4) Management parameters
### a) Annual detection probability
### b) Annual management adoption probability subsequent to detection
### c) Annual mortality probability under management
### d) Spread reduction when management adopted
###5) A binary vector of initially infested nodes OR
###   A vector of invasion risk probabilities and proportion of nodes initially invaded
###Key outputs are:
###3-dimensional arrays of invasion, management and detection status for each node in each timestep of each permuation
###2-dimensional array of invasion probability (i.e. proportion of permutations pest present) for each node in each timestep
###Line graphs summarising number of total population (as a proportion of carrying capacity) nodes infested, where infestations
###are detected, the proportion of infestations detected and nodes under management against time 
###########################################################################
###########################################################################

#######################################################################
###This version implements parallel processing with a PSOCK cluster and parallel::parLapply for the permutation loop.
###See file "ParallelSetup.r" for notes on steps for setting up parallel processing
###Allows for multiple land uses within a location, each with different management parameters 
###and carrying capacity
#######################################################################



local.dynamicsLU = function(sddprob = SDDprob, nodepropaguleproduction = NodePropaguleProduction,nodeenvestabprob = NodeEnvEstabProb,n=N0,
lddprob = LDDprob, lddrate = LDDrate,k_is_0 = K_is_0, nodeK = NodeK,nodepropaguleestablishment = NodePropaguleEstablishment,
nodespreadreduction = NodeSpreadReduction,nodefecundityreduction = NodeFecundityReduction,managing = Managing)
{
  # Keep the default LocalDynamics self-contained so it can be serialized to
  # PSOCK workers without relying on objects in the master's global environment.
  AllocateLandUseRecruits <- function(X)
    {
    Vect = vector(length=0)
    for(i in 2:length(X))
      Vect <- c(Vect,rep((i-1),times = X[i]))
    Sample <- sample(Vect,size = X[1],replace = F)
    Out <- vector(length = length(X)-1)
    for(j in 1:length(Out))
      Out[j] <- length(Sample[Sample==j])
    return(Out)
    }
FecundityContrib <- n * (1-nodefecundityreduction*managing)
Propagules <- rpois(nrow(sddprob), nodepropaguleproduction * rowSums(FecundityContrib))# propagules are produced after management effects on fecundity
      ###self-mediated spread
      Pin = 0
      Qin = 0
      Pout <- Propagules*(1-lddrate)
      if(sum(Pout)>0 ) 
       Pin <- t(rmultinom(1, size=sum(Pout*rowSums(sddprob)), prob=Pout %*% sddprob))  # propagules are dispersed
      
      ###human-mediated spread
      if (is.matrix(lddprob)==T) 
       {
       Pn = FecundityContrib/rowSums(FecundityContrib)
       Pn[is.na(Pn)==TRUE] = 0
       Qout  = rowSums(Propagules*lddrate*Pn *(1-nodespreadreduction*managing))       
       if(sum(Qout)>0)  
        Qin <- t(rmultinom(1, size=sum(Qout*rowSums(lddprob)), prob=Qout %*% lddprob))    # propagules are dispersed
       }
     
    # propagule success depends on availability of uninfested host plants
    Recruits <- rbinom(nrow(sddprob), rowSums(nodeK)-rowSums(n), 1 - exp(-nodepropaguleestablishment*nodeenvestabprob*(Pin+Qin)))  
    InVector = cbind(Recruits,nodeK-n)
    LUrecruits <- t(apply(InVector,1,FUN = AllocateLandUseRecruits))
    ###Possibly due to rounding error when K is large N0 can occasionally = K+1
    Nout <- ifelse(n + LUrecruits>nodeK,nodeK,n + LUrecruits)
return(Nout)
}


# -----------------------------------------------------------------------------
# Optional custom LocalDynamics arguments
# -----------------------------------------------------------------------------
# Ordinary LocalDynamicsArgs entries are passed unchanged. Wrap a vector,
# matrix, array, or list in INApestLocalDynamicsTimeArg() when its final
# dimension/list position is indexed by simulation timestep. This avoids
# guessing whether an arbitrary custom matrix is static or time-varying.
INApestLocalDynamicsTimeArg <- function(x) {
  structure(list(values = x), class = "INApestLocalDynamicsTimeArg")
}

.resolve_INApest_LocalDynamicsArgs <- function(LocalDynamicsArgs, timestep, Ntimesteps) {
  if (is.null(LocalDynamicsArgs)) LocalDynamicsArgs <- list()
  if (!is.list(LocalDynamicsArgs))
    stop("LocalDynamicsArgs must be a named list")
  if (!length(LocalDynamicsArgs)) return(list())
  if (is.null(names(LocalDynamicsArgs)) || any(!nzchar(names(LocalDynamicsArgs))))
    stop("Every LocalDynamicsArgs entry must have a non-empty name")
  if (anyDuplicated(names(LocalDynamicsArgs)))
    stop("LocalDynamicsArgs names must be unique")

  resolve_one <- function(x, name) {
    if (inherits(x, "INApestLocalDynamicsTimeArg")) {
      values <- x$values
      if (is.list(values)) {
        if (length(values) != Ntimesteps)
          stop(name, " wrapped list must have length Ntimesteps")
        return(values[[timestep]])
      }
      d <- dim(values)
      if (is.null(d)) {
        if (length(values) != Ntimesteps)
          stop(name, " wrapped vector must have length Ntimesteps")
        return(values[[timestep]])
      }
      if (tail(d, 1L) != Ntimesteps)
        stop(name, " wrapped array/matrix must have Ntimesteps in its final dimension")
      index <- lapply(d, seq_len)
      index[[length(d)]] <- timestep
      return(do.call(`[`, c(list(values), index, list(drop = TRUE))))
    }

    # Resolver functions are evaluated by the parent model. The returned
    # current value, not timestep itself, is passed to LocalDynamics.
    if (is.function(x)) {
      fm <- names(formals(x))
      call_args <- list(timestep = timestep, Ntimesteps = Ntimesteps)
      if (!is.null(fm) && !("..." %in% fm))
        call_args <- call_args[intersect(names(call_args), fm)]
      return(do.call(x, call_args))
    }

    x
  }

  out <- Map(resolve_one, LocalDynamicsArgs, names(LocalDynamicsArgs))
  names(out) <- names(LocalDynamicsArgs)
  out
}

###############################################################################
### Generic pathogen-state component for INApest abundance models
###############################################################################

INApestPathogen <- function(
    Model = c("SIS", "SIR", "SEIR", "Binary"),
    Beta = 0,
    RecoveryProb = 0,
    ProgressionProb = 1,
    PathogenMortalityProb = 0,
    ImmunityLossProb = 0,
    InitialInfected = 0,
    InitialExposed = 0,
    InitialRecovered = 0,
    IntroductionProb = 0,
    IntroductionNumber = 1,
    ContactMatrix = NULL,
    Transmission = c("frequency", "density"),
    DensityScale = 1,
    InitialPresent = 0,
    ClearanceProb = 0,
    TransmissionProb = NULL,
    PathogenHostExtinctionProb = 0,
    DetectionProb = 0,
    DetectionTriggersInfo = FALSE) {

  Model <- match.arg(Model)
  Transmission <- match.arg(Transmission)
  if (!is.logical(DetectionTriggersInfo) || length(DetectionTriggersInfo) != 1L || is.na(DetectionTriggersInfo))
    stop("DetectionTriggersInfo must be TRUE or FALSE")

  if (Model == "Binary") {
    if (is.null(TransmissionProb)) TransmissionProb <- Beta
    clip01 <- function(x) pmin(1, pmax(0, x))
    resolve_binary <- function(x, timestep, n_nodes, Ntimesteps, name) {
      if (is.function(x)) {
        fm <- names(formals(x)); a <- list(timestep=timestep,n_nodes=n_nodes,Ntimesteps=Ntimesteps)
        if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a),fm)]
        x <- do.call(x,a)
      }
      d <- dim(x)
      if (!is.null(d)) {
        if (length(d)==2L && all(d==c(n_nodes,Ntimesteps))) return(as.numeric(x[,timestep]))
        stop(name," must be scalar, length nodes, nodes x Ntimesteps, or resolver function")
      }
      x <- as.numeric(x)
      if (length(x)==1L) return(rep(x,n_nodes))
      if (length(x)==n_nodes) return(x)
      if (length(x)==Ntimesteps && Ntimesteps != n_nodes) return(rep(x[timestep],n_nodes))
      stop(name," has unsupported or ambiguous shape")
    }
    contact_binary <- function(timestep,n_nodes,Ntimesteps) {
      x <- ContactMatrix
      if (is.null(x)) return(diag(n_nodes))
      if (is.function(x)) {
        fm <- names(formals(x)); a <- list(timestep=timestep,n_nodes=n_nodes,Ntimesteps=Ntimesteps)
        if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a),fm)]
        x <- do.call(x,a)
      }
      d <- dim(x)
      if (length(d)==2L && all(d==c(n_nodes,n_nodes))) out <- x
      else if (length(d)==3L && all(d==c(n_nodes,n_nodes,Ntimesteps))) out <- x[,,timestep,drop=TRUE]
      else stop("Binary ContactMatrix must be nodes x nodes or nodes x nodes x Ntimesteps")
      if (any(!is.finite(out)) || any(out < 0)) stop("ContactMatrix entries must be finite and non-negative")
      as.matrix(out)
    }
    validate_binary <- function(n_nodes,Ntimesteps) {
      for (tt in seq_len(Ntimesteps)) {
        for (nm in c("RecoveryProb","PathogenMortalityProb","ImmunityLossProb")) {
          z <- resolve_binary(switch(nm, RecoveryProb=RecoveryProb, PathogenMortalityProb=PathogenMortalityProb, ImmunityLossProb=ImmunityLossProb),tt,n_nodes,Ntimesteps,nm)
          if (any(z != 0)) stop(nm," is not used by Model = 'Binary'; use ClearanceProb or PathogenHostExtinctionProb as appropriate")
        }
        for (nm in c("TransmissionProb","ClearanceProb","IntroductionProb","PathogenHostExtinctionProb","DetectionProb")) {
          x <- switch(nm, TransmissionProb=TransmissionProb, ClearanceProb=ClearanceProb, IntroductionProb=IntroductionProb, PathogenHostExtinctionProb=PathogenHostExtinctionProb, DetectionProb=DetectionProb)
          z <- resolve_binary(x,tt,n_nodes,Ntimesteps,nm)
          if (any(!is.finite(z)) || any(z < 0 | z > 1)) stop(nm," must resolve to [0,1]")
        }
        contact_binary(tt,n_nodes,Ntimesteps)
      }
      for (nm in c("InitialInfected","InitialExposed","InitialRecovered")) {
        z <- resolve_binary(switch(nm, InitialInfected=InitialInfected, InitialExposed=InitialExposed, InitialRecovered=InitialRecovered),1L,n_nodes,Ntimesteps,nm)
        if (any(z != 0)) stop(nm," is not used by Model = 'Binary'; use InitialPresent")
      }
      invisible(TRUE)
    }
    initial_binary <- function(Invaded,n_nodes=length(Invaded),Ntimesteps=1L) {
      validate_binary(n_nodes,Ntimesteps)
      z <- resolve_binary(InitialPresent,1L,n_nodes,Ntimesteps,"InitialPresent")
      if (any(!is.finite(z)) || any(z < 0 | z > 1)) stop("InitialPresent must resolve to [0,1]")
      if (all(z %in% c(0,1))) out <- as.integer(z) else out <- rbinom(n_nodes,1,z)
      as.integer(out * as.integer(Invaded > 0))
    }
    detect_binary <- function(PathogenPresent, timestep, Ntimesteps) {
      n_nodes <- length(PathogenPresent)
      p <- resolve_binary(DetectionProb, timestep, n_nodes, Ntimesteps, "DetectionProb")
      if (any(!is.finite(p)) || any(p < 0 | p > 1)) stop("DetectionProb must resolve to [0,1]")
      as.integer(stats::rbinom(n_nodes, 1L, clip01(p) * as.integer(PathogenPresent > 0)))
    }
    step_binary <- function(PathogenPresent,Invaded,timestep,Ntimesteps) {
      n_nodes <- length(Invaded); validate_binary(n_nodes,Ntimesteps)
      P <- as.integer(PathogenPresent > 0) * as.integer(Invaded > 0)
      tr <- resolve_binary(TransmissionProb,timestep,n_nodes,Ntimesteps,"TransmissionProb")
      clear <- resolve_binary(ClearanceProb,timestep,n_nodes,Ntimesteps,"ClearanceProb")
      intro <- resolve_binary(IntroductionProb,timestep,n_nodes,Ntimesteps,"IntroductionProb")
      ext <- resolve_binary(PathogenHostExtinctionProb,timestep,n_nodes,Ntimesteps,"PathogenHostExtinctionProb")
      C <- contact_binary(timestep,n_nodes,Ntimesteps)
      q <- sweep(C,1,tr*P,"*")
      q[q>1] <- 1
      ptrans <- 1-apply(1-q,2,prod)
      susceptible <- which(Invaded==1 & P==0)
      if (length(susceptible)) P[susceptible] <- rbinom(length(susceptible),1,clip01(ptrans[susceptible]))
      infected <- which(P==1)
      if (length(infected)) P[infected] <- P[infected] * rbinom(length(infected),1,1-clear[infected])
      candidates <- which(Invaded==1 & P==0)
      if (length(candidates)) P[candidates] <- pmax(P[candidates],rbinom(length(candidates),1,intro[candidates]))
      host_extinct <- integer(n_nodes); infected <- which(P==1)
      if (length(infected)) host_extinct[infected] <- rbinom(length(infected),1,ext[infected])
      Invaded[host_extinct==1] <- 0L; P[Invaded==0] <- 0L
      list(PathogenPresent=as.integer(P),Invaded=as.integer(Invaded),HostExtinction=host_extinct,NewPathogen=as.integer(P & !PathogenPresent))
    }
    out <- list(Model="Binary", States=c("Absent","Present"), TransmissionProb=TransmissionProb, ClearanceProb=ClearanceProb,
                InitialPresent=InitialPresent, IntroductionProb=IntroductionProb, ContactMatrix=ContactMatrix,
                PathogenHostExtinctionProb=PathogenHostExtinctionProb, DetectionProb=DetectionProb, DetectionTriggersInfo=DetectionTriggersInfo,
                binary_initial=initial_binary, binary_step=step_binary, binary_detect=detect_binary, binary_validate=validate_binary)
    class(out) <- c("INApestPathogen","list")
    return(out)
  }

  States <- switch(Model,
                   SIS = c("S", "I"),
                   SIR = c("S", "I", "R"),
                   SEIR = c("S", "E", "I", "R"))

  Spec <- list(
    Model = Model, Beta = Beta, RecoveryProb = RecoveryProb,
    ProgressionProb = ProgressionProb,
    PathogenMortalityProb = PathogenMortalityProb,
    ImmunityLossProb = ImmunityLossProb,
    InitialInfected = InitialInfected,
    InitialExposed = InitialExposed,
    InitialRecovered = InitialRecovered,
    IntroductionProb = IntroductionProb,
    IntroductionNumber = IntroductionNumber,
    ContactMatrix = ContactMatrix,
    Transmission = Transmission,
    DensityScale = DensityScale,
    DetectionProb = DetectionProb,
    DetectionTriggersInfo = DetectionTriggersInfo
  )

  clip01 <- function(x) pmin(1, pmax(0, x))

  context_dims <- function(context) {
    if (is.null(context$n_nodes) || is.null(context$Ntimesteps))
      stop("Pathogen context requires n_nodes and Ntimesteps")
    n_nodes <- as.integer(context$n_nodes)
    nts <- as.integer(context$Ntimesteps)
    n_lu <- if (is.null(context$n_landuses)) NULL else as.integer(context$n_landuses)
    list(n_nodes = n_nodes, Ntimesteps = nts, n_landuses = n_lu,
         n_units = if (is.null(n_lu)) n_nodes else n_nodes * n_lu)
  }

  resolve <- function(x, timestep, context, name) {
    dd <- context_dims(context)
    if (is.function(x)) {
      fm <- names(formals(x))
      args <- list(timestep = timestep, n_nodes = dd$n_nodes,
                   n_landuses = dd$n_landuses, Ntimesteps = dd$Ntimesteps)
      if (!is.null(fm) && !("..." %in% fm)) args <- args[intersect(names(args), fm)]
      x <- do.call(x, args)
    }

    d <- dim(x)
    if (is.null(dd$n_landuses)) {
      if (!is.null(d)) {
        if (length(d) == 2L && all(d == c(dd$n_nodes, dd$Ntimesteps)))
          return(as.numeric(x[, timestep]))
        stop(name, " matrix must have dimensions nodes x Ntimesteps")
      }
      x <- as.numeric(x)
      if (length(x) == 1L) return(rep(x, dd$n_units))
      if (length(x) == dd$n_nodes) return(x)
      if (length(x) == dd$Ntimesteps && dd$Ntimesteps != dd$n_nodes)
        return(rep(x[timestep], dd$n_units))
      stop(name, " must be scalar, length nodes, length Ntimesteps (when unambiguous), nodes x Ntimesteps, or a resolver function")
    }

    # Multiple-land-use context. Matrix/array forms are preferred when lengths
    # are ambiguous (e.g. number of nodes equals number of land uses).
    if (!is.null(d)) {
      if (length(d) == 2L && all(d == c(dd$n_nodes, dd$n_landuses)))
        return(as.numeric(x))
      if (length(d) == 3L && all(d == c(dd$n_nodes, dd$n_landuses, dd$Ntimesteps)))
        return(as.numeric(x[, , timestep]))
      stop(name, " for MLU must be nodes x land uses or nodes x land uses x Ntimesteps")
    }
    x <- as.numeric(x)
    if (length(x) == 1L) return(rep(x, dd$n_units))
    if (length(x) == dd$n_units) return(x)
    if (length(x) == dd$n_landuses && dd$n_landuses != dd$n_nodes)
      return(as.numeric(matrix(rep(x, each = dd$n_nodes), nrow = dd$n_nodes)))
    if (length(x) == dd$n_nodes && dd$n_nodes != dd$n_landuses)
      return(as.numeric(matrix(rep(x, dd$n_landuses), nrow = dd$n_nodes)))
    if (length(x) == dd$Ntimesteps && dd$Ntimesteps != dd$n_nodes && dd$Ntimesteps != dd$n_landuses)
      return(rep(x[timestep], dd$n_units))
    stop(name, " has an ambiguous or unsupported MLU shape; use a nodes x land uses matrix or nodes x land uses x Ntimesteps array")
  }

  contact_matrix <- function(timestep, context) {
    dd <- context_dims(context)
    x <- Spec$ContactMatrix
    if (is.null(x)) return(diag(dd$n_units))
    if (is.function(x)) {
      fm <- names(formals(x))
      args <- list(timestep = timestep, n_nodes = dd$n_nodes,
                   n_landuses = dd$n_landuses, Ntimesteps = dd$Ntimesteps)
      if (!is.null(fm) && !("..." %in% fm)) args <- args[intersect(names(args), fm)]
      x <- do.call(x, args)
    }
    d <- dim(x)
    if (is.null(d)) stop("ContactMatrix must be a square matrix, a 3-D source x target x timestep array, NULL, or a resolver function")
    if (length(d) == 2L) {
      if (!all(d == c(dd$n_units, dd$n_units)))
        stop("ContactMatrix must have dimensions pathogen units x pathogen units")
      out <- x
    } else if (length(d) == 3L) {
      if (!all(d == c(dd$n_units, dd$n_units, dd$Ntimesteps)))
        stop("Time-varying ContactMatrix must have dimensions pathogen units x pathogen units x Ntimesteps")
      out <- x[, , timestep, drop = TRUE]
    } else stop("ContactMatrix must be a matrix or 3-D array")
    if (any(!is.finite(out)) || any(out < 0))
      stop("ContactMatrix entries must be finite and non-negative")
    as.matrix(out)
  }

  validate_context <- function(context) {
    dd <- context_dims(context)
    prob_names <- c("RecoveryProb", "PathogenMortalityProb", "ImmunityLossProb", "IntroductionProb", "DetectionProb")
    if (Model == "SEIR") prob_names <- c(prob_names, "ProgressionProb")
    for (tt in seq_len(dd$Ntimesteps)) {
      contact_matrix(tt, context)
      for (nm in prob_names) {
        z <- resolve(Spec[[nm]], tt, context, nm)
        if (any(!is.finite(z)) || any(z < 0 | z > 1))
          stop(nm, " must resolve to values in [0,1] at every timestep")
      }
      rec <- resolve(Spec$RecoveryProb, tt, context, "RecoveryProb")
      mort <- resolve(Spec$PathogenMortalityProb, tt, context, "PathogenMortalityProb")
      if (any(rec + mort > 1 + 1e-12))
        stop("RecoveryProb + PathogenMortalityProb must not exceed 1 at any timestep")
      b <- resolve(Spec$Beta, tt, context, "Beta")
      if (any(!is.finite(b)) || any(b < 0))
        stop("Beta must resolve to finite non-negative values at every timestep")
      ds <- resolve(Spec$DensityScale, tt, context, "DensityScale")
      if (any(!is.finite(ds)) || any(ds <= 0))
        stop("DensityScale must resolve to positive values at every timestep")
      intro_n <- resolve(Spec$IntroductionNumber, tt, context, "IntroductionNumber")
      if (any(!is.finite(intro_n)) || any(intro_n < 0) || any(intro_n != floor(intro_n)))
        stop("IntroductionNumber must resolve to non-negative whole numbers at every timestep")
    }

    # Fail rather than silently ignore state-specific initial conditions or
    # processes that the selected compartment model cannot represent.
    if (!("E" %in% States)) {
      z <- resolve(Spec$InitialExposed, 1L, context, "InitialExposed")
      if (any(z != 0)) stop("InitialExposed must be zero unless Model = 'SEIR'")
    }
    if (!("R" %in% States)) {
      z <- resolve(Spec$InitialRecovered, 1L, context, "InitialRecovered")
      if (any(z != 0)) stop("InitialRecovered must be zero unless the model contains R")
      for (tt in seq_len(dd$Ntimesteps)) {
        z <- resolve(Spec$ImmunityLossProb, tt, context, "ImmunityLossProb")
        if (any(z != 0)) stop("ImmunityLossProb must be zero unless the model contains R")
      }
    }
    invisible(dd)
  }

  initial <- function(N, context) {
    validate_context(context)
    Nvec <- as.integer(as.numeric(N))
    resolve_count <- function(x, name) {
      z <- resolve(x, 1L, context, name)
      if (any(!is.finite(z)) || any(z < 0) || any(z != floor(z)))
        stop(name, " must resolve to non-negative whole-number host counts")
      as.integer(z)
    }
    I <- pmin(Nvec, resolve_count(Spec$InitialInfected, "InitialInfected"))
    E <- if ("E" %in% States) resolve_count(Spec$InitialExposed, "InitialExposed") else integer(length(Nvec))
    R <- if ("R" %in% States) resolve_count(Spec$InitialRecovered, "InitialRecovered") else integer(length(Nvec))
    if (any(I + E + R > Nvec)) stop("Initial pathogen-state counts exceed total host abundance")
    out <- matrix(0L, nrow = length(Nvec), ncol = length(States), dimnames = list(NULL, States))
    out[, "S"] <- Nvec - I - E - R
    out[, "I"] <- I
    if ("E" %in% States) out[, "E"] <- E
    if ("R" %in% States) out[, "R"] <- R
    out
  }

  thin_row <- function(counts, target) {
    counts <- as.integer(counts); target <- as.integer(target); total <- sum(counts)
    if (target >= total) return(counts)
    if (target <= 0L) return(integer(length(counts)))
    if (length(counts) == 1L) return(target)
    out <- integer(length(counts)); remain_draw <- target; remain_total <- total
    for (j in seq_len(length(counts) - 1L)) {
      out[j] <- rhyper(1L, counts[j], remain_total - counts[j], remain_draw)
      remain_draw <- remain_draw - out[j]
      remain_total <- remain_total - counts[j]
    }
    out[length(counts)] <- remain_draw
    out
  }

  reconcile <- function(State, N, context = NULL) {
    Nvec <- as.integer(as.numeric(N)); old <- rowSums(State)
    if (length(Nvec) != nrow(State)) stop("Host abundance and PathogenState lengths differ")
    for (i in seq_len(nrow(State))) {
      if (Nvec[i] < old[i]) State[i, ] <- thin_row(State[i, ], Nvec[i])
      else if (Nvec[i] > old[i]) State[i, "S"] <- State[i, "S"] + (Nvec[i] - old[i])
    }
    storage.mode(State) <- "integer"
    State
  }

  step <- function(State, N, timestep, context) {
    validate_context(context)
    n_units <- nrow(State); Nvec <- as.integer(as.numeric(N))
    beta <- pmax(0, resolve(Spec$Beta, timestep, context, "Beta"))
    rec <- resolve(Spec$RecoveryProb, timestep, context, "RecoveryProb")
    mort <- resolve(Spec$PathogenMortalityProb, timestep, context, "PathogenMortalityProb")
    prog <- resolve(Spec$ProgressionProb, timestep, context, "ProgressionProb")
    waning <- resolve(Spec$ImmunityLossProb, timestep, context, "ImmunityLossProb")
    intro_p <- resolve(Spec$IntroductionProb, timestep, context, "IntroductionProb")
    intro_n <- as.integer(resolve(Spec$IntroductionNumber, timestep, context, "IntroductionNumber"))
    density_scale <- resolve(Spec$DensityScale, timestep, context, "DensityScale")
    if (any(rec + mort > 1 + 1e-12)) stop("RecoveryProb + PathogenMortalityProb must not exceed 1")

    S0 <- as.integer(State[, "S"]); I0 <- as.integer(State[, "I"])
    E0 <- if ("E" %in% States) as.integer(State[, "E"]) else integer(n_units)
    R0 <- if ("R" %in% States) as.integer(State[, "R"]) else integer(n_units)
    live <- S0 + E0 + I0 + R0
    contact <- contact_matrix(timestep, context) # rows = infectious source, columns = recipient target
    infectious_pressure <- as.numeric(crossprod(I0, contact))
    if (Transmission == "frequency") {
      contact_population <- as.numeric(crossprod(live, contact))
      prevalence_pressure <- ifelse(contact_population > 0, infectious_pressure / contact_population, 0)
      foi <- beta * prevalence_pressure
    } else {
      foi <- beta * infectious_pressure / density_scale
    }
    p_inf <- clip01(-expm1(-pmax(0, foi)))
    new_inf <- rbinom(n_units, S0, p_inf)

    introduced <- integer(n_units)
    event <- rbinom(n_units, 1L, intro_p)
    idx <- which(event > 0 & S0 - new_inf > 0)
    if (length(idx)) introduced[idx] <- pmin(intro_n[idx], S0[idx] - new_inf[idx])

    progress <- if ("E" %in% States) rbinom(n_units, E0, prog) else integer(n_units)
    stay_i <- recover <- deaths <- integer(n_units)
    for (i in which(I0 > 0)) {
      z <- as.numeric(rmultinom(1L, I0[i], c(1 - rec[i] - mort[i], rec[i], mort[i])))
      stay_i[i] <- z[1]; recover[i] <- z[2]; deaths[i] <- z[3]
    }
    lose <- if ("R" %in% States) rbinom(n_units, R0, waning) else integer(n_units)

    State[, "S"] <- S0 - new_inf - introduced + lose
    if ("E" %in% States) {
      State[, "E"] <- E0 - progress + new_inf + introduced
      State[, "I"] <- stay_i + progress
    } else State[, "I"] <- stay_i + new_inf + introduced
    if ("R" %in% States) State[, "R"] <- R0 - lose + recover
    else State[, "S"] <- State[, "S"] + recover

    Nvec <- Nvec - deaths
    if (any(Nvec < 0) || any(rowSums(State) != Nvec)) stop("Internal pathogen-state conservation error")
    list(State = State, N = Nvec, Deaths = deaths,
         Introduced = introduced, NewInfections = new_inf)
  }

  Engine <- list(States = States, Validate = validate_context, Resolve = resolve,
                 ContactMatrix = contact_matrix, Initial = initial,
                 Reconcile = reconcile, Step = step)
  structure(c(Spec, list(States = States, Engine = Engine)), class = "INApestPathogen")
}

INApestPathogenContactMatrix <- function(NodeContact, LandUseMixing = NULL) {
  NodeContact <- as.matrix(NodeContact)
  if (nrow(NodeContact) != ncol(NodeContact) || any(!is.finite(NodeContact)) || any(NodeContact < 0))
    stop("NodeContact must be a finite non-negative square source x target matrix")
  if (is.null(LandUseMixing)) return(NodeContact)
  LandUseMixing <- as.matrix(LandUseMixing)
  if (nrow(LandUseMixing) != ncol(LandUseMixing) || any(!is.finite(LandUseMixing)) || any(LandUseMixing < 0))
    stop("LandUseMixing must be a finite non-negative square source x target matrix")
  # MLU state is flattened from nodes x land uses with node varying fastest.
  # Kronecker order therefore gives source/target blocks by land-use pair.
  kronecker(LandUseMixing, NodeContact)
}

INApestPathogenOutputs <- function(PathogenStateResults) {
  d <- dim(PathogenStateResults)
  if (length(d) == 4L) {
    st <- dimnames(PathogenStateResults)[[2]]
    get_state <- function(nm) if (nm %in% st) PathogenStateResults[, nm, , , drop = FALSE] else NULL
    I <- get_state("I")
    return(list(S = get_state("S"), E = get_state("E"), I = I, R = get_state("R"),
                PathogenPresence = I > 0,
                InfectiousNodeCount = apply(I > 0, c(3, 4), sum),
                TotalInfectious = apply(I, c(3, 4), sum)))
  }
  if (length(d) == 5L) {
    st <- dimnames(PathogenStateResults)[[3]]
    get_state <- function(nm) if (nm %in% st) PathogenStateResults[, , nm, , , drop = FALSE] else NULL
    I <- get_state("I")
    return(list(S = get_state("S"), E = get_state("E"), I = I, R = get_state("R"),
                PathogenPresence = I > 0,
                InfectiousNodeCount = apply(I > 0, c(4, 5), sum),
                TotalInfectious = apply(I, c(4, 5), sum)))
  }
  stop("Expected Meta or MLU PathogenStateResults array")
}


INApestMetaParallelMultipleLandUse = function(
ModelName, #Name for storing results to file 
Nperm,                  #Number of permutations per parameter combination
Ntimesteps,                 #Simulation duration timesteps can be any length of time
LocalDynamics = local.dynamicsLU, #Local population growth, dispersal and management function; user-defined functions are supported
LocalDynamicsArgs = list(), #Named custom arguments passed to LocalDynamics; wrap time-varying values with INApestLocalDynamicsTimeArg()
Pathogen = NULL, #Optional INApestPathogen() specification; N remains host abundance by node x land use
Nlanduses,
DetectionProb,          #Vector of probabilities per land use, matrix (nodes x land uses) or 3d array (land uses x nodes x timesteps). Must be between 0 and 1.
DetectionSD = NULL, #Option to provide standard deviation for management probability can be single number or vector (nodes)
ManageProb,             ##Vector of probabilities per land use, matrix (nodes x land uses) or 3d array (land uses x nodes x timesteps). Must be between 0 and 1.
ManageSD = NULL, #Option to provide standard deviation for management probability can be single number or vector (nodes)
MortalityProb,           ##Vector of probabilities per land use, matrix (nodes x land uses) or 3d array (land uses x nodes x timesteps). Must be between 0 and 1.
MortalitySD = NULL, #Option to provide standard deviation for management probability can be single number or vector (nodes)
FecundityReduction = 0, #Proportional reduction in per-capita fecundity under management: scalar, vector (land uses), matrix (nodes x land uses), or array (nodes x land uses x timesteps)
SpreadReduction,        ##Vector of probabilities per land use, matrix (nodes x land uses) or 3d array (land uses x nodes x timesteps). Must be between 0 and 1.
SpreadReductionSD = NULL, #Option to provide standard deviation for management probability can be single number or vector (nodes)
InitialPopulation = NA,        #Matrix (nodes x land uses) of population sizes at start of simulations
InitBioP = NA,		#Proportion of nodes infested at start of simulations
InvasionRisk = NA,           #Vector or matrix (nodes x timesteps) of probabilities of invasion from external sources
InitialInfo = NA,        #Vector or of nodes with information at start of simulations
InitInfoP = 0,		#Proportion of nodes with information at start of simulations
ExternalInfoProb = NA,           #Vector of probabilities of communication from external sources
InfoRetentionProb = 1,       #Probability that existing information is retained between timesteps. Can be single number, vector (nodes) or matrix (nodes x timesteps)
InfoPersistenceSteps = NA,    #Number of timesteps information persists after last known local presence. Can be single number, vector (nodes) or matrix (nodes x timesteps); NA uses InfoRetentionProb
EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
Survival = 1,           # local population survival probability. Set to 1 for no environmental limitation on survival. Can be single number, vector (nodes) or matrix (nodes x timesteps)
K,		       #Population carrying capacity - vector (nodes)
PropaguleProduction, #Propagules produced per individual, can be single value, vector (nodes) or matrix (nodes x years)
PropaguleEstablishment, #Propagules establishment probability. The likelihood of a dispersing propagule encountering a single
                        #host plant or establishment site within a node. Can be a ratio of search radius or patch size to node area
IncursionStartPop=NA,      #option to set population size for new incursions
SDDprob,                   #Natural dispersal probability matrix, or 3D array (nodes x nodes x timesteps)
SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
LDDprob = NA,         #Option to provide long distance (human-mediated) dispersal matrix or 3D array (nodes x nodes x timesteps) instead of distance-independent dispesal rate
			      #e.g. could be weighted by law of human visitation or data on stock movements
LDDrate = 0,         #Proportion of available propagules entering LDD
OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
OngoingExternalInfo = F,   ##Option to include ongoing communication from external sources
OutputDir = NA,		      #Directory for storing results
DoPlots = TRUE	     #Option to omit printing of line graphs.Default is to print.
)
{
if(!is.function(LocalDynamics))
  stop("LocalDynamics must be a function")
if (is.null(LocalDynamicsArgs)) LocalDynamicsArgs <- list()
if (!is.list(LocalDynamicsArgs))
  stop("LocalDynamicsArgs must be a named list")
if (length(LocalDynamicsArgs) &&
    (is.null(names(LocalDynamicsArgs)) || any(!nzchar(names(LocalDynamicsArgs)))))
  stop("Every LocalDynamicsArgs entry must have a non-empty name")
if (anyDuplicated(names(LocalDynamicsArgs)))
  stop("LocalDynamicsArgs names must be unique")
# Force the argument before any parallel worker closure is created. This keeps
# the selected default or user-supplied function as an explicit model input.
force(LocalDynamics)
force(LocalDynamicsArgs)
UserLocalDynamicsArgs <- LocalDynamicsArgs
PathogenOriginal <- Pathogen
if(!is.null(PathogenOriginal) && !inherits(PathogenOriginal, "INApestPathogen")) stop("Pathogen must be NULL or an object returned by INApestPathogen()")
UsePathogen <- !is.null(PathogenOriginal)
if(UsePathogen) { PathogenEngine <- PathogenOriginal$Engine; PathogenContext <- list(n_nodes = nrow(SDDprob), n_landuses = Nlanduses, Ntimesteps = Ntimesteps); PathogenEngine$Validate(PathogenContext); force(PathogenEngine); force(PathogenContext) }
###POTENTIAL ADDITIONS
###1) Make detection prob a function of population size. Could be based on individual detection prob so that DetectionProb = 1-(1-DPindividual)^N)
###   DPindividual could vary between nodes


###Allow SDD and LDD connectivity to vary through time
if(length(dim(SDDprob)) == 3 && (dim(SDDprob)[1] != dim(SDDprob)[2] || dim(SDDprob)[3] != Ntimesteps))
  stop("SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
if(length(dim(LDDprob)) == 3 && (dim(LDDprob)[1] != nrow(SDDprob) || dim(LDDprob)[2] != nrow(SDDprob) || dim(LDDprob)[3] != Ntimesteps))
  stop("LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")

###Allow management-induced fecundity reduction to vary by land use, node and time
FecundityReductionDims <- dim(FecundityReduction)
if(is.null(FecundityReductionDims)) {
  if(!(length(FecundityReduction) == 1 || length(FecundityReduction) == Nlanduses))
    stop("FecundityReduction must be scalar, length Nlanduses, nodes x land uses, or nodes x land uses x Ntimesteps")
} else if(length(FecundityReductionDims) == 2) {
  if(!all(FecundityReductionDims == c(nrow(SDDprob), Nlanduses)))
    stop("FecundityReduction matrix must have dimensions nodes x land uses")
} else if(length(FecundityReductionDims) == 3) {
  if(!all(FecundityReductionDims == c(nrow(SDDprob), Nlanduses, Ntimesteps)))
    stop("FecundityReduction array must have dimensions nodes x land uses x Ntimesteps")
} else {
  stop("FecundityReduction must be scalar, length Nlanduses, nodes x land uses, or nodes x land uses x Ntimesteps")
}
if(any(!is.finite(FecundityReduction)) || any(FecundityReduction < 0) || any(FecundityReduction > 1))
  stop("FecundityReduction values must be between 0 and 1")

ResolveFecundityReductionLU <- function(timestep) {
  if(is.null(dim(FecundityReduction))) {
    if(length(FecundityReduction) == 1)
      return(matrix(FecundityReduction, nrow = nrow(SDDprob), ncol = Nlanduses))
    return(matrix(rep(FecundityReduction, each = nrow(SDDprob)), nrow = nrow(SDDprob), ncol = Nlanduses))
  }
  if(length(dim(FecundityReduction)) == 2) return(FecundityReduction)
  FecundityReduction[,,timestep]
}

###Allow information retention to vary by node and through time
if(is.matrix(InfoRetentionProb) == T && (nrow(InfoRetentionProb) != nrow(SDDprob) || ncol(InfoRetentionProb) != Ntimesteps))
  stop("InfoRetentionProb matrix must have dimensions nodes x Ntimesteps")
if(is.matrix(InfoRetentionProb) == F && !(length(InfoRetentionProb) == 1 || length(InfoRetentionProb) == nrow(SDDprob)))
  stop("InfoRetentionProb must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
if(any(is.na(InfoRetentionProb)) || any(InfoRetentionProb < 0) || any(InfoRetentionProb > 1))
  stop("InfoRetentionProb values must be between 0 and 1")

###Allow programmed information persistence after last known local presence to vary by node and through time
if(is.matrix(InfoPersistenceSteps) == T && (nrow(InfoPersistenceSteps) != nrow(SDDprob) || ncol(InfoPersistenceSteps) != Ntimesteps))
  stop("InfoPersistenceSteps matrix must have dimensions nodes x Ntimesteps")
if(is.matrix(InfoPersistenceSteps) == F && !(length(InfoPersistenceSteps) == 1 || length(InfoPersistenceSteps) == nrow(SDDprob)))
  stop("InfoPersistenceSteps must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
if(any(!is.na(InfoPersistenceSteps) & (!is.finite(InfoPersistenceSteps) | InfoPersistenceSteps < 0 | InfoPersistenceSteps != floor(InfoPersistenceSteps))))
  stop("InfoPersistenceSteps values must be non-negative whole numbers or NA")
UseInfoPersistence = any(!is.na(InfoPersistenceSteps))
if(UseInfoPersistence == T && any(InfoRetentionProb < 1))
  warning("Both InfoPersistenceSteps and InfoRetentionProb specify information loss. Programmed stopping takes priority where InfoPersistenceSteps is not NA; InfoRetentionProb is only used where InfoPersistenceSteps is NA.",call. = F)

# pre-evaluate some variables for efficiency
if(length(dim(K)) <3)
{
K_is_0 <- rowSums(K)<=0
inv_K <- 1 / sum(colSums(K))
NodeK = K
Pk  = K/rowSums(K)
Pk[is.na(Pk)] = 0
}

if(length(dim(K)) ==3)
   {
    K_is_0 <- rowSums(K[,,1])<=0
    NodeK = K[,,1]
    Pk  = NodeK/rowSums(NodeK)
    Pk[is.na(Pk)] = 0
   }  
if(length(dim(PropaguleProduction)) < 3)
  NodePropaguleProduction = PropaguleProduction

if(length(dim(PropaguleEstablishment)) <3)
  NodePropaguleEstablishment = PropaguleEstablishment

if(length(dim(EnvEstabProb)) <3)
  NodeEnvEstabProb <- EnvEstabProb

if(length(dim(Survival)) < 3)
  NodeSurvival <- Survival


###Declare matrix for information spread simulations
if(is.matrix(SEAM) == T)
     {
     diag(SEAM) = 0
     RandSEAM <- matrix(NA,nrow = nrow(SDDprob),ncol=nrow(SDDprob))
     }

###Assign standard deviation value to management in no value provided
if(is.null(ManageSD) == T)
	ManageSD = ManageProb/10
if(is.null(SpreadReductionSD) == T)
	SpreadReductionSD = (1-SpreadReduction)/10
if(is.null(DetectionSD) == T)
	DetectionSD = DetectionProb/10
if(is.null(MortalitySD) == T)
	MortalitySD = MortalityProb/10


###########################################################
###Start of simulation
###########################################################
    
detected_cores <- parallel::detectCores()
if (is.na(detected_cores)) detected_cores <- 2L
n_cores <- max(1L, min(Nperm, detected_cores - 1L))

###Run one stochastic realisation. Function arguments and local helpers are
###captured in this closure, avoiding fragile manual worker export lists.
PermutationWorker <- function(i_perm)
  {
  ###Set initial dispersal connectivity
  NodeSDDprob = SDDprob
  if(length(dim(SDDprob)) == 3)
    NodeSDDprob = SDDprob[,,1]
  NodeLDDprob = LDDprob
  if(length(dim(LDDprob)) == 3)
    NodeLDDprob = LDDprob[,,1]
    ###Custom function for allocating exact number of cases to land uses
    ###Currently very clunky. Would be good to seek built in R approach rather than
    ###custom function
    SampleVector <- function(X)
      {
      Vect = vector(length=0)
      for(i  in 2:length(X))
        Vect <- c(Vect,rep((i-1),times = X[i]))
      Sample <- sample(Vect,size = X[1],replace = F)
      Out <- vector(length = length(X)-1)
      for(j in 1:length(Out))
        Out[j] <- length(Sample[Sample==j])
      return(Out)  
      }
    
    
  InvasionResultsLoop <- array(dim = c(nrow(SDDprob),Nlanduses,Ntimesteps))
  PopulationResultsLoop <- InvasionResultsLoop
  ManagingResultsLoop <- array(dim = c(nrow(SDDprob),Nlanduses,Ntimesteps))
  DetectedResultsLoop <- array(dim = c(nrow(SDDprob),Nlanduses,Ntimesteps))
  if(UsePathogen) {
    PathogenStates <- PathogenEngine$States
    PathogenStateResultsLoop <- array(0L, dim = c(nrow(SDDprob), Nlanduses, length(PathogenStates), Ntimesteps), dimnames = list(NULL, NULL, PathogenStates, NULL))
    PathogenDetectedResultsLoop <- matrix(0L,nrow=nrow(SDDprob),ncol=Ntimesteps)
  }
###Assign initial infestations according either to "InitialInvasion" binary vector OR
###"InvasionRisk" probabilities and/or initial proportion of nodes infested ("InitBioP") OR
###just "InitBioP" if neither "InitialInvasion" or "InvasionRisk" supplied by user
InitBio = matrix(ncol = Nlanduses, nrow = nrow(SDDprob))
InitBio[,] = 0
InintInfested = rep(0,times = nrow(SDDprob))

if(is.matrix(InitialPopulation) == T && nrow(InitialPopulation) == nrow(SDDprob) && ncol(InitialPopulation) == Nlanduses)
  InitBio = InitialPopulation

if(is.matrix(InitialPopulation) == F || nrow(InitialPopulation) != nrow(SDDprob) || ncol(InitialPopulation) != Nlanduses)
{
risk = NULL
if(is.matrix(InvasionRisk) == T && nrow(InvasionRisk) == nrow(SDDprob))
  risk = InvasionRisk[,1]
if(is.matrix(InvasionRisk) == F && length(InvasionRisk) == nrow(SDDprob))
  risk = InvasionRisk

if(is.na(InitBioP) == F)
  Infested = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitBioP),prob = risk)
if(is.na(InitBioP) == T && is.null(risk) == F)
  {
  Infested = rbinom(1:nrow(SDDprob),size = 1,prob = risk)
  Infested = which(Infested == 1)
  }
if(is.na(InitBioP) == T && is.null(risk) == T)
  Infested = integer(0)

if(is.na(IncursionStartPop) == T)
  InintInfested[Infested] = 1
if(is.na(IncursionStartPop) == F)
  InintInfested[Infested] = IncursionStartPop
###Find alternative to for loop
InVector = cbind(InintInfested,K)
InitialPopulation <- t(apply(InVector,1,FUN = SampleVector))
InitBio = InitialPopulation
}

###Ensure initial population not greater than carrying capacity
for(i in 1:Nlanduses)
  InitBio[,i] = apply(cbind(NodeK[,i], InitBio[,i]),MARGIN = 1,FUN = min)

# initialise the population
N <- InitBio
if(UsePathogen) { PathogenStateFlat <- PathogenEngine$Initial(N, PathogenContext) }
if(sum(N) == 0 && OngoingExternalInvasion == F)
  warning("No initial populations and no future external invasions")

###Select nodes with information at start of simulation  according either to "InitialInfo" binary vector OR
###"ExternalInfoProb" probabilities and/or initial proportion of nodes with information ("InitInfoP") OR
###just "InitInfoP" if neither "InitialInfo" or "ExternalInfoProb" supplied by user.
###If no initial info variables provided, no nodes have info at start of simulations
InitInfo = rep(0,times = nrow(SDDprob))
if(length(InitialInfo) == nrow(SDDprob) || (is.na(InitInfoP) == F && InitInfoP>0) || is.na(sum(ExternalInfoProb)) == F )
{
if(length(InitialInfo) != nrow(SDDprob))
  {
  if(length(ExternalInfoProb) == nrow(SDDprob))
    {
    if(is.na(InitInfoP) == F)
      Info = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitInfoP),prob = ExternalInfoProb)
    if(is.na(InitInfoP) == T)
      {
      Info = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb) 
      Info = which(Info == 1)
      } 
    }
  if(length(ExternalInfoProb) != nrow(SDDprob))
    {
    if(is.matrix(ExternalInfoProb) == F)
      Info = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitInfoP))
    if(is.matrix(ExternalInfoProb) == T)
      {
      Info = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb[,1])
      Info = which(Info == 1)
      }
    }
  InitInfo[Info] = 1
  
  }
if(length(InitialInfo) == nrow(SDDprob))
  InitInfo = InitialInfo  
}

###Randomly assign annual detection probability, based on mean and sd
###If DetectionProb given as vector (land uses) or matrix (nodes x land uses) 
if((is.matrix(DetectionProb)==FALSE &&(length(DetectionProb) == Nlanduses) ||(is.matrix(DetectionProb)==TRUE && nrow(DetectionProb) == nrow(SDDprob) )))
      {
      NodeDetectionProb = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
      for(lu in 1:Nlanduses)
        NodeDetectionProb[,lu] = rnorm(DetectionProb[lu],DetectionSD[lu],n = nrow(SDDprob))
      NodeDetectionProb[NodeDetectionProb<0] = 0
      NodeDetectionProb[NodeDetectionProb>1] = 1
      }

###If DetectionProb given as 3d array (nodes x land uses x timesteps) use values for first timestep to get initial detections
if(length(dim(DetectionProb))==3)
      {
      NodeDetectionProb = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
      for(lu in 1:Nlanduses)
        NodeDetectionProb[,lu] = rnorm(DetectionProb[,lu,1],DetectionSD[,lu,1],n = nrow(SDDprob))
      NodeDetectionProb[NodeDetectionProb<0] = 0
      NodeDetectionProb[NodeDetectionProb>1] = 1
      }


###Randomly assign probability of mangement adoption upon detection of infestation
###If ManageProb given as vector (land uses) or matrix (nodes x land uses) 
if((is.matrix(ManageProb)==FALSE &&(length(ManageProb) == Nlanduses) ||(is.matrix(ManageProb)==TRUE && nrow(ManageProb) == nrow(SDDprob) )))
  {
  NodeManageProb = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
  for(lu in 1:Nlanduses)
    NodeManageProb[,lu] = rnorm(ManageProb[lu],ManageSD[lu],n = nrow(SDDprob))
  NodeManageProb[NodeManageProb<0] = 0
  NodeManageProb[NodeManageProb>1] = 1
  }

###Randomly assign spread reduction factor when management adopted
###If SpreadReduction given as vector (land uses) or matrix (nodes x land uses) 
if((is.matrix(SpreadReduction)==FALSE &&(length(SpreadReduction) == Nlanduses) ||(is.matrix(SpreadReduction)==TRUE && nrow(SpreadReduction) == nrow(SDDprob) )))
  {
  NodeSpreadReduction = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
  for(lu in 1:Nlanduses)
    NodeSpreadReduction[,lu] = rnorm(SpreadReduction[lu],ManageSD[lu],n = nrow(SDDprob))
  NodeSpreadReduction[NodeSpreadReduction<0] = 0
  NodeSpreadReduction[NodeSpreadReduction>1] = 1
  }


###Randomly assign mortality probability when management applied
###If MortalityProb given as vector (land uses) or matrix (nodes x land uses) 
if((is.matrix(MortalityProb)==FALSE &&(length(MortalityProb) == Nlanduses) ||(is.matrix(MortalityProb)==TRUE && nrow(MortalityProb) == nrow(SDDprob) )))
  {
  NodeMortalityProb = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
  for(lu in 1:Nlanduses)
    NodeMortalityProb[,lu] = rnorm(MortalityProb[lu],MortalitySD[lu],n = nrow(SDDprob))
  NodeMortalityProb[NodeMortalityProb<0] = 0
  NodeMortalityProb[NodeMortalityProb>1] = 1
  }

###Populate invasion status vector ahead of timestep loop
Invaded = ifelse(InitBio>0,1,0) 


###Probability of info at start of simulation depends on
###Presence of pest and detection probability
###Select nodes that have detected infestation 

LUdetectionProb = 1-(1-NodeDetectionProb)^(InitBio)
InitDetection = rbinom(1:nrow(SDDprob),size = 1,prob = 1-apply(1-LUdetectionProb,1,prod))
InitInfo[InitInfo == 0] = InitDetection[InitInfo == 0]
###Populate information status vector ahead of timestep loop
HaveInfo = InitInfo
if(UsePathogen && isTRUE(PathogenOriginal$DetectionTriggersInfo))
  {
  PathogenDetectionP <- PathogenEngine$Resolve(PathogenOriginal$DetectionProb, 1L, PathogenContext, "DetectionProb")
  p_not <- (1-PathogenDetectionP)^PathogenStateFlat[,"I"]
  p_node <- 1-apply(matrix(p_not,nrow=nrow(SDDprob),ncol=Nlanduses),1,prod)
  InitPathogenDetection <- rbinom(nrow(SDDprob),1,p_node)
  HaveInfo[HaveInfo == 0] <- InitPathogenDetection[HaveInfo == 0]
  }

###Track the most recent timestep with known local presence
LastKnownPresence = rep(NA,nrow(SDDprob))
if(UseInfoPersistence == T)
  {
  InitialKnownPresence = which(InitDetection == 1)
  if(length(InitialKnownPresence) > 0)
    LastKnownPresence[InitialKnownPresence] = 0
  }

 
###Declare matrices outside loop
Managing = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
N0 = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
Propagules = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
Recruits = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))  
  # run simulation
for (timestep in 1:Ntimesteps) 
  { 
 
  ###Allow for variation in dispersal connectivity through time
  if(length(dim(SDDprob)) == 3)
    NodeSDDprob = SDDprob[,,timestep]
  if(length(dim(LDDprob)) == 3)
    NodeLDDprob = LDDprob[,,timestep]
  NodeFecundityReduction <- ResolveFecundityReductionLU(timestep)

  ###Allow for variation in establishment through time
  ###e.g.  climate change predictions
  ###Note: could be done outside loop, but would take heaps of memory to store 
  if(is.matrix(EnvEstabProb) == T)
    NodeEnvEstabProb <- EnvEstabProb[,timestep]
   
  if(is.matrix(Survival) == T)
    NodeSurvival <- Survival[,timestep]
    
    
  ###If carrying capacity provided as 3D array assign values for relevant timestep
  if(length(dim(K)) == 3)
    {
    K_is_0 <- K[,,timestep]<=0
    inv_K <- 1 / sum(K[,,timestep])
    NodeK = K[,,timestep] 
    }  

  ###If propagule production provided as matrix assign values for relevant timestep
  if(is.matrix(PropaguleProduction) == TRUE)
    NodePropaguleProduction = PropaguleProduction[,timestep] 
  
  if(is.matrix(PropaguleEstablishment) == TRUE)
    NodePropaguleEstablishment = PropaguleEstablishment[,timestep]
      
  ###Randomly assign annual detection probability, based on mean and sd
  ###If DetectionProb given as 3d array (nodes x land uses x timesteps)
  if(length(dim(DetectionProb))==3)
    {
    NodeDetectionProb = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
    for(lu in 1:Nlanduses)
      NodeDetectionProb[,lu] = rnorm(DetectionProb[,lu,timestep],DetectionSD[,lu,timestep],n = nrow(SDDprob))
    NodeDetectionProb[NodeDetectionProb<0] = 0
    NodeDetectionProb[NodeDetectionProb>1] = 1
    }

  ###Randomly assign probability of mangement adoption upon detection of infestation
  ###If ManageProb given as 3d array (nodes x land uses x timesteps)
  if(length(dim(ManageProb))==3)
    {
    NodeManageProb = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
    for(lu in 1:Nlanduses)
      NodeManageProb[,lu] = rnorm(ManageProb[,lu,timestep],ManageSD[,lu,timestep],n = nrow(SDDprob))
    NodeManageProb[NodeManageProb<0] = 0
    NodeManageProb[NodeManageProb>1] = 1
    }

  ###Randomly assign spread reduction factor when management adopted
  ###If DetectionProb given as 3d array (nodes x land uses x timesteps)
  if(length(dim(SpreadReduction))==3)
    {
    NodeSpreadReduction = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
    for(lu in 1:Nlanduses)
      NodeSpreadReduction[,lu] = rnorm(SpreadReduction[,lu,timestep],DetectionSD[,lu,timestep],n = nrow(SDDprob))
    NodeSpreadReduction[NodeSpreadReduction<0] = 0
    NodeSpreadReduction[NodeSpreadReduction>1] = 1
    }
  
  ###Randomly assign annual eradication probability when management applied
  ###If MortalityProb given as 3d array (nodes x land uses x timesteps)
  if(length(dim(MortalityProb))==3)
    {
    NodeMortalityProb = matrix(ncol = Nlanduses,nrow = nrow(SDDprob))
    for(lu in 1:Nlanduses)
      NodeMortalityProb[,lu] = rnorm(MortalityProb[,lu,timestep],MortalitySD[,lu,timestep],n = nrow(SDDprob))
    NodeMortalityProb[NodeMortalityProb<0] = 0
    NodeMortalityProb[NodeMortalityProb>1] = 1
    }


  ###Assign management status to nodes   
  ###Management is only applied to nodes which have information
  ###i.e. where pest has been detected or following communication of information
  ###from neighbouring infested farms 
  Managing[] = rbinom(Nlanduses*nrow(SDDprob),size = 1,prob = NodeManageProb*HaveInfo)
  ###Identify nodes with known extant infestations 
  Detected = Invaded*HaveInfo
  
  
  
  ###Adjust starting population for natural and managed mortality
  N0[] = rbinom(Nlanduses*nrow(SDDprob),size = c(N),prob = NodeSurvival*(1-NodeMortalityProb*Managing))
    ###Track known local presence from actual management mortality
    ###Condition on realised total deaths so the existing population draw is unchanged
    if(UseInfoPersistence == T)
      {
      ManagementMortality = NodeMortalityProb*Managing
      TotalMortalityProb = 1-NodeSurvival*(1-ManagementMortality)
      ConditionalManagementMortality = matrix(0,nrow = nrow(SDDprob),ncol = Nlanduses)
      ManagementMortalityCells = which((N-N0) > 0 & ManagementMortality > 0 & TotalMortalityProb > 0)
      if(length(ManagementMortalityCells) > 0)
        ConditionalManagementMortality[ManagementMortalityCells] = (NodeSurvival*ManagementMortality)[ManagementMortalityCells]/TotalMortalityProb[ManagementMortalityCells]
      ManagementKillProb = 1-apply((1-ConditionalManagementMortality)^(N-N0),1,prod)
      ManagementKillProb[ManagementKillProb < 0] = 0
      ManagementKillProb[ManagementKillProb > 1] = 1
      CertainManagementKillNodes = which(ManagementKillProb >= 1)
      if(length(CertainManagementKillNodes) > 0)
        LastKnownPresence[CertainManagementKillNodes] = timestep
      PotentialManagementKillNodes = which(ManagementKillProb > 0 & ManagementKillProb < 1)
      if(length(PotentialManagementKillNodes) > 0)
        {
        ManagementKilled = rbinom(n = length(PotentialManagementKillNodes),size = 1,prob = ManagementKillProb[PotentialManagementKillNodes])
        KnownPresence = PotentialManagementKillNodes[ManagementKilled == 1]
        if(length(KnownPresence) > 0)
          LastKnownPresence[KnownPresence] = timestep
        }
      }
  if(sum(N0)<=0 )
    N = N0
  Pin <-0
  Qin <- 0  
    # natural dispersal 
  if(sum(N0)>0 ) 
  {
  CoreLocalDynamicsArgs <- list(
    sddprob = NodeSDDprob,
    nodepropaguleproduction = NodePropaguleProduction,
    nodeenvestabprob = NodeEnvEstabProb,
    n = N0,
    lddprob = NodeLDDprob,
    lddrate = LDDrate,
    k_is_0 = K_is_0,
    nodeK = NodeK,
    nodepropaguleestablishment = NodePropaguleEstablishment,
    nodespreadreduction = NodeSpreadReduction,
    managing = Managing
  )
  LocalDynamicsFormals <- names(formals(LocalDynamics))
  LocalDynamicsAcceptsFecundityReduction <-
    "nodefecundityreduction" %in% LocalDynamicsFormals || "..." %in% LocalDynamicsFormals
  if (LocalDynamicsAcceptsFecundityReduction)
    CoreLocalDynamicsArgs$nodefecundityreduction <- NodeFecundityReduction
  else if (any(NodeFecundityReduction * Managing > 0))
    stop("Custom LocalDynamics must accept a 'nodefecundityreduction' argument (or ...) when FecundityReduction is active")
  ResolvedLocalDynamicsArgs <- .resolve_INApest_LocalDynamicsArgs(
      UserLocalDynamicsArgs, timestep = timestep, Ntimesteps = Ntimesteps
    )
    if (length(ResolvedLocalDynamicsArgs)) {
      duplicate_args <- intersect(names(ResolvedLocalDynamicsArgs), names(CoreLocalDynamicsArgs))
      if (length(duplicate_args))
        stop("LocalDynamicsArgs may not override INApest core LocalDynamics argument(s): ",
             paste(duplicate_args, collapse = ", "))
      LocalDynamicsFormals <- names(formals(LocalDynamics))
      unknown_args <- setdiff(names(ResolvedLocalDynamicsArgs), LocalDynamicsFormals)
      if (!("..." %in% LocalDynamicsFormals) && length(unknown_args))
        stop("Custom LocalDynamics does not accept LocalDynamicsArgs argument(s): ",
             paste(unknown_args, collapse = ", "))
      CoreLocalDynamicsArgs <- c(CoreLocalDynamicsArgs, ResolvedLocalDynamicsArgs)
    }
    N <- do.call(LocalDynamics, CoreLocalDynamicsArgs)
  } 
 ###Apply programmed stopping after last known local presence
NodeInfoPersistenceSteps = InfoPersistenceSteps
if(is.matrix(InfoPersistenceSteps) == T)
  NodeInfoPersistenceSteps = InfoPersistenceSteps[,timestep]
if(length(NodeInfoPersistenceSteps) == 1)
  NodeInfoPersistenceSteps = rep(NodeInfoPersistenceSteps,nrow(SDDprob))
ProgrammedInfoNodes = which(HaveInfo == 1 & !is.na(NodeInfoPersistenceSteps))
if(length(ProgrammedInfoNodes) > 0)
  {
  TimeSinceKnownPresence = timestep-LastKnownPresence
  InfoStopNodes = ProgrammedInfoNodes[is.na(LastKnownPresence[ProgrammedInfoNodes]) | TimeSinceKnownPresence[ProgrammedInfoNodes] >= NodeInfoPersistenceSteps[ProgrammedInfoNodes]]
  if(length(InfoStopNodes) > 0)
    HaveInfo[InfoStopNodes] = 0
  }

###Allow information to decay after management and spread where no programmed stop is supplied
NodeInfoRetentionProb = InfoRetentionProb
if(is.matrix(InfoRetentionProb) == T)
  NodeInfoRetentionProb = InfoRetentionProb[,timestep]
if(length(NodeInfoRetentionProb) == 1)
  NodeInfoRetentionProb = rep(NodeInfoRetentionProb,nrow(SDDprob))
InfoDecayNodes = which(HaveInfo == 1 & is.na(NodeInfoPersistenceSteps) & NodeInfoRetentionProb < 1)
if(length(InfoDecayNodes) > 0)
  HaveInfo[InfoDecayNodes] = rbinom(n = length(InfoDecayNodes),size = 1,prob = NodeInfoRetentionProb[InfoDecayNodes])

###Update info vector for any info spread (if SEAM supplied)
 ###Only zero values updated here so information can refresh nodes that lost information
 if(is.matrix(SEAM) == T)
  {
  RandSEAM[] <- rbinom(n=nrow(SDDprob)^2, size=1, prob = SEAM*apply(Detected,1,max))
  InfoTransferred = ifelse(colSums(RandSEAM)>0,1,0)
  HaveInfo[HaveInfo == 0] = InfoTransferred[HaveInfo == 0]
  }
 
 
 ###Add invasion resulting from colonisation from external sources
 if(OngoingExternalInvasion == T)
  {
  if(is.matrix(InvasionRisk) == F)
   ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
  if(is.matrix(InvasionRisk) == T)
   ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,timestep])
  Invaded[Invaded == 0] = ExternalInvasion[Invaded==0]
  
  ###If no start pop size all land uses invaded
  if(is.na(IncursionStartPop) == T) 
    N <- ifelse(N + ExternalInvasion>NodeK,NodeK,N + ExternalInvasion)
  
  ###If start pop size given, share individuals in proportion to available hosts/sites for each land use
  if(is.na(IncursionStartPop) == F) 
    {
    InVector = cbind(ExternalInvasion,NodeK-N0)
    LUExternalInvasion <- t(apply(InVector,1,FUN = SampleVector))
    N <- ifelse(N + LUExternalInvasion>NodeK,NodeK,N + LUExternalInvasion)
    }
  
  }

  if(UsePathogen)
    {
    PathogenStateFlat <- PathogenEngine$Reconcile(PathogenStateFlat, N, PathogenContext)
    PathogenStep <- PathogenEngine$Step(PathogenStateFlat, N, timestep, PathogenContext)
    PathogenStateFlat <- PathogenStep$State
    N[] <- PathogenStep$N
    for(ss in seq_along(PathogenStates)) PathogenStateResultsLoop[,,ss,timestep] <- matrix(PathogenStateFlat[,ss], nrow = nrow(SDDprob), ncol = Nlanduses)
    PathogenDetectionP <- PathogenEngine$Resolve(PathogenOriginal$DetectionProb, timestep, PathogenContext, "DetectionProb")
    p_not <- (1-PathogenDetectionP)^PathogenStateFlat[,"I"]
    p_node <- 1-apply(matrix(p_not,nrow=nrow(SDDprob),ncol=Nlanduses),1,prod)
    PathogenDetectedNow <- rbinom(nrow(SDDprob),1,p_node)
    PathogenDetectedResultsLoop[,timestep] <- PathogenDetectedNow
    if(isTRUE(PathogenOriginal$DetectionTriggersInfo))
      {
      if(UseInfoPersistence == T) LastKnownPresence[PathogenDetectedNow == 1] <- timestep
      HaveInfo[HaveInfo == 0] <- PathogenDetectedNow[HaveInfo == 0]
      }
    }
 
  ###Add nodes with information resulting from external sources
  if(OngoingExternalInfo == T)
    {
    if(is.matrix(ExternalInfoProb) == F)
      ExternalInfo = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb)
    if(is.matrix(ExternalInfoProb) == T)
      ExternalInfo = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb[,timestep])
    HaveInfo[HaveInfo == 0] = ExternalInfo[HaveInfo==0]
     }
 ###Update infestation vector
 Invaded = ifelse(N>0,1,0)
 head(Invaded)
 ###Record nodes adopting management
 ManagingResultsLoop[,,timestep] = Managing
  
 ###Record infested nodes
 InvasionResultsLoop[,,timestep] = Invaded

 ###Record populations
 PopulationResultsLoop[,,timestep] = N

 ###Select new nodes where infestation detected
 LUdetectionProb = 1-(1-NodeDetectionProb)^(N)
 NewHaveInfo = rbinom(1:nrow(SDDprob),size = 1,prob = 1-apply(1-LUdetectionProb,1,prod))
 
 
 ###Record newly detected infestations as known local presence
 if(UseInfoPersistence == T)
   {
   KnownPresence = which(NewHaveInfo == 1)
   if(length(KnownPresence) > 0)
     LastKnownPresence[KnownPresence] = timestep
   }
 
 ###Add newly detected infestations to info vector
 ###Only zero values updated here so information can refresh nodes that lost information
 HaveInfo[HaveInfo==0] = NewHaveInfo[HaveInfo==0]  
 
 ###Record detection status
 DetectedResultsLoop[,,timestep] = HaveInfo*Invaded 
 }
 if(UsePathogen) return(list(Invasion = InvasionResultsLoop, Population = PopulationResultsLoop, Managing = ManagingResultsLoop, Detected = DetectedResultsLoop, PathogenState = PathogenStateResultsLoop, PathogenDetected = PathogenDetectedResultsLoop))
 simplify2array(list(InvasionResultsLoop, PopulationResultsLoop, ManagingResultsLoop, DetectedResultsLoop), higher = TRUE)
}

###Use a common PSOCK/parLapply architecture across parallel INApest variants.
###Static scheduling and L'Ecuyer-CMRG worker streams support reproducible
###parallel simulations when the caller fixes the R seed.
if(n_cores == 1L)
  {
  PermutationResults <- lapply(seq_len(Nperm), PermutationWorker)
  } else {
  cluster <- parallel::makeCluster(n_cores, type = "PSOCK")
  on.exit(if(inherits(cluster, "cluster")) parallel::stopCluster(cluster), add = TRUE)
  parallel::clusterSetRNGStream(cluster)
  PermutationResults <- parallel::parLapply(cluster, seq_len(Nperm), PermutationWorker)
  parallel::stopCluster(cluster)
  cluster <- NULL
  }
if(!UsePathogen) {
  PermOut <- simplify2array(PermutationResults, higher = TRUE)
  if(length(dim(PermOut)) == 4) dim(PermOut) = c(dim(PermOut),1)
  InvasionResults <- array(PermOut[,,,1,,drop=FALSE],dim = c(nrow(SDDprob),Nlanduses,Ntimesteps,Nperm))
  PopulationResults <- array(PermOut[,,,2,,drop=FALSE],dim = c(nrow(SDDprob),Nlanduses,Ntimesteps,Nperm))
  ManagingResults <- array(PermOut[,,,3,,drop=FALSE],dim = c(nrow(SDDprob),Nlanduses,Ntimesteps,Nperm))
  DetectedResults <- array(PermOut[,,,4,,drop=FALSE],dim = c(nrow(SDDprob),Nlanduses,Ntimesteps,Nperm))
} else {
  InvasionResults <- simplify2array(lapply(PermutationResults, `[[`, "Invasion"), higher = TRUE)
  PopulationResults <- simplify2array(lapply(PermutationResults, `[[`, "Population"), higher = TRUE)
  ManagingResults <- simplify2array(lapply(PermutationResults, `[[`, "Managing"), higher = TRUE)
  DetectedResults <- simplify2array(lapply(PermutationResults, `[[`, "Detected"), higher = TRUE)
  if(Nperm == 1L) {
    dim(InvasionResults) <- c(nrow(SDDprob),Nlanduses,Ntimesteps,1L); dim(PopulationResults) <- c(nrow(SDDprob),Nlanduses,Ntimesteps,1L)
    dim(ManagingResults) <- c(nrow(SDDprob),Nlanduses,Ntimesteps,1L); dim(DetectedResults) <- c(nrow(SDDprob),Nlanduses,Ntimesteps,1L)
  }
  PathogenStates <- PathogenEngine$States
  PathogenStateResults <- array(0L, dim = c(nrow(SDDprob),Nlanduses,length(PathogenStates),Ntimesteps,Nperm), dimnames = list(NULL,NULL,PathogenStates,NULL,NULL))
  PathogenDetectedResults <- array(0L,dim=c(nrow(SDDprob),Ntimesteps,Nperm))
  for(pp in seq_len(Nperm)) { PathogenStateResults[,,,,pp] <- PermutationResults[[pp]]$PathogenState; PathogenDetectedResults[,,pp] <- PermutationResults[[pp]]$PathogenDetected }
}
###########################################################
###End of Simulation
###########################################################

###########################################################
###Save results for post-hoc analyses
###########################################################
###ModelName used to generate filenames
###Use standard format for ease of reading results to produce heat maps 
###and conduct post-hoc stats comparing managment settings/scenarios 
if(is.na(OutputDir) == T)
	OutputDir = ""
FileNameStem = paste0(OutputDir,ModelName)

###These are 3D arrays with dimensions (Nodes,Timesteps,Realisations)
saveRDS(ManagingResults, paste0(FileNameStem,"InfoLargeOut.rds"))
saveRDS(PopulationResults, paste0(FileNameStem,"PopulationLargeOut.rds"))
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))
if(UsePathogen) {
  saveRDS(PathogenStateResults, paste0(FileNameStem,"PathogenStateLargeOut.rds"))
  saveRDS(PathogenDetectedResults, paste0(FileNameStem,"PathogenDetectedLargeOut.rds"))
}

##########################################################
###Store annual node-level invasion probs for heat maps
###and estimation of invasion threat to other regions
##########################################################

InvasionProb = matrix(ncol = Ntimesteps, nrow = nrow(SDDprob))
for(timestep in 1:Ntimesteps)
{
TimestepData = InvasionResults[,,timestep,,drop=FALSE]
dim(TimestepData)
NodeInvaded <- apply(TimestepData, c(1,4), max)
if(is.null(dim(NodeInvaded))) NodeInvaded = matrix(NodeInvaded,nrow=nrow(SDDprob),ncol=Nperm)
InvasionProb[,timestep] = rowSums(NodeInvaded)/Nperm
}
saveRDS(InvasionProb, paste0(FileNameStem,"InvasionProb.rds"))
if(DoPlots == T)
{
###########################################################
###Produce summary figs when processing completed
###########################################################

Title = ModelName


###Change in total population with time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

PopulationSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(PopulationSummary) = c("Realisation",   "Timestep",  "NodesInfested")

#if(is.matrix(K) == TRUE)
#    inv_K <- 1 / colSums(K)
inv_K
for(perm in 1:Nperm)
{
PopulationData = PopulationResults[,,,perm]
dim(PopulationData)
NodesInfested = apply(PopulationData,3,sum)*inv_K
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesInfested)
PopulationSummary = rbind(PopulationSummary,Results)
}


Filename = paste0(FileNameStem,"PopulationRaw.png")
png(Filename)
plot(PopulationSummary$Timestep,PopulationSummary$NodesInfested,ylim = c(0,1),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Total population (proportion of K)", main = Title)

for(perm in 1:Nperm)
{
Sub = PopulationSummary[PopulationSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesInfested,col  = perm)
}
dev.off()

Quantiles = as.data.frame(aggregate(PopulationSummary$NodesInfested, by = list(PopulationSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"PopulationSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,1), xlab = "Time since incursion detected (timesteps)",
ylab = "Total population (proportion of K)", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()

###Change in number of nodes infested with time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

InvasionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(InvasionSummary) = c("Realisation",   "Timestep",  "NodesInfested")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,,perm]
NodesInfested = colSums(apply(InvasionData,c(1,3),max))

Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesInfested)
InvasionSummary = rbind(InvasionSummary,Results)
}


Filename = paste0(FileNameStem,"InvasionRaw.png")
png(Filename)
plot(InvasionSummary$Timestep,InvasionSummary$NodesInfested,ylim = c(0,max(InvasionSummary$NodesInfested)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Number of nodes infested", main = Title)

for(perm in 1:Nperm)
{
Sub = InvasionSummary[InvasionSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesInfested,col  = perm)
}
dev.off()

Quantiles = as.data.frame(aggregate(InvasionSummary$NodesInfested, by = list(InvasionSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"InvasionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Number of nodes infested", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


###Change in number of farms managing through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided

ManagingSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(ManagingSummary) = c("Realisation",   "Timestep",  "NodesManaging")

for(perm in 1:Nperm)
{
ManagingData = ManagingResults[,,,perm]
dim(ManagingData)
NodesManaging = apply(ManagingData,3,sum)/Nlanduses
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesManaging)
ManagingSummary = rbind(ManagingSummary,Results)
}



 
Filename = paste0(FileNameStem,"ManagingRaw.png")
png(Filename)
plot(ManagingSummary$Timestep,ManagingSummary$NodesManaging,ylim = c(0,max(ManagingSummary$NodesManaging)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes under management", main = Title)

for(perm in 1:Nperm)
{
Sub = ManagingSummary[ManagingSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesManaging,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(ManagingSummary$NodesManaging, by = list(ManagingSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"ManagingSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes under management", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


###Change in number of known extant infestations through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

DetectedSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedSummary) = c("Realisation",   "Timestep",  "NodesDetected")

for(perm in 1:Nperm)
{
DetectedData = DetectedResults[,,,perm]
dim(DetectedData)
NodesDetected = colSums(apply(DetectedData,c(1,3),max))
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,NodesDetected)
DetectedSummary = rbind(DetectedSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedRaw.png")
png(Filename)
plot(DetectedSummary$Timestep,DetectedSummary$NodesDetected,ylim = c(0,max(DetectedSummary$NodesDetected)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes pest detected", main = Title)

for(perm in 1:Nperm)
{
Sub = DetectedSummary[DetectedSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$NodesDetected,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(DetectedSummary$NodesDetected, by = list(DetectedSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"DetectedSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Nodes pest detected", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()


###Change in proportion of extant infestations detected through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided
 
DetectedProportionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedProportionSummary) = c("Realisation",   "Timestep",  "DetectedProportion")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,,perm]
DetectedData = DetectedResults[,,,perm]
NodesDetected = colSums(apply(DetectedData,c(1,3),max))
NodesInvaded = colSums(apply(InvasionData,c(1,3),max))
DetectedProportion = NodesDetected/NodesInvaded
DetectedProportion[is.na(DetectedProportion)==T] = 1
Realisation = perm 
Timestep = 1:Ntimesteps
Results = data.frame(Realisation,Timestep,DetectedProportion)
DetectedProportionSummary = rbind(DetectedProportionSummary,Results)
}


Filename = paste0(FileNameStem,"DetectedProportionRaw.png")
png(Filename)
plot(DetectedProportionSummary$Timestep,DetectedProportionSummary$DetectedProportion,ylim = c(0,max(DetectedProportionSummary$DetectedProportion)),pch = NA
, xlab = "Time since incursion detected (timesteps)",
ylab = "Proportion of infested nodes detected", main = Title)

for(perm in 1:Nperm)
{
Sub = DetectedProportionSummary[DetectedProportionSummary$Realisation == perm,]
lines(Sub$Timestep,Sub$DetectedProportion,col  = perm)
}
dev.off()
Quantiles = as.data.frame(aggregate(DetectedProportionSummary$DetectedProportion, by = list(DetectedProportionSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

Filename = paste0(FileNameStem,"DetectedProportionSummary.png")
png(Filename)
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (timesteps)",
ylab = "Proportion infested nodes detected", main = Title)
lines(Quantiles[,1],Yvals[,2],lwd = 3)
lines(Quantiles[,1],Yvals[,1],lwd = 3,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 3,col = 2)
dev.off()
}
}

################################################################
################################################################
###End of function
################################################################
################################################################
