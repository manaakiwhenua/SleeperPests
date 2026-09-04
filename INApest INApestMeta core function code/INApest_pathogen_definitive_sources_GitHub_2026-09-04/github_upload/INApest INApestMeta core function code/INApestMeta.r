###############################################################################
### INApestMeta -- abundance-by-location invasion simulation engine
###
### INApestMeta replaces binary occupancy with host/pest abundance at each
### spatial node. Local populations survive and reproduce, propagules move by
### short- and long-distance pathways, and establishment is limited by habitat
### and carrying capacity.
###
### The simulation keeps biological state separate from information and response:
### surveillance determines what is known, informed nodes may be managed, and
### management can alter mortality, fecundity and spread. An optional SIS/SIR/
### SEIR pathogen process can be coupled to the node-level host abundance.
###############################################################################

# Default abundance-node local growth, dispersal and establishment process.
local.dynamics <- function(
    sddprob = SDDprob,
    nodepropaguleproduction = NodePropaguleProduction,
    nodeenvestabprob = NodeEnvEstabProb,
    n = N0,
    lddprob = LDDprob,
    lddrate = LDDrate,
    k_is_0 = K_is_0,
    nodeK = NodeK,
    nodepropaguleestablishment = NodePropaguleEstablishment,
    nodespreadreduction = NodeSpreadReduction,
    nodefecundityreduction = 0,
    managing = Managing,
    maxinteger = MaxInteger
) {
  EffectiveReproductivePopulation <- n * (1 - nodefecundityreduction * managing)
  Propagules <- rpois(nrow(sddprob), nodepropaguleproduction * EffectiveReproductivePopulation)

  # Treat dispersal entries as per-propagule source-to-target probabilities.
  # LDDrate first splits the realised propagule count; management then thins the
  # LDD branch. Any residual
  # probability of 1 - rowSum represents propagules that leave the modelled system.
  n_nodes <- nrow(sddprob)
  # Allocate realised integer propagules among destinations and outside loss.
  AllocateDispersers <- function(counts, probmat, name) {
    incoming <- numeric(n_nodes)
    if(!is.matrix(probmat)) return(incoming)
    rs <- rowSums(probmat)
    if(any(!is.finite(probmat)) || any(probmat < 0) || any(rs > 1 + 1e-12))
      stop(name, " must contain finite non-negative probabilities with every source-row sum <= 1")
    active <- which(counts > 0 & rs > 0)
    if(!length(active)) return(incoming)
    for(ii in active) {
      if(counts[ii] < maxinteger) {
        probs <- c(probmat[ii,], max(0, 1-rs[ii]))
        z <- rmultinom(1, size=counts[ii], prob=probs)
        incoming <- incoming + z[seq_len(n_nodes),1]
      } else {
        incoming <- incoming + counts[ii] * probmat[ii,]
      }
    }
    incoming
  }

  if(lddrate <= 0) LDDcount <- integer(n_nodes)
  else if(lddrate >= 1) LDDcount <- Propagules
  else LDDcount <- rbinom(n_nodes, size=Propagules, prob=lddrate)
  SDDcount <- Propagules-LDDcount
  Pin <- AllocateDispersers(SDDcount, sddprob, "SDDprob")

  Qin <- numeric(n_nodes)
  if(is.matrix(lddprob)) {
    keep <- pmin(1,pmax(0,1-nodespreadreduction*managing))
    LDDkept <- if(all(keep >= 1)) LDDcount else rbinom(n_nodes,size=LDDcount,prob=keep)
    Qin <- AllocateDispersers(LDDkept, lddprob, "LDDprob")
  }

  # Recruitment is limited by unoccupied capacity in the receiving node.
  Nout <- ifelse(
    k_is_0,
    0,
    n + rbinom(
      nrow(sddprob),
      nodeK - n,
      1 - exp(-nodepropaguleestablishment * nodeenvestabprob * (Pin + Qin))
    )
  )

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

# Resolve custom LocalDynamics arguments for the current timestep.
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
    Model = c("SIS", "SIR", "SEIR", "Binary"),  # Pathogen state model: Binary, SIS, SIR or SEIR
    Beta = 0,                                   # Transmission-rate parameter
    RecoveryProb = 0,                           # Probability an infectious host recovers
    ProgressionProb = 1,                        # Probability an exposed host becomes infectious
    PathogenMortalityProb = 0,                  # Probability infection kills an infectious host
    ImmunityLossProb = 0,                       # Probability a recovered host becomes susceptible again
    InitialInfected = 0,                        # Starting infectious host count/specification
    InitialExposed = 0,                         # Starting exposed host count/specification
    InitialRecovered = 0,                       # Starting recovered host count/specification
    IntroductionProb = 0,                       # Probability an existing susceptible host is externally infected
    IntroductionNumber = 1,                     # Hosts infected when an external introduction occurs
    ContactMatrix = NULL,                       # Source-to-target pathogen contact matrix
    Transmission = c("frequency", "density"),   # Frequency- or density-dependent transmission
    DensityScale = 1,                           # Scale for density-dependent transmission
    InitialPresent = 0,                         # Starting binary pathogen occupancy
    ClearanceProb = 0,                          # Probability pathogen clears without host loss
    TransmissionProb = NULL,                    # Binary source-to-target transmission probability
    PathogenHostExtinctionProb = 0,             # Probability pathogen removes infected host occupancy
    DetectionProb = 0,                          # Background-surveillance detection probability
    DetectionTriggersInfo = FALSE) {            # Whether pathogen detection updates INApest information

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
    I <- resolve_count(Spec$InitialInfected, "InitialInfected")
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


INApestMeta = function(
ModelName,                                      # Model and output name
Nperm,                                          # Number of stochastic simulation runs
Ntimesteps,                                     # Timesteps in each simulation
LocalDynamics = local.dynamics,                 # Local population growth, movement and management function
LocalDynamicsArgs = list(),                     # Named extra arguments passed to LocalDynamics
Pathogen = NULL,                                # Optional pathogen process specification
DetectionProb,                                  # Background-surveillance detection probability
DetectionSD = NULL,                             # Variation in background detection probability
ManageProb,                                     # Management probability when information is available
ManageSD = NULL,                                # Variation in management probability
MortalityProb,                                  # Management-driven host mortality probability
MortalitySD = NULL,                             # Variation in management mortality probability
FecundityReduction = 0,                         # Proportional reduction in reproduction under management
SpreadReduction,                                # Proportional reduction in dispersal under management
SpreadReductionSD = NULL,                       # Variation in spread reduction
InitialPopulation = NA,                         # Starting host abundance by node
InitBioP = NA,                                  # Proportion of nodes initially invaded
InvasionRisk = NA,                              # External invasion probability or node weighting
InitialInfo = NA,                               # Starting information state
InitInfoP = NA,                                 # Proportion of nodes initially with information
ExternalInfoProb = NA,                          # Information arriving from outside the modelled system
InfoRetentionProb = 1,                          # Probability existing information persists one timestep
InfoPersistenceSteps = NA,                      # Programmed timesteps information persists after evidence
EnvEstabProb = 1,                               # Environmental establishment probability
Survival = 1,                                   # Host survival probability between timesteps
K,                                              # Host carrying capacity by node
PropaguleProduction,                            # Propagules produced per host
PropaguleEstablishment,                         # Establishment probability for arriving propagules
IncursionStartPop=NA,                           # Host abundance assigned to new external incursions
SDDprob,                                        # Short-distance source-to-target dispersal
SEAM = 0,                                       # Information-transfer adjacency between nodes
LDDprob = NA,                                   # Long-distance source-to-target dispersal
LDDrate = 0,                                    # Fraction of propagules entering long-distance dispersal
OngoingExternalInvasion = F,                    # Allow new host incursions after initialisation
OngoingExternalInfo = F,                        # Allow new external information after initialisation
OutputDir = NA,                                 # Directory for saved outputs
DoPlots = TRUE,                                 # Legacy plotting option; plotting is post-processing
ExternalPathogenStateProb = NULL,               # Pathogen-state distribution for external host arrivals
InformationAcquisition = NULL,                  # Local evidence used for information: host, pathogen or both
InfoTriggeredDetectionProb = 0,                 # Detection probability where information already exists
InfoTriggeredDetectionSD = NULL,                # Variation in information-triggered detection
ReturnResults = FALSE                           # Return the in-memory result object
)
{

# ---------------------------------------------------------------------------
# Set up and validate the abundance-node simulation.
# ---------------------------------------------------------------------------
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
LocalDynamicsHasPreMortalityN <- "pre_mortality_n" %in% names(formals(LocalDynamics))
PathogenOriginal <- Pathogen
if(!is.null(PathogenOriginal) && !inherits(PathogenOriginal, "INApestPathogen")) stop("Pathogen must be NULL or an object returned by INApestPathogen()")
UsePathogen <- !is.null(PathogenOriginal)
if(UsePathogen) {
  PathogenEngine <- PathogenOriginal$Engine
  PathogenContext <- list(n_nodes = nrow(SDDprob), Ntimesteps = Ntimesteps)
  PathogenEngine$Validate(PathogenContext)
  force(PathogenEngine); force(PathogenContext)
}

# Choose which direct local biological evidence can create/refresh HaveInfo.
# NULL is the backwards-compatible pathway: host detection always informs, and
# pathogen detection informs only when Pathogen$DetectionTriggersInfo is TRUE.
ExplicitInformationAcquisition <- !is.null(InformationAcquisition)
if(ExplicitInformationAcquisition) {
  if(length(InformationAcquisition) != 1L || is.na(InformationAcquisition))
    stop("InformationAcquisition must be one of: host, pathogen, both")
  InformationAcquisition <- match.arg(as.character(InformationAcquisition), c("host", "pathogen", "both"))
  if(InformationAcquisition == "pathogen" && !UsePathogen)
    stop("InformationAcquisition = 'pathogen' requires Pathogen to be supplied")
  HostInformationAcquisition <- InformationAcquisition %in% c("host", "both")
  PathogenInformationAcquisition <- InformationAcquisition %in% c("pathogen", "both")
} else {
  HostInformationAcquisition <- TRUE
  PathogenInformationAcquisition <- UsePathogen && isTRUE(PathogenOriginal$DetectionTriggersInfo)
  InformationAcquisition <- if(PathogenInformationAcquisition) "both" else "host"
}
force(HostInformationAcquisition); force(PathogenInformationAcquisition); force(InformationAcquisition)

# Optional pathogen-state composition of accepted external host immigrants.
# This changes only pathogen-state bookkeeping: the existing host invasion draw,
# IncursionStartPop and carrying-capacity clipping remain authoritative for N.
ExternalPathogenStateProbResolved <- NULL
if(!is.null(ExternalPathogenStateProb)) {
  if(!UsePathogen)
    stop("ExternalPathogenStateProb requires Pathogen to be supplied")
  if(!is.numeric(ExternalPathogenStateProb) || is.null(names(ExternalPathogenStateProb)) ||
     any(!nzchar(names(ExternalPathogenStateProb))) || anyDuplicated(names(ExternalPathogenStateProb)))
    stop("ExternalPathogenStateProb must be a uniquely named numeric vector of pathogen-state probabilities")
  if(any(!names(ExternalPathogenStateProb) %in% PathogenEngine$States))
    stop("ExternalPathogenStateProb names must be pathogen states: ", paste(PathogenEngine$States, collapse = ", "))
  if(any(!is.finite(ExternalPathogenStateProb)) || any(ExternalPathogenStateProb < 0))
    stop("ExternalPathogenStateProb values must be finite and non-negative")
  ExternalPathogenStateProbResolved <- setNames(rep(0, length(PathogenEngine$States)), PathogenEngine$States)
  ExternalPathogenStateProbResolved[names(ExternalPathogenStateProb)] <- ExternalPathogenStateProb
  if(abs(sum(ExternalPathogenStateProbResolved) - 1) > 1e-10)
    stop("ExternalPathogenStateProb must sum to 1")
}
force(ExternalPathogenStateProbResolved)


# Allow SDD and LDD connectivity to vary through time
# static connectivity matrices must preserve the
# their time-varying forms.
if(is.matrix(SDDprob) && nrow(SDDprob) != ncol(SDDprob))
  stop("SDDprob matrix must be square")
if(is.matrix(LDDprob) && (nrow(LDDprob) != nrow(SDDprob) || ncol(LDDprob) != nrow(SDDprob)))
  stop("LDDprob matrix must have dimensions nodes x nodes")
if(!is.numeric(LDDrate) || length(LDDrate) != 1L || !is.finite(LDDrate) || LDDrate < 0 || LDDrate > 1)
  stop("LDDrate must be one finite probability between 0 and 1")
if(length(dim(SDDprob)) == 3 && (dim(SDDprob)[1] != dim(SDDprob)[2] || dim(SDDprob)[3] != Ntimesteps))
  stop("SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
if(length(dim(LDDprob)) == 3 && (dim(LDDprob)[1] != nrow(SDDprob) || dim(LDDprob)[2] != nrow(SDDprob) || dim(LDDprob)[3] != Ntimesteps))
  stop("LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
# Allow management-induced fecundity reduction to vary by node and through time
if(is.matrix(FecundityReduction)) {
  if(nrow(FecundityReduction) != nrow(SDDprob) || ncol(FecundityReduction) != Ntimesteps)
    stop("FecundityReduction matrix must have dimensions nodes x Ntimesteps")
} else if(!(length(FecundityReduction) == 1 || length(FecundityReduction) == nrow(SDDprob))) {
  stop("FecundityReduction must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
}
if(any(!is.finite(FecundityReduction)) || any(FecundityReduction < 0) || any(FecundityReduction > 1))
  stop("FecundityReduction values must be between 0 and 1")

# Optional second surveillance stream activated only by information available
# before the current surveillance round. Contract mirrors DetectionProb.
if(is.matrix(InfoTriggeredDetectionProb)) {
  if(nrow(InfoTriggeredDetectionProb) != nrow(SDDprob) || ncol(InfoTriggeredDetectionProb) != Ntimesteps)
    stop("InfoTriggeredDetectionProb matrix must have dimensions nodes x Ntimesteps")
} else if(!(length(InfoTriggeredDetectionProb) == 1 || length(InfoTriggeredDetectionProb) == nrow(SDDprob))) {
  stop("InfoTriggeredDetectionProb must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
}
if(any(!is.finite(InfoTriggeredDetectionProb)) || any(InfoTriggeredDetectionProb < 0) || any(InfoTriggeredDetectionProb > 1))
  stop("InfoTriggeredDetectionProb values must be between 0 and 1")
if(!is.null(InfoTriggeredDetectionSD)) {
  if(!is.numeric(InfoTriggeredDetectionSD) || !(length(InfoTriggeredDetectionSD) == 1 || length(InfoTriggeredDetectionSD) == nrow(SDDprob)) ||
     any(!is.finite(InfoTriggeredDetectionSD)) || any(InfoTriggeredDetectionSD < 0))
    stop("InfoTriggeredDetectionSD must be NULL, a non-negative scalar, or a non-negative vector of length nodes")
}

# Allow information retention to vary by node and through time
if(is.matrix(InfoRetentionProb) == T && (nrow(InfoRetentionProb) != nrow(SDDprob) || ncol(InfoRetentionProb) != Ntimesteps))
  stop("InfoRetentionProb matrix must have dimensions nodes x Ntimesteps")
if(is.matrix(InfoRetentionProb) == F && !(length(InfoRetentionProb) == 1 || length(InfoRetentionProb) == nrow(SDDprob)))
  stop("InfoRetentionProb must be a single value, vector of length nodes, or matrix nodes x Ntimesteps")
if(any(is.na(InfoRetentionProb)) || any(InfoRetentionProb < 0) || any(InfoRetentionProb > 1))
  stop("InfoRetentionProb values must be between 0 and 1")

# Allow programmed information persistence after last known local presence to vary by node and through time
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
if(is.matrix(K) == FALSE)
{
K_is_0 <- K<=0
inv_K <- 1 / sum(K)
NodeK = K
}

if(is.matrix(PropaguleProduction) == FALSE)
 NodePropaguleProduction = PropaguleProduction

if(is.matrix(PropaguleEstablishment) == FALSE)
  NodePropaguleEstablishment = PropaguleEstablishment

if(is.matrix(EnvEstabProb) == F)
  NodeEnvEstabProb <- EnvEstabProb

if(is.matrix(Survival) == F)
  NodeSurvival <- Survival


# Validate the socioeconomic information-transfer network when supplied.
if(is.matrix(SEAM) == T)
     {
     diag(SEAM) = 0
     RandSEAM <- matrix(NA,nrow = nrow(SDDprob),ncol=nrow(SDDprob))
     }

# Assign standard deviation value to management in no value provided
if(is.null(ManageSD) == T)
	ManageSD = mean(ManageProb)/10
if(is.null(SpreadReductionSD) == T)
	SpreadReductionSD = (1-mean(SpreadReduction))/10
if(is.null(DetectionSD) == T)
	DetectionSD = mean(DetectionProb)/10
if(is.null(InfoTriggeredDetectionSD) == T)
	InfoTriggeredDetectionSD = mean(InfoTriggeredDetectionProb)/10
UseInfoTriggeredSurveillance = any(InfoTriggeredDetectionProb != 0) || any(InfoTriggeredDetectionSD != 0)
if(is.null(MortalitySD) == T)
	MortalitySD = mean(MortalityProb)/10


###########################################################
### Start of simulation
###########################################################
    
# Run one stochastic realisation. Function arguments and local helpers are
# captured in this closure, avoiding fragile manual worker export lists.

# ---------------------------------------------------------------------------
# Run one complete stochastic host-pathogen history.
# ---------------------------------------------------------------------------
PermutationWorker <- function(i_perm)
  {
  # Max integer for propagule dispersal using rmultinom
  MaxInteger <- .Machine$integer.max  

  # Set initial dispersal connectivity
  NodeSDDprob = SDDprob
  if(length(dim(SDDprob)) == 3)
    NodeSDDprob = SDDprob[,,1]
  NodeLDDprob = LDDprob
  if(length(dim(LDDprob)) == 3)
    NodeLDDprob = LDDprob[,,1]
  NodeFecundityReduction <- if(is.matrix(FecundityReduction)) FecundityReduction[,1] else FecundityReduction
  if(length(NodeFecundityReduction) == 1)
    NodeFecundityReduction <- rep(NodeFecundityReduction, nrow(SDDprob))
  
  InvasionResultsLoop <- array(dim = c(nrow(SDDprob),Ntimesteps))
  PopulationResultsLoop <- InvasionResultsLoop
  ManagingResultsLoop <- InvasionResultsLoop
  DetectedResultsLoop <- InvasionResultsLoop
  BackgroundDetectedResultsLoop <- matrix(0L, nrow = nrow(SDDprob), ncol = Ntimesteps)
  InfoTriggeredDetectedResultsLoop <- matrix(0L, nrow = nrow(SDDprob), ncol = Ntimesteps)
  BackgroundDetectionProbabilityResultsLoop <- matrix(0, nrow = nrow(SDDprob), ncol = Ntimesteps)
  InfoTriggeredDetectionProbabilityResultsLoop <- matrix(0, nrow = nrow(SDDprob), ncol = Ntimesteps)
  InformationStateBeforeSurveillanceResultsLoop <- matrix(0L, nrow = nrow(SDDprob), ncol = Ntimesteps)
  HaveInfoResultsLoop <- matrix(0L, nrow = nrow(SDDprob), ncol = Ntimesteps)
  if(UsePathogen) {
    PathogenStates <- PathogenEngine$States
    PathogenStateResultsLoop <- array(0L, dim = c(nrow(SDDprob), length(PathogenStates), Ntimesteps),
                                      dimnames = list(NULL, PathogenStates, NULL))
    PathogenDetectedResultsLoop <- matrix(0L, nrow = nrow(SDDprob), ncol = Ntimesteps)
  }

  
  # Use first-timestep carrying capacity when K varies through time.
  if(is.matrix(K) == TRUE)
    {
    K_is_0 <- K[,1]<=0
    inv_K <- 1 / sum(K[,1])
    NodeK = K[,1] 
    }      
  
  
# Initialise host abundance from explicit starting values or sampled invasion settings.
InitBio = rep(0,times = nrow(SDDprob))
# Treat scalar NA as the default sentinel, including in one-node models.
HasInitialPopulation <- length(InitialPopulation) == nrow(SDDprob) && !anyNA(InitialPopulation)
if(HasInitialPopulation)
  InitBio = InitialPopulation

if(!HasInitialPopulation)
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
  InitBio[Infested] = 1
if(is.na(IncursionStartPop) == F)
  InitBio[Infested] = IncursionStartPop
}

# Cap starting host abundance at local carrying capacity.
InitBio[InitBio > NodeK] = NodeK[InitBio > NodeK] 

# Set the working host abundance and initialise pathogen state when present.
N <- InitBio
if(UsePathogen)
  PathogenState <- PathogenEngine$Initial(N, PathogenContext)
if(sum(N) == 0 && OngoingExternalInvasion == F)
  warning("No initial populations and no future external invasions")

# Initialise response information independently of true host abundance.
# Without an initial information specification, all nodes start uninformed.
InitInfo = rep(0,times = nrow(SDDprob))
# Treat scalar NA as the InitialInfo sentinel, including in one-node models.
HasInitialInfo <- length(InitialInfo) == nrow(SDDprob) && !anyNA(InitialInfo)
if(HasInitialInfo || (is.na(InitInfoP) == F && InitInfoP>0) || is.na(sum(ExternalInfoProb)) == F )
{
if(!HasInitialInfo)
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
      # Use the first timestep of matrix ExternalInfoProb as initial information weights.
      if(is.na(InitInfoP) == F)
        Info = sample(1:nrow(SDDprob),size = ceiling(nrow(SDDprob)*InitInfoP),prob = ExternalInfoProb[,1])
      if(is.na(InitInfoP) == T)
        {
        Info = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb[,1])
        Info = which(Info == 1)
        }
      }
    }
  InitInfo[Info] = 1
  
  }
if(HasInitialInfo)
  InitInfo = InitialInfo  
}

# Draw node-level background-surveillance detection probabilities.
# Use static background-detection inputs when supplied by scalar or node.
if(is.matrix(DetectionProb)==FALSE &&(length(DetectionProb) == 1 ||length(DetectionProb) == nrow(SDDprob) ))
      {
      NodeDetectionProb = rnorm(DetectionProb,DetectionSD,n = nrow(SDDprob))
      NodeDetectionProb[NodeDetectionProb<0] = 0
      NodeDetectionProb[NodeDetectionProb>1] = 1
      }

# Use first-timestep detection inputs for initial surveillance when detection varies through time.
if(is.matrix(DetectionProb)==TRUE && nrow(DetectionProb) == nrow(SDDprob) && ncol(DetectionProb) == Ntimesteps)
      {
      NodeDetectionProb = rnorm(DetectionProb[,1],DetectionSD,n = nrow(SDDprob))
      NodeDetectionProb[NodeDetectionProb<0] = 0
      NodeDetectionProb[NodeDetectionProb>1] = 1
      }


# Resolve information-triggered per-individual detection probability without
# consuming any additional random numbers when the optional pathway is off.
if(UseInfoTriggeredSurveillance)
  {
  if(is.matrix(InfoTriggeredDetectionProb)==FALSE && (length(InfoTriggeredDetectionProb) == 1 || length(InfoTriggeredDetectionProb) == nrow(SDDprob)))
    {
    NodeInfoTriggeredDetectionProb = rnorm(InfoTriggeredDetectionProb,InfoTriggeredDetectionSD,n = nrow(SDDprob))
    NodeInfoTriggeredDetectionProb[NodeInfoTriggeredDetectionProb<0] = 0
    NodeInfoTriggeredDetectionProb[NodeInfoTriggeredDetectionProb>1] = 1
    }
  }

# Draw node-level management-adoption probabilities where information is available.
# Use static management-adoption inputs when supplied by scalar or node.
if(is.matrix(ManageProb)==FALSE &&(length(ManageProb) == 1 ||length(ManageProb) == nrow(SDDprob) ))
      {
      NodeManageProb = rnorm(ManageProb,ManageSD,n = nrow(SDDprob))
      NodeManageProb[NodeManageProb<0] = 0
      NodeManageProb[NodeManageProb>1] = 1
      }

# Draw the management effect on outgoing spread.
# Use static spread-reduction inputs when supplied by scalar or node.
if(is.matrix(SpreadReduction)==FALSE &&(length(SpreadReduction) == 1 ||length(SpreadReduction) == nrow(SDDprob) ))
      {
      NodeSpreadReduction = rnorm(SpreadReduction,SpreadReductionSD,n = nrow(SDDprob))
      NodeSpreadReduction[NodeSpreadReduction<0] = 0
      NodeSpreadReduction[NodeSpreadReduction>1] = 1
      }

# Draw management-driven host mortality probabilities.
# Use static management-mortality inputs when supplied by scalar or node.
if(is.matrix(MortalityProb)==FALSE &&(length(MortalityProb) == 1 ||length(MortalityProb) == nrow(SDDprob) ))
      {
      NodeMortalityProb = rnorm(MortalityProb,MortalitySD,n = nrow(SDDprob))
      NodeMortalityProb[NodeMortalityProb<0] = 0
      NodeMortalityProb[NodeMortalityProb>1] = 1
      }

# Record current host/pest presence from starting abundance.
Invaded = ifelse(InitBio>0,1,0) 

# Apply initial surveillance to the starting host population.
# Accepted detections can add response information.
InitDetection = rbinom(1:nrow(SDDprob),size = 1,prob = 1-(1-NodeDetectionProb)^InitBio)
if(HostInformationAcquisition)
  InitInfo[InitInfo == 0] = InitDetection[InitInfo == 0]
# Set the working information state for the simulation.
HaveInfo = InitInfo
InitPathogenDetection <- integer(nrow(SDDprob))
if(UsePathogen && PathogenInformationAcquisition)
  {
  PathogenDetectionP <- PathogenEngine$Resolve(PathogenOriginal$DetectionProb, 1L, PathogenContext, "DetectionProb")
  InitPathogenDetection <- rbinom(nrow(SDDprob), size = 1, prob = 1-(1-PathogenDetectionP)^PathogenState[,"I"])
  HaveInfo[HaveInfo == 0] <- InitPathogenDetection[HaveInfo == 0]
  }

# Track the most recent timestep with direct local evidence accepted by the
# selected information-acquisition rule. User/external information alone does
# not create a local-evidence clock.
LastKnownPresence = rep(NA,nrow(SDDprob))
if(UseInfoPersistence == T)
  {
  InitialKnownPresence = integer(0)
  if(HostInformationAcquisition)
    InitialKnownPresence = union(InitialKnownPresence, which(InitDetection == 1))
  if(PathogenInformationAcquisition && ExplicitInformationAcquisition)
    InitialKnownPresence = union(InitialKnownPresence, which(InitPathogenDetection == 1))
  # When InformationAcquisition is NULL, preserve the established initial pathogen
  # information-clock semantics; later pathogen detections still refresh the clock.
  if(length(InitialKnownPresence) > 0)
    LastKnownPresence[InitialKnownPresence] = 0
  }

    
  # run simulation

# ---------------------------------------------------------------------------
# Advance biology, information and response through time.
# ---------------------------------------------------------------------------
for(timestep in 1:Ntimesteps) 
  { 
 
  # Resolve short- and long-distance connectivity for the current timestep.
  if(length(dim(SDDprob)) == 3)
    NodeSDDprob = SDDprob[,,timestep]
  if(length(dim(LDDprob)) == 3)
    NodeLDDprob = LDDprob[,,timestep]
  if(is.matrix(FecundityReduction))
    NodeFecundityReduction <- FecundityReduction[,timestep]

  # Resolve time-varying environmental establishment inputs.
  if(is.matrix(EnvEstabProb) == T)
    NodeEnvEstabProb <- EnvEstabProb[,timestep]
   
  if(is.matrix(Survival) == T)
    NodeSurvival <- Survival[,timestep]
    
    
  # Resolve carrying capacity for the current timestep.
  if(is.matrix(K) == TRUE)
    {
    K_is_0 <- K[,timestep]<=0
    inv_K <- 1 / sum(K[,timestep])
    NodeK = K[,timestep] 
    }  

  # Resolve propagule production for the current timestep.
  if(is.matrix(PropaguleProduction) == TRUE)
    NodePropaguleProduction = PropaguleProduction[,timestep] 
  
  if(is.matrix(PropaguleEstablishment) == TRUE)
    NodePropaguleEstablishment = PropaguleEstablishment[,timestep]
      
  # Resolve time-varying background detection for the current timestep.
  if(is.matrix(DetectionProb)==TRUE && nrow(DetectionProb) == nrow(SDDprob) && ncol(DetectionProb) == Ntimesteps)
   	{	
   	NodeDetectionProb = rnorm(DetectionProb[,timestep],DetectionSD,n = nrow(SDDprob))
   	NodeDetectionProb[NodeDetectionProb<0] = 0
   	NodeDetectionProb[NodeDetectionProb>1] = 1
   	}

  # Draw information-triggered detection probabilities for informed nodes.
  # If information-triggered detection probability is time-varying, draw the
  # per-node value for this timestep. Static values were drawn once above.
  if(UseInfoTriggeredSurveillance && is.matrix(InfoTriggeredDetectionProb)==TRUE && nrow(InfoTriggeredDetectionProb) == nrow(SDDprob) && ncol(InfoTriggeredDetectionProb) == Ntimesteps)
    {
    NodeInfoTriggeredDetectionProb = rnorm(InfoTriggeredDetectionProb[,timestep],InfoTriggeredDetectionSD,n = nrow(SDDprob))
    NodeInfoTriggeredDetectionProb[NodeInfoTriggeredDetectionProb<0] = 0
    NodeInfoTriggeredDetectionProb[NodeInfoTriggeredDetectionProb>1] = 1
    }

  # Resolve time-varying management adoption for the current timestep.
  if(is.matrix(ManageProb)==TRUE && nrow(ManageProb) == nrow(SDDprob) && ncol(ManageProb) == Ntimesteps)
   	{	
   	NodeManageProb = rnorm(ManageProb[,timestep],ManageSD,n = nrow(SDDprob))
   	NodeManageProb[NodeManageProb<0] = 0
   	NodeManageProb[NodeManageProb>1] = 1
   	}

  # Resolve time-varying spread reduction for the current timestep.
  if(is.matrix(SpreadReduction)==TRUE && nrow(SpreadReduction) == nrow(SDDprob) && ncol(SpreadReduction) == Ntimesteps)
   	{	
   	NodeSpreadReduction = rnorm(SpreadReduction[,timestep],SpreadReductionSD,n = nrow(SDDprob))
   	NodeSpreadReduction[NodeSpreadReduction<0] = 0
   	NodeSpreadReduction[NodeSpreadReduction>1] = 1
   	}
  
  # Resolve time-varying management mortality for the current timestep.
  if(is.matrix(MortalityProb)==TRUE && nrow(MortalityProb) == nrow(SDDprob) && ncol(MortalityProb) == Ntimesteps)
      {
      NodeMortalityProb = rnorm(MortalityProb[,timestep],MortalitySD,n = nrow(SDDprob))
      NodeMortalityProb[NodeMortalityProb<0] = 0
      NodeMortalityProb[NodeMortalityProb>1] = 1
      }


  # Use current information to activate management at each node.
  Managing = rbinom(1:nrow(SDDprob),size = 1,prob = NodeManageProb*HaveInfo)
  
  # Identify occupied nodes currently known to the response system.
  Detected = Invaded*HaveInfo
  
  # Adjust starting population for natural and managed mortality.
  # The authoritative host draw is deliberately unchanged. When a pathogen is
  # active, realised deaths are conditionally separated into natural deaths and
  # management deaths, preserving exactly the same final N0 while exposing two
  # biologically meaningful pathogen-state event boundaries.
  NBeforeMortality <- N
  ManagementMortality = NodeMortalityProb*Managing
  N0 = rbinom(nrow(SDDprob),N,NodeSurvival*(1-ManagementMortality))

  ConditionalManagementMortality = rep(0,nrow(SDDprob))
  if(UsePathogen || UseInfoPersistence == T) {
    TotalMortalityProb = 1-NodeSurvival*(1-ManagementMortality)
    ManagementMortalityNodes = which((N-N0) > 0 & ManagementMortality > 0 & TotalMortalityProb > 0)
    if(length(ManagementMortalityNodes) > 0)
      ConditionalManagementMortality[ManagementMortalityNodes] =
        (NodeSurvival*ManagementMortality)[ManagementMortalityNodes]/TotalMortalityProb[ManagementMortalityNodes]
  }

  ManagementDeaths <- NaturalDeaths <- NULL
  if(UsePathogen) {
    TotalDeaths <- as.integer(N-N0)
    ManagementDeaths <- integer(nrow(SDDprob))
    NodesWithPossibleManagementDeaths <- which(TotalDeaths > 0 & ConditionalManagementMortality > 0)
    if(length(NodesWithPossibleManagementDeaths) > 0)
      ManagementDeaths[NodesWithPossibleManagementDeaths] <- rbinom(
        n = length(NodesWithPossibleManagementDeaths),
        size = TotalDeaths[NodesWithPossibleManagementDeaths],
        prob = ConditionalManagementMortality[NodesWithPossibleManagementDeaths]
      )
    NaturalDeaths <- TotalDeaths-ManagementDeaths
    NAfterNaturalMortality <- as.integer(N)-NaturalDeaths
    PathogenState <- PathogenEngine$Reconcile(PathogenState, NAfterNaturalMortality, PathogenContext)
    PathogenState <- PathogenEngine$Reconcile(PathogenState, N0, PathogenContext)
  }

  # Track known local presence from actual management mortality. Host-only runs
  # retain the previous conditional at-least-one-kill calculation exactly.
  if(UseInfoPersistence == T)
    {
    if(UsePathogen) {
      if(HostInformationAcquisition) {
        KnownPresence = which(ManagementDeaths > 0)
        if(length(KnownPresence) > 0)
          LastKnownPresence[KnownPresence] = timestep
      }
    } else {
      ManagementKillProb = 1-(1-ConditionalManagementMortality)^(N-N0)
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
    }
  if(sum(N0)<=0 )
    N = N0
  Pin <-0
  Qin <- 0  
    # natural dispersal
  if(sum(N0)>0 || (LocalDynamicsHasPreMortalityN && sum(NBeforeMortality)>0)) 
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
    managing = Managing,
    maxinteger = MaxInteger
  )
  LocalDynamicsFormals <- names(formals(LocalDynamics))
  # Optional read-only start-of-step host abundance for custom dynamics that
  # explicitly request it. It is not passed merely because a function has ...,
  # so existing custom/default LocalDynamics calls remain unchanged.
  if(LocalDynamicsHasPreMortalityN)
    CoreLocalDynamicsArgs$pre_mortality_n <- NBeforeMortality
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
    # Default and legacy LocalDynamics return total host abundance. Reconcile
    # at this boundary: positive net recruitment enters S; negative net change
    # is thinned from the current pathogen-state composition.
    if(UsePathogen)
      PathogenState <- PathogenEngine$Reconcile(PathogenState, N, PathogenContext)

  } 
 # Apply programmed stopping after last known local presence
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

# Allow information to decay after management and spread where no programmed stop is supplied
NodeInfoRetentionProb = InfoRetentionProb
if(is.matrix(InfoRetentionProb) == T)
  NodeInfoRetentionProb = InfoRetentionProb[,timestep]
if(length(NodeInfoRetentionProb) == 1)
  NodeInfoRetentionProb = rep(NodeInfoRetentionProb,nrow(SDDprob))
InfoDecayNodes = which(HaveInfo == 1 & is.na(NodeInfoPersistenceSteps) & NodeInfoRetentionProb < 1)
if(length(InfoDecayNodes) > 0)
  HaveInfo[InfoDecayNodes] = rbinom(n = length(InfoDecayNodes),size = 1,prob = NodeInfoRetentionProb[InfoDecayNodes])

# Update info vector for any info spread (if SEAM supplied)
 # Only zero values updated here so information can refresh nodes that lost information
 if(is.matrix(SEAM) == T)
  {
  RandSEAM[] <- rbinom(n=nrow(SDDprob)^2, size=1, prob = SEAM*Detected)
  InfoTransferred = ifelse(colSums(RandSEAM)>0,1,0)
  HaveInfo[HaveInfo == 0] = InfoTransferred[HaveInfo == 0]
  }
 
 
 # Add invasion resulting from colonisation from external sources.
 # The existing host draw and capacity clipping are unchanged.
 NBeforeExternalInvasion <- N
 if(OngoingExternalInvasion == T)
  {
  if(is.matrix(InvasionRisk) == F)
   ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
  if(is.matrix(InvasionRisk) == T)
   ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,timestep])
  Invaded[Invaded == 0] = ExternalInvasion[Invaded==0]
  if(is.na(IncursionStartPop) == T) 
	N = N+ExternalInvasion
  if(is.na(IncursionStartPop) == F) 
	N = N+ExternalInvasion*IncursionStartPop
  N[N > NodeK] = NodeK[N > NodeK] 
  }

  # Synchronise accepted external host immigrants with pathogen state. NULL
  # retains the historical assumption that all ordinary host gains are S.
  # A named state distribution assigns only the actually accepted increase
  # after the unchanged host carrying-capacity calculation.
  if(UsePathogen)
    {
    if(is.null(ExternalPathogenStateProbResolved)) {
      PathogenState <- PathogenEngine$Reconcile(PathogenState, N, PathogenContext)
    } else {
      ExternalBaseN <- pmin(as.integer(NBeforeExternalInvasion), as.integer(N))
      PathogenState <- PathogenEngine$Reconcile(PathogenState, ExternalBaseN, PathogenContext)
      ExternalAccepted <- pmax(0L, as.integer(N) - as.integer(NBeforeExternalInvasion))
      for(ii in which(ExternalAccepted > 0L)) {
        ExternalByState <- as.integer(rmultinom(1L, size = ExternalAccepted[ii], prob = ExternalPathogenStateProbResolved))
        PathogenState[ii, PathogenEngine$States] <-
          PathogenState[ii, PathogenEngine$States] + ExternalByState
      }
      storage.mode(PathogenState) <- "integer"
      if(any(rowSums(PathogenState) != as.integer(N)))
        stop("External host pathogen-state assignment violated S/E/I/R = N")
    }

    # Pathogen transmission/progression/recovery occurs once after host events.
    PathogenStep <- PathogenEngine$Step(PathogenState, N, timestep, PathogenContext)
    PathogenState <- PathogenStep$State
    N <- PathogenStep$N
    PathogenStateResultsLoop[,,timestep] <- PathogenState
    PathogenDetectionP <- PathogenEngine$Resolve(PathogenOriginal$DetectionProb, timestep, PathogenContext, "DetectionProb")
    PathogenDetectedNow <- rbinom(nrow(SDDprob), size = 1, prob = 1-(1-PathogenDetectionP)^PathogenState[,"I"])
    PathogenDetectedResultsLoop[,timestep] <- PathogenDetectedNow
    if(PathogenInformationAcquisition)
      {
      if(UseInfoPersistence == T) LastKnownPresence[PathogenDetectedNow == 1] <- timestep
      HaveInfo[HaveInfo == 0] <- PathogenDetectedNow[HaveInfo == 0]
      }
    }
 
  # Add nodes with information resulting from external sources
  if(OngoingExternalInfo == T)
    {
    if(is.matrix(ExternalInfoProb) == F)
      ExternalInfo = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb)
    if(is.matrix(ExternalInfoProb) == T)
      ExternalInfo = rbinom(1:nrow(SDDprob),size = 1,prob = ExternalInfoProb[,timestep])
    HaveInfo[HaveInfo == 0] = ExternalInfo[HaveInfo==0]
     }
 # Update host/pest presence from current abundance.
 Invaded = ifelse(N>0,1,0)
 
 # Record nodes adopting management
 ManagingResultsLoop[,timestep] = Managing
  
 # Record infested nodes
 InvasionResultsLoop[,timestep] = Invaded

 # Record populations
 PopulationResultsLoop[,timestep] = N

 # Two host-surveillance streams. Freeze the information state before either
 # stream so a background detection cannot activate targeted surveillance in
 # the same round.
 InfoBeforeSurveillance = as.integer(HaveInfo != 0)
 InformationStateBeforeSurveillanceResultsLoop[,timestep] = InfoBeforeSurveillance
 BackgroundDetectionProbabilityResultsLoop[,timestep] = NodeDetectionProb
 BackgroundDetection = rbinom(1:nrow(SDDprob),size = 1,prob = 1-(1-NodeDetectionProb)^N)
 InfoTriggeredDetection = integer(nrow(SDDprob))
 if(UseInfoTriggeredSurveillance)
   {
   InfoTriggeredDetectionProbabilityResultsLoop[,timestep] = NodeInfoTriggeredDetectionProb
   InfoTriggeredDetection = rbinom(1:nrow(SDDprob),size = 1,
     prob = (1-(1-NodeInfoTriggeredDetectionProb)^N) * InfoBeforeSurveillance)
   }
 BackgroundDetectedResultsLoop[,timestep] = BackgroundDetection
 InfoTriggeredDetectedResultsLoop[,timestep] = InfoTriggeredDetection
 HostDetectionEvidence = pmax(BackgroundDetection,InfoTriggeredDetection)
 
 # Record direct host evidence as known local presence only when host evidence
 # is an enabled information-acquisition route.
 if(HostInformationAcquisition && UseInfoPersistence == T)
   {
   KnownPresence = which(HostDetectionEvidence == 1)
   if(length(KnownPresence) > 0)
     LastKnownPresence[KnownPresence] = timestep
   }
 
 # Add new host evidence to info only when host evidence is enabled.
 if(HostInformationAcquisition)
   HaveInfo[HaveInfo==0] = HostDetectionEvidence[HaveInfo==0]
 HaveInfoResultsLoop[,timestep] = HaveInfo
 
 # Legacy DetectedResults remains the persistent known-present state.
 DetectedResultsLoop[,timestep] = HaveInfo*Invaded 
 }
 Result <- list(Invasion = InvasionResultsLoop, Population = PopulationResultsLoop,
                Managing = ManagingResultsLoop, Detected = DetectedResultsLoop,
                BackgroundDetected = BackgroundDetectedResultsLoop,
                InfoTriggeredDetected = InfoTriggeredDetectedResultsLoop,
                BackgroundDetectionProbability = BackgroundDetectionProbabilityResultsLoop,
                InfoTriggeredDetectionProbability = InfoTriggeredDetectionProbabilityResultsLoop,
                InformationStateBeforeSurveillance = InformationStateBeforeSurveillanceResultsLoop,
                HaveInfo = HaveInfoResultsLoop)
 if(UsePathogen)
   {
   Result$PathogenState <- PathogenStateResultsLoop
   Result$PathogenDetected <- PathogenDetectedResultsLoop
   }
 Result
}

# Run stochastic realisations serially in permutation order.
# This is deliberately the same worker body used by INApestMetaParallel,
# so model logic remains aligned between serial and parallel implementations.

# ---------------------------------------------------------------------------
# Run all independent stochastic permutations.
# ---------------------------------------------------------------------------
PermutationResults <- lapply(seq_len(Nperm), PermutationWorker)

# ---------------------------------------------------------------------------
# Combine permutation histories into standard result arrays.
# ---------------------------------------------------------------------------
BindMeta3 <- function(field)
  {
  prototype <- PermutationResults[[1L]][[field]]
  out <- array(vector(typeof(prototype), nrow(SDDprob) * Ntimesteps * Nperm),
               dim = c(nrow(SDDprob), Ntimesteps, Nperm))
  for(pp in seq_len(Nperm)) out[,,pp] <- PermutationResults[[pp]][[field]]
  out
  }
InvasionResults <- BindMeta3("Invasion")
PopulationResults <- BindMeta3("Population")
ManagingResults <- BindMeta3("Managing")
DetectedResults <- BindMeta3("Detected")
BackgroundDetectedResults <- BindMeta3("BackgroundDetected")
InfoTriggeredDetectedResults <- BindMeta3("InfoTriggeredDetected")
InformationStateBeforeSurveillanceResults <- BindMeta3("InformationStateBeforeSurveillance")
HaveInfoResults <- BindMeta3("HaveInfo")
BackgroundDetectionProbabilityResults <- BindMeta3("BackgroundDetectionProbability")
InfoTriggeredDetectionProbabilityResults <- BindMeta3("InfoTriggeredDetectionProbability")
if(UsePathogen)
  {
  PathogenStates <- PathogenEngine$States
  PathogenStateResults <- array(0L, dim = c(nrow(SDDprob), length(PathogenStates), Ntimesteps, Nperm),
                                dimnames = list(NULL, PathogenStates, NULL, NULL))
  PathogenDetectedResults <- array(0L, dim = c(nrow(SDDprob), Ntimesteps, Nperm))
  for(pp in seq_len(Nperm))
    {
    PathogenStateResults[,,,pp] <- PermutationResults[[pp]]$PathogenState
    PathogenDetectedResults[,,pp] <- PermutationResults[[pp]]$PathogenDetected
    }
  }
###########################################################
### End of Simulation
###########################################################

###########################################################
### Save results for post-hoc analyses
###########################################################
### ModelName used to generate filenames
# Use standard format for ease of reading results to produce heat maps
# Support post-processing comparisons among management scenarios.
if(is.na(OutputDir) == T)
	OutputDir = ""
FileNameStem = paste0(OutputDir,ModelName)

# These are 3D arrays with dimensions (Nodes,Timesteps,Realisations)
saveRDS(ManagingResults, paste0(FileNameStem,"InfoLargeOut.rds"))
saveRDS(PopulationResults, paste0(FileNameStem,"PopulationLargeOut.rds"))
saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))
saveRDS(BackgroundDetectedResults, paste0(FileNameStem,"BackgroundDetectedLargeOut.rds"))
saveRDS(InfoTriggeredDetectedResults, paste0(FileNameStem,"InfoTriggeredDetectedLargeOut.rds"))
saveRDS(InformationStateBeforeSurveillanceResults, paste0(FileNameStem,"InformationStateBeforeSurveillanceLargeOut.rds"))
saveRDS(HaveInfoResults, paste0(FileNameStem,"HaveInfoLargeOut.rds"))
saveRDS(BackgroundDetectionProbabilityResults, paste0(FileNameStem,"BackgroundDetectionProbabilityLargeOut.rds"))
saveRDS(InfoTriggeredDetectionProbabilityResults, paste0(FileNameStem,"InfoTriggeredDetectionProbabilityLargeOut.rds"))
if(UsePathogen) {
  saveRDS(PathogenStateResults, paste0(FileNameStem,"PathogenStateLargeOut.rds"))
  saveRDS(PathogenDetectedResults, paste0(FileNameStem,"PathogenDetectedLargeOut.rds"))
}

##########################################################
### Store node-level invasion probabilities for each timestep.
### and estimation of invasion threat to other regions
##########################################################

InvasionProb = matrix(ncol = Ntimesteps, nrow = nrow(SDDprob))
for(timestep in 1:Ntimesteps)
{
TimestepData = matrix(InvasionResults[,timestep,,drop=FALSE],nrow=nrow(SDDprob),ncol=Nperm)
InvasionProb[,timestep] = rowSums(TimestepData)/Nperm
}
saveRDS(InvasionProb, paste0(FileNameStem,"InvasionProb.rds"))
if(DoPlots == T)
{
###########################################################
### Produce summary figs when processing completed
###########################################################

Title = ModelName


# Change in total population with time
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

PopulationSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(PopulationSummary) = c("Realisation",   "Timestep",  "NodesInfested")

if(is.matrix(K) == TRUE)
    inv_K <- 1 / colSums(K)

for(perm in 1:Nperm)
{
PopulationData = PopulationResults[,,perm]
dim(PopulationData)
NodesInfested = colSums(PopulationData)*inv_K
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

# Change in number of nodes infested with time
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

InvasionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(InvasionSummary) = c("Realisation",   "Timestep",  "NodesInfested")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,perm]
dim(InvasionData)
NodesInfested = colSums(InvasionData)
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


# Summarise the number of nodes under management through time.
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

ManagingSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(ManagingSummary) = c("Realisation",   "Timestep",  "NodesManaging")

for(perm in 1:Nperm)
{
ManagingData = ManagingResults[,,perm]
dim(ManagingData)
NodesManaging = colSums(ManagingData)
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


# Summarise occupied nodes known to the response system through time.
# Plots of raw values for each realisation and summaries (median and 95% CI) provided

DetectedSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedSummary) = c("Realisation",   "Timestep",  "NodesDetected")

for(perm in 1:Nperm)
{
DetectedData = DetectedResults[,,perm]
dim(DetectedData)
NodesDetected = colSums(DetectedData)
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


# Summarise the proportion of occupied nodes known through time.
# Plots of raw values for each realisation and summaries (median and 95% CI) provided
 
DetectedProportionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(DetectedProportionSummary) = c("Realisation",   "Timestep",  "DetectedProportion")

for(perm in 1:Nperm)
{
InvasionData = InvasionResults[,,perm]
DetectedData = DetectedResults[,,perm]
NodesDetected = colSums(DetectedData)
NodesInvaded = colSums(InvasionData)
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
if(ReturnResults)
  {
  ResultObject <- list(ModelName=ModelName, InvasionResults=InvasionResults, PopulationResults=PopulationResults,
    ManagingResults=ManagingResults, DetectedResults=DetectedResults,
    BackgroundDetectedResults=BackgroundDetectedResults, InfoTriggeredDetectedResults=InfoTriggeredDetectedResults,
    BackgroundDetectionProbabilityResults=BackgroundDetectionProbabilityResults, InfoTriggeredDetectionProbabilityResults=InfoTriggeredDetectionProbabilityResults,
    InformationStateBeforeSurveillanceResults=InformationStateBeforeSurveillanceResults, HaveInfoResults=HaveInfoResults,
    InvasionProb=InvasionProb)
  if(UsePathogen)
    {
    ResultObject$PathogenStateResults <- PathogenStateResults
    ResultObject$PathogenDetectedResults <- PathogenDetectedResults
    }
  class(ResultObject) <- c("INApestMeta","list")
  return(invisible(ResultObject))
  }
invisible(NULL)
}


################################################################
################################################################
### End of function
################################################################
################################################################
