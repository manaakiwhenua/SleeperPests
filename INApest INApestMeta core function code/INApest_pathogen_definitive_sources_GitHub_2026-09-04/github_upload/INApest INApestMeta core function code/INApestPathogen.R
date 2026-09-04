###############################################################################
### INApestPathogen -- pathogen-state specification for INApest
###
### Defines the generic pathogen process coupled to INApest host simulations.
### Binary mode tracks pathogen presence/absence in occupied host nodes. SIS,
### SIR and SEIR modes track susceptible, exposed, infectious and recovered host
### counts as required by the selected model.
###
### Host abundance and pathogen state remain coherent: pathogen compartments are
### reconciled to the host population, pathogen mortality can remove hosts, and
### pathogen detection can optionally update the wider INApest information state.
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


