###############################################################################
### INApestPointPathogen -- pathogen interaction for explicit point models
###
### Attaches a pathogen state to each explicit host point. Candidate infectious-
### susceptible pairs are defined by realised point geometry and ContactRadius,
### then filtered by ContactProb and an optional distance/contact kernel.
###
### Infection, progression, recovery, pathogen mortality, introduction and
### detection are updated without creating new hosts. Host deaths are returned to
### the parent point engine so point identifiers remain coherent during updates.
###############################################################################

INApestPointPathogenInteraction <- function(
    Pathogen,                                   # Optional pathogen process specification
    ContactRadius = Inf,                        # Maximum distance for candidate infectious contacts
    ContactKernel = NULL,                       # Optional distance/contact multiplier for candidate pairs
    ContactProb = 1,                            # Probability a candidate point contact is realised
    InitialStateField = "pathogen_state",       # Point column storing pathogen state
    DefaultState = "S") {                       # Default pathogen state for points without an explicit state

  if (!inherits(Pathogen, "INApestPathogen"))
    stop("Pathogen must be created by INApestPathogen().")
  states <- Pathogen$States
  if (!(DefaultState %in% states))
    stop("DefaultState must be one of the states in Pathogen$States.")
  if (!is.numeric(ContactRadius) || length(ContactRadius) != 1L ||
      is.na(ContactRadius) || ContactRadius < 0)
    stop("ContactRadius must be a non-negative scalar or Inf.")
  if (!is.null(ContactKernel) && !is.function(ContactKernel))
    stop("ContactKernel must be NULL or a function(distance, source, target, timestep, perm).")

  clip01 <- function(x) pmin(1, pmax(0, x))

  resolve_point <- function(x, points, timestep, perm, name) {
    n <- nrow(points)
    if (!n) return(numeric(0))
    if (is.function(x)) {
      fm <- names(formals(x))
      args <- list(points = points, timestep = timestep, perm = perm)
      if (!is.null(fm) && !("..." %in% fm))
        args <- args[intersect(names(args), fm)]
      out <- do.call(x, args)
    } else if (length(x) == 1L) {
      out <- rep(as.numeric(x), n)
    } else {
      out <- as.numeric(x)
      if (length(out) != n)
        stop(name, " must be scalar, a function, or resolve to one value per point.")
    }
    if (length(out) == 1L) out <- rep(out, n)
    if (length(out) != n || any(!is.finite(out)))
      stop(name, " must resolve to one finite value per point.")
    out
  }

  ensure_state <- function(points, timestep, perm) {
    if (!InitialStateField %in% names(points))
      points[[InitialStateField]] <- DefaultState
    z <- as.character(points[[InitialStateField]])
    z[is.na(z) | !nzchar(z)] <- DefaultState
    bad <- setdiff(unique(z), states)
    if (length(bad))
      stop("Unknown point pathogen state(s): ", paste(bad, collapse = ", "))
    points[[InitialStateField]] <- z
    points
  }

  initialize_fun <- function(points, perm = 1L) {
    explicit <- InitialStateField %in% names(points) &&
      any(!is.na(points[[InitialStateField]]) & nzchar(as.character(points[[InitialStateField]])))
    points <- ensure_state(points, 0L, perm)
    if (explicit || !nrow(points)) return(points)

    get_count <- function(x, name) {
      if (is.function(x) || length(x) != 1L || !is.finite(x) || x < 0 || x != floor(x))
        stop(name, " must be a single non-negative whole-number count for point-model automatic seeding; alternatively provide explicit ", InitialStateField, " in InitialPoints.")
      as.integer(x)
    }
    ni <- get_count(Pathogen$InitialInfected, "InitialInfected")
    ne <- if ("E" %in% states) get_count(Pathogen$InitialExposed, "InitialExposed") else 0L
    nr <- if ("R" %in% states) get_count(Pathogen$InitialRecovered, "InitialRecovered") else 0L
    if (ni + ne + nr > nrow(points))
      stop("Initial point pathogen-state counts exceed the number of InitialPoints.")
    available <- seq_len(nrow(points))
    assign_state <- function(k, state) {
      if (!k) return(invisible(NULL))
      pick <- sample(available, k, replace = FALSE)
      points[[InitialStateField]][pick] <<- state
      available <<- setdiff(available, pick)
      invisible(NULL)
    }
    assign_state(ne, "E")
    assign_state(ni, "I")
    assign_state(nr, "R")
    points
  }

  contact_fun <- function(points, home_range = NULL, timestep, perm, context = list()) {
    points <- ensure_state(points, timestep, perm)
    n <- nrow(points)
    if (n < 2L)
      return(data.frame(id1 = integer(0), id2 = integer(0), distance = numeric(0)))

    infectious <- which(points[[InitialStateField]] == "I")
    susceptible <- which(points[[InitialStateField]] == "S")
    if (!length(infectious) || !length(susceptible))
      return(data.frame(id1 = integer(0), id2 = integer(0), distance = numeric(0)))

    pairs <- expand.grid(i = infectious, j = susceptible, KEEP.OUT.ATTRS = FALSE)
    dx <- points$x[pairs$i] - points$x[pairs$j]
    dy <- points$y[pairs$i] - points$y[pairs$j]
    d <- sqrt(dx^2 + dy^2)
    keep <- d <= ContactRadius
    pairs <- pairs[keep, , drop = FALSE]
    d <- d[keep]
    if (!nrow(pairs))
      return(data.frame(id1 = integer(0), id2 = integer(0), distance = numeric(0)))

    cp <- if (is.function(ContactProb)) {
      fm <- names(formals(ContactProb))
      args <- list(distance = d,
                   source = points[pairs$i, , drop = FALSE],
                   target = points[pairs$j, , drop = FALSE],
                   timestep = timestep, perm = perm)
      if (!is.null(fm) && !("..." %in% fm)) args <- args[intersect(names(args), fm)]
      do.call(ContactProb, args)
    } else ContactProb
    if (length(cp) == 1L) cp <- rep(cp, nrow(pairs))
    if (length(cp) != nrow(pairs) || any(!is.finite(cp)))
      stop("ContactProb must resolve to one finite probability per candidate pair.")

    if (!is.null(ContactKernel)) {
      fm <- names(formals(ContactKernel))
      args <- list(distance = d,
                   source = points[pairs$i, , drop = FALSE],
                   target = points[pairs$j, , drop = FALSE],
                   timestep = timestep, perm = perm)
      if (!is.null(fm) && !("..." %in% fm)) args <- args[intersect(names(args), fm)]
      km <- do.call(ContactKernel, args)
      if (length(km) == 1L) km <- rep(km, nrow(pairs))
      if (length(km) != nrow(pairs) || any(!is.finite(km)) || any(km < 0))
        stop("ContactKernel must return one finite non-negative multiplier per candidate pair.")
      cp <- cp * km
    }
    cp <- clip01(cp)
    realised <- stats::rbinom(length(cp), 1L, cp) == 1L
    if (!any(realised))
      return(data.frame(id1 = integer(0), id2 = integer(0), distance = numeric(0)))

    data.frame(id1 = points$id[pairs$i[realised]],
               id2 = points$id[pairs$j[realised]],
               distance = d[realised],
               stringsAsFactors = FALSE)
  }

  detect_fun <- function(points, timestep, perm) {
    points <- ensure_state(points, timestep, perm)
    n <- nrow(points)
    if (!n) return(logical(0))
    p <- clip01(resolve_point(Pathogen$DetectionProb, points, timestep, perm, "DetectionProb"))
    points[[InitialStateField]] == "I" & (stats::rbinom(n, 1L, p) == 1L)
  }

  update_fun <- function(points, contacts = NULL, home_range = NULL,
                         timestep, perm, context = list()) {
    points <- ensure_state(points, timestep, perm)
    n <- nrow(points)
    if (!n) return(list(points = points, events = data.frame()))

    state0 <- as.character(points[[InitialStateField]])
    beta <- resolve_point(Pathogen$Beta, points, timestep, perm, "Beta")
    rec <- resolve_point(Pathogen$RecoveryProb, points, timestep, perm, "RecoveryProb")
    mort <- resolve_point(Pathogen$PathogenMortalityProb, points, timestep, perm, "PathogenMortalityProb")
    prog <- resolve_point(Pathogen$ProgressionProb, points, timestep, perm, "ProgressionProb")
    waning <- resolve_point(Pathogen$ImmunityLossProb, points, timestep, perm, "ImmunityLossProb")
    intro_p <- resolve_point(Pathogen$IntroductionProb, points, timestep, perm, "IntroductionProb")
    if (any(rec < 0 | rec > 1) || any(mort < 0 | mort > 1) || any(prog < 0 | prog > 1) ||
        any(waning < 0 | waning > 1) || any(intro_p < 0 | intro_p > 1))
      stop("Point pathogen probabilities must resolve to [0,1].")
    if (any(rec + mort > 1 + 1e-12))
      stop("RecoveryProb + PathogenMortalityProb must not exceed 1.")
    if (any(beta < 0)) stop("Beta must be non-negative.")

    # Transmission uses realised I->S contacts. Multiple infectious contacts
    # compound independently for the target during this timestep.
    new_infection <- rep(FALSE, n)
    if (!is.null(contacts) && nrow(contacts)) {
      target_ids <- unique(contacts$id2)
      for (tid in target_ids) {
        j <- match(tid, points$id)
        if (is.na(j) || state0[j] != "S") next
        m <- sum(contacts$id2 == tid)
        if (Pathogen$Transmission == "frequency") {
          p <- clip01(1 - (1 - clip01(beta[j]))^m)
        } else {
          ds <- resolve_point(Pathogen$DensityScale, points, timestep, perm, "DensityScale")[j]
          if (!is.finite(ds) || ds <= 0) stop("DensityScale must be positive.")
          p <- clip01(-expm1(-beta[j] * m / ds))
        }
        new_infection[j] <- stats::rbinom(1L, 1L, p) == 1L
      }
    }

    # Repeated/pathway introductions convert susceptible points only; no hosts
    # are created by the pathogen component.
    intro <- state0 == "S" & !new_infection & (stats::rbinom(n, 1L, intro_p) == 1L)

    new_state <- state0
    if ("E" %in% states) {
      new_state[new_infection | intro] <- "E"
      eidx <- which(state0 == "E")
      if (length(eidx)) {
        move <- stats::rbinom(length(eidx), 1L, prog[eidx]) == 1L
        new_state[eidx[move]] <- "I"
      }
    } else {
      new_state[new_infection | intro] <- "I"
    }

    # Competing I exits. Deaths cannot be applied inside Interaction because the
    # point-core contract requires ids to be preserved. We therefore mark them
    # for the parent adapter, which removes them immediately after Interaction.
    pathogen_death <- rep(FALSE, n)
    iidx <- which(state0 == "I")
    if (length(iidx)) {
      for (j in iidx) {
        z <- sample.int(3L, 1L, prob = c(1 - rec[j] - mort[j], rec[j], mort[j]))
        if (z == 2L) new_state[j] <- if ("R" %in% states) "R" else "S"
        if (z == 3L) pathogen_death[j] <- TRUE
      }
    }

    if ("R" %in% states) {
      ridx <- which(state0 == "R")
      if (length(ridx)) {
        lose <- stats::rbinom(length(ridx), 1L, waning[ridx]) == 1L
        new_state[ridx[lose]] <- "S"
      }
    }

    points[[InitialStateField]] <- new_state
    points$.pathogen_death <- pathogen_death
    pathogen_detected <- detect_fun(points, timestep, perm) & !pathogen_death
    if (isTRUE(Pathogen$DetectionTriggersInfo) && any(pathogen_detected)) {
      if ("have_info" %in% names(points)) points$have_info[pathogen_detected] <- TRUE
      if ("last_known_timestep" %in% names(points)) points$last_known_timestep[pathogen_detected] <- timestep
    }

    changed <- which(new_state != state0 | pathogen_death | intro | new_infection | pathogen_detected)
    events <- if (length(changed)) data.frame(
      point_id = points$id[changed],
      from_state = state0[changed],
      to_state = new_state[changed],
      new_infection = new_infection[changed],
      introduced = intro[changed],
      pathogen_death = pathogen_death[changed],
      pathogen_detected = pathogen_detected[changed],
      stringsAsFactors = FALSE
    ) else data.frame()

    list(points = points, events = events)
  }

  structure(list(Initialize = initialize_fun,
                 Contact = contact_fun,
                 Update = update_fun,
                 Detect = detect_fun,
                 Pathogen = Pathogen,
                 PathogenStateField = InitialStateField),
            class = c("INApestPointPathogenInteraction", "list"))
}

INApestPointPathogenApplyDeaths <- function(points) {
  if (!is.data.frame(points)) stop("points must be a data.frame.")
  if (!".pathogen_death" %in% names(points)) return(points)
  keep <- is.na(points$.pathogen_death) | !points$.pathogen_death
  out <- points[keep, , drop = FALSE]
  out$.pathogen_death <- NULL
  rownames(out) <- NULL
  out
}

INApestPointPathogenOutputs <- function(PointHistory, state_field = "pathogen_state") {
  if (!is.data.frame(PointHistory)) stop("PointHistory must be a data.frame.")
  if (!(state_field %in% names(PointHistory)))
    stop("PointHistory does not contain ", state_field, ".")
  if (!all(c("perm", "timestep") %in% names(PointHistory)))
    stop("PointHistory must contain perm and timestep.")
  tab <- aggregate(
    rep(1L, nrow(PointHistory)),
    by = list(perm = PointHistory$perm,
              timestep = PointHistory$timestep,
              pathogen_state = PointHistory[[state_field]]),
    FUN = sum
  )
  names(tab)[4] <- "n"
  tab
}
