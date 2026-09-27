###############################################################################
### INApestBiocontrolPointAdapter.R
### Point-pest / node-biocontrol adapter for INApest
### Candidate v0.1 - 19 September 2026
###
### Design contract
### - Pest state remains explicit points with x/y coordinates and stable IDs.
### - Biocontrol agents remain node x agent-stage populations.
### - Existing INApestBiocontrolAgent()/INApestBiocontrol() definitions are
###   consumed; this file adds only point support and point-ID resolution.
### - Event order matches the validated non-point companion:
###     release -> attack -> transition -> recruitment -> movement.
### - Multiple agents contribute simultaneous hazards; one pest point can be
###   removed at most once in a timestep.
###############################################################################

.ibp_stop <- function(...) stop(..., call. = FALSE)

.ibp_field <- function(x, names, required = FALSE, default = NULL) {
  if (is.null(x)) {
    if (required) .ibp_stop("Missing required object while resolving: ", names[1L])
    return(default)
  }
  nms <- names(x)
  if (!is.null(nms)) {
    hit <- match(tolower(names), tolower(nms), nomatch = 0L)
    hit <- hit[hit > 0L]
    if (length(hit)) return(x[[hit[1L]]])
  }
  if (required) .ibp_stop("Required field not found: ", paste(names, collapse = "/"))
  default
}

.ibp_call <- function(fun, args) {
  if (!is.function(fun)) .ibp_stop("Expected a function.")
  fml <- names(formals(fun))
  if (is.null(fml) || "..." %in% fml) return(do.call(fun, args))
  do.call(fun, args[intersect(names(args), fml)])
}

.ibp_agents <- function(Biocontrol) {
  z <- .ibp_field(Biocontrol, c("Agents", "agents"), default = NULL)
  if (is.null(z)) {
    eng <- .ibp_field(Biocontrol, c("Engine", "engine"), default = NULL)
    z <- .ibp_field(eng, c("Agents", "agents"), default = NULL)
  }
  if (is.null(z) && is.list(Biocontrol) &&
      !is.null(.ibp_field(Biocontrol, c("Stages", "stages"), default = NULL))) {
    z <- list(Biocontrol)
  }
  if (is.null(z)) .ibp_stop(
    "Biocontrol does not expose an Agents list. Supply an object returned by ",
    "the definitive INApestBiocontrol() companion."
  )
  if (!is.list(z)) z <- list(z)
  if (!length(z)) .ibp_stop("Biocontrol Agents list is empty.")
  z
}

.ibp_dt <- function(Biocontrol) {
  z <- .ibp_field(Biocontrol, c("TimestepLength", "timestep_length", "dt"), default = 1)
  if (!is.numeric(z) || length(z) != 1L || !is.finite(z) || z <= 0)
    .ibp_stop("Biocontrol timestep length must be one finite positive value.")
  as.numeric(z)
}

INApestBiocontrolPointSupport <- function(
  xmin, xmax, ymin, ymax, nrow, ncol,
  row_origin = c("ymax", "ymin"),
  outside = c("exclude", "error"),
  crs = NA_character_
) {
  row_origin <- match.arg(row_origin)
  outside <- match.arg(outside)
  vals <- c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  if (any(!is.finite(vals)) || xmax <= xmin || ymax <= ymin)
    .ibp_stop("Point-support extent must be finite with xmax > xmin and ymax > ymin.")
  if (!is.numeric(nrow) || length(nrow) != 1L || nrow < 1L || nrow != floor(nrow) ||
      !is.numeric(ncol) || length(ncol) != 1L || ncol < 1L || ncol != floor(ncol))
    .ibp_stop("nrow and ncol must be positive integers.")
  nrow <- as.integer(nrow); ncol <- as.integer(ncol)
  structure(list(
    xmin = as.numeric(xmin), xmax = as.numeric(xmax),
    ymin = as.numeric(ymin), ymax = as.numeric(ymax),
    nrow = nrow, ncol = ncol,
    xres = (xmax - xmin) / ncol,
    yres = (ymax - ymin) / nrow,
    n_nodes = nrow * ncol,
    row_origin = row_origin, outside = outside,
    crs = as.character(crs)[1L]
  ), class = "INApestBiocontrolPointSupport")
}

print.INApestBiocontrolPointSupport <- function(x, ...) {
  cat("INApestBiocontrolPointSupport\n")
  cat("  nodes:", x$n_nodes, "(", x$nrow, "x", x$ncol, ")\n")
  cat("  extent:", x$xmin, x$xmax, x$ymin, x$ymax, "\n")
  cat("  row_origin:", x$row_origin, "\n")
  invisible(x)
}

INApestBiocontrolPointNode <- function(x, y, PointSupport) {
  if (!inherits(PointSupport, "INApestBiocontrolPointSupport"))
    .ibp_stop("PointSupport must be created by INApestBiocontrolPointSupport().")
  if (length(x) != length(y)) .ibp_stop("x and y must have equal length.")
  x <- as.numeric(x); y <- as.numeric(y)
  inside <- is.finite(x) & is.finite(y) &
    x >= PointSupport$xmin & x <= PointSupport$xmax &
    y >= PointSupport$ymin & y <= PointSupport$ymax
  if (PointSupport$outside == "error" && any(!inside))
    .ibp_stop("One or more pest points fall outside BiocontrolPointSupport.")
  out <- rep(NA_integer_, length(x))
  if (!any(inside)) return(out)
  xx <- x[inside]; yy <- y[inside]
  col <- as.integer(floor((xx - PointSupport$xmin) / PointSupport$xres)) + 1L
  col[xx == PointSupport$xmax] <- PointSupport$ncol
  if (PointSupport$row_origin == "ymax") {
    row <- as.integer(floor((PointSupport$ymax - yy) / PointSupport$yres)) + 1L
    row[yy == PointSupport$ymin] <- PointSupport$nrow
  } else {
    row <- as.integer(floor((yy - PointSupport$ymin) / PointSupport$yres)) + 1L
    row[yy == PointSupport$ymax] <- PointSupport$nrow
  }
  good <- row >= 1L & row <= PointSupport$nrow & col >= 1L & col <= PointSupport$ncol
  node <- rep(NA_integer_, length(row))
  node[good] <- as.integer((row[good] - 1L) * PointSupport$ncol + col[good])
  out[inside] <- node
  out
}

.ibp_as_matrix <- function(x, nr, nc, name) {
  if (is.data.frame(x)) x <- as.matrix(x)
  if (is.matrix(x)) {
    if (!all(dim(x) == c(nr, nc)))
      .ibp_stop(name, " must have dimensions n_nodes x n_agent_stages (", nr, " x ", nc, ").")
    z <- matrix(as.numeric(x), nr, nc, dimnames = dimnames(x))
  } else if (is.numeric(x) && length(x) == nr * nc) {
    z <- matrix(as.numeric(x), nr, nc)
  } else if (is.numeric(x) && nr == 1L && length(x) == nc) {
    z <- matrix(as.numeric(x), 1L, nc)
  } else {
    .ibp_stop(name, " must be a numeric n_nodes x n_agent_stages matrix.")
  }
  if (any(!is.finite(z)) || any(z < 0)) .ibp_stop(name, " must contain finite non-negative values.")
  if (any(abs(z - round(z)) > 1e-8)) .ibp_stop(name, " must contain integer counts for stochastic point coupling.")
  matrix(as.integer(round(z)), nr, nc, dimnames = dimnames(z))
}

.ibp_agent_name <- function(agent, i) {
  z <- .ibp_field(agent, c("Name", "name"), default = paste0("agent", i))
  as.character(z)[1L]
}

.ibp_agent_stages <- function(agent) {
  z <- .ibp_field(agent, c("Stages", "stages"), required = TRUE)
  z <- as.character(z)
  if (!length(z) || anyNA(z) || any(!nzchar(z)) || anyDuplicated(z))
    .ibp_stop("Every biocontrol agent must have unique, non-empty stage names.")
  z
}

.ibp_resolve_initial <- function(agent, context, perm) {
  st <- .ibp_agent_stages(agent)
  z <- .ibp_field(agent, c("InitialState", "initial_state", "Initial"), required = TRUE)
  if (is.function(z)) z <- .ibp_call(z, list(
    context = context, n_nodes = context$n_nodes, stages = st,
    perm = perm, Ntimesteps = context$Ntimesteps
  ))
  m <- .ibp_as_matrix(z, context$n_nodes, length(st), "Biocontrol InitialState")
  colnames(m) <- st
  m
}

.ibp_resolve_release <- function(agent, context, timestep, perm) {
  st <- .ibp_agent_stages(agent)
  z <- .ibp_field(agent, c("Release", "release"), default = NULL)
  if (is.null(z)) return(matrix(0L, context$n_nodes, length(st), dimnames = list(NULL, st)))
  if (is.function(z)) z <- .ibp_call(z, list(
    timestep = timestep, perm = perm, context = context,
    n_nodes = context$n_nodes, stages = st, Ntimesteps = context$Ntimesteps
  ))
  if (length(dim(z)) == 3L) {
    if (dim(z)[1L] != context$n_nodes || dim(z)[2L] != length(st) || timestep > dim(z)[3L])
      .ibp_stop("Biocontrol Release array must be n_nodes x n_agent_stages x time.")
    z <- z[, , timestep, drop = FALSE][, , 1L]
  }
  m <- .ibp_as_matrix(z, context$n_nodes, length(st), "Biocontrol Release")
  colnames(m) <- st
  m
}

.ibp_resolve_transition <- function(agent, timestep, perm, context) {
  st <- .ibp_agent_stages(agent); ns <- length(st)
  z <- .ibp_field(agent, c("Transition", "transition"), default = diag(ns))
  if (is.function(z)) z <- .ibp_call(z, list(
    timestep = timestep, perm = perm, context = context, stages = st
  ))
  if (length(dim(z)) == 3L) {
    if (dim(z)[1L] != ns || dim(z)[2L] != ns || timestep > dim(z)[3L])
      .ibp_stop("Biocontrol Transition array must be agent_stage x agent_stage x time.")
    z <- z[, , timestep]
  }
  if (!is.matrix(z) || !all(dim(z) == c(ns, ns)) || any(!is.finite(z)) || any(z < 0))
    .ibp_stop("Biocontrol Transition must be a finite non-negative square stage matrix.")
  if (any(colSums(z) > 1 + 1e-10))
    .ibp_stop("Biocontrol Transition column sums must not exceed 1 (columns are source stages).")
  z
}

.ibp_resolve_movement <- function(agent, timestep, perm, context) {
  nn <- context$n_nodes
  z <- .ibp_field(agent, c("Movement", "movement"), default = NULL)
  if (is.null(z)) return(diag(nn))
  if (is.function(z)) z <- .ibp_call(z, list(
    timestep = timestep, perm = perm, context = context, n_nodes = nn
  ))
  if (length(dim(z)) == 3L) {
    if (dim(z)[1L] != nn || dim(z)[2L] != nn || timestep > dim(z)[3L])
      .ibp_stop("Biocontrol Movement array must be node x node x time.")
    z <- z[, , timestep]
  }
  if (!is.matrix(z) || !all(dim(z) == c(nn, nn)) || any(!is.finite(z)) || any(z < 0))
    .ibp_stop("Biocontrol Movement must be a finite non-negative node x node matrix.")
  if (any(rowSums(z) > 1 + 1e-10))
    .ibp_stop("Biocontrol Movement row sums must not exceed 1 (rows are source nodes).")
  z
}

.ibp_resolve_rate <- function(agent, timestep, perm, context) {
  z <- .ibp_field(agent, c("AttackRate", "attack_rate"), required = TRUE)
  if (is.function(z)) z <- .ibp_call(z, list(
    timestep = timestep, perm = perm, context = context, n_nodes = context$n_nodes
  ))
  if (!is.numeric(z) || !(length(z) %in% c(1L, context$n_nodes)) || any(!is.finite(z)) || any(z < 0))
    .ibp_stop("Biocontrol AttackRate must be a finite non-negative scalar or n_nodes vector.")
  rep_len(as.numeric(z), context$n_nodes)
}

.ibp_attack_stage_index <- function(agent) {
  st <- .ibp_agent_stages(agent)
  z <- .ibp_field(agent, c("AttackStage", "attack_stage"), required = TRUE)
  if (is.numeric(z)) idx <- as.integer(z) else idx <- match(as.character(z), st)
  if (!length(idx) || anyNA(idx) || any(idx < 1L | idx > length(st)))
    .ibp_stop("AttackStage does not match the agent stage names.")
  unique(idx)
}

.ibp_target_mask <- function(points, agent, Architecture, HostStages) {
  z <- .ibp_field(agent, c("TargetStage", "target_stage"), default = NULL)
  if (is.null(z) || (is.character(z) && any(tolower(z) %in% c("all", "default", "any"))))
    return(rep(TRUE, nrow(points)))
  if (!"stage" %in% names(points))
    .ibp_stop("Stage-targeted biocontrol requires points$stage.")
  if (Architecture == "MetaPoint") {
    return(as.character(points$stage) %in% as.character(z))
  }
  ps <- points$stage
  if (is.numeric(z)) return(as.integer(ps) %in% as.integer(z))
  if (is.null(HostStages)) {
    # Permit character numerals when no stage-name map is supplied.
    zz <- suppressWarnings(as.integer(as.character(z)))
    if (anyNA(zz)) .ibp_stop("Character TargetStage in PointTransitionMatrix requires HostStages names.")
    return(as.integer(ps) %in% zz)
  }
  idx <- match(as.character(z), as.character(HostStages))
  if (anyNA(idx)) .ibp_stop("TargetStage does not match HostStages.")
  as.integer(ps) %in% idx
}

.ibp_offspring <- function(agent) {
  z <- .ibp_field(agent, c("OffspringPerAttack", "offspring_per_attack"), default = 1)
  if (!is.numeric(z) || length(z) != 1L || !is.finite(z) || z < 0 || abs(z - round(z)) > 1e-8)
    .ibp_stop("OffspringPerAttack must be one non-negative integer for point coupling.")
  as.integer(round(z))
}

.ibp_recruit_stage <- function(agent) {
  st <- .ibp_agent_stages(agent)
  z <- .ibp_field(agent, c("RecruitStage", "recruit_stage"), required = TRUE)
  idx <- if (is.numeric(z)) as.integer(z) else match(as.character(z)[1L], st)
  if (length(idx) != 1L || is.na(idx) || idx < 1L || idx > length(st))
    .ibp_stop("RecruitStage does not match the agent stage names.")
  idx
}

.ibp_multinom_counts <- function(n, probs) {
  if (n <= 0L) return(integer(length(probs)))
  probs <- pmax(0, as.numeric(probs))
  s <- sum(probs)
  if (s > 1 + 1e-10) .ibp_stop("Probability mass exceeds 1.")
  cprob <- c(probs, max(0, 1 - s))
  as.integer(stats::rmultinom(1L, size = as.integer(n), prob = cprob)[seq_along(probs), 1L])
}

.ibp_transition_state <- function(state, T) {
  nr <- nrow(state); ns <- ncol(state)
  out <- matrix(0L, nr, ns, dimnames = dimnames(state))
  for (node in seq_len(nr)) for (src in seq_len(ns)) {
    if (state[node, src] > 0L) out[node, ] <- out[node, ] + .ibp_multinom_counts(state[node, src], T[, src])
  }
  out
}

.ibp_move_state <- function(state, M) {
  nn <- nrow(state); ns <- ncol(state)
  out <- matrix(0L, nn, ns, dimnames = dimnames(state))
  for (stage in seq_len(ns)) for (src in seq_len(nn)) {
    if (state[src, stage] > 0L) out[, stage] <- out[, stage] + .ibp_multinom_counts(state[src, stage], M[src, ])
  }
  out
}

INApestPointBiocontrolInitial <- function(
  Biocontrol, PointSupport, Ntimesteps, perm = 1L,
  Architecture = c("MetaPoint", "PointTransitionMatrix"), HostStages = NULL
) {
  Architecture <- match.arg(Architecture)
  if (is.null(Biocontrol)) return(NULL)
  if (!inherits(PointSupport, "INApestBiocontrolPointSupport"))
    .ibp_stop("BiocontrolPointSupport is required when Biocontrol is active in a point engine.")
  agents <- .ibp_agents(Biocontrol)
  context <- list(
    Architecture = Architecture, architecture = Architecture,
    n_nodes = PointSupport$n_nodes, Ntimesteps = as.integer(Ntimesteps),
    host_stages = HostStages, PointSupport = PointSupport
  )
  states <- lapply(seq_along(agents), function(i) .ibp_resolve_initial(agents[[i]], context, perm))
  names(states) <- vapply(seq_along(agents), function(i) .ibp_agent_name(agents[[i]], i), character(1))
  structure(list(State = states, Context = context), class = "INApestPointBiocontrolState")
}

INApestPointBiocontrolStateHistory <- function(state, perm, timestep) {
  if (is.null(state)) return(data.frame())
  out <- list(); k <- 0L
  for (a in seq_along(state$State)) {
    z <- state$State[[a]]
    nm <- names(state$State)[a]
    for (s in seq_len(ncol(z))) {
      k <- k + 1L
      out[[k]] <- data.frame(
        perm = perm, timestep = timestep, agent = nm,
        node = seq_len(nrow(z)), stage = colnames(z)[s],
        abundance = as.integer(z[, s]), stringsAsFactors = FALSE
      )
    }
  }
  do.call(rbind, out)
}

INApestPointBiocontrolStep <- function(
  points, PointBiocontrolState, Biocontrol, PointSupport,
  timestep, perm = 1L,
  Architecture = c("MetaPoint", "PointTransitionMatrix"), HostStages = NULL
) {
  Architecture <- match.arg(Architecture)
  if (is.null(Biocontrol)) return(list(
    points = points, state = PointBiocontrolState,
    attacks = data.frame(), history = data.frame(), n_attacked = 0L
  ))
  if (!is.data.frame(points) || !all(c("id", "x", "y") %in% names(points)))
    .ibp_stop("Point biocontrol requires a points data.frame containing id, x and y.")
  if (is.null(PointBiocontrolState))
    .ibp_stop("PointBiocontrolState is NULL; call INApestPointBiocontrolInitial() first.")

  agents <- .ibp_agents(Biocontrol)
  context <- PointBiocontrolState$Context
  context$timestep <- timestep; context$perm <- perm
  dt <- .ibp_dt(Biocontrol)
  states <- PointBiocontrolState$State
  if (length(states) != length(agents)) .ibp_stop("Biocontrol state/agent count mismatch.")

  # 1. Scheduled release occurs before attack.
  for (a in seq_along(agents)) {
    rel <- .ibp_resolve_release(agents[[a]], context, timestep, perm)
    states[[a]] <- states[[a]] + rel
  }

  node <- INApestBiocontrolPointNode(points$x, points$y, PointSupport)
  np <- nrow(points); na <- length(agents)
  hazard <- matrix(0, np, na)
  agent_names <- vapply(seq_along(agents), function(i) .ibp_agent_name(agents[[i]], i), character(1))
  colnames(hazard) <- agent_names

  # Build simultaneous per-point hazards from every eligible agent.
  if (np) for (a in seq_along(agents)) {
    eligible <- .ibp_target_mask(points, agents[[a]], Architecture, HostStages) & !is.na(node)
    if (!any(eligible)) next
    atk_stage <- .ibp_attack_stage_index(agents[[a]])
    attackers <- rowSums(states[[a]][, atk_stage, drop = FALSE])
    rate <- .ibp_resolve_rate(agents[[a]], timestep, perm, context)
    ii <- which(eligible)
    hazard[ii, a] <- dt * rate[node[ii]] * attackers[node[ii]]
  }

  total_hazard <- if (np) rowSums(hazard) else numeric(0)
  attack_prob <- 1 - exp(-total_hazard)
  attacked <- if (np) stats::runif(np) < attack_prob else logical(0)
  responsible <- rep(NA_integer_, np)
  if (any(attacked)) {
    for (i in which(attacked)) {
      h <- hazard[i, ]
      if (sum(h) <= 0) next
      u <- stats::runif(1L) * sum(h)
      responsible[i] <- which(cumsum(h) >= u)[1L]
    }
  }
  attacked <- attacked & !is.na(responsible)

  attack_events <- data.frame()
  attack_count <- lapply(seq_along(agents), function(a) integer(PointSupport$n_nodes))
  if (any(attacked)) {
    idx <- which(attacked)
    ae <- points[idx, intersect(c("id", "parent_id", "x", "y", "stage"), names(points)), drop = FALSE]
    names(ae)[names(ae) == "id"] <- "point_id"
    if (!"parent_id" %in% names(ae)) ae$parent_id <- NA_integer_
    if (!"stage" %in% names(ae)) ae$stage <- NA
    ae$perm <- perm; ae$timestep <- timestep
    ae$node <- node[idx]
    ae$agent <- agent_names[responsible[idx]]
    ae$hazard <- total_hazard[idx]
    ae$attack_probability <- attack_prob[idx]
    ae$detail <- "biocontrol_mortality"
    attack_events <- ae[, c("perm", "timestep", "point_id", "parent_id", "x", "y", "stage", "node", "agent", "hazard", "attack_probability", "detail"), drop = FALSE]
    for (a in seq_along(agents)) {
      jj <- idx[responsible[idx] == a]
      if (length(jj)) attack_count[[a]] <- tabulate(node[jj], nbins = PointSupport$n_nodes)
    }
    points <- points[-idx, , drop = FALSE]
  }

  # 2 was attack. 3. Existing + released agents transition.
  for (a in seq_along(agents)) {
    T <- .ibp_resolve_transition(agents[[a]], timestep, perm, context)
    states[[a]] <- .ibp_transition_state(states[[a]], T)
  }

  # 4. Attack-derived recruits enter after transition.
  for (a in seq_along(agents)) {
    n_off <- .ibp_offspring(agents[[a]])
    if (n_off <= 0L || !any(attack_count[[a]])) next
    rs <- .ibp_recruit_stage(agents[[a]])
    states[[a]][, rs] <- states[[a]][, rs] + as.integer(attack_count[[a]] * n_off)
  }

  # 5. Agent movement occurs last.
  for (a in seq_along(agents)) {
    M <- .ibp_resolve_movement(agents[[a]], timestep, perm, context)
    states[[a]] <- .ibp_move_state(states[[a]], M)
  }

  PointBiocontrolState$State <- states
  hist <- INApestPointBiocontrolStateHistory(PointBiocontrolState, perm, timestep)
  list(
    points = points,
    state = PointBiocontrolState,
    attacks = attack_events,
    history = hist,
    n_attacked = as.integer(sum(attacked))
  )
}

###############################################################################
### End point adapter
###############################################################################
