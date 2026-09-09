###############################################################################
### INApestPointTransitionMatrix / INApestVertebratePoint
###
### Combines explicit continuous-space points with demographic/life stages. A
### point can change stage and, where configured, change x/y position during the
### same transition. Reproductive dispersal and stage-transition movement are
### therefore distinct spatial processes.
###
### The Vertebrate capability layer adds optional Birth, HomeRange, Control and
### Interaction modules without changing the simpler point-transition behaviour
### when Vertebrate = NULL. Pathogen interactions can be supplied through the
### generic point interaction interface.
###############################################################################

# Keep probability values within the valid range from 0 to 1.
.iv_clip01 <- function(x) pmin(1, pmax(0, x))

# Validate the optional Birth/HomeRange/Control/Interaction module list.
.iv_validate_vertebrate <- function(Vertebrate) {
  if (is.null(Vertebrate)) return(invisible(NULL))
  if (!is.list(Vertebrate)) stop("Vertebrate must be NULL or a named list.")
  if (length(Vertebrate) && (is.null(names(Vertebrate)) || any(!nzchar(names(Vertebrate)))))
    stop("Vertebrate must be a named list.")
  allowed <- c("Birth", "HomeRange", "Control", "Interaction")
  unknown <- setdiff(names(Vertebrate), allowed)
  if (length(unknown))
    stop("Unknown Vertebrate component(s): ", paste(unknown, collapse = ", "),
         ". Allowed components are Birth, HomeRange, Control and Interaction.")
  invisible(NULL)
}

# Return one named vertebrate module, or NULL when it is not supplied.
.iv_module <- function(Vertebrate, name) {
  if (is.null(Vertebrate) || is.null(Vertebrate[[name]])) return(NULL)
  Vertebrate[[name]]
}

# Call a user vertebrate hook with only the arguments it accepts.
.iv_call_hook <- function(fun, args, name = "vertebrate hook") {
  if (!is.function(fun)) stop(name, " must be a function.")
  fm <- names(formals(fun))
  if (is.null(fm) || "..." %in% fm) return(do.call(fun, args))
  do.call(fun, args[intersect(names(args), fm)])
}

# Row-bind event tables after filling missing columns with NA.
.iv_rbind_fill <- function(...) {
  xs <- list(...)
  if (length(xs) == 1L && is.list(xs[[1]]) && !is.data.frame(xs[[1]])) xs <- xs[[1]]
  xs <- xs[!vapply(xs, is.null, logical(1))]
  if (!length(xs)) return(data.frame())
  xs <- lapply(xs, function(x) {
    x <- as.data.frame(x, stringsAsFactors = FALSE)
    for (nm in names(x)) if (is.factor(x[[nm]])) x[[nm]] <- as.character(x[[nm]])
    x
  })
  all_names <- unique(unlist(lapply(xs, names), use.names = FALSE))
  xs <- lapply(xs, function(x) {
    miss <- setdiff(all_names, names(x))
    for (nm in miss) x[[nm]] <- rep(NA, nrow(x))
    x[, all_names, drop = FALSE]
  })
  do.call(rbind, xs)
}

.iv_point_reserved <- c(
  "x", "y", "stage", "id", "parent_id", "birth_timestep",
  "have_info", "detected", "managing", "last_known_timestep"
)

# Standardise starting point records and engine-owned state fields.
.iv_initial_point_state <- function(InitialPoints, initial_stage) {
  p <- as.data.frame(InitialPoints, stringsAsFactors = FALSE)
  for (nm in names(p)) if (is.factor(p[[nm]])) p[[nm]] <- as.character(p[[nm]])
  n <- nrow(p)
  p$x <- as.numeric(p$x)
  p$y <- as.numeric(p$y)
  p$stage <- as.integer(initial_stage)
  p$id <- seq_len(n)
  p$parent_id <- rep(NA_integer_, n)
  p$birth_timestep <- rep(0L, n)
  p$have_info <- rep(FALSE, n)
  p$detected <- rep(FALSE, n)
  p$managing <- rep(FALSE, n)
  p$last_known_timestep <- rep(NA_integer_, n)
  # Put engine-owned fields first and retain every other user field after them.
  extras <- setdiff(names(p), .iv_point_reserved)
  p[, c(.iv_point_reserved, extras), drop = FALSE]
}

# Standardise newly arriving external point records.
.iv_external_point_state <- function(ext, next_id, timestep, Nstages) {
  if (!is.data.frame(ext) || !all(c("x", "y") %in% names(ext)))
    stop("ExternalIncursionGenerator must return NULL/zero rows or a data.frame with x and y.")
  ep <- as.data.frame(ext, stringsAsFactors = FALSE)
  for (nm in names(ep)) if (is.factor(ep[[nm]])) ep[[nm]] <- as.character(ep[[nm]])
  if (!"stage" %in% names(ep)) ep$stage <- 1L
  ep$stage <- as.integer(ep$stage)
  if (any(is.na(ep$stage)) || any(ep$stage < 1L | ep$stage > Nstages))
    stop("External incursion stage must be 1..Nstages.")
  n <- nrow(ep)
  ep$x <- as.numeric(ep$x); ep$y <- as.numeric(ep$y)
  ep$id <- seq.int(next_id, length.out = n)
  ep$parent_id <- rep(NA_integer_, n)
  ep$birth_timestep <- rep(as.integer(timestep), n)
  ep$have_info <- rep(FALSE, n)
  ep$detected <- rep(FALSE, n)
  ep$managing <- rep(FALSE, n)
  ep$last_known_timestep <- rep(NA_integer_, n)
  extras <- setdiff(names(ep), .iv_point_reserved)
  ep[, c(.iv_point_reserved, extras), drop = FALSE]
}

# Resolve one home-range movement scale per explicit point.
.iv_home_range_point <- function(module, points, timestep, perm, context = list()) {
  if (is.null(module) || !nrow(points)) return(NULL)
  ans <- if (is.function(module)) {
    .iv_call_hook(module, list(
      points = points, state = points, timestep = timestep, perm = perm,
      context = context
    ), "Vertebrate$HomeRange")
  } else module

  if (is.data.frame(ans)) {
    if (!all(c("id", "sigma") %in% names(ans)))
      stop("Point HomeRange data.frame output must contain id and sigma.")
    sigma <- ans$sigma[match(points$id, ans$id)]
  } else {
    sigma <- as.numeric(ans)
    if (length(sigma) == 1L) sigma <- rep(sigma, nrow(points))
    if (length(sigma) != nrow(points))
      stop("Point HomeRange must resolve to one sigma per point (or a scalar).")
  }
  if (any(is.na(sigma)) || any(!is.finite(sigma)) || any(sigma < 0))
    stop("Home-range sigma values must be finite and non-negative.")
  as.numeric(sigma)
}

# Resolve home-range movement scales by node and life stage.
.iv_home_range_node <- function(module, population, timestep, perm, context = list()) {
  n_nodes <- nrow(population); n_stages <- ncol(population)
  if (is.null(module)) return(NULL)
  ans <- if (is.function(module)) {
    .iv_call_hook(module, list(
      population = population, state = population, timestep = timestep,
      perm = perm, context = context
    ), "Vertebrate$HomeRange")
  } else module
  if (is.matrix(ans)) {
    if (!identical(dim(ans), c(n_nodes, n_stages)))
      stop("Node HomeRange matrix must have dimensions nodes x stages.")
    sigma <- ans
  } else {
    ans <- as.numeric(ans)
    if (length(ans) == 1L) {
      sigma <- matrix(ans, n_nodes, n_stages)
    } else if (length(ans) == n_stages && length(ans) != n_nodes) {
      sigma <- matrix(rep(ans, each = n_nodes), n_nodes, n_stages)
    } else if (length(ans) == n_nodes && length(ans) != n_stages) {
      sigma <- matrix(rep(ans, n_stages), n_nodes, n_stages)
    } else if (length(ans) == n_nodes && length(ans) == n_stages) {
      stop("Node HomeRange vector is ambiguous because nodes equals stages; supply a nodes x stages matrix.")
    } else stop("Node HomeRange must be scalar, length nodes, length stages, or a nodes x stages matrix.")
  }
  if (any(is.na(sigma)) || any(!is.finite(sigma)) || any(sigma < 0))
    stop("Home-range sigma values must be finite and non-negative.")
  sigma
}

# Select control devices active in the current timestep.
.iv_active_devices <- function(devices, timestep) {
  if (is.null(devices)) return(NULL)
  if (!is.data.frame(devices)) stop("Vertebrate$Control$Devices must be a data.frame.")
  d <- as.data.frame(devices, stringsAsFactors = FALSE)
  if (!nrow(d)) return(d)
  start_name <- if ("start" %in% names(d)) "start" else if ("active_from" %in% names(d)) "active_from" else NULL
  end_name <- if ("end" %in% names(d)) "end" else if ("active_to" %in% names(d)) "active_to" else NULL
  keep <- rep(TRUE, nrow(d))
  if (!is.null(start_name)) keep <- keep & (is.na(d[[start_name]]) | d[[start_name]] <= timestep)
  if (!is.null(end_name)) keep <- keep & (is.na(d[[end_name]]) | d[[end_name]] >= timestep)
  d[keep, , drop = FALSE]
}

# Create neutral per-point control effects.
.iv_point_effects_empty <- function(points) {
  data.frame(
    id = points$id,
    kill_prob = rep(0, nrow(points)),
    detect_prob = rep(0, nrow(points)),
    fecundity_reduction = rep(0, nrow(points)),
    stringsAsFactors = FALSE
  )
}

# Normalise custom point-control output to the engine contract.
.iv_normalise_point_effects <- function(ans, points, name = "control model") {
  cost <- 0
  if (is.null(ans)) return(list(effects = .iv_point_effects_empty(points), cost = 0))
  if (is.list(ans) && !is.data.frame(ans) && !is.null(ans$Effects)) {
    cost <- if (is.null(ans$Cost)) 0 else as.numeric(ans$Cost)[1]
    ans <- ans$Effects
  }
  if (is.numeric(ans) && !is.data.frame(ans)) {
    v <- as.numeric(ans)
    if (length(v) == 1L) v <- rep(v, nrow(points))
    if (length(v) != nrow(points)) stop(name, " numeric output must be scalar or one kill probability per point.")
    ans <- data.frame(id = points$id, kill_prob = v)
  }
  if (!is.data.frame(ans)) stop(name, " must return a data.frame, numeric kill probability, or list(Effects=..., Cost=...).")
  if (!"id" %in% names(ans)) {
    if (nrow(ans) != nrow(points)) stop(name, " output without id must have one row per point.")
    ans$id <- points$id
  }
  if (anyDuplicated(ans$id)) stop(name, " output contains duplicate point ids.")
  idx <- match(points$id, ans$id)
  out <- .iv_point_effects_empty(points)
  for (nm in c("kill_prob", "detect_prob", "fecundity_reduction")) {
    if (nm %in% names(ans)) {
      z <- ans[[nm]][idx]
      z[is.na(z)] <- 0
      out[[nm]] <- as.numeric(z)
    }
    if (any(!is.finite(out[[nm]])) || any(out[[nm]] < 0 | out[[nm]] > 1))
      stop(name, " ", nm, " values must be between 0 and 1.")
  }
  if (!is.finite(cost) || cost < 0) stop(name, " Cost must be finite and non-negative.")
  list(effects = out, cost = cost)
}

# Combine point-control effects from multiple control sources.
.iv_combine_point_effects <- function(a, b) {
  out <- a
  out$kill_prob <- 1 - (1 - a$kill_prob) * (1 - b$kill_prob)
  out$detect_prob <- 1 - (1 - a$detect_prob) * (1 - b$detect_prob)
  out$fecundity_reduction <- 1 - (1 - a$fecundity_reduction) * (1 - b$fecundity_reduction)
  out
}

# Calculate default distance-based control effects for points.
.iv_point_device_default <- function(points, devices, home_range) {
  out <- .iv_point_effects_empty(points)
  if (is.null(devices) || !nrow(devices) || !nrow(points)) return(list(effects = out, cost = 0))
  if (!all(c("x", "y") %in% names(devices)))
    stop("Default point-device control requires device x and y columns.")

  d <- devices
  defaults <- list(g0 = 1, effort = 1, kill = 1, detect = 0, fecundity_reduction = 0, cost = 0)
  for (nm in names(defaults)) if (!nm %in% names(d)) d[[nm]] <- defaults[[nm]]
  for (nm in c("g0", "kill", "detect", "fecundity_reduction")) {
    z <- as.numeric(d[[nm]])
    if (any(!is.finite(z)) || any(z < 0 | z > 1)) stop("Device ", nm, " values must be in [0,1].")
    d[[nm]] <- z
  }
  d$effort <- as.numeric(d$effort)
  d$cost <- as.numeric(d$cost)
  if (any(!is.finite(d$effort)) || any(d$effort < 0)) stop("Device effort must be finite and non-negative.")
  if (any(!is.finite(d$cost)) || any(d$cost < 0)) stop("Device cost must be finite and non-negative.")

  for (j in seq_len(nrow(d))) {
    sigma <- if ("sigma" %in% names(d) && is.finite(as.numeric(d$sigma[j]))) {
      rep(as.numeric(d$sigma[j]), nrow(points))
    } else {
      if (is.null(home_range)) stop("Default point-device control requires Vertebrate$HomeRange or a device sigma column.")
      home_range
    }
    if (any(sigma < 0)) stop("Device/home-range sigma must be non-negative.")
    d2 <- (points$x - as.numeric(d$x[j]))^2 + (points$y - as.numeric(d$y[j]))^2
    base <- numeric(nrow(points))
    pos <- sigma > 0
    base[pos] <- d$g0[j] * exp(-d2[pos] / (2 * sigma[pos]^2))
    base[!pos & d2 == 0] <- d$g0[j]
    encounter <- 1 - (1 - .iv_clip01(base))^d$effort[j]
    out$kill_prob <- 1 - (1 - out$kill_prob) * (1 - encounter * d$kill[j])
    out$detect_prob <- 1 - (1 - out$detect_prob) * (1 - encounter * d$detect[j])
    out$fecundity_reduction <- 1 - (1 - out$fecundity_reduction) *
      (1 - encounter * d$fecundity_reduction[j])
  }
  list(effects = out, cost = sum(d$cost * d$effort))
}

# Calculate default area-based control effects for points.
.iv_point_area_default <- function(points, area, timestep, perm, context) {
  if (is.null(area)) return(list(effects = .iv_point_effects_empty(points), cost = 0))
  if (is.function(area)) {
    ans <- .iv_call_hook(area, list(points = points, state = points, timestep = timestep,
                                    perm = perm, context = context), "Vertebrate$Control$Area")
    return(.iv_normalise_point_effects(ans, points, "Vertebrate$Control$Area"))
  }
  if (is.numeric(area)) return(.iv_normalise_point_effects(area, points, "Vertebrate$Control$Area"))
  if (!is.list(area)) stop("Vertebrate$Control$Area must be NULL, numeric, a function, or a list.")
  out <- .iv_point_effects_empty(points)
  resolve <- function(x, label) {
    if (is.null(x)) return(rep(0, nrow(points)))
    if (is.function(x)) x <- .iv_call_hook(x, list(points = points, state = points,
      timestep = timestep, perm = perm, context = context), label)
    x <- as.numeric(x)
    if (length(x) == 1L) x <- rep(x, nrow(points))
    if (length(x) != nrow(points) || any(!is.finite(x)) || any(x < 0 | x > 1))
      stop(label, " must resolve to scalar or one [0,1] value per point.")
    x
  }
  out$kill_prob <- resolve(area$Kill, "Control$Area$Kill")
  out$detect_prob <- resolve(area$Detect, "Control$Area$Detect")
  out$fecundity_reduction <- resolve(area$FecundityReduction, "Control$Area$FecundityReduction")
  cost <- if (is.null(area$Cost)) 0 else {
    x <- if (is.function(area$Cost)) .iv_call_hook(area$Cost, list(points = points, state = points,
      timestep = timestep, perm = perm, context = context), "Control$Area$Cost") else area$Cost
    sum(as.numeric(x))
  }
  if (!is.finite(cost) || cost < 0) stop("Control$Area$Cost must be finite and non-negative.")
  list(effects = out, cost = cost)
}

# Resolve all point-level vertebrate control effects.
.iv_point_control <- function(module, points, home_range, timestep, perm, context = list()) {
  if (is.null(module) || !nrow(points)) return(list(effects = .iv_point_effects_empty(points), cost = 0))
  if (!is.list(module)) stop("Vertebrate$Control must be NULL or a list.")
  devices <- .iv_active_devices(module$Devices, timestep)
  dev <- if (is.function(module$Model)) {
    ans <- .iv_call_hook(module$Model, list(
      points = points, state = points, devices = devices, home_range = home_range,
      timestep = timestep, perm = perm, context = context
    ), "Vertebrate$Control$Model")
    .iv_normalise_point_effects(ans, points, "Vertebrate$Control$Model")
  } else .iv_point_device_default(points, devices, home_range)
  area <- .iv_point_area_default(points, module$Area, timestep, perm, context)
  list(effects = .iv_combine_point_effects(dev$effects, area$effects), cost = dev$cost + area$cost)
}

# Create neutral node x stage control effects.
.iv_node_effects_empty <- function(n_nodes, n_stages) {
  list(
    kill_prob = matrix(0, n_nodes, n_stages),
    detect_prob = matrix(0, n_nodes, n_stages),
    fecundity_reduction = matrix(0, n_nodes, n_stages),
    cost = 0
  )
}

# Resolve node/stage values to a probability matrix.
.iv_node_prob_matrix <- function(x, n_nodes, n_stages, name) {
  if (is.null(x)) return(matrix(0, n_nodes, n_stages))
  if (is.matrix(x)) {
    if (!identical(dim(x), c(n_nodes, n_stages))) stop(name, " matrix must be nodes x stages.")
    out <- x
  } else {
    x <- as.numeric(x)
    if (length(x) == 1L) out <- matrix(x, n_nodes, n_stages)
    else if (length(x) == n_stages && length(x) != n_nodes) out <- matrix(rep(x, each = n_nodes), n_nodes, n_stages)
    else if (length(x) == n_nodes && length(x) != n_stages) out <- matrix(rep(x, n_stages), n_nodes, n_stages)
    else if (length(x) == n_nodes && length(x) == n_stages)
      stop(name, " vector is ambiguous because nodes equals stages; supply a matrix.")
    else stop(name, " must be scalar, length nodes, length stages, or nodes x stages.")
  }
  out <- matrix(as.numeric(out), n_nodes, n_stages)
  if (any(!is.finite(out)) || any(out < 0 | out > 1)) stop(name, " values must be in [0,1].")
  out
}

# Normalise custom node-control output to the engine contract.
.iv_normalise_node_effects <- function(ans, population, name = "control model") {
  n_nodes <- nrow(population); n_stages <- ncol(population)
  out <- .iv_node_effects_empty(n_nodes, n_stages)
  if (is.null(ans)) return(out)
  if (is.numeric(ans) || is.matrix(ans)) {
    out$kill_prob <- .iv_node_prob_matrix(ans, n_nodes, n_stages, paste0(name, " kill_prob"))
    return(out)
  }
  if (!is.list(ans)) stop(name, " must return numeric/matrix kill probabilities or a list.")
  out$kill_prob <- .iv_node_prob_matrix(ans$kill_prob, n_nodes, n_stages, paste0(name, " kill_prob"))
  out$detect_prob <- .iv_node_prob_matrix(ans$detect_prob, n_nodes, n_stages, paste0(name, " detect_prob"))
  out$fecundity_reduction <- .iv_node_prob_matrix(ans$fecundity_reduction, n_nodes, n_stages,
                                                   paste0(name, " fecundity_reduction"))
  out$cost <- if (is.null(ans$cost)) if (is.null(ans$Cost)) 0 else sum(as.numeric(ans$Cost)) else sum(as.numeric(ans$cost))
  if (!is.finite(out$cost) || out$cost < 0) stop(name, " cost must be finite and non-negative.")
  out
}

# Combine node-control effects from multiple sources.
.iv_combine_node_effects <- function(a, b) {
  list(
    kill_prob = 1 - (1 - a$kill_prob) * (1 - b$kill_prob),
    detect_prob = 1 - (1 - a$detect_prob) * (1 - b$detect_prob),
    fecundity_reduction = 1 - (1 - a$fecundity_reduction) * (1 - b$fecundity_reduction),
    cost = a$cost + b$cost
  )
}

# Calculate default device-based node control effects.
.iv_node_device_default <- function(population, devices, home_range) {
  n_nodes <- nrow(population); n_stages <- ncol(population)
  out <- .iv_node_effects_empty(n_nodes, n_stages)
  if (is.null(devices) || !nrow(devices)) return(out)
  if (!"node" %in% names(devices))
    stop("Default node-device control requires a device/programme 'node' column. Use Control$Model for other geometries.")
  d <- devices
  defaults <- list(density = 1, g0 = 1, effort = 1, kill = 1, detect = 0,
                   fecundity_reduction = 0, cost = 0)
  for (nm in names(defaults)) if (!nm %in% names(d)) d[[nm]] <- defaults[[nm]]
  if (any(is.na(d$node)) || any(d$node < 1 | d$node > n_nodes)) stop("Device node values must be 1..number of nodes.")
  for (nm in c("g0", "kill", "detect", "fecundity_reduction")) {
    z <- as.numeric(d[[nm]])
    if (any(!is.finite(z)) || any(z < 0 | z > 1)) stop("Device ", nm, " values must be in [0,1].")
    d[[nm]] <- z
  }
  for (nm in c("density", "effort", "cost")) {
    z <- as.numeric(d[[nm]])
    if (any(!is.finite(z)) || any(z < 0)) stop("Device ", nm, " values must be finite and non-negative.")
    d[[nm]] <- z
  }
  for (j in seq_len(nrow(d))) {
    node <- as.integer(d$node[j])
    sigma <- if ("sigma" %in% names(d) && is.finite(as.numeric(d$sigma[j]))) {
      rep(as.numeric(d$sigma[j]), n_stages)
    } else {
      if (is.null(home_range)) stop("Default node-device control requires Vertebrate$HomeRange or a device sigma column.")
      home_range[node, ]
    }
    if (any(sigma < 0)) stop("Device/home-range sigma must be non-negative.")
    encounter <- 1 - exp(-2 * pi * d$g0[j] * sigma^2 * d$density[j] * d$effort[j])
    encounter <- .iv_clip01(encounter)
    out$kill_prob[node, ] <- 1 - (1 - out$kill_prob[node, ]) * (1 - encounter * d$kill[j])
    out$detect_prob[node, ] <- 1 - (1 - out$detect_prob[node, ]) * (1 - encounter * d$detect[j])
    out$fecundity_reduction[node, ] <- 1 - (1 - out$fecundity_reduction[node, ]) *
      (1 - encounter * d$fecundity_reduction[j])
  }
  out$cost <- sum(d$cost * d$effort)
  out
}

# Calculate default area-based node control effects.
.iv_node_area_default <- function(population, area, timestep, perm, context) {
  n_nodes <- nrow(population); n_stages <- ncol(population)
  if (is.null(area)) return(.iv_node_effects_empty(n_nodes, n_stages))
  if (is.function(area)) {
    ans <- .iv_call_hook(area, list(population = population, state = population, timestep = timestep,
                                    perm = perm, context = context), "Vertebrate$Control$Area")
    return(.iv_normalise_node_effects(ans, population, "Vertebrate$Control$Area"))
  }
  if (is.numeric(area) || is.matrix(area))
    return(.iv_normalise_node_effects(area, population, "Vertebrate$Control$Area"))
  if (!is.list(area)) stop("Vertebrate$Control$Area must be NULL, numeric/matrix, a function, or a list.")
  resolve <- function(x, nm) {
    if (is.function(x)) x <- .iv_call_hook(x, list(population = population, state = population,
      timestep = timestep, perm = perm, context = context), nm)
    .iv_node_prob_matrix(x, n_nodes, n_stages, nm)
  }
  out <- .iv_node_effects_empty(n_nodes, n_stages)
  out$kill_prob <- resolve(area$Kill, "Control$Area$Kill")
  out$detect_prob <- resolve(area$Detect, "Control$Area$Detect")
  out$fecundity_reduction <- resolve(area$FecundityReduction, "Control$Area$FecundityReduction")
  if (!is.null(area$Cost)) {
    z <- if (is.function(area$Cost)) .iv_call_hook(area$Cost, list(population = population, state = population,
      timestep = timestep, perm = perm, context = context), "Control$Area$Cost") else area$Cost
    out$cost <- sum(as.numeric(z))
    if (!is.finite(out$cost) || out$cost < 0) stop("Control$Area$Cost must be finite and non-negative.")
  }
  out
}

# Resolve all node-level vertebrate control effects.
.iv_node_control <- function(module, population, home_range, timestep, perm, context = list()) {
  n_nodes <- nrow(population); n_stages <- ncol(population)
  if (is.null(module)) return(.iv_node_effects_empty(n_nodes, n_stages))
  if (!is.list(module)) stop("Vertebrate$Control must be NULL or a list.")
  devices <- .iv_active_devices(module$Devices, timestep)
  dev <- if (is.function(module$Model)) {
    ans <- .iv_call_hook(module$Model, list(
      population = population, state = population, devices = devices,
      home_range = home_range, timestep = timestep, perm = perm, context = context
    ), "Vertebrate$Control$Model")
    .iv_normalise_node_effects(ans, population, "Vertebrate$Control$Model")
  } else .iv_node_device_default(population, devices, home_range)
  area <- .iv_node_area_default(population, module$Area, timestep, perm, context)
  .iv_combine_node_effects(dev, area)
}

# Apply a vertebrate Birth module to explicit points.
.iv_point_birth <- function(module, points, transition, timestep, perm, context = list()) {
  if (is.null(module)) return(NULL)
  ans <- .iv_call_hook(module, list(
    points = points, state = points, transition = transition,
    timestep = timestep, perm = perm, context = context
  ), "Vertebrate$Birth")
  if (is.null(ans)) return(data.frame(parent_id = integer(0), stage = integer(0)))
  if (is.numeric(ans) && !is.data.frame(ans)) {
    counts <- as.numeric(ans)
    if (length(counts) == 1L) counts <- rep(counts, nrow(points))
    if (length(counts) != nrow(points) || any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts)))
      stop("Point Vertebrate$Birth numeric output must be non-negative whole-number offspring counts per parent.")
    n_offspring <- sum(counts)
    if (n_offspring == 0L) return(data.frame(parent_id = integer(0), stage = integer(0)))
    return(data.frame(
      parent_id = rep(points$id, counts),
      stage = rep(1L, n_offspring),
      stringsAsFactors = FALSE
    ))
  }
  if (!is.data.frame(ans) || !"parent_id" %in% names(ans))
    stop("Point Vertebrate$Birth must return NULL, offspring counts per parent, or a data.frame containing parent_id.")
  out <- as.data.frame(ans, stringsAsFactors = FALSE)
  for (nm in names(out)) if (is.factor(out[[nm]])) out[[nm]] <- as.character(out[[nm]])
  if (!"stage" %in% names(out)) out$stage <- 1L
  out$stage <- as.integer(out$stage)
  if (any(!out$parent_id %in% points$id)) stop("Point Vertebrate$Birth returned a parent_id not present in the current population.")
  out
}

# Apply a vertebrate Birth module to node x stage populations.
.iv_node_birth <- function(module, population, transition, timestep, perm, context = list()) {
  if (is.null(module)) return(NULL)
  ans <- .iv_call_hook(module, list(
    population = population, state = population, transition = transition,
    timestep = timestep, perm = perm, context = context
  ), "Vertebrate$Birth")
  if (is.null(ans)) return(list(mean = rep(0, nrow(population)), mother_counts = NULL))
  if (is.list(ans)) {
    mean <- ans$mean
    mothers <- ans$mother_counts
  } else {
    mean <- ans; mothers <- NULL
  }
  mean <- as.numeric(mean)
  if (length(mean) == 1L) mean <- rep(mean, nrow(population))
  if (length(mean) != nrow(population) || any(!is.finite(mean)) || any(mean < 0))
    stop("Node Vertebrate$Birth must return non-negative expected births per node (or list(mean=..., mother_counts=...)).")
  if (!is.null(mothers)) {
    mothers <- as.numeric(mothers)
    if (length(mothers) == 1L) mothers <- rep(mothers, nrow(population))
    if (length(mothers) != nrow(population) || any(!is.finite(mothers)) || any(mothers < 0))
      stop("Node Vertebrate$Birth mother_counts must be non-negative and length nodes.")
  }
  list(mean = mean, mother_counts = mothers)
}

# Apply custom interaction/social processes to explicit points.
.iv_point_interaction <- function(module, points, home_range, timestep, perm, context = list()) {
  if (is.null(module) || !nrow(points)) return(list(points = points, contacts = data.frame(), events = NULL))
  if (is.function(module)) module <- list(Update = module)
  if (!is.list(module)) stop("Vertebrate$Interaction must be NULL, a function, or a list.")
  contacts <- data.frame()
  if (is.function(module$Contact)) {
    contacts <- .iv_call_hook(module$Contact, list(
      points = points, state = points, home_range = home_range,
      timestep = timestep, perm = perm, context = context
    ), "Vertebrate$Interaction$Contact")
    if (is.null(contacts)) contacts <- data.frame()
    if (!is.data.frame(contacts)) stop("Interaction$Contact must return a data.frame or NULL.")
    if (nrow(contacts) && !all(c("id1", "id2") %in% names(contacts)))
      stop("Interaction$Contact output must contain id1 and id2.")
  }
  updated <- points; extra_events <- NULL
  if (is.function(module$Update)) {
    ans <- .iv_call_hook(module$Update, list(
      points = points, state = points, contacts = contacts, home_range = home_range,
      timestep = timestep, perm = perm, context = context
    ), "Vertebrate$Interaction$Update")
    if (is.list(ans) && !is.data.frame(ans) && !is.null(ans$points)) {
      updated <- ans$points; extra_events <- ans$events
    } else updated <- ans
    if (!is.data.frame(updated)) stop("Interaction$Update must return a point data.frame or list(points=...).")
    if (!all(.iv_point_reserved %in% names(updated)))
      stop("Interaction$Update must retain all INApest reserved point fields.")
    if (anyDuplicated(updated$id) || !setequal(updated$id, points$id))
      stop("Interaction$Update must retain exactly the same point ids; births/deaths belong in other model processes.")
    updated <- updated[match(points$id, updated$id), , drop = FALSE]
  }
  list(points = updated, contacts = contacts, events = extra_events)
}

# Apply custom interaction/social processes to node populations.
.iv_node_interaction <- function(module, population, home_range, timestep, perm, context = list()) {
  if (is.null(module)) return(population)
  fun <- if (is.function(module)) module else if (is.list(module) && is.function(module$Update)) module$Update else NULL
  if (is.null(fun)) return(population)
  ans <- .iv_call_hook(fun, list(
    population = population, state = population, home_range = home_range,
    timestep = timestep, perm = perm, context = context
  ), "Vertebrate$Interaction")
  if (is.list(ans) && !is.null(ans$population)) ans <- ans$population
  if (!is.matrix(ans) || !identical(dim(ans), dim(population)))
    stop("Node Interaction must return a population matrix with unchanged nodes x stages dimensions.")
  if (any(!is.finite(ans)) || any(ans < 0)) stop("Node Interaction population values must be finite and non-negative.")
  floor(ans)
}

###############################################################################
### INApestPointTransitionMatrix
### Point-based analogue of INApestMetaTransitionMatrix for small-area incursions
###
### Each biological point represents exactly ONE individual (e.g. invasive tree)
### or ONE colony/nest (e.g. yellow-legged hornet nest).
###
### Main changes relative to INApestMeta:
# * no spatial nodes;
# * no SDD connectivity matrix;
# * propagules disperse from explicit x/y coordinates using kernels;
# * establishment can be moderated by a native base-R habitat grid;
# * active propagules can nudge from a provisional destination to the
# closest location with the highest better habitat inside a search radius;
# * SEAM-like information spread is represented by distance-based transfer
### from persistent information sites.
###
### Coordinates, habitat-grid extents, kernel distances, HabitatSearchRadius,
# KRadius and InfoRadius must use compatible linear units. GIS libraries may be
### used upstream to prepare the base-R habitat grid, but are not required here.
###############################################################################

# Keep probability values within the valid range from 0 to 1.
.ipp_clip01 <- function(x) pmin(1, pmax(0, x))


# Validate scalar, timestep or resolver probability inputs.
.ipp_validate_probability_schedule <- function(x,
                                               Ntimesteps,
                                               name,
                                               allow_null = FALSE) {

  if (is.null(x)) {
    if (allow_null)
      return(invisible(TRUE))
    stop(name, " cannot be NULL.")
  }

  if (!is.numeric(x) || is.function(x))
    stop(
      name,
      " must be a numeric scalar or numeric vector of length Ntimesteps."
    )

  if (!(length(x) == 1L || length(x) == Ntimesteps))
    stop(
      name,
      " must be a scalar or a numeric vector of length Ntimesteps (",
      Ntimesteps,
      ")."
    )

  if (any(!is.finite(x)))
    stop(name, " must contain only finite values.")

  if (any(x < 0 | x > 1))
    stop(name, " values must all be between 0 and 1 inclusive.")

  invisible(TRUE)
}


# Resolve a point-level input for the current points and timestep.
.ipp_resolve <- function(x, points, timestep, perm, Ntimesteps, name) {

  n <- nrow(points)
  if (n == 0L)
    return(numeric(0))

  if (is.function(x)) {

    out <- x(
      points = points,
      timestep = timestep,
      perm = perm
    )

  } else if (length(x) == 1L) {

    out <- rep(as.numeric(x), n)

  } else if (is.numeric(x) && length(x) == Ntimesteps) {

    out <- rep(as.numeric(x[timestep]), n)

  } else {

    stop(
      name,
      " must be a scalar, a numeric vector of length Ntimesteps, ",
      "or a function(points, timestep, perm)."
    )
  }

  if (length(out) == 1L)
    out <- rep(out, n)

  if (length(out) != n)
    stop(name, " must return one value per supplied point.")

  as.numeric(out)
}


# Resolve and optionally draw point-level probabilities.
.ipp_probability <- function(mean_value,
                             sd_value,
                             points,
                             timestep,
                             perm,
                             Ntimesteps,
                             name) {

  mu <- .ipp_resolve(
    mean_value,
    points,
    timestep,
    perm,
    Ntimesteps,
    name
  )

  if (is.null(sd_value))
    return(.ipp_clip01(mu))

  sd <- .ipp_resolve(
    sd_value,
    points,
    timestep,
    perm,
    Ntimesteps,
    paste0(name, "SD")
  )

  .ipp_clip01(
    stats::rnorm(
      length(mu),
      mean = mu,
      sd = pmax(0, sd)
    )
  )
}


###############################################################################
### Spatial-grid and habitat helpers
###############################################################################

### Generic native base-R spatial modifier grid.
###
### values may be either:
# * a numeric matrix for a static spatial surface; or
# * a numeric 3-D array [row, column, timestep] for a changing surface.
### For repeated grids across timestep ranges, use INApestSpatialSchedule().
###
### Values must be in [0, 1]. NA is permitted and is treated as 0 when the
# model queries the surface. This object is used for habitat suitability and
### for optional spatial multipliers on management parameters.
###
### By default row 1 is the northern/top row (adjacent to ymax), matching the
# usual raster-matrix convention. Set row_origin = "ymin" if row 1 is the
# southern/bottom row instead.
INApestSpatialGrid <- function(values,
                              xmin,             # Minimum x coordinate of the grid extent
                              xmax,             # Maximum x coordinate of the grid extent
                              ymin,             # Minimum y coordinate of the grid extent
                              ymax,             # Maximum y coordinate of the grid extent
                              row_origin = c("ymax", "ymin")) { # Whether row 1 starts at ymax or ymin

  row_origin <- match.arg(row_origin)

  d <- dim(values)
  if (is.null(d) || !(length(d) %in% c(2L, 3L)))
    stop("values must be a numeric matrix or a 3-D numeric array [row, column, timestep].")

  if (!is.numeric(values))
    stop("Spatial-grid values must be numeric.")

  if (d[1] < 1L || d[2] < 1L)
    stop("Spatial grid must contain at least one row and one column.")

  if (length(d) == 3L && d[3] < 1L)
    stop("A time-varying spatial grid must contain at least one timestep layer.")

  finite_values <- values[!is.na(values)]
  if (length(finite_values) && any(!is.finite(finite_values)))
    stop("Spatial-grid values must be finite or NA.")

  if (length(finite_values) && any(finite_values < 0 | finite_values > 1))
    stop("Spatial-grid values must all be between 0 and 1 inclusive.")

  extent <- c(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)
  if (!is.numeric(extent) || length(extent) != 4L || any(!is.finite(extent)))
    stop("xmin, xmax, ymin and ymax must be finite numeric scalars.")

  if (xmax <= xmin || ymax <= ymin)
    stop("Spatial-grid extent requires xmax > xmin and ymax > ymin.")

  nr <- d[1]
  nc <- d[2]

  structure(
    list(
      values = values,
      xmin = as.numeric(xmin),
      xmax = as.numeric(xmax),
      ymin = as.numeric(ymin),
      ymax = as.numeric(ymax),
      nrow = as.integer(nr),
      ncol = as.integer(nc),
      nlayer = if (length(d) == 3L) as.integer(d[3]) else 1L,
      xres = (as.numeric(xmax) - as.numeric(xmin)) / nc,
      yres = (as.numeric(ymax) - as.numeric(ymin)) / nr,
      row_origin = row_origin
    ),
    class = "INApestSpatialGrid"
  )
}


print.INApestSpatialGrid <- function(x, ...) {
  cat(
    "INApestSpatialGrid:",
    x$ncol, "columns x", x$nrow, "rows x", x$nlayer, "layer(s)\n",
    "extent:", x$xmin, x$xmax, x$ymin, x$ymax, "\n",
    "resolution:", x$xres, "x", x$yres, "\n",
    "row 1 origin:", x$row_origin, "\n"
  )
  invisible(x)
}


# Backwards-compatible habitat constructor. Habitat and management surfaces
# deliberately use the same spatial-grid architecture.
INApestHabitatGrid <- function(values,
                              xmin,             # Minimum x coordinate of the grid extent
                              xmax,             # Maximum x coordinate of the grid extent
                              ymin,             # Minimum y coordinate of the grid extent
                              ymax,             # Maximum y coordinate of the grid extent
                              row_origin = c("ymax", "ymin")) { # Whether row 1 starts at ymax or ymin

  out <- INApestSpatialGrid(
    values = values,
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    row_origin = row_origin
  )

  class(out) <- c("INApestHabitatGrid", "INApestSpatialGrid")
  out
}


print.INApestHabitatGrid <- function(x, ...) {
  cat(
    "INApestHabitatGrid:",
    x$ncol, "columns x", x$nrow, "rows x", x$nlayer, "layer(s)\n",
    "extent:", x$xmin, x$xmax, x$ymin, x$ymax, "\n",
    "resolution:", x$xres, "x", x$yres, "\n",
    "row 1 origin:", x$row_origin, "\n"
  )
  invisible(x)
}



# Human-friendly change-point schedule for reusing static spatial grids across
### timestep ranges. The simplest form is:
###
### INApestSpatialSchedule(
# early = grid_early,
# response = grid_response,
# intensive = grid_intensive,
# from = c(1, 6, 11)
### )
###
### which means grid_early applies at timesteps 1:5, grid_response at 6:10,
# and grid_intensive from timestep 11 onward. Use `to` only when explicit
# end timesteps are preferable. The same grid object can be supplied more
# than once if a management pattern returns later.
INApestSpatialSchedule <- function(...,
                                   from = 1,    # First timestep for each scheduled grid
                                   to = NULL,   # Last timestep for each scheduled grid
                                   grids = NULL) { # Optional list of scheduled grids

  dots <- list(...)

  if (!is.null(grids)) {
    if (length(dots) > 0L)
      stop("Supply spatial grids either directly in ... or through grids=, not both.")

    if (inherits(grids, "INApestSpatialGrid") || inherits(grids, "SpatRaster")) {
      dots <- list(grids)
    } else if (is.list(grids)) {
      dots <- grids
    } else {
      stop("grids= must be an INApestSpatialGrid, a SpatRaster, or a list of grids.")
    }
  }

  # Also allow INApestSpatialSchedule(list(early=g1, later=g2), from=...).
  if (length(dots) == 1L &&
      is.list(dots[[1]]) &&
      !inherits(dots[[1]], "INApestSpatialGrid") &&
      !inherits(dots[[1]], "SpatRaster")) {
    dots <- dots[[1]]
  }

  if (length(dots) < 1L)
    stop("Supply at least one spatial grid to INApestSpatialSchedule().")

  ngrid <- length(dots)

  if (!is.numeric(from) || length(from) != ngrid || any(!is.finite(from)))
    stop("from must contain one finite timestep for each supplied grid.")

  if (any(from < 1) || any(from != floor(from)))
    stop("from timesteps must be positive integers.")

  from <- as.integer(from)

  if (is.unsorted(from, strictly = TRUE))
    stop("from timesteps must be strictly increasing.")

  if (is.null(to)) {
    to <- c(from[-1L] - 1L, Inf)
  } else {
    if (!is.numeric(to) || length(to) != ngrid || any(is.na(to)))
      stop("to must contain one end timestep for each supplied grid.")

    finite_to <- is.finite(to)
    if (any(to[finite_to] < 1) || any(to[finite_to] != floor(to[finite_to])))
      stop("Finite to timesteps must be positive integers.")

    if (any(to < from))
      stop("Each to timestep must be greater than or equal to its from timestep.")

    if (ngrid > 1L && any(from[-1L] <= to[-ngrid]))
      stop("Spatial schedule ranges must not overlap.")
  }

  # A schedule deliberately contains static grids. If a unique grid is
  # required at every timestep, users can still supply a 3-D grid directly.
  for (i in seq_len(ngrid)) {
    g <- dots[[i]]

    if (inherits(g, "INApestSpatialGrid")) {
      if (g$nlayer != 1L)
        stop(
          "Grid ", i, " in INApestSpatialSchedule() has ", g$nlayer,
          " layers. Scheduled grids must be static (one layer); use its own ",
          "schedule periods or supply the 3-D grid directly instead."
        )
    } else if (inherits(g, "SpatRaster")) {
      if (requireNamespace("terra", quietly = TRUE) && terra::nlyr(g) != 1L)
        stop("Scheduled SpatRaster grids must contain one layer each.")
    } else {
      stop(
        "Every scheduled item must be an INApestSpatialGrid (or optionally a one-layer SpatRaster)."
      )
    }
  }

  # For native base-R grids, require identical geometry. This catches a very
  # easy parameterisation mistake while still allowing the cell values to
  # change between operational periods.
  native <- vapply(dots, inherits, logical(1), what = "INApestSpatialGrid")
  if (sum(native) > 1L) {
    ref <- dots[[which(native)[1L]]]
    for (i in which(native)[-1L]) {
      g <- dots[[i]]
      same_geometry <- isTRUE(all.equal(
        c(ref$xmin, ref$xmax, ref$ymin, ref$ymax, ref$nrow, ref$ncol),
        c(g$xmin, g$xmax, g$ymin, g$ymax, g$nrow, g$ncol),
        tolerance = 0
      )) && identical(ref$row_origin, g$row_origin)

      if (!same_geometry)
        stop("All INApestSpatialGrid objects in one schedule must have the same extent, dimensions and row_origin.")
    }
  }

  labels <- names(dots)
  if (is.null(labels)) labels <- rep("", ngrid)
  blank <- is.na(labels) | !nzchar(labels)
  labels[blank] <- paste0("grid", which(blank))
  names(dots) <- labels

  structure(
    list(
      grids = dots,
      from = from,
      to = as.numeric(to),
      labels = labels
    ),
    class = "INApestSpatialSchedule"
  )
}


print.INApestSpatialSchedule <- function(x, ...) {
  cat("INApestSpatialSchedule:\n")

  for (i in seq_along(x$grids)) {
    period <- if (is.infinite(x$to[i])) {
      paste0("timestep ", x$from[i], " onward")
    } else if (x$from[i] == x$to[i]) {
      paste0("timestep ", x$from[i])
    } else {
      paste0("timesteps ", x$from[i], "-", x$to[i])
    }

    cat("  ", x$labels[i], ": ", period, "\n", sep = "")
  }

  invisible(x)
}


.ipp_is_spatial_grid <- function(x) inherits(x, "INApestSpatialGrid")
.ipp_is_spatial_schedule <- function(x) inherits(x, "INApestSpatialSchedule")
.ipp_is_habitat_grid <- function(x) .ipp_is_spatial_grid(x)


.ipp_schedule_index <- function(schedule,
                                timestep,
                                name = "SpatialSurface") {

  if (length(timestep) != 1L || !is.finite(timestep) || timestep < 1)
    stop("timestep must be a positive finite scalar when a spatial schedule is queried.")

  idx <- which(
    timestep >= schedule$from &
      timestep <= schedule$to
  )

  if (length(idx) == 0L)
    stop(
      name, " has no spatial grid assigned to timestep ", timestep,
      ". Adjust the from/to values in INApestSpatialSchedule()."
    )

  if (length(idx) > 1L)
    stop(name, " assigns more than one spatial grid to timestep ", timestep, ".")

  idx
}


# Select the spatial grid/layer used in the current timestep.
.ipp_spatial_layer <- function(SpatialSurface, timestep, name = "SpatialSurface") {

  if (is.null(SpatialSurface) || is.function(SpatialSurface))
    return(SpatialSurface)

  if (.ipp_is_spatial_schedule(SpatialSurface)) {
    i <- .ipp_schedule_index(SpatialSurface, timestep, name)
    return(.ipp_spatial_layer(SpatialSurface$grids[[i]], timestep, name))
  }

  if (.ipp_is_spatial_grid(SpatialSurface)) {

    h <- SpatialSurface

    if (h$nlayer == 1L)
      return(h)

    if (timestep > h$nlayer)
      stop(
        "A time-varying ", name, " grid must have at least Ntimesteps layers."
      )

    h$values <- matrix(
      h$values[, , timestep],
      nrow = h$nrow,
      ncol = h$ncol
    )
    h$nlayer <- 1L
    return(h)
  }

  # Optional compatibility path only: terra is not required by the model.
  if (inherits(SpatialSurface, "SpatRaster")) {

    if (!requireNamespace("terra", quietly = TRUE))
      stop(
        "Package 'terra' is required only when ", name, " is supplied as a SpatRaster. ",
        "Use INApestSpatialGrid() for the native base-R path."
      )

    if (terra::nlyr(SpatialSurface) == 1L)
      return(SpatialSurface)

    if (timestep > terra::nlyr(SpatialSurface))
      stop(
        "A multi-layer ", name, " SpatRaster must have at least Ntimesteps layers."
      )

    return(SpatialSurface[[timestep]])
  }

  stop(
    name, " must be NULL, an INApestSpatialGrid (or INApestHabitatGrid), ",
    "an INApestSpatialSchedule, a function(x, y, timestep), or optionally a terra SpatRaster."
  )
}


# Validate a spatial surface across the simulation horizon.
.ipp_validate_spatial_surface <- function(SpatialSurface,
                                          Ntimesteps,
                                          name) {
  if (is.null(SpatialSurface))
    return(invisible(TRUE))

  if (.ipp_is_spatial_schedule(SpatialSurface)) {
    # Check every modelled timestep so schedules fail early with a useful
    # message if a range has accidentally been left uncovered.
    for (tt in seq_len(Ntimesteps))
      invisible(.ipp_spatial_layer(SpatialSurface, tt, name))
    return(invisible(TRUE))
  }

  invisible(.ipp_spatial_layer(SpatialSurface, Ntimesteps, name))
  invisible(TRUE)
}


# Return x/y coordinates for native grid-cell centres.
.ipp_grid_cell_centres <- function(h) {

  x <- h$xmin + (seq_len(h$ncol) - 0.5) * h$xres

  y <- if (identical(h$row_origin, "ymax")) {
    h$ymax - (seq_len(h$nrow) - 0.5) * h$yres
  } else {
    h$ymin + (seq_len(h$nrow) - 0.5) * h$yres
  }

  list(x = x, y = y)
}


# Read spatial-surface values at explicit point coordinates.
.ipp_spatial_value <- function(x,
                               y,
                               SpatialSurface,
                               timestep,
                               name = "SpatialSurface") {

  if (length(x) == 0L)
    return(numeric(0))

  if (length(x) != length(y))
    stop("x and y must have the same length when a spatial surface is queried.")

  # NULL is the neutral multiplier everywhere.
  if (is.null(SpatialSurface))
    return(rep(1, length(x)))

  h <- .ipp_spatial_layer(SpatialSurface, timestep, name)

  if (is.function(h)) {

    v <- h(
      x = x,
      y = y,
      timestep = timestep
    )

    if (length(v) == 1L)
      v <- rep(v, length(x))

    if (length(v) != length(x))
      stop(name, " function must return one value per coordinate pair.")

  } else if (.ipp_is_spatial_grid(h)) {

    # Outside the supplied spatial extent is zero by default. For habitat
    # that means unsuitable; for management it means no local coverage.
    v <- rep(0, length(x))

    inside <- is.finite(x) & is.finite(y) &
      x >= h$xmin & x <= h$xmax &
      y >= h$ymin & y <= h$ymax

    if (any(inside)) {

      col <- floor((x[inside] - h$xmin) / h$xres) + 1L
      col <- pmin(h$ncol, pmax(1L, col))

      row <- if (identical(h$row_origin, "ymax")) {
        floor((h$ymax - y[inside]) / h$yres) + 1L
      } else {
        floor((y[inside] - h$ymin) / h$yres) + 1L
      }
      row <- pmin(h$nrow, pmax(1L, row))

      v[inside] <- h$values[cbind(row, col)]
    }

  } else {

    z <- terra::extract(
      h,
      cbind(x, y)
    )

    value_name <- names(h)[1]

    if (!(value_name %in% names(z)))
      stop("Could not identify the spatial-value column returned by terra::extract.")

    v <- as.numeric(z[[value_name]])
  }

  v[is.na(v)] <- 0

  if (any(!is.finite(v)))
    stop(name, " returned non-finite values.")

  if (any(v < 0 | v > 1))
    stop(name, " values must all be between 0 and 1 inclusive.")

  as.numeric(v)
}


# Management probability/effect with an optional spatial multiplier.
# The spatial surface modifies the mean first. Standard deviation remains the
# supplied uncertainty/variation around that local mean. A spatial multiplier
# of exactly zero is forced to a realised value of exactly zero so stochastic
# SD cannot create management activity outside the covered area.
.ipp_probability_spatial <- function(mean_value,
                                     sd_value,
                                     SpatialSurface,
                                     points,
                                     timestep,
                                     perm,
                                     Ntimesteps,
                                     name) {

  mu <- .ipp_resolve(
    mean_value,
    points,
    timestep,
    perm,
    Ntimesteps,
    name
  )

  spatial <- .ipp_spatial_value(
    points$x,
    points$y,
    SpatialSurface,
    timestep,
    paste0(name, "Spatial")
  )

  mu <- .ipp_clip01(mu * spatial)

  if (is.null(sd_value))
    return(mu)

  sd <- .ipp_resolve(
    sd_value,
    points,
    timestep,
    perm,
    Ntimesteps,
    paste0(name, "SD")
  )

  out <- .ipp_clip01(
    stats::rnorm(
      length(mu),
      mean = mu,
      sd = pmax(0, sd)
    )
  )

  out[spatial <= 0] <- 0
  out
}


# Habitat-specific wrappers retained so existing model code and user scripts
# continue to work. Habitat uses exactly the same grid architecture.
.ipp_habitat_layer <- function(HabitatSuitability, timestep) {
  .ipp_spatial_layer(HabitatSuitability, timestep, "HabitatSuitability")
}


.ipp_habitat_value <- function(x, y, HabitatSuitability, timestep) {
  .ipp_spatial_value(x, y, HabitatSuitability, timestep, "HabitatSuitability")
}

# Search within a radius for a better habitat location.
.ipp_habitat_search <- function(x,
                                y,
                                HabitatSuitability,
                                HabitatSearchRadius,
                                timestep,
                                HabitatSearchCandidates = 128L) {

  n <- length(x)

  if (n == 0L) {
    return(
      data.frame(
        x = numeric(0),
        y = numeric(0),
        habitat = numeric(0),
        habitat_nudged = logical(0)
      )
    )
  }

  if (length(y) != n)
    stop("x and y must have the same length for habitat search.")

  current <- .ipp_habitat_value(
    x,
    y,
    HabitatSuitability,
    timestep
  )

  out_x <- x
  out_y <- y
  out_h <- current
  nudged <- rep(FALSE, n)

  if (is.null(HabitatSuitability) || HabitatSearchRadius <= 0) {
    return(
      data.frame(
        x = out_x,
        y = out_y,
        habitat = out_h,
        habitat_nudged = nudged
      )
    )
  }

  h <- .ipp_habitat_layer(
    HabitatSuitability,
    timestep
  )

  ###########################################################################
  ### Native base-R grid search
  ###
  ### Candidate habitat locations are grid-cell centres inside the biological
  # search radius around the provisional dispersal destination. The point
  # moves only if a cell has strictly better suitability. It chooses the
  # highest suitability available; ties are resolved by shortest distance.
  ### Consequently the realised nudge can never exceed HabitatSearchRadius.
  ###########################################################################

  if (.ipp_is_habitat_grid(h)) {

    centres <- .ipp_grid_cell_centres(h)
    tol <- sqrt(.Machine$double.eps)
    r2 <- HabitatSearchRadius^2

    for (i in seq_len(n)) {

      candidate_cols <- which(
        abs(centres$x - x[i]) <= HabitatSearchRadius + tol
      )

      candidate_rows <- which(
        abs(centres$y - y[i]) <= HabitatSearchRadius + tol
      )

      if (!length(candidate_cols) || !length(candidate_rows))
        next

      rr <- rep(candidate_rows, times = length(candidate_cols))
      cc <- rep(candidate_cols, each = length(candidate_rows))

      d2 <- (centres$x[cc] - x[i])^2 +
        (centres$y[rr] - y[i])^2

      in_radius <- d2 <= r2 + tol
      if (!any(in_radius))
        next

      rr <- rr[in_radius]
      cc <- cc[in_radius]
      d2 <- d2[in_radius]

      hv <- h$values[cbind(rr, cc)]
      hv[is.na(hv)] <- 0
      hv <- .ipp_clip01(hv)

      best <- max(hv)

      if (!is.finite(best) || best <= current[i] + tol)
        next

      best_rows <- which(abs(hv - best) <= tol)
      j <- best_rows[which.min(d2[best_rows])]

      out_x[i] <- centres$x[cc[j]]
      out_y[i] <- centres$y[rr[j]]
      out_h[i] <- hv[j]
      nudged[i] <- TRUE
    }

  } else if (inherits(h, "SpatRaster")) {

    # Optional terra compatibility. This path is not needed when the native
    # base-R grid is used.
    value_name <- names(h)[1]

    for (i in seq_len(n)) {

      p <- terra::vect(
        data.frame(x = x[i], y = y[i]),
        geom = c("x", "y"),
        crs = terra::crs(h)
      )

      search_area <- terra::buffer(
        p,
        width = HabitatSearchRadius
      )

      z <- terra::extract(
        h,
        search_area,
        cells = TRUE,
        xy = TRUE,
        ID = FALSE
      )

      if (is.null(z) ||
          nrow(z) == 0L ||
          !(value_name %in% names(z)))
        next

      z <- z[
        is.finite(z[[value_name]]),
        ,
        drop = FALSE
      ]

      if (nrow(z) == 0L)
        next

      # Enforce the biological radius on the destination itself. Raster cells
      # intersecting a buffer can have centres beyond its radius.
      d2 <- (z$x - x[i])^2 + (z$y - y[i])^2
      keep <- d2 <= HabitatSearchRadius^2 + sqrt(.Machine$double.eps)
      z <- z[keep, , drop = FALSE]
      d2 <- d2[keep]

      if (nrow(z) == 0L)
        next

      best <- max(z[[value_name]], na.rm = TRUE)

      if (!is.finite(best) ||
          best <= current[i] + sqrt(.Machine$double.eps))
        next

      best_rows <- which(
        abs(z[[value_name]] - best) <= sqrt(.Machine$double.eps)
      )

      j <- best_rows[which.min(d2[best_rows])]

      out_x[i] <- z$x[j]
      out_y[i] <- z$y[j]
      out_h[i] <- .ipp_clip01(z[[value_name]][j])
      nudged[i] <- TRUE
    }

  } else {

    # Generic habitat-function fallback. This approximates the best habitat
    # with candidate points sampled uniformly inside the search disc.
    HabitatSearchCandidates <- max(
      1L,
      as.integer(HabitatSearchCandidates)
    )

    for (i in seq_len(n)) {

      angle <- stats::runif(
        HabitatSearchCandidates,
        0,
        2 * pi
      )

      dist <- HabitatSearchRadius *
        sqrt(stats::runif(HabitatSearchCandidates))

      cx <- c(
        x[i],
        x[i] + dist * cos(angle)
      )

      cy <- c(
        y[i],
        y[i] + dist * sin(angle)
      )

      hv <- .ipp_habitat_value(
        cx,
        cy,
        h,
        timestep
      )

      best <- max(hv)

      if (best <= current[i] + sqrt(.Machine$double.eps))
        next

      best_rows <- which(
        abs(hv - best) <= sqrt(.Machine$double.eps)
      )

      d2 <- (cx[best_rows] - x[i])^2 +
        (cy[best_rows] - y[i])^2

      j <- best_rows[which.min(d2)]

      out_x[i] <- cx[j]
      out_y[i] <- cy[j]
      out_h[i] <- hv[j]
      nudged[i] <- TRUE
    }
  }

  data.frame(
    x = out_x,
    y = out_y,
    habitat = out_h,
    habitat_nudged = nudged
  )
}


###############################################################################
### Dispersal helpers
###############################################################################

# Draw x/y displacement from a point dispersal kernel.
.ipp_draw_displacement <- function(kernel,
                                   parents,
                                   timestep,
                                   perm,
                                   name) {

  n <- nrow(parents)

  if (n == 0L)
    return(
      data.frame(
        dx = numeric(0),
        dy = numeric(0)
      )
    )

  if (!is.function(kernel))
    stop(name, " must be a function.")

  z <- kernel(
    n = n,
    parents = parents,
    timestep = timestep,
    perm = perm
  )

  if (is.matrix(z))
    z <- as.data.frame(z)

  if (!is.data.frame(z))
    stop(
      name,
      " must return a data.frame or matrix with columns dx and dy."
    )

  if (!all(c("dx", "dy") %in% names(z)))
    stop(name, " output must contain columns dx and dy.")

  if (nrow(z) != n)
    stop(name, " must return one displacement per propagule.")

  if (any(!is.finite(z$dx)) ||
      any(!is.finite(z$dy)))
    stop(name, " returned non-finite displacements.")

  z[, c("dx", "dy"), drop = FALSE]
}


###############################################################################
### Information helpers
###############################################################################

# Create the empty information-site history used by point models.
.ipp_empty_info_sites <- function() {

  data.frame(
    info_id = integer(0),
    source_point_id = integer(0),
    x = numeric(0),
    y = numeric(0),
    created_timestep = integer(0),
    last_known_timestep = integer(0),
    active = logical(0),
    stringsAsFactors = FALSE
  )
}


# Add or refresh point locations that hold response information.
.ipp_add_or_refresh_info_sites <- function(info_sites,
                                           locations,
                                           timestep) {

  if (nrow(locations) == 0L)
    return(info_sites)

  if (!"id" %in% names(locations))
    locations$id <- NA_integer_

  for (i in seq_len(nrow(locations))) {

    source_id <- as.integer(locations$id[i])
    existing <- integer(0)

    if (!is.na(source_id) &&
        nrow(info_sites) > 0L) {

      existing <- which(
        !is.na(info_sites$source_point_id) &
          info_sites$source_point_id == source_id
      )
    }

    if (length(existing) > 0L) {

      j <- existing[1]

      info_sites$x[j] <- locations$x[i]
      info_sites$y[j] <- locations$y[i]
      info_sites$last_known_timestep[j] <- timestep
      info_sites$active[j] <- TRUE

    } else {

      new_id <- if (nrow(info_sites) == 0L)
        1L
      else
        max(info_sites$info_id) + 1L

      info_sites <- rbind(
        info_sites,
        data.frame(
          info_id = new_id,
          source_point_id = source_id,
          x = locations$x[i],
          y = locations$y[i],
          created_timestep = timestep,
          last_known_timestep = timestep,
          active = TRUE,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  info_sites
}


# Transfer information among nearby explicit points.
.ipp_transfer_information <- function(points,
                                      info_sites,
                                      timestep,
                                      InfoRadius,
                                      InfoTransferProb,
                                      InfoKernel,
                                      perm,
                                      Ntimesteps) {

  if (nrow(points) == 0L)
    return(points)

  # Information can be transmitted from either:
  # 1) extant biological points that currently have information; or
  # 2) persistent information sites left by detection/management history.
  # This mirrors the INApestMeta idea that an invaded, informed location can
  # communicate information, while still allowing a removed point to leave a
  # temporary information footprint.
  src_sites <- info_sites[
    info_sites$active,
    ,
    drop = FALSE
  ]

  src_points <- points[
    points$have_info,
    c("id", "x", "y"),
    drop = FALSE
  ]

  if (nrow(src_points) > 0L) {

    names(src_points)[names(src_points) == "id"] <- "source_point_id"
    src_points$info_id <- NA_integer_
    src_points$created_timestep <- NA_integer_
    src_points$last_known_timestep <- NA_integer_
    src_points$active <- TRUE

    # Avoid double counting an extant informed point that also has a
    # persistent detection site at the same source point ID.
    if (nrow(src_sites) > 0L) {
      src_sites <- src_sites[
        is.na(src_sites$source_point_id) |
          !(src_sites$source_point_id %in% src_points$source_point_id),
        ,
        drop = FALSE
      ]
    }
  }

  src <- rbind(
    src_sites[, c(
      "info_id",
      "source_point_id",
      "x",
      "y",
      "created_timestep",
      "last_known_timestep",
      "active"
    ), drop = FALSE],
    if (nrow(src_points) > 0L)
      src_points[, c(
        "info_id",
        "source_point_id",
        "x",
        "y",
        "created_timestep",
        "last_known_timestep",
        "active"
      ), drop = FALSE]
    else
      .ipp_empty_info_sites()
  )

  if (nrow(src) == 0L)
    return(points)

  targets <- which(!points$have_info)

  if (length(targets) == 0L)
    return(points)

  for (j in targets) {

    d <- sqrt(
      (src$x - points$x[j])^2 +
        (src$y - points$y[j])^2
    )

    if (!is.null(InfoKernel)) {

      p_each <- InfoKernel(
        distance = d,
        source = src,
        target = points[j, , drop = FALSE],
        timestep = timestep,
        perm = perm
      )

      if (length(p_each) == 1L)
        p_each <- rep(p_each, length(d))

      if (length(p_each) != length(d))
        stop(
          "InfoKernel must return one probability ",
          "per source-target distance."
        )

      p_each <- .ipp_clip01(p_each)

    } else {

      if (InfoRadius <= 0)
        next

      p0 <- .ipp_resolve(
        InfoTransferProb,
        points[j, , drop = FALSE],
        timestep,
        perm,
        Ntimesteps,
        "InfoTransferProb"
      )

      p_each <- ifelse(
        d <= InfoRadius,
        p0[1],
        0
      )
    }

    p_transfer <- 1 - prod(1 - p_each)

    if (stats::rbinom(
      1L,
      1L,
      p_transfer
    ) == 1L)
      points$have_info[j] <- TRUE
  }

  points
}


###############################################################################
### Optional local carrying-capacity analogue
###############################################################################

# Apply the point model local weighted-capacity limit.
.ipp_apply_local_capacity <- function(candidates,
                                      existing_points,
                                      LocalK,
                                      KRadius,
                                      timestep,
                                      perm,
                                      Ntimesteps) {

  if (nrow(candidates) == 0L)
    return(logical(0))

  if (!is.function(LocalK) &&
      length(LocalK) == 1L &&
      is.infinite(LocalK))
    return(rep(TRUE, nrow(candidates)))

  if (KRadius <= 0)
    stop(
      "Finite LocalK requires KRadius > 0 because points ",
      "do not have an intrinsic node area."
    )

  k <- .ipp_resolve(
    LocalK,
    candidates,
    timestep,
    perm,
    Ntimesteps,
    "LocalK"
  )

  k <- pmax(
    0,
    floor(k)
  )

  keep <- rep(
    FALSE,
    nrow(candidates)
  )

  # Random sequential acceptance avoids a fixed row-order advantage.
  candidate_order <- sample(
    seq_len(nrow(candidates))
  )

  accepted_xy <- existing_points[
    ,
    c("x", "y"),
    drop = FALSE
  ]

  for (ii in candidate_order) {

    if (k[ii] <= 0)
      next

    n_local <- 0L

    if (nrow(accepted_xy) > 0L) {

      d2 <- (accepted_xy$x - candidates$x[ii])^2 +
        (accepted_xy$y - candidates$y[ii])^2

      n_local <- sum(
        d2 <= KRadius^2
      )
    }

    if (n_local < k[ii]) {

      keep[ii] <- TRUE

      accepted_xy <- rbind(
        accepted_xy,
        candidates[ii, c("x", "y"), drop = FALSE]
      )
    }
  }

  keep
}


###############################################################################
### Convenience radial dispersal kernels
###############################################################################

# Create an isotropic exponential-distance displacement kernel.
INApestPointKernelExponential <- function(mean_distance) {

  force(mean_distance)

  function(n,
           parents = NULL,
           timestep = NULL,
           perm = NULL) {

    d <- stats::rexp(
      n,
      rate = 1 / mean_distance
    )

    a <- stats::runif(
      n,
      0,
      2 * pi
    )

    data.frame(
      dx = d * cos(a),
      dy = d * sin(a)
    )
  }
}


# Create an isotropic lognormal-distance displacement kernel.
INApestPointKernelLognormal <- function(meanlog,
                                       sdlog) {

  force(meanlog)
  force(sdlog)

  function(n,
           parents = NULL,
           timestep = NULL,
           perm = NULL) {

    d <- stats::rlnorm(
      n,
      meanlog = meanlog,
      sdlog = sdlog
    )

    a <- stats::runif(
      n,
      0,
      2 * pi
    )

    data.frame(
      dx = d * cos(a),
      dy = d * sin(a)
    )
  }
}


# Create an isotropic fixed-distance displacement kernel.
INApestPointKernelFixed <- function(distance) {

  force(distance)

  function(n,
           parents = NULL,
           timestep = NULL,
           perm = NULL) {

    a <- stats::runif(
      n,
      0,
      2 * pi
    )

    data.frame(
      dx = distance * cos(a),
      dy = distance * sin(a)
    )
  }
}


###############################################################################
### Event helper
###############################################################################

# Build a standard point-model event-log row set.
.ipp_event <- function(perm,
                       timestep,
                       event,
                       point_id,
                       parent_id,
                       x,
                       y,
                       detail,
                       provisional_x = NA_real_,
                       provisional_y = NA_real_,
                       habitat = NA_real_,
                       habitat_nudged = NA) {

  n <- length(point_id)

  data.frame(
    perm = rep_len(perm, n),
    timestep = rep_len(timestep, n),
    event = rep_len(event, n),
    point_id = point_id,
    parent_id = rep_len(parent_id, n),
    x = x,
    y = y,
    detail = rep_len(detail, n),
    provisional_x = rep_len(provisional_x, n),
    provisional_y = rep_len(provisional_y, n),
    habitat = rep_len(habitat, n),
    habitat_nudged = rep_len(habitat_nudged, n),
    stringsAsFactors = FALSE
  )
}


###############################################################################
### Main function
###############################################################################


###############################################################################
### Point-based transition-matrix helpers
###############################################################################

# Resolve the stage-transition matrix for one point and timestep.
.ipptm_get_transition <- function(Transition, timestep, perm, Nstages) {

  A <- if (is.function(Transition)) {
    Transition(timestep = timestep, perm = perm)
  } else if (is.matrix(Transition)) {
    Transition
  } else if (is.array(Transition) && length(dim(Transition)) == 3L) {
    if (!identical(dim(Transition)[1:2], c(Nstages, Nstages)))
      stop("Transition array must have dimensions Nstages x Nstages x Ntimesteps.")
    if (dim(Transition)[3] < timestep)
      stop("Transition array has fewer time slices than Ntimesteps.")
    Transition[, , timestep]
  } else {
    stop("Transition must be an Nstages x Nstages matrix, a 3D time array, or a function(timestep, perm).")
  }

  if (!is.matrix(A) || !identical(dim(A), c(Nstages, Nstages)))
    stop("Each resolved Transition must be an Nstages x Nstages matrix.")

  if (any(!is.finite(A)) || any(A < 0))
    stop("Transition entries must be finite and non-negative.")

  # The current INApestMetaTransitionMatrix convention is retained:
  # A[1, 2:S]       = fecundity into stage 1
  # diag(A)[1:S-1]  = stasis
  # A[s+1, s]       = progression from stage s to s+1
  # A[S, S]         = terminal-stage survival/stasis
  # Other entries are deliberately rejected so ignored matrix entries cannot
  # silently imply biological transitions that this point engine does not use.
  allowed <- matrix(FALSE, Nstages, Nstages)
  diag(allowed) <- TRUE
  if (Nstages > 1L) {
    allowed[1, 2:Nstages] <- TRUE
    allowed[cbind(2:Nstages, 1:(Nstages - 1L))] <- TRUE
  }
  if (any(A[!allowed] != 0))
    stop("Transition contains non-zero entries outside fecundity, diagonal stasis/survival, or adjacent progression positions.")

  if (Nstages > 1L) {
    for (s in seq_len(Nstages - 1L)) {
      if (A[s, s] > 1 || A[s + 1L, s] > 1 || A[s, s] + A[s + 1L, s] > 1 + sqrt(.Machine$double.eps))
        stop("For each non-terminal source stage, stasis + progression must be between 0 and 1.")
    }
  }
  if (A[Nstages, Nstages] > 1)
    stop("Terminal-stage survival/stasis probability must be between 0 and 1.")

  A
}


# Resolve a stage-specific probability/value schedule.
.ipptm_resolve_stage_schedule <- function(x,
                                          points,
                                          timestep,
                                          perm,
                                          Ntimesteps,
                                          Nstages,
                                          name) {
  n <- nrow(points)
  if (n == 0L) return(numeric(0))

  if (is.function(x)) {
    out <- x(points = points, timestep = timestep, perm = perm)
    if (length(out) == 1L) out <- rep(out, n)
    if (length(out) != n) stop(name, " function must return one value per point.")
    return(as.numeric(out))
  }

  if (!is.numeric(x))
    stop(name, " must be numeric or a function(points, timestep, perm).")

  if (length(x) == 1L)
    return(rep(as.numeric(x), n))

  # Matrix is the unambiguous way to specify stage x timestep schedules.
  if (is.matrix(x)) {
    if (!identical(dim(x), c(Nstages, Ntimesteps)))
      stop(name, " matrix must have dimensions Nstages x Ntimesteps.")
    return(as.numeric(x[points$stage, timestep]))
  }

  # Preserve transition-matrix style: a vector of Nstages means stage-specific.
  # If Nstages == Ntimesteps, this interpretation wins; use a matrix for a
  # timestep schedule in that ambiguous case.
  if (length(x) == Nstages)
    return(as.numeric(x[points$stage]))

  if (length(x) == Ntimesteps)
    return(rep(as.numeric(x[timestep]), n))

  stop(name, " must be scalar, length Nstages, length Ntimesteps, Nstages x Ntimesteps matrix, or a function.")
}


# Resolve and draw stage-specific point probabilities.
.ipptm_probability <- function(mean_value,
                               sd_value,
                               points,
                               timestep,
                               perm,
                               Ntimesteps,
                               Nstages,
                               name) {
  mu <- .ipptm_resolve_stage_schedule(
    mean_value, points, timestep, perm, Ntimesteps, Nstages, name
  )
  if (is.null(sd_value)) return(.ipp_clip01(mu))
  sd <- .ipptm_resolve_stage_schedule(
    sd_value, points, timestep, perm, Ntimesteps, Nstages, paste0(name, "SD")
  )
  .ipp_clip01(stats::rnorm(length(mu), mean = mu, sd = pmax(0, sd)))
}



# Combine stage-specific probability with a spatial multiplier.
.ipptm_probability_spatial <- function(mean_value,
                                       sd_value,
                                       SpatialSurface,
                                       points,
                                       timestep,
                                       perm,
                                       Ntimesteps,
                                       Nstages,
                                       name) {
  mu <- .ipptm_resolve_stage_schedule(
    mean_value, points, timestep, perm, Ntimesteps, Nstages, name
  )

  spatial <- .ipp_spatial_value(
    points$x, points$y, SpatialSurface, timestep, paste0(name, "Spatial")
  )

  mu <- .ipp_clip01(mu * spatial)
  if (is.null(sd_value)) return(mu)

  sd <- .ipptm_resolve_stage_schedule(
    sd_value, points, timestep, perm, Ntimesteps, Nstages, paste0(name, "SD")
  )

  out <- .ipp_clip01(stats::rnorm(length(mu), mean = mu, sd = pmax(0, sd)))
  out[spatial <= 0] <- 0
  out
}


# Validate probability-like point-transition inputs.
.ipptm_validate_prob_like <- function(x, Ntimesteps, Nstages, name, allow_null = FALSE) {
  if (is.null(x)) {
    if (allow_null) return(invisible(TRUE))
    stop(name, " cannot be NULL.")
  }
  if (is.function(x)) return(invisible(TRUE))
  if (!is.numeric(x)) stop(name, " must be numeric or a function.")
  if (is.matrix(x)) {
    if (!identical(dim(x), c(Nstages, Ntimesteps)))
      stop(name, " matrix must have dimensions Nstages x Ntimesteps.")
  } else if (!(length(x) %in% unique(c(1L, Nstages, Ntimesteps)))) {
    stop(name, " must be scalar, length Nstages, length Ntimesteps, or Nstages x Ntimesteps matrix.")
  }
  if (any(!is.finite(x))) stop(name, " must contain finite values.")
  if (any(x < 0 | x > 1)) stop(name, " values must be between 0 and 1 inclusive.")
  invisible(TRUE)
}


# Return the reproductive dispersal kernel for a source stage.
.ipptm_get_kernel <- function(kernels, source_stage, Nstages, name, allow_null = FALSE) {
  if (is.null(kernels)) {
    if (allow_null) return(NULL)
    stop(name, " cannot be NULL.")
  }
  if (is.function(kernels)) return(kernels)
  if (!is.list(kernels) || length(kernels) != Nstages)
    stop(name, " must be a function or a list of length Nstages.")
  k <- kernels[[source_stage]]
  if (is.null(k) && allow_null) return(NULL)
  if (!is.function(k)) stop(name, " list entries used by a source stage must be functions.")
  k
}


# Return the movement kernel for a progressing life stage.
.ipptm_get_transition_kernel <- function(TransitionKernels, source_stage, Nstages) {
  if (is.null(TransitionKernels)) return(NULL)
  if (!is.list(TransitionKernels) || length(TransitionKernels) != Nstages - 1L)
    stop("TransitionKernels must be NULL or a list of length Nstages - 1; NULL entries mean transition in place.")
  k <- TransitionKernels[[source_stage]]
  if (!is.null(k) && !is.function(k))
    stop("Every non-NULL TransitionKernels entry must be a dispersal-kernel function.")
  k
}


# Apply local weighted capacity after point establishment or movement.
.ipptm_capacity_keep <- function(candidates,
                                 existing_points,
                                 LocalK,
                                 KRadius,
                                 Weights,
                                 timestep,
                                 perm,
                                 Ntimesteps) {
  if (nrow(candidates) == 0L) return(logical(0))
  if (!is.function(LocalK) && length(LocalK) == 1L && is.infinite(LocalK))
    return(rep(TRUE, nrow(candidates)))
  if (KRadius <= 0) stop("Finite LocalK requires KRadius > 0.")

  k <- .ipp_resolve(LocalK, candidates, timestep, perm, Ntimesteps, "LocalK")
  if (any(k < 0)) stop("LocalK must be non-negative.")

  keep <- rep(FALSE, nrow(candidates))
  order_i <- sample(seq_len(nrow(candidates)))
  accepted <- existing_points

  for (ii in order_i) {
    d2 <- if (nrow(accepted))
      (accepted$x - candidates$x[ii])^2 + (accepted$y - candidates$y[ii])^2
    else numeric(0)
    local <- if (length(d2)) accepted[d2 <= KRadius^2, , drop = FALSE] else accepted[FALSE, , drop = FALSE]
    occupied <- if (nrow(local)) sum(Weights[local$stage]) else 0
    candidate_weight <- Weights[candidates$stage[ii]]
    if (occupied + candidate_weight <= k[ii] + sqrt(.Machine$double.eps)) {
      keep[ii] <- TRUE
      accepted <- rbind(accepted, candidates[ii, names(accepted), drop = FALSE])
    }
  }
  keep
}


# Build a standard point-transition event-log row set.
.ipptm_event <- function(perm,
                         timestep,
                         event,
                         point_id,
                         parent_id = NA_integer_,
                         x,
                         y,
                         stage_from = NA_integer_,
                         stage_to = NA_integer_,
                         detail = NA_character_,
                         provisional_x = NA_real_,
                         provisional_y = NA_real_,
                         habitat = NA_real_,
                         habitat_nudged = NA) {
  n <- length(point_id)
  data.frame(
    perm = rep_len(perm, n),
    timestep = rep_len(timestep, n),
    event = rep_len(event, n),
    point_id = point_id,
    parent_id = rep_len(parent_id, n),
    x = x,
    y = y,
    stage_from = rep_len(stage_from, n),
    stage_to = rep_len(stage_to, n),
    detail = rep_len(detail, n),
    provisional_x = rep_len(provisional_x, n),
    provisional_y = rep_len(provisional_y, n),
    habitat = rep_len(habitat, n),
    habitat_nudged = rep_len(habitat_nudged, n),
    stringsAsFactors = FALSE
  )
}


###############################################################################
### INApestPointTransitionMatrix
###############################################################################

INApestVertebratePoint <- function(

  ModelName = "INApestVertebratePoint",         # Model and output name
  Nperm,                                        # Number of stochastic simulation runs
  Ntimesteps,                                   # Timesteps in each simulation
  Nstages,                                      # Number of demographic/life stages
  Weights = rep(1, Nstages),                    # Stage weights used for population/capacity totals
  Transition,                                   # Stage-transition matrix or node/time-varying equivalent

  # One row per individual/nest; required x, y; optional stage.
  InitialPoints = NULL,                         # Fixed starting occurrences; NULL when a generator is used

  ###########################################################################
  ### Reproductive dispersal and establishment
  ###########################################################################

  SDDkernel,                                    # Short-distance displacement kernel
  LDDkernel = NULL,                             # Long-distance displacement kernel
  LDDrate = 0,                                  # Fraction of propagules entering long-distance dispersal
  PropaguleEstablishment = 1,                   # Establishment probability for arriving propagules
  EnvEstabProb = 1,                             # Environmental establishment probability

  ###########################################################################
  ### Stage progression movement
  ###########################################################################

  # List of length Nstages - 1. Entry s is used when stage s progresses
  # to s+1. NULL means progression occurs at the same coordinates.
  TransitionKernels = NULL,                     # Stage-specific movement kernels during progression

  # By default, transition movement can use habitat-seeking to choose its
  # final location but does not get a second establishment penalty: the
  # matrix progression probability already represents successful transition.
  TransitionHabitatSearch = FALSE,              # Allow moving transitions to seek suitable habitat
  ApplyHabitatToTransitions = FALSE,            # Apply habitat establishment penalty to stage movement
  TransitionEstablishment = 1,                  # Establishment probability after a moving stage transition
  BlockedTransitionMortality = 0,               # Mortality when stage progression is blocked

  ###########################################################################
  ### Habitat suitability and active habitat finding (native base R)
  ###########################################################################

  # Preferred input: INApestHabitatGrid() for one surface, or
  # INApestSpatialSchedule() to reuse static grids across timestep ranges.
  # NULL means suitability = 1 everywhere. A function(x, y, timestep) is
  # also supported; SpatRaster is retained only as optional compatibility.
  HabitatSuitability = NULL,                    # Spatial habitat/establishment suitability surface
  HabitatSearchRadius = 0,                      # Radius within which new points can seek better habitat
  HabitatSearchCandidates = 128,                # Candidate locations sampled during habitat search

  ###########################################################################
  ### Optional local weighted carrying capacity analogue
  ###########################################################################

  LocalK = Inf,                                 # Local weighted carrying-capacity limit
  KRadius = 0,                                  # Radius used to evaluate local carrying capacity

  ###########################################################################
  ### Detection and management
  ###########################################################################

  # Each *Spatial argument is an optional [0,1] multiplier using the same
  # INApestSpatialGrid architecture as habitat. INApestSpatialSchedule()
  # can reuse grids across timestep ranges. NULL means multiplier 1
  # everywhere. Spatial zero means that process is exactly absent locally.

  # Detection and mortality may be scalar, stage vector, timestep vector,
  # Nstages x Ntimesteps matrix, or function(points, timestep, perm).
  DetectionProb = 0,                            # Background-surveillance detection probability
  DetectionSD = NULL,                           # Variation in background detection probability
  DetectionSpatial = NULL,                      # Spatial multiplier for background detection

  ManageProb = 0,                               # Management probability when information is available
  ManageSD = NULL,                              # Variation in management probability
  ManageSpatial = NULL,                         # Spatial multiplier for management adoption

  MortalityProb = 0,                            # Management-driven host mortality probability
  MortalitySD = NULL,                           # Variation in management mortality probability
  MortalitySpatial = NULL,                      # Spatial multiplier for management mortality

  # Deliberately the same contract as INApestMetaPoint:
  # scalar or vector of length Ntimesteps, bounded [0,1].
  FecundityReduction = 0,                       # Proportional reduction in reproduction under management
  FecundityReductionSD = NULL,                  # Variation in fecundity reduction
  FecundityReductionSpatial = NULL,             # Spatial multiplier for fecundity reduction

  SpreadReduction = 0,                          # Proportional reduction in dispersal under management
  SpreadReductionSD = NULL,                     # Variation in spread reduction
  SpreadReductionSpatial = NULL,                # Spatial multiplier for spread reduction
  SpreadReductionAppliesTo = c("LDD", "all"),   # Whether spread reduction affects LDD only or all spread

  ###########################################################################
  ### Information spread / persistence
  ###########################################################################

  InitialInfo = NULL,                           # Starting information state
  InfoRadius = 0,                               # Distance over which information can transfer
  InfoTransferProb = 0,                         # Probability of information transfer between nearby points
  InfoKernel = NULL,                            # Optional distance kernel for information transfer
  InfoRetentionProb = 1,                        # Probability existing information persists one timestep
  InfoPersistenceSteps = NA,                    # Programmed timesteps information persists after evidence
  ExternalInfoProb = 0,                         # Information arriving from outside the modelled system
  OngoingExternalInfo = FALSE,                  # Allow new external information after initialisation

  ###########################################################################
  ### External incursion
  ###########################################################################

  ExternalIncursionGenerator = NULL,            # Function generating new external occurrence points
  OngoingExternalInvasion = FALSE,              # Allow new host incursions after initialisation

  ###########################################################################
  ### Vertebrate specialist modules
  ###########################################################################

  # One deliberately compact extension point. Supported named components:
  # Birth, HomeRange, Control and Interaction. NULL preserves the parent
  # INApestPointTransitionMatrix biological behaviour.
  Vertebrate = NULL,                            # Optional Birth/HomeRange/Control/Interaction modules

  ###########################################################################
  ### Output
  ###########################################################################

  OutputDir = NA,                               # Directory for saved outputs
  SaveResults = FALSE,                          # Save standard simulation outputs to disk
  DoProgress = TRUE,                            # Print simulation progress to the console
  Seed = NULL,                                  # Random seed for reproducible simulations
  InfoTriggeredDetectionProb = 0,               # Detection probability where information already exists
  InfoTriggeredDetectionSD = NULL,              # Variation in information-triggered detection
  InfoTriggeredDetectionSpatial = NULL,         # Spatial multiplier for information-triggered detection
  InitialPointGenerator = NULL                  # Optional function(perm) returning starting x/y points
) {


  # ---------------------------------------------------------------------------
  # Set up and validate the explicit point x life-stage simulation.
  # ---------------------------------------------------------------------------
  SpreadReductionAppliesTo <- match.arg(SpreadReductionAppliesTo)
  if (!exists(".iv_validate_vertebrate", mode = "function"))
    stop("INApestVertebrateHelpers.R must be sourced before INApestVertebratePoint().")
  .iv_validate_vertebrate(Vertebrate)
  BirthModule <- .iv_module(Vertebrate, "Birth")
  HomeRangeModule <- .iv_module(Vertebrate, "HomeRange")
  ControlModule <- .iv_module(Vertebrate, "Control")
  InteractionModule <- .iv_module(Vertebrate, "Interaction")
  if (!is.null(Seed)) set.seed(Seed)

  if (!is.numeric(Nstages) || length(Nstages) != 1L || Nstages < 2L || Nstages != floor(Nstages))
    stop("Nstages must be an integer >= 2.")
  Nstages <- as.integer(Nstages)

  if (length(Weights) != Nstages || any(!is.finite(Weights)) || any(Weights <= 0))
    stop("Weights must contain Nstages finite positive values.")
  Weights <- as.numeric(Weights)

  if (is.null(InitialPoints) && is.null(InitialPointGenerator))
    stop("Supply either InitialPoints or InitialPointGenerator.")
  if (!is.null(InitialPoints) && !is.null(InitialPointGenerator))
    stop("Supply InitialPoints or InitialPointGenerator, not both.")
  if (!is.null(InitialPointGenerator) && !is.function(InitialPointGenerator))
    stop("InitialPointGenerator must be NULL or a function(perm) returning a data.frame with x and y.")

  if (!is.null(InitialPoints)) {
    if (!is.data.frame(InitialPoints) || !all(c("x", "y") %in% names(InitialPoints)))
      stop("InitialPoints must be a data.frame containing x and y.")
    if (any(!is.finite(InitialPoints$x)) || any(!is.finite(InitialPoints$y)))
      stop("InitialPoints x and y must be finite.")
  }

  if (Nperm < 1 || Nperm != floor(Nperm)) stop("Nperm must be a positive integer.")
  if (Ntimesteps < 1 || Ntimesteps != floor(Ntimesteps)) stop("Ntimesteps must be a positive integer.")



  # Resolve once as an early validation, then re-resolve by timestep.
  invisible(.ipptm_get_transition(Transition, 1L, 1L, Nstages))

  # SDD/LDD can be a single function for every reproductive stage or a list
  # of Nstages stage-specific functions.
  if (!(is.function(SDDkernel) || (is.list(SDDkernel) && length(SDDkernel) == Nstages)))
    stop("SDDkernel must be a function or a list of length Nstages.")
  if (!is.null(LDDkernel) && !(is.function(LDDkernel) || (is.list(LDDkernel) && length(LDDkernel) == Nstages)))
    stop("LDDkernel must be NULL, a function, or a list of length Nstages.")
  if (!is.null(TransitionKernels) && (!is.list(TransitionKernels) || length(TransitionKernels) != Nstages - 1L))
    stop("TransitionKernels must be NULL or a list of length Nstages - 1.")

  .ipp_validate_probability_schedule(LDDrate, Ntimesteps, "LDDrate")
  .ipp_validate_probability_schedule(FecundityReduction, Ntimesteps, "FecundityReduction")
  .ipp_validate_probability_schedule(FecundityReductionSD, Ntimesteps, "FecundityReductionSD", allow_null = TRUE)

  .ipptm_validate_prob_like(DetectionProb, Ntimesteps, Nstages, "DetectionProb")
  .ipptm_validate_prob_like(DetectionSD, Ntimesteps, Nstages, "DetectionSD", allow_null = TRUE)
  .ipptm_validate_prob_like(InfoTriggeredDetectionProb, Ntimesteps, Nstages, "InfoTriggeredDetectionProb")
  .ipptm_validate_prob_like(InfoTriggeredDetectionSD, Ntimesteps, Nstages, "InfoTriggeredDetectionSD", allow_null = TRUE)
  .ipptm_validate_prob_like(MortalityProb, Ntimesteps, Nstages, "MortalityProb")
  .ipptm_validate_prob_like(MortalitySD, Ntimesteps, Nstages, "MortalitySD", allow_null = TRUE)

  btm <- if (length(BlockedTransitionMortality) == 1L) rep(BlockedTransitionMortality, Nstages - 1L) else BlockedTransitionMortality
  if (!is.numeric(btm) || length(btm) != Nstages - 1L || any(!is.finite(btm)) || any(btm < 0 | btm > 1))
    stop("BlockedTransitionMortality must be a scalar or length Nstages - 1 vector in [0,1].")

  if (HabitatSearchRadius < 0 || KRadius < 0 || InfoRadius < 0)
    stop("HabitatSearchRadius, KRadius and InfoRadius must be >= 0.")


  # Validate spatial surfaces and, for time-varying grids, layer availability.
  .ipp_validate_spatial_surface(HabitatSuitability, Ntimesteps, "HabitatSuitability")
  .ipp_validate_spatial_surface(DetectionSpatial, Ntimesteps, "DetectionSpatial")
  .ipp_validate_spatial_surface(InfoTriggeredDetectionSpatial, Ntimesteps, "InfoTriggeredDetectionSpatial")
  .ipp_validate_spatial_surface(ManageSpatial, Ntimesteps, "ManageSpatial")
  .ipp_validate_spatial_surface(MortalitySpatial, Ntimesteps, "MortalitySpatial")
  .ipp_validate_spatial_surface(FecundityReductionSpatial, Ntimesteps, "FecundityReductionSpatial")
  .ipp_validate_spatial_surface(SpreadReductionSpatial, Ntimesteps, "SpreadReductionSpatial")

  if (OngoingExternalInvasion && !is.function(ExternalIncursionGenerator))
    stop("OngoingExternalInvasion = TRUE requires ExternalIncursionGenerator(timestep, perm).")

  # INApestMeta-style stochastic defaults, with zero effects remaining exactly off.
  if (is.null(ManageSD) && !is.function(ManageProb)) ManageSD <- ManageProb / 10
  if (is.null(DetectionSD) && !is.function(DetectionProb)) DetectionSD <- DetectionProb / 10
  if (is.null(InfoTriggeredDetectionSD) && !is.function(InfoTriggeredDetectionProb)) InfoTriggeredDetectionSD <- InfoTriggeredDetectionProb / 10
  if (is.null(InfoTriggeredDetectionSD) && is.function(InfoTriggeredDetectionProb)) InfoTriggeredDetectionSD <- 0
  UseInfoTriggeredSurveillance <- is.function(InfoTriggeredDetectionProb) || is.function(InfoTriggeredDetectionSD) ||
    any(as.numeric(InfoTriggeredDetectionProb) != 0, na.rm = TRUE) ||
    any(as.numeric(InfoTriggeredDetectionSD) != 0, na.rm = TRUE)
  if (is.null(MortalitySD) && !is.function(MortalityProb)) MortalitySD <- MortalityProb / 10
  if (is.null(FecundityReductionSD)) FecundityReductionSD <- FecundityReduction / 10
  if (is.null(SpreadReductionSD) && !is.function(SpreadReduction)) SpreadReductionSD <- (1 - SpreadReduction) / 10

  .ipp_validate_probability_schedule(FecundityReductionSD, Ntimesteps, "FecundityReductionSD")

  snapshots <- list(); events <- list(); final_points <- list(); info_results <- list(); summaries <- list()
  contacts <- list(); interaction_events <- list(); control_observations <- list()
  si <- 0L; ei <- 0L; sumi <- 0L; ci <- 0L; iei <- 0L; coi <- 0L


  # ---------------------------------------------------------------------------
  # Run independent stochastic point-transition histories.
  # ---------------------------------------------------------------------------
  for (perm in seq_len(Nperm)) {

    # Draw uncertain initial locations for THIS permutation only.
    # The realised points are then passed through the established initial
    # point-state pathway, so all existing stage, information, pathogen
    # and vertebrate behaviour is retained.
    if (!is.null(InitialPointGenerator)) {
      InitialPoints <- InitialPointGenerator(perm = perm)
      if (!is.data.frame(InitialPoints) || !all(c("x", "y") %in% names(InitialPoints)))
        stop("InitialPointGenerator must return a data.frame containing x and y.")
      if (any(!is.finite(InitialPoints$x)) || any(!is.finite(InitialPoints$y)))
        stop("InitialPointGenerator output x and y must be finite.")
    }

    initial_stage <- if ("stage" %in% names(InitialPoints)) {
      as.integer(InitialPoints$stage)
    } else {
      rep(1L, nrow(InitialPoints))
    }
    if (any(is.na(initial_stage)) || any(initial_stage < 1L | initial_stage > Nstages))
      stop("InitialPoints$stage must contain integers from 1 to Nstages.")
    n_initial <- nrow(InitialPoints)
    # Preserve every user-supplied individual attribute. The engine owns only
    # the reserved INApest fields; sex, age, behaviour, group membership,
    # infection state, genotype and other atomic columns remain attached.
    points <- .iv_initial_point_state(InitialPoints, initial_stage)
    next_id <- nrow(points) + 1L

    if (!is.null(InitialInfo)) {
      if (is.logical(InitialInfo) && length(InitialInfo) == nrow(points)) {
        points$have_info <- InitialInfo
      } else if (is.numeric(InitialInfo)) {
        idx <- intersect(as.integer(InitialInfo), seq_len(nrow(points)))
        points$have_info[idx] <- TRUE
      } else stop("InitialInfo must be NULL, a logical vector per initial point, or row indices.")
    } else if ("have_info" %in% names(InitialPoints)) {
      points$have_info <- as.logical(InitialPoints$have_info)
      points$have_info[is.na(points$have_info)] <- FALSE
    }

# InitialInfo / InitialPoints$have_info represents information known at model time 0.
if (nrow(points) && any(points$have_info)) {
  points$last_known_timestep[points$have_info] <- 0L
}
info_sites <- .ipp_empty_info_sites()

    # Initial detection uses stage-specific per-individual probabilities.
    if (nrow(points)) {
      pdet0 <- .ipptm_probability_spatial(DetectionProb, DetectionSD, DetectionSpatial, points, 1L, perm, Ntimesteps, Nstages, "DetectionProb")
      d0 <- stats::rbinom(nrow(points), 1L, pdet0) == 1L
      points$detected[d0] <- TRUE
      points$have_info[d0] <- TRUE
      points$last_known_timestep[d0] <- 0L
      if (any(d0)) info_sites <- .ipp_add_or_refresh_info_sites(info_sites, points[d0, c("id", "x", "y"), drop = FALSE], 0L)
    }


    # ---------------------------------------------------------------------------
    # Advance stage transitions, movement, information and response through time.
    # ---------------------------------------------------------------------------
    for (timestep in seq_len(Ntimesteps)) {
      if (DoProgress) cat("\r", "Realisation", perm, "Timestep", timestep, "...")

      A <- .ipptm_get_transition(Transition, timestep, perm, Nstages)
      fecundity <- c(0, A[1, 2:Nstages])
      ldd_rate_t <- if (length(LDDrate) == 1L) LDDrate else LDDrate[timestep]

      n_start <- nrow(points)
      stage_start <- tabulate(points$stage, nbins = Nstages)
      n_management_deaths <- 0L
      n_stage_deaths <- 0L
      n_stage_transitions <- 0L
      n_propagules <- 0L
      n_established <- 0L
      n_external <- 0L
      n_new_detections <- 0L
      n_new_background_detections <- 0L
      n_new_info_triggered_detections <- 0L
      n_info_before_surveillance <- 0L
      InfoBeforeSurveillancePoint <- logical(0)
      n_managing <- 0L
      n_control_deaths <- 0L
      n_control_detections <- 0L
      control_cost <- 0
      control_fecundity_reduction <- setNames(numeric(nrow(points)), as.character(points$id))
      pending_control_detection <- data.frame(
        id = integer(0), parent_id = integer(0), x = numeric(0), y = numeric(0),
        stage = integer(0), stringsAsFactors = FALSE
      )

      #########################################################################
      ### 0. Pre-deployed vertebrate control (independent of pest information)
      #########################################################################
      if (nrow(points) && !is.null(ControlModule)) {
        control_context <- list(
          transition = A, ModelName = ModelName, Nstages = Nstages,
          Weights = Weights, phase = "pre_response_control"
        )
        hr_control <- .iv_home_range_point(
          HomeRangeModule, points, timestep, perm, control_context
        )
        ce <- .iv_point_control(
          ControlModule, points, hr_control, timestep, perm, control_context
        )
        control_cost <- ce$cost
        control_fecundity_reduction <- setNames(
          ce$effects$fecundity_reduction, as.character(points$id)
        )

        # Retain the exact pre-control detection opportunity for PoA. The
        # detection and lethal effects below share one uniform draw, so the
        # companion also needs the realised death outcome to evaluate the
        # joint no-detection likelihood without changing biological semantics.
        control_obs <- points[, c("id","parent_id","x","y","stage"), drop=FALSE]
        control_obs$control_detection_prob <- ce$effects$detect_prob
        control_obs$control_kill_prob <- ce$effects$kill_prob

        # A shared uniform couples detection and lethal effects sensibly for
        # simple devices (e.g. kill traps with detection on encounter). Custom
        # Control$Model functions can alter the marginal probabilities.
        u_control <- stats::runif(nrow(points))
        control_detected <- u_control < ce$effects$detect_prob
        control_dead <- u_control < ce$effects$kill_prob
        n_control_detections <- sum(control_detected)
        n_control_deaths <- sum(control_dead)
        control_obs$control_dead <- as.logical(control_dead)
        control_obs$perm <- perm
        control_obs$timestep <- timestep
        coi <- coi + 1L
        control_observations[[coi]] <- control_obs[, c("perm","timestep","id","parent_id","x","y","stage",
          "control_detection_prob","control_kill_prob","control_dead"), drop=FALSE]

        if (any(control_detected)) {
          pending_control_detection <- points[
            control_detected, c("id", "parent_id", "x", "y", "stage"), drop = FALSE
          ]
        }
        if (any(control_dead)) {
          ei <- ei + 1L
          events[[ei]] <- .ipptm_event(
            perm, timestep, "death", points$id[control_dead],
            points$parent_id[control_dead], points$x[control_dead],
            points$y[control_dead], points$stage[control_dead], NA_integer_,
            "vertebrate_predeployed_control"
          )
          points <- points[!control_dead, , drop = FALSE]
          control_fecundity_reduction <- control_fecundity_reduction[as.character(points$id)]
        }
      }

      #########################################################################
      ### 1. Information-triggered response management and mortality
      #########################################################################
      if (nrow(points)) {
        p_manage <- .ipp_probability_spatial(ManageProb, ManageSD, ManageSpatial, points, timestep, perm, Ntimesteps, "ManageProb")
        points$managing <- stats::rbinom(nrow(points), 1L, p_manage * as.numeric(points$have_info)) == 1L
        n_managing <- sum(points$managing)

        p_mort <- .ipptm_probability_spatial(MortalityProb, MortalitySD, MortalitySpatial, points, timestep, perm, Ntimesteps, Nstages, "MortalityProb")
        management_dead <- points$managing & (stats::runif(nrow(points)) < p_mort)
        n_management_deaths <- sum(management_dead)

        if (any(management_dead)) {
          ei <- ei + 1L
          events[[ei]] <- .ipptm_event(perm, timestep, "death", points$id[management_dead], points$parent_id[management_dead], points$x[management_dead], points$y[management_dead], points$stage[management_dead], NA_integer_, "management_mortality")
          info_sites <- .ipp_add_or_refresh_info_sites(info_sites, points[management_dead, c("id", "x", "y"), drop = FALSE], timestep)
          points <- points[!management_dead, , drop = FALSE]
        }
      }

      #########################################################################
      ### 2. Births: default matrix fecundity or a vertebrate Birth hook
      #########################################################################
      propagule_parents <- points[FALSE, , drop = FALSE]
      propagule_state <- data.frame(parent_id = integer(0), stage = integer(0))
      if (nrow(points)) {
        f_red <- .ipp_probability_spatial(
          FecundityReduction, FecundityReductionSD, FecundityReductionSpatial,
          points, timestep, perm, Ntimesteps, "FecundityReduction"
        )
        f_red <- .ipp_clip01(f_red)
        ctrl_fr <- as.numeric(control_fecundity_reduction[as.character(points$id)])
        ctrl_fr[is.na(ctrl_fr)] <- 0
        fec_mult <- (1 - f_red * as.numeric(points$managing)) * (1 - ctrl_fr)

        if (is.null(BirthModule)) {
          lambda <- pmax(0, fecundity[points$stage] * fec_mult)
          n_by_parent <- stats::rpois(nrow(points), lambda)
          n_propagules <- sum(n_by_parent)
          if (n_propagules > 0L) {
            parent_rows <- rep(seq_len(nrow(points)), n_by_parent)
            propagule_parents <- points[parent_rows, , drop = FALSE]
            propagule_state <- data.frame(
              parent_id = points$id[parent_rows], stage = 1L,
              stringsAsFactors = FALSE
            )
          }
        } else {
          birth_context <- list(
            fecundity = fecundity, fecundity_multiplier = fec_mult,
            response_fecundity_reduction = f_red,
            control_fecundity_reduction = ctrl_fr,
            managing = points$managing, home_range = .iv_home_range_point(
              HomeRangeModule, points, timestep, perm,
              list(transition = A, phase = "birth")
            ),
            ModelName = ModelName, Nstages = Nstages, Weights = Weights
          )
          # A custom Birth hook returns realised offspring (or integer counts).
          # It receives fecundity_multiplier and is responsible for applying it
          # to its own breeding logic.
          propagule_state <- .iv_point_birth(
            BirthModule, points, A, timestep, perm, birth_context
          )
          if (nrow(propagule_state)) {
            if (any(is.na(propagule_state$stage)) ||
                any(propagule_state$stage < 1L | propagule_state$stage > Nstages))
              stop("Vertebrate$Birth offspring stage must be 1..Nstages.")
            propagule_parents <- points[
              match(propagule_state$parent_id, points$id), , drop = FALSE
            ]
          }
          n_propagules <- nrow(propagule_state)
        }
      }

      #########################################################################
      ### 3. Stage survival/stasis/progression, with optional relocation kernel
      #########################################################################
      if (nrow(points)) {
        old <- points
        alive <- rep(TRUE, nrow(points))
        progress <- rep(FALSE, nrow(points))

        for (s in seq_len(Nstages)) {
          idx <- which(old$stage == s)
          if (!length(idx)) next
          if (s == Nstages) {
            surv <- stats::rbinom(length(idx), 1L, A[s, s]) == 1L
            alive[idx] <- surv
          } else {
            p_stay <- A[s, s]
            p_prog <- A[s + 1L, s]
            p_surv <- pmin(1, p_stay + p_prog)
            surv <- stats::rbinom(length(idx), 1L, p_surv) == 1L
            alive[idx] <- surv
            if (p_surv > 0 && any(surv)) {
              prog_cond <- p_prog / p_surv
              progress[idx[surv]] <- stats::rbinom(sum(surv), 1L, prog_cond) == 1L
            }
          }
        }

        stage_dead <- !alive
        n_stage_deaths <- sum(stage_dead)
        if (any(stage_dead)) {
          ei <- ei + 1L
          events[[ei]] <- .ipptm_event(perm, timestep, "death", old$id[stage_dead], old$parent_id[stage_dead], old$x[stage_dead], old$y[stage_dead], old$stage[stage_dead], NA_integer_, "transition_matrix_mortality")
        }

        points <- old[alive, , drop = FALSE]

        # Use immutable point IDs rather than row positions for progression.
        # This keeps multi-stage transitions correct even if failed transitions
        # remove some points while another source stage is still being processed.
        if (any(progress & alive)) {
          for (s in seq_len(Nstages - 1L)) {
            source_ids <- old$id[alive & progress & old$stage == s]
            if (!length(source_ids)) next

            ii <- match(source_ids, points$id)
            keep_present <- !is.na(ii)
            source_ids <- source_ids[keep_present]
            ii <- ii[keep_present]
            if (!length(ii)) next

            k <- .ipptm_get_transition_kernel(TransitionKernels, s, Nstages)
            provisional_x <- points$x[ii]
            provisional_y <- points$y[ii]
            if (!is.null(k)) {
              z <- .ipp_draw_displacement(
                k, points[ii, , drop = FALSE], timestep, perm,
                paste0("TransitionKernels[[", s, "]]")
              )
              provisional_x <- provisional_x + z$dx
              provisional_y <- provisional_y + z$dy
            }

            if (TransitionHabitatSearch && !is.null(HabitatSuitability) && HabitatSearchRadius > 0) {
              dest <- .ipp_habitat_search(
                provisional_x, provisional_y, HabitatSuitability,
                HabitatSearchRadius, timestep, HabitatSearchCandidates
              )
            } else {
              hv <- .ipp_habitat_value(
                provisional_x, provisional_y, HabitatSuitability, timestep
              )
              dest <- data.frame(
                x = provisional_x, y = provisional_y,
                habitat = hv, habitat_nudged = FALSE
              )
            }

            success <- rep(TRUE, length(source_ids))
            if (ApplyHabitatToTransitions) {
              cand <- data.frame(x = dest$x, y = dest$y, stage = s + 1L)
              pte <- .ipp_resolve(
                TransitionEstablishment, cand, timestep, perm,
                Ntimesteps, "TransitionEstablishment"
              )
              pee <- .ipp_resolve(
                EnvEstabProb, cand, timestep, perm,
                Ntimesteps, "EnvEstabProb"
              )
              p_success <- .ipp_clip01(pte * pee * dest$habitat)
              success <- stats::rbinom(length(source_ids), 1L, p_success) == 1L
            }

            failed <- !success
            if (any(failed)) {
              failed_ids <- source_ids[failed]
              die_blocked <- stats::rbinom(
                length(failed_ids), 1L, btm[s]
              ) == 1L

              if (any(die_blocked)) {
                doomed_ids <- failed_ids[die_blocked]
                doomed_rows <- match(doomed_ids, points$id)
                doomed_rows <- doomed_rows[!is.na(doomed_rows)]

                if (length(doomed_rows)) {
                  ei <- ei + 1L
                  events[[ei]] <- .ipptm_event(
                    perm, timestep, "death",
                    points$id[doomed_rows], points$parent_id[doomed_rows],
                    points$x[doomed_rows], points$y[doomed_rows],
                    s, NA_integer_, "blocked_transition_mortality"
                  )
                  points <- points[-doomed_rows, , drop = FALSE]
                }
              }
            }

            success_ids <- source_ids[success]
            jj <- match(success_ids, points$id)
            present <- !is.na(jj)
            success_ids <- success_ids[present]
            jj <- jj[present]

            if (length(jj)) {
              m <- match(success_ids, source_ids)
              stage_from <- points$stage[jj]
              points$stage[jj] <- s + 1L
              points$x[jj] <- dest$x[m]
              points$y[jj] <- dest$y[m]
              n_stage_transitions <- n_stage_transitions + length(jj)

              ei <- ei + 1L
              events[[ei]] <- .ipptm_event(
                perm, timestep, "stage_transition",
                points$id[jj], points$parent_id[jj],
                points$x[jj], points$y[jj],
                stage_from, s + 1L,
                if (is.null(k)) "in_place" else "transition_kernel",
                provisional_x[m], provisional_y[m],
                dest$habitat[m], dest$habitat_nudged[m]
              )
            }
          }
        }
      }

      #########################################################################
      ### 4. Reproductive propagule SDD/LDD, habitat search and establishment
      #########################################################################
      if (nrow(propagule_parents)) {
        ps <- propagule_state
        is_ldd <- stats::rbinom(nrow(propagule_parents), 1L, ldd_rate_t) == 1L
        spread_red <- .ipp_probability_spatial(SpreadReduction, SpreadReductionSD, SpreadReductionSpatial, propagule_parents, timestep, perm, Ntimesteps, "SpreadReduction")
        suppressed <- if (SpreadReductionAppliesTo == "LDD")
          is_ldd & propagule_parents$managing & (stats::runif(nrow(propagule_parents)) < spread_red)
        else
          propagule_parents$managing & (stats::runif(nrow(propagule_parents)) < spread_red)

        pp <- propagule_parents[!suppressed, , drop = FALSE]
        ps <- ps[!suppressed, , drop = FALSE]
        is_ldd <- is_ldd[!suppressed]

        if (nrow(pp)) {
          dx <- dy <- numeric(nrow(pp))
          for (s in seq_len(Nstages)) {
            sidx <- which(pp$stage == s & !is_ldd)
            if (length(sidx)) {
              k <- .ipptm_get_kernel(SDDkernel, s, Nstages, "SDDkernel")
              z <- .ipp_draw_displacement(k, pp[sidx, , drop = FALSE], timestep, perm, "SDDkernel")
              dx[sidx] <- z$dx; dy[sidx] <- z$dy
            }
            lidx <- which(pp$stage == s & is_ldd)
            if (length(lidx)) {
              k <- .ipptm_get_kernel(LDDkernel, s, Nstages, "LDDkernel")
              z <- .ipp_draw_displacement(k, pp[lidx, , drop = FALSE], timestep, perm, "LDDkernel")
              dx[lidx] <- z$dx; dy[lidx] <- z$dy
            }
          }

          px <- pp$x + dx; py <- pp$y + dy
          dest <- .ipp_habitat_search(px, py, HabitatSuitability, HabitatSearchRadius, timestep, HabitatSearchCandidates)
          cand <- data.frame(
            x = dest$x, y = dest$y, stage = as.integer(ps$stage), parent_id = pp$id,
            dispersal_type = ifelse(is_ldd, "LDD", "SDD"),
            provisional_x = px, provisional_y = py,
            habitat = dest$habitat, habitat_nudged = dest$habitat_nudged,
            stringsAsFactors = FALSE
          )

          pestab <- .ipp_resolve(PropaguleEstablishment, cand, timestep, perm, Ntimesteps, "PropaguleEstablishment")
          penv <- .ipp_resolve(EnvEstabProb, cand, timestep, perm, Ntimesteps, "EnvEstabProb")
          p_est <- .ipp_clip01(pestab * penv * cand$habitat)
          established <- stats::rbinom(nrow(cand), 1L, p_est) == 1L
          cand <- cand[established, , drop = FALSE]
          ps <- ps[established, , drop = FALSE]

          if (nrow(cand)) {
            capacity_keep <- .ipptm_capacity_keep(cand, points, LocalK, KRadius, Weights, timestep, perm, Ntimesteps)
            cand <- cand[capacity_keep, , drop = FALSE]
            ps <- ps[capacity_keep, , drop = FALSE]
          }

          if (nrow(cand)) {
            n_established <- nrow(cand)
            recruits <- data.frame(
              x = cand$x, y = cand$y, stage = as.integer(ps$stage),
              id = seq.int(next_id, length.out = n_established),
              parent_id = cand$parent_id,
              birth_timestep = timestep,
              have_info = FALSE, detected = FALSE, managing = FALSE,
              last_known_timestep = NA_integer_, stringsAsFactors = FALSE
            )
            extra_birth_cols <- setdiff(names(ps), .iv_point_reserved)
            extra_birth_cols <- setdiff(extra_birth_cols, c("parent_id", "stage"))
            if (length(extra_birth_cols))
              recruits <- cbind(recruits, ps[, extra_birth_cols, drop = FALSE])
            next_id <- next_id + n_established
            ei <- ei + 1L
            events[[ei]] <- .ipptm_event(
              perm, timestep, "establishment", recruits$id, recruits$parent_id,
              recruits$x, recruits$y, pp$stage[match(recruits$parent_id, pp$id)], recruits$stage,
              cand$dispersal_type, cand$provisional_x, cand$provisional_y,
              cand$habitat, cand$habitat_nudged
            )
            points <- .iv_rbind_fill(points, recruits)
          }
        }
      }

      #########################################################################
      ### 5. External incursions
      #########################################################################
      if (OngoingExternalInvasion) {
        ext <- ExternalIncursionGenerator(timestep = timestep, perm = perm)
        if (!is.null(ext) && nrow(ext)) {
          if (!is.data.frame(ext) || !all(c("x", "y") %in% names(ext)))
            stop("ExternalIncursionGenerator must return NULL/zero rows or a data.frame with x and y.")
          if (!"stage" %in% names(ext)) ext$stage <- 1L
          if (any(ext$stage < 1L | ext$stage > Nstages)) stop("External incursion stage must be 1..Nstages.")
          n_external <- nrow(ext)
          ep <- .iv_external_point_state(ext, next_id, timestep, Nstages)
          next_id <- next_id + n_external
          points <- .iv_rbind_fill(points, ep)
          ei <- ei + 1L
          events[[ei]] <- .ipptm_event(perm, timestep, "external_incursion", ep$id, NA_integer_, ep$x, ep$y, NA_integer_, ep$stage, "external")
        }
      }

      #########################################################################
      ### 5b. Social/contact interaction hook
      #########################################################################
      if (nrow(points) && !is.null(InteractionModule)) {
        interaction_context <- list(
          transition = A, ModelName = ModelName, Nstages = Nstages,
          Weights = Weights, phase = "interaction"
        )
        hr_interaction <- .iv_home_range_point(
          HomeRangeModule, points, timestep, perm, interaction_context
        )
        iz <- .iv_point_interaction(
          InteractionModule, points, hr_interaction, timestep, perm,
          interaction_context
        )
        points <- iz$points
        # Generic point-pathogen adapters mark pathogen-associated deaths inside
        # Interaction so the identity-preserving hook can complete safely. The
        # parent engine owns actual point removal immediately afterwards.
        if (exists("INApestPointPathogenApplyDeaths", mode = "function") &&
            ".pathogen_death" %in% names(points)) {
          n_pathogen_deaths <- sum(points$.pathogen_death %in% TRUE, na.rm = TRUE)
          points <- INApestPointPathogenApplyDeaths(points)
          if (n_pathogen_deaths > 0L) {
            iei <- iei + 1L
            interaction_events[[iei]] <- data.frame(
              perm = perm, timestep = timestep,
              event = "pathogen_mortality", n = n_pathogen_deaths,
              stringsAsFactors = FALSE
            )
          }
        }
        if (nrow(iz$contacts)) {
          ci <- ci + 1L
          cc <- iz$contacts; cc$perm <- perm; cc$timestep <- timestep
          contacts[[ci]] <- cc[, c("perm", "timestep", setdiff(names(cc), c("perm", "timestep"))), drop = FALSE]
        }
        if (!is.null(iz$events)) {
          ee <- as.data.frame(iz$events, stringsAsFactors = FALSE)
          if (nrow(ee) > 0L) {
            iei <- iei + 1L
            ee$perm <- rep(perm, nrow(ee)); ee$timestep <- rep(timestep, nrow(ee))
            interaction_events[[iei]] <- ee[, c("perm", "timestep", setdiff(names(ee), c("perm", "timestep"))), drop = FALSE]
          }
        }
      }

      #########################################################################
      ### 6. Information persistence, transfer and external information
      #########################################################################
      if (nrow(points)) {
        informed <- which(points$have_info)
        if (length(informed)) {
          persistence <- .ipp_resolve(InfoPersistenceSteps, points[informed, , drop = FALSE], timestep, perm, Ntimesteps, "InfoPersistenceSteps")
          programmed <- !is.na(persistence)
          if (any(programmed)) {
            last_known <- points$last_known_timestep[informed]
            age <- timestep - last_known
            stop_info <- programmed & (is.na(last_known) | age >= persistence)
            if (any(stop_info)) points$have_info[informed[stop_info]] <- FALSE
          }
          retention_group <- informed[is.na(persistence) & points$have_info[informed]]
          if (length(retention_group)) {
            pret <- .ipp_resolve(InfoRetentionProb, points[retention_group, , drop = FALSE], timestep, perm, Ntimesteps, "InfoRetentionProb")
            retain <- stats::rbinom(length(retention_group), 1L, .ipp_clip01(pret)) == 1L
            points$have_info[retention_group[!retain]] <- FALSE
          }
        }
      }

      if (nrow(info_sites)) {
        active <- which(info_sites$active)
        if (length(active)) {
          site_points <- data.frame(x = info_sites$x[active], y = info_sites$y[active])
          persistence <- .ipp_resolve(InfoPersistenceSteps, site_points, timestep, perm, Ntimesteps, "InfoPersistenceSteps")
          programmed <- !is.na(persistence)
          if (any(programmed)) {
            age <- timestep - info_sites$last_known_timestep[active]
            stop_sites <- active[programmed & age >= persistence]
            if (length(stop_sites)) info_sites$active[stop_sites] <- FALSE
          }
          retention_group <- active[is.na(persistence) & info_sites$active[active]]
          if (length(retention_group)) {
            sp <- data.frame(x = info_sites$x[retention_group], y = info_sites$y[retention_group])
            pret <- .ipp_resolve(InfoRetentionProb, sp, timestep, perm, Ntimesteps, "InfoRetentionProb")
            retain <- stats::rbinom(length(retention_group), 1L, .ipp_clip01(pret)) == 1L
            info_sites$active[retention_group[!retain]] <- FALSE
          }
        }
      }

      if (nrow(points)) {
        points <- .ipp_transfer_information(points, info_sites, timestep, InfoRadius, InfoTransferProb, InfoKernel, perm, Ntimesteps)
        if (OngoingExternalInfo) {
          pe <- .ipp_resolve(ExternalInfoProb, points, timestep, perm, Ntimesteps, "ExternalInfoProb")
          extinf <- stats::rbinom(nrow(points), 1L, .ipp_clip01(pe)) == 1L
          points$have_info[extinf] <- TRUE
        }
      }

      #########################################################################
      ### Freeze information available before this round's surveillance.
      # Same-round control/background detections must not retrospectively
      ### activate information-triggered surveillance.
      #########################################################################
      if (nrow(points)) {
        InfoBeforeSurveillancePoint <- as.logical(points$have_info)
        n_info_before_surveillance <- sum(InfoBeforeSurveillancePoint)
      }

      #########################################################################
      ### 6b. Register detections generated by pre-deployed control
      # Detection is deliberately registered here, after information transfer,
      ### so it can trigger response management only from the next timestep.
      #########################################################################
      if (nrow(pending_control_detection)) {
        info_sites <- .ipp_add_or_refresh_info_sites(
          info_sites, pending_control_detection[, c("id", "x", "y"), drop = FALSE],
          timestep
        )
        jj <- match(pending_control_detection$id, points$id)
        surviving <- !is.na(jj)
        if (any(surviving)) {
          kk <- jj[surviving]
          points$detected[kk] <- TRUE
          points$have_info[kk] <- TRUE
          points$last_known_timestep[kk] <- timestep
        }
        ei <- ei + 1L
        events[[ei]] <- .ipptm_event(
          perm, timestep, "detection", pending_control_detection$id,
          pending_control_detection$parent_id, pending_control_detection$x,
          pending_control_detection$y, pending_control_detection$stage,
          pending_control_detection$stage, "vertebrate_control_detection"
        )
      }

      #########################################################################
      ### 7. Detection after population dynamics
      #########################################################################
      if (nrow(points)) {
        pdet <- .ipptm_probability_spatial(DetectionProb, DetectionSD, DetectionSpatial, points, timestep, perm, Ntimesteps, Nstages, "DetectionProb")
        background_detected_now <- stats::rbinom(nrow(points), 1L, pdet) == 1L
        pdet_info <- rep(0, nrow(points))
        info_triggered_detected_now <- rep(FALSE, nrow(points))
        if (UseInfoTriggeredSurveillance) {
          pdet_info <- .ipptm_probability_spatial(
            InfoTriggeredDetectionProb, InfoTriggeredDetectionSD, InfoTriggeredDetectionSpatial,
            points, timestep, perm, Ntimesteps, Nstages, "InfoTriggeredDetectionProb"
          )
          info_triggered_detected_now <- stats::rbinom(
            nrow(points), 1L, pdet_info * as.numeric(InfoBeforeSurveillancePoint)
          ) == 1L
        }
        detected_now <- background_detected_now | info_triggered_detected_now
        new_detection <- detected_now & !points$detected
        n_new_background_detections <- sum(background_detected_now)
        n_new_info_triggered_detections <- sum(info_triggered_detected_now)
        n_new_detections <- sum(new_detection)
        if (any(background_detected_now)) {
          ei <- ei + 1L
          events[[ei]] <- .ipptm_event(perm, timestep, "background_detection", points$id[background_detected_now], points$parent_id[background_detected_now], points$x[background_detected_now], points$y[background_detected_now], points$stage[background_detected_now], points$stage[background_detected_now], "background")
        }
        if (any(info_triggered_detected_now)) {
          ei <- ei + 1L
          events[[ei]] <- .ipptm_event(perm, timestep, "info_triggered_detection", points$id[info_triggered_detected_now], points$parent_id[info_triggered_detected_now], points$x[info_triggered_detected_now], points$y[info_triggered_detected_now], points$stage[info_triggered_detected_now], points$stage[info_triggered_detected_now], "info_triggered")
        }
        points$detected[detected_now] <- TRUE
        points$have_info[detected_now] <- TRUE
        points$last_known_timestep[detected_now] <- timestep
        if (any(detected_now)) info_sites <- .ipp_add_or_refresh_info_sites(info_sites, points[detected_now, c("id", "x", "y"), drop = FALSE], timestep)
        if (any(new_detection)) {
          ei <- ei + 1L
          events[[ei]] <- .ipptm_event(perm, timestep, "detection", points$id[new_detection], points$parent_id[new_detection], points$x[new_detection], points$y[new_detection], points$stage[new_detection], points$stage[new_detection], "detected")
        }
      }

      #########################################################################
      ### 8. Store outputs
      #########################################################################
      if (nrow(points)) {
        si <- si + 1L
        snap <- points
        snap$have_info_before_surveillance <- InfoBeforeSurveillancePoint
        snap$background_detection_prob <- pdet
        snap$info_triggered_detection_prob <- pdet_info
        snap$perm <- perm; snap$timestep <- timestep
        front <- c("perm", "timestep", .iv_point_reserved)
        snapshots[[si]] <- snap[, c(front, setdiff(names(snap), front)), drop = FALSE]
      }

      stage_end <- tabulate(points$stage, nbins = Nstages)
      sumi <- sumi + 1L
      row <- data.frame(
        perm = perm, timestep = timestep, n_start = n_start,
        n_control_deaths = n_control_deaths,
        n_control_detections = n_control_detections,
        control_cost = control_cost,
        n_management_deaths = n_management_deaths,
        n_stage_deaths = n_stage_deaths,
        n_stage_transitions = n_stage_transitions,
        n_propagules = n_propagules, n_established = n_established,
        n_external = n_external, n_end = nrow(points),
        weighted_population = if (nrow(points)) sum(Weights[points$stage]) else 0,
        n_detected = if (nrow(points)) sum(points$detected) else 0L,
        n_new_detections = n_new_detections,
        n_new_background_detections = n_new_background_detections,
        n_new_info_triggered_detections = n_new_info_triggered_detections,
        n_info_before_surveillance = n_info_before_surveillance,
        n_have_info = if (nrow(points)) sum(points$have_info) else 0L,
        n_managing = n_managing,
        stringsAsFactors = FALSE
      )
      for (s in seq_len(Nstages)) {
        row[[paste0("stage", s, "_start")]] <- stage_start[s]
        row[[paste0("stage", s, "_end")]] <- stage_end[s]
      }
      summaries[[sumi]] <- row
    }

    fp <- points; fp$perm <- rep(perm, nrow(fp)); final_points[[perm]] <- fp
    ip <- info_sites; ip$perm <- rep(perm, nrow(ip)); info_results[[perm]] <- ip
  }

  if (DoProgress) cat("\n")

  out <- list(
    ModelName = ModelName,
    PointHistory = if (length(snapshots)) .iv_rbind_fill(snapshots) else data.frame(),
    EventLog = if (length(events)) do.call(rbind, events) else data.frame(),
    FinalPoints = if (length(final_points)) .iv_rbind_fill(final_points) else data.frame(),
    InfoSites = if (length(info_results)) .iv_rbind_fill(info_results) else data.frame(),
    ContactHistory = if (length(contacts)) .iv_rbind_fill(contacts) else data.frame(),
    InteractionEvents = if (length(interaction_events)) .iv_rbind_fill(interaction_events) else data.frame(),
    ControlObservationHistory = if (length(control_observations)) .iv_rbind_fill(control_observations) else data.frame(),
    Summary = if (length(summaries)) do.call(rbind, summaries) else data.frame()
  )
  class(out) <- c("INApestVertebratePoint", "INApestPointTransitionMatrix", "list")

  if (SaveResults) {
    if (is.na(OutputDir)) OutputDir <- ""
    if (nzchar(OutputDir) && !dir.exists(OutputDir)) dir.create(OutputDir, recursive = TRUE)
    saveRDS(out, file.path(OutputDir, paste0(ModelName, "_VertebratePointResults.rds")))
  }
  out
}


###############################################################################
### Generic INApest point-transition pathogen facade
###
### The engine above is the validated standalone extension of the current
# INApestPointTransitionMatrix engine. With Vertebrate = NULL its biological
# population dynamics follow the parent point-transition model. The facade
### below exposes only the generic pathogen interaction to ordinary INApest use.
###############################################################################

.INApestPointTransitionMatrix_engine <- INApestVertebratePoint

# This public facade accepts the fully documented INApestVertebratePoint argument
# set through ... and adds only the optional point-pathogen interaction.

INApestPointTransitionMatrix <- function(..., Pathogen = NULL) {
  if (!is.null(Pathogen) && !inherits(Pathogen, "INApestPointPathogenInteraction"))
    stop("Pathogen must be NULL or created by INApestPointPathogenInteraction().")

  vertebrate_spec <- if (is.null(Pathogen)) NULL else list(Interaction = Pathogen)
  out <- .INApestPointTransitionMatrix_engine(..., Vertebrate = vertebrate_spec)

  if (!is.null(Pathogen)) {
    out$PathogenEvents <- out$InteractionEvents
  } else {
    out$PathogenEvents <- data.frame()
  }
  class(out) <- unique(c("INApestPointTransitionMatrix", class(out)))
  out
}