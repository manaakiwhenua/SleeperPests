###############################################################################
### INApestVertebrateNode -- stage-structured vertebrate node simulation engine
###
### Extends the node x life-stage architecture with animal-specific mechanisms.
### The optional Vertebrate list supplies Birth, HomeRange, Control and Interaction
### components while retaining the underlying Transition Matrix population,
### dispersal, surveillance, information and management machinery.
###
### Use this architecture when animal-specific spatial response, social structure
### or control changes the question; Vertebrate = NULL preserves parent behaviour.
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
    return(data.frame(parent_id = rep(points$id, counts), stage = 1L, stringsAsFactors = FALSE))
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

# Default stage-structured node dynamics used by the vertebrate extension.
.iv_local_dynamics_transition_matrix <- function(
    nodetransition = NodeTransition,
    weights = Weights,
    sddprob = SDDprob,
    nodeenvestabprob = NodeEnvEstabProb,
    n0 = N0,
    lddprob = LDDprob,
    lddrate = LDDrate,
    transition_sddprob = NULL,
    transition_lddprob = NULL,
    transition_lddrate = 0,
    nodeK = NodeK,
    node.seedbankK = NodeSeedbankK,
    nodepropaguleestablishment = NodePropaguleEstablishment,
    nodespreadreduction = NodeSpreadReduction,
    nodefecundityreduction = NodeFecundityReduction,
    nodecontrolfecundityreduction = 0,
    nodebirthmean = NULL,
    nodebirthmothers = NULL,
    managing = Managing,
    MaxInteger = MaxInteger,
    ApplyFootprintToTransitions = FALSE,
    BlockedTransitionMortality = 0,
    DispersalDensityFactor = 0
) {
  n_pops <- nrow(n0)
  S <- ncol(n0)
  n <- t(n0)

  if(is.null(dim(nodefecundityreduction))) {
    if(length(nodefecundityreduction) > 1 && n_pops == S && length(nodefecundityreduction) == n_pops)
      stop("nodefecundityreduction vector is ambiguous because n_pops equals S; supply an n_pops x S matrix")
    if(length(nodefecundityreduction) == 1) {
      nodefecundityreduction <- matrix(nodefecundityreduction, nrow = n_pops, ncol = S)
    } else if(length(nodefecundityreduction) == n_pops) {
      nodefecundityreduction <- matrix(rep(nodefecundityreduction, S), nrow = n_pops, ncol = S)
    } else if(length(nodefecundityreduction) == S) {
      nodefecundityreduction <- matrix(rep(nodefecundityreduction, each = n_pops), nrow = n_pops, ncol = S)
    } else {
      stop("nodefecundityreduction must be scalar, length n_pops, length S, or an n_pops x S matrix")
    }
  } else if(!identical(dim(nodefecundityreduction), c(n_pops, S))) {
    stop("nodefecundityreduction matrix must have dimensions n_pops x S")
  }
  if(any(!is.finite(nodefecundityreduction)) || any(nodefecundityreduction < 0) || any(nodefecundityreduction > 1))
    stop("nodefecundityreduction values must be between 0 and 1")
  # Independent, pre-deployed vertebrate control can also reduce fecundity.
  # It is multiplicative with information-triggered INApest management.
  if(is.null(dim(nodecontrolfecundityreduction))) {
    if(length(nodecontrolfecundityreduction) > 1 && n_pops == S &&
       length(nodecontrolfecundityreduction) == n_pops)
      stop("nodecontrolfecundityreduction vector is ambiguous because n_pops equals S; supply an n_pops x S matrix")
    if(length(nodecontrolfecundityreduction) == 1) {
      nodecontrolfecundityreduction <- matrix(nodecontrolfecundityreduction, nrow = n_pops, ncol = S)
    } else if(length(nodecontrolfecundityreduction) == n_pops) {
      nodecontrolfecundityreduction <- matrix(rep(nodecontrolfecundityreduction, S), nrow = n_pops, ncol = S)
    } else if(length(nodecontrolfecundityreduction) == S) {
      nodecontrolfecundityreduction <- matrix(rep(nodecontrolfecundityreduction, each = n_pops), nrow = n_pops, ncol = S)
    } else {
      stop("nodecontrolfecundityreduction must be scalar, length n_pops, length S, or an n_pops x S matrix")
    }
  } else if(!identical(dim(nodecontrolfecundityreduction), c(n_pops, S))) {
    stop("nodecontrolfecundityreduction matrix must have dimensions n_pops x S")
  }
  if(any(!is.finite(nodecontrolfecundityreduction)) ||
     any(nodecontrolfecundityreduction < 0) || any(nodecontrolfecundityreduction > 1))
    stop("nodecontrolfecundityreduction values must be between 0 and 1")

  EffectiveFecundityMultiplier <-
    (1 - nodefecundityreduction * managing) *
    (1 - nodecontrolfecundityreduction)

  # Stage weights may be supplied in the original form as one value per
  # stage, or as an n_pops x S matrix for node-specific stage weights.
  # A vector is expanded to an identical row for every population, exactly
  # preserving previous scale-invariance and density-dispersal applications.
  if (is.null(weights)) {
    w <- matrix(1, nrow = n_pops, ncol = S)
  } else if (is.matrix(weights)) {
    if (!identical(dim(weights), c(n_pops, S))) {
      stop("A weights matrix must have n_pops rows and S columns")
    }
    w <- weights
  } else {
    if (length(weights) != S) {
      stop("A weights vector must contain one value per stage")
    }
    w <- matrix(
      rep(weights, each = n_pops),
      nrow = n_pops,
      ncol = S
    )
  }
  
  if (any(!is.finite(w)) || any(w <= 0)) {
    stop("All stage weights must be finite and greater than zero")
  }
  
  # Recycle node-level inputs once rather than repeatedly inside loops.
  nodeK <- rep_len(nodeK, n_pops)
  node.seedbankK <- rep_len(node.seedbankK, n_pops)
  nodeenvestabprob <- rep_len(nodeenvestabprob, n_pops)
  nodepropaguleestablishment <- rep_len(nodepropaguleestablishment, n_pops)
  
  # Probability that a progression candidate dies when it fails to secure
  # a target-stage slot. A scalar applies to every source stage and node; a
  # vector of length S - 1 supplies one probability for each source stage;
  # and an n_pops x (S - 1) matrix permits node-by-stage values. The default
  # zero exactly preserves the former treatment of blocking as stasis.
  if (is.matrix(BlockedTransitionMortality)) {
    if (!identical(dim(BlockedTransitionMortality), c(n_pops, S - 1L))) {
      stop(
        "BlockedTransitionMortality must be scalar, length S - 1, or an n_pops x (S - 1) matrix"
      )
    }
    blocked_transition_mortality <- BlockedTransitionMortality
  } else if (length(BlockedTransitionMortality) == 1L) {
    blocked_transition_mortality <- matrix(
      BlockedTransitionMortality,
      nrow = n_pops,
      ncol = S - 1L
    )
  } else if (length(BlockedTransitionMortality) == S - 1L) {
    blocked_transition_mortality <- matrix(
      rep(BlockedTransitionMortality, each = n_pops),
      nrow = n_pops,
      ncol = S - 1L
    )
  } else {
    stop(
      "BlockedTransitionMortality must be scalar, length S - 1, or an n_pops x (S - 1) matrix"
    )
  }
  
  if (any(!is.finite(blocked_transition_mortality)) ||
      any(blocked_transition_mortality < 0 |
          blocked_transition_mortality > 1)) {
    stop("BlockedTransitionMortality probabilities must be between 0 and 1")
  }

  # ------------------------------------------------------------
  # Optional spatial movement during stage progression
  # ------------------------------------------------------------
  # A bare n_pops x n_pops matrix applies to every progression transition.
  # A list of length S - 1 permits transition-specific matrices; NULL list
  # elements retain the historical local-only transition for that source stage.
  # The top-level INApestMetaTransitionMatrix function resolves any time-varying
  # node x node x timestep arrays before calling LocalDynamics.
  .normalise_transition_movement <- function(x, name) {
    if (is.null(x)) {
      return(rep(list(NULL), S - 1L))
    }

    if (is.list(x)) {
      if (length(x) != S - 1L) {
        stop(name, " list must have length S - 1 (one element per progression transition)")
      }
      out <- x
    } else {
      out <- rep(list(x), S - 1L)
    }

    for (k in seq_len(S - 1L)) {
      P <- out[[k]]
      if (is.null(P)) next
      if (!is.matrix(P) || !identical(dim(P), c(n_pops, n_pops))) {
        stop(name, " entries must be NULL or n_pops x n_pops matrices")
      }
      if (any(!is.finite(P)) || any(P < 0)) {
        stop(name, " matrices must contain finite non-negative probabilities")
      }
      rs <- rowSums(P)
      if (any(rs > 1 + 1e-10)) {
        stop(name, " matrix row sums must not exceed 1; unused row mass represents movement outside the modelled landscape")
      }
      # Remove tiny floating-point excess so the residual outside probability
      # used by rmultinom is always non-negative.
      if (any(rs > 1)) {
        P <- P / pmax(1, rs)
        out[[k]] <- P
      }
    }
    out
  }

  transition_sdd <- .normalise_transition_movement(
    transition_sddprob,
    "transition_sddprob"
  )
  transition_ldd <- .normalise_transition_movement(
    transition_lddprob,
    "transition_lddprob"
  )

  if (!(length(transition_lddrate) %in% c(1L, S - 1L))) {
    stop("transition_lddrate must be scalar or length S - 1")
  }
  transition_lddrate <- rep_len(transition_lddrate, S - 1L)
  if (any(!is.finite(transition_lddrate)) ||
      any(transition_lddrate < 0 | transition_lddrate > 1)) {
    stop("transition_lddrate values must be between 0 and 1")
  }

  transition_movement_active <- vapply(
    seq_len(S - 1L),
    function(k) !is.null(transition_sdd[[k]]) || !is.null(transition_ldd[[k]]),
    logical(1)
  )

  # Allocate transitioning individuals from each source node to internal
  # destinations plus an implicit outside-landscape destination. As for the
  # existing transition-model dispersal semantics, missing row mass is export.
  .transition_flow_matrix <- function(source_counts, P) {
    flows <- matrix(0, nrow = n_pops, ncol = n_pops)
    exported <- numeric(n_pops)

    active_sources <- which(source_counts > 0)
    if (length(active_sources) == 0L) {
      return(list(flows = flows, exported = exported))
    }

    for (i in active_sources) {
      n_i <- floor(source_counts[i])
      if (n_i <= 0) next
      p_internal <- pmax(0, P[i, ])
      p_outside <- max(0, 1 - sum(p_internal))
      probs <- c(p_internal, p_outside)
      prob_total <- sum(probs)

      if (prob_total <= 0) next
      probs <- probs / prob_total

      if (n_i < MaxInteger) {
        allocation <- as.numeric(rmultinom(1, size = n_i, prob = probs))
      } else {
        # Preserve integer-valued stage abundances when counts exceed the
        # rmultinom integer limit. Rounding remainder is assigned to export.
        allocation <- floor(n_i * probs)
        remainder <- n_i - sum(allocation)
        allocation[length(allocation)] <- allocation[length(allocation)] + remainder
      }

      flows[i, ] <- allocation[seq_len(n_pops)]
      exported[i] <- allocation[n_pops + 1L]
    }

    list(flows = flows, exported = exported)
  }

  # Uniformly select accepted transitioners from source-specific arrivals at a
  # destination without allowing any source to contribute more than arrived.
  .take_without_replacement <- function(counts, n_take) {
    counts <- pmax(0, floor(counts))
    n_take <- min(floor(n_take), sum(counts))
    out <- numeric(length(counts))
    if (n_take <= 0 || sum(counts) <= 0) return(out)
    if (n_take >= sum(counts)) return(counts)

    if (sum(counts) >= MaxInteger) {
      expected <- n_take * counts / sum(counts)
      out <- floor(expected)
      remainder <- n_take - sum(out)
      if (remainder > 0) {
        spare <- counts - out
        ord <- order(expected - out, decreasing = TRUE)
        for (idx in ord) {
          if (remainder <= 0) break
          add <- min(spare[idx], remainder)
          if (add > 0) {
            out[idx] <- out[idx] + add
            remainder <- remainder - add
          }
        }
      }
      return(out)
    }

    remaining_take <- n_take
    remaining_total <- sum(counts)
    last_positive <- max(which(counts > 0))

    for (i in seq_along(counts)) {
      if (remaining_take <= 0) break
      if (counts[i] <= 0) {
        remaining_total <- remaining_total - counts[i]
        next
      }
      if (i == last_positive) {
        out[i] <- remaining_take
        break
      }
      out[i] <- rhyper(
        1,
        m = counts[i],
        n = remaining_total - counts[i],
        k = remaining_take
      )
      remaining_take <- remaining_take - out[i]
      remaining_total <- remaining_total - counts[i]
    }
    out
  }
  
  # ------------------------------------------------------------
  # STEP 1: Propagule production and transition extraction
  # ------------------------------------------------------------
  
  transition_is_list <- is.list(nodetransition)
  
  if (transition_is_list) {
    if (length(nodetransition) != n_pops) {
      stop("A nodetransition list must contain one matrix per population")
    }
    
    fecundity <- matrix(
      vapply(
        nodetransition,
        function(A) A[1, -1],
        numeric(S - 1)
      ),
      nrow = S - 1L,
      ncol = n_pops
    )
    
    transition_probabilities <- matrix(
      vapply(
        nodetransition,
        function(A) A[cbind(2:S, 1:(S - 1))],
        numeric(S - 1)
      ),
      nrow = S - 1L,
      ncol = n_pops
    )
    
    stasis_probabilities <- matrix(
      vapply(
        nodetransition,
        function(A) diag(A)[1:(S - 1)],
        numeric(S - 1)
      ),
      nrow = S - 1L,
      ncol = n_pops
    )
    
    terminal_survival_prob <- vapply(
      nodetransition,
      function(A) A[S, S],
      numeric(1)
    )
    
    EffectiveFecundity <- fecundity * t(EffectiveFecundityMultiplier[, 2:S, drop = FALSE])
    fec_means <- colSums(EffectiveFecundity * n[-1, , drop = FALSE])
    mother_counts <- colSums(
      n[-1, , drop = FALSE] * (fecundity > 0)
    )
  } else {
    fecundity <- nodetransition[1, -1]
    transition_probabilities <-
      nodetransition[cbind(2:S, 1:(S - 1))]
    stasis_probabilities <- diag(nodetransition)[1:(S - 1)]
    terminal_survival_prob <- nodetransition[S, S]
    
    EffectiveFecundity <- matrix(rep(fecundity, n_pops), nrow = S - 1L, ncol = n_pops) *
      t(EffectiveFecundityMultiplier[, 2:S, drop = FALSE])
    fec_means <- colSums(EffectiveFecundity * n[-1, , drop = FALSE])
    
    reproductive_stages <- which(fecundity > 0) + 1L
    mother_counts <- if (length(reproductive_stages)) {
      colSums(n[reproductive_stages, , drop = FALSE])
    } else {
      numeric(n_pops)
    }
  }
  
  # A custom vertebrate Birth hook may replace the transition-matrix
  # fecundity calculation with expected births per node. It is intentionally
  # still passed through the same stochastic propagule, dispersal and
  # establishment machinery. The hook receives management multipliers in its
  # context and is responsible for applying them to its own breeding logic.
  if (!is.null(nodebirthmean)) {
    nodebirthmean <- as.numeric(nodebirthmean)
    if (length(nodebirthmean) == 1L) nodebirthmean <- rep(nodebirthmean, n_pops)
    if (length(nodebirthmean) != n_pops || any(!is.finite(nodebirthmean)) || any(nodebirthmean < 0))
      stop("nodebirthmean must contain one finite non-negative expected birth count per node")
    fec_means <- nodebirthmean
  }
  if (!is.null(nodebirthmothers)) {
    nodebirthmothers <- as.numeric(nodebirthmothers)
    if (length(nodebirthmothers) == 1L) nodebirthmothers <- rep(nodebirthmothers, n_pops)
    if (length(nodebirthmothers) != n_pops || any(!is.finite(nodebirthmothers)) || any(nodebirthmothers < 0))
      stop("nodebirthmothers must contain one finite non-negative value per node")
    mother_counts <- nodebirthmothers
  }

  if (any(!is.finite(fec_means)) || any(fec_means < 0)) {
    stop("Propagule-production means must be finite and non-negative")
  }

  propagules <- rpois(n_pops, fec_means)
  
  # ------------------------------------------------------------
  # STEP 1b: Optional density-dependent animal dispersal
  # ------------------------------------------------------------
  #
  # Reweight dispersal before mortality and stage transitions. Propagules
  # produced in STEP 1 therefore respond to the same pre-transition census
  # that produced them. For the hornet model, gynes can consequently respond
  # to main nests that are present during gyne production, even though those
  # annual nests subsequently die during STEP 2.
  #
  # NA or 0 leaves the original dispersal matrix unchanged. For any positive
  # value, DispersalDensityFactor is the exponent alpha. The source-cell
  # diagonal is treated like every other destination. Probability absent from
  # an original row remains an outside-landscape option, permitting export.
  
  if (length(DispersalDensityFactor) != 1L) {
    stop("DispersalDensityFactor must be a single value")
  }
  
  density_dependent_dispersal <-
    !is.na(DispersalDensityFactor) &&
    DispersalDensityFactor != 0
  
  if (density_dependent_dispersal) {
    if (!is.finite(DispersalDensityFactor) ||
        DispersalDensityFactor <= 0) {
      stop(
        "DispersalDensityFactor must be finite and positive, or 0/NA to disable it"
      )
    }
    
    if (any(propagules > 0)) {
      # Weighted occupancy of stages 2:S in the pre-transition census.
      pre_transition_population <- rowSums(
        n0[, 2:S, drop = FALSE] * w[, 2:S, drop = FALSE]
      )
      
      # Cells with K <= 0 cannot attract dispersers. Calculating the ratio
      # only for positive capacities also avoids division by zero.
relative_occupancy <- rep(1, n_pops)
      positive_capacity <- nodeK > 0
      relative_occupancy[positive_capacity] <- pmin(
        1,
        pmax(
          0,
          pre_transition_population[positive_capacity] /
            nodeK[positive_capacity]
        )
      )
      
      # Zero occupancy gives attractiveness one; occupancy at or above K
      # gives attractiveness zero.
      cell_attractiveness <-
        pmax(0, 1 - relative_occupancy) ^ DispersalDensityFactor
      
      # The outside landscape retains constant attractiveness one.
      base_export_probability <- pmax(0, 1 - rowSums(sddprob))
      
      # R stores matrices column-wise. Repeating each destination value
      # n_pops times scales columns without sweep().
      cell_choice_weights <-
        sddprob * rep(cell_attractiveness, each = n_pops)
      
      total_choice_weights <-
        rowSums(cell_choice_weights) + base_export_probability
      
      # A row-length vector is recycled down every column, normalising the
      # destination and export choices without sweep() or matrix subsetting.
      row_multiplier <- numeric(n_pops)
      valid_rows <- total_choice_weights > 0
      row_multiplier[valid_rows] <- 1 / total_choice_weights[valid_rows]
      
      # If all reachable cells are full and no original export probability
      # exists, the adjusted row remains zero and propagules find no modelled
      # destination.
      sddprob <- cell_choice_weights * row_multiplier
    }
  }
  
  # ------------------------------------------------------------
  # STEP 2: Stage transitions
  # ------------------------------------------------------------
  
  if (any(!is.finite(terminal_survival_prob)) ||
      any(terminal_survival_prob < 0 | terminal_survival_prob > 1)) {
    stop("Terminal-stage survival probabilities must be between 0 and 1")
  }
  
  n[S, ] <- rbinom(
    n = n_pops,
    size = n[S, ],
    prob = terminal_survival_prob
  )
  
  transition_footprint <- if (ApplyFootprintToTransitions) {
    pmin(1, pmax(0, nodepropaguleestablishment))
  } else {
    1
  }
  
  # Weighted population in stages already processed above the target stage.
  capacity_above <- numeric(n_pops)
  
  # The descending order prevents multiple stage transitions per timestep.
  for (s in S:2) {
    N_prev <- n[s - 1, ]
    
    trans_prob <- if (transition_is_list) {
      transition_probabilities[s - 1, ]
    } else {
      transition_probabilities[s - 1]
    }
    
    stay_prob <- if (transition_is_list) {
      stasis_probabilities[s - 1, ]
    } else {
      stasis_probabilities[s - 1]
    }
    
    surv_total_prob <- pmin(trans_prob + stay_prob, 1)
    
    N_surv <- rbinom(
      n = n_pops,
      size = N_prev,
      prob = surv_total_prob
    )
    
    prog_cond_prob <- ifelse(
      surv_total_prob > 0,
      trans_prob / surv_total_prob,
      0
    )
    
    n_trans_candidates <- rbinom(
      n = n_pops,
      size = N_surv,
      prob = prog_cond_prob
    )
    
    # STEP 2c: target stage and all higher stages consume capacity.
    # stage_weight is node-specific when weights was supplied as a matrix.
    stage_weight <- w[, s]
    total_pop <- capacity_above + n[s, ] * stage_weight
    slots_available <- pmax(
      0,
      floor((nodeK - total_pop) / stage_weight)
    )
    
    max_slots <- nodeK / stage_weight
    free_fraction <- ifelse(
      max_slots > 0,
      slots_available / max_slots,
      0
    )
    
    candidate_prob <- pmin(
      1,
      pmax(0, transition_footprint * free_fraction)
    )
    
    # A transition with no transition-specific movement matrices follows the
    # historical local-only path exactly, including its random-number sequence.
    if (!transition_movement_active[s - 1L]) {
      n_trans_actual <- rbinom(
        n = n_pops,
        size = n_trans_candidates,
        prob = candidate_prob
      )
      n_trans_actual <- pmin(n_trans_actual, slots_available)
      n_blocked <- n_trans_candidates - n_trans_actual
      transition_additions <- n_trans_actual
    } else {
      P_sdd <- transition_sdd[[s - 1L]]
      P_ldd <- transition_ldd[[s - 1L]]

      # If both movement modes are supplied, transition_lddrate gives the
      # fraction entering LDD. If only one matrix is supplied, every moving
      # transition candidate uses that matrix and the rate is immaterial.
      if (!is.null(P_sdd) && !is.null(P_ldd)) {
        r_ldd <- transition_lddrate[s - 1L]
        if (r_ldd <= 0) {
          ldd_candidates <- numeric(n_pops)
          sdd_candidates <- n_trans_candidates
        } else if (r_ldd >= 1) {
          ldd_candidates <- n_trans_candidates
          sdd_candidates <- numeric(n_pops)
        } else {
          ldd_candidates <- rbinom(
            n = n_pops,
            size = n_trans_candidates,
            prob = r_ldd
          )
          sdd_candidates <- n_trans_candidates - ldd_candidates
        }
      } else if (!is.null(P_sdd)) {
        sdd_candidates <- n_trans_candidates
        ldd_candidates <- numeric(n_pops)
      } else {
        sdd_candidates <- numeric(n_pops)
        ldd_candidates <- n_trans_candidates
      }

      sdd_flow <- if (!is.null(P_sdd)) {
        .transition_flow_matrix(sdd_candidates, P_sdd)
      } else {
        list(flows = matrix(0, n_pops, n_pops), exported = numeric(n_pops))
      }
      ldd_flow <- if (!is.null(P_ldd)) {
        .transition_flow_matrix(ldd_candidates, P_ldd)
      } else {
        list(flows = matrix(0, n_pops, n_pops), exported = numeric(n_pops))
      }

      internal_flows <- sdd_flow$flows + ldd_flow$flows
      exported_by_source <- sdd_flow$exported + ldd_flow$exported
      arrivals <- colSums(internal_flows)

      # Capacity and optional transition footprint are evaluated at the
      # destination of the moving transition. Candidates exported beyond the
      # modelled landscape are successful departures and do not compete for
      # internal slots.
      accepted_by_destination <- rbinom(
        n = n_pops,
        size = arrivals,
        prob = candidate_prob
      )
      accepted_by_destination <- pmin(
        accepted_by_destination,
        slots_available,
        arrivals
      )

      accepted_by_source <- numeric(n_pops)
      destination_nodes <- which(arrivals > 0 & accepted_by_destination > 0)
      for (j in destination_nodes) {
        accepted_sources_j <- .take_without_replacement(
          internal_flows[, j],
          accepted_by_destination[j]
        )
        accepted_by_source <- accepted_by_source + accepted_sources_j
      }

      internal_by_source <- rowSums(internal_flows)
      n_blocked <- pmax(0, internal_by_source - accepted_by_source)
      transition_additions <- accepted_by_destination

      # Sanity check: every progression candidate is either accepted into an
      # internal destination, blocked internally, or exported from the model.
      if (any(abs(
        n_trans_candidates -
          (accepted_by_source + n_blocked + exported_by_source)
      ) > 1e-8)) {
        stop("Internal error while allocating movement during stage transition")
      }
    }
    
    # STEP 2e: resolve stasis and mortality of blocked candidates.
    # Individuals allocated to stasis are unaffected. Only candidates that
    # attempted progression, remained inside the modelled landscape, and did
    # not obtain a target-stage slot are exposed to BlockedTransitionMortality.
    # Blocked survivors return to their source node and source stage. Movement
    # outside the modelled landscape is treated as export, not blocking.
    n_stay <- N_surv - n_trans_candidates
    blocked_mortality_prob <- blocked_transition_mortality[, s - 1L]
    n_blocked_survivors <- n_blocked
    
    # Avoid unnecessary random draws at the exact 0 and 1 limits. In
    # particular, BlockedTransitionMortality = 0 preserves the previous RNG
    # sequence as well as the previous biological behaviour.
    certain_death <- n_blocked > 0 & blocked_mortality_prob == 1
    n_blocked_survivors[certain_death] <- 0
    
    stochastic_blocking <-
      n_blocked > 0 &
      blocked_mortality_prob > 0 &
      blocked_mortality_prob < 1
    
    if (any(stochastic_blocking)) {
      n_blocked_survivors[stochastic_blocking] <- rbinom(
        n = sum(stochastic_blocking),
        size = n_blocked[stochastic_blocking],
        prob = 1 - blocked_mortality_prob[stochastic_blocking]
      )
    }
    
    n[s, ] <- n[s, ] + transition_additions
    n[s - 1, ] <- n_stay + n_blocked_survivors
    
    # This now equals the weighted population in stages s:S at each
    # destination node. For local transitions transition_additions is the
    # historical n_trans_actual vector; for moving transitions it is the
    # accepted destination total.
    capacity_above <- total_pop + transition_additions * stage_weight
  }
  
  # ------------------------------------------------------------
  # STEP 3: Propagule dispersal
  # ------------------------------------------------------------
  # Route realised integer propagules, rather than converting an expected
  # fractional number of movers into rmultinom(size). This is the same
  # correction validated for INApestMeta: first assign each realised
  # propagule to SDD/LDD, then allocate it from its source row to an internal
  # destination or the residual outside-landscape category.
  Pin <- numeric(n_pops)
  Qin <- numeric(n_pops)

  # Route realised reproductive propagules through a transition-matrix dispersal surface.
  .route_integer_propagules_tm <- function(source_counts, P, label) {
    arrivals <- numeric(n_pops)
    if(!is.matrix(P)) return(arrivals)
    if(!identical(dim(P),c(n_pops,n_pops))) stop(label," must be an n_pops x n_pops matrix")
    active_sources <- which(source_counts > 0)
    if(!length(active_sources)) return(arrivals)
    for(i in active_sources) {
      n_i <- floor(source_counts[i])
      if(n_i <= 0) next
      p_internal <- pmax(0,as.numeric(P[i,]))
      rs <- sum(p_internal)
      if(!is.finite(rs) || rs > 1 + 1e-10) stop(label," source-row probabilities must be finite, non-negative, and sum to at most 1")
      if(rs > 1) p_internal <- p_internal / rs
      probs <- c(p_internal,max(0,1-sum(p_internal)))
      if(sum(probs) <= 0) next
      if(n_i < MaxInteger) {
        allocation <- as.numeric(rmultinom(1,size=n_i,prob=probs))
      } else {
        # Retain the previous large-count fallback without fractional rmultinom sizes.
        allocation <- floor(n_i * probs)
        allocation[length(allocation)] <- allocation[length(allocation)] + (n_i-sum(allocation))
      }
      arrivals <- arrivals + allocation[seq_len(n_pops)]
    }
    arrivals
  }

  if(any(propagules > 0)) {
    # Integer SDD/LDD branching gives the intended Poisson-thinning
    # marginals when propagule production itself is Poisson.
    if(lddrate <= 0) {
      ldd_sources <- numeric(n_pops)
      sdd_sources <- propagules
    } else if(lddrate >= 1) {
      ldd_sources <- propagules
      sdd_sources <- numeric(n_pops)
    } else {
      ldd_sources <- rbinom(n_pops,size=propagules,prob=lddrate)
      sdd_sources <- propagules - ldd_sources
    }

    Pin <- .route_integer_propagules_tm(sdd_sources,sddprob,"sddprob")

    if(is.matrix(lddprob) && any(ldd_sources > 0)) {
      spread_reduction <- pmin(1,pmax(0,rep_len(nodespreadreduction,n_pops) * rep_len(managing,n_pops)))
      keep_prob <- 1-spread_reduction
      keep_sources <- ldd_sources
      stochastic_keep <- ldd_sources > 0 & keep_prob > 0 & keep_prob < 1
      keep_sources[keep_prob <= 0] <- 0
      if(any(stochastic_keep))
        keep_sources[stochastic_keep] <- rbinom(sum(stochastic_keep),size=ldd_sources[stochastic_keep],prob=keep_prob[stochastic_keep])
      Qin <- .route_integer_propagules_tm(keep_sources,lddprob,"lddprob")
    }
  }

  # ------------------------------------------------------------
  # STEP 4: Recruitment into seedbank
  # ------------------------------------------------------------
  
  seedbank_slots <- pmax(0, floor(node.seedbankK - n[1, ]))
  Pin <- pmax(0, floor(Pin))
  Qin <- pmax(0, floor(Qin))
  env_prob <- pmin(1, pmax(0, nodeenvestabprob))
  
  # A value below one retains the passive-propagule maternal-footprint model.
  # Setting every value to one removes mother dependence: any naturally
  # dispersed propagule that reaches a cell can search all its free slots.
  unrestricted_natural_search <- all(
    is.finite(nodepropaguleestablishment) &
      nodepropaguleestablishment >= 1
  )
  
  if (unrestricted_natural_search) {
    accessible_fraction <- as.numeric(Pin > 0)
    coverage_pressure <- accessible_fraction
  } else {
    footprint_area <- nodepropaguleestablishment * nodeK
    coverage_pressure <- numeric(n_pops)
    
    if (any(Pin > 0) && any(mother_counts > 0)) {
      coverage_numerator <- as.numeric(
        (mother_counts * footprint_area) %*% sddprob
      )
      
      positive_capacity <- nodeK > 0
      coverage_pressure[positive_capacity] <-
        coverage_numerator[positive_capacity] /
        nodeK[positive_capacity]
    }
    
    accessible_fraction <- pmin(
      1,
      pmax(0, -expm1(-coverage_pressure))
    )
  }
  
  natural_slots <- floor(seedbank_slots * accessible_fraction)
  
  # Retain rare, non-zero colonisation links without creating slots when
  # no naturally dispersed propagule arrived.
  natural_slots <- ifelse(
    Pin > 0 & accessible_fraction > 0,
    pmax(1, natural_slots),
    0
  )
  
  natural_slots <- pmin(natural_slots, seedbank_slots)
  other_slots <- pmax(0, seedbank_slots - natural_slots)
  
  lambda_P <- ifelse(
    natural_slots > 0,
    env_prob * Pin / pmax(natural_slots, 1),
    0
  )
  
  lambda_Q <- ifelse(
    seedbank_slots > 0,
    env_prob * Qin / pmax(seedbank_slots, 1),
    0
  )
  
  # Qin contributes inside and outside the natural footprint; Pin only inside.
  p_natural_slots <- -expm1(-(lambda_P + lambda_Q))
  p_other_slots <- -expm1(-lambda_Q)
  
  recruits <-
    rbinom(n_pops, natural_slots, p_natural_slots) +
    rbinom(n_pops, other_slots, p_other_slots)
  
  recruits <- pmin(recruits, Pin + Qin)
  n[1, ] <- n[1, ] + recruits
  
  t(n)
}

INApestVertebrateNode = function(
ModelName = "INApestVertebrateNode",            # Model and output name
Nperm,                                          # Number of stochastic simulation runs
Ntimesteps,                                     # Timesteps in each simulation
Nstages,                                        # Number of demographic/life stages
Weights,                                        # Stage weights used for population/capacity totals
Transition,                                     # Stage-transition matrix or node/time-varying equivalent
LocalDynamics = .iv_local_dynamics_transition_matrix, # Local population growth, movement and management function
DetectionProb,                                  # Background-surveillance detection probability
DetectionSD = NULL,                             # Variation in background detection probability
ManageProb,                                     # Management probability when information is available
ManageSD = NULL,                                # Variation in management probability
MortalityProb,                                  # Management-driven host mortality probability
MortalitySD = NULL,                             # Variation in management mortality probability
FecundityReduction = 0,                         # Proportional reduction in reproduction under management
SpreadReduction,                                # Proportional reduction in dispersal under management
SpreadReductionSD = NULL,                       # Variation in spread reduction
InitialPopulation = NA,                         # Starting vertebrate abundance by node x stage
InitBioP = NA,                                  # Proportion of nodes initially invaded
InvasionRisk = NA,                              # External invasion probability or node weighting
InitialInfo = NA,                               # Starting information state
InitInfoP = NA,                                 # Proportion of nodes initially with information
ExternalInfoProb = 0.0,                         # Information arriving from outside the modelled system
InfoRetentionProb = 1,                          # Probability existing information persists one timestep
InfoPersistenceSteps = NA,                      # Programmed timesteps information persists after evidence
EnvEstabProb = 1,                               # Environmental establishment probability
K,                                              # Host carrying capacity by node
SeedbankK,                                      # Stage-1 / seedbank carrying capacity by node
PropaguleEstablishment,                         # Establishment probability for arriving propagules
IncursionStartPop=NA,                           # Host abundance assigned to new external incursions
SDDprob,                                        # Short-distance source-to-target dispersal
SEAM = 0,                                       # Information-transfer adjacency between nodes
LDDprob = NA,                                   # Long-distance source-to-target dispersal
LDDrate = 0,                                    # Fraction of propagules entering long-distance dispersal
TransitionSDDprob = NULL,                       # Short-distance movement during stage progression
TransitionLDDprob = NULL,                       # Long-distance movement during stage progression
TransitionLDDrate = 0,                          # LDD fraction for moving stage transitions
DispersalDensityFactor = 0,                     # Strength of density-dependent dispersal redistribution
BlockedTransitionMortality = 0,                 # Mortality when stage progression is blocked
OngoingExternalInvasion = F,                    # Allow new host incursions after initialisation
OngoingExternalInfo = F,                        # Allow new external information after initialisation
Vertebrate = NULL,                              # Optional Birth/HomeRange/Control/Interaction modules
OutputDir = NA,                                 # Directory for saved outputs
DoPlots = TRUE,                                 # Legacy plotting option; plotting is post-processing
InfoTriggeredDetectionProb = 0,                 # Detection probability where information already exists
InfoTriggeredDetectionSD = NULL,                # Variation in information-triggered detection
SaveResults = TRUE,                             # Save standard simulation outputs to disk
DoProgress = TRUE                               # Print simulation progress to the console
)
{

# ---------------------------------------------------------------------------
# Set up and validate the vertebrate node simulation.
# ---------------------------------------------------------------------------
if (!exists(".iv_validate_vertebrate", mode = "function"))
  stop("INApestVertebrateHelpers.R must be sourced before INApestVertebrateNode().")
.iv_validate_vertebrate(Vertebrate)
BirthModule <- .iv_module(Vertebrate, "Birth")
HomeRangeModule <- .iv_module(Vertebrate, "HomeRange")
ControlModule <- .iv_module(Vertebrate, "Control")
InteractionModule <- .iv_module(Vertebrate, "Interaction")
if(!is.function(LocalDynamics))
  stop("LocalDynamics must be a function")
# Force the argument before any parallel worker closure is created. This keeps
# the selected default or user-supplied function as an explicit model input.
force(LocalDynamics)
# Validate the host-engine contract before starting stochastic simulation.
n_nodes <- nrow(SDDprob)
if(!is.numeric(Nperm) || length(Nperm)!=1L || !is.finite(Nperm) || Nperm<1 || Nperm!=floor(Nperm)) stop("Nperm must be a positive integer")
if(!is.numeric(Ntimesteps) || length(Ntimesteps)!=1L || !is.finite(Ntimesteps) || Ntimesteps<1 || Ntimesteps!=floor(Ntimesteps)) stop("Ntimesteps must be a positive integer")
if(!is.numeric(Nstages) || length(Nstages)!=1L || !is.finite(Nstages) || Nstages<2 || Nstages!=floor(Nstages)) stop("Nstages must be an integer of at least 2")
if(is.null(dim(SDDprob))) stop("SDDprob must be a nodes x nodes matrix or nodes x nodes x Ntimesteps array")
n_nodes <- dim(SDDprob)[1L]
# Validate fixed or time-varying node connectivity.
.ValidateConnectivityTM <- function(x,name,allow_disabled_scalar=FALSE) {
  if(allow_disabled_scalar && length(x)==1L && (is.na(x) || identical(as.numeric(x),0))) return(invisible(TRUE))
  d <- dim(x)
  if(is.null(d)) stop(name," must be a nodes x nodes matrix or nodes x nodes x Ntimesteps array")
  if(length(d)==2L) {
    if(!all(d==c(n_nodes,n_nodes))) stop(name," matrix must have dimensions nodes x nodes")
    mats <- list(x)
  } else if(length(d)==3L) {
    if(!all(d==c(n_nodes,n_nodes,Ntimesteps))) stop(name," 3D array must have dimensions nodes x nodes x Ntimesteps")
    mats <- lapply(seq_len(Ntimesteps),function(tt) matrix(x[,,tt],nrow=n_nodes,ncol=n_nodes))
  } else stop(name," must be a matrix or 3D array")
  for(M in mats) {
    if(any(!is.finite(M)) || any(M<0)) stop(name," must contain finite non-negative probabilities")
    if(any(rowSums(M)>1+1e-10)) stop(name," source-row sums must not exceed 1")
  }
  invisible(TRUE)
}
.ValidateConnectivityTM(SDDprob,"SDDprob",FALSE)
.ValidateConnectivityTM(LDDprob,"LDDprob",TRUE)
# Return the connectivity matrix used in the current timestep.
.SliceConnectivityTM <- function(x,timestep) {
  if(length(x)==1L && is.na(x)) return(x)
  if(length(dim(x))==3L) return(matrix(x[,,timestep],nrow=n_nodes,ncol=n_nodes))
  x
}
# Validate supported demographic transition-matrix forms.
.ValidateTransitionTM <- function(x) {
  if(is.matrix(x)) {
    if(!all(dim(x)==c(Nstages,Nstages))) stop("Transition matrix must have dimensions Nstages x Nstages")
    if(any(!is.finite(x)) || any(x<0)) stop("Transition must contain finite non-negative values")
    return(invisible(TRUE))
  }
  if(is.list(x)) {
    if(length(x)!=n_nodes) stop("Transition list must contain one Nstages x Nstages matrix per node")
    for(ii in seq_len(n_nodes)) {
      if(!is.matrix(x[[ii]]) || !all(dim(x[[ii]])==c(Nstages,Nstages))) stop("Every Transition list entry must be an Nstages x Nstages matrix")
      if(any(!is.finite(x[[ii]])) || any(x[[ii]]<0)) stop("Transition must contain finite non-negative values")
    }
    return(invisible(TRUE))
  }
  d <- dim(x)
  if(length(d)==3L) {
    if(!all(d==c(Nstages,Nstages,Ntimesteps))) stop("3D Transition array must have dimensions Nstages x Nstages x Ntimesteps")
  } else if(length(d)==4L) {
    if(!all(d==c(Nstages,Nstages,n_nodes,Ntimesteps))) stop("4D Transition array must have dimensions Nstages x Nstages x nodes x Ntimesteps")
  } else stop("Transition must be a matrix, a per-node list, a Nstages x Nstages x Ntimesteps array, or a Nstages x Nstages x nodes x Ntimesteps array")
  if(any(!is.finite(x)) || any(x<0)) stop("Transition must contain finite non-negative values")
  invisible(TRUE)
}
# Resolve demographic transition matrices for the current timestep.
.ResolveTransitionTM <- function(x,timestep) {
  if(is.matrix(x) || is.list(x)) return(x)
  d <- dim(x)
  if(length(d)==3L) return(matrix(x[,,timestep],nrow=Nstages,ncol=Nstages))
  lapply(seq_len(n_nodes),function(ii) matrix(x[,,ii,timestep],nrow=Nstages,ncol=Nstages))
}
.ValidateTransitionTM(Transition)
if(is.matrix(Weights)) {
  if(!all(dim(Weights)==c(n_nodes,Nstages))) stop("Weights matrix must have dimensions nodes x Nstages")
} else if(length(Weights)!=Nstages) stop("Weights must be length Nstages or a nodes x Nstages matrix")
if(any(!is.finite(Weights)) || any(Weights<=0)) stop("Weights must contain finite values greater than zero")
# Validate node-level scalar/vector/time-varying inputs.
.ValidateNodeTimeTM <- function(x,name,unit_interval=FALSE,allow_na_scalar=FALSE,whole=FALSE) {
  if(allow_na_scalar && length(x)==1L && is.na(x)) return(invisible(TRUE))
  d <- dim(x); ok <- FALSE
  if(is.null(d)) ok <- length(x) %in% c(1L,n_nodes)
  if(!is.null(d) && length(d)==2L) ok <- all(d==c(n_nodes,Ntimesteps))
  if(!ok) stop(name," must be scalar, length nodes, or nodes x Ntimesteps")
  if(any(!is.finite(x))) stop(name," must contain finite values")
  if(any(x<0)) stop(name," must contain non-negative values")
  if(unit_interval && any(x>1)) stop(name," values must be between 0 and 1")
  if(whole && any(x!=floor(x))) stop(name," must contain whole numbers")
  invisible(TRUE)
}
.ValidateNodeTimeTM(K,"K",FALSE,FALSE,FALSE)
.ValidateNodeTimeTM(SeedbankK,"SeedbankK",FALSE,FALSE,FALSE)
.ValidateNodeTimeTM(EnvEstabProb,"EnvEstabProb",TRUE,FALSE,FALSE)
.ValidateNodeTimeTM(PropaguleEstablishment,"PropaguleEstablishment",FALSE,FALSE,FALSE)
.ValidateNodeTimeTM(ManageProb,"ManageProb",TRUE,FALSE,FALSE)
.ValidateNodeTimeTM(SpreadReduction,"SpreadReduction",TRUE,FALSE,FALSE)
.ValidateNodeTimeTM(ExternalInfoProb,"ExternalInfoProb",TRUE,TRUE,FALSE)
.ValidateNodeTimeTM(InvasionRisk,"InvasionRisk",TRUE,TRUE,FALSE)
if(!is.na(InitBioP) && (!is.numeric(InitBioP) || length(InitBioP)!=1L || !is.finite(InitBioP) || InitBioP<0 || InitBioP>1)) stop("InitBioP must be NA or one probability between 0 and 1")
if(!(length(IncursionStartPop)==1L && is.na(IncursionStartPop)) && (!is.numeric(IncursionStartPop) || length(IncursionStartPop)!=1L || !is.finite(IncursionStartPop) || IncursionStartPop<0 || IncursionStartPop!=floor(IncursionStartPop))) stop("IncursionStartPop must be NA or one non-negative whole number")
# Validate stage-specific probability inputs.
.ValidateStageProbTM <- function(x,name) {
  d <- dim(x); ok <- FALSE
  if(is.null(d)) ok <- length(x) %in% c(1L,Nstages)
  if(!is.null(d) && length(d)==2L) ok <- all(d==c(n_nodes,Nstages))
  if(!is.null(d) && length(d)==3L) ok <- all(d==c(n_nodes,Nstages,Ntimesteps))
  if(!ok) stop(name," must be scalar, length Nstages, nodes x Nstages, or nodes x Nstages x Ntimesteps")
  if(any(!is.finite(x)) || any(x<0 | x>1)) stop(name," values must be probabilities between 0 and 1")
  invisible(TRUE)
}
.ValidateStageProbTM(DetectionProb,"DetectionProb")
.ValidateStageProbTM(MortalityProb,"MortalityProb")
if(!(length(InitialPopulation)==1L && is.na(InitialPopulation))) {
  if(!is.matrix(InitialPopulation) || !all(dim(InitialPopulation)==c(n_nodes,Nstages))) stop("InitialPopulation must be NA or a nodes x Nstages matrix")
  if(any(!is.finite(InitialPopulation)) || any(InitialPopulation<0) || any(InitialPopulation!=floor(InitialPopulation))) stop("InitialPopulation must contain non-negative whole numbers")
}
if(isTRUE(OngoingExternalInfo) && length(ExternalInfoProb)==1L && is.na(ExternalInfoProb)) stop("OngoingExternalInfo requires finite ExternalInfoProb")
if(isTRUE(OngoingExternalInvasion) && length(InvasionRisk)==1L && is.na(InvasionRisk)) stop("OngoingExternalInvasion requires finite InvasionRisk")
if(!is.numeric(LDDrate) || length(LDDrate)!=1L || !is.finite(LDDrate) || LDDrate<0 || LDDrate>1) stop("LDDrate must be one probability between 0 and 1")
if(!is.numeric(TransitionLDDrate) || any(!is.finite(TransitionLDDrate)) || any(TransitionLDDrate<0 | TransitionLDDrate>1)) stop("TransitionLDDrate values must be probabilities between 0 and 1")
if(!(length(DispersalDensityFactor)==1L && (is.na(DispersalDensityFactor) || (is.finite(DispersalDensityFactor) && DispersalDensityFactor>=0)))) stop("DispersalDensityFactor must be NA, zero, or one finite positive value")
# Host-engine contract validation is complete.

  
# Max integer for propagule dispersal using rmultinom
MaxInteger <- .Machine$integer.max  
  
# Connectivity was validated with the other host-engine inputs above.

# Optional movement while individuals progress between life-history stages.
# A bare matrix/3D array applies to every progression transition. A list of
# length Nstages-1 allows each source-stage transition to have its own matrix,
# time-varying array, or NULL (local-only transition).
n_nodes <- nrow(SDDprob)
# Validate stage-transition movement inputs.
ValidateTransitionMovement <- function(x, name) {
  validate_one <- function(z, label) {
    if(is.null(z)) return(invisible(NULL))
    dz <- dim(z)
    if(is.null(dz))
      stop(label, " must be a nodes x nodes matrix or nodes x nodes x Ntimesteps array")
    if(length(dz) == 2L) {
      if(!all(dz == c(n_nodes, n_nodes)))
        stop(label, " matrix must have dimensions nodes x nodes")
    } else if(length(dz) == 3L) {
      if(!all(dz == c(n_nodes, n_nodes, Ntimesteps)))
        stop(label, " 3D array must have dimensions nodes x nodes x Ntimesteps")
    } else {
      stop(label, " must be a matrix or 3D array")
    }
    if(any(!is.finite(z)) || any(z < 0))
      stop(label, " must contain finite non-negative probabilities")
    if(length(dz) == 2L) {
      if(any(rowSums(z) > 1 + 1e-10))
        stop(label, " row sums must not exceed 1")
    } else {
      for(tt in seq_len(Ntimesteps))
        if(any(rowSums(matrix(z[,,tt],nrow=n_nodes,ncol=n_nodes)) > 1 + 1e-10))
          stop(label, " row sums must not exceed 1 in any timestep")
    }
    invisible(NULL)
  }

  if(is.null(x)) return(invisible(NULL))
  if(is.list(x)) {
    if(length(x) != Nstages - 1L)
      stop(name, " list must have length Nstages - 1")
    for(k in seq_len(Nstages - 1L))
      validate_one(x[[k]], paste0(name, "[[", k, "]]"))
  } else {
    validate_one(x, name)
  }
  invisible(NULL)
}

# Resolve stage-transition movement for the current timestep.
ResolveTransitionMovement <- function(x, timestep) {
  resolve_one <- function(z) {
    if(is.null(z)) return(NULL)
    if(length(dim(z)) == 3L) return(matrix(z[,,timestep],nrow=n_nodes,ncol=n_nodes))
    z
  }
  if(is.null(x)) return(NULL)
  if(is.list(x)) return(lapply(x, resolve_one))
  resolve_one(x)
}

# Check whether stage-transition movement is active.
HasTransitionMovement <- function(x) {
  if(is.null(x)) return(FALSE)
  if(is.list(x)) return(any(!vapply(x, is.null, logical(1))))
  TRUE
}

ValidateTransitionMovement(TransitionSDDprob, "TransitionSDDprob")
ValidateTransitionMovement(TransitionLDDprob, "TransitionLDDprob")
if(!(length(TransitionLDDrate) %in% c(1L, Nstages - 1L)))
  stop("TransitionLDDrate must be scalar or length Nstages - 1")
if(any(!is.finite(TransitionLDDrate)) || any(TransitionLDDrate < 0) || any(TransitionLDDrate > 1))
  stop("TransitionLDDrate values must be between 0 and 1")
TransitionMovementConfigured <- HasTransitionMovement(TransitionSDDprob) || HasTransitionMovement(TransitionLDDprob)

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

# Use first-timestep carrying capacity when K varies through time.
if(is.matrix(K) == TRUE)
  {
  K_is_0 <- K[,1]<=0
  inv_K <- 1 / sum(K[,1])
  NodeK = K[,1] 
}    

# --- Seedbank carrying capacity ---
if (is.matrix(SeedbankK) == FALSE)
{
  NodeSeedbankK <- SeedbankK
}

if (is.matrix(SeedbankK) == TRUE)
{
  NodeSeedbankK <- SeedbankK[, 1]
}

  
if(is.matrix(PropaguleEstablishment) == FALSE)
  NodePropaguleEstablishment = PropaguleEstablishment

if(is.matrix(EnvEstabProb) == F)
  NodeEnvEstabProb <- EnvEstabProb

# Allocate host-abundance histories.
PopulationResults = array(dim = c(nrow(SDDprob),Ntimesteps,Nperm))

# Allocate host-abundance histories.
PopulationStageResults = array(dim = c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))

# Allocate host/pest presence histories.
InvasionResults = array(dim = c(nrow(SDDprob),Ntimesteps,Nperm))


# Allocate known-presence and surveillance histories.
DetectedResults = InvasionResults

# PoA observation contract. Routine pre-deployed vertebrate control is a
# Background observation opportunity because it operates regardless of HaveInfo.
# Ordinary Background and InfoTriggered surveillance remain separate draws.
BackgroundDetectedResults = array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
BackgroundSurveillanceDetectedResults = array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
InfoTriggeredDetectedResults = array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
InformationStateBeforeSurveillanceResults = array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
HaveInfoResults = array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
BackgroundDetectionProbabilityResults = array(0, dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))
InfoTriggeredDetectionProbabilityResults = array(0, dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))
# Pre-control state/probability are retained so the PoA companion can calculate
# the no-detection likelihood outside the biological engine.
RoutineControlObservationAbundanceResults = array(0, dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))
RoutineControlDetectionProbabilityResults = array(0, dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))

# Allocate management histories.
ManagingResults = InvasionResults

# Vertebrate specialist control diagnostics
ControlDeathResults = array(0, dim = c(nrow(SDDprob), Nstages, Ntimesteps, Nperm))
ControlDetectionResults = array(0, dim = c(nrow(SDDprob), Ntimesteps, Nperm))
ControlCostResults = matrix(0, nrow = Ntimesteps, ncol = Nperm)

# Validate the socioeconomic information-transfer network when supplied.
if(is.matrix(SEAM) == T)
     {
     diag(SEAM) = 0
     RandSEAM <- matrix(NA,nrow = nrow(SDDprob),ncol=nrow(SDDprob))
     }

# Validate management-induced fecundity reduction
n_nodes <- nrow(SDDprob)
FecundityReductionDims <- dim(FecundityReduction)
if(is.null(FecundityReductionDims)) {
  valid_lengths <- unique(c(1L, n_nodes, Nstages))
  if(!(length(FecundityReduction) %in% valid_lengths))
    stop("FecundityReduction must be scalar, length nodes, length Nstages, nodes x Ntimesteps, or nodes x Nstages x Ntimesteps")
  if(length(FecundityReduction) > 1 && n_nodes == Nstages && length(FecundityReduction) == n_nodes)
    stop("FecundityReduction vector is ambiguous because number of nodes equals Nstages; use a nodes x Ntimesteps matrix for node-specific values or a nodes x Nstages x Ntimesteps array for stage-specific values")
} else if(length(FecundityReductionDims) == 2) {
  if(!all(FecundityReductionDims == c(n_nodes, Ntimesteps)))
    stop("FecundityReduction matrix must have dimensions nodes x Ntimesteps")
} else if(length(FecundityReductionDims) == 3) {
  if(!all(FecundityReductionDims == c(n_nodes, Nstages, Ntimesteps)))
    stop("FecundityReduction array must have dimensions nodes x Nstages x Ntimesteps")
} else {
  stop("FecundityReduction must be scalar, length nodes, length Nstages, nodes x Ntimesteps, or nodes x Nstages x Ntimesteps")
}
if(any(!is.finite(FecundityReduction)) || any(FecundityReduction < 0) || any(FecundityReduction > 1))
  stop("FecundityReduction values must be between 0 and 1")

# Resolve fecundity reduction to node x stage values.
ResolveFecundityReductionTM <- function(timestep) {
  if(is.null(dim(FecundityReduction))) {
    if(length(FecundityReduction) == 1)
      return(matrix(FecundityReduction, nrow = n_nodes, ncol = Nstages))
    if(length(FecundityReduction) == n_nodes)
      return(matrix(rep(FecundityReduction, Nstages), nrow = n_nodes, ncol = Nstages))
    return(matrix(rep(FecundityReduction, each = n_nodes), nrow = n_nodes, ncol = Nstages))
  }
  if(length(dim(FecundityReduction)) == 2)
    return(matrix(rep(FecundityReduction[,timestep], Nstages), nrow = n_nodes, ncol = Nstages))
  FecundityReduction[,,timestep]
}

# Assign standard deviation value to management in no value provided
if(is.null(ManageSD) == T)
	ManageSD = mean(ManageProb)/10
if(is.null(SpreadReductionSD) == T)
	SpreadReductionSD = (1-mean(SpreadReduction))/10

n_nodes <- nrow(SDDprob)


# --- Detection SD ---
if (is.null(DetectionSD)) {
  if (is.matrix(DetectionProb)) {
    # Calculate mean per stage (column)
    stage_means <- colMeans(DetectionProb, na.rm = TRUE)
    # Replicate per node
    DetectionSD <- matrix(stage_means / 10, nrow = n_nodes, ncol = ncol(DetectionProb), byrow = TRUE)
  } else {
    # If scalar or vector, repeat across all nodes
    DetectionSD <- matrix(mean(DetectionProb, na.rm = TRUE) / 10,
                          nrow = n_nodes, ncol = Nstages)
  }
}

# --- Information-triggered detection SD / shape validation ---
.InfoTriggeredNodeShapeOK <- function(x, allow_time=TRUE) {
  d <- dim(x)
  if(is.null(d)) return(length(x) == 1L || length(x) == Nstages)
  if(length(d) == 2L) return(all(d == c(n_nodes,Nstages)))
  if(allow_time && length(d) == 3L) return(all(d == c(n_nodes,Nstages,Ntimesteps)))
  FALSE
}
if(!.InfoTriggeredNodeShapeOK(InfoTriggeredDetectionProb, TRUE))
  stop("InfoTriggeredDetectionProb must be scalar, length Nstages, nodes x stages, or nodes x stages x Ntimesteps")
if(any(!is.finite(InfoTriggeredDetectionProb)) || any(InfoTriggeredDetectionProb < 0) || any(InfoTriggeredDetectionProb > 1))
  stop("InfoTriggeredDetectionProb values must be between 0 and 1")
if(is.null(InfoTriggeredDetectionSD)) {
  if(is.matrix(InfoTriggeredDetectionProb)) {
    stage_means <- colMeans(InfoTriggeredDetectionProb, na.rm=TRUE)
    InfoTriggeredDetectionSD <- matrix(stage_means/10, nrow=n_nodes, ncol=Nstages, byrow=TRUE)
  } else {
    InfoTriggeredDetectionSD <- matrix(mean(InfoTriggeredDetectionProb, na.rm=TRUE)/10, nrow=n_nodes, ncol=Nstages)
  }
}
if(!.InfoTriggeredNodeShapeOK(InfoTriggeredDetectionSD, FALSE))
  stop("InfoTriggeredDetectionSD must be scalar, length Nstages, or nodes x stages")
if(any(!is.finite(InfoTriggeredDetectionSD)) || any(InfoTriggeredDetectionSD < 0))
  stop("InfoTriggeredDetectionSD must contain finite non-negative values")
UseInfoTriggeredSurveillance <- any(InfoTriggeredDetectionProb != 0) || any(InfoTriggeredDetectionSD != 0)

.ResolveInfoTriggeredNodeMatrix <- function(x, timestep, name) {
  d <- dim(x)
  if(is.null(d)) {
    if(length(x) == 1L) return(matrix(x, nrow=n_nodes, ncol=Nstages))
    if(length(x) == Nstages) return(matrix(rep(x, each=n_nodes), nrow=n_nodes, ncol=Nstages))
    stop(name, " must be scalar or length Nstages")
  }
  if(length(d) == 2L) return(matrix(x, nrow=n_nodes, ncol=Nstages))
  if(length(d) == 3L) return(matrix(x[,,timestep], nrow=n_nodes, ncol=Nstages))
  stop(name, " has unsupported dimensions")
}

# --- Mortality SD ---
if (is.null(MortalitySD)) {
  if (is.matrix(MortalityProb)) {
    # Calculate mean per stage (column)
    stage_means <- colMeans(MortalityProb, na.rm = TRUE)
    # Replicate per node
    MortalitySD <- matrix(stage_means / 10, nrow = n_nodes, ncol = ncol(MortalityProb), byrow = TRUE)
  } else {
    # If scalar or vector, repeat across all nodes
    MortalitySD <- matrix(mean(MortalityProb, na.rm = TRUE) / 10,
                          nrow = n_nodes, ncol = Nstages)
  }
}




###########################################################
### Start of simulation
###########################################################
    

# ---------------------------------------------------------------------------
# Run independent stochastic vertebrate histories.
# ---------------------------------------------------------------------------
for (perm in 1:Nperm) 
{ 
# Initialise host abundance from explicit starting values or sampled invasion settings.
n_nodes <- nrow(SDDprob)
InitBio <- matrix(0, n_nodes, Nstages)

if (is.matrix(InitialPopulation) &&
    nrow(InitialPopulation) == n_nodes &&
    ncol(InitialPopulation) == Nstages) {
  
  InitBio <- InitialPopulation
  
} else {
  
  risk <- if (is.matrix(InvasionRisk) && nrow(InvasionRisk) == n_nodes)
    InvasionRisk[,1] else if (length(InvasionRisk) == n_nodes) InvasionRisk else NULL
  
  if (!is.na(InitBioP)) {
    Infested <- sample.int(n_nodes, ceiling(n_nodes * InitBioP), prob = risk)
  } else if (!is.null(risk)) {
    Infested <- which(rbinom(n_nodes, 1, risk) == 1)
  } else {
    Infested <- integer(0)
  }
  
  InitBio[Infested,1] <- if (is.na(IncursionStartPop)) 1 else IncursionStartPop
}

# --- Ensure initial population (weighted) does not exceed carrying capacity ---

if (exists("Weights") && length(Weights) == Nstages) {
  # Calculate weighted population per node
  weighted_pop <- InitBio %*% Weights  # (n_nodes x 1)
  
  # Identify nodes exceeding K
  overcap <- which(weighted_pop > NodeK)
  
  if (length(overcap) > 0) {
    # Scale down all stage values proportionally
    scale_factor <- NodeK[overcap] / weighted_pop[overcap]
    InitBio[overcap, ] <- InitBio[overcap, , drop = FALSE] * scale_factor
  }
  
} else {
  # Fallback: simple elementwise comparison (if Weights missing)
  for (i in seq_len(nrow(InitBio))) {
    if (any(InitBio[i, ] > NodeK[i])) {
      # Clip all stages proportionally to match K
      total_pop <- sum(InitBio[i, ])
      if (total_pop > 0) {
        InitBio[i, ] <- InitBio[i, ] * (NodeK[i] / total_pop)
      }
    }
  }
}

# Ensure all stage populations are integers
InitBio <- floor(InitBio)

# Set the working host abundance and initialise pathogen state when present.
N <- InitBio
if(sum(N) == 0 && OngoingExternalInvasion == F)
  warning("No initial populations and no future external invasions")

# Initialise response information independently of true host abundance.
# Without an initial information specification, all nodes start uninformed.
# Robust InitInfo setup
n_nodes <- nrow(SDDprob)
HasExplicitInitialInfo <- !(length(InitialInfo)==1L && is.na(InitialInfo))
if(HasExplicitInitialInfo) {
  if(!is.numeric(InitialInfo) || length(InitialInfo)!=n_nodes || any(!is.finite(InitialInfo)) || any(!InitialInfo %in% c(0,1)))
    stop("InitialInfo must be scalar NA or a binary vector of length nodes")
  InitInfo <- as.integer(InitialInfo)
} else {
  InitInfo <- integer(n_nodes)
  ExternalInfoMissing <- length(ExternalInfoProb)==1L && is.na(ExternalInfoProb)
  ext0 <- if(ExternalInfoMissing) { if(!is.na(InitInfoP)) rep(1,n_nodes) else rep(0,n_nodes) } else if(is.matrix(ExternalInfoProb)) ExternalInfoProb[,1] else rep_len(ExternalInfoProb,n_nodes)
  if(any(!is.finite(ext0)) || any(ext0<0 | ext0>1)) stop("ExternalInfoProb values must be probabilities between 0 and 1")
  if(!is.na(InitInfoP)) {
    if(!is.numeric(InitInfoP) || length(InitInfoP)!=1L || !is.finite(InitInfoP) || InitInfoP<0 || InitInfoP>1) stop("InitInfoP must be NA or one probability between 0 and 1")
    n_sample <- min(n_nodes,ceiling(n_nodes*InitInfoP))
    if(n_sample>0L) {
      weights <- ext0
      if(sum(weights)<=0) weights <- rep(1,n_nodes)
      Info <- sample.int(n_nodes,size=n_sample,prob=weights)
      InitInfo[Info] <- 1L
    }
  } else if(any(ext0>0)) {
    InitInfo <- rbinom(n_nodes,size=1,prob=ext0)
  }
}


# --- Determine dimensions ---
n_nodes <- nrow(SDDprob)

# Randomly assign detection probabilities
# --- Case 1: Single value or vector (stages only) ---
if (!is.array(DetectionProb) && (length(DetectionProb) == 1 || length(DetectionProb) == Nstages)) {
  
  NodeDetectionProb <- array(NA, dim = c(n_nodes, Nstages, Ntimesteps))
  
  for (s in 1:Nstages) {
    for (t in 1:Ntimesteps) {
      NodeDetectionProb[, s, t] <- pmax(0, pmin(1,
                                                rnorm(n_nodes,
                                                      mean = if (length(DetectionProb) == 1) DetectionProb else DetectionProb[s],
                                                      sd   = if (is.matrix(DetectionSD)) DetectionSD[, s]   # <--- per node & stage
                                                      else if (length(DetectionSD) == 1) DetectionSD
                                                      else DetectionSD[s]
                                                )
      ))
    }
  }
}


# --- Case 2: Matrix (nodes × stages) ---
if (is.matrix(DetectionProb) && all(dim(DetectionProb) == c(n_nodes, Nstages))) {
  
  NodeDetectionProb <- array(NA, dim = c(n_nodes, Nstages, Ntimesteps))
  
  for (t in 1:Ntimesteps) {
    NodeDetectionProb[, , t] <- pmax(0, pmin(1,
                                             matrix(
                                               rnorm(n_nodes * Nstages,
                                                     mean = as.vector(DetectionProb),
                                                     sd   = as.vector(DetectionSD)),
                                               nrow = n_nodes, ncol = Nstages)
    ))
  }
}

# --- Case 3: 3D array (nodes × stages × timesteps) ---
if (length(dim(DetectionProb)) == 3 && all(dim(DetectionProb)[1:2] == c(n_nodes, Nstages))) {
  
  NodeDetectionProb <- array(NA, dim = dim(DetectionProb))
  
  for (t in 1:Ntimesteps) {
    NodeDetectionProb[, , t] <- pmax(0, pmin(1,
                                             matrix(
                                               rnorm(n_nodes * Nstages,
                                                     mean = as.vector(DetectionProb[, , t]),
                                                     sd   = as.vector(DetectionSD)),
                                               nrow = n_nodes, ncol = Nstages)
    ))
  }
}

# Realised InfoTriggered detection probabilities for this permutation.
NodeInfoTriggeredDetectionProb <- array(0, dim=c(n_nodes,Nstages,Ntimesteps))
if(UseInfoTriggeredSurveillance) {
  for(t in seq_len(Ntimesteps)) {
    mu <- .ResolveInfoTriggeredNodeMatrix(InfoTriggeredDetectionProb, t, "InfoTriggeredDetectionProb")
    sdmat <- .ResolveInfoTriggeredNodeMatrix(InfoTriggeredDetectionSD, t, "InfoTriggeredDetectionSD")
    NodeInfoTriggeredDetectionProb[,,t] <- pmax(0,pmin(1,
      matrix(rnorm(n_nodes*Nstages, mean=as.vector(mu), sd=as.vector(sdmat)),
             nrow=n_nodes,ncol=Nstages)))
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
# --- Determine dimensions ---
n_nodes <- nrow(SDDprob)

# --- Case 1: Single value or vector (stages only) ---
if (!is.array(MortalityProb) && (length(MortalityProb) == 1 || length(MortalityProb) == Nstages)) {
  
  NodeMortalityProb <- array(NA, dim = c(n_nodes, Nstages, Ntimesteps))
  
  for (s in 1:Nstages) {
    for (t in 1:Ntimesteps) {
      NodeMortalityProb[, s, t] <- pmax(0, pmin(1,
                                                rnorm(n_nodes,
                                                      mean = if (length(MortalityProb) == 1) MortalityProb else MortalityProb[s],
                                                      sd   = if (is.matrix(MortalitySD)) MortalitySD[, s]   # <--- per node & stage
                                                      else if (length(MortalitySD) == 1) MortalitySD
                                                      else MortalitySD[s]
                                                )
      ))
    }
  }
}


# --- Case 2: Matrix (nodes × stages) ---
if (is.matrix(MortalityProb) && all(dim(MortalityProb) == c(n_nodes, Nstages))) {
  
  NodeMortalityProb <- array(NA, dim = c(n_nodes, Nstages, Ntimesteps))
  
  for (t in 1:Ntimesteps) {
    NodeMortalityProb[, , t] <- pmax(0, pmin(1,
                                             matrix(
                                               rnorm(n_nodes * Nstages,
                                                     mean = as.vector(MortalityProb),
                                                     sd   = as.vector(MortalitySD)),
                                               nrow = n_nodes, ncol = Nstages)
    ))
  }
}

# --- Case 3: 3D array (nodes × stages × timesteps) ---
if (length(dim(MortalityProb)) == 3 && all(dim(MortalityProb)[1:2] == c(n_nodes, Nstages))) {
  
  NodeMortalityProb <- array(NA, dim = dim(MortalityProb))
  
  for (t in 1:Ntimesteps) {
    NodeMortalityProb[, , t] <- pmax(0, pmin(1,
                                             matrix(
                                               rnorm(n_nodes * Nstages,
                                                     mean = as.vector(MortalityProb[, , t]),
                                                     sd   = as.vector(MortalitySD)),
                                               nrow = n_nodes, ncol = Nstages)
    ))
  }
}






# Record current host/pest presence from starting abundance.
Invaded <- ifelse(rowSums(InitBio) > 0, 1, 0)

# Apply initial surveillance to the starting host population.
# Accepted detections can add response information.




# Calculate probability of detection per node
prob_detect <- 1 - apply((1 - matrix(NodeDetectionProb[,,1],nrow=n_nodes,ncol=Nstages))^InitBio, 1, prod)

# Draw initial detection (0/1) per node
InitDetection <- rbinom(n = nrow(InitBio), size = 1, prob = prob_detect)


# Add detections to nodes which already have info (e.g. pre-emptive control and hygiene measures)
InitInfo[InitInfo == 0] = InitDetection[InitInfo == 0]

# Set the working information state for the simulation.
HaveInfo = InitInfo

# Track the most recent timestep with known local presence
LastKnownPresence = rep(NA,nrow(SDDprob))
if(UseInfoPersistence == T)
  {
  InitialKnownPresence = which(HaveInfo == 1)
  if(length(InitialKnownPresence) > 0)
    LastKnownPresence[InitialKnownPresence] = 0
  }


  # run simulation

  # ---------------------------------------------------------------------------
  # Advance stage dynamics and vertebrate response processes through time.
  # ---------------------------------------------------------------------------
  for (timestep in 1:Ntimesteps) 
    { 
    # Print progress
    if(DoProgress) cat("\r", "Realisation ", perm, "Timestep ", timestep, "...")
 
    # Resolve short- and long-distance connectivity for the current timestep.
    NodeSDDprob <- .SliceConnectivityTM(SDDprob,timestep)
    NodeLDDprob <- .SliceConnectivityTM(LDDprob,timestep)
    NodeTransitionSDDprob <- ResolveTransitionMovement(TransitionSDDprob,timestep)
    NodeTransitionLDDprob <- ResolveTransitionMovement(TransitionLDDprob,timestep)
    NodeFecundityReduction <- ResolveFecundityReductionTM(timestep)
# Allow for variation in recruit establishment through time
    if(is.matrix(EnvEstabProb) == T)
      NodeEnvEstabProb <- EnvEstabProb[,timestep]
    
    NodeTransition <- .ResolveTransitionTM(Transition,timestep)

    
    
    # Resolve carrying capacity for the current timestep.
    if(is.matrix(K) == TRUE)
      {
      K_is_0 <- K[,timestep]<=0
      inv_K <- 1 / sum(K[,timestep])
      NodeK = K[,timestep] 
      }  
    
    # If seedbank carrying capacity provided as matrix assign values for relevant timestep
    if (is.matrix(SeedbankK) == TRUE)
      {
      NodeSeedbankK <- SeedbankK[, timestep]
      }
    
    
    if(is.matrix(PropaguleEstablishment) == TRUE)
      NodePropaguleEstablishment = PropaguleEstablishment[,timestep]

  ###########################################################################
  ### Pre-deployed vertebrate control, independent of information status
  ###########################################################################
  NodeControlFecundityReduction <- matrix(0, nrow = nrow(N), ncol = Nstages)
  ControlDetectedNodes <- rep(FALSE, nrow(N))
  # Retain the exact pre-control population exposed to routine devices/area control.
  RoutineControlObservationAbundanceResults[,,timestep,perm] <- N
  if(!is.null(ControlModule) && sum(N) > 0)
    {
    control_context <- list(
      transition = NodeTransition, ModelName = ModelName, Nstages = Nstages,
      Weights = Weights, K = NodeK, phase = "pre_response_control"
    )
    NodeHomeRange <- .iv_home_range_node(
      HomeRangeModule, N, timestep, perm, control_context
    )
    ce <- .iv_node_control(
      ControlModule, N, NodeHomeRange, timestep, perm, control_context
    )
    NodeControlFecundityReduction <- ce$fecundity_reduction
    ControlCostResults[timestep,perm] <- ce$cost
    RoutineControlDetectionProbabilityResults[,,timestep,perm] <- ce$detect_prob

    # Whole-node detection probability implied by per-individual stage exposure.
    ControlProbDetectNode <- 1 - apply((1 - ce$detect_prob)^N, 1, prod)
    ControlDetectedNodes <- rbinom(nrow(N), size = 1, prob = ControlProbDetectNode) == 1
    ControlDetectionResults[,timestep,perm] <- as.integer(ControlDetectedNodes)

    NpreControl <- matrix(0, nrow = nrow(N), ncol = Nstages)
    for(s in seq_len(Nstages))
      NpreControl[,s] <- rbinom(nrow(N), size = N[,s], prob = 1 - ce$kill_prob[,s])
    ControlDeathResults[,,timestep,perm] <- N - NpreControl
    N <- NpreControl
    }

  # Draw node-level management-adoption probabilities where information is available.
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
  
  # Use current information to activate management at each node.
  Managing = rbinom(1:nrow(SDDprob),size = 1,prob = NodeManageProb*HaveInfo)
  
  Managing = Managing*HaveInfo
  
  # Identify occupied nodes currently known to the response system.
  Detected = Invaded*HaveInfo
  
  # Apply management
  # N: current population matrix (nodes x stages)
  # Managing: vector of management adoption per node (0/1) for current timestep
  # NodeMortalityProb: array (nodes x stages x timesteps)
  # timestep: current timestep index
  
  # Apply management-driven mortality
  N0 <- matrix(NA, nrow = nrow(N), ncol = ncol(N))
  for (s in 1:Nstages) {
    N0[, s] <- rbinom(n = nrow(N), size = N[, s], prob = 1 - (NodeMortalityProb[, s,timestep] * Managing))
  }
  
  # Track known local presence from actual management mortality
  if(UseInfoPersistence == T)
    {
    KnownPresence = which(rowSums(N-N0) > 0)
    if(length(KnownPresence) > 0)
      LastKnownPresence[KnownPresence] = timestep
    }
  
  # Update population matrix
  if (sum(N0) <= 0) {
    N <- N0
  } else if (sum(N0) > 0) {
    N <- N0
  }
  
  if(sum(N0)>0 ) 
  {
  NodeBirthSpec <- NULL
  if(!is.null(BirthModule))
    {
    response_fec_mult <- 1 - NodeFecundityReduction * Managing
    birth_context <- list(
      fecundity_multiplier = response_fec_mult * (1 - NodeControlFecundityReduction),
      response_fecundity_reduction = NodeFecundityReduction,
      control_fecundity_reduction = NodeControlFecundityReduction,
      managing = Managing,
      home_range = .iv_home_range_node(
        HomeRangeModule, N0, timestep, perm,
        list(transition = NodeTransition, K = NodeK, phase = "birth")
      ),
      K = NodeK, ModelName = ModelName, Nstages = Nstages, Weights = Weights
    )
    NodeBirthSpec <- .iv_node_birth(
      BirthModule, N0, NodeTransition, timestep, perm, birth_context
    )
    }
      
  LocalDynamicsArgs <- list(nodetransition = NodeTransition, weights = Weights, sddprob = NodeSDDprob,
                            nodeenvestabprob = NodeEnvEstabProb, n0 = N0, lddprob = NodeLDDprob,
                            lddrate = LDDrate, nodeK = NodeK, node.seedbankK = NodeSeedbankK,
                            nodepropaguleestablishment = NodePropaguleEstablishment,
                            nodespreadreduction = NodeSpreadReduction, managing = Managing,
                            MaxInteger = MaxInteger, BlockedTransitionMortality = BlockedTransitionMortality,
                            DispersalDensityFactor = DispersalDensityFactor)
  LocalDynamicsFormals <- names(formals(LocalDynamics))
  LocalDynamicsAcceptsFecundityReduction <- "nodefecundityreduction" %in% LocalDynamicsFormals || "..." %in% LocalDynamicsFormals
  if(LocalDynamicsAcceptsFecundityReduction)
    LocalDynamicsArgs$nodefecundityreduction <- NodeFecundityReduction
  else if(any(NodeFecundityReduction * Managing > 0))
    stop("Custom LocalDynamics must accept a 'nodefecundityreduction' argument (or ...) when FecundityReduction is active")

  LocalDynamicsAcceptsControlFecundity <-
    "nodecontrolfecundityreduction" %in% LocalDynamicsFormals || "..." %in% LocalDynamicsFormals
  if(LocalDynamicsAcceptsControlFecundity)
    LocalDynamicsArgs$nodecontrolfecundityreduction <- NodeControlFecundityReduction
  else if(any(NodeControlFecundityReduction > 0))
    stop("Custom LocalDynamics must accept 'nodecontrolfecundityreduction' (or ...) when vertebrate control reduces fecundity")

  if(!is.null(NodeBirthSpec))
    {
    LocalDynamicsAcceptsBirth <- "nodebirthmean" %in% LocalDynamicsFormals || "..." %in% LocalDynamicsFormals
    if(!LocalDynamicsAcceptsBirth)
      stop("Custom LocalDynamics must accept 'nodebirthmean' (or ...) when Vertebrate$Birth is active")
    LocalDynamicsArgs$nodebirthmean <- NodeBirthSpec$mean
    if(!is.null(NodeBirthSpec$mother_counts))
      {
      if(!("nodebirthmothers" %in% LocalDynamicsFormals || "..." %in% LocalDynamicsFormals))
        stop("Custom LocalDynamics must accept 'nodebirthmothers' (or ...) when Birth returns mother_counts")
      LocalDynamicsArgs$nodebirthmothers <- NodeBirthSpec$mother_counts
      }
    }

  LocalDynamicsAcceptsTransitionMovement <-
    "..." %in% LocalDynamicsFormals ||
    all(c("transition_sddprob", "transition_lddprob", "transition_lddrate") %in% LocalDynamicsFormals)
  if(LocalDynamicsAcceptsTransitionMovement) {
    LocalDynamicsArgs$transition_sddprob <- NodeTransitionSDDprob
    LocalDynamicsArgs$transition_lddprob <- NodeTransitionLDDprob
    LocalDynamicsArgs$transition_lddrate <- TransitionLDDrate
  } else if(TransitionMovementConfigured) {
    stop("Custom LocalDynamics must accept 'transition_sddprob', 'transition_lddprob' and 'transition_lddrate' arguments (or ...) when transition movement is active")
  }
  N <- do.call(LocalDynamics, LocalDynamicsArgs)
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
 
 
 # Add invasion resulting from colonisation from external sources
 if(OngoingExternalInvasion == T)
  {
  if(is.matrix(InvasionRisk) == F)
   ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk)
  if(is.matrix(InvasionRisk) == T)
   ExternalInvasion = rbinom(1:nrow(SDDprob),size = 1,prob = InvasionRisk[,timestep])
  Invaded[Invaded == 0] = ExternalInvasion[Invaded==0]
  if(is.na(IncursionStartPop) == T) 
	  N[,1] = N[,1]+ExternalInvasion
  if(is.na(IncursionStartPop) == F) 
	  N[,1] = N[,1]+ExternalInvasion*IncursionStartPop
 
  }

  # Optional aggregate social/contact/disease-state interaction. The hook can
  # redistribute the existing node x class population but cannot change array
  # dimensions; births and deaths belong in their dedicated processes.
  if(!is.null(InteractionModule) && sum(N) > 0)
    {
    interaction_context <- list(
      transition = NodeTransition, K = NodeK, ModelName = ModelName,
      Nstages = Nstages, Weights = Weights, phase = "interaction"
    )
    InteractionHomeRange <- .iv_home_range_node(
      HomeRangeModule, N, timestep, perm, interaction_context
    )
    N <- .iv_node_interaction(
      InteractionModule, N, InteractionHomeRange, timestep, perm, interaction_context
    )
    }
  
  # Ensure all stage populations are integers
  N <- floor(N)
 
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
 Invaded = ifelse(rowSums(N)>0,1,0)
 
 # Record nodes adopting management
 ManagingResults[,timestep,perm] = Managing
  
 # Record infested nodes
 InvasionResults[,timestep,perm] = Invaded

 # Record populations
weighted_population <- if (is.matrix(Weights)) {
  rowSums(N[, 2:Nstages, drop = FALSE] *
            Weights[, 2:Nstages, drop = FALSE])
} else {
  as.numeric(N[, 2:Nstages, drop = FALSE] %*% Weights[2:Nstages])
}
PopulationResults[, timestep, perm] <- weighted_population
 
 # Record stage populations
 PopulationStageResults[,,timestep,perm] = N

 # Freeze information available before this round's surveillance.
 # Routine control detections occurred earlier but are deliberately not registered
 # until after this snapshot, so they cannot activate InfoTriggered surveillance
 # retrospectively in the same timestep.
 InfoBeforeSurveillance <- as.integer(HaveInfo != 0)
 InformationStateBeforeSurveillanceResults[,timestep,perm] <- InfoBeforeSurveillance

 # Record the realised per-individual observation probabilities used this round.
 BackgroundDetectionProbabilityResults[,,timestep,perm] <- matrix(NodeDetectionProb[,,timestep],nrow=n_nodes,ncol=Nstages)
 InfoTriggeredDetectionProbabilityResults[,,timestep,perm] <- matrix(NodeInfoTriggeredDetectionProb[,,timestep],nrow=n_nodes,ncol=Nstages)

 # Ordinary Background surveillance after population dynamics.
 DetectionProbPerStage <- 1 - (1 - matrix(NodeDetectionProb[,,timestep],nrow=n_nodes,ncol=Nstages))^N
 ProbDetectNode <- 1 - apply(1 - DetectionProbPerStage, 1, prod)
 BackgroundDetection <- rbinom(n=nrow(SDDprob), size=1, prob=ProbDetectNode)

 # Additional targeted surveillance is gated only by information that existed
 # before the current observation results.
 InfoTriggeredDetection <- integer(nrow(SDDprob))
 if(UseInfoTriggeredSurveillance) {
   InfoDetectionProbPerStage <- 1 - (1 - matrix(NodeInfoTriggeredDetectionProb[,,timestep],nrow=n_nodes,ncol=Nstages))^N
   ProbInfoDetectNode <- 1 - apply(1 - InfoDetectionProbPerStage, 1, prod)
   InfoTriggeredDetection <- rbinom(
     n=nrow(SDDprob), size=1, prob=ProbInfoDetectNode*InfoBeforeSurveillance)
 }

 # Routine control and ordinary survey detections are both Background evidence.
 # The combined field stores event counts, so a node detected by both pathways
 # can contribute two observation events in one timestep.
 BackgroundSurveillanceDetectedResults[,timestep,perm] <- BackgroundDetection
 BackgroundDetectedResults[,timestep,perm] <- as.integer(ControlDetectedNodes) + BackgroundDetection
 InfoTriggeredDetectedResults[,timestep,perm] <- InfoTriggeredDetection
 HostDetectionEvidence <- pmax(as.integer(ControlDetectedNodes), BackgroundDetection, InfoTriggeredDetection)

 # All observation pathways update response information only after surveillance.
 if(UseInfoPersistence == T) {
   KnownPresence <- which(HostDetectionEvidence == 1)
   if(length(KnownPresence) > 0) LastKnownPresence[KnownPresence] <- timestep
 }
 HaveInfo[HaveInfo == 0] <- HostDetectionEvidence[HaveInfo == 0]
 HaveInfoResults[,timestep,perm] <- HaveInfo

 # Legacy DetectedResults remains the persistent known-present state.
 DetectedResults[,timestep,perm] = HaveInfo*Invaded 
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
if(SaveResults) {
  saveRDS(ManagingResults, paste0(FileNameStem,"InfoLargeOut.rds"))
  saveRDS(PopulationResults, paste0(FileNameStem,"PopulationLargeOut.rds"))
  saveRDS(PopulationStageResults, paste0(FileNameStem,"PopulationStageLargeOut.rds"))
  saveRDS(InvasionResults, paste0(FileNameStem,"InvasionLargeOut.rds"))
  saveRDS(DetectedResults, paste0(FileNameStem,"DetectedLargeOut.rds"))
  saveRDS(BackgroundDetectedResults, paste0(FileNameStem,"BackgroundDetectedLargeOut.rds"))
  saveRDS(BackgroundSurveillanceDetectedResults, paste0(FileNameStem,"BackgroundSurveillanceDetectedLargeOut.rds"))
  saveRDS(InfoTriggeredDetectedResults, paste0(FileNameStem,"InfoTriggeredDetectedLargeOut.rds"))
  saveRDS(InformationStateBeforeSurveillanceResults, paste0(FileNameStem,"InformationStateBeforeSurveillanceLargeOut.rds"))
  saveRDS(HaveInfoResults, paste0(FileNameStem,"HaveInfoLargeOut.rds"))
  saveRDS(BackgroundDetectionProbabilityResults, paste0(FileNameStem,"BackgroundDetectionProbabilityLargeOut.rds"))
  saveRDS(InfoTriggeredDetectionProbabilityResults, paste0(FileNameStem,"InfoTriggeredDetectionProbabilityLargeOut.rds"))
  saveRDS(RoutineControlObservationAbundanceResults, paste0(FileNameStem,"RoutineControlObservationAbundanceLargeOut.rds"))
  saveRDS(RoutineControlDetectionProbabilityResults, paste0(FileNameStem,"RoutineControlDetectionProbabilityLargeOut.rds"))
  saveRDS(ControlDeathResults, paste0(FileNameStem,"VertebrateControlDeaths.rds"))
  saveRDS(ControlDetectionResults, paste0(FileNameStem,"VertebrateControlDetections.rds"))
  saveRDS(ControlCostResults, paste0(FileNameStem,"VertebrateControlCost.rds"))
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
if(SaveResults) saveRDS(InvasionProb, paste0(FileNameStem,"InvasionProb.rds"))
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

# Return a structured object as well as retaining the parent function's RDS
# outputs. This is additive and does not alter the legacy saved results.
out <- list(
  ModelName = ModelName,
  # Legacy vertebrate result aliases
  Population = PopulationResults,
  PopulationStage = PopulationStageResults,
  Invasion = InvasionResults,
  Detection = DetectedResults,
  Managing = ManagingResults,
  ControlDeaths = ControlDeathResults,
  ControlDetections = ControlDetectionResults,
  ControlCost = ControlCostResults,
  InvasionProbability = InvasionProb,
  # Canonical INApest / PoA result contract
  PopulationResults = PopulationResults,
  PopulationStageResults = PopulationStageResults,
  InvasionResults = InvasionResults,
  DetectedResults = DetectedResults,
  ManagingResults = ManagingResults,
  BackgroundDetectedResults = BackgroundDetectedResults,
  BackgroundSurveillanceDetectedResults = BackgroundSurveillanceDetectedResults,
  InfoTriggeredDetectedResults = InfoTriggeredDetectedResults,
  InformationStateBeforeSurveillanceResults = InformationStateBeforeSurveillanceResults,
  HaveInfoResults = HaveInfoResults,
  BackgroundDetectionProbabilityResults = BackgroundDetectionProbabilityResults,
  InfoTriggeredDetectionProbabilityResults = InfoTriggeredDetectionProbabilityResults,
  RoutineControlObservationAbundanceResults = RoutineControlObservationAbundanceResults,
  RoutineControlDetectionProbabilityResults = RoutineControlDetectionProbabilityResults,
  ControlDetectionResults = ControlDetectionResults,
  ControlDeathResults = ControlDeathResults,
  ControlCostResults = ControlCostResults
)
class(out) <- c("INApestVertebrateNode", "list")
invisible(out)
}


################################################################
################################################################
### End of function
################################################################
################################################################
