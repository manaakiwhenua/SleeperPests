###############################################################################
### INApestVertebrateNode
### Standalone vertebrate-specialist node extension of the current
### INApestMetaTransitionMatrix engine. Base R; no package dependency.
### Public extension is one argument: Vertebrate = list(Birth, HomeRange,
### Control, Interaction). NULL preserves parent biological semantics.
###############################################################################

###############################################################################
### Shared helpers for INApest vertebrate extensions
###
### Public vertebrate functions add one top-level argument only:
###   Vertebrate = list(
###     Birth = NULL,
###     HomeRange = NULL,
###     Control = NULL,
###     Interaction = NULL
###   )
###
### These helpers are deliberately private (dot-prefixed). They use base R only.
###############################################################################

.iv_clip01 <- function(x) pmin(1, pmax(0, x))

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

.iv_module <- function(Vertebrate, name) {
  if (is.null(Vertebrate) || is.null(Vertebrate[[name]])) return(NULL)
  Vertebrate[[name]]
}

.iv_call_hook <- function(fun, args, name = "vertebrate hook") {
  if (!is.function(fun)) stop(name, " must be a function.")
  fm <- names(formals(fun))
  if (is.null(fm) || "..." %in% fm) return(do.call(fun, args))
  do.call(fun, args[intersect(names(args), fm)])
}

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

.iv_point_effects_empty <- function(points) {
  data.frame(
    id = points$id,
    kill_prob = rep(0, nrow(points)),
    detect_prob = rep(0, nrow(points)),
    fecundity_reduction = rep(0, nrow(points)),
    stringsAsFactors = FALSE
  )
}

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

.iv_combine_point_effects <- function(a, b) {
  out <- a
  out$kill_prob <- 1 - (1 - a$kill_prob) * (1 - b$kill_prob)
  out$detect_prob <- 1 - (1 - a$detect_prob) * (1 - b$detect_prob)
  out$fecundity_reduction <- 1 - (1 - a$fecundity_reduction) * (1 - b$fecundity_reduction)
  out
}

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

.iv_node_effects_empty <- function(n_nodes, n_stages) {
  list(
    kill_prob = matrix(0, n_nodes, n_stages),
    detect_prob = matrix(0, n_nodes, n_stages),
    fecundity_reduction = matrix(0, n_nodes, n_stages),
    cost = 0
  )
}

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

.iv_combine_node_effects <- function(a, b) {
  list(
    kill_prob = 1 - (1 - a$kill_prob) * (1 - b$kill_prob),
    detect_prob = 1 - (1 - a$detect_prob) * (1 - b$detect_prob),
    fecundity_reduction = 1 - (1 - a$fecundity_reduction) * (1 - b$fecundity_reduction),
    cost = a$cost + b$cost
  )
}

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
  
  Pin <- numeric(n_pops)
  Qin <- numeric(n_pops)
  
  if (any(propagules > 0)) {
    # Self-mediated spread.
    sdd_sources <- propagules * (1 - lddrate)
    sdd_destination_weights <- as.numeric(sdd_sources %*% sddprob)
    total_p <- sum(sdd_destination_weights)
    
    if (total_p > 0) {
      Pin <- if (floor(total_p) < MaxInteger) {
        as.numeric(
          rmultinom(
            n = 1,
            size = floor(total_p),
            prob = sdd_destination_weights
          )
        )
      } else {
        # Algebraically equivalent to colSums(sweep(...)), without
        # allocating another n_pops x n_pops matrix.
        as.numeric(floor(sdd_sources) %*% sddprob)
      }
    }
    
    # Human-mediated spread. Management is applied at each source before
    # calculating both the destination weights and the multinomial size.
    if (is.matrix(lddprob)) {
      spread_reduction <- pmin(
        1,
        pmax(
          0,
          rep_len(nodespreadreduction, n_pops) *
            rep_len(managing, n_pops)
        )
      )
      
      ldd_sources <- propagules * lddrate * (1 - spread_reduction)
      ldd_destination_weights <- as.numeric(ldd_sources %*% lddprob)
      total_q <- sum(ldd_destination_weights)
      
      if (total_q > 0) {
        Qin <- if (floor(total_q) < MaxInteger) {
          as.numeric(
            rmultinom(
              n = 1,
              size = floor(total_q),
              prob = ldd_destination_weights
            )
          )
        } else {
          as.numeric(floor(ldd_sources) %*% lddprob)
        }
      }
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
ModelName = "INApestVertebrateNode", #Name for storing results to file 
Nperm,                  #Number of permutations per parameter combination
Ntimesteps,                 #Simulation duration timesteps can be any length of time
Nstages,               #Number of stages in transition matrix
Weights,               #Weight for converting stage populations to total populations
Transition,            #Transition matrix (N stages x N stages), list of matrices (length = N nodes)
                       #or 4D array (Nstages x N stages x N nodes x N timesteps)
LocalDynamics = .iv_local_dynamics_transition_matrix, #Local population growth, dispersal and management function; user-defined functions are supported
DetectionProb,          #Vector of Per-individual detection probability for each stage, or matrix of probabilities per stage per node (e.g. farm) 
                        #or 3D array of probabilities per stage per node per year (must be between 0 and 1)
DetectionSD = NULL, #Option to provide standard deviation for detection probability. Can be single number or vector (nodes)
ManageProb,             #Probability or vector of probabilities vector length nrow(SDDprob)of node adopting management upon detection
ManageSD = NULL, #Option to provide standard deviation for management probability. Can be single number or vector (nodes)
MortalityProb,           ##Vector of Per-individual mortality probability for each stage, or matrix of probabilities per stage per node (e.g. farm) 
                         #or 3D array of probabilities per stage per node per year (must be between 0 and 1)
MortalitySD = NULL, #Option to provide standard deviation for mortality probability. Can be single number or vector (nodes)
FecundityReduction = 0, #Proportional reduction in fecundity under management: scalar; vector (nodes or stages); matrix (nodes x timesteps); or array (nodes x stages x timesteps)
SpreadReduction,        #Reduction in dispersal probability when management adopted. Must be between 0 (no spread reduction) and 1 (complete prevention of spread). Can be single value or vector length nrow(SDDprob)
SpreadReductionSD = NULL, #Option to provide standard deviation for spread reduction probability can be single number or vector (nodes)
InitialPopulation = NA,        #matrix (nodes x stages) of population sizes at start of simulations
InitBioP = NA,		#Proportion of nodes infested at start of simulations
InvasionRisk = NA,           #Vector of probabilities of invasion from external sources
InitialInfo = NA,        #Vector or of nodes with information at start of simulations
InitInfoP = NA,		#Proportion of nodes with information at start of simulations
ExternalInfoProb = 0.0,           #Vector of probabilities of communication from external sources
InfoRetentionProb = 1,       #Probability that existing information is retained between timesteps. Can be single number, vector (nodes) or matrix (nodes x timesteps)
InfoPersistenceSteps = NA,    #Number of timesteps information persists after last known local presence. Can be single number, vector (nodes) or matrix (nodes x timesteps); NA uses InfoRetentionProb
EnvEstabProb = 1,           #Environmentally determined establishment probability. Can be single value, vector (nodes) or matrix (nodes x timesteps)
K,		       #Population carrying capacity - vector (nodes) or matrix (nodes x timesteps)
SeedbankK,        # Seedbank carrying capacity - vector (nodes) or matrix (nodes x timesteps)
PropaguleEstablishment, #Propagules establishment probability
IncursionStartPop=NA,      #option to set population size for new incursions
SDDprob,                   #Natural dispersal probability matrix, or 3D array (nodes x nodes x timesteps)
SEAM = 0,			#Option to provide socioeconomic adjacency matrix for information spread
LDDprob = NA,         #Option to provide long distance (human-mediated) dispersal matrix or 3D array (nodes x nodes x timesteps) instead of distance-independent dispesal rate
			      #e.g. could be weighted by law of human visitation or data on stock movements
LDDrate = 0,         #Proportion of available propagules entering LDD
TransitionSDDprob = NULL, #Optional movement during stage progression. Matrix applies to all transitions; list length Nstages-1 allows stage-specific NULL/matrix/3D-array entries
TransitionLDDprob = NULL, #Optional LDD movement during stage progression, with the same structure as TransitionSDDprob
TransitionLDDrate = 0,    #LDD fraction for a dispersing transition when both transition SDD and LDD matrices are supplied; scalar or length Nstages-1
DispersalDensityFactor = 0,
BlockedTransitionMortality = 0,
OngoingExternalInvasion = F,   ##Option to include ongoing invasion from external sources
OngoingExternalInfo = F,   ##Option to include ongoing communication from external sources
Vertebrate = NULL,         #Optional list(Birth, HomeRange, Control, Interaction)
OutputDir = NA,		      #Directory for storing results
DoPlots = TRUE,	     #Option to omit printing of line graphs.Default is to print.
InfoTriggeredDetectionProb = 0, #Additional per-individual stage-specific detection under pre-existing information
InfoTriggeredDetectionSD = NULL, #Optional uncertainty for information-triggered detection
SaveResults = TRUE, #Retain historical RDS side effects by default; PoA wrapper disables them
DoProgress = TRUE #Retain historical console progress by default; PoA wrapper disables it
)
{
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
###POTENTIAL ADDITIONS
###1) Make detection prob a function of population size. Could be based on individual detection prob so that DetectionProb = 1-(1-DPindividual)^N)
###   DPindividual could vary between nodes
###2) Allow provision of natural mortality rate to permit extinction of local populations (may happen in climates where R0 is very low?)
  
###Max integer for propagule dispersal using rmultinom
MaxInteger <- .Machine$integer.max  
  
###Allow SDD and LDD connectivity to vary through time
if(length(dim(SDDprob)) == 3 && (dim(SDDprob)[1] != dim(SDDprob)[2] || dim(SDDprob)[3] != Ntimesteps))
  stop("SDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")
if(length(dim(LDDprob)) == 3 && (dim(LDDprob)[1] != nrow(SDDprob) || dim(LDDprob)[2] != nrow(SDDprob) || dim(LDDprob)[3] != Ntimesteps))
  stop("LDDprob 3D array must have dimensions nodes x nodes x Ntimesteps")

###Optional movement while individuals progress between life-history stages.
###A bare matrix/3D array applies to every progression transition. A list of
###length Nstages-1 allows each source-stage transition to have its own matrix,
###time-varying array, or NULL (local-only transition).
n_nodes <- nrow(SDDprob)
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
        if(any(rowSums(z[,,tt,drop=FALSE][,,1]) > 1 + 1e-10))
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

ResolveTransitionMovement <- function(x, timestep) {
  resolve_one <- function(z) {
    if(is.null(z)) return(NULL)
    if(length(dim(z)) == 3L) return(z[,,timestep])
    z
  }
  if(is.null(x)) return(NULL)
  if(is.list(x)) return(lapply(x, resolve_one))
  resolve_one(x)
}

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
if(is.matrix(K) == FALSE)
{
K_is_0 <- K<=0
inv_K <- 1 / sum(K)
NodeK = K
}

###If carrying capacity provided as matrix assign values from first timestep for population initialisation
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

###Declare array tracking population size
###of individual nodes in each timestep of each realisation
PopulationResults = array(dim = c(nrow(SDDprob),Ntimesteps,Nperm))

###Declare array tracking population size
###of individual nodes in each timestep of each realisation
PopulationStageResults = array(dim = c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))

###Declare array tracking invasion status
###of individual nodes in each timestep of each realisation
InvasionResults = array(dim = c(nrow(SDDprob),Ntimesteps,Nperm))


###Declare array tracking detection status 
###of individual nodes in each timestep of each realisation
DetectedResults = InvasionResults

###PoA observation contract. Routine pre-deployed vertebrate control is a
###Background observation opportunity because it operates regardless of HaveInfo.
###Ordinary Background and InfoTriggered surveillance remain separate draws.
BackgroundDetectedResults = array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
BackgroundSurveillanceDetectedResults = array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
InfoTriggeredDetectedResults = array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
InformationStateBeforeSurveillanceResults = array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
HaveInfoResults = array(0L, dim=c(nrow(SDDprob),Ntimesteps,Nperm))
BackgroundDetectionProbabilityResults = array(0, dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))
InfoTriggeredDetectionProbabilityResults = array(0, dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))
###Pre-control state/probability are retained so the PoA companion can calculate
###the no-detection likelihood outside the biological engine.
RoutineControlObservationAbundanceResults = array(0, dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))
RoutineControlDetectionProbabilityResults = array(0, dim=c(nrow(SDDprob),Nstages,Ntimesteps,Nperm))

###Declare array for tracking management adoption status 
###of individual nodes in each timestep of each realisation
###This is a measure of potential disruption to farm businesses
###or ongoing management burden (surveillance and removal)
###for publicly-owned lands
ManagingResults = InvasionResults

###Vertebrate specialist control diagnostics
ControlDeathResults = array(0, dim = c(nrow(SDDprob), Nstages, Ntimesteps, Nperm))
ControlDetectionResults = array(0, dim = c(nrow(SDDprob), Ntimesteps, Nperm))
ControlCostResults = matrix(0, nrow = Ntimesteps, ncol = Nperm)

###Declare matrix for information spread simulations
if(is.matrix(SEAM) == T)
     {
     diag(SEAM) = 0
     RandSEAM <- matrix(NA,nrow = nrow(SDDprob),ncol=nrow(SDDprob))
     }

###Validate management-induced fecundity reduction
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

###Assign standard deviation value to management in no value provided
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
###Start of simulation
###########################################################
    
for (perm in 1:Nperm) 
{ 
###Assign initial infestations according either to "InitialInvasion" binary vector OR
###"InvasionRisk" probabilities and/or initial proportion of nodes infested ("InitBioP") OR
###just "InitBioP" if neither "InitialInvasion" or "InvasionRisk" supplied by user
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

# initialise the population
N <- InitBio
if(sum(N) == 0 && OngoingExternalInvasion == F)
  warning("No initial populations and no future external invasions")

###Select nodes with information at start of simulation  according either to "InitialInfo" binary vector OR
###"ExternalInfoProb" probabilities and/or initial proportion of nodes with information ("InitInfoP") OR
###just "InitInfoP" if neither "InitialInfo" or "ExternalInfoProb" supplied by user.
###If no initial info variables provided, no nodes have info at start of simulations
# Robust InitInfo setup
n_nodes <- nrow(SDDprob)  # number of nodes
InitInfo <- rep(0, n_nodes)  # default: no nodes initially informed

# Only proceed if any initial info or probabilities are provided
if(!all(is.na(InitialInfo)) || !is.na(InitInfoP) || !all(is.na(ExternalInfoProb))) {
  
  # If InitialInfo is valid length, just use it
  if(length(InitialInfo) == n_nodes) {
    InitInfo <- InitialInfo
    
  } else {
    # Ensure InitInfoP is numeric and in [0,1]
    if(is.na(InitInfoP)) InitInfoP <- 0
    
    # Ensure ExternalInfoProb is numeric of length n_nodes
    if(is.na(sum(ExternalInfoProb))) ExternalInfoProb <- rep(1, n_nodes)
    if(length(ExternalInfoProb) != n_nodes) ExternalInfoProb <- rep(1, n_nodes)
    
    # Determine number of nodes to sample
    n_sample <- ceiling(n_nodes * InitInfoP)
    
    # Only sample if n_sample > 0
    if(n_sample > 0) {
      Info <- sample(1:n_nodes, size = n_sample, prob = ExternalInfoProb)
      InitInfo[Info] <- 1
    }
  }
}





# --- Determine dimensions ---
n_nodes <- nrow(SDDprob)

###Randomly assign detection probabilities
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

###Realised InfoTriggered detection probabilities for this permutation.
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


###Randomly assign probability of mangement adoption upon detection of infestation
###If ManageProb given as single value or vector (nodes)
if(is.matrix(ManageProb)==FALSE &&(length(ManageProb) == 1 ||length(ManageProb) == nrow(SDDprob) ))
      {
      NodeManageProb = rnorm(ManageProb,ManageSD,n = nrow(SDDprob))
      NodeManageProb[NodeManageProb<0] = 0
      NodeManageProb[NodeManageProb>1] = 1
      }

###Randomly assign spread reduction factor when management adopted
###If SpreadReduction given as single value or vector (nodes)
if(is.matrix(SpreadReduction)==FALSE &&(length(SpreadReduction) == 1 ||length(SpreadReduction) == nrow(SDDprob) ))
      {
      NodeSpreadReduction = rnorm(SpreadReduction,SpreadReductionSD,n = nrow(SDDprob))
      NodeSpreadReduction[NodeSpreadReduction<0] = 0
      NodeSpreadReduction[NodeSpreadReduction>1] = 1
      }

###Randomly assign mortality probability when management applied
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






###Populate invasion status vector ahead of timestep loop
Invaded <- ifelse(rowSums(InitBio) > 0, 1, 0)

###Probability of info at start of simulation depends on
###Presence of pest and detection probability
###Select nodes that have detected infestation 




# Calculate probability of detection per node
prob_detect <- 1 - apply((1 - NodeDetectionProb[,,1])^InitBio, 1, prod)

# Draw initial detection (0/1) per node
InitDetection <- rbinom(n = nrow(InitBio), size = 1, prob = prob_detect)


###Add detections to nodes which already have info (e.g. pre-emptive control and hygiene measures)
InitInfo[InitInfo == 0] = InitDetection[InitInfo == 0]

###Populate information status vector ahead of timestep loop
HaveInfo = InitInfo

###Track the most recent timestep with known local presence
LastKnownPresence = rep(NA,nrow(SDDprob))
if(UseInfoPersistence == T)
  {
  InitialKnownPresence = which(InitDetection == 1)
  if(length(InitialKnownPresence) > 0)
    LastKnownPresence[InitialKnownPresence] = 0
  }


  # run simulation
  for (timestep in 1:Ntimesteps) 
    { 
    ###Print progress
    if(DoProgress) cat("\r", "Realisation ", perm, "Timestep ", timestep, "...")
 
    ###Allow for variation in dispersal connectivity through time
    NodeSDDprob = SDDprob
    if(length(dim(SDDprob)) == 3)
      NodeSDDprob = SDDprob[,,timestep]
    NodeLDDprob = LDDprob
    if(length(dim(LDDprob)) == 3)
      NodeLDDprob = LDDprob[,,timestep]
    NodeTransitionSDDprob <- ResolveTransitionMovement(TransitionSDDprob, timestep)
    NodeTransitionLDDprob <- ResolveTransitionMovement(TransitionLDDprob, timestep)
    NodeFecundityReduction <- ResolveFecundityReductionTM(timestep)
###Allow for variation in recruit establishment through time
    ###e.g.  climate change predictions
    ###Note: could be done outside loop, but would take heaps of memory to store 
    if(is.matrix(EnvEstabProb) == T)
      NodeEnvEstabProb <- EnvEstabProb[,timestep]
    
    if (is.array(Transition) && length(dim(Transition)) == 3) {
      NodeTransition <- Transition[,,timestep]
      } else if (is.list(Transition)) {
        NodeTransition <- Transition
      } else {
      NodeTransition <- Transition
      }
    
    
    ###If carrying capacity provided as matrix assign values for relevant timestep
    if(is.matrix(K) == TRUE)
      {
      K_is_0 <- K[,timestep]<=0
      inv_K <- 1 / sum(K[,timestep])
      NodeK = K[,timestep] 
      }  
    
    ###If seedbank carrying capacity provided as matrix assign values for relevant timestep
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
  ###Retain the exact pre-control population exposed to routine devices/area control.
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

    ###Whole-node detection probability implied by per-individual stage exposure.
    ControlProbDetectNode <- 1 - apply((1 - ce$detect_prob)^N, 1, prod)
    ControlDetectedNodes <- rbinom(nrow(N), size = 1, prob = ControlProbDetectNode) == 1
    ControlDetectionResults[,timestep,perm] <- as.integer(ControlDetectedNodes)

    NpreControl <- matrix(0, nrow = nrow(N), ncol = Nstages)
    for(s in seq_len(Nstages))
      NpreControl[,s] <- rbinom(nrow(N), size = N[,s], prob = 1 - ce$kill_prob[,s])
    ControlDeathResults[,,timestep,perm] <- N - NpreControl
    N <- NpreControl
    }

  ###Randomly assign probability of mangement adoption upon detection of infestation
  ###If ManageProb given as matrix (nodes x timesteps)
  if(is.matrix(ManageProb)==TRUE && nrow(ManageProb) == nrow(SDDprob) && ncol(ManageProb) == Ntimesteps)
   	{	
   	NodeManageProb = rnorm(ManageProb[,timestep],ManageSD,n = nrow(SDDprob))
   	NodeManageProb[NodeManageProb<0] = 0
   	NodeManageProb[NodeManageProb>1] = 1
   	}

  ###Randomly assign spread reduction factor when management adopted
  ###If SpreadReduction given as matrix (nodes x timesteps)
  if(is.matrix(SpreadReduction)==TRUE && nrow(SpreadReduction) == nrow(SDDprob) && ncol(SpreadReduction) == Ntimesteps)
   	{	
   	NodeSpreadReduction = rnorm(SpreadReduction[,timestep],ManageSD,n = nrow(SDDprob))
   	NodeSpreadReduction[NodeSpreadReduction<0] = 0
   	NodeSpreadReduction[NodeSpreadReduction>1] = 1
   	}
  
  ###Assign management status to nodes   
  ###Management is only applied to nodes which have information
  ###i.e. where pest has been detected or following communication of information
  ###from neighbouring infested farms 
  Managing = rbinom(1:nrow(SDDprob),size = 1,prob = NodeManageProb*HaveInfo)
  
  ###Management is only applied to nodes which have information
  ###i.e. where pest has been detected or following communication of information
  ###from neighbouring infested farms 
  Managing = Managing*HaveInfo
  
  ###Identify nodes with known extant infestations 
  Detected = Invaded*HaveInfo
  
  ###Apply management
  # N: current population matrix (nodes x stages)
  # Managing: vector of management adoption per node (0/1) for current timestep
  # NodeMortalityProb: array (nodes x stages x timesteps)
  # timestep: current timestep index
  
  # Apply management-driven mortality
  N0 <- matrix(NA, nrow = nrow(N), ncol = ncol(N))
  for (s in 1:Nstages) {
    N0[, s] <- rbinom(n = nrow(N), size = N[, s], prob = 1 - (NodeMortalityProb[, s,timestep] * Managing))
  }
  
  ###Track known local presence from actual management mortality
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
  RandSEAM[] <- rbinom(n=nrow(SDDprob)^2, size=1, prob = SEAM*Detected)
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
  if(is.na(IncursionStartPop) == T) 
	  N[,1] = N[,1]+ExternalInvasion
  if(is.na(IncursionStartPop) == F) 
	  N[,1] = N[,1]+ExternalInvasion*IncursionStartPop
 
  }

  ###Optional aggregate social/contact/disease-state interaction. The hook can
  ###redistribute the existing node x class population but cannot change array
  ###dimensions; births and deaths belong in their dedicated processes.
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
 Invaded = ifelse(rowSums(N)>0,1,0)
 
 ###Record nodes adopting management
 ManagingResults[,timestep,perm] = Managing
  
 ###Record infested nodes
 InvasionResults[,timestep,perm] = Invaded

 ###Record populations
weighted_population <- if (is.matrix(Weights)) {
  rowSums(N[, 2:Nstages, drop = FALSE] *
            Weights[, 2:Nstages, drop = FALSE])
} else {
  as.numeric(N[, 2:Nstages, drop = FALSE] %*% Weights[2:Nstages])
}
PopulationResults[, timestep, perm] <- weighted_population
 
 ###Record stage populations
 PopulationStageResults[,,timestep,perm] = N

 ###Freeze information available before this round's surveillance.
 ###Routine control detections occurred earlier but are deliberately not registered
 ###until after this snapshot, so they cannot activate InfoTriggered surveillance
 ###retrospectively in the same timestep.
 InfoBeforeSurveillance <- as.integer(HaveInfo != 0)
 InformationStateBeforeSurveillanceResults[,timestep,perm] <- InfoBeforeSurveillance

 ###Record the realised per-individual observation probabilities used this round.
 BackgroundDetectionProbabilityResults[,,timestep,perm] <- NodeDetectionProb[,,timestep]
 InfoTriggeredDetectionProbabilityResults[,,timestep,perm] <- NodeInfoTriggeredDetectionProb[,,timestep]

 ###Ordinary Background surveillance after population dynamics.
 DetectionProbPerStage <- 1 - (1 - NodeDetectionProb[,,timestep])^N
 ProbDetectNode <- 1 - apply(1 - DetectionProbPerStage, 1, prod)
 BackgroundDetection <- rbinom(n=nrow(SDDprob), size=1, prob=ProbDetectNode)

 ###Additional targeted surveillance is gated only by information that existed
 ###before the current observation results.
 InfoTriggeredDetection <- integer(nrow(SDDprob))
 if(UseInfoTriggeredSurveillance) {
   InfoDetectionProbPerStage <- 1 - (1 - NodeInfoTriggeredDetectionProb[,,timestep])^N
   ProbInfoDetectNode <- 1 - apply(1 - InfoDetectionProbPerStage, 1, prod)
   InfoTriggeredDetection <- rbinom(
     n=nrow(SDDprob), size=1, prob=ProbInfoDetectNode*InfoBeforeSurveillance)
 }

 ###Routine control and ordinary survey detections are both Background evidence.
 ###The combined field stores event counts, so a node detected by both pathways
 ###can contribute two observation events in one timestep.
 BackgroundSurveillanceDetectedResults[,timestep,perm] <- BackgroundDetection
 BackgroundDetectedResults[,timestep,perm] <- as.integer(ControlDetectedNodes) + BackgroundDetection
 InfoTriggeredDetectedResults[,timestep,perm] <- InfoTriggeredDetection
 HostDetectionEvidence <- pmax(as.integer(ControlDetectedNodes), BackgroundDetection, InfoTriggeredDetection)

 ###All observation pathways update response information only after surveillance.
 if(UseInfoPersistence == T) {
   KnownPresence <- which(HostDetectionEvidence == 1)
   if(length(KnownPresence) > 0) LastKnownPresence[KnownPresence] <- timestep
 }
 HaveInfo[HaveInfo == 0] <- HostDetectionEvidence[HaveInfo == 0]
 HaveInfoResults[,timestep,perm] <- HaveInfo

 ###Legacy DetectedResults remains the persistent known-present state.
 DetectedResults[,timestep,perm] = HaveInfo*Invaded 
 }
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
###Store annual node-level invasion probs for heat maps
###and estimation of invasion threat to other regions
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
###Produce summary figs when processing completed
###########################################################

Title = ModelName


###Change in total population with time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

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

###Change in number of nodes infested with time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

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


###Change in number of farms managing through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided

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


###Change in number of known extant infestations through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided 

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


###Change in proportion of extant infestations detected through time
###Plots of raw values for each realisation and summaries (median and 95% CI) provided
 
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

###Return a structured object as well as retaining the parent function's RDS
###outputs. This is additive and does not alter the legacy saved results.
out <- list(
  ModelName = ModelName,
  ###Legacy vertebrate result aliases
  Population = PopulationResults,
  PopulationStage = PopulationStageResults,
  Invasion = InvasionResults,
  Detection = DetectedResults,
  Managing = ManagingResults,
  ControlDeaths = ControlDeathResults,
  ControlDetections = ControlDetectionResults,
  ControlCost = ControlCostResults,
  InvasionProbability = InvasionProb,
  ###Canonical INApest / PoA result contract
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
###End of function
################################################################
################################################################
