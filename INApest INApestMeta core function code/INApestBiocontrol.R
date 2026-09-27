###############################################################################
### INApestBiocontrol -- non-pathogen biocontrol companion for INApest
### Definitive commented source, validated 2026-09-18.
###
### The companion represents one or more independent biocontrol-agent
### populations alongside an INApest pest/host population.  It is deliberately
### kept separate from the pest state so the same agent model can be coupled to
### several non-point INApest architectures without duplicating biology inside
### each parent simulation engine.
###
### Supported parent architectures in this validated release
### ---------------------------------------------------------
###   * Binary node occupancy
###   * Meta abundance by node
###   * Multiple Land Use abundance by node x land-use class
###   * Transition Matrix abundance by node x host stage
###   * Vertebrate Node abundance by node x host stage
###
### Agent state is node based.  Each agent can have its own life stages,
### releases, stage transition matrix, spatial movement matrix, attacking stage,
### host target stage(s), attack rate and offspring produced per successful
### attack.  Several agents can act in the same timestep.  Competing agents are
### evaluated as simultaneous hazards so results do not depend on arbitrary
### agent ordering in the input list.
###
### The companion uses the same Transition Matrix orientation as INApest:
### columns are source stages and rows are destination stages.  Residual column
### probability is mortality.  Agent movement matrices use source rows and
### destination columns; residual row probability is export from the modelled
### landscape.
###
### Validated event order within one biocontrol step
### ------------------------------------------------
###   1. Add scheduled releases.
###   2. Calculate simultaneous attack hazards and remove attacked hosts.
###   3. Advance existing biocontrol stages through their transition matrices.
###   4. Add attack-derived recruits to the configured recruit stage.
###   5. Move configured agent stages among nodes.
###
### This source is a comment-only documentation pass over the validated v0.2
### executable code.  No executable statement was changed by this source pass;
### COMMENT_ONLY_VERIFICATION.txt in the definitive bundle records that check.
###############################################################################

# Keep probabilities inside [0, 1] after numerical calculations.
# This mirrors the defensive probability clipping used throughout INApest.
.inabc_clip01 <- function(x) pmin(1, pmax(0, x))


# Evaluate a time-varying user resolver using only arguments it declares.
# Non-function inputs pass through unchanged.  This lets schedules be supplied
# either directly or as functions without forcing one rigid function signature.
.inabc_call_resolver <- function(x, timestep, context, agent = NULL) {
  if (!is.function(x)) return(x)
  fm <- names(formals(x))
  a <- list(timestep = timestep, context = context, agent = agent)
  if (!is.null(fm) && !("..." %in% fm)) a <- a[intersect(names(a), fm)]
  do.call(x, a)
}


# Resolve a scalar/node/timestep attack-rate style input to one value per node
# for the current timestep.  Ambiguous shapes are rejected rather than guessed.
.inabc_resolve_node <- function(x, timestep, context, name, agent = NULL) {
  x <- .inabc_call_resolver(x, timestep, context, agent)
  n <- context$n_nodes
  nt <- context$Ntimesteps
  d <- dim(x)
  if (is.null(d)) {
    x <- as.numeric(x)
    if (length(x) == 1L) return(rep(x, n))
    if (length(x) == n) return(x)
    if (length(x) == nt && nt != n) return(rep(x[timestep], n))
    stop(name, " must be scalar, length n_nodes, nodes x Ntimesteps, or a resolver function")
  }
  if (length(d) == 2L && identical(d, c(n, nt))) return(as.numeric(x[, timestep]))
  stop(name, " must be scalar, length n_nodes, nodes x Ntimesteps, or a resolver function")
}


# Resolve the physical duration represented by the current INApest timestep.
# Attack rates are continuous-time hazards, so this value controls the mapping
# p(attack) = 1 - exp(-hazard * dt).
.inabc_resolve_dt <- function(x, timestep, context) {
  z <- .inabc_call_resolver(x, timestep, context, NULL)
  if (length(z) == 1L) out <- as.numeric(z)
  else if (length(z) == context$Ntimesteps) out <- as.numeric(z[timestep])
  else stop("TimestepLength must be scalar, length Ntimesteps, or a resolver function")
  if (!is.finite(out) || out <= 0) stop("TimestepLength must resolve to a finite value > 0")
  out
}


# Resolve an agent life-stage transition matrix for the current timestep.
# INApest convention is retained: matrix columns are source stages and rows are
# destination stages; unused column mass represents mortality.
.inabc_resolve_transition <- function(x, timestep, context, agent) {
  x <- .inabc_call_resolver(x, timestep, context, agent)
  s <- length(agent$Stages)
  d <- dim(x)
  if (length(d) == 2L && identical(d, c(s, s))) out <- x
  else if (length(d) == 3L && identical(d, c(s, s, context$Ntimesteps))) out <- x[, , timestep]
  else stop("Agent Transition must be stages x stages, stages x stages x Ntimesteps, or a resolver function")
  out <- as.matrix(out)
  if (any(!is.finite(out)) || any(out < 0)) stop("Agent Transition entries must be finite and non-negative")
  # Match the parent INApest Transition Matrix convention: columns are source
  # stages and rows are destination stages. Residual column mass is mortality.
  if (any(colSums(out) > 1 + 1e-12)) stop("Agent Transition column sums must be <= 1; residual mass is mortality")
  out
}


# Resolve a source-node x destination-node movement matrix.  Missing row mass
# is treated as movement outside the modelled landscape (export).
.inabc_resolve_movement <- function(x, timestep, context, agent) {
  if (is.null(x)) return(NULL)
  x <- .inabc_call_resolver(x, timestep, context, agent)
  n <- context$n_nodes
  d <- dim(x)
  if (length(d) == 2L && identical(d, c(n, n))) out <- x
  else if (length(d) == 3L && identical(d, c(n, n, context$Ntimesteps))) out <- x[, , timestep]
  else stop("Agent Movement must be nodes x nodes, nodes x nodes x Ntimesteps, or a resolver function")
  out <- as.matrix(out)
  if (any(!is.finite(out)) || any(out < 0)) stop("Agent Movement entries must be finite and non-negative")
  if (any(rowSums(out) > 1 + 1e-12)) stop("Agent Movement row sums must be <= 1; residual mass is export")
  out
}


# Resolve starting state or scheduled releases to a node x agent-stage matrix.
# Scalar/node-vector inputs are placed in ReleaseStage; explicit matrices can
# seed or release several agent stages simultaneously.
.inabc_resolve_state_matrix <- function(x, timestep, context, agent, name) {
  x <- .inabc_call_resolver(x, timestep, context, agent)
  n <- context$n_nodes
  s <- length(agent$Stages)
  d <- dim(x)
  if (is.null(d)) {
    x <- as.numeric(x)
    if (length(x) == 1L) {
      out <- matrix(0, n, s, dimnames = list(NULL, agent$Stages))
      out[, agent$ReleaseStage] <- x
      return(out)
    }
    if (length(x) == n) {
      out <- matrix(0, n, s, dimnames = list(NULL, agent$Stages))
      out[, agent$ReleaseStage] <- x
      return(out)
    }
    stop(name, " must be scalar, length n_nodes, nodes x stages, nodes x stages x Ntimesteps, or a resolver function")
  }
  if (length(d) == 2L && identical(d, c(n, s))) return(matrix(as.numeric(x), n, s, dimnames = list(NULL, agent$Stages)))
  if (length(d) == 3L && identical(d, c(n, s, context$Ntimesteps))) {
    return(matrix(as.numeric(x[, , timestep]), n, s, dimnames = list(NULL, agent$Stages)))
  }
  stop(name, " must be scalar, length n_nodes, nodes x stages, nodes x stages x Ntimesteps, or a resolver function")
}


# Allocate an integer count among competing destinations/effects in one draw.
# The helper is used for stage transitions, movement and competition among
# biocontrol agents so total integer counts are conserved exactly.
.inabc_allocate_multinomial <- function(n, prob) {
  n <- as.integer(max(0, floor(n)))
  if (n == 0L) return(integer(length(prob)))
  p <- pmax(0, as.numeric(prob))
  if (sum(p) <= 0) return(integer(length(prob)))
  as.integer(stats::rmultinom(1L, n, p)[, 1L])
}


# Advance one agent population through its life-stage transition matrix.
# Every source-stage individual is allocated once among destination stages plus
# an implicit mortality category.
.inabc_transition_agent <- function(state, P) {
  n <- nrow(state); s <- ncol(state)
  out <- matrix(0L, n, s, dimnames = dimnames(state))
  for (i in seq_len(n)) {
    for (from in seq_len(s)) {
      count <- as.integer(max(0, floor(state[i, from])))
      if (!count) next
      # INApest stage-transition convention: column = source, row = destination.
      probs <- c(P[, from], max(0, 1 - sum(P[, from])))
      z <- .inabc_allocate_multinomial(count, probs)
      out[i, ] <- out[i, ] + z[seq_len(s)]
    }
  }
  out
}


# Move selected life stages among nodes after local demography and recruitment.
# Individuals assigned to residual row probability leave the modelled system.
.inabc_move_agent <- function(state, P, stages) {
  if (is.null(P)) return(state)
  n <- nrow(state); s <- ncol(state)
  out <- state
  for (st in stages) {
    moved <- integer(n)
    exported <- integer(n)
    for (i in seq_len(n)) {
      count <- as.integer(max(0, floor(state[i, st])))
      if (!count) next
      probs <- c(P[i, ], max(0, 1 - sum(P[i, ])))
      z <- .inabc_allocate_multinomial(count, probs)
      moved <- moved + z[seq_len(n)]
      exported[i] <- z[n + 1L]
    }
    out[, st] <- moved
  }
  out
}


###############################################################################
### Define one biocontrol agent
###############################################################################
# This constructor stores the biological rules for one independent agent.
# Values may be constant or, where resolved by the engine, time varying.
INApestBiocontrolAgent <- function(
    Name,                         # Unique agent name used in histories and outputs
    Stages = c("juvenile", "adult"), # Ordered agent life-stage names
    InitialState = 0,              # Starting agent abundance: scalar/node vector/node x stage
    Release = 0,                   # Scheduled release abundance in each timestep
    ReleaseStage = length(Stages),  # Stage receiving scalar/vector releases
    Transition = diag(length(Stages)), # Agent transition matrix; columns source, rows destination
    Movement = NULL,               # Optional source-node x destination-node movement matrix
    MovementStages = seq_along(Stages), # Agent stages to which Movement is applied
    AttackStage = length(Stages),   # Agent stage responsible for host attack
    TargetStage = NULL,             # Host stage(s) attacked in stage-structured parent models
    AttackRate = 0,                # Continuous-time per-agent attack-rate parameter by node
    RecruitStage = 1L,              # Agent stage receiving offspring from successful attacks
    OffspringPerAttack = 1L) {      # Integer offspring/recruits produced per successful attack

  Stages <- as.character(Stages)
  if (!length(Stages) || any(!nzchar(Stages)) || anyDuplicated(Stages)) stop("Stages must be unique non-empty names")
  s <- length(Stages)
  idx <- function(z, nm) {
    if (is.character(z)) z <- match(z, Stages)
    z <- as.integer(z)
    if (length(z) != 1L || is.na(z) || z < 1L || z > s) stop(nm, " must identify one agent stage")
    z
  }
  ReleaseStage <- idx(ReleaseStage, "ReleaseStage")
  AttackStage <- idx(AttackStage, "AttackStage")
  RecruitStage <- idx(RecruitStage, "RecruitStage")
  if (is.character(MovementStages)) MovementStages <- match(MovementStages, Stages)
  MovementStages <- as.integer(MovementStages)
  if (anyNA(MovementStages) || any(MovementStages < 1L | MovementStages > s)) stop("MovementStages contains an invalid stage")
  OffspringPerAttack <- as.integer(OffspringPerAttack)
  if (length(OffspringPerAttack) != 1L || is.na(OffspringPerAttack) || OffspringPerAttack < 0L)
    stop("OffspringPerAttack must be a non-negative whole number")

  structure(list(
    Name = as.character(Name)[1L], Stages = Stages,
    InitialState = InitialState, Release = Release, ReleaseStage = ReleaseStage,
    Transition = Transition, Movement = Movement, MovementStages = MovementStages,
    AttackStage = AttackStage, TargetStage = TargetStage, AttackRate = AttackRate,
    RecruitStage = RecruitStage, OffspringPerAttack = OffspringPerAttack
  ), class = "INApestBiocontrolAgent")
}


###############################################################################
### Describe the parent INApest host architecture
###############################################################################
# Parent engines create this small context object once, allowing the companion
# to validate target dimensions without depending on parent-engine internals.
INApestBiocontrolContext <- function(
    Architecture = c("binary", "meta", "mlu", "transition", "vertebrate_node"), # Parent host-state architecture
    n_nodes,                       # Number of spatial nodes
    Ntimesteps,                    # Number of simulation timesteps
    n_landuses = NULL,             # Number of land-use classes for MLU only
    host_stages = NULL) {          # Host stage labels/indices for stage-structured models
  Architecture <- match.arg(Architecture)
  n_nodes <- as.integer(n_nodes)
  Ntimesteps <- as.integer(Ntimesteps)
  if (length(n_nodes) != 1L || is.na(n_nodes) || n_nodes < 1L) stop("n_nodes must be a positive integer")
  if (length(Ntimesteps) != 1L || is.na(Ntimesteps) || Ntimesteps < 1L) stop("Ntimesteps must be a positive integer")
  if (Architecture == "mlu") {
    n_landuses <- as.integer(n_landuses)
    if (length(n_landuses) != 1L || is.na(n_landuses) || n_landuses < 1L) stop("MLU context requires positive n_landuses")
  }
  if (Architecture %in% c("transition", "vertebrate_node")) {
    if (is.null(host_stages)) stop("Stage-structured context requires host_stages")
    if (length(host_stages) == 1L && is.numeric(host_stages)) host_stages <- seq_len(as.integer(host_stages))
  }
  structure(list(
    Architecture = Architecture, n_nodes = n_nodes, Ntimesteps = Ntimesteps,
    n_landuses = n_landuses, host_stages = host_stages
  ), class = "INApestBiocontrolContext")
}


###############################################################################
### Assemble one or more agents into an INApest biocontrol companion
###############################################################################
# The returned object exposes Validate, Initial and Step methods through Engine.
# Parent INApest sources call only this stable companion contract.
INApestBiocontrol <- function(Agents, TimestepLength = 1) { # Agents plus physical duration of each simulation timestep
  if (inherits(Agents, "INApestBiocontrolAgent")) Agents <- list(Agents)
  if (!is.list(Agents) || !length(Agents) || !all(vapply(Agents, inherits, logical(1), "INApestBiocontrolAgent")))
    stop("Agents must be one INApestBiocontrolAgent or a list of them")
  nm <- vapply(Agents, `[[`, character(1), "Name")
  if (any(!nzchar(nm)) || anyDuplicated(nm)) stop("Biocontrol agent names must be unique and non-empty")
  names(Agents) <- nm

  # Check all agent schedules against the parent model before simulation starts.
  # Up-front validation avoids discovering malformed time-varying inputs part-way
  # through a long stochastic or parallel run.
  validate <- function(context) {
    if (!inherits(context, "INApestBiocontrolContext")) stop("context must be INApestBiocontrolContext()")
    for (a in Agents) {
      if (context$Architecture %in% c("transition", "vertebrate_node")) {
        ts <- a$TargetStage
        if (is.character(ts)) ts <- match(ts, context$host_stages)
        ts <- as.integer(ts)
        if (!length(ts) || anyNA(ts) || any(ts < 1L | ts > length(context$host_stages)))
          stop("Every agent in a stage-structured host model must specify valid TargetStage value(s)")
      } else if (!is.null(a$TargetStage)) {
        stop("TargetStage is only used with stage-structured host architectures")
      }

      # Validate every scheduled value up front so a malformed time-varying
      # specification fails before a long stochastic run begins.
      for (tt in seq_len(context$Ntimesteps)) {
        .inabc_resolve_dt(TimestepLength, tt, context)
        .inabc_resolve_transition(a$Transition, tt, context, a)
        .inabc_resolve_movement(a$Movement, tt, context, a)
        ar <- .inabc_resolve_node(a$AttackRate, tt, context, "AttackRate", a)
        if (any(!is.finite(ar)) || any(ar < 0))
          stop("AttackRate must resolve to finite non-negative values")
        rel <- .inabc_resolve_state_matrix(a$Release, tt, context, a, "Release")
        if (any(!is.finite(rel)) || any(rel < 0))
          stop("Release must resolve to finite non-negative counts")
      }
    }
    invisible(TRUE)
  }


  # Build the initial node x agent-stage state for each configured agent.
  # Target is accepted for a stable parent-engine contract even where initial
  # agent abundance does not depend on current host abundance.
  initial <- function(Target, context) {
    validate(context)
    out <- lapply(Agents, function(a) {
      z <- .inabc_resolve_state_matrix(a$InitialState, 1L, context, a, "InitialState")
      if (any(!is.finite(z)) || any(z < 0)) stop("InitialState must contain finite non-negative counts")
      matrix(as.integer(round(z)), nrow = context$n_nodes, ncol = length(a$Stages), dimnames = list(NULL, a$Stages))
    })
    names(out) <- names(Agents)
    out
  }


  # Apply one complete biocontrol timestep to the supplied host Target.
  # The function returns the updated host target, updated independent agent
  # populations, and explicit attack/recruitment impacts for diagnostics.
  step <- function(Target, State, timestep, context) {
    if (!inherits(context, "INApestBiocontrolContext")) stop("context must be INApestBiocontrolContext()")
    if (!is.list(State) || !identical(names(State), names(Agents))) stop("Biocontrol State does not match configured agents")
    dt <- .inabc_resolve_dt(TimestepLength, timestep, context)

    # Add scheduled releases before attack so released adults can act immediately.
    current <- vector("list", length(Agents)); names(current) <- names(Agents)
    for (j in seq_along(Agents)) {
      a <- Agents[[j]]
      z <- State[[j]]
      rel <- .inabc_resolve_state_matrix(a$Release, timestep, context, a, "Release")
      if (any(!is.finite(rel)) || any(rel < 0)) stop("Release must contain finite non-negative counts")
      current[[j]] <- z + matrix(as.integer(round(rel)), nrow = nrow(z), ncol = ncol(z), dimnames = dimnames(z))
    }

    # Present the parent host state as comparable node-level components.
    # Binary/Meta have one component; MLU has one per land-use class; TM and
    # Vertebrate Node have one per host stage.  This common representation keeps
    # attack biology independent of the parent host storage format.
    arch <- context$Architecture
    TargetOut <- Target
    if (arch %in% c("binary", "meta")) {
      if (length(Target) != context$n_nodes) stop("Binary/Meta target must have length n_nodes")
      comps <- list(as.numeric(Target))
      comp_keys <- "all"
    } else if (arch == "mlu") {
      if (!is.matrix(Target) || !identical(dim(Target), c(context$n_nodes, context$n_landuses)))
        stop("MLU target must be n_nodes x n_landuses")
      comps <- lapply(seq_len(context$n_landuses), function(k) as.numeric(Target[, k]))
      comp_keys <- paste0("landuse_", seq_len(context$n_landuses))
    } else {
      # Transition and Vertebrate Node hosts share the node x host-stage contract.
      s_host <- length(context$host_stages)
      if (!is.matrix(Target) || !identical(dim(Target), c(context$n_nodes, s_host)))
        stop("Transition target must be n_nodes x host stages")
      comps <- lapply(seq_len(s_host), function(k) as.numeric(Target[, k]))
      comp_keys <- as.character(context$host_stages)
    }

    attacks_by_agent <- lapply(Agents, function(a) matrix(0L, nrow = context$n_nodes, ncol = length(comps), dimnames = list(NULL, comp_keys)))
    names(attacks_by_agent) <- names(Agents)
    recruits <- setNames(integer(length(Agents)), names(Agents))
    recruits_node <- lapply(Agents, function(a) integer(context$n_nodes)); names(recruits_node) <- names(Agents)

    # Apply attacks simultaneously within each host component.  For agent j,
    # hazard_j = AttackRate_j * abundance of its attacking stage.  Total host
    # attack probability is 1-exp(-sum(hazard_j)*dt).  Realised host losses are
    # drawn once and then allocated among agents in proportion to their hazards,
    # avoiding arbitrary sequential competition among agents.
    for (cc in seq_along(comps)) {
      H <- pmax(0L, as.integer(floor(comps[[cc]])))
      eligible <- vapply(Agents, function(a) {
        if (!(arch %in% c("transition", "vertebrate_node"))) return(TRUE)
        ts <- a$TargetStage
        if (is.character(ts)) ts <- match(ts, context$host_stages)
        cc %in% as.integer(ts)
      }, logical(1))
      jj <- which(eligible)
      if (!length(jj)) next

      hazards <- matrix(0, nrow = context$n_nodes, ncol = length(jj))
      for (kk in seq_along(jj)) {
        a <- Agents[[jj[kk]]]
        rate <- .inabc_resolve_node(a$AttackRate, timestep, context, "AttackRate", a)
        active <- current[[jj[kk]]][, a$AttackStage]
        hazards[, kk] <- pmax(0, rate * active)
      }
      total_h <- rowSums(hazards)
      p_attack <- 1 - exp(-total_h * dt)
      p_attack <- .inabc_clip01(p_attack)
      total_kill <- stats::rbinom(context$n_nodes, size = H, prob = p_attack)

      for (i in which(total_kill > 0L)) {
        w <- hazards[i, ]
        if (sum(w) <= 0) next
        alloc <- .inabc_allocate_multinomial(total_kill[i], w)
        for (kk in seq_along(jj)) {
          j <- jj[kk]
          attacks_by_agent[[j]][i, cc] <- alloc[kk]
          recruits_node[[j]][i] <- recruits_node[[j]][i] + alloc[kk] * Agents[[j]]$OffspringPerAttack
        }
      }
      comps[[cc]] <- H - total_kill
    }

    # Reconstruct exactly the host-state shape expected by the parent engine.
    # Binary architecture is reduced back to occupancy after host attack.
    if (arch %in% c("binary", "meta")) TargetOut <- comps[[1L]]
    else if (arch == "mlu") {
      TargetOut <- do.call(cbind, comps)
      dimnames(TargetOut) <- dimnames(Target)
    } else {
      TargetOut <- do.call(cbind, comps)
      dimnames(TargetOut) <- dimnames(Target)
    }
    if (arch == "binary") TargetOut <- as.integer(TargetOut > 0)

    # Advance the independent agent populations after host attack.  Existing
    # agents transition first; new attack-derived offspring then enter RecruitStage
    # and therefore cannot mature or die until the next timestep.  Configured
    # spatial movement is applied last.
    next_state <- vector("list", length(Agents)); names(next_state) <- names(Agents)
    for (j in seq_along(Agents)) {
      a <- Agents[[j]]
      P <- .inabc_resolve_transition(a$Transition, timestep, context, a)
      z <- .inabc_transition_agent(current[[j]], P)
      z[, a$RecruitStage] <- z[, a$RecruitStage] + recruits_node[[j]]
      M <- .inabc_resolve_movement(a$Movement, timestep, context, a)
      z <- .inabc_move_agent(z, M, a$MovementStages)
      next_state[[j]] <- z
      recruits[j] <- sum(recruits_node[[j]])
    }

    list(
      Target = TargetOut,
      State = next_state,
      Impact = list(AttacksByAgent = attacks_by_agent, RecruitsByAgent = recruits,
                    RecruitsByAgentNode = recruits_node)
    )
  }

  Engine <- list(Validate = validate, Initial = initial, Step = step)
  structure(list(Agents = Agents, TimestepLength = TimestepLength, Engine = Engine), class = "INApestBiocontrol")
}


###############################################################################
### Common result-history helpers used by parent INApest engines
###############################################################################


###############################################################################
### Allocate biocontrol result histories
###############################################################################
# Histories remain separate from legacy pest/host outputs so Biocontrol = NULL
# does not alter existing result structures or random-number use.
INApestBiocontrolHistory <- function(Biocontrol, Context, Nperm = 1L) {
  if (!inherits(Biocontrol, "INApestBiocontrol"))
    stop("Biocontrol must be created by INApestBiocontrol()")
  if (!inherits(Context, "INApestBiocontrolContext"))
    stop("Context must be created by INApestBiocontrolContext()")
  Nperm <- as.integer(Nperm)
  if (length(Nperm) != 1L || is.na(Nperm) || Nperm < 1L)
    stop("Nperm must be a positive integer")

  n_nodes <- Context$n_nodes
  nt <- Context$Ntimesteps
  n_components <- switch(
    Context$Architecture,
    binary = 1L,
    meta = 1L,
    mlu = Context$n_landuses,
    transition = length(Context$host_stages),
    vertebrate_node = length(Context$host_stages),
    stop("Unsupported biocontrol architecture")
  )
  component_names <- switch(
    Context$Architecture,
    binary = "all",
    meta = "all",
    mlu = paste0("landuse_", seq_len(Context$n_landuses)),
    transition = as.character(Context$host_stages),
    vertebrate_node = as.character(Context$host_stages)
  )

  state <- lapply(Biocontrol$Agents, function(a) {
    array(0L, dim = c(n_nodes, length(a$Stages), nt, Nperm),
          dimnames = list(NULL, a$Stages, NULL, NULL))
  })
  attacks <- lapply(Biocontrol$Agents, function(a) {
    array(0L, dim = c(n_nodes, n_components, nt, Nperm),
          dimnames = list(NULL, component_names, NULL, NULL))
  })
  recruits <- lapply(Biocontrol$Agents, function(a) {
    array(0L, dim = c(n_nodes, nt, Nperm))
  })
  names(state) <- names(attacks) <- names(recruits) <- names(Biocontrol$Agents)
  list(State = state, Attacks = attacks, Recruits = recruits)
}


# Record one timestep of agent state, attacks and attack-derived recruitment.
# The same helper is used by serial and worker-local parallel simulations.
INApestBiocontrolRecord <- function(History, State, Impact, timestep, perm = 1L) {
  if (is.null(History)) return(NULL)
  for (nm in names(State)) {
    History$State[[nm]][, , timestep, perm] <- State[[nm]]
    if (!is.null(Impact$AttacksByAgent[[nm]]))
      History$Attacks[[nm]][, , timestep, perm] <- Impact$AttacksByAgent[[nm]]
    if (!is.null(Impact$RecruitsByAgentNode[[nm]]))
      History$Recruits[[nm]][, timestep, perm] <- Impact$RecruitsByAgentNode[[nm]]
  }
  History
}


# Combine worker-local biocontrol histories across stochastic permutations.
# Semantic dimnames (agent stages and host components) are retained so combined
# PSOCK output can be indexed by the same labels as serial output.
INApestBiocontrolBindWorkerHistories <- function(worker_results, field = "BiocontrolHistory") {
  if (!length(worker_results)) return(NULL)
  first <- worker_results[[1L]][[field]]
  if (is.null(first)) return(NULL)
  nperm <- length(worker_results)
  out <- first
  # Worker histories are nodes x ... x timesteps with no permutation dimension,
  # or nodes x ... x timesteps x 1. Normalise to a final permutation dimension.
  bind_one <- function(parts) {
    first_part <- parts[[1L]]
    d <- dim(first_part)
    if (is.null(d)) stop("Biocontrol worker history must be an array")

    # Preserve semantic dimension labels (especially agent stages and host
    # components) when adding the worker/permutation dimension.  The previous
    # binder rebuilt arrays from dimensions alone, which silently dropped
    # dimnames and made otherwise-correct histories impossible to index by
    # stage name after serial-worker or PSOCK combination.
    dn <- dimnames(first_part)
    if (length(d) >= 1L && tail(d, 1L) == 1L) {
      dbase <- head(d, -1L)
      dnbase <- if (is.null(dn)) NULL else head(dn, -1L)
    } else {
      dbase <- d
      dnbase <- dn
    }
    ans_dn <- if (is.null(dnbase)) NULL else c(dnbase, list(NULL))
    ans <- array(0L, dim = c(dbase, nperm), dimnames = ans_dn)

    for (pp in seq_len(nperm)) {
      z <- parts[[pp]]
      dz <- dim(z)
      if (length(dz) == length(dbase) + 1L && tail(dz, 1L) == 1L)
        z <- array(z, dim = dbase)
      idx <- c(rep(list(TRUE), length(dbase)), list(pp))
      ans <- do.call(`[<-`, c(list(ans), idx, list(value = z)))
    }
    ans
  }
  for (section in c("State", "Attacks", "Recruits")) {
    for (nm in names(first[[section]])) {
      out[[section]][[nm]] <- bind_one(lapply(worker_results, function(x) x[[field]][[section]][[nm]]))
    }
  }
  out
}
