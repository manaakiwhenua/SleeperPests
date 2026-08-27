###############################################################################
# INApest pathogen functionality - consolidated regression validation runner
# Release target: 25 August 2026 definitive pathogen-capable source bundle
#
# PURPOSE
#   1. Run the focused regression tests supplied with the definitive bundle.
#   2. Exercise those tests against the CURRENT src/ files, including the
#      serial Meta/MLU integration blocks that were originally written against
#      development filenames.
#   3. Add explicit current-release checks for Pathogen=NULL equivalence,
#      pathogen-detection -> information, public APIs, and serial/parallel Meta
#      parity where the R runtime supports it.
#   4. Write simple machine-readable and human-readable result files that can
#      be uploaded back into the pathogen report workflow.
#
# DEPENDENCIES
#   Base R only. The optional parallel checks use the base/recommended
#   'parallel' package through the INApest parallel wrappers themselves.
#
# HOW TO RUN
#   Easiest: put this script beside either
#     INApest_definitive_pathogen_information_trigger_2026-08-25.zip
#   or the extracted folder
#     INApest_definitive_pathogen_functions_2026-08-25/
#   then run:
#     source("INApest_pathogen_regression_validation_2026-08-25.R")
#   or:
#     Rscript INApest_pathogen_regression_validation_2026-08-25.R
#
#   Optional command-line arguments:
#     --root=/path/to/INApest_definitive_pathogen_functions_2026-08-25
#     --zip=/path/to/INApest_definitive_pathogen_information_trigger_2026-08-25.zip
#     --out=/path/to/output_directory
#
# OUTPUTS
#   pathogen_regression_results.csv
#   pathogen_regression_summary.txt
#   pathogen_regression_sessionInfo.txt
#   logs/<test_id>.log
#
# INTERPRETATION
#   PASS = test executed and all assertions passed.
#   FAIL = test executed and an assertion/runtime error occurred.
#   SKIP = test could not be meaningfully exercised in this runtime (for
#          example, true multi-worker PSOCK parallelism is unavailable).
###############################################################################

options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x) || !nzchar(x)) y else x

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix) {
  hit <- grep(paste0("^", prefix, "="), args, value = TRUE)
  if (!length(hit)) return(NULL)
  sub(paste0("^", prefix, "="), "", hit[1L])
}

user_root <- arg_value("--root")
user_zip  <- arg_value("--zip")
user_out  <- arg_value("--out")

script_dir <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile, mustWork = FALSE))
}, error = function(e) getwd())

release_folder_name <- "INApest_definitive_pathogen_functions_2026-08-25"
release_zip_name <- "INApest_definitive_pathogen_information_trigger_2026-08-25.zip"

find_release_root <- function() {
  if (!is.null(user_root) && dir.exists(user_root)) {
    return(normalizePath(user_root))
  }

  candidates <- unique(c(
    file.path(getwd(), release_folder_name),
    file.path(script_dir, release_folder_name),
    getwd(), script_dir
  ))
  for (x in candidates) {
    if (dir.exists(file.path(x, "src")) && dir.exists(file.path(x, "tests"))) {
      return(normalizePath(x))
    }
  }

  zip_candidates <- unique(c(
    user_zip,
    file.path(getwd(), release_zip_name),
    file.path(script_dir, release_zip_name)
  ))
  zip_candidates <- zip_candidates[!is.na(zip_candidates) & nzchar(zip_candidates)]
  zip_candidates <- zip_candidates[file.exists(zip_candidates)]
  if (length(zip_candidates)) {
    ex <- file.path(tempdir(), paste0("inapest_pathogen_release_", Sys.getpid()))
    dir.create(ex, recursive = TRUE, showWarnings = FALSE)
    utils::unzip(zip_candidates[1L], exdir = ex)
    r <- file.path(ex, release_folder_name)
    if (dir.exists(file.path(r, "src")) && dir.exists(file.path(r, "tests"))) {
      return(normalizePath(r))
    }
  }

  stop(
    "Could not locate the definitive pathogen release. Put this runner beside ",
    release_zip_name, " or the extracted ", release_folder_name,
    " folder, or provide --root= or --zip=."
  )
}

root <- find_release_root()
src_dir <- file.path(root, "src")
test_dir <- file.path(root, "tests")

stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- user_out %||% file.path(getwd(), paste0("INApest_pathogen_validation_", stamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
log_dir <- file.path(out_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

cat("INApest pathogen regression validation\n")
cat("Release root: ", root, "\n", sep = "")
cat("Output dir:   ", normalizePath(out_dir, mustWork = FALSE), "\n\n", sep = "")

# -----------------------------------------------------------------------------
# Stage current sources under the development aliases used by some supplied
# tests. This makes the integration sections exercise the current release
# rather than silently skipping because an older development filename is absent.
# -----------------------------------------------------------------------------
stage_root <- file.path(tempdir(), paste0("inapest_pathogen_stage_", Sys.getpid()))
stage_src <- file.path(stage_root, "src")
stage_tests <- file.path(stage_root, "tests")
dir.create(stage_src, recursive = TRUE, showWarnings = FALSE)
dir.create(stage_tests, recursive = TRUE, showWarnings = FALSE)

src_files <- list.files(src_dir, full.names = TRUE)
test_files <- list.files(test_dir, full.names = TRUE)
file.copy(src_files, stage_src, overwrite = TRUE)
file.copy(test_files, stage_tests, overwrite = TRUE)

aliases <- c(
  "INApestMetaPathogen.r" = "INApestMeta.r",
  "INApestMetaMultipleLandUsePathogen.r" = "INApestMetaMultipleLandUse.r",
  "INApestMetaParallelPathogen.r" = "INApestMetaParallel.r",
  "INApestMetaParallelMultipleLandUsePathogen.r" = "INApestMetaParallelMultipleLandUse.r"
)
for (alias in names(aliases)) {
  from <- file.path(stage_src, aliases[[alias]])
  to <- file.path(stage_src, alias)
  if (file.exists(from)) file.copy(from, to, overwrite = TRUE)
}

results <- data.frame(
  test_id = character(),
  category = character(),
  description = character(),
  status = character(),
  elapsed_seconds = numeric(),
  message = character(),
  log_file = character(),
  stringsAsFactors = FALSE
)

add_result <- function(test_id, category, description, status, elapsed, message = "", log_file = "") {
  results <<- rbind(results, data.frame(
    test_id = test_id,
    category = category,
    description = description,
    status = status,
    elapsed_seconds = round(as.numeric(elapsed), 3),
    message = message,
    log_file = log_file,
    stringsAsFactors = FALSE
  ))
}

run_block <- function(test_id, category, description, code, wd = stage_src) {
  t0 <- proc.time()[[3L]]
  log_file <- file.path(log_dir, paste0(test_id, ".log"))
  oldwd <- getwd()
  on.exit(setwd(oldwd), add = TRUE)
  setwd(wd)

  warnings_seen <- character()
  error_msg <- NULL
  output <- character()

  ok <- tryCatch({
    output <- capture.output(
      withCallingHandlers(
        force(code),
        warning = function(w) {
          warnings_seen <<- c(warnings_seen, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      ),
      type = "output"
    )
    TRUE
  }, error = function(e) {
    error_msg <<- conditionMessage(e)
    FALSE
  })

  elapsed <- proc.time()[[3L]] - t0
  lines <- c(
    paste("TEST:", test_id),
    paste("CATEGORY:", category),
    paste("DESCRIPTION:", description),
    paste("STATUS:", if (ok) "PASS" else "FAIL"),
    if (length(warnings_seen)) c("WARNINGS:", paste0("- ", warnings_seen)) else NULL,
    if (!is.null(error_msg)) c("ERROR:", error_msg) else NULL,
    "OUTPUT:", output
  )
  writeLines(lines, log_file)

  add_result(
    test_id, category, description,
    if (ok) "PASS" else "FAIL", elapsed,
    if (ok) paste(unique(warnings_seen), collapse = " | ") else error_msg,
    file.path("logs", basename(log_file))
  )

  cat(sprintf("%-42s %s\n", test_id, if (ok) "PASS" else "FAIL"))
  invisible(ok)
}

skip_block <- function(test_id, category, description, message) {
  add_result(test_id, category, description, "SKIP", 0, message, "")
  cat(sprintf("%-42s SKIP - %s\n", test_id, message))
  invisible(FALSE)
}

# -----------------------------------------------------------------------------
# Supplied current-release mechanism/integration tests
# -----------------------------------------------------------------------------
run_block(
  "01_binary_mechanisms", "Binary INApest",
  "Binary occupancy constraint, directed transmission, clearance, pathogen introduction, host extinction, time-varying transmission and input-contract checks.",
  sys.source(file.path(stage_tests, "test_binary_pathogen.R"), envir = new.env(parent = .GlobalEnv)),
  wd = stage_src
)

run_block(
  "02_meta_mlu_mechanisms_integration", "INApestMeta / MultipleLandUse",
  "SIS/SIR/SEIR state conservation, recruitment, thinning, progression, recovery, waning, pathogen mortality, pathogen introduction, land-use/time resolution, spatial contact, serial Meta/MLU integration and parallel API checks.",
  sys.source(file.path(stage_tests, "test_INApestPathogen_Meta_MLU.R"), envir = new.env(parent = .GlobalEnv)),
  wd = stage_src
)

run_block(
  "03_point_mechanisms", "INApestMetaPoint",
  "Individual point-state initialization, distance/contact transmission, susceptible default for recruits, progression, recovery, waning, pathogen mortality and pathogen introduction.",
  sys.source(file.path(stage_tests, "test_INApestPointPathogen.R"), envir = new.env(parent = .GlobalEnv)),
  wd = stage_src
)

run_block(
  "04_point_transition_mechanisms", "INApestPointTransitionMatrix",
  "Independence/persistence of demographic stage and pathogen state, stage-specific contact, transmission without demographic corruption, susceptible recruits and pathogen mortality.",
  sys.source(file.path(stage_tests, "test_INApestPointTransitionPathogen.R"), envir = new.env(parent = .GlobalEnv)),
  wd = stage_src
)

run_block(
  "05_node_transition_mechanisms", "INApestMetaTransitionMatrix",
  "Stage x pathogen product-state conservation, demographic progression preserving infection state, cross-stage transmission, management reconciliation and stage-transition movement preserving pathogen state.",
  sys.source(file.path(stage_tests, "test_INApestPathogenTransitionMatrix.R"), envir = new.env(parent = .GlobalEnv)),
  wd = stage_src
)

run_block(
  "06_detection_information_trigger", "Shared response behaviour",
  "Certain pathogen detection is recorded in both trigger modes, but creates INApest information only when DetectionTriggersInfo=TRUE; detection is absent where pathogen is absent.",
  sys.source(file.path(stage_tests, "test_pathogen_detection_triggers_information.R"), envir = new.env(parent = .GlobalEnv)),
  wd = stage_tests
)

# -----------------------------------------------------------------------------
# Extra current-release runtime API/source checks. These are intentionally
# simple: they verify that the exact consolidated parent functions can be
# sourced and expose the common Pathogen argument.
# -----------------------------------------------------------------------------
run_block(
  "07_parent_function_api", "Current consolidated core",
  "Source all pathogen-capable serial parent functions and verify their public Pathogen argument.",
  {
    source(file.path(stage_src, "INApestPathogen.R"))
    source(file.path(stage_src, "INApest.R"))
    source(file.path(stage_src, "INApestMeta.r"))
    source(file.path(stage_src, "INApestMetaMultipleLandUse.r"))
    source(file.path(stage_src, "INApestMetaTransitionMatrix.r"))
    source(file.path(stage_src, "INApestMetaPoint.R"))
    source(file.path(stage_src, "INApestPointTransitionMatrix.R"))

    fns <- c(
      "INApest", "INApestMeta", "INApestMetaMultipleLandUse",
      "INApestMetaTransitionMatrix", "INApestMetaPoint",
      "INApestPointTransitionMatrix"
    )
    missing <- fns[!vapply(fns, exists, logical(1), mode = "function")]
    if (length(missing)) stop("Missing parent functions: ", paste(missing, collapse = ", "))
    bad <- fns[!vapply(fns, function(nm) "Pathogen" %in% names(formals(get(nm))), logical(1))]
    if (length(bad)) stop("Pathogen argument missing from: ", paste(bad, collapse = ", "))
    cat("PASS: all serial pathogen-capable parent functions source and expose Pathogen\n")
  },
  wd = stage_src
)

run_block(
  "08_parallel_function_api", "Parallel wrappers",
  "Source all current parallel wrappers and verify that each exposes a Pathogen argument.",
  {
    source(file.path(stage_src, "INApestPathogen.R"))
    source(file.path(stage_src, "INApestMetaParallel.r"))
    source(file.path(stage_src, "INApestMetaParallelMultipleLandUse.r"))
    source(file.path(stage_src, "INApestMetaTransitionMatrixParallel.r"))
    source(file.path(stage_src, "INApestMetaPointParallel.R"))
    source(file.path(stage_src, "INApestPointTransitionMatrixParallel.R"))

    fns <- c(
      "INApestMetaParallel", "INApestMetaParallelMultipleLandUse",
      "INApestMetaTransitionMatrixParallel", "INApestMetaPointParallel",
      "INApestPointTransitionMatrixParallel"
    )
    missing <- fns[!vapply(fns, exists, logical(1), mode = "function")]
    if (length(missing)) stop("Missing parallel functions: ", paste(missing, collapse = ", "))
    bad <- fns[!vapply(fns, function(nm) "Pathogen" %in% names(formals(get(nm))), logical(1))]
    if (length(bad)) stop("Pathogen argument missing from: ", paste(bad, collapse = ", "))
    cat("PASS: all parallel pathogen-capable wrappers source and expose Pathogen\n")
  },
  wd = stage_src
)

# -----------------------------------------------------------------------------
# Current Pathogen=NULL explicit-vs-omitted equivalence.
# This is a current-release regression check, NOT a comparison to a historical
# pre-pathogen source baseline.
# -----------------------------------------------------------------------------
run_block(
  "09_meta_null_equivalence", "Backward-compatible optionality",
  "Within the current INApestMeta release, explicitly supplying Pathogen=NULL gives the same ordinary host outputs as omitting Pathogen under the same seed.",
  {
    source(file.path(stage_src, "INApestMeta.r"))

    run_meta_null <- function(outdir, explicit_null) {
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      a <- list(
        ModelName = "null_", Nperm = 2, Ntimesteps = 3,
        DetectionProb = 0.05, DetectionSD = 0,
        ManageProb = 0.2, ManageSD = 0,
        MortalityProb = 0.1, MortalitySD = 0,
        FecundityReduction = 0,
        SpreadReduction = 0.1, SpreadReductionSD = 0,
        InitialPopulation = c(10,0,0), InitBioP = NA,
        InvasionRisk = c(0,0,0), InitialInfo = c(0,0,0), InitInfoP = NA,
        ExternalInfoProb = c(0,0,0), InfoRetentionProb = 1,
        InfoPersistenceSteps = NA, EnvEstabProb = 1, Survival = 0.95,
        K = c(100,100,100), PropaguleProduction = 0.5,
        PropaguleEstablishment = 0.1, IncursionStartPop = 1,
        SDDprob = diag(3), SEAM = 0, LDDprob = diag(3), LDDrate = 0,
        OngoingExternalInvasion = FALSE, OngoingExternalInfo = FALSE,
        OutputDir = paste0(outdir, .Platform$file.sep), DoPlots = FALSE
      )
      if (explicit_null) a$Pathogen <- NULL
      set.seed(1234)
      do.call(INApestMeta, a)
      suffixes <- c("PopulationLargeOut.rds", "InvasionLargeOut.rds", "InfoLargeOut.rds", "DetectedLargeOut.rds")
      setNames(lapply(suffixes, function(s) readRDS(file.path(outdir, paste0("null_", s)))), suffixes)
    }

    d1 <- tempfile("meta_omit_")
    d2 <- tempfile("meta_null_")
    x <- run_meta_null(d1, FALSE)
    y <- run_meta_null(d2, TRUE)
    for (nm in names(x)) {
      if (!isTRUE(all.equal(x[[nm]], y[[nm]], check.attributes = TRUE))) {
        stop("Pathogen=NULL explicit-vs-omitted mismatch for ", nm)
      }
    }
    cat("PASS: current INApestMeta Pathogen=NULL explicit-vs-omitted equivalence\n")
  },
  wd = stage_src
)

# -----------------------------------------------------------------------------
# Meta serial vs parallel wrapper parity under deterministic pathogen settings.
# In runtimes with one detected core the parallel wrapper legitimately falls
# back to lapply; we record that fact in the log/summary and do NOT call it a
# true multi-worker PSOCK validation.
# -----------------------------------------------------------------------------
run_block(
  "10_meta_serial_parallel_wrapper_parity", "Parallel wrappers",
  "Compare serial INApestMeta with INApestMetaParallel under deterministic pathogen settings. Runtime core count is reported separately so true PSOCK multi-worker evidence can be distinguished from wrapper fallback.",
  {
    source(file.path(stage_src, "INApestPathogen.R"))
    source(file.path(stage_src, "INApestMeta.r"))
    serial_fun <- INApestMeta
    source(file.path(stage_src, "INApestMetaParallel.r"))
    parallel_fun <- INApestMetaParallel

    d1 <- tempfile("meta_serial_")
    d2 <- tempfile("meta_parallel_")
    dir.create(d1); dir.create(d2)

    path <- INApestPathogen("SIR", Beta = 0, RecoveryProb = 0,
                            InitialInfected = c(2,0,0))
    common <- list(
      ModelName = "par_", Nperm = 2, Ntimesteps = 2, Pathogen = path,
      DetectionProb = 0, DetectionSD = 0,
      ManageProb = 0, ManageSD = 0,
      MortalityProb = 0, MortalitySD = 0,
      FecundityReduction = 0,
      SpreadReduction = 0, SpreadReductionSD = 0,
      InitialPopulation = c(10,0,0), InitBioP = NA,
      InvasionRisk = c(0,0,0), InitialInfo = c(0,0,0), InitInfoP = NA,
      ExternalInfoProb = c(0,0,0), InfoRetentionProb = 1,
      InfoPersistenceSteps = NA, EnvEstabProb = 0, Survival = 1,
      K = c(100,100,100), PropaguleProduction = 0,
      PropaguleEstablishment = 0, IncursionStartPop = 1,
      SDDprob = diag(3), SEAM = 0, LDDprob = diag(3), LDDrate = 0,
      OngoingExternalInvasion = FALSE, OngoingExternalInfo = FALSE,
      DoPlots = FALSE
    )

    a1 <- common; a1$OutputDir <- paste0(d1, .Platform$file.sep)
    a2 <- common; a2$OutputDir <- paste0(d2, .Platform$file.sep)
    set.seed(777); do.call(serial_fun, a1)
    set.seed(777); do.call(parallel_fun, a2)

    suffixes <- c("PopulationLargeOut.rds", "PathogenStateLargeOut.rds", "InvasionLargeOut.rds")
    for (s in suffixes) {
      x <- readRDS(file.path(d1, paste0("par_", s)))
      y <- readRDS(file.path(d2, paste0("par_", s)))
      if (!isTRUE(all.equal(x, y, check.attributes = TRUE))) {
        stop("Serial/parallel wrapper mismatch for ", s)
      }
    }

    dc <- tryCatch(parallel::detectCores(), error = function(e) NA_integer_)
    nc <- if (is.na(dc)) 1L else max(1L, min(2L, dc - 1L))
    cat("Detected cores: ", dc, "\n", sep = "")
    cat("INApestMetaParallel effective worker count for Nperm=2: ", nc, "\n", sep = "")
    if (nc > 1L) {
      cat("PASS: this run exercised a true multi-worker PSOCK path\n")
    } else {
      cat("PASS: wrapper parity passed, but this runtime used the single-worker fallback\n")
    }
  },
  wd = stage_src
)

# -----------------------------------------------------------------------------
# Emit results and session information.
# -----------------------------------------------------------------------------
write.csv(results, file.path(out_dir, "pathogen_regression_results.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "pathogen_regression_sessionInfo.txt"))

n_pass <- sum(results$status == "PASS")
n_fail <- sum(results$status == "FAIL")
n_skip <- sum(results$status == "SKIP")

summary_lines <- c(
  "INApest pathogen regression validation summary",
  paste0("Date/time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Release root: ", root),
  paste0("R version: ", R.version.string),
  paste0("Platform: ", R.version$platform),
  paste0("Tests: ", nrow(results), " total; ", n_pass, " PASS; ", n_fail, " FAIL; ", n_skip, " SKIP"),
  "",
  "Important evidence boundary:",
  "- Test 09 checks Pathogen=NULL equivalence within the current release. It is not a comparison to a historical pre-pathogen source baseline.",
  "- Test 10 records whether a true multi-worker PSOCK path was available. A single-worker fallback can validate wrapper logic but is not evidence of multi-worker parallel execution.",
  "- The supplied historical test_Pathogen_NULL_and_parallel.R is not run directly because the release bundle does not contain its baseline/INApestMeta.r dependency and it refers to older development filenames.",
  "",
  "Per-test results:",
  paste0(sprintf("%-42s", results$test_id), "  ", results$status, ifelse(nzchar(results$message), paste0("  ", results$message), ""))
)
writeLines(summary_lines, file.path(out_dir, "pathogen_regression_summary.txt"))

cat("\n", paste(rep("=", 72), collapse = ""), "\n", sep = "")
cat("Validation complete\n")
cat("PASS: ", n_pass, "  FAIL: ", n_fail, "  SKIP: ", n_skip, "\n", sep = "")
cat("Results: ", normalizePath(out_dir, mustWork = FALSE), "\n", sep = "")
cat(paste(rep("=", 72), collapse = ""), "\n", sep = "")

if (n_fail > 0L) {
  stop("One or more pathogen regression validations FAILED. See pathogen_regression_summary.txt and logs/.")
}
