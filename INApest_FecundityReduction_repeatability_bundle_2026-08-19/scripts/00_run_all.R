# Entry point for the FecundityReduction repeatability bundle.
# Run from the bundle root: Rscript scripts/00_run_all.R

bundle_root <- r"[C:\Users\MasonN\OneDrive - MWLR\MPI sleeper pests\BSI SSIF\INApestDevelopment\INApest_FecundityReduction_repeatability_bundle_2026-08-19]"
Sys.setenv(INAPEST_BUNDLE_ROOT=bundle_root)
cores <- parallel::detectCores()

if (is.na(cores) || cores <= 1L) {
  # This reproduces the archived validation path, including the parallel
  # functions' one-core fallback, and should allow exact archived comparisons.
  rerun_dir <- file.path(bundle_root, "results", "rerun")
  dir.create(rerun_dir, recursive=TRUE, showWarnings=FALSE)
  Sys.setenv(INAPEST_TEST_OUT=rerun_dir)
  source(file.path(bundle_root, "scripts", "01_validate_fecundity_reduction_current.R"))
  source(file.path(bundle_root, "scripts", "02_current_meta_serial_parallel_smoke.R"))
  source(file.path(bundle_root, "scripts", "03_compare_archived_and_rerun.R"))
  source(file.path(bundle_root, "scripts", "04_plot_results.R"))
  cat("\nPASS: one-core archived validation workflow completed.\n")
  cat("Rerun results:", rerun_dir, "\n")
} else {
  # On a normal multi-core machine, parallel RNG streams are intentionally
  # different from the archived one-core run. Reproduce the scientific tables
  # with serial current-source functions, then validate true PSOCK separately.
  rerun_dir <- file.path(bundle_root, "results", "native_serial_reproduction")
  dir.create(rerun_dir, recursive=TRUE, showWarnings=FALSE)
  Sys.setenv(INAPEST_TEST_OUT=rerun_dir)
  source(file.path(bundle_root, "scripts", "06_reproduce_scientific_results_serial.R"))
  source(file.path(bundle_root, "scripts", "04_plot_results.R"))
  cat("\nPASS: native serial scientific reproduction completed.\n")
  cat("Scientific rerun results:", rerun_dir, "\n")
  cat("Now run scripts/05_native_psock_validation_current.R separately to force and test the multi-worker PSOCK branch.\n")
}
