bundle_root <- normalizePath(Sys.getenv("INAPEST_BUNDLE_ROOT", unset = "."), mustWork = TRUE)
archived <- file.path(bundle_root, "results", "archived")
rerun <- Sys.getenv("INAPEST_TEST_OUT", unset = file.path(bundle_root, "results", "rerun"))
files <- c(
  "equivalence_results.csv",
  "parameterisation_smoke_results.csv",
  "parameterisation_parallel_equivalence.csv",
  "mortality_by_fecundity_results.csv",
  "landuse_mortality_fecundity_results.csv",
  "validation_environment.csv"
)
res <- do.call(rbind, lapply(files, function(f) {
  a <- readLines(file.path(archived, f), warn=FALSE)
  b <- readLines(file.path(rerun, f), warn=FALSE)
  data.frame(file=f, exact_text_match=identical(a,b), stringsAsFactors=FALSE)
}))
write.csv(res, file.path(rerun, "archived_vs_rerun_R_comparison.csv"), row.names=FALSE)
print(res)
if (!all(res$exact_text_match)) stop("At least one rerun CSV differs from the archived result.")
