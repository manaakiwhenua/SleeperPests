###############################################################################
### Regression checks for pathogen PoF integration patch
###############################################################################

root <- normalizePath(file.path(getwd()), mustWork = TRUE)
source(file.path(root, "src", "INApestPathogen.R"))
source(file.path(root, "src", "INApestPointPathogen.R"))
source(file.path(root, "src", "INApestPathogenPoF.R"))

assert <- function(x, msg) if (!isTRUE(x)) stop(msg, call. = FALSE)

### 1. Saved-output filename stem and in-memory object must agree.
p <- INApestPathogen(Model = "SEIR", Beta = 0, RecoveryProb = 0,
                     ProgressionProb = 0, DetectionProb = 0.5)
state <- array(0, dim = c(2, 4, 1, 2),
               dimnames = list(NULL, c("S","E","I","R"), NULL, NULL))
state[,"S",1,] <- 1
state[1,"S",1,2] <- 0
state[1,"I",1,2] <- 1
obj <- list(PathogenStateResults = state)

td <- tempfile("ipof_saved_"); dir.create(td)
stem <- file.path(td, "demo")
saveRDS(state, paste0(stem, "PathogenStateLargeOut.rds"))

mem <- INApestPathogenPoF(obj, p, ObservationHistory = list(0))
disk <- INApestPathogenPoF(stem, p, ObservationHistory = list(0))
disk2 <- INApestPathogenPoF(list(OutputDir = td, ModelName = "demo"), p,
                            ObservationHistory = list(0))
assert(isTRUE(all.equal(mem$Summary, disk$Summary, tolerance = 0)),
       "Saved-output stem differs from in-memory PoF result.")
assert(isTRUE(all.equal(mem$Summary, disk2$Summary, tolerance = 0)),
       "OutputDir/ModelName input differs from in-memory PoF result.")

### 2. Point PoF must use the exact point helper resolver.
p2 <- INApestPathogen(Model = "SIR", Beta = 0, RecoveryProb = 0,
                      DetectionProb = function(points, timestep, perm)
                        ifelse(points$x < 0.5, 0.9, 0.2))
pi <- INApestPointPathogenInteraction(p2)
assert(is.function(pi$ResolvePoint), "Point interaction does not expose ResolvePoint.")
h <- data.frame(perm = c(1,1), timestep = c(1,1), id = c(1,2),
                x = c(0,1), y = c(0,0), pathogen_state = c("I","I"))
po <- list(PointHistory = h)
L <- INApestPathogenObservationLikelihood(po, pi, 1, 0)
assert(abs(L - 0.08) < 1e-12, "Point resolver likelihood should equal 0.08.")

### 3. Serial node transition source must persist pathogen detections.
tm <- paste(readLines(file.path(root, "src", "INApestMetaTransitionMatrix.r"), warn = FALSE), collapse = "\n")
assert(grepl("PathogenDetectedLargeOut\\.rds", tm),
       "Serial transition-matrix source does not save PathogenDetectedLargeOut.rds.")

cat("Pathogen PoF integration patch regression checks PASS\n")
