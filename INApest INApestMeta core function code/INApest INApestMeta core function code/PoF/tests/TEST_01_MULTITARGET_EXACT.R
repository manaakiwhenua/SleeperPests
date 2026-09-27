options(stringsAsFactors = FALSE)
script_path <- function() {
  z <- grep("^--file=", commandArgs(FALSE), value=TRUE)
  if (length(z)) return(normalizePath(sub("^--file=", "", z[1L]), winslash="/", mustWork=TRUE))
  if (!is.null(sys.frame(1)$ofile)) return(normalizePath(sys.frame(1)$ofile, winslash="/", mustWork=TRUE))
  stop("Could not determine script path.")
}
test_dir <- dirname(script_path())
root <- normalizePath(file.path(test_dir, ".."), winslash="/", mustWork=TRUE)
source(file.path(root, "src", "INApestProofOfFreedom.R"))

assert_close <- function(x, y, tol = 1e-12, label = "") {
  if (!isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol, check.attributes = FALSE)))
    stop("FAIL - ", label, ": got ", paste(x, collapse=","), " expected ", paste(y, collapse=","), call. = FALSE)
  cat("PASS -", label, "\n")
}
assert_true <- function(x, label) {
  if (!isTRUE(x)) stop("FAIL - ", label, call. = FALSE)
  cat("PASS -", label, "\n")
}

# Four correlated particles:
# 1 both free; 2 pathogen present only; 3 biocontrol present only; 4 both present.
path_present <- array(c(0,1,0,1), dim = c(1,1,4))
q_state <- array(0L, dim = c(1,2,1,4),
                 dimnames = list(NULL, c("juvenile","adult"), NULL, NULL))
q_state[1,"adult",1,c(3,4)] <- 1L
model <- list(
  PathogenPresentResults = path_present,
  BiocontrolHistory = list(State=list(Q=q_state), Attacks=list(), Recruits=list())
)
pathogen <- structure(list(Model="Binary", DetectionProb=0.5,
                           DetectionTriggersInfo=FALSE), class="INApestPathogen")
surv <- INApestBiocontrolSurveillance(DetectionProb=0.25, DetectStages="adult")
prior <- c(0.4,0.1,0.1,0.4)

fit <- INApestPoF(
  model, Pathogen=pathogen, Targets=c("pathogen","biocontrol"),
  ObservationHistory=list(pathogen=0, biocontrol=0),
  PriorWeights=prior, BiocontrolSurveillance=surv
)

# Conditional no-detection likelihoods:
# pathogen = [1, .5, 1, .5]; biocontrol = [1, 1, .75, .75].
# Combined = [1, .5, .75, .375]; evidence = .675.
expected_w <- c(0.4, 0.05, 0.075, 0.15) / 0.675
assert_close(fit$Summary$ObservationEvidence, 0.675, label="joint observation evidence")
assert_close(fit$PosteriorWeights, expected_w, label="joint posterior weights")
assert_close(fit$Summary$PosteriorPoF_Pathogen, sum(expected_w[c(1,3)]), label="pathogen marginal PoF")
assert_close(fit$Summary$PosteriorPoF_Biocontrol, sum(expected_w[c(1,2)]), label="biocontrol marginal PoF")
assert_close(fit$Summary$PosteriorPoF_Joint, expected_w[1], label="joint pathogen+biocontrol PoF")
assert_true(abs(fit$Summary$PosteriorPoF_Joint -
                fit$Summary$PosteriorPoF_Pathogen * fit$Summary$PosteriorPoF_Biocontrol) > 0.1,
            "joint PoF is not product of correlated marginals")

lp <- c(1,0.5,1,0.5); lq <- c(1,1,0.75,0.75)
w_pq <- .ipof_norm_weights(.ipof_norm_weights(prior * lp) * lq)
w_qp <- .ipof_norm_weights(.ipof_norm_weights(prior * lq) * lp)
assert_close(w_pq, w_qp, label="same-timestep update order independence")
assert_close(w_pq, fit$PosteriorWeights, label="combined-likelihood equivalence")

# Hidden-stage benchmark: juveniles count against freedom even when only adults
# are detectable. A juvenile-only particle therefore has no-detection likelihood
# one but is biologically not free.
q2 <- array(0L, dim=c(1,2,1,2), dimnames=list(NULL,c("juvenile","adult"),NULL,NULL))
q2[1,"juvenile",1,2] <- 1L
m2 <- list(BiocontrolHistory=list(State=list(Q=q2),Attacks=list(),Recruits=list()))
fs2 <- INApestBiocontrolFreedomState(m2)
l2 <- INApestBiocontrolObservationLikelihood(m2, surv, 1, 0)
assert_true(identical(as.logical(fs2$Freedom[1,]), c(TRUE,FALSE)), "juvenile-only state is not free")
assert_close(l2, c(1,1), label="juvenile-only state is invisible to adult-only surveillance")

cat("\nOVERALL: PASS\n")
