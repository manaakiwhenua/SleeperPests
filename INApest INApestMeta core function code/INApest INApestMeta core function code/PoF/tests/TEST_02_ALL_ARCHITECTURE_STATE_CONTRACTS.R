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

pass <- function(label) cat("PASS -", label, "\n")
check <- function(ok, label) { if (!isTRUE(ok)) stop("FAIL - ", label, call.=FALSE); pass(label) }

expected <- matrix(c(TRUE,FALSE,TRUE,
                     TRUE,TRUE,FALSE), nrow=2, byrow=TRUE)

make_node <- function() {
  x <- array(0L, dim=c(2,2,2,3), dimnames=list(NULL,c("juvenile","adult"),NULL,NULL))
  x[1,"adult",1,2] <- 1L
  x[2,"juvenile",2,3] <- 2L
  list(BiocontrolHistory=list(State=list(Q=x),Attacks=list(),Recruits=list()))
}

make_point <- function() {
  h <- expand.grid(perm=1:3,timestep=1:2,agent="Q",node=1:2,
                   stage=c("juvenile","adult"), KEEP.OUT.ATTRS=FALSE,
                   stringsAsFactors=FALSE)
  h$abundance <- 0L
  h$abundance[h$perm==2 & h$timestep==1 & h$node==1 & h$stage=="adult"] <- 1L
  h$abundance[h$perm==3 & h$timestep==2 & h$node==2 & h$stage=="juvenile"] <- 2L
  s <- expand.grid(perm=1:3,timestep=1:2)
  list(BiocontrolHistory=h, Summary=s)
}

node_arch <- c("Binary", "Meta", "Multiple Land Use", "Transition Matrix", "Vertebrate Node")
point_arch <- c("MetaPoint", "Point Transition Matrix", "Vertebrate Point")

for (arch in node_arch) {
  m <- make_node()
  fs <- INApestBiocontrolFreedomState(m)
  check(identical(unname(fs$Freedom), expected), paste(arch, "biocontrol freedom extraction"))
  surv <- INApestBiocontrolSurveillance(0.5, DetectStages="adult")
  l1 <- INApestBiocontrolObservationLikelihood(m, surv, 1, 0)
  check(isTRUE(all.equal(l1, c(1,0.5,1), tolerance=1e-12)), paste(arch, "adult observation likelihood"))
  l2 <- INApestBiocontrolObservationLikelihood(m, surv, 2, 0)
  check(isTRUE(all.equal(l2, c(1,1,1), tolerance=1e-12)), paste(arch, "hidden juvenile observation"))
}

for (arch in point_arch) {
  m <- make_point()
  fs <- INApestBiocontrolFreedomState(m)
  check(identical(unname(fs$Freedom), expected), paste(arch, "biocontrol freedom extraction"))
  surv <- INApestBiocontrolSurveillance(0.5, DetectStages="adult")
  l1 <- INApestBiocontrolObservationLikelihood(m, surv, 1, 0)
  check(isTRUE(all.equal(l1, c(1,0.5,1), tolerance=1e-12)), paste(arch, "adult observation likelihood"))
  l2 <- INApestBiocontrolObservationLikelihood(m, surv, 2, 0)
  check(isTRUE(all.equal(l2, c(1,1,1), tolerance=1e-12)), paste(arch, "hidden juvenile observation"))
}

cat("\nArchitectures checked:", length(node_arch)+length(point_arch), "\n")
cat("OVERALL: PASS\n")
