###############################################################################
# INApestMLURK4LocalDynamics.R
# Host-only stochastic continuous-time LocalDynamics adaptor for MLU
# v0.6 -- 30 September 2026
###############################################################################
if (!exists("INApestStochasticLocalDynamics", mode="function"))
  stop("Source INApestRK4.R and INApestStochasticFluxBridge.R first")

INApestMLURK4LocalDynamics <- function(
    HostFluxFunction,
    TimestepLength=1,
    RKMaxStep=0.025,
    Parameters=NULL,
    StartTime=0,
    DispersalDynamics=NULL,
    RunWhenEmpty=FALSE,
    CapacityTolerance=1e-9) {
  if(!is.function(HostFluxFunction)) stop("HostFluxFunction must be a function")
  if(is.null(DispersalDynamics)) {
    if(!exists("local.dynamicsLU",mode="function")) stop("Source current INApestMetaMultipleLandUse.r first")
    DispersalDynamics <- get("local.dynamicsLU",mode="function")
  }
  if(!is.function(DispersalDynamics)) stop("DispersalDynamics must be a function")
  if(!is.logical(RunWhenEmpty)||length(RunWhenEmpty)!=1L||is.na(RunWhenEmpty)) stop("RunWhenEmpty must be TRUE/FALSE")
  cap_tol <- as.numeric(CapacityTolerance)[1L]
  if(!is.finite(cap_tol)||cap_tol<0) stop("CapacityTolerance must be finite and >=0")
  stochastic_local <- INApestStochasticLocalDynamics(HostFluxFunction,TimestepLength,RKMaxStep,Parameters,StartTime)
  disperse <- DispersalDynamics
  f <- function(sddprob,nodepropaguleproduction,nodeenvestabprob,n,lddprob,lddrate,
                k_is_0,nodeK,nodepropaguleestablishment,nodespreadreduction,
                nodefecundityreduction=0,managing,timestep=NULL,Ntimesteps=NULL) {
    if(!is.matrix(n)||!is.matrix(nodeK)||!identical(dim(n),dim(nodeK))) stop("MLU host state and nodeK must be matching node x land-use matrices")
    n_after <- stochastic_local(n0=n,timestep=timestep,Ntimesteps=Ntimesteps,
      nodeK=nodeK,k_is_0=k_is_0,nodeenvestabprob=nodeenvestabprob,
      nodepropaguleproduction=nodepropaguleproduction,nodepropaguleestablishment=nodepropaguleestablishment,
      nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,
      managing=managing,sddprob=sddprob,lddprob=lddprob,lddrate=lddrate)
    if(!is.matrix(n_after)||!identical(dim(n_after),dim(n))) stop("MLU RK host bridge did not preserve node x land-use shape")
    over <- which(n_after > nodeK + cap_tol*pmax(1,abs(nodeK)),arr.ind=TRUE)
    if(nrow(over)) stop("RK host state exceeded MLU nodeK before dispersal at cell(s): ",paste(apply(over,1,paste,collapse=":"),collapse=", "))
    out <- disperse(sddprob=sddprob,nodepropaguleproduction=nodepropaguleproduction,
      nodeenvestabprob=nodeenvestabprob,n=n_after,lddprob=lddprob,lddrate=lddrate,
      k_is_0=k_is_0,nodeK=nodeK,nodepropaguleestablishment=nodepropaguleestablishment,
      nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,
      managing=managing)
    attr(out,"INApestMLURK4") <- list(StateAfterRK=n_after,Flux=attr(n_after,"INApestStochasticFlux"),Timestep=timestep,Ntimesteps=Ntimesteps)
    out
  }
  attr(f,"INApestMLURK4LocalDynamics") <- list(Version="0.6",StateMode="integer",TimestepLength=TimestepLength,RKMaxStep=RKMaxStep,RunWhenEmpty=RunWhenEmpty)
  attr(f,"INApestRunWhenEmpty") <- RunWhenEmpty
  class(f) <- c("INApestMLURK4LocalDynamics","function")
  f
}
