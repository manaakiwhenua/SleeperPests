###############################################################################
# INApestMetaRK4HostPathogenLocalDynamics.R
# Coupled stochastic continuous-time host + pathogen LocalDynamics for Meta
# v0.4 -- 30 September 2026
###############################################################################

if (!exists("INApestStochasticCompartmentIntegrate", mode="function"))
  stop("Source INApestRK4.R and INApestStochasticCompartmentBridge.R first")

INApestContinuousPathogenRates <- function(
    RecoveryRate = 0,
    ProgressionRate = 0,
    PathogenMortalityRate = 0,
    ImmunityLossRate = 0) {
  out <- list(RecoveryRate=RecoveryRate, ProgressionRate=ProgressionRate,
              PathogenMortalityRate=PathogenMortalityRate,
              ImmunityLossRate=ImmunityLossRate)
  class(out) <- c("INApestContinuousPathogenRates","list")
  out
}

.inapest_hpr_host_flux <- function(fun,t,N,Parameters=NULL,runtime=list()) {
  ans <- .inapest_rk4_call_supported(
    fun,c(list(t=t,time=t,state=N,State=N,pars=Parameters,Parameters=Parameters),runtime)
  )
  if (!is.list(ans)) stop("HostFluxFunction must return GainRate and Hazards")
  gain <- ans$GainRate; if (is.null(gain)) gain <- ans$Gains
  if (is.null(gain)) stop("HostFluxFunction must return GainRate")
  gain <- as.numeric(gain)
  if (length(gain)==1L) gain <- rep(gain,length(N))
  if (length(gain)!=length(N) || any(!is.finite(gain)) || any(gain<0))
    stop("Host GainRate must resolve to finite non-negative values per node")
  H <- ans$Hazards
  if (is.null(H) || (is.list(H) && length(H)==0L))
    return(list(GainRate=gain,Hazards=matrix(numeric(),length(N),0L)))
  if (is.list(H)) {
    if (is.null(names(H)) || any(!nzchar(names(H))) || anyDuplicated(names(H)))
      stop("Host hazard list must have unique non-empty names")
    mat <- vapply(H,function(z) {
      z <- as.numeric(z); if(length(z)==1L) z <- rep(z,length(N))
      if(length(z)!=length(N)) stop("Each host hazard must be scalar or length nodes")
      z
    },numeric(length(N)))
    if (is.null(dim(mat))) mat <- matrix(mat,ncol=1L,dimnames=list(NULL,names(H)))
    colnames(mat) <- names(H)
  } else {
    mat <- as.matrix(H)
    if(nrow(mat)==1L && length(N)>1L) mat <- mat[rep(1L,length(N)),,drop=FALSE]
    if(nrow(mat)!=length(N)) stop("Host Hazards must have one row per node")
    if(is.null(colnames(mat))) colnames(mat) <- paste0("loss",seq_len(ncol(mat)))
  }
  if(any(!is.finite(mat)) || any(mat<0)) stop("Host Hazards must be finite and non-negative")
  list(GainRate=gain,Hazards=mat)
}

.inapest_hpr_rate_model <- function(t, State, HostFluxFunction, HostParameters,
                                    PathogenRates, Pathogen, PathogenEngine,
                                    PathogenContext, timestep, Ntimesteps,
                                    runtime=list()) {
  comps <- colnames(State); n <- nrow(State); k <- ncol(State)
  if (!all(c("S","I") %in% comps)) stop("Coupled pathogen state must contain S and I")
  N <- rowSums(State)
  hf <- .inapest_hpr_host_flux(HostFluxFunction,t,N,HostParameters,runtime)
  gain <- matrix(0,n,k,dimnames=dimnames(State)); gain[,"S"] <- hf$GainRate
  trans <- array(0,c(n,k,k),dimnames=list(rownames(State),comps,comps))

  resolve_rate <- function(x,name) {
    z <- PathogenEngine$Resolve(x,timestep,PathogenContext,name)
    z <- as.numeric(z)
    if(length(z)!=n || any(!is.finite(z)) || any(z<0)) stop(name," must resolve to finite non-negative rates per node")
    z
  }
  rec <- resolve_rate(PathogenRates$RecoveryRate,"RecoveryRate")
  prog <- resolve_rate(PathogenRates$ProgressionRate,"ProgressionRate")
  pmort <- resolve_rate(PathogenRates$PathogenMortalityRate,"PathogenMortalityRate")
  wan <- resolve_rate(PathogenRates$ImmunityLossRate,"ImmunityLossRate")
  beta <- as.numeric(PathogenEngine$Resolve(Pathogen$Beta,timestep,PathogenContext,"Beta"))
  if(length(beta)!=n || any(!is.finite(beta)) || any(beta<0)) stop("Beta must resolve to finite non-negative rates per node")
  density_scale <- as.numeric(PathogenEngine$Resolve(Pathogen$DensityScale,timestep,PathogenContext,"DensityScale"))
  if(length(density_scale)!=n || any(!is.finite(density_scale)) || any(density_scale<=0)) stop("DensityScale must resolve positive values per node")
  C <- PathogenEngine$ContactMatrix(timestep,PathogenContext)
  I <- State[,"I"]
  infectious_pressure <- as.numeric(crossprod(I,C))
  if (Pathogen$Transmission == "frequency") {
    contact_population <- as.numeric(crossprod(N,C))
    pressure <- ifelse(contact_population>0,infectious_pressure/contact_population,0)
    lambda <- beta*pressure
  } else lambda <- beta*infectious_pressure/density_scale
  lambda <- pmax(0,lambda)

  if ("E" %in% comps) {
    trans[,"S","E"] <- lambda
    trans[,"E","I"] <- prog
  } else {
    if(any(prog>0)) stop("ProgressionRate must be zero for SIS/SIR")
    trans[,"S","I"] <- lambda
  }
  if (Pathogen$Model == "SIS") {
    trans[,"I","S"] <- rec
    if(any(wan>0)) stop("ImmunityLossRate must be zero for SIS")
  } else {
    if (!("R" %in% comps)) stop("SIR/SEIR coupled state requires R")
    trans[,"I","R"] <- rec
    trans[,"R","S"] <- wan
  }

  host_causes <- colnames(hf$Hazards)
  cause_names <- c(if(length(host_causes)) paste0("host:",host_causes) else character(),
                   if(any(pmort>0)) "pathogen" else character())
  exits <- array(0,c(n,k,length(cause_names)),dimnames=list(rownames(State),comps,cause_names))
  if(length(host_causes)) for(j in seq_along(host_causes))
    exits[,,paste0("host:",host_causes[j])] <- matrix(rep(hf$Hazards[,j],k),n,k)
  if(any(pmort>0)) exits[,"I","pathogen"] <- pmort

  list(GainRate=gain,TransitionHazards=trans,ExitHazards=exits)
}

INApestMetaRK4HostPathogenLocalDynamics <- function(
    HostFluxFunction,
    PathogenRates = INApestContinuousPathogenRates(),
    TimestepLength = 1,
    RKMaxStep = 0.025,
    HostParameters = NULL,
    StartTime = 0,
    DispersalDynamics = NULL,
    RunWhenEmpty = FALSE,
    CapacityTolerance = 1e-9) {
  if(!is.function(HostFluxFunction)) stop("HostFluxFunction must be a function")
  if(!inherits(PathogenRates,"INApestContinuousPathogenRates")) stop("PathogenRates must come from INApestContinuousPathogenRates()")
  TimestepLength <- as.numeric(TimestepLength)[1L]; RKMaxStep <- as.numeric(RKMaxStep)[1L]; StartTime <- as.numeric(StartTime)[1L]
  if(!is.finite(TimestepLength)||TimestepLength<=0) stop("TimestepLength must be finite and >0")
  if(!is.finite(RKMaxStep)||RKMaxStep<=0) stop("RKMaxStep must be finite and >0")
  if(!is.finite(StartTime)) stop("StartTime must be finite")
  if(is.null(DispersalDynamics)) {
    if(!exists("local.dynamics",mode="function")) stop("Source current INApestMeta.r before constructing this adaptor")
    DispersalDynamics <- get("local.dynamics",mode="function")
  }
  if(!is.function(DispersalDynamics)) stop("DispersalDynamics must be a function")
  if(!is.logical(RunWhenEmpty)||length(RunWhenEmpty)!=1L||is.na(RunWhenEmpty)) stop("RunWhenEmpty must be TRUE/FALSE")
  CapacityTolerance <- as.numeric(CapacityTolerance)[1L]
  if(!is.finite(CapacityTolerance)||CapacityTolerance<0) stop("CapacityTolerance must be finite and >=0")

  host_fun <- HostFluxFunction; rates <- PathogenRates; pars <- HostParameters
  dt <- TimestepLength; hmax <- RKMaxStep; t0 <- StartTime; disperse <- DispersalDynamics; cap_tol <- CapacityTolerance

  f <- function(sddprob,nodepropaguleproduction,nodeenvestabprob,n,lddprob,lddrate,
                k_is_0,nodeK,nodepropaguleestablishment,nodespreadreduction,
                nodefecundityreduction=0,managing,maxinteger,
                pathogen_state,pathogen,pathogen_engine,pathogen_context,
                timestep=NULL,Ntimesteps=NULL) {
    if(is.null(pathogen) || !inherits(pathogen,"INApestPathogen")) stop("Coupled RK LocalDynamics requires Pathogen")
    if(pathogen$Model=="Binary") stop("Binary pathogen mode is not supported by continuous compartment coupling")
    if(!is.matrix(pathogen_state) || any(rowSums(pathogen_state)!=as.integer(n)))
      stop("pathogen_state must be a compartment matrix whose row sums equal host abundance entering coupled RK")
    step_index <- if(is.null(timestep)) 1L else as.integer(timestep)[1L]
    if(is.na(step_index)||step_index<1L) stop("timestep must be >=1")
    nts <- if(is.null(Ntimesteps)) 1L else as.integer(Ntimesteps)[1L]
    intro <- pathogen_engine$Resolve(pathogen$IntroductionProb,step_index,pathogen_context,"IntroductionProb")
    if(any(intro!=0)) stop("Continuous Meta host-pathogen coupling v0.4 requires Pathogen$IntroductionProb = 0; pathogen introduction is deferred to a later extension")
    runtime <- list(nodeK=nodeK,k_is_0=k_is_0,nodeenvestabprob=nodeenvestabprob,
                    nodepropaguleproduction=nodepropaguleproduction,
                    nodepropaguleestablishment=nodepropaguleestablishment,
                    nodespreadreduction=nodespreadreduction,
                    nodefecundityreduction=nodefecundityreduction,managing=managing,
                    sddprob=sddprob,lddprob=lddprob,lddrate=lddrate,
                    timestep=step_index,Ntimesteps=nts)
    rate_fun <- function(t,state,...) .inapest_hpr_rate_model(
      t,state,host_fun,pars,rates,pathogen,pathogen_engine,pathogen_context,
      step_index,nts,runtime)
    z <- INApestStochasticCompartmentIntegrate(
      State=pathogen_state,CompartmentRateFunction=rate_fun,Duration=dt,
      MaxStep=hmax,StartTime=t0+(step_index-1L)*dt)
    p_after <- z$State; n_after <- rowSums(p_after)

    Kvec <- as.numeric(nodeK); if(length(Kvec)==1L) Kvec <- rep(Kvec,length(n_after))
    if(length(Kvec)!=length(n_after)||any(!is.finite(Kvec))||any(Kvec<0)) stop("nodeK must resolve finite non-negative values per Meta node")
    over <- which(n_after > Kvec + cap_tol*pmax(1,abs(Kvec)))
    if(length(over)) stop("Coupled RK host state exceeded Meta nodeK before dispersal at node(s): ",paste(over,collapse=", "),". The adaptor does not silently clip.")

    n_out <- disperse(sddprob=sddprob,nodepropaguleproduction=nodepropaguleproduction,
      nodeenvestabprob=nodeenvestabprob,n=n_after,lddprob=lddprob,lddrate=lddrate,
      k_is_0=k_is_0,nodeK=nodeK,nodepropaguleestablishment=nodepropaguleestablishment,
      nodespreadreduction=nodespreadreduction,nodefecundityreduction=nodefecundityreduction,
      managing=managing,maxinteger=maxinteger)
    p_out <- pathogen_engine$Reconcile(p_after,n_out,pathogen_context)
    if(any(rowSums(p_out)!=as.integer(n_out))) stop("Coupled RK dispersal reconciliation violated pathogen-state accounting")
    list(N=as.integer(n_out),PathogenState=p_out,
         CoupledRK=list(PreDispersalHost=as.integer(n_after),PreDispersalPathogen=p_after,
                        Flux=z,Timestep=step_index,Ntimesteps=nts))
  }
  attr(f,"INApestMetaRK4HostPathogenLocalDynamics") <- list(
    Version="0.4",StateMode="integer-compartment",TimestepLength=dt,RKMaxStep=hmax,
    RunWhenEmpty=RunWhenEmpty,Ordering=c("parent survival/management mortality",
      "joint RK stochastic host/pathogen biology","Meta dispersal/establishment",
      "external invasion","biocontrol/surveillance downstream"))
  attr(f,"INApestCouplesPathogen") <- TRUE
  attr(f,"INApestRunWhenEmpty") <- RunWhenEmpty
  class(f) <- c("INApestMetaRK4HostPathogenLocalDynamics","function")
  f
}
