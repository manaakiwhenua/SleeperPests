###############################################################################
# INApest pathogen functionality - behavioural validation / report demonstrations
# Release target: corrected 25 August 2026 pathogen-capable source bundle
# Base R only.
###############################################################################
options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(prefix) {
  hit <- grep(paste0('^', prefix, '='), args, value = TRUE)
  if (!length(hit)) return(NULL)
  sub(paste0('^', prefix, '='), '', hit[1L])
}
user_root <- arg_value('--root')
user_out <- arg_value('--out')
script_dir <- tryCatch({
  ofile <- sys.frame(1)$ofile
  if (is.null(ofile)) getwd() else dirname(normalizePath(ofile, mustWork = FALSE))
}, error = function(e) getwd())
release_name <- 'INApest_definitive_pathogen_functions_2026-08-25'
find_root <- function() {
  if (!is.null(user_root) && dir.exists(user_root)) return(normalizePath(user_root))
  cand <- unique(c(file.path(getwd(), release_name), file.path(script_dir, release_name), getwd(), script_dir))
  for (x in cand) if (file.exists(file.path(x,'src','INApestPathogen.R'))) return(normalizePath(x))
  stop('Could not locate ', release_name, '. Use --root=/path/to/release.')
}
root <- find_root()
src <- file.path(root,'src')
out <- if (is.null(user_out)) file.path(getwd(), paste0('INApest_pathogen_behavioral_validation_',format(Sys.time(),'%Y%m%d_%H%M%S'))) else user_out
dir.create(out, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out,'figures'), recursive=TRUE, showWarnings=FALSE)
source(file.path(src,'INApestPathogen.R'))
source(file.path(src,'INApestMeta.r'))

set.seed(20260825)
rows <- list(); row_i <- 0L
add_row <- function(group, parameter, value, metric, estimate, se=NA_real_, expected=NA_real_, nrep=NA_integer_, note='') {
  row_i <<- row_i + 1L
  rows[[row_i]] <<- data.frame(group=group, parameter=parameter, value=as.numeric(value), metric=metric,
                               estimate=as.numeric(estimate), se=as.numeric(se), expected=as.numeric(expected),
                               nrep=as.integer(nrep), note=note, stringsAsFactors=FALSE)
}
mean_se <- function(x) c(mean=mean(x), se=sd(x)/sqrt(length(x)))

run_engine <- function(pathogen, N0, steps=20L, nrep=1000L, establish_threshold=10L) {
  ctx <- list(n_nodes=length(N0), Ntimesteps=steps)
  cum <- peak <- finalN <- numeric(nrep)
  establish <- numeric(nrep)
  ever_node2 <- numeric(nrep)
  newinf_rate_t1 <- numeric(nrep)
  for (r in seq_len(nrep)) {
    N <- as.numeric(N0)
    st <- pathogen$Engine$Initial(N, ctx)
    init_nonS <- sum(N - st[,'S'])
    total_events <- init_nonS
    peak[r] <- sum(st[,'I'])
    if (length(N)>=2L) ever_node2[r] <- as.integer(st[2,'I']>0)
    for (tt in seq_len(steps)) {
      z <- pathogen$Engine$Step(st,N,tt,ctx)
      if (tt==1L) {
        denom <- sum(st[,'S'])
        newinf_rate_t1[r] <- if (denom>0) sum(z$NewInfections)/denom else 0
      }
      st <- z$State; N <- z$N
      total_events <- total_events + sum(z$NewInfections) + sum(z$Introduced)
      peak[r] <- max(peak[r], sum(st[,'I']))
      if (length(N)>=2L) ever_node2[r] <- max(ever_node2[r], as.integer(st[2,'I']>0))
    }
    cum[r] <- total_events
    establish[r] <- as.integer(total_events >= establish_threshold)
    finalN[r] <- sum(N)
  }
  list(cumulative=cum, peakI=peak, establishment=establish, finalN=finalN,
       ever_node2=ever_node2, newinf_rate_t1=newinf_rate_t1)
}

# 1. Binary one-step transmission: observed target infection probability should
# match TransmissionProb when a single infected source has a unit directed link.
bin_p <- c(0,.05,.1,.25,.5,.75,1); nbin <- 3000L
for (p in bin_p) {
  spec <- INApestPathogen('Binary', TransmissionProb=p, InitialPresent=c(1,0),
                          ContactMatrix=matrix(c(1,1,0,1),2,2,byrow=TRUE),
                          ClearanceProb=0, IntroductionProb=0)
  y <- numeric(nbin)
  for (r in seq_len(nbin)) y[r] <- spec$binary_step(c(1,0),c(1,1),1,1)$PathogenPresent[2]
  z <- mean_se(y); add_row('binary_transmission','TransmissionProb',p,'target_infection_probability',z['mean'],z['se'],p,nbin,
                           'One infected source; one susceptible occupied target; unit directed contact.')
}

# 2. Abundance SIR: stronger beta should increase epidemic establishment and burden.
betas <- c(0,.25,.5,1,2); nr <- 1000L
for (b in betas) {
  p <- INApestPathogen('SIR',Beta=b,RecoveryProb=.25,InitialInfected=1)
  z <- run_engine(p,100,steps=20,nrep=nr,establish_threshold=10)
  a <- mean_se(z$establishment); c1 <- mean_se(z$cumulative); pk <- mean_se(z$peakI)
  add_row('beta_response','Beta',b,'establishment_probability',a['mean'],a['se'],NA,nr,'Establishment defined here as >=10 infection events including the initial infection.')
  add_row('beta_response','Beta',b,'mean_cumulative_infection_events',c1['mean'],c1['se'],NA,nr)
  add_row('beta_response','Beta',b,'mean_peak_infectious_hosts',pk['mean'],pk['se'],NA,nr)
}

# 3. Recovery probability: faster recovery should reduce epidemic burden.
recovery <- c(.05,.1,.25,.5,1)
for (q in recovery) {
  p <- INApestPathogen('SIR',Beta=1.2,RecoveryProb=q,InitialInfected=1)
  z <- run_engine(p,100,steps=20,nrep=nr,establish_threshold=10)
  c1 <- mean_se(z$cumulative); pk <- mean_se(z$peakI)
  add_row('recovery_response','RecoveryProb',q,'mean_cumulative_infection_events',c1['mean'],c1['se'],NA,nr)
  add_row('recovery_response','RecoveryProb',q,'mean_peak_infectious_hosts',pk['mean'],pk['se'],NA,nr)
}

# 4. Pathogen mortality: greater disease-associated mortality should increase host loss.
mortality <- c(0,.05,.1,.25,.5)
for (q in mortality) {
  p <- INApestPathogen('SIR',Beta=1.2,RecoveryProb=.1,PathogenMortalityProb=q,InitialInfected=5)
  z <- run_engine(p,100,steps=20,nrep=nr,establish_threshold=10)
  loss <- 100-z$finalN; m <- mean_se(loss); c1 <- mean_se(z$cumulative)
  add_row('mortality_response','PathogenMortalityProb',q,'mean_host_loss',m['mean'],m['se'],NA,nr)
  add_row('mortality_response','PathogenMortalityProb',q,'mean_cumulative_infection_events',c1['mean'],c1['se'],NA,nr,
          'Mortality also shortens infectious residence, so infection burden need not change linearly.')
}

# 5. Spatial contact: stronger source->target connectivity should increase spread.
weights <- c(0,.01,.05,.2,.5,1)
for (w in weights) {
  C <- matrix(c(1,w,0,1),2,2,byrow=TRUE)
  p <- INApestPathogen('SIS',Beta=1,RecoveryProb=.2,InitialInfected=c(5,0),ContactMatrix=C)
  z <- run_engine(p,c(50,50),steps=10,nrep=nr,establish_threshold=10)
  m <- mean_se(z$ever_node2)
  add_row('connectivity_response','source_to_target_contact',w,'probability_target_node_ever_infected',m['mean'],m['se'],NA,nr)
}

# 6. Repeated pathogen introductions: exact finite-time benchmark.
intro <- c(0,.01,.05,.1,.25,.5); nintro <- 3000L; Tintro <- 10L
for (q in intro) {
  p <- INApestPathogen('SIS',Beta=0,RecoveryProb=0,IntroductionProb=q,IntroductionNumber=1,InitialInfected=0)
  z <- run_engine(p,50,steps=Tintro,nrep=nintro,establish_threshold=1)
  obs <- mean_se(z$establishment); expct <- 1-(1-q)^Tintro
  add_row('introduction_response','IntroductionProb',q,'probability_any_introduction',obs['mean'],obs['se'],expct,nintro,
          'With beta=0 and recovery=0, this is an exact Bernoulli repeated-introduction benchmark.')
}

# 7. Waning immunity: under repeated exposure, faster immunity loss should increase reinfection burden.
waning <- c(0,.02,.05,.1,.25,.5)
for (q in waning) {
  p <- INApestPathogen('SIR',Beta=1.2,RecoveryProb=.3,ImmunityLossProb=q,
                       InitialInfected=5,IntroductionProb=.03,IntroductionNumber=1)
  z <- run_engine(p,100,steps=40,nrep=nr,establish_threshold=10)
  c1 <- mean_se(z$cumulative)
  add_row('waning_response','ImmunityLossProb',q,'mean_cumulative_infection_events',c1['mean'],c1['se'],NA,nr,
          'Repeated introductions prevent a single early stochastic fade-out from dominating the comparison.')
}

# 8. Frequency vs density dependence at the same initial prevalence but different population sizes.
# The frequency-dependent per-susceptible risk should be approximately invariant to N;
# the density-dependent risk should rise with absolute infectious abundance.
for (mode in c('frequency','density')) {
  for (N0 in c(50,200)) {
    I0 <- as.integer(.1*N0)
    p <- INApestPathogen('SIS',Beta=.5,RecoveryProb=0,InitialInfected=I0,
                         Transmission=mode,DensityScale=100)
    z <- run_engine(p,N0,steps=1,nrep=3000,establish_threshold=1)
    m <- mean_se(z$newinf_rate_t1)
    add_row('transmission_formulation','population_size',N0,paste0(mode,'_new_infection_probability_per_susceptible'),m['mean'],m['se'],NA,3000,
            'Initial prevalence is fixed at 10%. DensityScale=100 for density-dependent transmission.')
  }
}

# 9. Detection aggregation through the actual INApestMeta parent function.
# For persistent I=5, expected node detection is 1-(1-p)^5.
detect_p <- c(0,.01,.05,.1,.25,.5); ndet <- 500L
for (q in detect_p) {
  td <- tempfile('det_'); dir.create(td)
  pp <- INApestPathogen('SIR',Beta=0,RecoveryProb=0,InitialInfected=5,DetectionProb=q,DetectionTriggersInfo=FALSE)
  INApestMeta(ModelName='det_',Nperm=ndet,Ntimesteps=1,Pathogen=pp,
    DetectionProb=0,DetectionSD=0,ManageProb=0,ManageSD=0,MortalityProb=0,MortalitySD=0,
    FecundityReduction=0,SpreadReduction=0,SpreadReductionSD=0,
    InitialPopulation=10,InitBioP=NA,InvasionRisk=0,InitialInfo=0,InitInfoP=NA,
    ExternalInfoProb=0,InfoRetentionProb=1,InfoPersistenceSteps=NA,
    EnvEstabProb=0,Survival=1,K=100,PropaguleProduction=0,PropaguleEstablishment=0,
    IncursionStartPop=1,SDDprob=matrix(1,1,1),SEAM=0,LDDprob=matrix(1,1,1),LDDrate=0,
    OngoingExternalInvasion=FALSE,OngoingExternalInfo=FALSE,OutputDir=paste0(td,'/'),DoPlots=FALSE)
  x <- readRDS(file.path(td,'det_PathogenDetectedLargeOut.rds'))
  obs <- as.numeric(x[1,1,]); m <- mean_se(obs); expct <- 1-(1-q)^5
  add_row('detection_response','DetectionProb',q,'node_detection_probability_I5',m['mean'],m['se'],expct,ndet,
          'Observed through INApestMeta PathogenDetected output, with five persistent infectious hosts.')
}

results <- do.call(rbind,rows)
write.csv(results,file.path(out,'pathogen_behavioral_results.csv'),row.names=FALSE)

# Assertions are deliberately about qualitative/benchmark behaviour rather than
# exact stochastic trajectories.
get_est <- function(g,m) { z<-results[results$group==g & results$metric==m,]; z[order(z$value),] }
assertions <- data.frame(test=character(),pass=logical(),detail=character(),stringsAsFactors=FALSE)
add_assert <- function(test,pass,detail) assertions <<- rbind(assertions,data.frame(test=test,pass=isTRUE(pass),detail=detail,stringsAsFactors=FALSE))

z <- get_est('binary_transmission','target_infection_probability')
add_assert('binary transmission matches one-step probability', max(abs(z$estimate-z$expected))<.04, paste('max abs error =',round(max(abs(z$estimate-z$expected)),4)))
z <- get_est('beta_response','establishment_probability')
add_assert('higher beta increases establishment', cor(z$value,z$estimate,method='spearman')>.9 && tail(z$estimate,1)>z$estimate[1]+.5, paste('Spearman =',round(cor(z$value,z$estimate,method='spearman'),3)))
z <- get_est('recovery_response','mean_cumulative_infection_events')
add_assert('higher recovery reduces burden', cor(z$value,z$estimate,method='spearman')<-.9, paste('Spearman =',round(cor(z$value,z$estimate,method='spearman'),3)))
z <- get_est('mortality_response','mean_host_loss')
add_assert('higher pathogen mortality increases host loss', cor(z$value,z$estimate,method='spearman')>.9, paste('Spearman =',round(cor(z$value,z$estimate,method='spearman'),3)))
z <- get_est('connectivity_response','probability_target_node_ever_infected')
add_assert('stronger connectivity increases spatial spread', cor(z$value,z$estimate,method='spearman')>.9 && z$estimate[1]<.02, paste('Spearman =',round(cor(z$value,z$estimate,method='spearman'),3),'; zero-link estimate =',round(z$estimate[1],3)))
z <- get_est('introduction_response','probability_any_introduction')
add_assert('repeated introduction matches exact benchmark', max(abs(z$estimate-z$expected))<.04, paste('max abs error =',round(max(abs(z$estimate-z$expected)),4)))
z <- get_est('waning_response','mean_cumulative_infection_events')
add_assert('waning immunity increases recurrent burden', cor(z$value,z$estimate,method='spearman')>.8, paste('Spearman =',round(cor(z$value,z$estimate,method='spearman'),3)))
zf <- get_est('transmission_formulation','frequency_new_infection_probability_per_susceptible')
zd <- get_est('transmission_formulation','density_new_infection_probability_per_susceptible')
add_assert('frequency dependence is approximately size invariant', abs(diff(zf$estimate))<.03, paste('difference =',round(diff(zf$estimate),4)))
add_assert('density dependence rises with absolute infectious abundance', diff(zd$estimate)>.03, paste('difference =',round(diff(zd$estimate),4)))
z <- get_est('detection_response','node_detection_probability_I5')
add_assert('Meta pathogen detection matches 1-(1-p)^I', max(abs(z$estimate-z$expected))<.06, paste('max abs error =',round(max(abs(z$estimate-z$expected)),4)))
write.csv(assertions,file.path(out,'pathogen_behavioral_assertions.csv'),row.names=FALSE)

# Simple report-ready SVGs. These are descriptive demonstrations, not fitted models.
plot_group <- function(group,metric,xlab,ylab,file,show_expected=FALSE) {
  z <- get_est(group,metric)
  grDevices::svg(file.path(out,'figures',file),width=7,height=4.5)
  plot(z$value,z$estimate,type='b',pch=16,xlab=xlab,ylab=ylab,main='')
  if (show_expected && any(is.finite(z$expected))) lines(z$value,z$expected,lty=2)
  if (all(is.finite(z$se))) arrows(z$value,z$estimate-1.96*z$se,z$value,z$estimate+1.96*z$se,angle=90,code=3,length=.04)
  box(); grDevices::dev.off()
}
try({
  plot_group('binary_transmission','target_infection_probability','Transmission probability','Target infection probability','01_binary_transmission.svg',TRUE)
  plot_group('beta_response','establishment_probability','Beta','Establishment probability','02_beta_establishment.svg')
  plot_group('recovery_response','mean_cumulative_infection_events','Recovery probability','Mean cumulative infection events','03_recovery_burden.svg')
  plot_group('mortality_response','mean_host_loss','Pathogen mortality probability','Mean host loss','04_pathogen_mortality_host_loss.svg')
  plot_group('connectivity_response','probability_target_node_ever_infected','Source-to-target contact weight','Probability target node ever infected','05_connectivity_spread.svg')
  plot_group('introduction_response','probability_any_introduction','Introduction probability per timestep','Probability of any introduction','06_repeated_introduction.svg',TRUE)
  plot_group('waning_response','mean_cumulative_infection_events','Immunity loss probability','Mean cumulative infection events','07_waning_burden.svg')
  plot_group('detection_response','node_detection_probability_I5','Per-infectious-host detection probability','Node detection probability','08_detection_probability.svg',TRUE)
}, silent=TRUE)

summary <- c(
  'INApest pathogen behavioural validation summary',
  paste0('Date/time: ',format(Sys.time(),'%Y-%m-%d %H:%M:%S %Z')),
  paste0('Release root: ',root),
  paste0('R version: ',R.version.string),
  paste0('Platform: ',R.version$platform),
  paste0('Assertions: ',sum(assertions$pass),' PASS / ',sum(!assertions$pass),' FAIL'),
  '',
  paste0(ifelse(assertions$pass,'PASS','FAIL'),' - ',assertions$test,' (',assertions$detail,')'),
  '',
  'Interpretation boundary:',
  '- These are generic behavioural validation scenarios designed to exercise expected epidemiological directions and exact limiting benchmarks.',
  '- They are not calibration to a particular wildlife pathogen or host system.',
  '- Regression tests separately verify the function-specific implementation contracts.'
)
writeLines(summary,file.path(out,'pathogen_behavioral_summary.txt'))
writeLines(capture.output(sessionInfo()),file.path(out,'pathogen_behavioral_sessionInfo.txt'))
cat(paste(summary,collapse='\n'),'\n')
if(any(!assertions$pass)) stop('One or more behavioural assertions failed; inspect output tables before interpreting.')
