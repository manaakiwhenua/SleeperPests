options(stringsAsFactors = FALSE)
script_path <- function() {
  z <- grep('^--file=', commandArgs(FALSE), value=TRUE)
  if (length(z)) return(normalizePath(sub('^--file=','',z[1L]),winslash='/',mustWork=TRUE))
  if (!is.null(sys.frame(1)$ofile)) return(normalizePath(sys.frame(1)$ofile,winslash='/',mustWork=TRUE))
  stop('Could not determine script path.')
}
test_dir <- dirname(script_path()); root <- normalizePath(file.path(test_dir,'..'),winslash='/',mustWork=TRUE)
e <- new.env(parent=.GlobalEnv)
for (f in c('INApestBiocontrol.R','INApestBiocontrolPointAdapter.R','INApestPointTransitionMatrix_Biocontrol_v0.1.R','INApestVertebratePointParallel.R'))
  sys.source(file.path(root,'engine',f),envir=e)
source(file.path(root,'src','INApestProofOfFreedom.R'))
check <- function(ok,label){if(!isTRUE(ok))stop('FAIL - ',label,call.=FALSE);cat('PASS -',label,'\n')}
zero_kernel <- function(n,parents=NULL,timestep=NULL,perm=NULL)data.frame(dx=rep(0,n),dy=rep(0,n))
support <- e$INApestBiocontrolPointSupport(xmin=-1,xmax=3,ymin=-1,ymax=1,nrow=1L,ncol=1L)
init_q <- matrix(c(0L,1L),1L,2L,dimnames=list(NULL,c('juvenile','adult')))
q <- e$INApestBiocontrolAgent(Name='Q',Stages=c('juvenile','adult'),InitialState=init_q,Release=NULL,
  Transition=diag(2),Movement=NULL,AttackStage='adult',TargetStage=1L,AttackRate=0,
  RecruitStage='juvenile',OffspringPerAttack=1L)
bc <- e$INApestBiocontrol(q)
point_interaction <- function(points,timestep,...){points$x <- points$x + 0.1;points}
args <- list(ModelName='POF_BC_VP_PSOCK',Nperm=2L,Ntimesteps=2L,Nstages=2L,Weights=c(1,1),
  Transition=diag(2),InitialPoints=data.frame(x=0.5,y=0.5,stage=1L),SDDkernel=zero_kernel,
  LDDkernel=NULL,LDDrate=0,PropaguleEstablishment=1,EnvEstabProb=1,DetectionProb=0,
  ManageProb=0,MortalityProb=0,FecundityReduction=0,SpreadReduction=0,InfoRadius=0,
  InfoTransferProb=0,ExternalInfoProb=0,OngoingExternalInfo=FALSE,OngoingExternalInvasion=FALSE,
  Vertebrate=list(Interaction=list(Update=point_interaction)),OutputDir=tempdir(),SaveResults=FALSE,
  DoProgress=FALSE,Biocontrol=bc,BiocontrolPointSupport=support)
serial <- do.call(e$INApestVertebratePoint,c(args,list(Seed=20260924L)))
parallel <- do.call(e$INApestVertebratePointParallel,c(args,list(Cores=2L,Backend='psock',Seed=20260924L)))
check(is.data.frame(parallel$BiocontrolHistory) && nrow(parallel$BiocontrolHistory)>0L,
      'patched Vertebrate Point PSOCK preserves BiocontrolHistory')
check(is.data.frame(parallel$BiocontrolPointEvents),
      'patched Vertebrate Point PSOCK preserves BiocontrolPointEvents')
ord <- function(z){if(!nrow(z))return(z); z[do.call(order,z[c('perm','timestep','agent','node','stage')]),,drop=FALSE]}
a <- ord(serial$BiocontrolHistory); b <- ord(parallel$BiocontrolHistory); rownames(a)<-NULL;rownames(b)<-NULL
check(isTRUE(all.equal(a,b,check.attributes=FALSE,tolerance=0)),
      'Vertebrate Point serial/PSOCK biocontrol-history parity')
fs <- INApestBiocontrolFreedomState(parallel,'Q')
check(identical(dim(fs$Freedom),c(2L,2L)) && all(!fs$Freedom),
      'PoF extractor reads patched Vertebrate Point PSOCK output')
cat('\nOVERALL: PASS\n')
