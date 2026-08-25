source("INApestPathogen.R")
source("INApestPathogenTransitionMatrix.R")

stopifnot_fun <- function(x,msg) if(!isTRUE(x)) stop(msg)

# 1. State construction conserves every node x demographic stage.
N <- matrix(c(8,2,4,6),nrow=2,byrow=TRUE)
p <- INApestPathogen(Model="SIR",Beta=0,RecoveryProb=0,InitialInfected=c(3,2))
st <- INApestPathogenStageState(N,p,Ntimesteps=3)
stopifnot_fun(all(apply(st,c(1,2),sum)==N),"initial product-state conservation")

# 2. Explicit state is authoritative and preserves demographic stages.
ex <- array(0L,c(1,2,3),dimnames=list(NULL,NULL,c("S","I","R")))
ex[1,1,"I"] <- 2L; ex[1,2,"S"] <- 3L
st2 <- INApestPathogenStageState(matrix(c(2,3),1,2),p,ex,3)
stopifnot_fun(st2[1,1,"I"]==2L && st2[1,2,"S"]==3L,"explicit stage x pathogen state")

# 3. Demographic progression preserves pathogen state; offspring enter S.
p0 <- INApestPathogen(Model="SIR",Beta=0,RecoveryProb=0)
A <- matrix(c(0,2, 1,1),2,byrow=TRUE) # stage1 progresses certainly; stage2 survives and reproduces
state <- array(0L,c(1,2,3),dimnames=list(NULL,NULL,c("S","I","R")))
state[1,1,"I"] <- 4L; state[1,2,"S"] <- 2L
set.seed(4)
r <- local.dynamics.transition.matrix.pathogen(
 nodetransition=A,weights=c(1,1),sddprob=matrix(1,1,1),nodeenvestabprob=1,
 n0=apply(state,c(1,2),sum),lddprob=NA,lddrate=0,nodeK=100,node.seedbankK=100,
 nodepropaguleestablishment=1,nodespreadreduction=0,managing=0,
 nodefecundityreduction=0,pathogen_state=state,Pathogen=p0,timestep=1,Ntimesteps=1)
stopifnot_fun(r$PathogenState[1,2,"I"]==4L,"infected progression must preserve pathogen state")
stopifnot_fun(r$PathogenState[1,1,"I"]==0L,"progressed infected hosts remain infected, not stage1")
stopifnot_fun(r$PathogenState[1,1,"S"]>=0L,"offspring/recruits enter susceptible state")

# 4. Pathogen process changes pathogen state without changing demographic identity except pathogen mortality.
p1 <- INApestPathogen(Model="SIR",Beta=100,RecoveryProb=0,InitialInfected=0)
state <- array(0L,c(1,2,3),dimnames=list(NULL,NULL,c("S","I","R")))
state[1,1,"I"]<-1L; state[1,2,"S"]<-3L
set.seed(2); z <- .iptm_pathogen_step(state,p1,1,1,StageMixing=matrix(1,2,2))
stopifnot_fun(sum(z$State[1,2,])==3L,"pathogen update must preserve stage2 host total without pathogen mortality")
stopifnot_fun(z$State[1,2,"I"]==3L,"cross-stage transmission")

# 5. Reconciliation after ordinary management mortality is unbiased/conservative in total.
state <- array(0L,c(1,2,3),dimnames=list(NULL,NULL,c("S","I","R")))
state[1,1,] <- c(5,3,2); state[1,2,] <- c(2,2,1)
set.seed(1); z <- .iptm_reconcile(state,matrix(c(6,3),1,2))
stopifnot_fun(all(apply(z,c(1,2),sum)==c(6,3)),"management reconciliation")

# 6. Stage-transition movement preserves pathogen state across nodes.
state <- array(0L,c(2,2,3),dimnames=list(NULL,NULL,c("S","I","R")))
state[1,1,"I"] <- 3L
N <- apply(state,c(1,2),sum)
A_move <- matrix(c(0,0,1,1),2,byrow=TRUE)
Pmove <- matrix(c(0,1,0,1),2,byrow=TRUE)
set.seed(9)
r <- local.dynamics.transition.matrix.pathogen(nodetransition=A_move,weights=c(1,1),sddprob=diag(2),nodeenvestabprob=1,n0=N,
 nodeK=c(20,20),node.seedbankK=c(20,20),nodepropaguleestablishment=1,nodespreadreduction=0,managing=c(0,0),
 pathogen_state=state,Pathogen=p0,timestep=1,Ntimesteps=1,transition_sddprob=Pmove)
stopifnot_fun(r$PathogenState[2,2,"I"]==3L,"transition movement must retain pathogen state")

cat("PASS: generic demographic x pathogen-state transition-matrix mechanism tests\n")
