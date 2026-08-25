source("INApestPathogen.R")
source("INApestPointPathogen.R")

check <- function(x, msg) if (!isTRUE(x)) stop(msg)
set.seed(1)

p <- INApestPathogen(Model="SIR", Beta=1, RecoveryProb=0,
                     PathogenMortalityProb=0, InitialInfected=0,
                     IntroductionProb=0)
mod <- INApestPointPathogenInteraction(p, ContactRadius=2, ContactProb=1)
pts <- data.frame(id=1:3, x=c(0,1,10), y=0,
                  pathogen_state=c("I","S","S"),
                  have_info=FALSE, detected=FALSE, managing=FALSE,
                  last_known_timestep=NA_integer_)
ct <- mod$Contact(pts, timestep=1, perm=1, context=list())
check(nrow(ct)==1L && ct$id1==1L && ct$id2==2L, "distance contact")
z <- mod$Update(pts, ct, timestep=1, perm=1, context=list())
check(z$points$pathogen_state[2]=="I", "contact transmission")
check(z$points$pathogen_state[3]=="S", "no distant transmission")

# Missing state on a recruit defaults to S.
rec <- pts[1:2,]; rec$pathogen_state[2] <- NA
z <- mod$Update(rec, data.frame(id1=integer(),id2=integer()), timestep=1, perm=1, context=list())
check(z$points$pathogen_state[2]=="S", "new recruit defaults susceptible")

# Certain progression.
p2 <- INApestPathogen(Model="SEIR", Beta=0, RecoveryProb=0,
                      ProgressionProb=1, InitialInfected=0)
m2 <- INApestPointPathogenInteraction(p2)
e <- pts[1,]; e$pathogen_state <- "E"
z <- m2$Update(e, data.frame(id1=integer(),id2=integer()), timestep=1, perm=1, context=list())
check(z$points$pathogen_state=="I", "certain progression")

# Certain recovery and waning operate on the state present at start of step.
p3 <- INApestPathogen(Model="SIR", Beta=0, RecoveryProb=1, ImmunityLossProb=1,
                      InitialInfected=0)
m3 <- INApestPointPathogenInteraction(p3)
i <- pts[1,]; i$pathogen_state <- "I"
z <- m3$Update(i, data.frame(id1=integer(),id2=integer()), timestep=1, perm=1, context=list())
check(z$points$pathogen_state=="R", "newly recovered does not wane same step")
r <- pts[1,]; r$pathogen_state <- "R"
z <- m3$Update(r, data.frame(id1=integer(),id2=integer()), timestep=1, perm=1, context=list())
check(z$points$pathogen_state=="S", "waning")

# Pathogen-associated mortality is marked then applied by the parent adapter.
p4 <- INApestPathogen(Model="SIR", Beta=0, RecoveryProb=0,
                      PathogenMortalityProb=1, InitialInfected=0)
m4 <- INApestPointPathogenInteraction(p4)
i <- pts[1:2,]; i$pathogen_state <- c("I","S")
z <- m4$Update(i, data.frame(id1=integer(),id2=integer()), timestep=1, perm=1, context=list())
check(z$points$.pathogen_death[1] && !z$points$.pathogen_death[2], "mortality marker")
a <- INApestPointPathogenApplyDeaths(z$points)
check(nrow(a)==1L && a$id==2L, "apply mortality")

# Introduction changes state without creating a point.
p5 <- INApestPathogen(Model="SIS", Beta=0, RecoveryProb=0,
                      IntroductionProb=1, IntroductionNumber=1,
                      InitialInfected=0)
m5 <- INApestPointPathogenInteraction(p5)
s <- pts[2:3,]; s$pathogen_state <- "S"
z <- m5$Update(s, data.frame(id1=integer(),id2=integer()), timestep=1, perm=1, context=list())
check(nrow(z$points)==2L && all(z$points$pathogen_state=="I"), "introduction no host creation")

cat("PASS: generic point pathogen interaction mechanism suite\n")

# Automatic initial seeding uses global counts only when no explicit state field is supplied.
p6 <- INApestPathogen(Model="SIR", Beta=0, RecoveryProb=0,
                      InitialInfected=2, InitialRecovered=1)
m6 <- INApestPointPathogenInteraction(p6)
seedpts <- data.frame(id=1:5, x=1:5, y=0)
set.seed(2)
seeded <- m6$Initialize(seedpts, perm=1)
check(sum(seeded$pathogen_state=="I")==2L && sum(seeded$pathogen_state=="R")==1L,
      "automatic initial point seeding")
explicit <- seedpts; explicit$pathogen_state <- c("I","S","S","S","S")
seeded2 <- m6$Initialize(explicit, perm=1)
check(sum(seeded2$pathogen_state=="I")==1L, "explicit initial state authoritative")
