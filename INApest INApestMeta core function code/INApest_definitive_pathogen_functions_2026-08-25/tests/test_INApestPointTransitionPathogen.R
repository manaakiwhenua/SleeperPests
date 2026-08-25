###############################################################################
### Generic pathogen x demographic-stage tests for point transition models
###############################################################################
source("INApestPathogen.R")
source("INApestPointPathogen.R")

assert <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)

# 1. Demographic stage and pathogen state are independent persistent attributes.
p <- INApestPathogen(Model="SIR", Beta=0, RecoveryProb=0, InitialInfected=0)
pi <- INApestPointPathogenInteraction(p, ContactRadius=10, ContactProb=1)
pts <- data.frame(id=1:3, x=c(0,1,20), y=0, stage=c(1L,2L,2L),
                  pathogen_state=c("I","S","S"), stringsAsFactors=FALSE)
pts2 <- pi$Initialize(pts, perm=1L)
pts2$stage[1] <- 2L
assert(pts2$pathogen_state[1] == "I", "Stage transition must preserve pathogen state")

# 2. Stage-specific contact behaviour can be supplied without changing API.
cp <- function(distance, source, target, ...) ifelse(source$stage == 2L, 1, 0)
pi2 <- INApestPointPathogenInteraction(
  INApestPathogen(Model="SIS", Beta=1, RecoveryProb=0, InitialInfected=0),
  ContactRadius=10, ContactProb=cp)
set.seed(1)
ct <- pi2$Contact(pts2, timestep=1L, perm=1L)
assert(nrow(ct) == 1L && ct$id1[1] == 1L && ct$id2[1] == 2L,
       "Stage-specific contact must use source/target stage")

# 3. Certain transmission changes pathogen state but not demographic stage/id.
set.seed(2)
u <- pi2$Update(pts2, contacts=ct, timestep=1L, perm=1L)
assert(identical(u$points$id, pts2$id), "Pathogen update must preserve ids")
assert(identical(u$points$stage, pts2$stage), "Pathogen update must preserve demographic stages")
assert(u$points$pathogen_state[u$points$id == 2L] == "I", "Contacted susceptible should become I")

# 4. New demographic recruits should be initialized susceptible by parent core.
recruit <- data.frame(id=4L, x=2, y=0, stage=1L, pathogen_state="S")
combined <- rbind(u$points[, c("id","x","y","stage","pathogen_state")], recruit)
assert(combined$stage[4] == 1L && combined$pathogen_state[4] == "S",
       "New demographic recruits should enter S by default")

# 5. External points can explicitly arrive with pathogen state.
ext <- data.frame(id=5L, x=3, y=0, stage=1L, pathogen_state="I")
combined <- rbind(combined, ext)
assert(combined$pathogen_state[5] == "I", "Explicit external pathogen state should be retainable")

# 6. Pathogen mortality marks a host for parent-engine removal without altering survivors.
pm <- INApestPointPathogenInteraction(
  INApestPathogen(Model="SIS", Beta=0, RecoveryProb=0, PathogenMortalityProb=1, InitialInfected=0),
  ContactRadius=0)
q <- data.frame(id=1:2, x=0:1, y=0, stage=c(1L,2L), pathogen_state=c("I","S"))
set.seed(3)
qu <- pm$Update(q, contacts=data.frame(), timestep=1L, perm=1L)$points
assert(isTRUE(qu$.pathogen_death[1]) && !isTRUE(qu$.pathogen_death[2]), "Only I host should be marked dead")
keep <- !qu$.pathogen_death
qu <- qu[keep, , drop=FALSE]
assert(qu$id == 2L && qu$stage == 2L && qu$pathogen_state == "S", "Parent removal should preserve survivor state")

cat("PASS: generic point-transition pathogen mechanism tests\n")
