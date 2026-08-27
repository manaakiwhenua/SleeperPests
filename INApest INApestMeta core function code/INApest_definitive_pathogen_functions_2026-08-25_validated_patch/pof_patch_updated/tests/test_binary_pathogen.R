source("INApestPathogen.R")

expect <- function(x,msg) if(!isTRUE(x)) stop(msg)

# Initial pathogen presence is conditional on host presence.
p <- INApestPathogen(Model="Binary", TransmissionProb=0, ClearanceProb=0, InitialPresent=c(1,1,0))
z <- p$binary_initial(c(1,0,1),3,3)
expect(identical(as.integer(z),c(1L,0L,0L)),"initial occupancy constraint failed")

# Certain directed transmission from source 1 to target 2.
C <- matrix(0,3,3); C[1,2] <- 1
p <- INApestPathogen(Model="Binary", TransmissionProb=1, ClearanceProb=0, ContactMatrix=C)
z <- p$binary_step(c(1,0,0),c(1,1,0),1,2)
expect(identical(z$PathogenPresent,c(1L,1L,0L)),"directed transmission failed")
expect(all(z$PathogenPresent <= z$Invaded),"P <= N constraint failed")

# Certain clearance removes pathogen but not host.
p <- INApestPathogen(Model="Binary", TransmissionProb=0, ClearanceProb=1)
z <- p$binary_step(c(1,0),c(1,1),1,1)
expect(identical(z$PathogenPresent,c(0L,0L)),"clearance failed")
expect(identical(z$Invaded,c(1L,1L)),"clearance altered hosts")

# Introduction changes pathogen occupancy without creating a host.
p <- INApestPathogen(Model="Binary", TransmissionProb=0, ClearanceProb=0, IntroductionProb=1)
z <- p$binary_step(c(0,0),c(1,0),1,1)
expect(identical(z$PathogenPresent,c(1L,0L)),"introduction failed")
expect(identical(z$Invaded,c(1L,0L)),"introduction created host")

# Pathogen-driven local extinction removes both local host and pathogen.
p <- INApestPathogen(Model="Binary", TransmissionProb=0, ClearanceProb=0, PathogenHostExtinctionProb=1)
z <- p$binary_step(c(1,0),c(1,1),1,1)
expect(identical(z$Invaded,c(0L,1L)),"host extinction failed")
expect(identical(z$PathogenPresent,c(0L,0L)),"pathogen persisted without host")
expect(identical(z$HostExtinction,c(1L,0L)),"host extinction event not recorded")

# Time-varying transmission.
Tp <- matrix(c(0,0, 1,0),nrow=2,ncol=2)
C <- matrix(c(0,1,0,0),2,2,byrow=TRUE)
p <- INApestPathogen(Model="Binary", TransmissionProb=Tp, ContactMatrix=C)
z1 <- p$binary_step(c(1,0),c(1,1),1,2)
z2 <- p$binary_step(c(1,0),c(1,1),2,2)
expect(z1$PathogenPresent[2]==0L,"time schedule t1 failed")
expect(z2$PathogenPresent[2]==1L,"time schedule t2 failed")

# Binary model rejects compartment-only inputs.
bad <- try(INApestPathogen(Model="Binary",TransmissionProb=0,RecoveryProb=0.2)$binary_validate(2,1),silent=TRUE)
expect(inherits(bad,"try-error"),"unused compartment parameter was not rejected")

cat("PASS: binary INApestPathogen mechanism tests\n")
