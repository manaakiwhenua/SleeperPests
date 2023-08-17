###Set root directory where testing data and scripts stored
main.dir = r"[C:\Users\MasonN\OneDrive - MWLR\Documents\GitHub\SleeperPests\Chilean needle grass\Test Code and Data]"
dir.create(main.dir,showWarnings = F)
setwd(main.dir)
ClimexVars = c("EI","EI_Niwa_204","EI_Niwa209")
ClimateScenarios = c("Current","Future_2040","Future_2090")
EI_Prob_CurveType = "Logit"
#EI_Prob_CurveType = "SplitLinear"

###Install required packages
#install.packages("actuar")
#install.packages("magick")
#install.packages("sf")


###Run test source code
for(cs in 1:3)
source("Scripts/INA__BlindRiver_TestingApril2023.R")

###Run test source code
for(cs in 1:3)
source("Scripts/INA__BlindRiver_TestingApril2023v2.R")

###Run test source code implementing parallel processing
###This uses rapid custom function to generate heat maps/gifs 
###for visualisation of results
###Maps not intended for publication
for(cs in 1:3)
  source("Scripts/INA__BlindRiver_TestingApril2023Parallel.R")

###Run test source code implementing simulations
###without calling INAscene
for(cs in 1:3)
  source("Scripts/INA__BlindRiver_TestingApril2023ParallelNoINAscene.R")

###Run test source code illustrating effect of long distance dispersal (LDD)
for(cs in 1:3)
  source("Scripts/INA__BlindRiver_NoLDDParallel.R")

###Run test source code illustrating effect of agency-led communication vs. 
###"Passive" information spread through social networks 
for(cs in 1:3)
  source("Scripts/INA__BlindRiver_SocialNetworksParallel.R")


###Print heat maps illustrating effect long-distance dispersal 
###WARNING: Only runs in R version 4.1.3 on Manaaki Whenua spatial lab system
###Requires various spatial layers and is sensitive to changes in version of 
###some geospatial packages
for(cs in 1:3)
  source("Scripts/NoLDDHeatMaps.R")


###Print heat maps illustrating effect of agency-led communication vs. 
###"Passive" information spread through social networks 
###WARNING: Only runs in R version 4.1.3 on Manaaki Whenua spatial lab system
###Requires various spatial layers and is sensitive to changes in version of 
###some geospatial packages
for(cs in 1:3)
  source("Scripts/SocialNetworkHeatMaps.R")


