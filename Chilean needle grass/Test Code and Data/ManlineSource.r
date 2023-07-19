#rm(list = ls())
val.path <- paste0( "N:/Projects/BaseData/RPackages/R-4.1.3")
if (dir.exists(val.path))
{
  .libPaths(val.path)
  .libPaths()
} else
{
  cat(val.path, "\r\n")
  stop( "Cannot find the lib path specified")
}

###Set root directory where testing data and scripts stored
main.dir = r"[C:\Users\MasonN\OneDrive - MWLR\Documents\GitHub\SleeperPests\Chilean needle grass\Test Code and Data]"
dir.create(main.dir,showWarnings = F)
setwd(main.dir)
memory.limit(size = 600000)
ClimexVars = c("EI","EI_Niwa_204","EI_Niwa209")
ClimateScenarios = c("Current","Future_2040","Future_2090")
EI_Prob_CurveType = "Logit"
#EI_Prob_CurveType = "SplitLinear"


###Run test source code
for(cs in 1:3)
source("Scripts/INA__BlindRiver_TestingApril2023.R")

###Run test source code
for(cs in 1:3)
source("Scripts/INA__BlindRiver_TestingApril2023v2.R")