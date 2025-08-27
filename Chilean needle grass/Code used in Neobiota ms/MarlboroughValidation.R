###Read in at risk farms
MARLatRisk <- sf::st_read(r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\Current\Inputs\MARLSheepBeefAtRiskCurrent.wgs84.shp]")
FarmNames <- MARLatRisk$farm_id

###Read in infestation at 1987
InputDir = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\Inputs\]"
HistoricInfestation = read.csv(paste0(InputDir,"BlindRiverInfestedFarms.txt"),header = T, as.is =T)

###Calculate number of farms infested at 1987
Ninfested1987 <- length(FarmNames[FarmNames %in% HistoricInfestation$farm_id])
####Read in latest infestation data held my Marlborough district council
Infestation2020 <- sf::st_read(r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\MARLoccurrences\Chilean_Needlegrass_Distribution.shp]")
Infestation2020 <- sf::st_transform(Infestation2020, "+proj=longlat +datum=WGS84")
Infestation2020 <- sf::st_make_valid(Infestation2020)

###Join infestation polygons to at risk farms
InfestedFarms <-sf::st_join(Infestation2020,MARLatRisk)

###Calculate number of infested farms at 2020
InfestedFarms2020 <- unique(InfestedFarms$farm_id)
Ninfested2020<- length(InfestedFarms2020[is.na(InfestedFarms2020)==F])

CurrentDir = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CurrentFutureClimate\HistoricExamplesLDDMatrixInfoSpread_v2\Current\MARL\]"
InvasionFile = paste0(CurrentDir,"DetProb_0_eradProb_0.05_SpreadReduction_0.999_InvasionLargeOut.rds")
InvasionResults = readRDS(InvasionFile)
InvasionSummary = as.data.frame(matrix(ncol = 3, nrow = 0))
colnames(InvasionSummary) = c("Realisation",   "Timestep",  "NodesInfested")
Dim = dim(InvasionResults)
for(perm in 1:Dim[3])
{
  InvasionData = InvasionResults[,,perm]
  NodesInfested = colSums(InvasionData)
  Realisation = perm 
  Timestep = 1:Dim[2]
  Results = data.frame(Realisation,Timestep,NodesInfested)
  InvasionSummary = rbind(InvasionSummary,Results)
}


Quantiles = as.data.frame(aggregate(InvasionSummary$NodesInfested, by = list(InvasionSummary$Timestep),quantile,prob = c(0.025,0.5,0.975)))
Yvals = as.data.frame(Quantiles[,2])

ImageDir = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CNG paper Images\]"
dir.create(ImageDir)
Filename = paste0(ImageDir,"Figure 2.png")
png(Filename, height = 3200, width = 3200)
par(cex.lab = 7,cex.axis = 6.4, mar = c(18,20,16,12),mgp = c(14,6,0),#mfrow = c(2,2),
    cex.main = 7.6,bty = "n")
plot(Quantiles[,1],Yvals[,1], pch = NA, ylim = c(0,max(Yvals)), xlab = "Time since incursion detected (years)",
     ylab = "Number of farms infested")
axis(side = 1,lwd = 12,at = c(0,20,40,60,80), col = "gray80")
axis(side = 2,lwd = 12,at = c(0,50,100,150,200), col = "gray80")
lines(Quantiles[,1],Yvals[,2],lwd = 12)
lines(Quantiles[,1],Yvals[,1],lwd = 12,col = 2)
lines(Quantiles[,1],Yvals[,3],lwd = 12,col = 2)

###############################
###Plot points for Marlborough historic data
##############################
###Starting infestation in year 1987
points(0,15,col = 1,cex = 15,pch = 19) 
###Historic N farms at year 2000 from Bell 2006 https://nzpps.org/_journal/index.php/nzpp/article/view/4417/4245
points(13,56*0.67,col = 1,cex = 15,pch = 19) 
###67% of infested paddocks under grazing at 2005
points(18,(96-20)*0.67,col = 1,cex = 15,pch = 19) ###Subtract 20 here as 20 new incursions due to subdivisions 
###N infested at 2020 - from Marlborough District Council data
points(33,Ninfested2020,col = 1,cex = 15,pch = 19)
dev.off()
