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
Filename = paste0(ImageDir,"InvasionSummary.png")
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

###Option to plot points for Marlborough historic data
points(13,56*0.67,col = 1,cex = 15,pch = 19) ###Historic N farms from Bell 2006 https://nzpps.org/_journal/index.php/nzpp/article/view/4417/4245
###67% of infested paddocks under grazing
points(18,(96-20)*0.67,col = 1,cex = 15,pch = 19) ###Subtract 20 here as 20 new incursions due to subdivisions 
dev.off()