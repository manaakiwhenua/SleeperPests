InfestationImpactSummary = function(
      Nyears,
      InvasionProbResults,
      ImageDir,
      FarmPolygons,
      TitleStem,
      FileNameStem)
{
###Area of infested farms
Filename = paste0(ImageDir,"/",FileNameStem,"_AreaAffected.png")
png(Filename, width = 1600,height = 1400, res = 300)
par(cex.axis = 0.8,cex.lab = 1, cex.main = 1.2,mar = c(5,5,5,1),mgp = c(3,1,0))
AreaAffected = colSums(sweep(InvasionProbResults,1,FarmPolygons$"size_ha","*"))/1000
Y = AreaAffected
plot(1:Nyears,Y,
xlab = "Time (years)", ylab = "Area affected (ha x 1000)", main = paste0(TitleStem,"\nArea affected"),
las = 1,pch = NA)
lines(1:Nyears,Y,lwd = 3)
dev.off()

###Area on infested farms as a proportion of total for all climatically suitable farms
Filename = paste0(ImageDir,"/",FileNameStem,"_AreaAffectedPtotal.png")
png(Filename, width = 1600,height = 1400, res = 300)
par(cex.axis = 0.8,cex.lab = 1, cex.main = 1.2,mar = c(5,5,5,1),mgp = c(3,1,0))
Y = AreaAffected/(sum(FarmPolygons$"size_ha")/1000)
plot(1:Nyears,Y,
xlab = "Time (years)", ylab = "Area affected (proportion of total)", main = paste0(TitleStem,"\nArea affected (proportion of total)"),
las = 1,pch = NA, ylim = c(0,1))
lines(1:Nyears,Y,lwd = 3)
dev.off()

###Total stock units on infested farms
Filename = paste0(ImageDir,"/",FileNameStem,"_StockUnitsAffected.png")
png(Filename, width = 1600,height = 1400, res = 300)
par(cex.axis = 0.8,cex.lab = 1, cex.main = 1.2,mar = c(5,5,5,1),mgp = c(3,1,0))
FarmPolygons$StockUnits = FarmPolygons$"bef_nos"*5.5 +   FarmPolygons$"shp_nos"*0.9 
StockUnitsAffected = colSums(sweep(InvasionProbResults,1,FarmPolygons$StockUnits,"*"))/1000
Y = StockUnitsAffected
plot(1:Nyears,Y,
xlab = "Time (years)", ylab = "Stock units affected (ewe equivalents x 1000)", main = paste0(TitleStem,"\nStock units affected"),
las = 1,pch = NA)
lines(1:Nyears,Y,lwd = 3)
dev.off()


###Total stock units on infested farms as a proportion of total for all climatically suitable farms
Filename = paste0(ImageDir,"/",FileNameStem,"_StockUnitsAffectedPtotal.png")
png(Filename, width = 1600,height = 1400, res = 300)
par(cex.axis = 0.8,cex.lab = 1, cex.main = 1.2,mar = c(5,5,5,1),mgp = c(3,1,0))
FarmPolygons$StockUnits = FarmPolygons$"bef_nos"*5.5 +   FarmPolygons$"shp_nos"*0.9 
Y = StockUnitsAffected/(sum(FarmPolygons$StockUnits)/1000)
plot(1:Nyears,Y,
xlab = "Time (years)", ylab = "Stock units affected (proportion of total)", main = paste0(TitleStem,"\nStock units affected (proportion of total)"),
las = 1,pch = NA, ylim = c(0,1))
lines(1:Nyears,Y,lwd = 3)
dev.off()


###Reduction in stock units on infested farms
Filename = paste0(ImageDir,"/",FileNameStem,"_StockUnitsReduction.png")
png(Filename, width = 1600,height = 1400, res = 300)
par(cex.axis = 0.8,cex.lab = 1, cex.main = 1.2,mar = c(5,5,5,1),mgp = c(3,1,0))
FarmPolygons$StockUnits = FarmPolygons$"bef_nos"*5.5 +   FarmPolygons$"shp_nos"*0.9 
StockUnitsReduction = 0.25*colSums(sweep(InvasionProbResults,1,FarmPolygons$StockUnits,"*"))/1000
Y = StockUnitsReduction
plot(1:Nyears,Y,
xlab = "Time (years)", ylab = "Stock units reduction (ewe equivalents x 1000)", main = paste0(TitleStem,"\nStock units reduction"),
las = 1,pch = NA)
lines(1:Nyears,Y,lwd = 3)
dev.off()

###Reduction in stock units on infested farms as a proportion of total for all climatically suitable farms
Filename = paste0(ImageDir,"/",FileNameStem,"_StockUnitsReductionPtotal.png")
png(Filename, width = 1600,height = 1400, res = 300)
par(cex.axis = 0.8,cex.lab = 1, cex.main = 1.2,mar = c(5,5,5,1),mgp = c(3,1,0))
FarmPolygons$StockUnits = FarmPolygons$"bef_nos"*5.5 +   FarmPolygons$"shp_nos"*0.9 
Y = StockUnitsReduction/(sum(FarmPolygons$StockUnits)/1000)
plot(1:Nyears,Y,
xlab = "Time (years)", ylab = "Stock units reduction (proportion of total)", main = paste0(TitleStem,"\nStock units reduction (proportion of total)"),
las = 1,pch = NA, ylim = c(0,0.25))
lines(1:Nyears,Y,lwd = 3)
dev.off()

}





