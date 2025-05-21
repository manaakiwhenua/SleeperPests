##################################################
###Simple boxplot of predicted probabilities for 
###continuous responses.
###Designed for easy inspection of xgb.cv predictions
##################################################
xgbm.cv.fit.scatterplot = function(pred,CVtrain_y,path)
  {
  Cor = round(cor(pred,CVtrain_y),digits = 3)
  Title = paste0("Cor = ",Cor)
  Filename = paste0(path,"FitScatterplot.png")
  png(Filename, height = 1600,width = 1600)
  par(mar = c(10,12,12,2), cex.main = 4,cex.lab = 3.6,cex.axis = 3.4,mgp = c(7,3.5,0))
  plot(pred~CVtrain_y, main = Title,
          xlab = paste0("Observed response"),ylab = paste0("Fitted response"))
  abline(0,1,col = 2,lwd = 3)
  dev.off()
  }
