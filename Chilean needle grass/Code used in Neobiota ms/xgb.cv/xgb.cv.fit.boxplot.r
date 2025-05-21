##################################################
###Simple boxplot of predicted probabilities for 
###binary and multiclass responses
###Separate plots fitted for each level of multiclass responses
###Designed for easy inspection of xgb.cv predictions
##################################################

xgbm.cv.fit.boxplot.multi = function(pred, ###$pred from xgb.cv output
                                     CVtrain_yFactor,Classes,path)
{
  for(class in 1:length(Classes))
    {  
    Class = Classes[class]
    Y = ifelse(CVtrain_yFactor==Classes[class],1,0)
    ###Reorder Pred and Y to ensure "success" is second level of Y
    Pred = pred[order(Y),class]
    Y=Y[order(Y)]
    ClassROC = roc_auc_vec(estimate = Pred,truth = as.factor(Y),event_level = "second")
    Title = paste0(Class, " ROC = ",round(ClassROC,digits = 3))
    Filename = paste0(path,Class,"_FitBoxplot.png")
    png(Filename, height = 1600,width = 1600)
    par(mar = c(10,12,12,2), cex.main = 4,cex.lab = 3.6,cex.axis = 3.4,mgp = c(7,2,0))
    boxplot(Pred~Y, ylim = c(0,1),main = Title,
            xlab = paste0(Class," observed"),ylab = paste0(Class," fitted probability"))
    dev.off()
    }
}

xgbm.cv.fit.boxplot.logistic = function(pred,###$pred from xgb.cv output
                                        CVtrain_y,ROC,path)
  {
  Y = CVtrain_y
  Pred = pred[order(Y)]
  Y=Y[order(Y)]
  Title = paste0("ROC = ",round(ROC,digits = 3))
  Filename = paste0(path,"FitBoxplot.png")
  png(Filename, height = 1600,width = 1600)
  par(mar = c(10,12,12,2), cex.main = 4,cex.lab = 3.6,cex.axis = 3.4,mgp = c(7,2,0))
  boxplot(Pred~Y, ylim = c(0,1),main = Title,
          xlab = paste0("Observed success"),ylab = paste0("Fitted probability"))
  dev.off()
  }
