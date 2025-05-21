################################################################
################################################################
###Functions calling partial plot function partial_dep()
###from hstats package:
###https://cran.r-project.org/web/packages/hstats/index.html
###for each fold model in a xgb.cv output
###and plotting results for each fold models as well as mean accross folds 
###Options provided for multiclass and binary/continuous responses
###discrete (binary) predictors supported
###Author Norman Mason with support from Michael Mayer
###https://github.com/mayer79/
################################################################
################################################################

library(hstats)
library(xgboost)
###Multiclass resonses
xgbm.cv.partial.multiclass = function
(cv, #output from xgb.cv
Nfolds,#number of fold models
CVtrain_x,#predictor data (as matrix)
          #can be same as data used to train model or a synthetic dataset
var,#predictor variable (integer between 1 and number of predictors)
path,#folder for saving image file
Prevalence = NA, #the proportion of observations for the class in question
                #this is a proxy for "zero effect" of the predictor
ResponseLevel = 1,
Classes
)
{
###predict(xgboost) stacks predictions for different levels of response
###in a single vector so need to apply reshape option in partial_dep function calls
Reshape = FALSE
if(is.null(dim(cv$pred))==FALSE)
  Reshape = TRUE
Predictors = colnames(CVtrain_x)

###Loop through the fold models to obtain partial dependence estimates
PartialPreds = as.data.frame(matrix(nrow = 0,ncol = 3))
colnames(PartialPreds) = c("x","y","Fold")
for(fold in 1:Nfolds)
  {
  Model = cv$models[[fold]]
  ###call partial_dep() function
  Partial = partial_dep(object = Model, v = Predictors[var], X = CVtrain_x,reshape = Reshape)
  Partial = Partial$data[,c(1,(ResponseLevel+1))]
  Partial$Fold = fold
  colnames(Partial) = colnames(PartialPreds)
  PartialPreds = rbind(PartialPreds,Partial)
  }
###Plot partial depency for each fold model and mean across fold models
PartialPredsMean = aggregate(PartialPreds[,2],by = list(PartialPreds[,1]),FUN = mean)
Title = paste0(Classes[ResponseLevel],"  ",Predictors[var])
Filename = paste0(path,Classes[ResponseLevel],".",Predictors[var],".partial.png")
png(Filename,width = 1800, height = 1600)
par(mfcol = c(1,1), cex.main = 4, cex.lab = 3, cex.axis = 3, mar = c(10,8,8,2), 
    mgp = c(4,1,0),font.lab = 2,oma = c(0, 0, 9, 0))
plot(PartialPreds[,1],PartialPreds[,2], ylim = c(0,max(PartialPreds[,2])),pch = NA,
     xlab = Predictors[var], ylab = paste0(Classes[ResponseLevel]),main = Title)
for(fold in 1:Nfolds)
  lines(PartialPreds[PartialPreds$Fold == fold,1],PartialPreds[PartialPreds$Fold == fold,2],lty=3,col = fold,lwd = 3)
lines(PartialPredsMean[,1],PartialPredsMean[,2],lty=1,col = 1,lwd = 10)
###Prevalence is proxy for neutral partial effect of predictors
###Acts as a check the partial effects are estimated correctly
if(is.na(Prevalence) == T)
  abline(h=1/length(Classes),lwd = 8, lty = 2,col = 2)
if(is.na(Prevalence) == F)
  abline(h=Prevalence,lwd = 8, lty = 2,col = 2)
dev.off()
}

###Binary and continuous responses
xgbm.cv.partial = function(cv,Nfolds,CVtrain_x,var,path,CVtrain_y,ResponseName = "yhat")
  {
  MeanY = mean(CVtrain_y)
  Predictors = colnames(CVtrain_x)
  PartialPreds = as.data.frame(matrix(nrow = 0,ncol = 3))
  colnames(PartialPreds) = c("x","y","Fold")
  for(fold in 1:Nfolds)
    {
    Model = cv$models[[fold]]
    ###call partialPlot to get partial effect predictions
    Partial = partial_dep(object = Model, v = Predictors[var], X = CVtrain_x)
    Partial = Partial$data
    Partial$Fold = fold
    colnames(Partial) = colnames(PartialPreds)
    PartialPreds = rbind(PartialPreds,Partial)
    }
  ###Plot partial depency for each fold model and mean across fold models
  PartialPredsMean = aggregate(PartialPreds[,2],by = list(PartialPreds[,1]),FUN = mean)
  Title = Predictors[var]
  Filename = paste0(path,ResponseName,".",Predictors[var],".partial.png")
  png(Filename,width = 1800, height = 1600)
  par(mfcol = c(1,1), cex.main = 4, cex.lab = 3, cex.axis = 3, mar = c(10,8,8,2), 
      mgp = c(4,1,0),font.lab = 2,oma = c(0, 0, 9, 0))
  plot(PartialPreds[,1],PartialPreds[,2], ylim = c(min(CVtrain_y),max(CVtrain_y)),pch = NA,
       xlab = Predictors[var], ylab = ResponseName,main = Title)
  for(fold in 1:Nfolds)
    lines(PartialPreds[PartialPreds$Fold == fold,1],PartialPreds[PartialPreds$Fold == fold,2],lty=3,col = fold,lwd = 3)
  lines(PartialPredsMean[,1],PartialPredsMean[,2],lty=1,col = 1,lwd = 10)
  ###Mean of response is proxy for neutral patrial effect of predictors
  ###Acts as a check the partial effects are estimated correctly
  abline(h=MeanY,lwd = 8, lty = 2,col = 2)
  dev.off()
  }