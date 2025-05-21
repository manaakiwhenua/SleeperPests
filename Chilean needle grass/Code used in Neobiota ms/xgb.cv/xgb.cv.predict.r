###Return predictions on to new data for each fold model
###Predictions are not summarised across fold models to allow flexibility 
###In expressing results
###Predictors must be identical to predictors names used to fit models in xgb.cv
xgb.cv.predict = function(cv,PredData, Predictors = Predictors,Nfolds)
{
PredData = as.matrix(PredData[,colnames(PredData) %in% Predictors])
Preds = vector(length = 0)
Fold = vector(length = 0)
for(fold in 1:Nfolds)
 {
 Model = xgb.Booster.complete(cv$models[[fold]])
 Preds = c(Preds,predict(Model, newdata = PredData))
 Fold = c(Fold, rep(fold, times = nrow(PredData)))
 }
return(cbind(Fold,Preds))  
}

