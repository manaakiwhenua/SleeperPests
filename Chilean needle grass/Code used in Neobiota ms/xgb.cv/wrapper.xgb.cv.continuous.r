library(caret)
library(ggplot2)
library(pdp)
library(xgboost)
library(boot)
library(yardstick)
source(paste0(ScriptDir,"xgb.cv.importance.plot.R"))
source(paste0(ScriptDir,"xgb.cv.partial.r"))
source(paste0(ScriptDir,"xgb.cv.fit.boxplot.r"))
source(paste0(ScriptDir,"xgb.cv.fit.scatterplot.r"))
source(paste0(ScriptDir,"xgb.cv.interaction.r"))
source(paste0(ScriptDir,"xgb.cv.makefolds.R"))

xgb.cv.continuous = function(Data,Predictors,Response,Objective = "reg:absoluteerror",Metric = "mae",path,Nfolds = 10,Nrounds = 10000,LearningRate=0.001
                           ,Nthread = 2,MaxDepth=3,save = TRUE, Folds = NULL, Monotone = NULL,DoInteraction = TRUE)
{
CVtrain_x = as.matrix(Data[, colnames(Data) %in% Predictors])
CVtrain_y = Data[,colnames(Data) == Response]

###Convert fold vector (if supplied) to list of obsrvations in each fold
###Assumes length of fold vector = nrow(Data)
K = Nfolds
FoldList = NULL
if(is.null(Folds)==FALSE)
  {
  K = min(Nfolds,length(unique(Folds)))
  FoldList <- xgb.cv.makefolds(as.factor(Folds), K)
  }
Nfolds = K

if(is.null(Monotone)==TRUE)
  Monotone = rep(0,times = ncol(CVtrain_x))

cv <- xgb.cv(data = CVtrain_x, stratified = TRUE,label = CVtrain_y,nrounds = Nrounds, nthread = Nthread, nfold = Nfolds,folds = FoldList,monotone_constraints =Monotone,
             max_depth = MaxDepth, eta = min(50,Nrounds),objective = Objective,metric = Metric,prediction = TRUE,print_every_n = 50,learning_rate = LearningRate,
             save_models = TRUE,early_stopping_rounds = 50,callbacks = list(cb.cv.predict(save_models = TRUE)))
Nfolds = length(cv$models)
#mean(CVtrain_y)
#mean(cv$pred)
Cor = round(cor(CVtrain_y,cv$pred),digits = 3)
if(save==TRUE)
  saveRDS(cv,paste0(path,"xgb.cv.continuous.rds"))

###Print scatter plot of predicted against observed values of response
xgbm.cv.fit.scatterplot(cv$pred,Data[, colnames(Data) == Response],path)

####Use custom function to generate predictor importance bar plots
Filename = paste0(path,"PredictorImportance.png")
Names = colnames(CVtrain_x)
Filename = paste0(path,"PredictorImportance.png")
Importance <- xgb.cv.importance.plot(cv, #ouput from xgb.cv. Be sure to use callback to save cv models
                       Nfolds, #number of fold models used in cross-validaton
                       Predictors= Names[Names%in% Predictors],#names of predictor variables
                                                              #this ensures names in right order    
                                                              #for importance function
                       Filename)#location to print bar plot 



####Use custom function to generate partial dependency plots
PartialDir = paste0(path,"PartialDependencePlots/")
dir.create(PartialDir,showWarnings = FALSE)
for(var in 1:length(Predictors))
    xgbm.cv.partial(cv,Nfolds = Nfolds,na.omit(CVtrain_x),var,path = PartialDir,CVtrain_y=CVtrain_y,ResponseName = Response)

###Do interaction last as hstats changes model predictions somehow in partial plots
if(DoInteraction == TRUE)
  Interaction = xgb.cv.interaction(cv,na.omit(CVtrain_x),Predictors,Nfolds)

OutList = list()
Key = "Model"
OutList[[Key]] = cv
Key = "Correlation"
OutList[[Key]] = Cor
Key = "Predictor importance"
OutList[[Key]]= Importance
Key = "Interaction"
if(DoInteraction == TRUE)
  OutList[[Key]] = Interaction
return(OutList)
}