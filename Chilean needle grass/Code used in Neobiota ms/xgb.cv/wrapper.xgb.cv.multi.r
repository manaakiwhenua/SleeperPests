###load required libraries
library(caret)
library(ggplot2)
library(pdp)
library(xgboost)
library(boot)
library(yardstick)
library(hstats)

###load required scripts
source(paste0(ScriptDir,"xgb.cv.importance.plot.R"))
source(paste0(ScriptDir,"xgb.cv.partial.r"))
source(paste0(ScriptDir,"xgb.cv.fit.boxplot.r"))
source(paste0(ScriptDir,"xgb.cv.interaction.r"))
source(paste0(ScriptDir,"xgb.cv.makefolds.R"))

xgb.cv.multi = function(Data,Predictors,Response,Objective = "multi:softprob",Metric = "mlogloss",path,Nfolds = 10,Nrounds = 10000,LearningRate = 0.001,
                        Nthread = 2,MaxDepth=3,save = TRUE,Folds = NULL, DoInteraction = TRUE)
{

###Need to provide number of classes to xgb.cv
Nclasses = length(unique(Data[,colnames(Data) == Response]))

###Get the data into the right format for xgb.cv
Data[,colnames(Data) == Response] = as.factor(Data[,colnames(Data) == Response])
Classes = levels(Data[,colnames(Data) == Response])  
CVtrain_x = as.matrix(Data[, colnames(Data) %in% Predictors])
CVtrain_y = as.numeric(Data[,colnames(Data) == Response])-1

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

###Call xgb.cv
cv <- xgb.cv(data = CVtrain_x, stratified = TRUE,label = CVtrain_y,nrounds = Nrounds, nthread = Nthread, nfold = Nfolds,folds = FoldList,
             max_depth = MaxDepth, eta = min(50,Nrounds), objective = Objective,metric = Metric,num_class = Nclasses,prediction = TRUE,print_every_n = 50,learning_rate = LearningRate,
             save_models = TRUE,early_stopping_rounds = 50,callbacks = list(cb.cv.predict(save_models = TRUE)))

Nfolds = length(cv$models)
###Save the output to disk
if(save==TRUE)
  saveRDS(cv,paste0(path,"xgb.cv.multi.rds"))

###Test accuracy of predictions
PredClass = apply(cv$pred,1,which.max)-1
Confusion = confusionMatrix(as.factor(PredClass),as.factor(CVtrain_y))

###Calculate ROC
###Currently throws a warning on first use
ROC_DF = data.frame(cv$pred,Data[,colnames(Data) == Response])
colnames(ROC_DF) = c(Classes,"Observed")
ROC = roc_auc(
  ROC_DF,
  "Observed",
  Classes
)


###Print box plots of predicted probabilities against observed occurrences for each class 
xgbm.cv.fit.boxplot.multi(cv$pred,Data[, colnames(Data) == Response],Classes,path)

###Use custom function to generate predictor importance bar plots
###and save to disk
Names = colnames(CVtrain_x)
Filename = paste0(path,"PredictorImportance.png")
Importance <- xgb.cv.importance.plot(cv, #ouput from xgb.cv. Be sure to use callback to save cv models
                       Nfolds, #number of fold models used in cross-validaton
                       Predictors= Names[Names%in% Predictors],#names of predictor variables
                                                               #this ensures names in right order    
                                                               #for importance function
                       Filename)#location to print bar plot 

###Use custom function to generate partial dependency plots
###and save to disk
PartialDir = paste0(path,"PartialDependencePlots/")
dir.create(PartialDir,showWarnings = FALSE)
Classes = levels(as.factor(Data[,colnames(Data) == Response]))
for(class in 1:length(Classes))
  for(var in 1:length(Predictors))
    xgbm.cv.partial.multiclass(cv = cv,Nfolds = Nfolds,na.omit(CVtrain_x),var=var,path = PartialDir,Prevalence = nrow(Data[Data[,colnames(Data) == Response]==Classes[class],])/nrow(Data),
                               ResponseLevel = class,Classes = Classes)


###Do interaction last as thas can be most time-consuming step
###when there a many variables
InteractionList = list()
if(DoInteraction == TRUE)
{
for(rl in 1:ncol(cv$pred))
  {
  Interaction = xgb.cv.interaction(cv,na.omit(CVtrain_x),Predictors,Nfolds,ResponseLevel = rl)
  Key = paste0("ResponseLevel_",rl)
  InteractionList[[Key]] = Interaction
  }
}
###Built list object for output
OutList = list()
Key = "Model"
OutList[[Key]] = cv
Key = "ROC"
OutList[[Key]] = ROC
Key = "ConfusionMatrix"
OutList[[Key]] = Confusion
Key = "Predictor importance"
OutList[[Key]]= Importance
Key = "Interaction"
if(DoInteraction == TRUE)
  OutList[[Key]] = InteractionList

###Return output to environment
return(c(OutList))
}