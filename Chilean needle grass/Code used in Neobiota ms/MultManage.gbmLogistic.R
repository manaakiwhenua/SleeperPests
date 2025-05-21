rm(list = ls())
library(hstats)
library(graphics)
rm(list = ls())
ScriptDir = r"[C:\Users\MasonN\OneDrive - MWLR\Documents\GitHub\SleeperPests\Chilean needle grass\Code used in Neobiota ms\xgb.cv\]"
ImageDir = r"[R:\Projects\SL2116_MPI_SleeperPests\Analysis\Chilean needle grass\Updated Analyses 2022\CNG paper Images\]"
dir.create(ImageDir)
setwd(ImageDir)
source(paste0(ScriptDir,"xgb.cv.interaction.r"))
cs = 1
SpreadReductions = c(0.05,0.1,0.3,0.5,0.9,0.99,0.999)
EradicationProbs = c(0,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.75)
CommProbs = c(0,0.025,0.05,0.1,0.15,0.2,0.25,0.3,0.5,0.75)
DetectionProbs = c(0,0.025,0.05,0.075,0.1,0.15,0.2,0.25,0.3,0.4,0.5,0.75) 

X = SpreadReductions
Y = EradicationProbs


Filename = paste0(ImageDir,"MultiManagementSummaryLong_ClimateScenario_",1,".csv")
CurrentData = read.csv(Filename,as.is = T, header = T)
head(CurrentData)
CurrentData$Climate = rep("Current", times = nrow(CurrentData))

cs = 3
Filename = paste0(ImageDir,"MultiManagementSummaryLong_ClimateScenario_",3,".csv")
FutureData = read.csv(Filename,as.is = T, header = T)
FutureData$Climate = rep("Future", times = nrow(FutureData))

AllData = rbind.data.frame(CurrentData,FutureData)
AllData = AllData[AllData$DetectionProb>0,]
head(AllData)
AllData$Invasion = AllData$Invasion*100
AllData$Success = ifelse(AllData$Invasion<1,1,0)
Response = "Success"
head(AllData)
AllData$Climate = ifelse(AllData$Climate == "Current",0,1)
###make sure in same order as in data frame
Predictors = c("CommunicationProb","SpreadReduction","EradicationProb","DetectionProb","Climate")
source(paste0(ScriptDir,"wrapper.xgb.cv.logistic.R"))

path = paste0(ImageDir,"MultManageGBM_Logistic/")
dir.create(path)
Nthread = 4
MaxDepth = 3
Nfolds = 10
Nrounds = 10000
LearningRate = 0.1
CV <- xgb.cv.logistic(Data = AllData,Predictors = Predictors,Response=Response,path=path,Nrounds = Nrounds,LearningRate = LearningRate,
                Nthread = 2,MaxDepth=MaxDepth,save = TRUE)
saveRDS(CV,paste0(path,"Success.xgb.cv.out.rds"))

CVout = readRDS(paste0(path,"Success.xgb.cv.out.rds"))
write.csv(CVout$Interaction, paste0(path,"Interactions.csv"),row.names = F)
write.csv(CVout$`Predictor importance`[1], paste0(path,"Gain.csv"),row.names = F)
write.csv(CVout$`Predictor importance`[2], paste0(path,"Cover.csv"),row.names = F)
write.csv(CVout$`Predictor importance`[3], paste0(path,"Frequency.csv"),row.names = F)

###########################################
###Inspect interactions
###########################################
###load the xgb.cv output
ModelFile = paste0(path,"xgb.cv.logisitc.rds")
cv = readRDS(ModelFile)
CVtrain_x = AllData[,colnames(AllData) %in% Predictors]
Names = colnames(CVtrain_x)
Filename = paste0(path,"PredictorImportance.png")
Importance <- xgb.cv.importance.plot(cv, #ouput from xgb.cv. Be sure to use callback to save cv models
                                     Nfolds, #number of fold models used in cross-validaton
                                     Predictors= Names[Names%in% Predictors],#names of predictor variables
                                     #this ensures names in right order    
                                     #for importance function
                                     Filename)

###Produce perspective plots for interactions between
###continuous predictors
Interactions = xgb.cv.interaction(cv,AllData[,colnames(AllData) %in% Predictors],Predictors,Nfolds)

###First set up predictor data frame 
Xpred = AllData[1:50,colnames(AllData) %in% Predictors]
for(i in 1:ncol(Xpred))
  Xpred[,i] = seq(0,max(AllData[,colnames(AllData)==colnames(Xpred[i])]),length = 50)
Xpred = as.matrix(Xpred)


Var2 = "EradicationProb"
Var1 = "DetectionProb"
xgb.cv.perspective(cv,Nfolds,CVtrain_x=Xpred,Var1,Var2,path,Response = Response,ResponseName = "Probability of success",Var1Lab = "Detection", Var2Lab = "Eradication")

Var1 = "CommunicationProb"
Var2 = "DetectionProb"
xgb.cv.perspective(cv,Nfolds,CVtrain_x=Xpred,Var1,Var2,path,Response = Response,ResponseName = "Probability of success",Var1Lab = "Communication", Var2Lab = "Detection")

Var1 = "CommunicationProb"
Var2 = "EradicationProb"
xgb.cv.perspective(cv,Nfolds,CVtrain_x=Xpred,Var1,Var2,path,Response = Response,ResponseName = "Probability of success",Var2Lab = "Eradication", Var1Lab = "Communication")

Var1 = "SpreadReduction"
Var2 = "EradicationProb"
xgb.cv.perspective(cv,Nfolds,CVtrain_x=Xpred,Var1,Var2,path,Response = Response,ResponseName = "Probability of success",Var2Lab = "Eradication", Var1Lab = "Spread reduction")

