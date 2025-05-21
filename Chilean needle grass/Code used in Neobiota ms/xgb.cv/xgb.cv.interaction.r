##########################################################
##########################################################
###Functions for calculating and visualising two-way interactions
###Response level argument permits handling of mulitclass responses
###This is set to 1 by default for binary and continuous responses
###Author Norman Mason with support from Michael Mayer
###https://github.com/mayer79/
###Functions call hstats() and partial_dep()
###from hstats package:
###https://cran.r-project.org/web/packages/hstats/index.html
##########################################################
##########################################################

library(hstats)
xgb.cv.interaction = function(cv,CVtrain_x,Predictors,Nfolds,ResponseLevel = 1, verbose = TRUE)
  {
  ###predict(xgboost) stacks predictions for different levels of response
  ###in a single vector so need to apply reshape option in hstats function calls
  ###for multiclass responses
  Reshape = FALSE
  if(is.null(dim(cv$pred))==FALSE)
    Reshape = TRUE
  ###Uses rbind to combine results from fold models in
  ###a single dataframe
  AllInt = as.data.frame(matrix(nrow = 0,ncol = 2))
  colnames(AllInt) = c("Int","IntVars")
  ###hstats subsamples the predictor matrix
  ###set the same seed acrsoss folds to isolate
  ###effect of variation across folds
  SetSeed = runif(n=1,min = 1,max = 10000000)
  for(fold in 1:Nfolds)
    {
    ###xgb.Booster.complete() required if loading model from disk
    Model = xgb.Booster.complete(cv$models[[fold]])
    
    ###Estimate two-way interactions (independent of main effects)
    ###all possible predictor pairs explored
    ###threeway_m = 0 exclude three-way interactions
    set.seed(SetSeed)
    ##pairwise_m = length(Predictors) to ensure all pairwise interactions covered
    s <- hstats(Model, X = as.matrix(na.omit(CVtrain_x)),reshape = Reshape,verbose = TRUE,pairwise_m = length(Predictors),threeway_m = 0)
    
    ###extract the pairwise interactions
    ###normalize = FALSE does not standardise for size on main effects
    Int = (h2_pairwise(s, normalize = FALSE))
     ###Where there are no interactions Int may be empty.
    if(is.null(Int)==FALSE)
     {
     Int$IntVars = as.character(row.names(Int))
     Int = Int[[1]]
     Int <- as.data.frame(Int[,ResponseLevel])
     colnames(Int) = "Int"
     Int$IntVars = as.character(row.names(Int))
     AllInt = rbind(AllInt,Int)
     }
    }
  
  ###Assign null value in case of no interactions
  IntOut = NULL
  
  ###Calculate mean and standard deviation across fold models
  if(nrow(AllInt)>0)
   {
   AllInt[is.na(AllInt)==TRUE] = 0
   MeanInt = aggregate(AllInt$Int, by = list(AllInt$IntVars),FUN = mean)
   SDInt = aggregate(AllInt$Int, by = list(AllInt$IntVars),FUN = sd)
   IntOut = cbind(MeanInt,SDInt[,2])
   IntOut = IntOut[order(IntOut[,2],decreasing = T),]
   colnames(IntOut) = c("Interaction","Mean","sd")
   }
  return(IntOut)
  }

xgb.cv.perspective = function(cv,Nfolds,CVtrain_x,Var1,Var2,path,Response,ResponseLab = "Response",Var1Lab = NA, Var2Lab = NA,ResponseLevel=1,theta = NA)
  {
  if(is.na(Var1Lab)==TRUE)
    Var1Lab = Var1
  if(is.na(Var2Lab)==TRUE)
    Var2Lab = Var2
  ###predict(xgboost) stacks predictions for different levels of response
  ###in a single vector so need to apply reshape option in partial_dep function calls
  Reshape = FALSE
  if(is.null(dim(cv$pred))==FALSE)
    Reshape = TRUE
  
  ###Set up a predictor grid for 3d perspective plots
  ###constrain values to min/max to avoid extrapolation
  ###may still extrapolate to combnations of predictor values
  ###that don't occur in the data
  Xpred = CVtrain_x[1:50,]
  for(i in 1:ncol(Xpred))
    Xpred[,i] = seq(min(na.omit(CVtrain_x[,i])),max(na.omit(CVtrain_x[,i])),length = 50)
  Xpred = as.matrix(Xpred) 
  Allpd = as.data.frame(matrix(ncol=3,nrow = 0))
  colnames(Allpd) = c(Var1,Var2,"y")
  ###Get partial predictions for each fold model
  for(fold in 1:Nfolds)
    {
    ###xgb.Booster.complete() required if loading model from disk
    Model = xgb.Booster.complete(cv$models[[fold]])
    Foldpd <- partial_dep(object = Model, v = c(Var1,Var2), X = Xpred,
                          grid_size = nrow(Xpred)^2,reshape = Reshape)
    Foldpd <- Foldpd$data[,c(1:2,(2+ResponseLevel))]
    colnames(Foldpd) <- colnames(Allpd)
    Allpd = rbind(Allpd,Foldpd)
    }
  ###Calculate mean predictions across fold models
  MeanIntPD = aggregate(Allpd[,3],by = list(Allpd[,1],Allpd[,2]),FUN = mean)
  
  ###Extract values for X,Y and Z used in persp() function call
  X = sort(unique(MeanIntPD[,1]))
  Y=sort(unique(MeanIntPD[,2]))
  Z = matrix(ncol = length(Y),nrow = length(X))
  for(x in 1:length(X))
    Z[x,] = MeanIntPD[MeanIntPD[,1]==X[x],3]
  
  if(is.na(theta == TRUE))
  {
  ###Attempt to improve visualisation based on predictions
  A = max(Z[1,])
  B = max(Z[nrow(Z),])
  C = max(Z[,1])
  D = max(Z[,ncol(Z)])
  if(A<=B && C <=D)
    theta = 315
  if(A<=B && C>D)
    theta = 225
  if(A>B && C<=D)
    theta = 45
  if(A>B && C>D)
    theta = 135
  }
  ###print to file
  Filename = paste0(path,Var1,".",Var2,".",Response,".","Perspective.png")
  png(Filename,width = 3200,height=3200)
  par(cex.lab = 7,cex.axis = 5, mar = c(18,20,16,12),mgp = c(14,6,0),
      cex.main = 7.6,bty = "n")
  persp(X,Y,Z,theta=theta, phi=40, r = sqrt(10), d = 3,               # viewing pars
        shade = 0.5,ticktype = "detailed", xlab = Var1Lab,ylab = Var2Lab,zlab = ResponseLab)
  dev.off()
  }

