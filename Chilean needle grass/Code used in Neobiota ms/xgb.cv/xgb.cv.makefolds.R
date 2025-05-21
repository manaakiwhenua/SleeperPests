###Uses grouping vector to create k folds
###where all observations in a group are in the same fold
library(groupdata2)
xgb.cv.makefolds <- function(Folds, K) 
  {
  Fold = fold(data = as.data.frame(Folds),k=K,id_col = "Folds")
  List = list()
  Index = 1:length(Folds)
  for(i in 1:K)
    {
    Key = i
    List[[Key]] = Index[Fold$.folds == i]
    } 
  return(List)  
  }