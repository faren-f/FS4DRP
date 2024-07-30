
library(randomForest)

RandomForest = function(ytrain,Xtrain,Xtest,ntree=150,mtry=100){

  RF = randomForest(y = ytrain,x = Xtrain, ntree=ntree,mtry=mtry)
  y_pred = predict(RF, newdata=Xtest)
  
  return(y_pred)
}