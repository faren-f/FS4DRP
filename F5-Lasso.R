library(glmnet)
library(caret)

Lasso = function(ytrain,Xtrain,Xtest){
  
  control = trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 5,
                         verboseIter = FALSE)
  
  tune = expand.grid(alpha = 1,lambda = seq(.000001,.001,.00001))

  model = caret::train(y= ytrain,
                       x = Xtrain,
                       method = "glmnet",
                       metric="RMSE",
                       allowParallel = TRUE,
                       tuneGrid = tune,
                       trControl = control)
  
  y_pred = predict(model,Xtest)
  model$bestTune
  plot(model$results$RMSE)
  return(y_pred)
  
}

