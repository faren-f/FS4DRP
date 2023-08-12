library(glmnet)
library(caret)

Ridge = function(ytrain, Xtrain, Xtest, weight = NULL){
  
  train_data = cbind(Xtrain,ytrain)
  control = trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 5,
                         verboseIter = FALSE)

  tune = expand.grid(alpha = 0,lambda = seq(0.01,5,by = 0.1))

  model = caret::train(x = Xtrain, y = ytrain,
                       method = "glmnet",
                       weights = weight,
                       metric="RMSE",
                       allowParallel = TRUE,
                       tuneGrid = tune,
                       trControl = control)
  y_pred = predict(model,Xtest)
  return(y_pred)
  
}
