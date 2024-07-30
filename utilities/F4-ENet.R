library(glmnet)
library(caret)

ElasticNet = function(ytrain,Xtrain,Xtest){
  train_data = cbind(Xtrain,ytrain)
  control = trainControl(method = "repeatedcv",
                         number = 5,
                         repeats = 5,
                         verboseIter = FALSE)
  
  tune = expand.grid(alpha = seq(.05, 1, length = 10),
                     lambda = seq(.000001,0.0001,.000001))


  model = caret::train(Xtrain, ytrain,
                         method = "glmnet",
                         metric="RMSE",
                         allowParallel = TRUE,
                         tuneGrid = tune,
                         trControl = control)
  
  y_pred = predict(model,Xtest)
  return(y_pred)
  
}
