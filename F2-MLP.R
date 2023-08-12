library(keras)
library(tensorflow)

MLP = function(ytrain,Xtrain,Xtest){
  model = keras_model_sequential()
  
  # Add layers to the model
  model %>% 
    layer_dense(units = 150, activation = 'sigmoid', input_shape = ncol(Xtrain),
    kernel_regularizer = regularizer_l2(0.5)) %>% 
    
    layer_dense(units = 100, activation = 'sigmoid')%>%  
    
    layer_dense(units = 1, activation = 'linear')
  
  
  # Compile the model
  
  model %>% compile(
    loss = 'mse',
    optimizer = optimizer_adam(learning_rate = 0.0001))
  
  
  callbacks = list(callback_early_stopping(monitor = "val_loss", patience = 5, 
                                           restore_best_weights = TRUE))
  
  # Fit the model 
  model %>% fit(
    Xtrain, 
    ytrain, 
    epochs = 150, 
    batch_size = 5, 
    validation_split = 0.2,
    callbacks = callbacks)
  
  y_pred = predict(model, Xtest)
  
  return(y_pred)
  
}
