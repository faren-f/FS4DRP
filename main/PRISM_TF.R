rm(list=ls())

library(caTools)
source("../utilities/F1-Ridge.R")
source("../utilities/F2-MLP.R")
source("../utilities/F3-RandomForest.R")
source("../utilities/F4-ENet.R")
source("../utilities/F5-Lasso.R")

sen = readRDS("../Data/sensitivity_matrix_AUC.rds")
TF = read.table("../Data/TF_gsea_PRISM.csv", sep = ",",header = TRUE, row.names = 1)

Models = c("RandomForest", "ElasticNet", "Lasso", "Ridge", "MLP")
Result = c()

for (i in 1:ncol(sen)){            # drug loop
  print(paste0("The drug number is: ", as.character(i)))
  Mean_Corr = c()
  
  for (M in Models){             # model loop
    model = get(M)
    Corr = c()
    
    for (j in 1:100){           # repeat loop
      print(paste0("The repeat number is: ", as.character(j)))
      
      X = TF[!is.na(sen[,i]),]
      y = sen[!is.na(sen[,i]),i]
      
      sample = sample.split(y, SplitRatio = .8)
      
      Xtrain = subset(X, sample == TRUE)
      Xtest  = subset(X, sample == FALSE)
      ytrain = subset(y, sample == TRUE)
      ytest  = subset(y, sample == FALSE)
      
      # Normalization
      Mu = apply(Xtrain, 2, mean)
      SD = apply(Xtrain, 2, sd)
      Xtrain = scale(Xtrain)
      
      for(t in 1:ncol(Xtest)){
        Xtest[,t] = (Xtest[,t] - Mu[t]) / SD[t]
      }
      
      Mu_y = mean(ytrain)
      SD_y = sd(ytrain)
      ytrain = scale(ytrain)
      ytest = (ytest-Mu_y)/SD_y
      
      
      # Models
      y_pred = model(ytrain = ytrain, Xtrain = Xtrain, Xtest = Xtest)
      
      # Evaluation
      corr = cor(ytest,y_pred)
      Corr = c(Corr, corr)
      
    }
    Mean_Corr = c(Mean_Corr, mean(Corr))
    
  }
  Result = rbind(Result, Mean_Corr)
}


