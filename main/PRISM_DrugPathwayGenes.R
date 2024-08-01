rm(list=ls())

library(caTools)
source("../utilities/F1-Ridge.R")
source("../utilities/F2-MLP.R")
source("../utilities/F3-RandomForest.R")
source("../utilities/F4-ENet.R")
source("../utilities/F5-Lasso.R")
source("../utilities/F6-Drug_Pathway_Level_genes.R")

sen = readRDS("../Data/sensitivity_matrix_AUC.rds")
GE = readRDS("../Data/expresion_matrix.rds")

Models = c("RandomForest","ElasticNet", "Lasso","Ridge","MLP")

Result = c()
N_genes = c()
for (i in 1:ncol(sen)){            # drug loop
  print(paste0("The drug number is: ", as.character(i)))
  Mean_Corr = c()
  
  drug = colnames(sen)[i]
  pathway_gene_set = Drug_Pathway_gene_set(drug, level=1)
  
  if(!isEmpty(pathway_gene_set)){
    
    I = intersect(colnames(GE), pathway_gene_set)
    GE_I = GE[,I]
    N_genes = c(N_genes, length(I))
  
    for (M in Models){             # model loop
      model = get(M)
      
      Corr = c()
      for (j in 1:100){           # repeat loop
        print(paste0("The repeat number is: ", as.character(j)))
        
        X = GE_I[!is.na(sen[,i]),]
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
    result = c(Mean_Corr, N_genes)

  }else{
    result = 0
  }
  
  Result = rbind(Result, result)
}


