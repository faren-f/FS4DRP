rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

sen_PRISM = readRDS("../Data/Sen_PRISM.rds")
res_TCGA = readRDS("../Data/Res_TCGA.rds")

GE_PRISM = readRDS("../Data/expresion_matrix.rds")
GE_TCGA = readRDS("../Data/expresion_matrix_TCGA.rds")

# Remove genes whose Q3 is zero
q3_genes = apply(GE_TCGA,2,quantile,prob=0.75)
if(sum(q3_genes==0)>0){
  GE_TCGA = GE_TCGA[,-which(q3_genes==0)]
  GE_PRISM = GE_PRISM[,-which(q3_genes==0)]
}

Models = c("RandomForest","ElasticNet", "Lasso","Ridge","MLP")

clusterExport(cl, c("GE_PRISM","GE_TCGA","sen_PRISM","res_TCGA","Models"))
clusterEvalQ(cl, c(source("F1-Ridge.R"), 
                   source("F2-MLP.R"),
                   source("F3-RandomForest.R"),
                   source("F4-ENet.R"),
                   source("F5-Lasso.R"), 
                   source("F10-Combat_Normalization.R")))

DrugLoop = function(i){
  
  Xtrain = GE_PRISM[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  Xtest = GE_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
  
  X_Normalization = Combat_Scale(Xtrain,Xtest)
  
  Xtrain = X_Normalization[[1]]
  Xtest = X_Normalization[[2]]
  
  source("F9-Feature_Selection_PRISM_TCGA.R")
  selected_features = c("Pathway_activity")
  Omics_List = Feature_Selection_PRISM_TCGA(selected_features, Xtrain = Xtrain, Xtest = Xtest)
  Xtrain = Omics_List[[1]]
  index = Omics_List[[2]]
  Xtest = Omics_List[[3]]
  
  # Ytrain normalization
  ytrain = scale(ytrain)
  ytrain = ytrain[,1]
  
  # Models
  result = c()
  for(M in Models){
    model = get(M)
    y_pred = model(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
    corr = cor(ytest,y_pred)
    Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
    result = rbind(result, c(corr, Ranksum))
  }
  
  return(result)
}

N_drug = ncol(sen_PRISM)
result = parLapply(cl, sapply(1:N_drug, list), DrugLoop) 

Result = data.frame()
for (k in 1:N_drug){
  Result = rbind(Result, result[[k]])
}

N_Models = length(Models)
Result_each_Model = list()
for(m in 1:N_Models){
  R = Result[m,]
  for(d in seq(N_Models,(N_Models*N_drug)-1,N_Models)){
    R = rbind(R, Result[m+d,])
  }
  Result_each_Model[m] = list(R)
}

stopCluster(cl)

saveRDS(Result[,1], "../Data/Pathway_activity_RF.rds")
saveRDS(Result[,2], "../Data/Pathway_activity_ENet.rds")
saveRDS(Result[,3], "../Data/Pathway_activity_Lasso.rds")
saveRDS(Result[,4], "../Data/Pathway_activity_Ridge.rds")
saveRDS(Result[,5], "../Data/Pathway_activity_MLP.rds")

