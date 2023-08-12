rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

sen_PRISM = readRDS("../Data/Sen_PRISM.rds")
res_TCGA = readRDS("../Data/Res_TCGA.rds")

dR_PRISM = read.table("../Dara/gsea_PRISM.csv",sep = ",",header = TRUE, row.names = 1)
dR_TCGA = read.table("../Data/gsea_TCGA.csv",sep = ",",header = TRUE, row.names = 1)

# Remove genes whose Q3 is zero
q3_genes = apply(dR_TCGA,2,quantile,prob=0.75)
if(sum(q3_genes==0)>0){
  dR_TCGA = dR_TCGA[,-which(q3_genes==0)]
  dR_PRISM = dR_PRISM[,-which(q3_genes==0)]
}

Models = c("LinearcRegresion", "RandomForest","ElasticNet", "Lasso","Ridge","MLP")

clusterExport(cl, c("dR_PRISM","dR_TCGA","sen_PRISM","res_TCGA","Models"))
clusterEvalQ(cl, c(source("F1-Ridge.R"), 
                   source("F2-MLP.R"),
                   source("F3-RandomForest.R"),
                   source("F4-ENet.R"),
                   source("F5-Lasso.R"), 
                   source("F10-Combat_Normalization.R")))

DrugLoop = function(i){
  
  Xtrain = dR_PRISM[!is.na(sen_PRISM[,i]),]
  ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
  
  Xtest = dR_TCGA[!is.na(res_TCGA[,i]),]
  ytest = res_TCGA[!is.na(res_TCGA[,i]),i]

  X_Normalization = Combat_Scale(Xtrain,Xtest)
  
  Xtrain = X_Normalization[[1]]
  Xtest = X_Normalization[[2]]
  
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

saveRDS(Result[,1], "../Data/TF_activity_RF.rds")
saveRDS(Result[,2], "../Data/TF_activity_ENet.rds")
saveRDS(Result[,3], "../Data/TF_activity_Lasso.rds")
saveRDS(Result[,4], "../Data/TF_activity_Ridge.rds")
saveRDS(Result[,5], "../Data/TF_activity_MLP.rds")


