rm(list=ls())

library(parallel)
no_cores = detectCores()
cl = makeCluster(no_cores-2)

sen_PRISM = readRDS("../Data/Sen_PRISM2.rds")
res_TCGA = readRDS("../Data/Res_TCGA.rds")

GE_PRISM = readRDS("../Data/expresion_matrix_PRISM2.rds")
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
                   source("F7-Drug_Pathway_gene_set.R"),
                   source("F9-Drug_Pathway_Level_genes_eachTarget.R"),
                   source("F10-Combat_Normalization.R")))


DrugLoop = function(i){
  
  pathway_gene_set = Drug_Pathway_gene_set_eachTarget(drug = i, level=1)[["all"]]
  result = c()
  if(!isEmpty(pathway_gene_set)){
    
    I = intersect(colnames(GE_PRISM),pathway_gene_set)
    X = GE_PRISM[,I]
    X_TCGA = GE_TCGA[,I]
    N_genes = length(I)
    
    Xtrain = X[!is.na(sen_PRISM[,i]),]
    ytrain = sen_PRISM[!is.na(sen_PRISM[,i]),i]
    
    Xtest = X_TCGA[!is.na(res_TCGA[,i]),]
    ytest = res_TCGA[!is.na(res_TCGA[,i]),i]
    
    length(ytest)
    if(length(ytest)>10){
      
      X_Normalization = Combat_Scale(Xtrain,Xtest)
      
      Xtrain = X_Normalization[[1]]
      Xtest = X_Normalization[[2]]
      
      # Ytrain normalization
      ytrain = scale(ytrain)
      ytrain = ytrain[,1]
      
      # Models
      for(M in Models){
        model = get(M)
        y_pred = model(ytrain = ytrain ,Xtrain = Xtrain,Xtest = Xtest)
        corr = cor(ytest,y_pred)
        Ranksum = wilcox.test(y_pred[ytest==1], y_pred[ytest==2], alternative ="greater")$p.value
        result = rbind(result, c(corr, Ranksum, N_genes))
      }
      
    }else{
      corr = 0
      Ranksum = 1
      N_genes = 0
      
      for(rep in 1:length(Models)){
        result = rbind(result, c(corr, Ranksum, N_genes))
      }
    }
  }else{
    corr = 0
    Ranksum = 1
    N_genes = 0
    for(rep in 1:length(Models)){
      result = rbind(result, c(corr, Ranksum, N_genes))
    }
  }
  
  return(result)
}


drugs = colnames(sen_PRISM)
result = parLapply(cl, sapply(drugs, list), DrugLoop) 


Result = data.frame()
for (k in drugs){
  Result = rbind(Result, result[[k]])
}

stopCluster(cl)

saveRDS(Result[,1], "../Data/Drug_Pathway_genes_RF.rds")
saveRDS(Result[,2], "../Data/Drug_Pathway_genes_ENet.rds")
saveRDS(Result[,3], "../Data/Drug_Pathway_genes_Lasso.rds")
saveRDS(Result[,4], "../Data/Drug_Pathway_genes_Ridge.rds")
saveRDS(Result[,5], "../Data/Drug_Pathway_genes_MLP.rds")

