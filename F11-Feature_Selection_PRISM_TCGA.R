
source("F10-PW_activies.R")
Feature_Selection_PRISM_TCGA = function(selected_features,Xtrain,Xtest){
  if (prod(selected_features == "")){
    writeLines("selected_features is empty!\nEnter your desired features")
    omics_train = c()
    omics_test = c()
    index = c()
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == "all_genes")){
    omics_train = Xtrain
    omics_test = Xtest
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == "Landmark_genes")){
    l1000_genes = readRDS("Data/Landmark_genes.rds")
    I_G = intersect(l1000_genes,colnames(Xtest))
    
    omics_train = Xtrain[,I_G]
    omics_test = Xtest[,I_G]
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
    
  }else if (prod(selected_features == "DrugPathway_genes")){
    omics_train = Pathway_genes(X = Xtrain, drug=)
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index)
    
  }else if (prod(selected_features == "Pathway_activity")){
    omics_train = Progeny_pw_act(Xtrain)
    omics_test = Progeny_pw_act(Xtest)
    
    index = rep(1,ncol(omics_train))
    Omics = list(omics_train,index,omics_test)
  }
  return(Omics)
}
