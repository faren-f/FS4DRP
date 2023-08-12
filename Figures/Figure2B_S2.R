rm(list=ls())

library(ggplot2)

sen = readRDS("../Data/sensitivity_matrix_AUC.rds")
N_drugs = 1448

#Read Landmark results
RF_Landmark = c()
ENet_Landmark = c()
Lasso_Landmark = c()
Ridge_Landmark = c()
MLP_Landmark = c()

for(i in 1:N_drugs){
  print(i)
  R_Landmark = readRDS(paste0("../Data/Results_Landmark_All_Models/Result_",as.character(i),".rds"))
  Ridge_Landmark = rbind(Ridge_Landmark, R_Landmark[4,])
  MLP_Landmark = rbind(MLP_Landmark, R_Landmark[5,])
  Lasso_Landmark = rbind(Lasso_Landmark, R_Landmark[3,])
  ENet_Landmark = rbind(ENet_Landmark, R_Landmark[2,])
  RF_Landmark = rbind(RF_Landmark, R_Landmark[1,])
}

#Read Whole Genes results
RF_WG = c()
ENet_WG = c()
Lasso_WG = c()
Ridge_WG = c()
MLP_WG = c()

for(i in 1:N_drugs){
  print(i)
  R_WG = readRDS(paste0("Processed_from_SLURM/Results_WholeGenes_All_Models/Result_",as.character(i),".rds"))
  Ridge_WG = rbind(Ridge_WG, R_WG[4,])
  MLP_WG = rbind(MLP_WG, R_WG[5,])
  Lasso_WG = rbind(Lasso_WG, R_WG[3,])
  ENet_WG = rbind(ENet_WG, R_WG[2,])
  RF_WG = rbind(RF_WG, R_WG[1,])
}

#Read Drug Pathway results
RF_PW = c()
ENet_PW = c()
Lasso_PW = c()
Ridge_PW = c()
MLP_PW = c()
I_zeros = c()
for(i in 1:N_drugs){
  print(i)
  R_PW = readRDS(paste0("Processed_from_SLURM/Results_DrugPathways_All_Models/Result_",as.character(i),".rds"))
  if(is.null(nrow(R_PW))){
    I_zeros = c(I_zeros,i)
  }else{
    Ridge_PW = rbind(Ridge_PW, R_PW[4,])
    MLP_PW = rbind(MLP_PW, R_PW[5,])
    Lasso_PW = rbind(Lasso_PW, R_PW[3,])
    ENet_PW = rbind(ENet_PW, R_PW[2,])
    RF_PW = rbind(RF_PW, R_PW[1,])
  }
}

#Read Progeny Pathway activities results
RF_PA = c()
ENet_PA = c()
Lasso_PA = c()
Ridge_PA = c()
MLP_PA = c()
I = c()
for(i in 1:N_drugs){
  print(i)
  if(file.exists(paste0("Processed_from_SLURM/Results_Progeny_Pathway_Activities_All_Models/Result_",as.character(i),".rds"))){
    R_PA = readRDS(paste0("Processed_from_SLURM/Results_Progeny_Pathway_Activities_All_Models/Result_",as.character(i),".rds"))
    Ridge_PA = rbind(Ridge_PA, R_PA[4,])
    MLP_PA = rbind(MLP_PA, R_PA[5,])
    Lasso_PA = rbind(Lasso_PA, R_PA[3,])
    ENet_PA = rbind(ENet_PA, R_PA[2,])
    RF_PA = rbind(RF_PA, R_PA[1,])
  }else{
    I = c(I,i)
  }
}


#Read decoupleR Transcription factor activities results
RF_TF = c()
ENet_TF = c()
Lasso_TF = c()
Ridge_TF = c()
MLP_TF = c()

for(i in 1:N_drugs){
  print(i)
  R_TF = readRDS(paste0("Processed_from_SLURM/Results_decoupleR_TF_Activities_All_Models/Result_",as.character(i),".rds"))
  Ridge_TF = rbind(Ridge_TF, R_TF[4,])
  MLP_TF = rbind(MLP_TF, R_TF[5,])
  Lasso_TF = rbind(Lasso_TF, R_TF[3,])
  ENet_TF = rbind(ENet_TF, R_TF[2,])
  RF_TF = rbind(RF_TF, R_TF[1,])
}


# 1) Ridge
boxplot(cbind(Ridge_PW[,1], Ridge_Landmark[,1],
              Ridge_WG[,1], Ridge_TF[,1], Ridge_PA[,1]), 
        names=NA, cex=.5, main="MeanCorr", 
        col = c("#FCC0C5","#A0E0E4","#F2E6D6","#F582A8","#A49393"), 
        ylim = c(-0.3,0.7))

# 2) MLP
boxplot(cbind(MLP_PW[,1], MLP_Landmark[,1], 
              MLP_WG[,1], MLP_TF[,1], 
              MLP_PA[,1]), names=NA, cex=.5, 
        main="MeanCorr", col = c("#FCC0C5","#A0E0E4", "#F2E6D6", "#F582A8", "#A49393"), 
        ylim = c(-0.3,0.7))

# 3) Lasso
boxplot(cbind(Lasso_WG[,1], Lasso_PW[,1], Lasso_PA[,1],
              Lasso_Landmark[,1],  Lasso_TF[,1]), names=NA, cex=.5, 
        main="MeanCorr", col = c("#F2E6D6","#FCC0C5","#A49393", "#A0E0E4", "#F582A8"), 
        ylim = c(-0.3,0.7))


# 4) ENet
boxplot(cbind(ENet_WG[,1], ENet_PW[,1], 
              ENet_Landmark[,1], ENet_PA[,1], 
              ENet_TF[,1]), names=NA, cex=.5, 
        main="MeanCorr", col = c("#F2E6D6","#FCC0C5", "#A0E0E4","#A49393","#F582A8"), 
        ylim = c(-0.3,0.7))


# 5) Random Forest
boxplot(cbind(RF_PW[,1], RF_WG[,1], 
              RF_Landmark[,1], RF_TF[,1], RF_PA[,1]), names=NA, cex=.5, 
        main="MeanCorr", col = c("#FCC0C5","#F2E6D6","#A0E0E4","#F582A8","#A49393"), 
        ylim = c(-0.3,0.7))


