rm(list=ls())
library(ggplot2)

sen = readRDS("../Data/sensitivity_matrix_AUC.rds")
GE = readRDS("../Data/expresion_matrix.rds")
N_drugs = 1448

#Landmark
l1000_genes = readRDS("../Data/Landmark_genes.rds")

#TF activity
TF = read.table("../Data/TF_gsea_PRISM.csv",
                sep = ",",header = TRUE, row.names = 1)

# Pathway activity
PA = readRDS("../Data/PW_activity.rds")

#PW
#Read Drug Pathway results
Ridge_PW = c()
I_zeros = c()
for(i in 1:N_drugs){
  print(i)
  R_PW = readRDS(paste0("../Data/Results_DrugPathways_All_Models/Result_",as.character(i),".rds"))
  if(is.null(nrow(R_PW))){
    I_zeros = c(I_zeros,i)
  }else{
    Ridge_PW = rbind(Ridge_PW, R_PW[4,])
  }
}


data = data.frame(
  name=letters[1:5],
  value = c(ncol(GE), mean(Ridge_PW[,5]), length(l1000_genes),
             ncol(TF), ncol(PA)),
  sd=c(0, sd(Ridge_PW[,5]),0, 0, 0))

ggplot(data) +
  geom_bar(aes(x=name, y=value), colour="black", stat="identity", size = 0.3, 
           fill= c("#F2E6D6","#FCC0C5","#A0E0E4","#F582A8","#A49393"), alpha=1, 
           width = 0.5) + scale_y_log10()+
  geom_errorbar(aes(x=name, ymin=value-sd, ymax=value+sd), width=0.2,
                colour="black", alpha = 0.8, size = 0.6)+ theme_classic()




