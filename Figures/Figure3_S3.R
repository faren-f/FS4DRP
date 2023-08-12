#Run after TrainPRISM_TestTCGA scripts
rm(list=ls())

library(ggplot2)
WholeGenes_Ridge = readRDS("all_genes_ridge.rds")
Landmark_Ridge = readRDS("Landmark_Ridge.rds")
Drug_Pathways_Ridge = readRDS("Drug_Pathway_genes_ridge.rds")
TF_activity_Ridge = readRDS("TF_activity_ridge.rds")
Pathway_activity_Ridge = readRDS("Pathway_activity_ridge.rds")

WholeGenes_Lasso = readRDS("all_genes_Lasso.rds")
Landmark_Lasso = readRDS("Landmark_Lasso.rds")
Drug_Pathways_Lasso = readRDS("Drug_Pathway_genes_Lasso.rds")
TF_activity_Lasso = readRDS("TF_activity_Lasso.rds")
Pathway_activity_Lasso = readRDS("Pathway_activity_Lasso.rds")

WholeGenes_ENet = readRDS("all_genes_ENet.rds")
Landmark_ENet = readRDS("Landmark_ENet.rds")
Drug_Pathways_ENet = readRDS("Drug_Pathway_genes_ENet.rds")
TF_activity_ENet = readRDS("TF_activity_ENet.rds")
Pathway_activity_ENet = readRDS("Pathway_activity_ENet.rds")

WholeGenes_RF = readRDS("all_genes_RF.rds")
Landmark_RF = readRDS("Landmark_RF.rds")
Drug_Pathways_RF = readRDS("Drug_Pathway_genes_RF.rds")
TF_activity_RF = readRDS("TF_activity_RF.rds")
Pathway_activity_RF = readRDS("Pathway_activity_RF.rds")

WholeGenes_MLP = readRDS("all_genes_MLP.rds")
Landmark_MLP = readRDS("Landmark_MLP.rds")
Drug_Pathways_MLP = readRDS("Drug_Pathway_genes_MLP.rds")
TF_activity_MLP = readRDS("TF_activity_MLP.rds")
Pathway_activity_MLP = readRDS("Pathway_activity_MLP.rds")


#Grouped barplot for Ranksum Test
WholeGenes_Ridge_ranksum = sum(WholeGenes_Ridge$Ranksum<0.05)
Landmark_Ridge_ranksum = sum(Landmark_Ridge$Ranksum<0.05)
Drug_Pathways_Ridge_ranksum = sum(Drug_Pathways_Ridge$Ranksum<0.05)
TF_activity_Ridge_ranksum = sum(TF_activity_Ridge$Ranksum<0.05)
Pathway_activity_Ridge_ranksum = sum(Pathway_activity_Ridge$Ranksum<0.05)

WholeGenes_Lasso_ranksum = sum(WholeGenes_Lasso$Ranksum<0.05)
Landmark_Lasso_ranksum = sum(Landmark_Lasso$Ranksum<0.05)
Drug_Pathways_Lasso_ranksum = sum(Drug_Pathways_Lasso$Ranksum<0.05)
TF_activity_Lasso_ranksum = sum(TF_activity_Lasso$Ranksum<0.05)
Pathway_activity_Lasso_ranksum = sum(Pathway_activity_Lasso$Ranksum<0.05)

WholeGenes_ENet_ranksum = sum(WholeGenes_ENet$Ranksum<0.05)
Landmark_ENet_ranksum = sum(Landmark_ENet$Ranksum<0.05)
Drug_Pathways_ENet_ranksum = sum(Drug_Pathways_ENet$Ranksum<0.05)
TF_activity_ENet_ranksum = sum(TF_activity_ENet$Ranksum<0.05)
Pathway_activity_ENet_ranksum = sum(Pathway_activity_ENet$Ranksum<0.05)

WholeGenes_RF_ranksum = sum(WholeGenes_RF$Ranksum<0.05)
Landmark_RF_ranksum = sum(Landmark_RF$Ranksum<0.05)
Drug_Pathways_RF_ranksum = sum(Drug_Pathways_RF$Ranksum<0.05)
TF_activity_RF_ranksum = sum(TF_activity_RF$Ranksum<0.05)
Pathway_activity_RF_ranksum = sum(Pathway_activity_RF$Ranksum<0.05)

WholeGenes_MLP_ranksum = sum(WholeGenes_MLP$Ranksum<0.05)
Landmark_MLP_ranksum = sum(Landmark_MLP$Ranksum<0.05)
Drug_Pathways_MLP_ranksum = sum(Drug_Pathways_MLP$Ranksum<0.05)
TF_activity_MLP_ranksum = sum(TF_activity_MLP$Ranksum<0.05)
Pathway_activity_MLP_ranksum = sum(Pathway_activity_MLP$Ranksum<0.05)

ML_methods = factor(rep(c("Ridge","Lasso","ENet","RF","MLP"),5),
                       levels = c("Ridge","MLP","RF","ENet","Lasso"))
FR_methods = factor(rep(c("WG","LM","PW","PA","TF"),c(5,5,5,5,5)),
                        levels = c("WG","LM","PW","PA","TF"))
N_sig_drugs = c(WholeGenes_Ridge_ranksum,WholeGenes_Lasso_ranksum,
  WholeGenes_ENet_ranksum,WholeGenes_RF_ranksum,WholeGenes_MLP_ranksum,
  Landmark_Ridge_ranksum,Landmark_Lasso_ranksum,Landmark_ENet_ranksum,
  Landmark_RF_ranksum,Landmark_MLP_ranksum,Drug_Pathways_Ridge_ranksum,
  Drug_Pathways_Lasso_ranksum,Drug_Pathways_ENet_ranksum,Drug_Pathways_RF_ranksum,
  Drug_Pathways_MLP_ranksum,Pathway_activity_Ridge_ranksum,Pathway_activity_Lasso_ranksum,
  Pathway_activity_ENet_ranksum,Pathway_activity_RF_ranksum,Pathway_activity_MLP_ranksum,
  TF_activity_Ridge_ranksum,TF_activity_Lasso_ranksum,TF_activity_ENet_ranksum,
  TF_activity_RF_ranksum,TF_activity_MLP_ranksum)

data = data.frame(FR_methods,ML_methods,N_sig_drugs)
ggplot(data, aes(fill=FR_methods, y=N_sig_drugs, x=ML_methods)) + 
  geom_bar(position="dodge", stat="identity",width = 0.5, color = "black", size = 0.2)+
  scale_fill_manual(values=c("#F2E6D6","#A0E0E4","#FCC0C5","#A49393","#F582A8")) + 
  theme_classic()+
  scale_x_discrete(guide = guide_axis(angle = 60))


#Grouped barplot for AUC
WholeGenes_Ridge_AUC = mean(WholeGenes_Ridge$AUC)
Landmark_Ridge_AUC = mean(Landmark_Ridge$AUC)
Drug_Pathways_Ridge_AUC = mean(Drug_Pathways_Ridge$AUC)
TF_activity_Ridge_AUC = mean(TF_activity_Ridge$AUC)
Pathway_activity_Ridge_AUC = mean(Pathway_activity_Ridge$AUC)

WholeGenes_Lasso_AUC = mean(WholeGenes_Lasso$AUC)
Landmark_Lasso_AUC = mean(Landmark_Lasso$AUC)
Drug_Pathways_Lasso_AUC = mean(Drug_Pathways_Lasso$AUC)
TF_activity_Lasso_AUC = mean(TF_activity_Lasso$AUC)
Pathway_activity_Lasso_AUC = mean(Pathway_activity_Lasso$AUC)

WholeGenes_ENet_AUC = mean(WholeGenes_ENet$AUC)
Landmark_ENet_AUC = mean(Landmark_ENet$AUC)
Drug_Pathways_ENet_AUC = mean(Drug_Pathways_ENet$AUC)
TF_activity_ENet_AUC = mean(TF_activity_ENet$AUC)
Pathway_activity_ENet_AUC = mean(Pathway_activity_ENet$AUC)

WholeGenes_RF_AUC = mean(WholeGenes_RF$AUC)
Landmark_RF_AUC = mean(Landmark_RF$AUC)
Drug_Pathways_RF_AUC = mean(Drug_Pathways_RF$AUC)
TF_activity_RF_AUC = mean(TF_activity_RF$AUC)
Pathway_activity_RF_AUC = mean(Pathway_activity_RF$AUC)

WholeGenes_MLP_AUC = mean(WholeGenes_MLP$AUC)
Landmark_MLP_AUC = mean(Landmark_MLP$AUC)
Drug_Pathways_MLP_AUC = mean(Drug_Pathways_MLP$AUC)
TF_activity_MLP_AUC = mean(TF_activity_MLP$AUC)
Pathway_activity_MLP_AUC = mean(Pathway_activity_MLP$AUC)

ML_methods = factor(rep(c("Ridge","Lasso","ENet","RF","MLP"),5),
                    levels = c("Ridge","MLP","Lasso","ENet","RF"))
FR_methods = factor(rep(c("WG","LM","PW","PA","TF"),c(5,5,5,5,5)),
                    levels = c("WG","LM","PW","PA","TF"))
AUC_ave = c(WholeGenes_Ridge_AUC,WholeGenes_Lasso_AUC,
                WholeGenes_ENet_AUC,WholeGenes_RF_AUC,WholeGenes_MLP_AUC,
                Landmark_Ridge_AUC,Landmark_Lasso_AUC,Landmark_ENet_AUC,
                Landmark_RF_AUC,Landmark_MLP_AUC,Drug_Pathways_Ridge_AUC,
                Drug_Pathways_Lasso_AUC,Drug_Pathways_ENet_AUC,Drug_Pathways_RF_AUC,
                Drug_Pathways_MLP_AUC,Pathway_activity_Ridge_AUC,Pathway_activity_Lasso_AUC,
                Pathway_activity_ENet_AUC,Pathway_activity_RF_AUC,Pathway_activity_MLP_AUC,
                TF_activity_Ridge_AUC,TF_activity_Lasso_AUC,TF_activity_ENet_AUC,
                TF_activity_RF_AUC,TF_activity_MLP_AUC)

data = data.frame(FR_methods,ML_methods,AUC_ave)
ggplot(data, aes(fill=FR_methods, y=AUC_ave, x=ML_methods)) + 
  geom_bar(position="dodge", stat="identity",width = 0.5, color = "black", size = 0.2)+
  scale_fill_manual(values=c("#F2E6D6","#A0E0E4","#FCC0C5","#A49393","#F582A8")) + 
  theme_minimal()


