rm(list=ls())

res_TCGA = readRDS("data/Res_TCGA_temp2.rds")
sen = readRDS("data/sensitivity_matrix_AUC.rds")

GE_PRISM = readRDS("data/expresion_matrix.rds")
GE_TCGA = readRDS("data/expresion_matrix_TCGA_temp.rds")

I_GE = intersect(colnames(GE_PRISM),colnames(GE_TCGA))
GE_PRISM = GE_PRISM[,I_GE]
GE_TCGA = GE_TCGA[,I_GE]

# Intercept of good drugs in PRISM with TCGA
I_D_G = intersect(colnames(sen_G),drugs)
sen_PRISM_G = sen_G[,I_D_G]
res_TCGA_G = res_TCGA[,I_D_G]

# Intercept of all drugs in PRISM with TCGA
I_D = intersect(colnames(sen),drugs)
sen_PRISM = sen[,I_D]
res_TCGA = res_TCGA[,I_D]


# Save data
saveRDS(GE_PRISM,"data/expresion_matrix_PRISM2.rds")
saveRDS(GE_TCGA,"data/expresion_matrix_TCGA.rds")

saveRDS(sen_PRISM,"data/Sen_PRISM2.rds")
saveRDS(res_TCGA,"data/Res_TCGA.rds")

