rm(list=ls())

res_TCGA = readRDS("Data/Res_TCGA_temp2.rds")
sen = readRDS("Data/Sen_PRISM_temp1.rds")

GE_PRISM = readRDS("Data/expresion_matrix_temp1.rds")
GE_TCGA = readRDS("Data/expresion_matrix_TCGA_temp.rds")

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
saveRDS(GE_PRISM,"Data/expresion_matrix.rds")
saveRDS(GE_TCGA,"Data/expresion_matrix_TCGA.rds")

saveRDS(sen_PRISM,"Data/Sen_PRISM.rds")
saveRDS(res_TCGA,"Data/Res_TCGA.rds")

