rm(list=ls())
library("readxl")

# Read TCGA clinical data from table (2016-Ding)
Response_table = read_excel("data/raw_data/bioinfo16_supplementary_tables.xlsx",
                      sheet = 3,na = "---")
colnames(Response_table) = Response_table[2,]
Response_table = Response_table[c(-1,-2,-3),]
Response_table$drug_name = tolower(Response_table$drug_name) 
Response_table = data.frame(Response_table)

# Drug names that are available in TCGA
drugs = unique(Response_table$drug_name)

# drugs that are common between TCGA and PRISM
sen = readRDS("data/Sen_PRISM.rds")
TCGA_PRISM_drugs = intersect(drugs,colnames(sen))

# Patients & Cancer type in TCGA
Cancer_Patient = Response_table[!duplicated(Response_table$bcr_patient_barcode),1:2]

# Drug response measures
res_measure = unique(Response_table$measure_of_response)
res_binarized = rep(0, nrow(Response_table))

for (i in 1:length(res_measure)){
    res_binarized[which(Response_table$measure_of_response == res_measure[i])] = i
}
Response_table$response_binarized = res_binarized

# Patients-Response matrix
response_mat = matrix(0,nrow(Cancer_Patient),length(TCGA_PRISM_drugs))
rownames(response_mat) = Cancer_Patient[,2]
colnames(response_mat) = TCGA_PRISM_drugs
save_name = c()
for(i in Cancer_Patient[,2]){
  for (j in TCGA_PRISM_drugs){
    I = Response_table$bcr_patient_barcode == i & Response_table$drug_name == j
    if(sum(I)==1){
      response_mat[i,j] = Response_table[I,15]
      
    }else if(sum(I)==2){
      if(Response_table[which(I)[1],15]==Response_table[which(I)[2],15]){
        response_mat[i,j] = Response_table[which(I)[1],15]
      }else{
        save_name = rbind(save_name,c(i,j))
      }
    }else{
      response_mat[i,j] = NA
    }
  }
}
response_mat["TCGA-IB-7645","gemcitabine"] = 1

#save data
saveRDS(response_mat,"data/Res_TCGA_temp1.rds")

