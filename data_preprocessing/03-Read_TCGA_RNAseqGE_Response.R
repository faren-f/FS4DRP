rm(list=ls())
library("readxl")

# Read clinical data from table (2016-Ding)
res_TCGA = readRDS("data/Res_TCGA_temp1.rds")
Response_table = read_excel("data/raw_data/bioinfo16_supplementary_tables.xlsx",
                      sheet = 3,na = "---")
colnames(Response_table) = Response_table[2,]
Response_table = Response_table[c(-1,-2,-3),]
Response_table$drug_name = tolower(Response_table$drug_name) 
Response_table = data.frame(Response_table)

Cancer_types = unique(Response_table$Cancer)
Cancer = sapply(strsplit(Cancer_types,"\\("), FUN = function(x){return(x[2])})
Cancer = sapply(strsplit(Cancer,"\\)"), FUN = function(x){return(x[1])})
S = c()
for(i in Cancer_types){
  S = c(S,sum(Response_table$Cancer == i))
}

Cancer_type = rep(Cancer,times = S)
Response_table$Cancer_type = Cancer_type
TCGA_Patients = cbind(Cancer_type,Response_table$bcr_patient_barcode,
                    Response_table$drug_name, Response_table$measure_of_response,
                    Response_table$Cancer)
## Read RNAseq data from (http://firebrowse.org) 

ACC = read.table("data/raw_data/ACC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
ACC = ACC[,1:(ncol(ACC)-1)]
BLCA = read.table("data/raw_data/BLCA.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
BRCA = read.table("data/raw_data/BRCA.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
BRCA = BRCA[,1:(ncol(BRCA)-1)]
CESC = read.table("data/raw_data/CESC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
COAD = read.table("data/raw_data/COAD.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
ESCA = read.table("data/raw_data/ESCA.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
HNSC = read.table("data/raw_data/HNSC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
KIRC = read.table("data/raw_data/KIRC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
KIRP = read.table("data/raw_data/KIRP.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
LGG = read.table("data/raw_data/LGG.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
LIHC = read.table("data/raw_data/LIHC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
LUAD = read.table("data/raw_data/LUAD.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
LUSC = read.table("data/raw_data/LUSC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
MESO = read.table("data/raw_data/MESO.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
OV = read.table("data/raw_data/OV.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
PAAD = read.table("data/raw_data/PAAD.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
PAAD = PAAD[,1:(ncol(PAAD)-1)]
PCPG = read.table("data/raw_data/PCPG.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",fill = TRUE, header=FALSE)
PRAD = read.table("data/raw_data/PRAD.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt",fill = TRUE, header=FALSE)
READ = read.table("data/raw_data/READ.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
SARC = read.table("data/raw_data/SARC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
SKCM = read.table("data/raw_data/SKCM.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
STAD = read.table("data/raw_data/STAD.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
TGCT = read.table("data/raw_data/TGCT.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
THCA = read.table("data/raw_data/THCA.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
UCEC = read.table("data/raw_data/UCEC.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)
UCS = read.table("data/raw_data/UCS.rnaseqv2_RSEM_genes_normalized.txt",fill = TRUE, header=FALSE)

Cancers = c("ACC","BLCA","BRCA","CESC","COAD","ESCA",
            "HNSC","KIRC","KIRP","LGG","LIHC","LUAD",
            "LUSC","MESO","OV","PAAD","PCPG","PRAD","READ",
            "SARC","SKCM","STAD","TGCT","THCA","UCEC","UCS")

TCGA_GE = c()
TT = c()
for (i in Cancers){
  print(i)
  Cancer_i = get(i)
  Cancer_i = Cancer_i[-2,]
  patient_all = Cancer_i[1,]
  patient_all = substr(Cancer_i[1,],1,12)
  Cancer_i[1,] = patient_all
  
  patient_with_response = Response_table[which(Response_table$Cancer_type == i),2]
  Intersect = intersect(patient_all, patient_with_response)
  
  patient_with_response = Response_table[Response_table$bcr_patient_barcode %in% Intersect,]
  Cancer_i = Cancer_i[,c(1,which(patient_all %in% Intersect))]
  dup = which(duplicated(as.character(Cancer_i[1,])))
  
  colnames(Cancer_i) = Cancer_i[1,]
  Cancer_i = Cancer_i[-1,]
  rownames(Cancer_i) = Cancer_i[,1]
  Cancer_i = Cancer_i[,c(-1)]
  Cancer_i = t(Cancer_i)
  Cancer_i_2 = apply(Cancer_i,2,as.numeric)
  
  if(length(dup)>0){
    dup = dup-1
    for (k in dup){
      rep = apply(Cancer_i_2[c(k-1,k),],2,mean)
      Cancer_i_2[k-1,] = rep
    }
    Cancer_i_2 = Cancer_i_2[-dup,]
    Cancer_i = Cancer_i[-dup,]
  }
  rownames(Cancer_i_2) = rownames(Cancer_i)
  colnames(Cancer_i_2) = colnames(Cancer_i)
  Cancer_i = Cancer_i_2
  
  col = colnames(Cancer_i)
  col = sapply(strsplit(col,"\\|"),FUN = function(x){return(x[1])})
  colnames(Cancer_i) = col
  Cancer_i = log2(Cancer_i + 1)
  
  GE = readRDS("data/expresion_matrix.rds")
  intersect_genes = intersect(colnames(GE),colnames(Cancer_i))
  
  GE = GE[,intersect_genes]
  Cancer_i = Cancer_i[,intersect_genes]
  TT = c(TT,rep(i,nrow(Cancer_i)))
  TCGA_GE = rbind(TCGA_GE,Cancer_i)
}

I_samples = intersect(rownames(TCGA_GE),rownames(res_TCGA))
TCGA_GE = TCGA_GE[I_samples,]
res_TCGA = res_TCGA[I_samples,]

saveRDS(TCGA_GE,"data/expresion_matrix_TCGA_temp.rds")
saveRDS(res_TCGA,"data/Res_TCGA_temp2.rds")


