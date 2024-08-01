rm(list = ls())

# Library -----------------------------------------------------------------
library('rtracklayer')
# Read Data ---------------------------------------------------------------

response = read.csv("data/raw_data/secondary-screen-dose-response-curve-parameters.csv")
RNAseq = read.table("data/raw_data/CCLE_RNAseq_rsem_transcripts_tpm_20180929.txt.gz",
                    header = TRUE, check.names = FALSE)

gene_transfer = import("Data/raw_data/gencode.v19.genes.v7_model.patched_contigs.gtf.gz")
gene_transfer = data.frame(gene_transfer)

# Finding drug targets for all drugs
drug_targets = response[!duplicated(response$name),c(12,14)]
rownames(drug_targets) = 1:nrow(drug_targets)
# Pre-processing ----------------------------------------------------------
## Log normalization of genes
expr_raw = RNAseq[,c(-1,-2)]
expr = log2(expr_raw + 1)


## Remove cell lines that do not exist in response from expression
ccle_name_intersect = colnames(expr) %in% response$ccle_name
expr = expr[,ccle_name_intersect]

expr = cbind(RNAseq[,1],expr)

## Remove expressions with low mean 
mean_expr = apply(expr[,-1], 1, mean)
hist(mean_expr,100,xlim = c(0,5))
abline(v = 0.2, col ="red")
expr = expr[mean_expr > 0.2,]

rownames(expr) = 1:nrow(expr)

## Remove cell lines that do not exist in response from expression
expr1 = expr[,-1]

## Convert ccle_name to depmap_id in expression matrix
depid_name = response[,c(2,3)]
depid_name = depid_name[!duplicated(depid_name[,2]),]
depid_name = depid_name[!is.na(depid_name[,2]),]
rownames(depid_name) = depid_name[,2]
dep_id = depid_name[colnames(expr1),1] 
colnames(expr1) = dep_id
expr = cbind(expr[,1],expr1)
rm(expr1)

## Find duplicated genes and calculate the average of them  
dup_ind = which(duplicated(expr[,1]))
dup_ENS = unique(expr[dup_ind,1])

ind = c()
for (i in dup_ENS){
  ind_i = which(expr[,1] == i)
  ave_dup_ENS_i = apply(expr[ind_i,-1],2,mean)
  expr[ind_i[1],-1] = ave_dup_ENS_i
  ind = c(ind,ind_i[-1])
}
Expr = expr[-ind,]

#'@Build_expression_matrix_[sample*genes].......................................
rownames(Expr) = Expr[,1]
Expr = Expr[,-1]
Expr = t(Expr)

## gene_ids common between gene_transfer & expression matrix
gene_transfer1 = gene_transfer[,c("gene_id","gene_name")]
gene_transfer1 = gene_transfer1[!duplicated(gene_transfer1[,1]),]

intersect_gene_id = intersect(gene_transfer1$gene_id, colnames(Expr))
Expr = Expr[,intersect_gene_id]
colnames(Expr) = gene_transfer1[gene_transfer1$gene_id %in% intersect_gene_id,2]
# Removing repeatative gene-symboles 
col_Exp = colnames(Expr)
dup_gene_symbs = unique(col_Exp[duplicated(colnames(Expr))])

ind_extra = c()
for (k in dup_gene_symbs){
  ind_rep = which(colnames(Expr)==k)
  Expr[,ind_rep] = apply(Expr[,ind_rep],1,mean)
  ind_extra = c(ind_extra,ind_rep[-1])
}
Expr = Expr[,-ind_extra]
#'@Build_response_matrix_[sample*Drug]..........................................

cell_id = rownames(Expr)
drug_name = unique(response$name)
AUC = matrix(0,length(cell_id),length(drug_name))

rownames(AUC) = cell_id
colnames(AUC) = drug_name

response = response[!is.na(response$depmap_id),]


c = 0
for (i in cell_id){
  c = c+1
  print(c)
  
  for (j in drug_name){
    cell_i_drug_j = response$depmap_id == i & response$name == j
    
    if (sum(cell_i_drug_j) == 0){
      AUC[i,j] = NA

    }else if (sum(cell_i_drug_j) == 1){
             AUC[i,j] = response[cell_i_drug_j,"auc"]
             
    }else{
      AUC[i,j] = mean(response[cell_i_drug_j,"auc"])
      
      }
   }
}

# Save Data ---------------------------------------------------------------
saveRDS(AUC, file = "data/sensitivity_matrix_AUC.rds")
saveRDS(Expr, file = "data/expresion_matrix.rds")






