library(reactome.db)

Drug_Pathway_gene_set = function(drug){
  Drug_Pathways = readRDS("Processed_data/S26/Drug_Pathways.rds")
  conv_table = readRDS("Processed_data/S7/biomart_conversion_table.rds")
  
  path2gene = as.list(reactomePATHID2EXTID)

  Drug_Pathways_i = Drug_Pathways[[drug]]
  Drug_Pathways_i = unlist(Drug_Pathways_i)
  
  pw = c()
  pathway_gene_set_entrez = c()
  
  for(j in Drug_Pathways_i){
    entrez_j = path2gene[[j]]
    l = length(entrez_j)
    if(l<50){
      pw = c(pw,j)
      pathway_gene_set_entrez = c(pathway_gene_set_entrez,entrez_j)
    }
  }
  pathway_gene_set_entrez = unique(pathway_gene_set_entrez)
  pathway_gene_set = conv_table[which(conv_table$entrezgene_id %in% pathway_gene_set_entrez),c(3,4)]
  
  return(pathway_gene_set)
}
