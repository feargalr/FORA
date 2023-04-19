
FORA <- function(gene_list_to_test,vector_all_detected_genes,filename="FORA_results",write_to_file=FALSE){
  #When creating objects for the GO look up I made one for detected GO terms/genes anone for all.
  #This is probably a somewhat lazy work around...

  EG2GO_detected = EG2GO[EG2GO$ensembl_gene_id %in% vector_all_detected_genes,]
  vector_all_detected_genes = unique(EG2GO_detected$ensembl_gene_id)

  # Making look up tables for Ensembl to GO term and vice versa.
  ens2GO <- by(EG2GO$go_id,
               EG2GO$ensembl_gene_id,
               function(x) as.character(x))

  Go2ens_all <- by(EG2GO$ensembl_gene_id,
                   EG2GO$go_id,
                   function(x) as.character(x))


  Go2ENS_detected <- by(EG2GO_detected$ensembl_gene_id,
                        EG2GO_detected$go_id,
                        function(x) as.character(x))

  Go2ENS_detected = Go2ENS_detected[lapply(Go2ENS_detected,length)>=4]
  Go2ENS_detected = Go2ENS_detected[!names(Go2ENS_detected) %in% c(CCterms,MFterms)]


  enr_results = list()
  enr_results_reactome = list()

  genes_in_go_term = list()
  genes_in_reactome = list()



  gene_list = gene_list_to_test

  gene_list_entrez = ens2entrez[gene_list,"entrezgene_id"]
  gene_list_entrez = gene_list_entrez[!is.na(gene_list_entrez)]

  detected_entrez = ens2entrez[vector_all_detected_genes,"entrezgene_id"]
  detected_entrez = detected_entrez[!is.na(detected_entrez)]

  for (pathway in names(Reactome2gene)){
    genes1 = Reactome2gene[[pathway]]
    genes2 = genes1[genes1 %in% detected_entrez]
    Reactome2gene[[pathway]] = genes2

  }
  Reactome2gene = Reactome2gene[lapply(Reactome2gene, length)>=4]


  gene_list_entrez_hsa = ens2entrez_nona[gene_list,"entrezgene_id"]
  gene_list_entrez_hsa = gene_list_entrez_hsa[!is.na(gene_list_entrez_hsa)]

  gene_list_entrez_hsa_detecting = ens2entrez_nona[vector_all_detected_genes,"entrezgene_id"]
  gene_list_entrez_hsa_detecting = gene_list_entrez_hsa_detecting[!is.na(gene_list_entrez_hsa_detecting)]


  genes_in_go_term = list()


  #### Over representation analysis for GO Terms ####
  hypergeo_res = c()
  for (goterm in names(Go2ENS_detected)){
    go_genes = Go2ENS_detected[goterm][[1]]


    genes_in_go_term[[goterm]] = gene_list[gene_list %in% go_genes]

    n1 = length(gene_list[gene_list %in% go_genes]) #Number of genes overlapping

    n2 = length(go_genes) #Number of genes in go term

    n3 = length(ens2GO) - n2 # Number of genes not in go term

    n4 = length(gene_list) - n1 # Number of genes correlated that aren't in go term

    hypergeo_res = c(hypergeo_res,phyper(n1-1,n2,n3,n4,lower.tail = FALSE))

  }
  names(hypergeo_res) = names(Go2ENS_detected)
  adjusted_p_GO = p.adjust(hypergeo_res,method="fdr")
  #adjusted_p_GO = adjusted_p_GO[adjusted_p_GO<0.1]
  enr_results = adjusted_p_GO



  #### Over representation analysis for Reactome Pathways ####
  hypergeo_res = c()
  genes_in_reactome = list()


  for (pathway in names(Reactome2gene)){
    go_genes = Reactome2gene[pathway][[1]]


    genes_in_reactome[[pathway]] = gene_list_entrez[gene_list_entrez %in% go_genes]

    n1 = length(gene_list_entrez[gene_list_entrez %in% go_genes]) #Number of genes overlapping

    n2 = length(go_genes) #Number of genes in go term

    n3 = length(Gene2Reactome) - n2 # Number of genes not in go term

    n4 = length(gene_list_entrez) - n1 # Number of genes correlated that aren't in go term

    hypergeo_res = c(hypergeo_res,phyper(n1-1,n2,n3,n4,lower.tail = FALSE))

  }




  names(hypergeo_res) = names(Reactome2gene)
  adjusted_p_reactome = p.adjust(hypergeo_res,method="fdr")
  #adjusted_p_reactome = adjusted_p_reactome[adjusted_p_reactome<0.1]
  enr_results_reactome = adjusted_p_reactome



  ########## BTMs ##########
  my_btms = btms
  hypergeo_res = c()
  genes_in_btm = list()
  symbol_list = ens2symbol[gene_list,2]

  for (pathway in names(btms)){
    btm_genes = btms[pathway][[1]]


    my_btms[[pathway]] = symbol_list[symbol_list %in% btm_genes]

    n1 = length(symbol_list[symbol_list %in% btm_genes]) #Number of genes overlapping

    n2 = length(btm_genes) #Number of genes in go term

    n3 = length(unique(unlist(btms))) - n2 # Number of genes not in BTM

    n4 = length(symbol_list) - n1 # Number of genes correlated that aren't in go term

    hypergeo_res = c(hypergeo_res,phyper(n1-1,n2,n3,n4,lower.tail = FALSE))
  }
  names(hypergeo_res) = names(btms)
  adjusted_p_btm = p.adjust(hypergeo_res,method="fdr")
  #adjusted_p_btm = adjusted_p_btm[adjusted_p_btm<0.1]
  enr_results_btm = adjusted_p_btm








  #### Filtering & formatting results
  frame_to_output = data.frame()


  pvals = enr_results
  gos = names(enr_results)

  l1 = genes_in_go_term
  l2 = l1[names(pvals)]

  new_data.df = data.frame(Pathway_ID=names(pvals),
                           Description = as.character(goterms[names(pvals)]),
                           Category = paste("GO -",go2cat[names(pvals),"Category"]),
                           FDR=pvals)
  new_data.df$Genes = l2
  new_data.df$Genes = as.vector( new_data.df$Genes)
  frame_to_output = rbind(frame_to_output,new_data.df)

  if (length(enr_results_btm) > 0){

    pvals = enr_results_btm
    gos = names(enr_results_btm)

    l1 = my_btms
    l2 = l1[names(pvals)]

    new_data.df = data.frame(Pathway_ID=names(pvals),
                             Description = as.character(btm_ann[names(pvals),3]),
                             Category = paste("BTM -",btm_ann[names(pvals),"Module category"]),
                             FDR=pvals)
    new_data.df$Genes = l2
    new_data.df$Genes = as.vector( new_data.df$Genes)
    frame_to_output = rbind(frame_to_output,new_data.df)
  }



  pvals = enr_results_reactome
  gos = names(enr_results_reactome)

  l1 = genes_in_reactome
  l2 = l1[names(pvals)]

  new_data.df = data.frame(Pathway_ID=names(pvals),
                           Description = as.character(unlist(Reactome2name[names(pvals)])),
                           Category = rep("REACTOME",length(pvals)),
                           FDR=pvals)
  new_data.df$Genes = l2
  new_data.df$Genes = as.vector( new_data.df$Genes)


  frame_to_output = rbind(frame_to_output,new_data.df)


  frame_to_output2 <- transform(frame_to_output, Genes = gsub("[c\"()]", "", as.character(Genes)))
  colnames(frame_to_output2)[5] = "Genes"
  frame_to_output3 = frame_to_output2

  ## Creating a column which shows the genes that contribute to the ORA signal
  for (row in 1:nrow(frame_to_output3)){
  x = frame_to_output3[row,5]
  x = strsplit(x," ")
  x = unlist(x)
  x = gsub(",","",x)
  
  if(all(grepl("ENS",x))){
    y = ens2symbol[x,"external_gene_name"]
  } 
  else if (all(grepl("hsa:",x))){
    y = hsa2symbol[x,"external_gene_name"]
  } 
  else if (str_detect(string = x[1],pattern="[A-Za-z]"))
    y = x
  else {
    y = entrez2symbol[x,"external_gene_name"]
  }
  y = paste(y,collapse=", ")
  frame_to_output3[row,5] = y

}



  ## Ordering and subsetting
  frame_to_output3 = frame_to_output3[order(frame_to_output3$FDR,decreasing = FALSE),]
  frame_to_output3 = frame_to_output3[frame_to_output3$FDR < 0.1,]
  frame_to_output3 = frame_to_output3[!is.na(frame_to_output3$Description),]
  frame_to_output3$Description <- gsub("\r", "", frame_to_output3$Description)



  if (nrow(frame_to_output3) == 0 ) {
    message("Sorry, no significant pathways")
    stop() }
  else if (nrow(frame_to_output3) > 0){
    return(frame_to_output3)
  }

  if (write_to_file){
    write.csv(frame_to_output3, paste0(filename, ".csv"), row.names = FALSE)
  }


}
