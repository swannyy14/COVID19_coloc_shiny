get_gene_position_rank_single <- function(gene, gene_pos_df) {
  # gene_pos_df contains genes and index which is ordered depending on the TSS (chromosome:position)
  # return the associated minimum rank for given gene
  # return 2e6 if gene is character NA
  # return 1e6 if the gene is not found
  
  if (gene == "NA") return(2e6)
  
  gene <- str_split(gene, ";")[[1]]
  gene_rank <- gene_pos_df %>%
    filter(
      gene_name %in% gene | gene_id_nover %in% gene
    ) %>% 
    pull(index) %>% 
    .[1]
  
  if (all(is.na(gene_rank))) {
    return(1e6)
  } else {
    return(gene_rank)
  }
}

get_gene_position_rank_vector <- Vectorize(get_gene_position_rank_single, vectorize.args = "gene")

order_genes <- function(gene, gene_pos_df) {
  gene_index <- get_gene_position_rank_vector(gene, gene_pos_df)
  return(gene[order(gene_index)])
}
