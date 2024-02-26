plot_signif_qtl_map_for_gene <- function(qtl_summary, gene, gene_summary_level, covid_phenotypes, qtl_type) {
  
  # make sure gene exists in the summary data subset (not performing this check results in error in the output - shiny still runs though)
  if (!(gene %in% qtl_summary[[gene_summary_level]])) {
    return(NULL)
  }
  
  # plot significant QTL map across all QTL map for chosen gene
  plot_title <- paste(
    paste0("Significant ", qtl_type, " Map for Gene ", gene)
  )
  
  qtl_subset <- qtl_summary %>% filter(.[[gene_summary_level]] == gene)
  clumps <- paste(unique(qtl_subset$GWAS_clump), collapse = ", ")
  
  # Add extra text to information that comes up when hovering over each point
  extra_hovertext_col <- ifelse(gene_summary_level == "Gene", "MP", "Gene")
  qtl_subset <- qtl_subset %>%
    mutate(
      extra_hovertext = 
        paste(
          paste("QTL Type:", QTL_type),
          paste("Clump:", GWAS_clump),
          paste0(extra_hovertext_col, ": ", .[[extra_hovertext_col]]),
          paste("Population:", population),
          paste("Top Coloc Variant:", top_coloc_variant),
          paste("Top Coloc Variant PP4:", top_coloc_variant_PP4),
          sep = "\n"
        )
    )
  
  # scatter plot of PP4 for each QTL map for a fixed Gene / MP
  if (qtl_type == "All") {
    g <- ggplot(qtl_subset, aes(color = QTL_type))
  } else {
    g <- ggplot(qtl_subset, aes())
  }
  
  g +
    facet_grid(
      cols = c(vars(GWAS)),
      labeller = labeller(
        GWAS = function(x) covid_phenotypes[x]
      )
    ) +
    # following line will produce warning, but "text" aesthetics is needed for hovertext
    geom_point(aes(x = QTL_map, y = PP4, text = extra_hovertext)) + 
    scale_x_discrete(labels = function(x) str_trunc(x, 30, side = "center")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)) +
    ggtitle(plot_title, clumps)
}

plotly_signif_qtl_map_for_gene <- function(qtl_summary, gene, gene_summary_level, covid_phenotypes, qtl_type) {
  if (!(gene %in% qtl_summary[[gene_summary_level]])) {
    return(NULL)
  }
  
  # plot significant QTL map across all QTL map for chosen gene
  g <- plot_signif_qtl_map_for_gene(qtl_summary, gene, gene_summary_level, covid_phenotypes, qtl_type)
  
  ggplotly(g)
}

# deprecated heatmap new heatmap has additional functionality such as interactivity
plot_heatmap_for_clump <- function(qtl_summary, clumps, clumps_df) {
  # plot heatmap of PP4 for each Gene and QTL map within the same clump
  qtl_clump <- qtl_summary %>%
    filter(GWAS_clump %in% clumps)
  plot_title <- "PP4 by Gene and QTL Map"

  ggplot(qtl_clump) +
    geom_tile(aes(x = Gene, y = QTL_map, fill = PP4), color = "grey") +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)
    ) +
    scale_y_discrete(labels = function(x) str_trunc(x, 30, side = "center")) +
    ggtitle(plot_title,paste(clumps_df$Clump[clumps_df$Clump_name %in% clumps]))
}

plot_heatmaply_for_clump <- function(qtl_summary, clumps, clumps_df, gene_pos_df, height_scale = 1, choose_max_pp4_per_gene = FALSE) {
  # display heatmap of PP4 for Tissue vs Gene in selected clump
  # qtl_summary is the data summary file filtered to selected qtl type
  # clumps is the full clump name
  # clumps_df (created in 1.qtl_summary_data.R) containing extracted information from clump
  # height_scale changes the scale of the heatmap output (controlled by ui)
  # choose_max_pp4_per_gene: if TRUE, then for Genes with duplicate QTL map, entry with max PP4 is used. Otherwise,
  #   the QTL maps are renamed (e.g. <QTL_type>.1, <QTL_type>.2, ...)
  
  require(heatmaply)
  
  if (!any(clumps %in% qtl_summary$GWAS_clump)) {
    return(NULL)
  }
  
  # get the clump number (used in subtitle)
  clump <- clumps_df$Clump[clumps_df$Clump_name %in% clumps]
  clump <- paste(clump, collapse = ", ")
  
  # filter coloc summary to selected clump
  qtl_clump <- qtl_summary %>% filter(GWAS_clump %in% clumps)
  
  # Treat NA as another character
  qtl_clump$Gene <- ifelse(is.na(qtl_clump$Gene), "NA", qtl_clump$Gene)
  
  # append Predominant_cell metadata to QTL_map
  append_data <- "Predominant_Cell"
  if (!is.null(append_data)) {
    qtl_clump$QTL_map <- ifelse(is.na(qtl_clump[[append_data]]), qtl_clump$QTL_map,  paste0(qtl_clump$QTL_map, ".", qtl_clump[[append_data]]))
  }
  
  # check what to do with duplicate QTL map within each gene
  if (choose_max_pp4_per_gene) {
    # count number of Molecular Phenotype collapsed
    qtl_clump_MP_count <- qtl_clump %>%
      group_by(Gene, QTL_map) %>%
      summarise(mp_count = n(), .groups = "drop")
    
    # choose observation with max PP4 for each gene-qtl_map group
    qtl_clump <- qtl_clump %>% 
      group_by(Gene, QTL_map) %>% 
      filter(PP4 == max(PP4)) %>%
      ungroup()
    
    # in the rare case that PP4 is the same for multiple MP in the same gene-qtl_map,
    # sample just one
    qtl_clump <- qtl_clump %>% slice_sample(n = 1, by = c("Gene", "QTL_map"))
    
    qtl_clump <- qtl_clump %>%
      left_join(qtl_clump_MP_count, by = c("Gene", "QTL_map"))
    
    qtl_clump$QTL_map_unique <- qtl_clump$QTL_map
  } else {
    # identify QTL maps that are not unique (within each gene) and append numbers
    qtl_clump <- qtl_clump %>%
      group_by(Gene) %>%
      mutate(QTL_map_unique = append_suffix_to_duplicate(QTL_map))
  }
  
  # convert data to a format where each row represents the QTL map and column represents Gene
  qtl_clump_wide <- qtl_clump %>% pivot_wider(id_cols = "QTL_map_unique", names_from = "Gene", values_from = PP4)
  qtl_clump_wide <- as.data.frame(qtl_clump_wide)
  rownames(qtl_clump_wide) <- qtl_clump_wide$QTL_map_unique
  qtl_clump_wide$QTL_map_unique <- NULL
  
  # reorder genes based on TSS
  genes_ordered <- order_genes(colnames(qtl_clump_wide), gene_pos_df = gene_pos_df)
  qtl_clump_wide <- qtl_clump_wide[,genes_ordered,drop=FALSE]
  
  # create a matrix where each entry is a hovertext information
  hover_info <- matrix(data=NA, nrow = nrow(qtl_clump_wide), ncol = ncol(qtl_clump_wide))
  rownames(hover_info) <- rownames(qtl_clump_wide); colnames(hover_info) <- colnames(qtl_clump_wide)
  
  for (i in 1:nrow(qtl_clump)) {
    i_gene <- qtl_clump$Gene[i]
    i_qtl_map <- qtl_clump$QTL_map_unique[i]
    
    hover_info_text <- 
      paste0(
        qtl_clump$QTL_map[i], "\n",
        "QTL Type: ", qtl_clump$QTL_type[i], "\n",
        "Gene: ", qtl_clump$Gene[i], "\n",
        "MP: ", qtl_clump$MP[i], "\n",
        "GWAS: ", qtl_clump$GWAS[i], "\n",
        "Population: ", qtl_clump$population[i], "\n",
        "PP4: ", qtl_clump$PP4[i]
      )
    
    if (choose_max_pp4_per_gene) {
      hover_info_text <- 
        paste0(hover_info_text, "\n", "# Collapsed MP: ", qtl_clump$mp_count[i])
    }
    
    hover_info[i_qtl_map,i_gene] <- hover_info_text
  }
  
  heatmaply(
    qtl_clump_wide,
    Rowv = FALSE,
    Colv = FALSE,
    color = colorRampPalette(c("lightpink","red")),
    #colors = colorspace::heat_hcl(256, alpha = 1),
    custom_hovertext = hover_info,
    plot_method = "plotly",
    main = paste0("PP4 Heatmap of Genes in ", clump),
    height = 400*height_scale
  )
  
}

append_suffix_to_duplicate <- function(x) {
  duplicates <- unique(x[duplicated(x)])
  for (dup in duplicates) {
    i <- which(x == dup)
    x[i] <- paste0(x[i], c("", paste0(".",seq_len(length(i))[-1])))
  }
  return(x)
}

