library(tidyverse)

# coloc_summary: summary statistics from colocalization
# clumps_df: meta data about the clumps


coloc_summary <- read.table("summary.per_colocalized_instance.November2023.v2.txt", header=TRUE, sep="\t")

covid_phenotypes <- c(
  "A2" = "A2: Confirmed COVID-19 with severe Respiratory Symptoms vs Population",
  "B2" = "B2: Hospitalized Individuals due to COVID-19 vs Population",
  "C2" = "C2: Confirmed COVID-19 vs Population"
)

# extract clump summary
unique_clumps <- unique(coloc_summary$GWAS_clump)
clumps_df <- as.data.frame(str_split(unique_clumps, "\\.", simplify = TRUE))
colnames(clumps_df) <- c("Clump", "Chrom", "Start", "End")
clumps_df$Clump_name <- unique_clumps
clumps_df$Clump_number <- as.integer(str_extract(clumps_df$Clump, "[[:digit:]]+$"))
clumps_df <- clumps_df[order(clumps_df$Clump_number),]

# extract eQTL summary
eqtl <- coloc_summary %>% filter(QTL_type == "eQTL")
eqtl_clumps_df <- clumps_df[clumps_df$Clump_name %in% eqtl$GWAS_clump,]

# genes available in eqtl summary
eqtl_genes <- sort(unique(eqtl$Gene))

# gene position data frame
gene_pos_df <- read.csv("gene_positions.csv")
gene_pos_df$chr <- factor(gene_pos_df$chr, paste0("chr", c(1:22, "X", "Y", "M")))
gene_pos_df <- gene_pos_df[order(gene_pos_df$chr, gene_pos_df$start),]
gene_pos_df$index <- 1:nrow(gene_pos_df)

# qtl metadata
qtl_metadata <- read.table("Included.QTL.maps.metadata.lite.tsv", sep = "\t", header=TRUE)

# add qtl metadata to coloc summary
coloc_summary <- coloc_summary %>% left_join(qtl_metadata[,c("Map", "Predominant_Cell")], by = join_by(QTL_map == Map))
