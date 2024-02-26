
setwd("/projects/abv/GRC/leejs/COVID19_coloc_paper/shiny/") # magnus path
# setwd("/Users/leejs/OneDrive - AbbVie Inc (O365)/Documents/CGProjects/COVID-19 Coloc/covid19_coloc")

library(tidyverse)

#coloc_summary <- read.table("summary.per_colocalized_instance.April2023.txt", header=TRUE, sep="\t")
coloc_summary <- read.table("summary.per_colocalized_instance.November2023.txt", header=TRUE, sep="\t")

covid_phenotypes <- c(
  "A2" = "A2: Confirmed COVID-19 with severe Respiratory Symptoms vs Population",
  "B2" = "B2: Hospitalized Individuals due to COVID-19 vs Population",
  "C2" = "C2: Confirmed COVID-19 vs Population"
)

# QTL Type
table(coloc_summary$QTL_type)
# eQTL mQTL pQTL sQTL 
#  679  669   14 1329 

table(coloc_summary$QTL_type, coloc_summary$MP_level)
#       CpG gene protein transcript
# eQTL    0  679       0          0
# mQTL  669    0       0          0
# pQTL    0    0      14          0
# sQTL    0    0       0       1329

# Clump
table(coloc_summary$GWAS_clump)
length(unique(coloc_summary$GWAS_clump))

clumps_df <- distinct(as.data.frame(str_split(coloc_summary$GWAS_clump, "\\.", simplify = TRUE)))
colnames(clumps_df) <- c("Clump", "Chrom", "Start", "End")

# explore eqtl first
eqtl <- coloc_summary %>% filter(QTL_type == "eQTL")

# select gene

# plot 1: plot significant QTL map 
g <- sample(eqtl$Gene, 1)

plot_title <- paste(
  "Significant QTL Map for Gene", g
)

ggplot(eqtl %>% filter(Gene == g)) +
  facet_grid(
    cols = c(vars(GWAS), vars(GWAS_clump)),
    labeller = labeller(
      GWAS = function(x) covid_phenotypes[x]
    )
  ) +
  geom_point(aes(x = QTL_map, y = PP4)) +
  scale_x_discrete(labels = function(x) str_trunc(x, 20, side = "center")) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ggtitle(plot_title)
  
# plot 2: heatmap of the posteriors of gene in the same clump vs qtl map
clumps <- unique(eqtl$GWAS_clump[eqtl$Gene == g])
eqtl_clump <- eqtl[eqtl$GWAS_clump %in% clumps,]

ggplot(eqtl_clump) +
  geom_tile(aes(x = Gene, y = QTL_map, fill = PP4), color = "grey")

plot_heatmap_for_clump(eqtl, clumps, clumps_df)

# plot heatmap of PP4 for each Gene and QTL map within the same clump
eqtl_clump <- eqtl %>%
  filter(GWAS_clump %in% clumps)
plot_title <- "PP4 by Gene and QTL Map"

ggplot(eqtl_clump) +
  geom_tile(aes(x = Gene, y = QTL_map, fill = PP4), color = "grey") +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10)
  ) +
  scale_y_discrete(labels = function(x) str_trunc(x, 30, side = "center")) +
  ggtitle(plot_title,paste(clumps_df$Clump[clumps_df$Clump_name %in% clumps]))

## Advance the heatmap?
source("R/1.qtl_summary_data.R")

myclump <- clumps_df$Clump_name[clumps_df$Clump_number == 11]

qtl_clump <- coloc_summary %>% filter(GWAS_clump == myclump, QTL_type == "eQTL")

# identify QTL maps that are not unique and append numbers
qtl_clump$QTL_map_unique <- qtl_clump$QTL_map

non_unique <- table(qtl_clump$QTL_map)
non_unique <- names(non_unique)[non_unique > 1]

if (length(non_unique) > 0) {
  for (qtl_map in non_unique) {
    i <- which(qtl_clump$QTL_map_unique == qtl_map)
    qtl_clump$QTL_map_unique[i] <- paste(qtl_clump$QTL_map_unique[i], seq_len(length(i)), sep = ".")
  }
}

library(heatmaply)

qtl_clump_wide <- qtl_clump %>% pivot_wider(id_cols = "QTL_map_unique", names_from = "Gene", values_from = PP4)
qtl_clump_wide <- as.data.frame(qtl_clump_wide)
rownames(qtl_clump_wide) <- qtl_clump_wide$QTL_map_unique
qtl_clump_wide$QTL_map_unique <- NULL

hover_info <- matrix(data=NA, nrow = nrow(qtl_clump_wide), ncol = ncol(qtl_clump_wide))
rownames(hover_info) <- rownames(qtl_clump_wide); colnames(hover_info) <- colnames(qtl_clump_wide)

for (i in 1:nrow(qtl_clump)) {
  i_gene <- qtl_clump$Gene[i]
  i_qtl_map <- qtl_clump$QTL_map_unique[i]
  hover_info[i_qtl_map,i_gene] <- 
    paste0(
      qtl_clump$QTL_map[i], "\n",
      "Gene: ", qtl_clump$Gene[i], "\n",
      "MP: ", qtl_clump$MP[i], "\n",
      "GWAS: ", qtl_clump$GWAS[i], "\n",
      "Population: ", qtl_clump$population[i], "\n",
      "PP4: ", qtl_clump$PP4[i]
    )
}

heatmaply(
  qtl_clump_wide,
  Rowv = FALSE,
  Colv = FALSE,
  colors = plasma(256, alpha = 1, begin = 0, end = 1, direction = 1), custom_hovertext = hover_info,
  plot_method = "plotly",
  main = paste0("PP4 Heatmap of Genes in ", myclump)
)

## 

filter_function <- function(x, cn, val) {
  x %>% filter(.[[cn]] == val)
}

filter_function(qtl_clump, "Gene", "ABO")

## 
datf <- data.frame(
  Gene = c("A", "A", "A", "A", "B", "B", "B", "C", "C", "C"),
  QTL_map = c('a','a','b','c','a','b','b','b','b','c')
)

datf$QTL_map_unique <- datf$QTL_map

non_unique <- table(datf$QTL_map)
non_unique <- names(non_unique)[non_unique > 1]

if (length(non_unique) > 0) {
  for (qtl_map in non_unique) {
    i <- which(datf$QTL_map_unique == qtl_map)
    datf$QTL_map_unique[i] <- paste(datf$QTL_map_unique[i], seq_len(length(i)), sep = ".")
  }
}

qtl_clump %>% group_by(Gene) %>% mutate(QTL_map_unique = append_suffix_to_duplicate(QTL_map)) %>% View

append_suffix_to_duplicate <- function(x) {
  duplicates <- unique(x[duplicated(x)])
  for (dup in duplicates) {
    i <- which(x == dup)
    x[i] <- paste0(x[i], c("", paste0(".",seq_len(length(i))[-1])))
  }
  return(x)
}

