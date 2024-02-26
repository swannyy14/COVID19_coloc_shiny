setwd("/projects/abv/GRC/leejs/COVID19_coloc_paper/shiny")

# read and clean GTF
gtf <- read.table("/projects/abv/users/degner/Project/COVID_HG/for_Reza_and_John/gencode.v38.annotation.gene_level.gtf", sep = "\t")

colnames(gtf) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

gtf$gene_id <- str_match(gtf$attribute, "^gene_id (ENSG[[:digit:].]+)")[,2]
gtf$gene_id_nover <- str_remove(gtf$gene_id, "\\.[[:digit:]]+")
gtf$gene_name <- str_match(gtf$attribute, "gene_name ([\"[:alnum:]:\\-._/]+)")[,2]

gtf <- gtf[,c("chr", "start", "end", "gene_name", "gene_id", "gene_id_nover")]

# read summary file
#coloc <- read.table("summary.per_colocalized_instance.April2023.txt", header=TRUE, sep="\t")
coloc <- read.table("summary.per_colocalized_instance.November2023.txt", header=TRUE, sep="\t")

# compare genes in the coloc summary file versus gtf file
coloc_genes <- unique(coloc$Gene)
coloc_genes <- unique(unlist(str_split(coloc_genes, ";")))

# number of overlapping genes when searched by gene name
overlapping <- coloc_genes %in% gtf$gene_name
table(overlapping)

# many genes are in ENSG format
gene_no_overlap <- coloc_genes[!overlapping]
print(gene_no_overlap)

# try to map using ENSG ID
overlapping2 <- gene_no_overlap %in% gtf$gene_id_nover
table(overlapping2)

# genes that do not map with gene id or ensemble gene id
gene_no_overlap2 <- gene_no_overlap[!overlapping2]
print(gene_no_overlap2)

# any duplicate gene IDs?
gtf2 <- gtf %>% filter(gene_name %in% coloc_genes | gene_id_nover %in% coloc_genes)

duplicate_genes <- table(gtf2$gene_name)
duplicate_genes <- names(duplicate_genes)[duplicate_genes > 1]
# there are no duplicate genes in the set that overlaps with the coloc summary
# however there are many duplicate gene names in the entire gtf

write.csv(gtf, file = "gene_positions.csv", row.names = FALSE)

