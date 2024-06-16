
library(SCORPIUS)
library(slingshot)
library(TSCAN)
library(tradeSeq)
library(mclust, quietly = TRUE)
library(RColorBrewer)
library(gam)
library(SingleCellExperiment)
library(ggfortify)
library(ggplot2)
library(gridExtra)
library(sva)
library(readr)
library(stringr)
library(dyno)
library(tidyverse)
library(biomaRt)
library(edgeR)
library(dplyr)
library(Matrix)
library(limpca)
library(clusterProfiler)
library(org.Hs.eg.db)
library(UpSetR)
library(GGally)
library(ggplot2)   
library(corrplot)
library(grDevices)

###############################################################
#                                                        ######
# import data                                            ######
#                                                        ######
###############################################################

setwd('C:/Users/lukfa/Desktop/Thesis/Res_Tables')

dat.reg <- readRDS("dat_reg.RDS")
pheno <- readRDS("pheno.RDS")
pheno$Time_point <- paste0('t_',pheno$Time_point)

pca_resultall <- prcomp(t(dat.reg), scale. = TRUE)
autoplot(pca_resultall, data = pheno, colour = 'subgroup', label = F, label.size = 3) + theme_minimal()

gene_ids <- rownames(dat.reg)

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                   filters = 'ensembl_gene_id', 
                   values = gene_ids, 
                   mart = mart)

missed_ids <- setdiff(gene_ids, gene_info$ensembl_gene_id)
if (length(missed_ids) > 0) {
  cat("No symbols found for the following IDs:", paste(missed_ids, collapse=", "), "\n")
}
# Notify if any IDs were not found
if (length(missed_ids) > 0) {
  cat("No symbols found for the following IDs:", paste(missed_ids, collapse=", "), "\n")
}

# Create a map from Ensembl IDs to HGNC symbols
names_map <- setNames(gene_info$hgnc_symbol, gene_info$ensembl_gene_id)

# Preparing a comprehensive dataframe with all gene IDs and their corresponding symbols
complete_gene_info <- data.frame(ensembl_gene_id = gene_ids)
complete_gene_info$hgnc_symbol <- names_map[complete_gene_info$ensembl_gene_id]

# Replace NA with a placeholder or the original Ensembl ID in 'hgnc_symbol' for unmatched IDs
complete_gene_info$hgnc_symbol[is.na(complete_gene_info$hgnc_symbol)] <- complete_gene_info$ensembl_gene_id[is.na(complete_gene_info$hgnc_symbol)]
complete_gene_info <- na.omit(complete_gene_info)

complete_gene_info <- complete_gene_info %>%
  filter(!is.na(hgnc_symbol) & hgnc_symbol != "")

dat.reg <- as.data.frame(dat.reg)
dat.reg$ensembl_gene_id <- rownames(dat.reg)
dat.reg <- left_join(dat.reg, complete_gene_info, by = 'ensembl_gene_id')
sum(duplicated(dat.reg$hgnc_symbol))
dat.reg <- distinct(dat.reg, hgnc_symbol, .keep_all = TRUE)
missing_values <- is.na(dat.reg$hgnc_symbol)
sum(missing_values)
dat.reg <- dat.reg[dat.reg$ensembl_gene_id != 'ENSG00000249456',]
rownames(dat.reg) <- dat.reg$hgnc_symbol 

# Remove unnecessary columns
dat.reg <- dat.reg %>% dplyr::select(-ensembl_gene_id, -hgnc_symbol)


match_order <- match(colnames(dat.reg), rownames(pheno))

pheno_pseudo <- pheno[match_order,]

identical(rownames(pheno),colnames(dat.reg))

pheno <- pheno[,c('Time_point','PATNO','subgroup')]
pheno$Time_point <- paste0('t_',pheno$Time_point)

df <- cbind(pheno, t(dat.reg))

colnames(df)[colnames(df) == 'PATNO'] <- 'ID'

colnames(df)[colnames(df) == 'subgroup'] <- 'GROUP'

df$GROUP <- relevel(df$GROUP, ref = 'Healthy Control')

df$GROUP <- as.character(df$GROUP)
df$Time_point <- as.character(df$Time_point)
df$ID <- as.character(df$ID)

library(reshape2)

# Melt the dataframe to long format
df_long <- melt(df, id.vars = c("ID", "GROUP","Time_point"), variable.name = "variable", value.name = "value")

saveRDS(df_long, file = "ALASCA.rds")
