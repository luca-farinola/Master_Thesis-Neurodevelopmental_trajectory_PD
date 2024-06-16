library(ALASCA)
library(biomaRt)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(pheatmap)
library(pheatmap)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)

setwd("C:/Users/lukfa/Desktop/Thesis/Res_Tables")

df<- readRDS(file = "ALASCA.rds")

df$GROUP[df$GROUP == "Healthy Control"] <- "Control"
df <- df[df$GROUP == "Control"|df$GROUP == "GBA"|df$GROUP == "LRRK2",]

# ALASCA Analysis ----
ASCA_model <- ALASCA(
  df,
  formula = value~Time_point*GROUP + (1|ID),
  separate_effects = TRUE,
  scale_function = "sdt1",
  save = TRUE,
  n_validation_runs = 100,
  validate = TRUE,reduce_dimensions = TRUE, use_Rfast = FALSE
)

# Plot ----

plot(ASCA_model,
     effect = c(1,2),
     component = 1,
     flip_axis = FALSE, type = 'effect'
)

plot(ASCA_model,
     effect = c(1,2),
     component = 2,
     flip_axis = FALSE, type = 'effect'
)


plot(ASCA_model,
     effect = 2,
     component = 2,
     flip_axis = FALSE, type = 'histogram'
)

plot(ASCA_model,
     effect = 2,
     component = 1,
     flip_axis = FALSE, type = 'histogram'
)


plot(ASCA_model,
     effect = 2,
     component = 2, type = 'prediction'
)

# Top absolute value loadings 

df_loading.top <- get_loadings(ASCA_model, component = 2, effect = 2, n_limit = 2)[[1]]

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

genes <-  getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                         filters = "external_gene_name",
                         values = df_loading.top$covars,
                         mart = ensembl)

genes <- genes[genes$ensembl_gene_id != "ENSG00000277909",]
genes <- genes[genes$ensembl_gene_id != "ENSG00000096093",]
genes <- genes[genes$ensembl_gene_id != "ENSG00000276838",]
genes <- genes[genes$ensembl_gene_id != "ENSG00000281690",]

match_order <- match(colnames(dat.reg), rownames(pheno))

pheno <- pheno[match_order,]

identical(rownames(pheno),colnames(dat.reg))

for (x in 1:dim(genes)[2]) {
  
  subset <- as.data.frame(dat.reg[genes$ensembl_gene_id[x],])
  
  colnames(subset) <- 'expression'
  pheno <- pheno[,c('subgroup','Time_point')]
  
  subset <- cbind(pheno,subset)
  
  subset <- subset[subset$subgroup != 'Sporadic',]
  
  subset$subgroup <- factor(subset$subgroup, levels = c("Healthy Control", "GBA", "LRRK2"))
  
  time_points_colors <- c("t_0" = "coral1", "t_25" = "seagreen2", "t_65" = "cyan3")
  
  ggplot(subset, aes(x=subgroup, y= expression, fill=Time_point)) + 
    scale_fill_manual(values = time_points_colors) +
    geom_boxplot() + theme_minimal()
  
  ggsave(
    filename = paste0('Daje/', genes$external_gene_name[x], '.png'), 
    dpi = 300
  )
  print(genes$external_gene_name[x])
}


toploadings_counts <- dat.reg[genes$ensembl_gene_id,]

identical(rownames(pheno),colnames(toploadings_counts))
merged_data <- cbind(t(toploadings_counts), pheno)

top <- get_loadings(ASCA_model, component = 2, effect = 2, n_limit = 10)[[1]]


ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

conversion_top <-  getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                                  filters = "external_gene_name",
                                  values =top$covars ,
                                  mart = ensembl)

phenopheatmap <- pheno %>% filter( subgroup == 'Healthy Control'| subgroup == 'GBA' | subgroup == 'LRRK2' )

dat <- dat.reg[, rownames(phenopheatmap)]
custom_colors = c("Healthy Control" = "gold", "GBA" = "cadetblue", "LRRK2" = "darkorange")

pca <- dat.reg[rownames(dat.reg) %in% conversion_top$ensembl_gene_id, rownames(phenopheatmap)]
pca_result <- prcomp(t(pca), scale. = TRUE)
autoplot(pca_result, data = phenopheatmap, colour = 'subgroup', label = F, label.size = 3) +
  scale_color_manual(values = custom_colors) +
  theme_minimal()

pca_result <- prcomp(t(pca), scale. = TRUE)
autoplot(pca_result, data = phenopheatmap, colour = 'Day', label = F, label.size = 3) +
  theme_minimal()

pca <- dat.reg[, rownames(phenopheatmap)]
pca_result <- prcomp(t(pca), scale. = TRUE)
autoplot(pca_result, data = phenopheatmap, colour = 'subgroup', label = F, label.size = 3) +
  scale_color_manual(values = custom_colors) +
  theme_minimal()

pca_result <- prcomp(t(pca), scale. = TRUE)
autoplot(pca_result, data = phenopheatmap, colour = 'Day', label = F, label.size = 3) +
  theme_minimal()








