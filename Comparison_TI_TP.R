library(readr)
library(tibble)
library(UpSetR)
library(ggvenn)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(lme4)
library(pheatmap)

setwd('C:/Users/lukfa/Desktop/Thesis/Res_Tables')

dat.reg <- readRDS("dat_reg.RDS")
pheno <- readRDS("pheno.RDS")
pheno$Time_point <- paste0('t_',pheno$Time_point)

pheno <- pheno[rownames(pheno) != c('RNAB_PPMI3452_7426_da65_v1'),]
pheno <- pheno[rownames(pheno) != c('RNAB_PPMI3475_1575_da65_v2'),]
pheno <- pheno[rownames(pheno) != c('RNAB_PPMI51440_4204_da65_v1'),]

dat.reg <- dat.reg[,colnames(dat.reg) != c('RNAB_PPMI3452_7426_da65_v1')]
dat.reg <- dat.reg[,colnames(dat.reg) != c('RNAB_PPMI3475_1575_da65_v2')]
dat.reg <- dat.reg[,colnames(dat.reg) != c('RNAB_PPMI51440_4204_da65_v1')]

Res_TI <- read_csv('all_in_one_with_comparison_mixmod_Res_TI.csv')

Res_TP <- read_csv('all_in_one_with_comparison_mixmod_Res_TP.csv')

Res_TI <- Res_TI %>% 
  column_to_rownames("...1")
Res_TP <- Res_TP %>% 
  column_to_rownames("...1")


#Res_TI$FDR <- p.adjust(Res_TI$`LRRK2-GBA`, method = 'BH')
#Res_TP$FDR <- p.adjust(Res_TP$`LRRK2-GBA`, method = 'BH')

Res_TP$FDR <- p.adjust(Res_TP$`LRRK2-GBA`, method = 'BH')
Res_TI$FDR <- p.adjust(Res_TI$`LRRK2-GBA`, method = 'BH')

# Trajectory Inference Genes ....


sorted_TI <- Res_TI[Res_TI$`LRRK2-GBA` < 0.05,]

result_df_sorted_TI <- sorted_TI[order(sorted_TI$`LRRK2-GBA`),]

ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")


identical(rownames(pheno),colnames(dat.reg))

for (x in 1:20) {
  subset <- as.data.frame(dat.reg[rownames(result_df_sorted_TI[x,]),])
  
  conversion_map <-  getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                           filters = "ensembl_gene_id",
                           values = rownames(result_df_sorted_TI[x,]),
                           mart = ensembl)
  
  gene_name <- rownames(result_df_sorted_TI[x,])
  
  if (dim(conversion_map)[1] != 0 & !(rownames(result_df_sorted_TI[x,]) %in% c("ENSG00000279114","ENSG00000272661","ENSG00000228818","ENSG00000234203","ENSG00000256591", "ENSG00000283312", "ENSG00000276259")))  { 
    
    if (conversion_map$external_gene_name != '') {
      
      gene_name <- conversion_map$external_gene_name
      
    } else {
      
      gene_name <- conversion_map$ensembl_gene_id
      
    } }
  
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
    paste0('TI/',gene_name, '.png'), dpi = 300
  )
  print(conversion_map$ensembl_gene_id)
  print(result_df_sorted_TI[x,'LRRK2-GBA'])
}

# Time Points Genes ....

sorted_TP <- Res_TP[Res_TP$comparison_pval< 0.05,]

result_df_sorted_TP <- sorted_TP[order(sorted_TP$comparison_pval),]

for (x in 1:20) {
  subset <- as.data.frame(dat.reg[rownames(result_df_sorted_TP[x,]),])
  
  conversion_map <-  getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                           filters = "ensembl_gene_id",
                           values = rownames(result_df_sorted_TP[x,]),
                           mart = ensembl)
  
  gene_name <- rownames(result_df_sorted_TP[x,])
  
  if (dim(conversion_map)[1] != 0 & !(rownames(result_df_sorted_TP[x,]) %in% c("ENSG00000228818","ENSG00000234203","ENSG00000256591", "ENSG00000283312")))  { 
    
  if (conversion_map$external_gene_name != '') {
    
    gene_name <- conversion_map$external_gene_name
    
  } else {
    
    gene_name <- conversion_map$ensembl_gene_id
    
  } }
  
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
    paste0('TP/',gene_name, '.png'), dpi = 300
  )
  print(gene_name)
  print(result_df_sorted_TI[x,'LRRK2-GBA'])
}


x <- list(
  TP = rownames(sorted_TP), 
  TI = rownames(sorted_TI)
)

ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF"),
  stroke_size = 0.5, set_name_size = 4
)

all_genes <- unique(c(rownames(sorted_TP), rownames(sorted_TI)))

# Initialize dataframe
gene_classification <- data.frame(
  Gene = all_genes,
  Classification = rep(NA, length(all_genes)),
  stringsAsFactors = FALSE
)

# Classify genes
gene_classification$Classification <- sapply(gene_classification$Gene, function(gene) {
  if (gene %in% rownames(sorted_TP) & gene %in% rownames(sorted_TI)) {
    return("Time Points & Trajectory Inference")
  } else if (gene %in% rownames(sorted_TP)) {
    return("Time Points")
  } else if (gene %in% rownames(sorted_TI)) {
    return("Trajectory Inference")
  } else {
    return(NA)
  }
})


rownames(gene_classification) <- gene_classification$Gene

#Strong Time Effect 

breaks <- seq(-1, 1, length.out = 101)

phenopheatmap <- pheno %>% filter(subgroup == 'Healthy Control' | subgroup == 'GBA' | subgroup == 'LRRK2' )
pseudotime_cond <- as.data.frame(pseudotime)
pseudotime_cond <- pseudotime[rownames(pseudotime_cond) %in% rownames(phenopheatmap), , drop = FALSE]

datepheatmap <- dat.reg[rownames(dat.reg) %in% rownames(sorted_TP), rownames(pseudotime_cond)]

identical(rownames(pseudotime_cond),colnames(datepheatmap))

match_order <- match(colnames(datepheatmap), rownames(phenopheatmap))

phenopheatmap <- phenopheatmap[match_order,]

identical(rownames(phenopheatmap),colnames(datepheatmap))

match_order <- match(rownames(datepheatmap), rownames(gene_classification))

gene_classification <- gene_classification[match_order,]

identical(rownames(gene_classification),rownames(datepheatmap))

annotation_row <- data.frame(LMMS = gene_classification$Classification)
rownames(annotation_row) <- rownames(datepheatmap)
identical(rownames(annotation_row),rownames(datepheatmap))

annotation_col <- data.frame(Pseudotime = pseudotime_cond$`traj$time`, Time_Points = phenopheatmap$Time_point, Phenotype = phenopheatmap$subgroup)
rownames(annotation_col) <- rownames(pseudotime_cond)
identical(rownames(annotation_col),colnames(datepheatmap))

time_points_colors <- c("t_0" = "coral1", "t_25" = "seagreen2", "t_65" = "cyan3")
# You can customize the color mapping as per your requirements.

# Combine annotation colors
annotation_colors <- list(
  Pseudotime = colorRampPalette(c("aliceblue", "blue"))(100),  # Example for Pseudotime if needed
  Time_Points = time_points_colors, 
  LMMS = c("Time Points" = "skyblue1", "Trajectory Inference" = "indianred", "Time Points & Trajectory Inference" = "palevioletred1"), 
  Phenotype = c("Healthy Control" = "gold","GBA" = "thistle", "LRRK2" = "darkorange")
)


pheatmap(datepheatmap, 
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         show_colnames = FALSE, 
         show_rownames = FALSE, 
         main = 'Significant Genes', breaks = breaks)

# Not ordering the phenotypes 

phenopheatmap <- pheno %>% filter((subgroup == 'GBA' | subgroup == 'LRRK2'| subgroup == 'Healthy Control'))
datepheatmap <- dat.reg[rownames(dat.reg) %in% rownames(sorted_TP) & rownames(dat.reg) %in% rownames(sorted_TI) , rownames(phenopheatmap)]
annotation <- data.frame(phenotype = phenopheatmap$subgroup)
rownames(annotation) <- rownames(phenopheatmap)
identical(rownames(annotation),colnames(datepheatmap))
annotation$phenotype <- factor(annotation$phenotype, levels = c("Healthy Control", "GBA", "LRRK2"))
annotation <- annotation[order(annotation$phenotype), , drop = FALSE]

match_order <- match(rownames(annotation),colnames(datepheatmap))
datepheatmap <- datepheatmap[,match_order]
identical(rownames(annotation),colnames(datepheatmap))

pheatmap(datepheatmap, 
         annotation_col = annotation, 
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE,
         show_colnames = FALSE, 
         show_rownames = FALSE, 
         main = paste0("Trajectory Inference "), breaks = breaks)



# trying to sketch patterns in conditions ... 
condition <- 'Healthy Control'
breaks <- seq(-1, 1, length.out = 101)

gene_onlyTIorTP <- gene_classification[gene_classification$Classification == 'Time Points',]#| gene_classification$Classification == 'Trajectory Inference' ,]

phenopheatmap <- pheno %>% filter((subgroup == condition))
pseudotime_cond <- as.data.frame(pseudotime)
pseudotime_cond <- pseudotime[rownames(pseudotime_cond) %in% rownames(phenopheatmap), , drop = FALSE]

datepheatmap <- dat.reg[rownames(dat.reg) %in% rownames(gene_onlyTIorTP), rownames(pseudotime_cond)]

mart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
gene_info <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'), 
                   filters = 'ensembl_gene_id', 
                   values = rownames(datepheatmap), 
                   mart = mart)

datepheatmap <- as.data.frame(datepheatmap)
datepheatmap$ensembl_gene_id <- rownames(datepheatmap)
merged_df <- full_join(datepheatmap, gene_info, by = "ensembl_gene_id")

# Replace NA values in hgnc_symbol with the ensembl_gene_id
merged_df$hgnc_symbol <- ifelse(is.na(merged_df$hgnc_symbol) | merged_df$hgnc_symbol == "", 
                                merged_df$ensembl_gene_id, 
                                merged_df$hgnc_symbol)

# Set hgnc_symbol as row names
rownames(merged_df) <- merged_df$hgnc_symbol

datepheatmap <- merged_df


pseudotime_cond <- pseudotime_cond[order(pseudotime_cond$`traj$time`), , drop = FALSE]
match_order <- match(rownames(pseudotime_cond),colnames(datepheatmap))
datepheatmap <- datepheatmap[,match_order]

match_order <- match(rownames(datepheatmap), rownames(gene_classification))
gene_classification <- gene_classification[match_order,]
rownames(gene_classification) <- rownames(datepheatmap)

annotation_row <- data.frame(LMMS = gene_classification$Classification)
rownames(annotation_row) <- rownames(datepheatmap)
identical(rownames(annotation_row),rownames(datepheatmap))

annotation_col <- data.frame(Pseudotime = pseudotime_cond$`traj$time`)
rownames(annotation_col) <- rownames(pseudotime_cond)
identical(rownames(annotation_col),colnames(datepheatmap))


annotation_colors <- list(
  Pseudotime = colorRampPalette(c("aliceblue", "blue"))(100)
  #LMMS = c("Time Points & Trajectory Inference" = "indianred")
)


pheatmap(datepheatmap, 
         annotation_col = annotation_col, 
         #annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE,
         show_colnames = FALSE, 
         show_rownames = TRUE, 
         main = condition, breaks = breaks)

###########################################################################



breaks <- seq(-1, 1, length.out = 101)

phenopheatmap <- pheno %>% filter(subgroup == 'Healthy Control' | subgroup == 'GBA' | subgroup == 'LRRK2' )
pseudotime_cond <- as.data.frame(pseudotime)
pseudotime_cond <- pseudotime[rownames(pseudotime_cond) %in% rownames(phenopheatmap), , drop = FALSE]

datepheatmap <- dat.reg[rownames(dat.reg) %in% rownames(sorted_TP) | rownames(dat.reg) %in% rownames(sorted_TI), rownames(pseudotime_cond)]
identical(rownames(pseudotime_cond),colnames(datepheatmap))

pseudotime_cond <- pseudotime_cond[order(pseudotime_cond$`traj$time`), , drop = FALSE]
match_order <- match(rownames(pseudotime_cond),colnames(datepheatmap))
datepheatmap <- datepheatmap[,match_order]
identical(rownames(pseudotime_cond),colnames(datepheatmap))

match_order <- match(colnames(datepheatmap), rownames(phenopheatmap))
phenopheatmap <- phenopheatmap[match_order,]
identical(rownames(phenopheatmap),colnames(datepheatmap))

match_order <- match(rownames(datepheatmap), rownames(gene_classification))
gene_classification <- gene_classification[match_order,]
identical(rownames(gene_classification),rownames(datepheatmap))

annotation_row <- data.frame(LMMS = gene_classification$Classification)
rownames(annotation_row) <- rownames(datepheatmap)
identical(rownames(annotation_row),rownames(datepheatmap))

annotation_col <- data.frame(Pseudotime = pseudotime_cond$`traj$time`, Time_Points = phenopheatmap$Time_point, Phenotype = phenopheatmap$subgroup)
rownames(annotation_col) <- rownames(pseudotime_cond)
identical(rownames(annotation_col),colnames(datepheatmap))

time_points_colors <- c("t_0" = "coral1", "t_25" = "seagreen2", "t_65" = "cyan3")
# You can customize the color mapping as per your requirements.

# Combine annotation colors
annotation_colors <- list(
  Pseudotime = colorRampPalette(c("aliceblue", "blue"))(100),  # Example for Pseudotime if needed
  Time_Points = time_points_colors, 
  LMMS = c("Time Points" = "skyblue1", "Trajectory Inference" = "indianred", "Both" = "palevioletred1"), 
  Phenotype = c("Healthy Control" = "yellow", "GBA" = "cadetblue", "LRRK2" = "darkorange")
)


pheatmap(datepheatmap, 
         annotation_col = annotation_col, 
         annotation_row = annotation_row,
         annotation_colors = annotation_colors,
         cluster_cols = FALSE,
         show_colnames = FALSE, 
         show_rownames = FALSE, 
         main = 'Significant Genes', breaks = breaks)




Genes <- c('FRA10AC1', 'LINC00115', 'NEFM', 'YBX1P1', 'SFRP2', 'IMMP2L','PPP4C','RFX3','AGGF1')
#non_coding_genes <-  c('PPP4C','RFX3','AGGF1')

explore <-  getBM(attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                         filters = "external_gene_name",
                         values = Genes,
                         mart = ensembl)








subset <- as.data.frame(dat.reg['lnc-CD96-4',])
  
gene_name <- 'lnc-CD96-4'
  
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
    paste0('Daje/',gene_name, '.png'), dpi = 300
  )
  print(gene_name)
  print(result_df_sorted_TI[x,'LRRK2-GBA'])





















