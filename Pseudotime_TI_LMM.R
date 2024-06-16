rm(list=ls())

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
library(lme4)
library(flexplot)
library(lmerTest)
library(multcomp)


###############################################################
#                                                        ######
# import data                                            ######
#                                                        ######
###############################################################

setwd('C:/Users/lukfa/Desktop/Thesis')

counts<-data.table::fread("aggregated_expression/countTable.tsv")
genes_ID <- counts$Geneid 
genes_ID <- sub("\\.\\d+$", "", genes_ID)
rownames(counts)<- genes_ID
counts$Geneid<-NULL
rownames(counts)

meta_pheno <- read_csv('aggregated_expression/pheno_ipsc.csv')

columns <- colnames(counts)
columns <- columns[-1]

ppmi_numbers <- as.integer(str_extract(columns, "(?<=PPMI)\\d+"))
da_numbers <- as.integer(str_extract(columns, "(?<=da)\\d+"))

metadata_all <- data.frame(PATNO = ppmi_numbers, Time_point = da_numbers, patient = columns)

metadata_all <- merge(metadata_all, meta_pheno, by = 'PATNO')
rownames(metadata_all) <- metadata_all$patient
metadata_all$patient <- NULL

pheno<-metadata_all[metadata_all$subgroup=="Healthy Control"|
                      metadata_all$subgroup=="GBA" | metadata_all$subgroup=="LRRK2" | metadata_all$subgroup=="Sporadic" ,]

counts1<-data.frame(counts)
rownames(counts1)<-rownames(counts)
counts<-counts1[,rownames(pheno)]
identical(rownames(pheno), colnames(counts))

###############################################################
#                                                        ######
# Normalizzation                                         ######
#                                                        ######
###############################################################

# Filter out lowly expressed genes 

keep = edgeR::filterByExpr(counts,min.count = 10, min.prop = 0.80)
counts = counts[keep,]
counts<-as.matrix(counts)

dge <- DGEList(counts)

# Normalize the data using TMM method

dge <- calcNormFactors(dge,method =c("TMM"))

# Calculate log2-transformed CPM values

log2.cpm <- edgeR::cpm(dge, log=TRUE)
identical(colnames(log2.cpm), rownames(pheno))

pheno$Day<-as.factor(pheno$Time_point)
pheno$SITE<-as.factor(pheno$SITE)
pheno$SEX<-as.factor(pheno$SEX)
pheno$subgroup<-as.factor(pheno$subgroup)

###############################################################
#                                                        ######
# Surrogate Variable Analysis                            ######
#                                                        ######
###############################################################


mod1 <- model.matrix(~ SEX + age + Day + subgroup, data=pheno)
mod0 <- model.matrix(~1, data=pheno)
  
set.seed(1234)
dat <- as.matrix(log2.cpm)
dat <- na.omit(dat)
  
# Perform SVA
svseq <- sva::sva(dat, mod1, mod0)$sv
svseq <- data.frame(svseq)
colnames(svseq) <- paste0("sv", seq(1, ncol(svseq)))
  
phen <- cbind(pheno, svseq)
  
resid <- function(row, sv1, sv2, sv3, sv4, sv5, SEX, age) {
    fit <- try(lm(row ~ sv1 + sv2 + sv3 + sv4 + sv5 + SEX + age), silent = TRUE)
    if(inherits(fit, 'try-error')) return(rep(NA, length(row)))
    return(fit$residuals)
  }
  
dat.reg <- {
    sv1 <- svseq$sv1
    sv2 <- svseq$sv2
    sv3 <- svseq$sv3
    sv4 <- svseq$sv4
    sv5 <- svseq$sv5
    SEX <- as.factor(phen$SEX)
    age <- as.numeric(phen$age)
    t(apply(dat, 1, resid, sv1, sv2, sv3,sv4,sv5, SEX, age))
  }

pca_resultall <- prcomp(t(dat.reg), scale. = TRUE)
pca_resultall <- prcomp(t(log2.cpm), scale. = TRUE)
plot <- autoplot(pca_resultall, data = pheno, colour = 'Day', label = F, label.size = 3) + theme_minimal()

#dens.reg <- plot(density(dat.reg))
#dens.log2 <- plot(density(log2(counts+0.0001)))


# Save the first density plot
pdf("dens_reg.pdf", width = 8, height = 6)
plot(density(dat.reg))
dev.off()

# Save the second density plot
pdf("dens_log2.pdf", width = 8, height = 6)
plot(density(log2(counts + 0.0001)))
dev.off()

#ggsave("Immages/pca_after_sva.png", plot, dpi = 300)
#ggsave("Immages/dens_after_reg.png", p1, dpi = 300)
#ggsave("Immages/dens_before_reg.png", p2, dpi = 300)


dat.regHC <- dat.reg[,pheno$subgroup == 'Healthy Control']
phenoHC <- pheno[pheno$subgroup == 'Healthy Control',]

dat.regGBA <- dat.reg[,pheno$subgroup == 'GBA']
phenoGBA <- pheno[pheno$subgroup == 'GBA',]

dat.regLRRK2 <- dat.reg[,pheno$subgroup == 'LRRK2']
phenoLRRK2 <- pheno[pheno$subgroup == 'LRRK2',]

dat.regSporadic <- dat.reg[,pheno$subgroup == 'Sporadic']
phenoSporadic <- pheno[pheno$subgroup == 'Sporadic',]

#Removing Outliers

dat.regHC <- dat.regHC[,colnames(dat.regHC) != c('RNAB_PPMI3452_7426_da65_v1')]
phenoHC <- phenoHC[rownames(phenoHC) != c('RNAB_PPMI3452_7426_da65_v1'),]

dat.regSporadic <- dat.regSporadic[,colnames(dat.regSporadic) != c('RNAB_PPMI3475_1575_da65_v2')]
phenoSporadic <- phenoSporadic[rownames(phenoSporadic) != c('RNAB_PPMI3475_1575_da65_v2'),]

dat.regLRRK2 <- dat.regLRRK2[,colnames(dat.regLRRK2) != c('RNAB_PPMI51440_4204_da65_v1')]
phenoLRRK2 <- phenoLRRK2[rownames(phenoLRRK2) != c('RNAB_PPMI51440_4204_da65_v1'),]

dat.reg <- dat.reg[,colnames(dat.reg) != c('RNAB_PPMI3452_7426_da65_v1')]
dat.reg <- dat.reg[,colnames(dat.reg) != c('RNAB_PPMI3475_1575_da65_v2')]
dat.reg <- dat.reg[,colnames(dat.reg) != c('RNAB_PPMI51440_4204_da65_v1')]

pheno <- pheno[rownames(pheno) != c('RNAB_PPMI3452_7426_da65_v1'),]
pheno <- pheno[rownames(pheno) != c('RNAB_PPMI3475_1575_da65_v2'),]
pheno <- pheno[rownames(pheno) != c('RNAB_PPMI51440_4204_da65_v1'),]



#PCA 

pca_resultHC <- prcomp(t(dat.regHC), scale. = TRUE)
autoplot(pca_resultHC, data = phenoHC, colour = 'subgroup', label = F, label.size = 3) +
  scale_color_manual(values = c("Healthy Control" = "#00CFCF")) + theme_minimal()

###############################################################
#                                                        ######
# Projecting GBA                                         ######
#                                                        ######
###############################################################

# Project on HC 

GBA_PCA <- as.matrix(predict(pca_resultHC, newdata = t(dat.regGBA)))
LRRK2_PCA <- as.matrix(predict(pca_resultHC, newdata = t(dat.regLRRK2)))
Sporadic_PCA <- as.matrix(predict(pca_resultHC, newdata = t(dat.regSporadic)))
HC_PCA <- as.matrix(pca_resultHC$x)
projectedPCA_GBA <- rbind(GBA_PCA[,1:2], HC_PCA[,1:2])
projectedPCA_LRRK2 <- rbind(LRRK2_PCA[,1:2], HC_PCA[,1:2])
projectedPCA_Sporadic <- rbind(Sporadic_PCA[,1:2], HC_PCA[,1:2])

#plotting GBA and HC and removing Outliers 

phenoHC_GBA <- pheno[pheno$subgroup=="Healthy Control" | pheno$subgroup=="GBA",]
phenoHC_LRRK2 <- pheno[pheno$subgroup=="Healthy Control" | pheno$subgroup=="LRRK2",]
phenoHC_Sporadic <- pheno[pheno$subgroup=="Healthy Control" | pheno$subgroup=="Sporadic",]
phenoHC_GBA <- phenoHC_GBA[rownames(phenoHC_GBA) != c('RNAB_PPMI3452_7426_da65_v1'),] 
phenoHC_LRRK2 <- phenoHC_LRRK2[rownames(phenoHC_LRRK2) != c('RNAB_PPMI3452_7426_da65_v1'),] 
phenoHC_Sporadic <- phenoHC_Sporadic[rownames(phenoHC_Sporadic) != c('RNAB_PPMI3452_7426_da65_v1'),] 
phenoHC_Sporadic <- phenoHC_Sporadic[rownames(phenoHC_Sporadic) != c('RNAB_PPMI3475_1575_da65_v2'),] 

# Add metadata information 

combined_df_GBA <- merge(projectedPCA_GBA, phenoHC_GBA, by = 0, all.x = TRUE)
combined_df_LRRK2 <- merge(projectedPCA_LRRK2, phenoHC_LRRK2, by = 0, all.x = TRUE)
combined_df_Sporadic <- merge(projectedPCA_Sporadic, phenoHC_Sporadic, by = 0, all.x = TRUE)
rownames(combined_df_GBA) <- combined_df_GBA$Row.names
rownames(combined_df_LRRK2) <- combined_df_LRRK2$Row.names
rownames(combined_df_Sporadic) <- combined_df_Sporadic$Row.names
combined_df_GBA$Row.names <- NULL
combined_df_LRRK2$Row.names <- NULL
combined_df_Sporadic$Row.names <- NULL


# Trajectory With Scorpius 

trajectory_analysis <- function(pca, phenotype) {
  rd <- data.frame(pca[, 1:2])
  traj <- SCORPIUS::infer_trajectory(rd)
  plot <- SCORPIUS::draw_trajectory_plot(rd, as.factor(phenotype$Day), traj$path)
  pseudotime <- as.data.frame(traj$time)
  
  return(list(time = pseudotime,plot = plot))
}

set.seed(1234)
pseudo_HC <- trajectory_analysis(HC_PCA, phenoHC)
pseudo_GBA <- trajectory_analysis(GBA_PCA, phenoGBA)
pseudo_LRRK2 <- trajectory_analysis(LRRK2_PCA, phenoLRRK2)
pseudo_Sporadic <- trajectory_analysis(Sporadic_PCA, phenoSporadic)

#pdf("LRRK2_traj.pdf")
#pseudo_LRRK2$plot
#dev.off()

pseudo_HC$plot
pseudo_GBA$plot
pseudo_LRRK2$plot
pseudo_Sporadic$plot

pheno_pseudo <- function(phenoHC_, pseudo_){
  phenomine <- phenoHC_
  pseudotime <- rbind(as.data.frame(pseudo_HC$time),as.data.frame(pseudo_$time))
  colnames(pseudotime) <- c('pseudotime')
  phenomine <- phenomine[,c('PATNO','SEX','Time_point','subgroup','age')]
  identical(rownames(phenoHC_), rownames(pseudotime))
  phenoforlme <- merge(phenoHC_, pseudotime, by = 'row.names')
  rownames(phenoforlme) <- phenoforlme$Row.names
  phenoforlme$Row.names <- NULL
  phenoforlme$Time_point <- as.factor(phenoforlme$Time_point)
  phenoforlme$PATNO <- as.factor(phenoforlme$PATNO)
  return(phenoforlme[,c('PATNO','SEX','Time_point','subgroup','age','pseudotime')])
  
}

pheno_pseudo_GBA <- pheno_pseudo(phenoHC_GBA, pseudo_GBA)
pheno_pseudo_GBA$pseudotime <- 1 - pheno_pseudo_GBA$pseudotime

#pseudo_LRRK2$time$`traj$time` <- 1 - pseudo_LRRK2$time$`traj$time`
pheno_pseudo_LRRK2 <- pheno_pseudo(phenoHC_LRRK2, pseudo_LRRK2)
pheno_pseudo_LRRK2$pseudotime <- 1 - pheno_pseudo_LRRK2$pseudotime

#pseudo_Sporadic$time$Lineage1 <- max(pseudo_Sporadic_Sli$time$Lineage1) - pseudo_Sporadic$time$Lineage1
pseudo_Sporadic$time$`traj$time` <- 1 - pseudo_Sporadic$time$`traj$time`
pheno_pseudo_Sporadic <- pheno_pseudo(phenoHC_Sporadic, pseudo_Sporadic)
pheno_pseudo_Sporadic$pseudotime <- 1 - pheno_pseudo_Sporadic$pseudotime

condition <- 'LRRK2'

if (condition == 'GBA') {
  pheno_pseudo = pheno_pseudo_GBA
  combined_df = combined_df_GBA
  dat.regcond = dat.regGBA
} else if (condition == 'LRRK2') {
  pheno_pseudo = pheno_pseudo_LRRK2
  combined_df = combined_df_LRRK2
  dat.regcond = dat.regLRRK2
} else {pheno_pseudo = pheno_pseudo_Sporadic
        combined_df = combined_df_Sporadic
        dat.regcond = dat.regSporadic }

pheno_pseudo$subgroup <- factor(pheno_pseudo$subgroup, levels = c(condition,"Healthy Control"))
combined_df$subgroup <- factor(combined_df$subgroup, levels = c(condition,"Healthy Control"))

plot1 <- ggplot(data = combined_df, aes(x = PC1, y = PC2, color = subgroup)) +
  geom_point() +
  labs(x = "Principal Component 1", y = "Principal Component 2", title = condition) +
  theme_minimal()

plot2 <- ggplot(data = pheno_pseudo, aes(x = pseudotime, y = Time_point, color = subgroup )) +
  geom_jitter() +
  labs(x = "Pseudotime", y = "Time Points") +
  theme_minimal()

plot3 <- ggplot(pheno_pseudo, aes(x = Time_point, y = pseudotime, fill = subgroup)) + 
  geom_boxplot()+
  labs(x = "Pseudotime", y = "Time Points") +
  theme_minimal()

plot4 <- ggplot(data = pheno_pseudo, aes(x = Time_point,y = pseudotime, group = PATNO, color = subgroup)) +
  geom_line() +
  labs(x = "Pseudotime", y = "Time Points") +
  theme_minimal()

grid_arrange_plot <-grid.arrange(plot1,plot2,plot3,plot4, ncol = 2)



###############################################################
#                                                        ######
# Linear Mixed Eff Model                                 ######
#                                                        ######
###############################################################

# match pseudotime to use it as categorical variable 

pseudotime <- rbind(as.data.frame(pseudo_HC$time),as.data.frame(pseudo_GBA$time),as.data.frame(pseudo_LRRK2$time),as.data.frame(pseudo_Sporadic$time))

pseudotime$`traj$time` <- 1 - pseudotime$`traj$time`

pheno_pseudo <- cbind(pheno[,c('subgroup','PATNO')],pseudotime)

pheno_pseudo <-merge(pheno[,c('subgroup','PATNO')], pseudotime, 
                       by = 'row.names', all = TRUE) 

rownames(pheno_pseudo) <- pheno_pseudo$Row.names
pheno_pseudo$Row.names <- NULL

pheno_pseudo$PATNO <- as.character(pheno_pseudo$PATNO)
pheno_pseudo$subgroup <- factor(pheno_pseudo$subgroup, levels = c('Healthy Control','GBA','LRRK2','Sporadic'))
pheno_pseudo$subgroup <- relevel(pheno_pseudo$subgroup, ref = 'Healthy Control')
colnames(pheno_pseudo) <- c('subgroup','ID','pseudotime')


pheno$subgroup <- factor(pheno$subgroup, levels = c('Healthy Control','GBA','LRRK2','Sporadic'))
pheno$subgroup <- relevel(pheno$subgroup, ref = 'Healthy Control')

match_order <- match(colnames(dat.reg), rownames(pheno_pseudo))

pheno_pseudo <- pheno_pseudo[match_order,]

identical(rownames(pheno),colnames(dat.reg))

# substitute t.Time_point with t.Pseudotime and pheno with pheno_pseudo

mixm <- apply(dat.reg, 1, function(y){
  d1 <- data.frame(y = y, t = pheno)
  mixmodel <- lmer(y ~ t.Time_point * t.subgroup + (1|t.PATNO), data=d1)
  comparison <- glht(mixmodel, linfct = mcp(t.subgroup = "Tukey"))
  return(list(model = summary(mixmodel), comparison = summary(comparison)))
})

results_df <- data.frame(
  interaction_effect_GBA_pval = numeric(0), 
  interaction_effect_LRRK2_pval = numeric(0), 
  interaction_effect_Sporadic_pval = numeric(0),
  interaction_effect_GBA_fold_change = numeric(0),
  interaction_effect_LRRK2_fold_change = numeric(0),
  interaction_effect_Sporadic_fold_change = numeric(0),
  comparison_pval = numeric(0)
)

for (i in 1:length(mixm)) { 
  gene_name <- unique(rownames(dat.reg))[i]
  model_summary <- mixm[[i]]
  
  # Extract p-value and fold change for the interaction term in GBA
  interaction_effect_GBA_pval <- coef(summary(model_summary$model))["t.Time_point:t.subgroupGBA", 'Pr(>|t|)']
  interaction_effect_GBA_fold_change <- coef(summary(model_summary$model))["t.Time_point:t.subgroupGBA", 'Estimate']
  
  # Extract p-value and fold change for the interaction term in LRRK2
  interaction_effect_LRRK2_pval <- coef(summary(model_summary$model))["t.Time_point:t.subgroupLRRK2", 'Pr(>|t|)']
  interaction_effect_LRRK2_fold_change <- coef(summary(model_summary$model))["t.Time_point:t.subgroupLRRK2", 'Estimate']
  
  # Extract p-value and fold change for the interaction term in Sporadic
  interaction_effect_Sporadic_pval <- coef(summary(model_summary$model))["t.Time_point:t.subgroupSporadic", 'Pr(>|t|)']
  interaction_effect_Sporadic_fold_change <- coef(summary(model_summary$model))["t.Time_point:t.subgroupSporadic", 'Estimate']
  
  # Make the contrast
  comparison_pval <- model_summary$comparison$test$pvalues[4] # 4 is LRRK2-GBA
  
  # Create a row for the current gene
  row_gene_i <- data.frame(
    interaction_effect_GBA_pval = interaction_effect_GBA_pval,
    interaction_effect_LRRK2_pval = interaction_effect_LRRK2_pval,
    interaction_effect_Sporadic_pval = interaction_effect_Sporadic_pval,
    interaction_effect_GBA_fold_change = interaction_effect_GBA_fold_change,
    interaction_effect_LRRK2_fold_change = interaction_effect_LRRK2_fold_change,
    interaction_effect_Sporadic_fold_change = interaction_effect_Sporadic_fold_change,
    comparison_pval = comparison_pval
  )
  rownames(row_gene_i) <- gene_name
  
  # Append to the results data frame
  results_df <- rbind(results_df, row_gene_i)
}

write.csv(results_df,'Res_Tables/all_in_one_with_comparison_mixmod_Res_TI.csv')

# If use pheno_pseudo 
#write.csv(results_df,'Res_Tables/all_in_one_with_comparison_mixmod_Res_TP.csv')













