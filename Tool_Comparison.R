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

###############################################################
#                                                        ######
# import data                                            ######
#                                                        ######
###############################################################
setwd('C:/Users/lukfa/Desktop/Thesis/Res_Tables')

dat.reg <- readRDS("dat_reg.RDS")
pheno <- readRDS("pheno.RDS")
pheno$Time_point <- paste0('t_',pheno$Time_point)

pca_result <- prcomp(t(dat.reg), scale. = TRUE)
p1 <- autoplot(pca_result, data = pheno, colour = 'SEX', label = F, label.size = 3)
p2 <- autoplot(pca_result, data = pheno, colour = 'subgroup', label = F, label.size = 3)
p3 <- autoplot(pca_result, data = pheno, colour = 'Day', label = F, label.size = 3)
p4 <- autoplot(pca_result, data = pheno, colour = 'age', label = F, label.size = 3)

grid_arrange_plot <-grid.arrange(p1,p2,p3,p4, ncol = 2)

###############################################################
#                                                        ######
# Trajectory Analysis                                    ######
#                                                        ######
###############################################################

phenoHC<-pheno[pheno$subgroup=="GBA",]

phenoHC <- phenoHC[!rownames(phenoHC) == "RNAB_PPMI56680_3169_da25_v1", ] # GBA outlier
#phenoHC <- phenoHC[!rownames(phenoHC) == "RNAB_PPMI3452_7426_da25_v1", ] # HC outlier
dat.regHC<-dat.reg[,rownames(phenoHC)]
pca_result <- prcomp(t(dat.regHC), scale. = TRUE)
autoplot(pca_result, data = phenoHC, colour = 'Day', label = F, label.size = 3) +  theme_minimal()

## Scorpius 

# Trajectory 
rd1 <- data.frame(pca_result$x[,1:2])
traj <- SCORPIUS::infer_trajectory(rd1)
draw_trajectory_plot(rd1, as.factor(phenoHC$Day), traj$path)

#Pseudotime 
scorpius_pseudotime <- as.data.frame(traj$time)
colnames(scorpius_pseudotime) <- c('scorpius')

# Gene Selection 
gimp <- gene_importances(
  t(dat.regHC), 
  traj$time, 
  num_permutations = 100,
  ) 
gimp$qvalue <- p.adjust(gimp$pvalue, method = 'fdr')
gene_sel_scorpius <- gimp$gene[gimp$pvalue < 0.05]

# Tscan 

#Trajectory 
lpsmclust <- exprmclust(dat.regHC,clusternum	= 3)
lpsorder <- TSCANorder(lpsmclust)
lpsorder <- TSCANorder(lpsmclust)
plotmclust(lpsmclust,show_cell_names = FALSE)

#Pseudotime 
Tscan_pseudotime <- as.data.frame(TSCANorder(lpsmclust,flip=TRUE,orderonly=FALSE))
Tscan_pseudotime  <- subset(Tscan_pseudotime , select = c('Pseudotime'))
Tscan_pseudotime_inv  <- data.frame(lapply(Tscan_pseudotime, rev))
rownames(Tscan_pseudotime_inv) <- rownames(Tscan_pseudotime)

#Gene Selection
diffval <- difftest(dat.regHC,lpsorder)
diffval$qvalue <- p.adjust(diffval$pval, method = 'bonferroni')
gene_sel_tscan <- row.names(diffval[diffval$qvalue < 0.05,])

# Slingshot 

#Trajectory
sce <- SingleCellExperiment(dat.regHC)
reducedDims(sce) <- SimpleList(PCA = as.matrix(rd1))

kmeans_result <- kmeans(rd1, centers = 3)
cl1 <- kmeans_result$cluster
#cl1 <- Mclust(rd1)$classification

colData(sce)$GMM <- cl1
sce <- slingshot(sce, clusterLabels = 'GMM', reducedDim = 'PCA')
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(sce)$PCA, col = plotcol, pch=16, asp = 1,xlim = c(-200, 200))
lines(SlingshotDataSet(sce), lwd=2, col='black',xlim = c(-200, 200))
plot(reducedDims(sce)$PCA, col = brewer.pal(9,'Set1')[sce$GMM], pch=16, asp = 1,xlim = c(-200, 200))
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

#Pseudotime
slingshot_pseudotime <- as.data.frame(slingPseudotime(sce))

#Gene Selection 

min_val <- min(dat.regHC)
dat.regHC_positive <- dat.regHC + abs(min_val)

sce <- fitGAM(dat.regHC_positive,sds = SlingshotDataSet(sce))

ATres <- associationTest(sce)
ATres <- na.omit(ATres)
ATres$qvalue <- p.adjust(ATres$pvalue, method = 'fdr')
gene_sel_slingshot <- row.names(ATres)[ATres$qvalue < 0.05]


#COMPARISON    

joined_pseudotime1 <- merge(slingshot_pseudotime, Tscan_pseudotime_inv, by = "row.names")
joined_pseudotime2 <- merge(scorpius_pseudotime, Tscan_pseudotime_inv, by = "row.names")
joined_pseudotime <- merge(joined_pseudotime1, joined_pseudotime2, by = "Row.names")
rownames(joined_pseudotime) <- joined_pseudotime$Row.names

joined_pseudotime  <- subset(joined_pseudotime , select = c('Lineage1','Pseudotime.x','scorpius'))
colnames(joined_pseudotime) <- c('Slingshot','TSCAN','Scorpius')


# Create ggplot for Slingshot vs TSCAN
plot1 <- ggplot(joined_pseudotime, aes(x = Slingshot, y = TSCAN)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Slingshot vs TSCAN", x = "Slingshot", y = "TSCAN")

# Create ggplot for Slingshot vs Scorpius
plot2 <- ggplot(joined_pseudotime, aes(x = Slingshot, y = Scorpius)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "Slingshot vs Scorpius", x = "Slingshot", y = "Scorpius")

# Create ggplot for TSCAN vs Scorpius
plot3 <- ggplot(joined_pseudotime, aes(x = TSCAN, y = Scorpius)) +
  geom_point() +
  geom_smooth(method = "lm", se = TRUE) +
  labs(title = "TSCAN vs Scorpius", x = "TSCAN", y = "Scorpius")

grid_arrange_plot <-grid.arrange(plot1,plot2,plot3, ncol = 3)

# Calculate correlation matrix
cor_matrix <- cor(joined_pseudotime[, c("Slingshot", "TSCAN", "Scorpius")])
corrplot(cor_matrix, method = "number")

###############################################################
#                                                        ######
# GAM + bonferroni correction                            ######
#                                                        ######
###############################################################

t <- as.numeric(traj$time)
Y <- data.frame(dat.regHC)

gam.P <- apply(Y,1,function(z){
  d1 <- data.frame(z=z, t=t)
  tmp <- gam(z ~ s(t, bs = "cr", k = 3), data=d1)
  P <- summary(tmp)
  P[4][[1]][1]
})

gam.results<-data.frame(gam.P)
rownames(gam.results)<-names(gam.P)
p.bonf<-p.adjust(gam.results[,1], "bonferroni")
gam.results<-cbind(gam.results,p.bonf)
gam.results.bonf<-gam.results[which(gam.results$p.bonf< 0.005),]
gam.results.bonf<-gam.results.bonf[order(gam.results.bonf$gam.P),] 

topgenes <- rownames(gam.results.bonf)

###############################################################
#                                                        ######
# Plotting regressions                                   ######
#                                                        ######
###############################################################

# Cubic vs local regresion

layout_matrix <- matrix(c(1, 2, 3, 4,5,6), nrow = 2, byrow = TRUE)
layout(layout_matrix)

for (i in head(row.names(diffval)[diffval$qval < 0.05])) {
  row <- Y[i,]
  
  d1 <- data.frame(z = as.numeric(row),
                   t = as.numeric(t))
  
  cubic_regression_spline <- gam(z ~ s(t, bs = "cr", k = 3), data = d1)
  local_regression_smooth <- gam(z ~ lo(t), data = d1)
  
  plot(d1$t, d1$z, xlab = "Time", ylab = "Response")
  curve(predict(local_regression_smooth, newdata = data.frame(t = x)), col = "red", lwd = 2, add = TRUE)
  curve(predict(cubic_regression_spline, newdata = data.frame(t = x)), col = "blue", lwd = 2, add = TRUE)
  
}

plot_list <- list() 

for (i in head(row.names(diffval)[diffval$qval < 0.05])) {
  plot_object <- singlegeneplot(dat.regHC[i,], TSCANorder(lpsmclust, flip = TRUE, orderonly = FALSE))
  
  grob_object <- ggplotGrob(plot_object)
  
  plot_list[[i]] <- grob_objects
}

grid.arrange(grobs = plot_list, ncol = 2)

layout(1)


row <- Y[100,]

d1 <- data.frame(z = as.numeric(row),
                 t = as.numeric(t))

cubic_regression_spline <- gam(z ~ s(t, bs = "cr", k = 3), data = d1)
local_regression_smooth <- gam(z ~ lo(t), data = d1)

# Plot the first two curves

plot(d1$t, d1$z, xlab = "Time", ylab = "Response")
curve(predict(local_regression_smooth, newdata = data.frame(t = x)), col = "red", lwd = 2, add = TRUE)
curve(predict(cubic_regression_spline, newdata = data.frame(t = x)), col = "blue", lwd = 2, add = TRUE)
singlegeneplot(dat.regHC[100,], TSCANorder(lpsmclust, flip = TRUE, orderonly = FALSE))












