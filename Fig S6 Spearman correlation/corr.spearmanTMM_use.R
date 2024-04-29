################################################################################
#
#                       Mold growth metatranscriptomics
#                      
#       Balasubrahmaniam et al 2024. Moving Beyond Species: Fungal function 
#       provides novel targets for potential indicators of mold growth in homes.
#
################################################################################

# Script name: corr.spearmanTMM_use
# Created on: 21 November 2023
# Last updated: 28 April 2024
# Author: Neeraja Balasubrahmaniam
# Description: This script plots a correlation plot for all genes that are 
# differentially expressed (log2FC >= 2, FDR-adjusted p <= 0.001).

# ---- Setup -------------------------------------------------------------------
# Install necessary packages
install.packages("corrplot")

# Load packages
library(corrplot)

# Package versions
packageVersion("corrplot")

# Reporting packages versions 
# > packageVersion("corrplot")
# [1] ‘0.92’

# Set working directory to read data
setwd("")

# Read gene expression matrix: TMM normalized expression values for all genes
# differentially expressed with log2FC >= 2, FDR-adjusted p <= 0.001

data = read.table("diffExpr.P0.001_C2.matrix", header=T, com='', row.names=1, check.names=F, sep='\t')
# or read .RDS file commented below
# data = readRDS("diffExpr.P0.001_C2.matrix.rds")

data = log((data+1)) #log transformation

data = as.matrix(data) # convert to matrix

# Find distance between samples, used for hierarchical clustering in the correlation plot.
sample_dist = dist(t(data), method='euclidean') #Euclidean distances for hierarchical clustering
hc_samples = hclust(sample_dist, method='complete')

# Perform correlation
sample_cor = cor(data, method='spearman', use='pairwise.complete.obs') 
sample_cor <- sample_cor[hc_samples$labels[c(hc_samples$order)], 
                         hc_samples$labels[c(hc_samples$order)]] #row and col names


# Compute p-values based on correlation
cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

# matrix of the p-value of the correlation
p.mat <- cor.mtest(data)
p.mat = p.mat[hc_samples$labels[c(hc_samples$order)], hc_samples$labels[c(hc_samples$order)]]

# Adjust p-values for multiple comparisons based on FDR
p_adj = matrix(p.adjust(as.vector(as.matrix(p.mat)), method='fdr'),ncol=27)

colnames(p_adj) = colnames(p.mat) # set row names
rownames(p_adj) = rownames(p.mat) # set column names

# Plotting
corrplot(sample_cor, method="circle", diag=TRUE, outline = FALSE, 
         col=(COL2("RdBu", 200)), order = 'hclust', hclust.method = 'complete',
         cl.lim = c(0, 1), is.corr = TRUE, addrect = 3,
         tl.col="black", tl.srt=90, tl.cex = 1.5, cl.cex = 1.5,
         p.mat = p_adj, sig.level = 0.05, insig = "blank")

# Save to .pdf file
pdf(file = "Spearmancorrplot_allTMM.pdf", height = 8, width = 8)

corrplot(sample_cor, method="circle", diag=TRUE, outline = FALSE, 
         col=(COL2("RdBu", 200)), order = 'hclust', hclust.method = 'complete',
         cl.lim = c(0, 1), is.corr = TRUE, addrect = 3,
         tl.col="black", tl.srt=90, tl.cex = 1.5, cl.cex = 1.5,
         p.mat = p_adj, sig.level = 0.05, insig = "blank")

dev.off()
