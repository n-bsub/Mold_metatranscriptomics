################################################################################
#
#                       Mold growth metatranscriptomics
#                      
#       Balasubrahmaniam et al 2024. Moving Beyond Species: Fungal function 
#       provides novel targets for potential indicators of mold growth in homes.
#
#                       By: Neeraja Balasubrahmaniam
#
################################################################################

# Script name: Target_genes_heatmap_MOLD
# Created on: 13 October 2023
# Last updated: 22 April 2024
# Author: Neeraja Balasubrahmaniam
# Description: This script plots the expression of target genes (100% 
# upregulated & 100% and 85% upregulated target genes, Figure 6) on a heatmap. 
# Inkscape [v.1.3] (https://inkscape.org/) was used to add additional labels 
# for the figure. TMM normalized CPM (Counts Per Million) values that are
# log2 transformed and row mean-centered were used to plot the heatmap.

# ---- Setup -------------------------------------------------------------------

# Install packages
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ComplexHeatmap")

# Load packages
library(ComplexHeatmap)

# Package version
packageVersion("ComplexHeatmap")

# Report package version
# > packageVersion("ComplexHeatmap")
# [1] ‘2.15.1’

# Set working directory
setwd("")

#---- Read data ----------------------------------------------------------------
# TMM-normalized gene expression values for target genes 100% and 100%&85% upregulated genes combined
data=read.table(file="group1&2_tmm.txt", header = FALSE, 
                sep='	');  
# TMM-normalized gene expression values for target genes 85% upregulated genes
data2=read.table(file="group3_tmm.txt", header = FALSE, 
                sep='	');  

# Note: TMM values of all genes present and their fungal swissprot annotations
# are provided in Github within the heatmap folder as a .RDS file.

# These transformations are done for both data and data2 from above
df <- data.frame(data[,-c(1,3)])
df2 <- data.frame(data2[,-c(1,3)])

colnames(df) <-  make.unique(as.character(df[1,]))
rownames(df) <- df[,1]
df <- df[,-1]
df <- df[-1,]

colnames(df2) <-  make.unique(as.character(df2[1,]))
rownames(df2) <- df2[,1]
df2 <- df2[,-1]
df2 <- df2[-1,]

df.matrix = as.matrix(df)
df.numeric = matrix(as.numeric(df.matrix),    # Convert to numeric matrix
                               ncol = ncol(df.matrix))

colnames(df.numeric) <- colnames(df)
rownames(df.numeric) <- rownames(df)

df.numeric = log2(df.numeric+1) # log transformation
# Centering rows
df.numeric = t(scale(t(df.numeric), scale=F))

df2.matrix = as.matrix(df2)
df2.numeric = matrix(as.numeric(df2.matrix),    # Convert to numeric matrix
                    ncol = ncol(df2.matrix))

colnames(df2.numeric) <- colnames(df2)
rownames(df2.numeric) <- rownames(df2)

df2.numeric = log2(df2.numeric+1) # log transformation
# Centering rows
df2.numeric = t(scale(t(df2.numeric), scale=F))

# Color palette for both heatmap plots col_fun and col_fun2
library(circlize)
col_fun = colorRamp2(c(min(df.numeric), median(df.numeric), max(df.numeric)), 
                     c("red","white","darkblue"))
col_fun(seq(-3, 3))

col_fun2 = colorRamp2(c(min(df2.numeric), median(df2.numeric), max(df2.numeric)), 
                     c("red","white","darkblue"))
col_fun2(seq(-3, 3))

# 
final_colnames = as.character(data[1,4:30])
column_labels = structure(final_colnames, names = colnames(df.numeric))

final_colnames2 = as.character(data2[1,3:29])
column_labels2 = structure(final_colnames2, names = colnames(df2.numeric))

#---- Plotting------------------------------------------------------------------


# Plot heatmap for 85% upregulated target genes; Manuscript Fig 6
p1 = Heatmap(df.numeric, name = " ", 
        cluster_columns = FALSE, cluster_rows = FALSE, 
        show_row_dend = FALSE, row_names_side = "left",
        col = col_fun, 
        column_labels = column_labels[colnames(df.numeric)],
        row_names_gp = gpar(fontsize = 16, fontfamily = "sans"),
        column_names_gp = gpar(fontsize = 15, fontfamily = "sans"),
        heatmap_legend_param = list(labels_gp = gpar(fontsize = 12, 
                                                     fontfamily = "sans")),
        )
p1

# Plot heatmap for 85% upregulated target genes; SI Fig S13
p2 = Heatmap(df2.numeric, name = " ", 
             cluster_columns = FALSE, cluster_rows = FALSE, 
             show_row_dend = FALSE, row_names_side = "left",
             col = col_fun2, 
             column_labels = column_labels2[colnames(df2.numeric)],
             row_names_gp = gpar(fontsize = 16, fontfamily = "sans"),
             column_names_gp = gpar(fontsize = 15, fontfamily = "sans"),
             heatmap_legend_param = list(labels_gp = gpar(fontsize = 12, 
                                                          fontfamily = "sans")),
)
p2

#---- To save figure as PDF-----------------------------------------------------
pdf(file = "Group1&2_heatmap.pdf", height = 8, width = 8)
p1   # or p2 if plotting p2
dev.off()
