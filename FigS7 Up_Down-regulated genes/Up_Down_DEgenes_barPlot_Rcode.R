################################################################################
#
#                       Mold growth metatranscriptomics
#                      
#       Balasubrahmaniam et al 2024. Moving Beyond Species: Fungal function 
#       provides novel targets for potential indicators of mold growth in homes.
#
################################################################################

# Script name: Up_Down_DEgenes_barPlot_Rcode
# Created on: 12 November 2023
# Last updated: 28 April 2024
# Author: Neeraja Balasubrahmaniam
# Description: This script plots a bar plot for genes that are upregulated and 
# downregulated in each ERH pairwise comparison, including fungal annotated genes.
# Inkscape was used to finalize the figure and add additional labels.

# ---- Setup -------------------------------------------------------------------
# Install necessary packages
install.packages("ggplot2")

# Load packages
library(ggplot2)

# Package versions
packageVersion("ggplot2")

# Reporting packages versions 
# > packageVersion("ggplot2")
# [1] ‘3.4.3’

# Set working directory to read data
setwd("")

# Read data for number of DE genes
# All DE genes are reported in the tab separated file: All_DE_genes_log2fc_tmm_annots.txt 
# present within the directory
data=read.table(file="Data_for_2sided_barPlot_all&fungal.txt", header = TRUE, sep='	');
df=as.data.frame(data)

DE = ifelse(df$up_down == "up", df$DE.genes, -df$DE.genes)

# Plotting
p1 = ggplot(df,aes(x=reorder(RH,order), y = DE, fill=all_or_fungal))+ 
  geom_bar(stat="identity", position="dodge", width = 0.5)+
  scale_y_continuous(breaks = pretty(DE, n = 6)) +
  theme(axis.text.x = element_text(angle = 0, size = 20))+
  theme(axis.text.y   = element_text(size=20, color = "black", family = "sans"),
        axis.text.x   = element_text(size=20, color = "black", family = "sans"),
        axis.title.y  = element_text(size=20, color = "black", family = "sans"),
        axis.title.x  = element_text(size=20, color = "black", family = "sans"),
        legend.text = element_text(size=20, color = "black", family = "sans"),
        legend.title = element_text(size=20, color = "black", family = "sans"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  geom_abline(aes(intercept = 0, slope = 0)) +
  labs(x="RH comparison", y="Number of differentially expressed genes", fill=NULL) 

p1

# Save as .pdf figure
pdf("Up_Down_DEgenes_barPlot_all&fungalgenes_1.pdf", height = 8, width = 10)

p1

dev.off()

