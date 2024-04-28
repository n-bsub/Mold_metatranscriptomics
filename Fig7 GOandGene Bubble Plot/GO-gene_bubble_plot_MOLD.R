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

# Script name: GO-gene_bubble_plot_MOLD
# Created on: 8 February 2024
# Last updated: 23 April 2024
# Author: Neeraja Balasubrahmaniam
# Description: This script plots log2FC and functional categories of target 
# genes in bubble plots. Inkscape [v.1.3] (https://inkscape.org/) was 
# used for further modifications of sizes of plots on the .pdf plot  
# that this script produced.

# ---- Setup -------------------------------------------------------------------
# Install necessary packages
install.packages(tidyverse)
install.packages(patchwork)
install.packages(ggplot2)

# Load packages
library(tidyverse)
library(patchwork)
library(ggplot2)

# Package versions
packageVersion("tidyverse")
packageVersion("patchwork")
packageVersion("ggplot2")

# Reporting packages versions 
# > packageVersion("tidyverse")
# [1] ‘2.0.0’
# > packageVersion("patchwork")
# [1] ‘1.1.2’
# > packageVersion("ggplot2")
# [1] ‘3.4.3’

#---- Code to create bubble plots of log2FC data--------------------------------
# log2FC and FDR-adjusted p values of all DE genes, their fungal annotations
# are provided in a separate .txt file as well as in SI Table S10 

#---- Set working directory to read data
setwd("")

#---- Read data for 100% upregulated genes (here, group1)-----------------------
data1=read.table(file="group1_go_gene_longer.txt", header = TRUE, sep='	');

df1 <- data.frame(data1)

df1 <- data.frame(data1)
df1$log2FC = as.numeric(df1$log2FC)

df1 <- df1[ order(df1$gene, decreasing = TRUE), ]

df1$term.group=factor(df1$term.group,levels = c("Morphological", 
                                                "Secondary metabolism",
                                                "Stress response",
                                                "Mitochondria"))
df1$gene = factor(df1$gene, levels=rev(sort(unique(df1$gene))))


#---- Read data for 100% and 85% upregulated genes (here, group2)---------------
data2=read.table(file="group2_go_gene_longer.txt", header = TRUE, sep='	');

df2 <- data.frame(data2)
df2$log2FC = as.numeric(df2$log2FC)

df2 <- df2[ order(df2$gene, decreasing = TRUE), ]

df2$term.group=factor(df2$term.group,levels = c("Morphological", 
                                                "Secondary metabolism",
                                                "Stress response",
                                                "Mitochondria"))
df2$gene = factor(df2$gene, levels=rev(sort(unique(df2$gene))))

#---- Read data for 85% upregulated genes (here, group3)------------------------
data3=read.table(file="group3_go_gene_longer.txt", header = TRUE, sep='	');

df3 <- data.frame(data3)
df3$log2FC = as.numeric(df3$log2FC)

df3 <- df3[ order(df3$gene, decreasing = TRUE), ]

df3$term.group=factor(df3$term.group,levels = c("Morphological", 
                                                "Secondary metabolism",
                                                "Stress response",
                                                "Mitochondria"))
df3$gene = factor(df3$gene, levels=rev(sort(unique(df3$gene))))

#---- Plotting group1 as bubble plot p1-----------------------------------------
p1 = ggplot(df1, aes(x=term.group , y=gene)) +
  geom_point(aes(size = log2FC, fill = term.group), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(0.000001, 1000), range = c(8,30), 
                        breaks = c(10,100,500,750)) +   
  labs( x= "", y = "", size = "log2FC", fill = "")  + 
  theme_bw() +
  theme(axis.text.x = element_text(size=30, color = "black", family = "sans", 
                                   angle = 30, vjust = 1, hjust=0.9),
        axis.text.y = element_text(size=30, color = "black", family = "sans"),
        legend.text = element_text(size=30, color = "black", family = "sans"),
        legend.title = element_text(size=30, color = "black", family = "sans"),
        panel.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25),
        panel.grid.minor = element_line(linewidth = 0.25),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  guides(fill = guide_legend(override.aes = list(size = 10)))

#---- Plotting group2 as bubble plot p2-----------------------------------------
p2 = ggplot(df2, aes(x=term.group , y=gene)) +
  geom_point(aes(size = log2FC, fill = term.group), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(0.000001, 1000), range = c(8,30), 
                        breaks = c(10,100,500,750)) +   
  labs( x= "", y = "", size = "log2FC", fill = "")  + 
  theme_bw() +
  theme(axis.text.x = element_text(size=30, color = "black", family = "sans", 
                                   angle = 30, vjust = 1, hjust=0.9),
        axis.text.y = element_text(size=30, color = "black", family = "sans"),
        legend.text = element_text(size=30, color = "black", family = "sans"),
        legend.title = element_text(size=30, color = "black", family = "sans"),
        panel.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25),
        panel.grid.minor = element_line(linewidth = 0.25),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  guides(fill = guide_legend(override.aes = list(size = 10)))

#---- Plotting group3 as bubble plot p3-----------------------------------------
p3 = ggplot(df3, aes(x=term.group , y=gene)) +
  geom_point(aes(size = log2FC, fill = term.group), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(0.000001, 1000), range = c(8,30), 
                        breaks = c(10,100,500,750)) +   
  labs( x= "", y = "", size = "log2FC", fill = "")  + 
  theme_bw() +
  theme(axis.text.x = element_text(size=30, color = "black", family = "sans", 
                                   angle = 30, vjust = 1, hjust=0.9),
        axis.text.y = element_text(size=30, color = "black", family = "sans"),
        legend.text = element_text(size=30, color = "black", family = "sans"),
        legend.title = element_text(size=30, color = "black", family = "sans"),
        panel.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25),
        panel.grid.minor = element_line(linewidth = 0.25),
        axis.line = element_line(colour = "black"),
        legend.position = "right") +
  guides(fill = guide_legend(override.aes = list(size = 10)))

#------ Use patchwork to combine p1, p2, p3 and plot as p4----------------------
p4 = p1 + p2 + p3 + plot_layout(guides = 'collect')

p4

#---- Save combined plot p4 as PDF----------------------------------------------
pdf(file = "GO-gene-bubble.pdf", height = 14, width = 20)
p4
dev.off()