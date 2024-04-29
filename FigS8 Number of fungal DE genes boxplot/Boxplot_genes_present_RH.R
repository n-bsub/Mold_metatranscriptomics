################################################################################
#
#                       Mold growth metatranscriptomics
#                      
#       Balasubrahmaniam et al 2024. Moving Beyond Species: Fungal function 
#       provides novel targets for potential indicators of mold growth in homes.
#
################################################################################

# Script name: Boxplot_genes_present_RH
# Created on: 24 February 2024
# Last updated: 26 April 2024
# Author: Neeraja Balasubrahmaniam
# Description: This script plots boxplots for the number of genes with a fungal
# annotation present across samples in each ERH condition (50%, 85% and 100%).
# All genes present in each sample are provided separately along with 
# fungal annotations and raw counts in 3 .xlsx files (50%, 85% and 100% RH)


# ---- Setup -------------------------------------------------------------------
# Install necessary packages
install.packages("tidyverse")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("scales") # to access break formatting functions


# Load packages
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(scales) # to access break formatting functions


# Package versions
packageVersion("tidyverse")
packageVersion("dplyr")
packageVersion("ggplot2")
packageVersion("ggpubr")
packageVersion("scales")


# Reporting packages versions 
# > packageVersion("tidyverse")
# [1] ‘2.0.0’
# > packageVersion("dplyr")
# [1] ‘1.1.3’
# > packageVersion("ggplot2")
# [1] ‘3.4.3’
# > packageVersion("ggpubr")
# [1] ‘0.6.0’
# > packageVersion("scales")
# [1] ‘1.2.1’


# Set working directory to read data
setwd("") # path of directory in ""

# Read data for number of genes present in each sample
data=read.table(file="Fungal_gene_number_RH.txt", header = TRUE, sep='	');

df=as.data.frame(data)
colnames(df2) = c("Site", "ERH", "Conc.")

df$ERH = factor(x = df$ERH, 
                levels = c("50%", "85%", 
                           "100%"))

# Plotting
p1 = df %>%
  ggplot(mapping = aes(x = ERH, y = 	
                         Number.of.fungal.annotated.genes, fill = `ERH`)) +
  stat_boxplot(geom = "errorbar", width = 0.2) + 
  geom_boxplot(colour = "black", linewidth = 0.8) +
  scale_y_continuous(trans='log10') +
  labs(x = "ERH", y = "Number of fungal annotated genes") +
  theme(axis.title.y = element_text(size = 14, color = "black", margin=margin(r=12))) +
  theme(axis.title.x = element_text(size = 14, color = "black", margin=margin(t=12))) +
  theme(axis.text.y = element_text(size = 14, color = "black")) +
  theme(axis.text.x = element_text(size = 14, color = "black")) +
  theme(axis.ticks.y = element_line(linewidth = 0.7, linetype = "solid", color = "black")) +
  theme(axis.ticks.x = element_line(linewidth = 0.7, linetype = "solid", color = "black")) +
  theme(axis.ticks.length.y = unit(0.1, "cm")) +
  theme(axis.ticks.length.x = unit(0.1, "cm")) +
  theme(axis.line.y = element_line(linewidth = 0.7, linetype = "solid", color = "black")) +
  theme(axis.line.x = element_line(linewidth = 0.7, linetype = "solid", color = "black")) +
  theme(legend.title = element_text(size = 12, hjust = 0, color = "black")) +
  theme(legend.text = element_text(size = 12, color="black")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, linewidth = 0.7, linetype = "solid", color = "black")) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.major.x = element_blank())
p1

# Save plot as .pdf
pdf("Boxplot_fungal_gene_number.pdf", height = 5, width = 5)
p1
dev.off()

