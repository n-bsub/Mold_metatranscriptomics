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

# Script name: GO_Bubble_plot_MOLD
# Created on: 26 October 2023
# Last updated: 22 April 2024
# Author: Neeraja Balasubrahmaniam
# Description: This script plots the enriched growth GO terms on a bubble plot.

# ---- Setup -------------------------------------------------------------------

# Install packages
install.packages("ggplot2")
install.packages("dplyr")
install.packages("viridis")
install.packages("ggh4x")

# Load packages
library(ggplot2)
library(dplyr)
library(viridis)
library(ggh4x)

# Package versions
packageVersion("ggplot2")
packageVersion("dplyr")
packageVersion("viridis")
packageVersion("ggh4x")

# Reporting package versions
# > packageVersion("ggplot2")
# [1] ‘3.4.3’
# > packageVersion("dplyr")
# [1] ‘1.1.3’
# > packageVersion("viridis")
# [1] ‘0.6.3’
# > packageVersion("ggh4x")
# [1] ‘0.2.8’

#---- Set working directory-----------------------------------------------------
setwd("")

#---- Read GO dataset-----------------------------------------------------------
GO_all=read.table("Bubble_plot_All_GOs_USE.txt",header=T,row.names=NULL,sep="\t")
colnames(GO_all)=c("category","GO","GO_order", "RH", "RHgroup", "-log10FDR", 
                   "Count") # set colnames
# Or read .RDS file instead of above
#GO_all <- readRDS("GO_All.rds")

#---- Plotting------------------------------------------------------------------

p1 = ggplot(GO_all, aes(x=factor(RH, levels = c("100% vs 85%", 
                                           "100% vs 50%",
                                           "85% vs 50%",
                                           "85% vs 100%",
                                           "50% vs 100%",
                                           "50% vs 85%")) , 
                   y=reorder(GO,-GO_order), group = RHgroup, 
                   color = `-log10FDR`, size = Count)) +
  geom_point(alpha = 0.8) +
  scale_size(range = c(3, 20)) +
  xlab(NULL) +
  ylab(NULL) +
  scale_color_viridis(option="viridis", limits = c(1, 15), direction = -1, 
                      oob = scales::squish) +
  theme_bw() +
  theme(axis.text.x = element_text(size=15, color = "black", family = "sans", 
                                   angle = 30, vjust = 1, hjust=0.9),
        axis.text.y = element_text(size=15, color = "black", family = "sans"),
        legend.text = element_text(size=15, color = "black", family = "sans"),
        legend.title = element_text(size=15, color = "black", family = "sans"),
        panel.background = element_blank(),
        panel.grid.major = element_line(linewidth = 0.25),
        #panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
        #panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  

  facet_grid2(~ factor(RHgroup, levels=c('Saturated', 'Elevated', 'Low')), 
              scales = "free_x", space='free',
              strip = strip_themed(
                background_x = list(element_rect(fill = "#cdc1d9"),
                                    element_rect(fill = "#e3dfe7"),
                                    element_rect(fill = "#f2fffb")))) +
  theme(strip.background = element_blank(),
        strip.placement = "inside",
        strip.text = element_text(size=18, face = "bold", color = "black")) +
  scale_x_discrete(breaks = c("100% vs 85%", 
                              "100% vs 50%",
                              "85% vs 50%",
                              "85% vs 100%",
                              "50% vs 100%",
                              "50% vs 85%"))

p1

#---- Save plot as PDF----------------------------------------------------------
pdf(file = "GO_bubble_combo.pdf", width = 14, height = 8); # to save as PDF file

p1

dev.off()


