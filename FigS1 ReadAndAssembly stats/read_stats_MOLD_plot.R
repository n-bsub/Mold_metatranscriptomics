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

# Script name: read_stats_MOLD_plot
# Created on: 6 February 2024
# Last updated: 22 April 2024
# Author: Neeraja Balasubrahmaniam
# Description: This script plots read statistics for the metatranscriptome assembly. 
# Inkscape [v.1.3] (https://inkscape.org/) was used to combine the plot generated
# in this script and in assembly_stats_MOLD_Plot.R.

# ---- Setup -------------------------------------------------------------------

# Install packages
install.packages("ggplot2")

# Load packages
library(ggplot2)

# Package versions
packageVersion("ggplot2")

# Reporting packages versions 
# > packageVersion("ggplot2")
# [1] ‘3.4.3’

# Set working directory
setwd("")

#---- Read data-----------------------------------------------------------------
data=read.table(file="read_stats_toPlot.txt", header = TRUE, sep='	');
df=as.data.frame(data)
df = df %>%
  pivot_longer(!c(Sample, Order), names_to = "type", values_to = "read")



#---- Plotting------------------------------------------------------------------
p1 = ggplot(df,aes(reorder(Sample,-Order), read/(10^7), 
              fill = factor(type, levels = c("Survived.Trimmomatic.filtering", 
                                             "Survived.k.mer.based.filtering", "Initial.read.count"))))+ 
  geom_bar(stat="identity", position="dodge", color = "black", width = 0.7)+
  scale_fill_viridis_d(option = "magma", direction = -1, 
                       breaks=c('Initial.read.count', 
                                'Survived.k.mer.based.filtering', 'Survived.Trimmomatic.filtering'),
                       labels=c("Initial read count","Survived k-mer based filtering",
                                "Survived Trimmomatic filtering")) +
  theme(axis.text.x = element_text(angle = 0, size = 24))+
  coord_flip() +
  theme(axis.text.y   = element_text(size=24, color = "black", family = "sans"),
        axis.text.x   = element_text(size=24, color = "black", family = "sans"),
        axis.title.y  = element_text(size=24, color = "black", family = "sans"),
        axis.title.x  = element_text(size=24, color = "black", family = "sans"),
        legend.text = element_text(size=24, color = "black", family = "sans"),
        legend.title = element_text(size=24, color = "black", family = "sans"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "bottom"
  ) +
  guides(fill=guide_legend(ncol=1))+
  labs(x="Sample", y = expression(paste("Read count (x", 10^7, ")")), fill=NULL) 

p1

#---- Save the plot p1 as a PDF-------------------------------------------------

pdf("Read_stats_plot_MOLD1.pdf", height = 14, width = 8)
p1
dev.off()

