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

# Script name: assembly_stats_MOLD_plot
# Created on: 6 February 2024
# Last updated: 22 April 2024
# Author: Neeraja Balasubrahmaniam
# Description: This script plots read statistics for the metatranscriptome assembly. 
# Inkscape [v.1.3] (https://inkscape.org/) was used to combine the plot generated
# in this script and in read_stats_MOLD_Plot.R.

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
data=read.table(file="assembly_stats_toPlot.txt", header = TRUE, sep='	');
df=as.data.frame(data)                            # Total Trinity transcripts/contigs = 1983474
                                                  # Total Trinity genes = 1023948

#---- Plotting------------------------------------------------------------------
p2 = ggplot(df,aes(factor(Type, levels = c("KO annotation", "GO mapping", "SP annotation",
                                     "CD-HIT-EST clusters", "Trinity contigs")), 
              Contig.count/(10^6), fill = Type))+ 
  geom_bar(stat="identity", position="identity", color = "black", width = 0.8)+
  
  scale_fill_manual(values = c("#600000", "#d62e2e", "#e39533", "#b4a12d","#4c5a23"),
                    breaks=c('Trinity contigs', 'CD-HIT-EST clusters', 'SP annotation',
                             'GO mapping', 'KO annotation'),
                    labels=c('Trinity contigs', 'CD-HIT-EST clusters', 'SP annotation',
                             'GO mapping', 'KO annotation')) +
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
        axis.line = element_line(colour = "black", linewidth = 1),
        #legend.position = "bottom"
        #panel.border = element_rect(colour = "black", fill=NA, size=1)
  ) +
  labs(x=NULL, y = expression(paste("Contig count (x", 10^6, ")")), fill=NULL) 

p2

#---- Save the plot p2 as a PDF-------------------------------------------------

pdf("Assembly_stats_plot_MOLD.pdf", height = 6, width = 12)
p2
dev.off()
