################################################################################
#
#                       Mold growth metatranscriptomics
#                      
#       Balasubrahmaniam et al 2024. Moving Beyond Species: Fungal function 
#       provides novel targets for potential indicators of mold growth in homes.
#
################################################################################

# Script name: Rcode_PCoA_Relative_Taxa_Aitchison
# Created on: 24 February 2024
# Last updated: 26 April 2024
# Author: Neeraja Balasubrahmaniam and Jon C. King
# Description: This script performs PCA on gene expression and and PCoA on 
# relative abundance of taxa (ASVs, species, genus) data using Aitchison distances.

# ---- Setup -------------------------------------------------------------------
# Install necessary packages
install.packages("tidyverse")
install.packages("vegan")
install.packages("glue")
install.packages("ape")
install.packages("patchwork")
install.packages("dplyr")
install.packages("ggplot2")

# Load packages
library(tidyverse)
library(vegan)
library(glue)
library(ape)
library(patchwork)
library(dplyr)
library(ggplot2)

# Package versions
packageVersion("tidyverse")
packageVersion("vegan")
packageVersion("glue")
packageVersion("ape")
packageVersion("patchwork")
packageVersion("dplyr")
packageVersion("ggplot2")

# Reporting packages versions 
# > packageVersion("tidyverse")
# [1] ‘2.0.0’
# > packageVersion("vegan")
# [1] ‘2.6.4’
# > packageVersion("glue")
# [1] ‘1.6.2’
# > packageVersion("ape")
# [1] ‘5.7.1’
# > packageVersion("patchwork")
# [1] ‘1.1.2’
# > packageVersion("dplyr")
# [1] ‘1.1.3’
# > packageVersion("ggplot2")
# [1] ‘3.4.3’

# Set working directory to read data
setwd("") # path of directory in ""

# import data
data <- readRDS("tblREL.rds")
data = tblREL

# prepare data
data_ASV <- data
data_Gen <- data %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum))
data_Spec <- data %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), sum))

d_list <- list(data_ASV, data_Spec, data_Gen) # combine ASV, species, genus into list
d_list2 <- list() # initialize
d_list3 <- list() # initialize
for (i in seq(1:3)) {
  d_list2[[i]] <- d_list[[i]] %>%
    select(where(is.numeric)) %>%
    t() %>%
    data.frame() %>%
    mutate(ERH = case_when(str_detect(rownames(.), "50") ~ "50%",
                           str_detect(rownames(.), "85") ~ "85%",
                           str_detect(rownames(.), "100") ~ "100%"),
           Site = case_when(str_detect(rownames(.), "CO") ~ "CO",
                             str_detect(rownames(.), "S1") ~ "OH.1",
                             str_detect(rownames(.), "S2") ~ "OH.2",
                             str_detect(rownames(.), "JB") ~ "OH.3",
                             str_detect(rownames(.), "PA") ~ "PA",
                             str_detect(rownames(.), "KS") ~ "KS",
                             str_detect(rownames(.), "SF") ~ "CA",
                             str_detect(rownames(.), "WA") ~ "WA",
                             str_detect(rownames(.), "MI") ~ "MI")) %>%
    relocate(c("Site","ERH"), .before = everything())

clr <- decostand(d_list2[[i]][,-(1:2)], method = "clr",
                 MARGIN = 2, pseudocount = 1);          # center log ratio

# center log ratio transformed data with Euclidean distances give Aitchison distances

d_list3[[i]] <- vegdist(clr, method = "euclidean");  # distance matrix

}

names(d_list3) <- c("ASVs", "Species", "Genus") 


# PERMANOVA
for (i in seq(1:3)) {
  print(names(d_list3)[i]);
  #list_dna_rel_clr_eucl[[i]] <- vegdist(clr, method = "euclidean");  # Euclidean
  #list_dna_rel_rclr_eucl[[i]] <- vegdist(rclr, method = "euclidean")
  print(adonis2(d_list3[[i]] ~ ERH, d_list2[[i]], permutations = 10000, 
                method = "euclidean")) #method = "aitchison", pseudocount=1
}

# Perform PCoA and plot (grouped by ERH)

letters1 = c("A", "C", "E") # initialize for plot labels

d_list4 <- list() # initialize
for (i in seq(1:3)) {
  df <- pcoa(d_list3[[i]])$vectors %>%  # PCoA
    data.frame() %>%
    bind_cols(select(d_list2[[i]], ERH)) %>%
    relocate(c("ERH"), .before = everything())
  df$ERH <- factor(df$ERH, levels = c("100%", "85%", "50%")) 
  d_list4[[i]] <- ggplot(df, aes(x = Axis.1,                       # plotting
                   y = Axis.2, color = ERH, shape = ERH)) +
    geom_point(size=3) +
    scale_shape_manual(values=c(15, 17, 19)) +
    scale_color_manual(values=c('#b91cd7','#009cff', '#00f2aa'))+
    stat_ellipse(geom = "polygon",
                 aes(fill = ERH), 
                 alpha = 0, linewidth = 1) +
    labs(title = paste0("(", letters1[i],") ", "Taxa (", names(d_list3)[i], ")"),      
         x = paste0("PC1 (", 
                    format(round(pcoa(d_list3[[i]])$values$Relative_eig[1],3)*100, nsmall=1),"%)"),
         y = paste0("PC2 (", 
                    format(round(pcoa(d_list3[[i]])$values$Relative_eig[2],3)*100, nsmall=1),"%)")) +
    xlim(-42, 61) +    #setting plot X axis limits
    ylim(-56, 46) +    #setting plot Y axis limits
    theme(axis.text.y   = element_text(size=18, color = "black", family = "sans"),
          axis.text.x   = element_text(size=18, color = "black", family = "sans"),
          axis.title.y  = element_text(size=18, color = "black", family = "sans"),
          axis.title.x  = element_text(size=18, color = "black", family = "sans"),
          legend.text = element_text(size=18, color = "black", family = "sans"),
          legend.title = element_text(size=18, color = "black", family = "sans"),
          plot.title = element_text(size=18, color = "black", family = "sans", hjust = 0),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    ) +
    coord_fixed(ratio = 1)
}
lapply(d_list4, function(x) x)  # show plots; list of 3 plots- ASVs, species, genus

#---- Set working directory to read gene expression raw counts data
setwd("")

# ---- Load gene expression data for PCA
data = read.table("RSEM.gene.counts.matrix",               # gene counts table
                   header=T, com='', row.names=1, check.names=F, sep='\t')
# Or load .rds file instead of gene counts table from previous line
#data = readRDS("RSEM.gene.counts.matrix.rds")

data = data[rowSums(data)>=10,] # remove genes which when summed across samples < 10
cs = colSums(data) 
data = t( t(data)/cs) * 1e6; # Counts per Million (CPM)
data = log2(data+1) # log transformation of CPM values

primary_data_PCA = as.matrix(primary_data_PCA) # convert dataframe to matrix

prin_comp_data = primary_data_PCA

#---- Perform PCA using prcomp

pca = prcomp(prin_comp_data, center = FALSE, scale. = FALSE)

pc_pct_variance = (pca$sdev^2)/sum(pca$sdev^2) # these are %variance on the axes

PCA.loadings=pca$x
PCA.scores = pca$rotation # these are the points on the graph (use only PC1, PC2)

#---- Plotting------------------------------------------------------------------

#---- Get a samples-to-RH file in the same order as colnames of prin_comp_data (samples)
Samples_to_RH=read.table(file="sampletoRH&Site_gene.txt", header = FALSE, sep='	'); 
Samples_to_RH = data.frame(Samples_to_RH)

#---- make a combined dataframe of the columns RH,PC1,PC2
PCA.df1 <- as.data.frame(PCA.scores[,c(1,2)]) # only PC1 vs PC2
df1 = cbind(Samples_to_RH,PCA.df1) # combine RH, Site, RH_level, PC1, PC2
colnames(df1) <- c("SampleID", "ERH", "RH.level", "Site", "PC1", "PC2") #colnames

# Format label % variance values for plot
pretty_pe <- format(round(pc_pct_variance*100, digits =1), nsmall=1, trim=TRUE)

labels <- c(glue("PC1 ({pretty_pe[1]}%)"),
            glue("PC2 ({pretty_pe[2]}%)"))

# Factor ERH to 100%, 85%, 50% order for plot
df1$ERH <- factor(df1$ERH, levels = c("100%", "85%", "50%"))

#---- Plotting PCA grouped by ERH-------------------------------------
letters2 = c("D", "E", "F")  # change based on plot order preference

d_list5 <- list()
for (j in seq(1:3)) {
  d_list5[[j]] = ggplot(df1, aes(x = PC1, 
                                 y = PC2, color = ERH, shape = ERH)) +
    geom_point(size=3) +
    scale_shape_manual(values=c(15, 17, 19)) +
    scale_color_manual(values=c('#b91cd7','#009cff', '#00f2aa'))+
    stat_ellipse(geom = "polygon",
                 aes(fill = ERH), 
                 alpha = 0, linewidth = 1) +
    labs(title = paste0("(", letters2[j],") ", "Gene Expression"))+
    xlab(labels[1]) +  
    ylab(labels[2]) +    
    xlim(-0.414,0.4) +    #setting plot X axis limits
    ylim(-0.4,0.4) +    #setting plot Y axis limits
    theme(axis.text.y   = element_text(size=18, color = "black", family = "sans"),
          axis.text.x   = element_text(size=18, color = "black", family = "sans"),
          axis.title.y  = element_text(size=18, color = "black", family = "sans"),
          axis.title.x  = element_text(size=18, color = "black", family = "sans"),
          legend.text = element_text(size=18, color = "black", family = "sans"),
          legend.title = element_text(size=18, color = "black", family = "sans"),
          plot.title = element_text(size=18, color = "black", family = "sans", hjust = 0),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    ) +
    coord_fixed(ratio = 1)
}

# plot gene expression 3 times to compare each with 3 of taxa PCoA plots
lapply(d_list5, function(x) x) 


# PCoA-grouped by Site

letters1 = c("B", "D", "F") # chnage based on prefernce of plot order

my_cols = c("#FFD92F", "#FC8D62", "#8DA0CB","#E78AC3","#A6D854","#66C2A5"
            ,"#00bb00", "#B3B3B3",  "#E5C494")   # color palette for plot

d_list6 <- list()
for (i in seq(1:3)) {
  df <- pcoa(d_list3[[i]])$vectors %>%
    data.frame() %>%
    bind_cols(select(d_list2[[i]], Site, ERH)) %>%
    relocate(c("Site","ERH"), .before = everything())
  df$ERH <- factor(df$ERH, levels = c("100%", "85%", "50%"))
  d_list6[[i]] <- ggplot(df, aes(x = Axis.1, 
                   y = Axis.2, color = Site, shape = ERH)) +
    geom_point(size=3) +
    scale_shape_manual(values=c(15, 17, 19)) +
    scale_color_manual(values=my_cols)+
    labs(title = paste0("(", letters1[i],") ", "Taxa (", names(d_list3)[i], ")"),      
         x = paste0("PC1 (", 
                    format(round(pcoa(d_list3[[i]])$values$Relative_eig[1],3)*100, nsmall=1),"%)"),
         y = paste0("PC2 (", 
                    format(round(pcoa(d_list3[[i]])$values$Relative_eig[2],3)*100, nsmall=1),"%)")) +
    xlim(-42, 61) +    #setting plot X axis limits-use for combined plot
    ylim(-56, 46) + #setting plot Y axis limits-use for combined plot
    #xlim(-30,50) +    #older limits if plotting this as standalone
    #ylim(-46, 28) +    #older limits if plotting this as standalone
    theme(axis.text.y   = element_text(size=18, color = "black", family = "sans"),
          axis.text.x   = element_text(size=18, color = "black", family = "sans"),
          axis.title.y  = element_text(size=18, color = "black", family = "sans"),
          axis.title.x  = element_text(size=18, color = "black", family = "sans"),
          legend.text = element_text(size=18, color = "black", family = "sans"),
          legend.title = element_text(size=18, color = "black", family = "sans"),
          plot.title = element_text(size=18, color = "black", family = "sans", hjust = 0),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    ) +
    coord_fixed(ratio = 1)
}
lapply(d_list6, function(x) x)


#---- Plotting PCA grouped by Site-------------------------------------
letters2 = c("A", "C", "E")
d_list7 <- list()
for (j in seq(1:3)) {
  d_list7[[j]] = ggplot(df1, aes(x = PC1, 
                                 y = PC2, color = Site, shape = ERH)) +
    geom_point(size=3) +
    scale_shape_manual(values=c(15, 17, 19)) +
    scale_color_manual(values=my_cols)+
    labs(title = paste0("(", letters2[j],") ", "Gene Expression"))+
    xlab(labels[1]) +  
    ylab(labels[2]) +    
    xlim(-0.414,0.4) +    #setting plot X axis limits
    ylim(-0.4,0.4) +    #setting plot Y axis limits
    theme(axis.text.y   = element_text(size=18, color = "black", family = "sans"),
          axis.text.x   = element_text(size=18, color = "black", family = "sans"),
          axis.title.y  = element_text(size=18, color = "black", family = "sans"),
          axis.title.x  = element_text(size=18, color = "black", family = "sans"),
          legend.text = element_text(size=18, color = "black", family = "sans"),
          legend.title = element_text(size=18, color = "black", family = "sans"),
          plot.title = element_text(size=18, color = "black", family = "sans", hjust = 0),
          panel.background = element_blank(),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=1)
    ) +
    coord_fixed(ratio = 1)
}
lapply(d_list7, function(x) x)

# To plot figure S3
p6 = ggdraw() +
  draw_plot(d_list4[[1]], x = 0, y = 0.67, width = .5, height = .33) +
  
  draw_plot(d_list4[[2]], x = 0, y = 0.33, width = .5, height = 0.33) +
  
  draw_plot(d_list4[[3]], x = 0, y = 0, width = .5, height = 0.33) +
  
  draw_plot(d_list6[[1]], x = .5, y = 0.67, width = .5, height = .33) +
  
  draw_plot(d_list6[[2]], x = .5, y = 0.33, width = .5, height = .33) +
  
  draw_plot(d_list6[[3]], x = .5, y = 0, width = .5, height = 0.33) 

p6

# To Save as pdf, provide path within quotes

#pdf(file = "", height = 14, width = 12) 

#p6

#dev.off()

