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

# Script name: Statistical Testing
# Created on: 14 November 2023
# Last updated: 22 April 2024
# Author: Neeraja Balasubrahmaniam
# Description: This script performs all statistical testing reported in the 
# manuscript. This includes adonis2 performed on relative abundance distance 
# matrices to check significance of different ERH groupings. This also includes
# testing to see if groups based on number of fungal annotated genes
# and fungal concentration differ by ERH condition.

# ---- Setup -------------------------------------------------------------------
# Install necessary packages
install.packages(vegan)
install.packages('devtools')
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
install.packages("tidyverse")
install.packages("readxl")

# Load packages
library(vegan)
library(pairwiseAdonis)
library(tidyverse)
library(readxl)

# Package versions
packageVersion("vegan")
packageVersion("pairwiseAdonis")
packageVersion("tidyverse")
packageVersion("readxl")

# Reporting packages versions 
# > packageVersion("vegan")
# [1] ‘2.6.4’
# > packageVersion("pairwiseAdonis")
# [1] ‘0.4.1’
# > packageVersion("tidyverse")
# [1] ‘2.0.0’
# > packageVersion("readxl")
# [1] ‘1.4.3’


#---- Test if ERH groups are significantly different based
# on gene expression data --------------

#---- Set working directory to read gene expression raw counts data
setwd("c:/Users/neerj/OneDrive - The Ohio State University/CAREER project_osuonedrive/Writing/Draft_MOLD/Significance_testing_all")

# ---- Load gene expression data for PCA

data_gene = readRDS("RSEM.gene.counts.matrix.rds") # load gene counts matrix
data_gene = data_gene[rowSums(data_gene)>=10,] # remove genes which when summed across samples < 10
cs = colSums(data_gene) 
data_gene = t( t(data_gene)/cs) * 1e6; # Counts per Million (CPM)
data_gene = log2(data_gene+1) # log transformation of CPM values
data_gene = t(data_gene)

#---- Get Euclidean distances---------------------------------------------------
euclid <- vegdist(data_gene, method = "euclidean")    # Euclidean distance

#---- Get a samples-to-RH file in the same sample order as the primary_data
meta_gene=read.table(file="sampletoRH&Site_gene.txt", header = FALSE, 
                     sep='	'); 
colnames(meta_gene)=c("SampleID","RH", "RH_level", "Site")

euclid <- vegdist(data_gene, method = "euclidean")      # Euclidean distance
euclid = euclid/1000            # scaling values between 0 and 1: does not change downstream analysis
euclid.matrix = as.matrix(euclid)               # Euclidean distance as matrix
euclid.df = as.data.frame(euclid.matrix)
euclid.df <- cbind(rownames(euclid.df), data.frame(euclid.df, row.names=NULL))
colnames(euclid.df)[1]="SampleID"

#---- Create a combined matrix with RH and distance data------------------------
euclid.df.meta = inner_join(meta_gene, euclid.df, by="SampleID")

#Are groups different based on RH: Permanova testing using adonis2
PERMANOVA_gene <- adonis2(euclid ~ RH, data = euclid.df.meta, 
                        permutations = 10000)
p_gene = PERMANOVA_gene$`Pr(>F)`[1]
R2_gene = PERMANOVA_gene$R2[1]

#---- Perform pairwise adonis as post-hoc---------------------------------------
# for gene expression
pairwise.p_all.dist = pairwise.adonis2(euclid~RH, data = euclid.df.meta, 
                                       nperm=10000)

pairwise_p_final = numeric()

pairwise_p_final["RH_50vs85_p"] = pairwise.p_all.dist[["50%_vs_85%"]][["Pr(>F)"]][1]
pairwise_p_final["RH_50vs100_p"] = pairwise.p_all.dist[["50%_vs_100%"]][["Pr(>F)"]][1]
pairwise_p_final["RH_85vs100_p"] = pairwise.p_all.dist[["85%_vs_100%"]][["Pr(>F)"]][1]

pairwise_p_adjust_final = p.adjust(pairwise_p_final, method = "BH")
pairwise_p_adjust_final = as.data.frame(pairwise_p_adjust_final)

print(pairwise_p_adjust_final)


#---- Test if sample grouping by ERH condition are significantly different based 
#---- on taxa (ASVs, species, genus) relative abundance (Bray-Curtis)

# Set working directory to read data
setwd("") # path of directory in ""

# import data
data <- readRDS("tblREL.rds") #load relative abundance of taxa

# prepare data
data_ASV <- data
data_Gen <- data %>%
  group_by(Genus) %>%
  summarise(across(where(is.numeric), sum))
data_Spec <- data %>%
  group_by(Species) %>%
  summarise(across(where(is.numeric), sum))
d_list4 <- list(data_ASV, data_Spec, data_Gen)
d_list5 <- list()
d_list6 <- list()
for (i in seq(1:3)) {
  d_list5[[i]] <- d_list4[[i]] %>%
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
  
}

# distance matrix: Bray-Curtis
d_list6 <- lapply(d_list5, function(x) vegdist(x[,-c(1,2)], method = "bray")) 
names(d_list6) <- c("ASVs", "Species", "Genus")
d_list_bray <- list()

# PERMANOVA
for (i in seq(1:3)) {
  d_list_bray[[i]] = adonis2(d_list6[[i]] ~ ERH, d_list5[[i]], permutations = 10000, 
                method = "bray")
}
names(d_list_bray) = c("ASVs_taxa_Bray", "Species_taxa_Bray", "Genus_taxa_Bray")

# The R2 and Pr(>F) values provide R2 and p-values respectively

#---- Perform pairwise adonis2 as post-hoc---------------------------------------

pairwise.p_bray = list()

for (i in seq(1:3)) {
  pairwise.p_bray[[i]] = pairwise.adonis2(d_list6[[i]] ~ ERH, data = d_list5[[i]], nperm = 10000)
}

pairwise.p_bray = pairwise.adonis2(d_list6[[1]]~ ERH, data = d_list5[[i]], 
                                       nperm=10000)

pairwise_p_final_bray = numeric()
pairwise_p_final_bray["P_100vs50_Bray_ASV"] = pairwise.p_bray[[1]][["100%_vs_50%"]][["Pr(>F)"]][1]
pairwise_p_final_bray["P_100vs85_Bray_ASV"] = pairwise.p_bray[[1]][["100%_vs_85%"]][["Pr(>F)"]][1]
pairwise_p_final_bray["P_185vs50_Bray_ASV"] = pairwise.p_bray[[1]][["50%_vs_85%"]][["Pr(>F)"]][1]

pairwise_p_final_bray["P_100vs50_Bray_Species"] = pairwise.p_bray[[2]][["100%_vs_50%"]][["Pr(>F)"]][1]
pairwise_p_final_bray["P_100vs85_Bray_Species"] = pairwise.p_bray[[2]][["100%_vs_85%"]][["Pr(>F)"]][1]
pairwise_p_final_bray["P_185vs50_Bray_Species"] = pairwise.p_bray[[2]][["50%_vs_85%"]][["Pr(>F)"]][1]

pairwise_p_final_bray["P_100vs50_Bray_Genus"] = pairwise.p_bray[[3]][["100%_vs_50%"]][["Pr(>F)"]][1]
pairwise_p_final_bray["P_100vs85_Bray_Genus"] = pairwise.p_bray[[3]][["100%_vs_85%"]][["Pr(>F)"]][1]
pairwise_p_final_bray["P_185vs50_Bray_Genus"] = pairwise.p_bray[[3]][["50%_vs_85%"]][["Pr(>F)"]][1]


pairwise_p_adjust_final_bray = p.adjust(pairwise_p_final_bray, method = "BH")
pairwise_p_adjust_final_bray = as.data.frame(pairwise_p_adjust_final_bray)

# Print FDR-adjusted p-values for all comparisons
print(pairwise_p_adjust_final_bray) 


#---- Test if sample grouping by ERH condition are significantly different based 
#---- on taxa (ASVs, species, genus) relative abundance (Aitchison)

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
d_list_aitch <- list()


# PERMANOVA
for (i in seq(1:3)) {
  d_list_aitch[[i]] = adonis2(d_list3[[i]] ~ ERH, d_list2[[i]], permutations = 10000, 
                method = "euclidean") 
}

names(d_list_aitch) = c("ASVs_taxa_Aitchison", "Species_taxa_Aitchison", "Genus_taxa_Aitchison")

# The R2 and Pr(>F) values provide R2 and p-values respectively

#---- Perform pairwise adonis as post-hoc---------------------------------------

pairwise.p_Aitch = list()

for (i in seq(1:3)) {
  pairwise.p_Aitch[[i]] = pairwise.adonis2(d_list3[[i]] ~ ERH, data = d_list2[[i]], nperm = 10000)
}


pairwise_p_final_Aitch = numeric()
pairwise_p_final_Aitch["P_100vs50_Aitch_ASV"] = pairwise.p_Aitch[[1]][["100%_vs_50%"]][["Pr(>F)"]][1]
pairwise_p_final_Aitch["P_100vs85_Aitch_ASV"] = pairwise.p_Aitch[[1]][["100%_vs_85%"]][["Pr(>F)"]][1]
pairwise_p_final_Aitch["P_85vs50_Aitch_ASV"] = pairwise.p_Aitch[[1]][["50%_vs_85%"]][["Pr(>F)"]][1]

pairwise_p_final_Aitch["P_100vs50_Aitch_Species"] = pairwise.p_Aitch[[2]][["100%_vs_50%"]][["Pr(>F)"]][1]
pairwise_p_final_Aitch["P_100vs85_Aitch_Species"] = pairwise.p_Aitch[[2]][["100%_vs_85%"]][["Pr(>F)"]][1]
pairwise_p_final_Aitch["P_85vs50_Aitch_Species"] = pairwise.p_Aitch[[2]][["50%_vs_85%"]][["Pr(>F)"]][1]

pairwise_p_final_Aitch["P_100vs50_Aitch_Genus"] = pairwise.p_Aitch[[3]][["100%_vs_50%"]][["Pr(>F)"]][1]
pairwise_p_final_Aitch["P_100vs85_Aitch_Genus"] = pairwise.p_Aitch[[3]][["100%_vs_85%"]][["Pr(>F)"]][1]
pairwise_p_final_Aitch["P_85vs50_Aitch_Genus"] = pairwise.p_Aitch[[3]][["50%_vs_85%"]][["Pr(>F)"]][1]


pairwise_p_adjust_final_Aitch = p.adjust(pairwise_p_final_Aitch, method = "BH")
pairwise_p_adjust_final_Aitch = as.data.frame(pairwise_p_adjust_final_Aitch)

# Print FDR-adjusted p-values for all comparisons
print(pairwise_p_adjust_final_Aitch) 

#---- Significance testing on number of fungal annotated genes at 
#---- each RH condition: Kruskal Wallis followed by pairwise Wilcoxon-----------

df = read.table("Fungal_gene_number_RH.txt", header=T, com='', row.names=NULL, check.names=F, sep='\t')
colnames(df) = c("Site", "RH", "DE.gene.number")

#Shapiro-wilk to test to check if normally distributed

shapiro.test(df$DE.gene.number)

#Not normally distributed, hence, use Krusla Wallis non-parametric test

# Kruskal Wallis test
kw_test = kruskal.test(df$DE.gene.number ~ df$RH, data = df)

print(kw_test[["p.value"]])

# p values is significant, therefore do Wilcoxon pairwsie tests 

Wilcoxon_RH <- pairwise.wilcox.test(df$DE.gene.number, df$RH, p.adjust.method = "BH")

print(Wilcoxon_RH[["p.value"]])

# 100% vs 50% and 100% vs 85% are significantly different.

#---- Code to run significance testing on fungal concentration at each
#---- RH condition: Kruskal Wallis followed by pairwise Wilcoxon

metadata <- read_xlsx("fungal_conc_nb.xlsx")

df = read.table("Fungal_abs_abun_RH.txt", header=T, com='', row.names=NULL, check.names=F, sep='\t')

# format metadata
colnames(metadata) <- colnames(metadata) %>%                                    
  str_remove("X7004.") %>%
  str_remove(".MSITS2a_R1.fastq") %>%
  str_replace(".100", "_100") %>%
  str_replace(".85", "_85") %>%
  str_replace(".50", "_50")
metadata$ID <- metadata$ID %>% 
  str_replace("OH_S1", "OH.1") %>%
  str_replace("OH_S2", "OH.2") %>%
  str_replace("OH_JB", "OH.3") %>%
  str_replace("SF", "CA")
metadata$Sample <- metadata$Sample %>% 
  str_replace("OH_S1", "OH.1") %>%
  str_replace("OH_S2", "OH.2") %>%
  str_replace("OH_JB", "OH.3") %>%
  str_replace("CA_SF", "CA")

# format metadata
meta1 <- metadata[,c("ID", "Sample", "RH", "Fungi")] %>%
  data.frame() %>%
  rename(`Relative humidity` = RH)
meta1$`Relative humidity` <- factor(meta1$`Relative humidity`, 
                                    levels = c("0.5", "0.85", "1"), labels = c("50%","85%","100%"))

colnames(meta1)[2] = "Site"

meta <- meta1

fungal_conc = meta[,(2:4)]
colnames(fungal_conc) = c("Site", "ERH", "Conc.")



#Shapiro-wilk to test to check if normally distributed

shapiro.test(fungal_conc$Conc.)

#Not normally distributed, hence, use Krusla Wallis non-parametric test

# Kruskal Wallis test
kw_test = kruskal.test(fungal_conc$Conc. ~ fungal_conc$ERH, data = fungal_conc)

print(kw_test[["p.value"]])

# p values is significant, therefore do Wilcoxon pairwise tests 

Wilcoxon_RH <- pairwise.wilcox.test(fungal_conc$Conc., fungal_conc$ERH, p.adjust.method = "BH")

print(Wilcoxon_RH[["p.value"]])

# 100% vs 50% and 100% vs 85% are significantly different.


