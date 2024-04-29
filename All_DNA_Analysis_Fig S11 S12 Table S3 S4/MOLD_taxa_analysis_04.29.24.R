################################################################################
#
#                       Mold growth metatranscriptomics
#                      
#       Balasubrahmaniam et al 2024. Moving Beyond Species: Fungal function 
#       provides novel targets for potential indicators of mold growth in homes.
#
################################################################################

# Script name: MOLD_taxa_analysis
# Created on: July/August 2023
# Last updated: 04/29/2024
# Author: Jon C. King and Neeraja Balasubrahmaniam
# Purpose: This script takes 1) output from the DADA2 ITS sequencing pipeline
# and 2) corresponding data from QPCR and combines them. Measures of the fungal
# community (e.g., concentration and diversity) are then computed and analyzed.


# ---- Setup -------------------------------------------------------------------
setwd("") # set working directory


# install packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq")
install.packages("readxl")
install.packages("tidyverse")
install.packages("vegan")
install.packages("writexl")

# load packages
library(phyloseq)
library(readxl)
library(tidyverse)
library(vegan)
library(writexl)

# package version
packageVersion("phyloseq")
packageVersion("readxl")
packageVersion("tidyverse")
packageVersion("vegan")
packageVersion("writexl")

# Report package versions
# > packageVersion("phyloseq")
# [1] ‘1.42.0’
# > packageVersion("readxl")
# [1] ‘1.4.3’
# > packageVersion("tidyverse")
# [1] ‘2.0.0’
# > packageVersion("vegan")
# [1] ‘2.6.4’
# > packageVersion("writexl")
# [1] ‘1.4.2’

# ---- Import data -------------------------------------------------------------

# import data
seqtab.nochim <- readRDS("seqtab.nochim.rds")
taxa <- readRDS("taxa.rds")
metadata <- read_xlsx("fungal_conc_nb.xlsx")


# ---- Prepare data (1st round) ------------------------------------------------

# This section formats, re-arranges, and repackages outputs from DADA2 and qPCR. 
# Tables of ASV relative and absolute abundance are created first followed by 
# computation of alpha diversity measures. 

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


# format ASV table
ASV.rel <- rownames_to_column(data.frame(t(seqtab.nochim)), var = "ASV")
colnames(ASV.rel) <- colnames(ASV.rel) %>%                                    
  str_remove("X7004.") %>%
  str_remove(".MSITS2a_R1.fastq") %>%
  str_replace(".100", "_100") %>%
  str_replace(".85", "_85") %>%
  str_replace(".50", "_50") %>%
  str_replace("S1", "OH.1") %>%
  str_replace("S2", "OH.2") %>%
  str_replace("JB", "OH.3") %>%
  str_replace("SF", "CA")
# format taxonomy table
ASV.tax <- rownames_to_column(taxa, var = "ASV")
# combine ASV and taxonomy
table1 <- full_join(ASV.tax, ASV.rel, by = "ASV") %>% 
  select(-(ASV)) %>%
  mutate(Name = case_when(is.na(Species) ~ NA,                                  
                          !is.na(Species) ~ paste(str_remove(Genus, "g__"),
                                                  str_remove(Species, "s__")))) %>%
  relocate(Name, .before = Species) %>%
  select(-(Species)) %>%
  rename(Species = Name)


# check positive controls (PCs)
pc.check <- select(table1, c(Species, PC1, PC2))
pc.check %>%
  filter(PC1 > 0) %>%
  arrange(desc(PC1))
# remove positive controls (PCs) and ASVs unique to them
table2 <- table1 %>%
  select(!c("PC1", "PC2"))
sum(rowSums(select(table2, where(is.numeric))) == 0)                            # count PC-only ASVs
rm.list <- table2 %>%                                                           # find PC-only ASVs
  rowwise() %>%
  summarise(sum = sum(c_across(CO_100:WA_85))) %>%
  pull()
tblREL <- table2 %>%                                                            # remove PC-only ASVs
  slice(which(rm.list > 0)) %>%
  select(-c(Kingdom))
sum(rowSums(select(tblREL, where(is.numeric))) == 0)                            # check


# calculate absolute abundance of taxa 
read.no <- tblREL %>%                                                           # reads per sample
  summarise_if(is.numeric, sum) %>%
  unlist()
fungi <- pull(metadata, Fungi, name = ID)                                       # fungal conc
matrABS <- as.matrix(tblREL[,sapply(tblREL, is.numeric)]) %>%                   # convert rel to abs
  sweep(2, read.no, FUN = "/") %>%
  sweep(2, fungi, FUN = "*")
tblABS <- cbind(tblREL[,c("Phylum", "Class", "Order",                           # recompose table
                          "Family", "Genus", "Species")], data.frame(matrABS))
# abundance by taxonomic rank
ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")            # ranks
listREL <- list(); listABS <- list()                                            # empty list
for (x in ranks) {
  listREL[[x]] <- tblREL %>%                                                    # abundance by rank
    group_by(tblREL[[x]]) %>%
    summarise_if(is.numeric, sum)
  listABS[[x]] <- tblABS %>%
    group_by(tblABS[[x]]) %>%
    summarise_if(is.numeric, sum)
}
db.Phy <- listABS[["Phylum"]]; db.Cla <- listABS[["Class"]];                    # assign
db.Ord <- listABS[["Order"]]; db.Fam <- listABS[["Family"]];
db.Gen <- listABS[["Genus"]]; db.Spe <- listABS[["Species"]]



# rarefy for alpha diversity
sort(read.no)                                                                   # min = 43891
tblRAR <- rrarefy(t(tblREL[,sapply(tblREL, is.numeric)]), 43891) %>%            # rarefy
  t() %>%
  data.frame()
sapply(tblRAR, sum)                                                             # check


# alpha diversity table
tblALP <- data.frame(reads = sapply(tblREL[,sapply(tblREL, is.numeric)], sum),
                     reads.r = sapply(tblRAR, sum),
                     richness = sapply(tblRAR, function(x) specnumber(x, MARGIN = 1)),
                     shannon = sapply(tblRAR, function(x) diversity(x, index = "shannon"))) %>%
  rownames_to_column(var = "ID")


# format metadata
meta1 <- metadata[,c("ID", "Sample", "RH", "Fungi")] %>%
  data.frame() %>%
  rename(`Relative humidity` = RH)
meta1$`Relative humidity` <- factor(meta1$`Relative humidity`, 
                                    levels = c("0.5", "0.85", "1"), labels = c("50%","85%","100%"))

colnames(meta1)[2] = "Site"

meta <- meta1

# join metadata and alpha diversity
meta <- left_join(meta, tblALP, by = "ID")
metaCont <- meta
metaCont$`Relative humidity` <- (as.character(metaCont$`Relative humidity`)) %>%
  replace(metaCont$`Relative humidity` == "50%", "0.5") %>%
  replace(metaCont$`Relative humidity` == "85%", "0.85") %>%
  replace(metaCont$`Relative humidity` == "100%", "1") %>%
  as.numeric()


# distance matrix
ASV <- data.frame(t(tblREL[,!colnames(tblREL) %in% c("Phylum", "Class",
                                                     "Order", "Family",
                                                     "Genus", "Species")])) %>%
  rownames_to_column(var = "ID")
input.beta <- left_join(meta, ASV, by = "ID")
colnames(input.beta) <- c(colnames(input.beta[,1:8]), 
                          paste("ASV", seq(1:(2205-8)), sep = ""))

#-----------------------------------------------------------------------------------

# clean environment
rm(list = c("ASV", "ASV.rel", "ASV.tax", "db.Cla", "db.Fam", "db.Ord", "fungi", 
            "G.name", "inputASV", "inputMeta", "inputTax", "listABS", "listREL", 
            "matrABS", "meta1", "metadata", "pc.check", "P.name", "ranks", 
            "read.no", "rm.list", "rm.list2", "seqtab.nochim", "S.name", 
            "table1", "table2", "taxa", "tblABS", "tblALP", "tblRAR", "tblREL", 
            "temp1", "x"))

# -----Plotting fungal concentration with ERH -------------------------------------
# This is SI Fig 11

fungal_conc = meta[,(2:4)]
colnames(fungal_conc) = c("Site", "ERH", "Conc.")

fungal_conc$ERH = factor(x = fungal_conc$ERH, levels = c("50%", "85%", "100%"))

p1 = fungal_conc %>%
  ggplot(mapping = aes(x = ERH, y = 	
                         Conc., fill = `ERH`)) +
  stat_boxplot(geom = "errorbar", width = 0.2) + 
  geom_boxplot(colour = "black", linewidth = 0.8) +
  scale_y_continuous(trans='log10') +
  labs(x = "ERH", y = expression(paste("Fungal concentration (spore equivalents  ", mg^-1, ")"))) +
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
pdf("Boxplot_fungal_conc.pdf", height = 5, width = 5)
p1
dev.off()

# -----Plotting phyla composition in individual samples-------------------------------------

# prepare for plotting: Phyla
P.name <- c("ID", "ERH", "Site", pull(db.Phy, `tblABS[[x]]`))                # get names
P.name[13] <- "Unknown"                                                         # define 'NA'
temp1 <- select(db.Phy, -(`tblABS[[x]]`)) %>%                                   # re-arrange
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "ID") %>%
  left_join(meta[,c("ID", "Site", "Relative humidity")], by = "ID") %>%
  relocate(c(`Relative humidity`, Site), .after = ID) 
colnames(temp1) <- P.name                                                       # apply names
phy <- pivot_longer(temp1, cols = colnames(temp1)[4:13], names_to = "Phylum",
                    values_to = "conc")

# Plotting phyla composition in samples (SI Fig S12)

phy_mod = phy
phy_mod = phy_mod %>%
  mutate(across('Phylum', str_replace, 'p__Ascomycota', 'Ascomycota')) %>%
  mutate(across('Phylum', str_replace, "p__Basidiomycota", "Basidiomycota")) %>%
  mutate(across('Phylum', str_replace, 'p__Blastocladiomycota', 'Blastocladiomycota')) %>%
  mutate(across('Phylum', str_replace, 'p__Chytridiomycota', 'Chytridiomycota')) %>%
  mutate(across('Phylum', str_replace, 'p__Fungi_phy_Incertae_sedis', 'Fungi Incertae sedis')) %>%
  mutate(across('Phylum', str_replace, "p__Monoblepharomycota", "Monoblepharomycota")) %>%
  mutate(across('Phylum', str_replace, "p__Mortierellomycota", "Mortierellomycota")) %>%
  mutate(across('Phylum', str_replace, 'p__Olpidiomycota', 'Olpidiomycota')) %>%
  mutate(across('Phylum', str_replace, 'p__Rozellomycota', 'Rozellomycota'))

phy_mod$Site <- factor(phy_mod$Site, levels=c("CA", "WA", "CO", "KS", "MI", "OH.1", "OH.2", "OH.3", "PA"))
phy_mod$ERH <- factor(phy_mod$ERH, levels=c("100%", "85%", "50%"))


p1 = ggplot(data = phy_mod) +
  geom_col(mapping = aes(x = Site, y = conc, fill = Phylum), colour="black", position = "fill") +
  labs(x = "Site", y = "Relative abundance") +
  scale_fill_viridis_d(option = "viridis", direction = -1) +
  theme(axis.title.y = element_text(family = "sans", face = "plain", size = 16, color = "black", margin=margin(r=12))) +
  theme(axis.title.x = element_text(family = "sans", face = "plain", size = 16, color = "black", margin=margin(t=12))) +
  theme(axis.text.y = element_text(family = "sans", face = "plain", size = 14, color = "black")) +
  theme(axis.text.x = element_text(family = "sans", face = "plain", size = 14, color = "black", angle = 45, hjust = 1)) +
  theme(axis.ticks.y = element_line(linewidth = 0.7, linetype = "solid", color = "black")) +
  theme(axis.ticks.x = element_line(linewidth = 0.7, linetype = "solid", color = "black")) +
  theme(axis.ticks.length.y = unit(0.1, "cm")) +
  theme(axis.ticks.length.x = unit(0.1, "cm")) +
  theme(axis.line.y = element_line(linewidth = 0.7, linetype = "solid", color = "black")) +
  theme(axis.line.x = element_line(linewidth = 0.7, linetype = "solid", color = "black")) +
  theme(legend.title = element_text(family = "sans", face = "plain", size = 14, hjust = 0, color = "black")) +
  theme(legend.text = element_text(family = "sans", face = "plain", size = 14, color="black")) +
  theme(panel.background = element_blank()) +
  theme(panel.border = element_rect(fill = NA, linewidth = 0.7, linetype = "solid", color = "black")) +
  theme(panel.grid.major.y = element_blank()) +
  theme(panel.grid.major.x = element_blank()) +
  facet_wrap(vars(ERH))+
  theme(strip.background = element_rect(fill="white", linewidth = 1, linetype = "solid", color = "black"))+
  theme(strip.text = element_text(colour = 'black', size = 16, face = "bold", family = "sans"))

p1

pdf(file = "Phyla_Composition_all_RH.pdf", height = 8, width = 14)
p1
dev.off()

# -----Differential abundance testing of fungal taxa (species and genus)---------------------------

# prepare for differential abundance testing

G.name <- c("ID", "RH", pull(db.Gen, `tblABS[[x]]`))                            # get names
S.name <- c("ID", "RH", pull(db.Spe, `tblABS[[x]]`))
gen <- select(db.Gen, -(`tblABS[[x]]`)) %>%                                     # re-arrange
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "ID") %>%
  left_join(metaCont[,c("ID", "Relative humidity")], by = "ID") %>%
  relocate(`Relative humidity`, .after = ID)
spec <- select(db.Spe, -(`tblABS[[x]]`)) %>%
  t() %>%
  data.frame() %>%
  rownames_to_column(var = "ID") %>%
  left_join(metaCont[,c("ID", "Relative humidity")], by = "ID") %>%
  relocate(`Relative humidity`, .after = ID)
colnames(gen) <- G.name; colnames(spec) <- S.name                               # apply names


# subset for summary stats and diff abundance testing
gen50 <- gen[gen$RH == 0.5,]
gen85 <- gen[gen$RH == 0.85,]
gen100 <- gen[gen$RH == 1.0,]
spec50 <- spec[spec$RH == 0.5,]
spec85 <- spec[spec$RH == 0.85,]
spec100 <- spec[spec$RH == 1.0,]


# ---- Differential abundance of taxa testing and statistics -------------------------------

# summary statistics for genera abundance at humidity levels
RH50freq <- vector(); RH50p25 <- vector(); RH50p50 <- vector();                 # empty vectors
RH50p75 <- vector(); RH85freq <- vector(); RH85p25 <- vector();
RH85p50 <- vector(); RH85p75 <- vector(); RH100freq <- vector();
RH100p25 <- vector(); RH100p50 <- vector(); RH100p75 <- vector()
for (i in 1:ncol(gen50[,-(1:2)])) {                                             # freq & quantiles
  RH50freq[i] <- sum(gen50[[i+2]] > 0);
  RH50p25[i] <- quantile(gen50[[i+2]], probs = seq(0, 1, 0.25))[2];
  RH50p50[i] <- quantile(gen50[[i+2]], probs = seq(0, 1, 0.25))[3];
  RH50p75[i] <- quantile(gen50[[i+2]], probs = seq(0, 1, 0.25))[4];
}
for (i in 1:ncol(gen85[,-(1:2)])) {
  RH85freq[i] <- sum(gen85[[i+2]] > 0);
  RH85p25[i] <- quantile(gen85[[i+2]], probs = seq(0, 1, 0.25))[2];
  RH85p50[i] <- quantile(gen85[[i+2]], probs = seq(0, 1, 0.25))[3];
  RH85p75[i] <- quantile(gen85[[i+2]], probs = seq(0, 1, 0.25))[4];
}
for (i in 1:ncol(gen100[,-(1:2)])) {
  RH100freq[i] <- sum(gen100[[i+2]] > 0);
  RH100p25[i] <- quantile(gen100[[i+2]], probs = seq(0, 1, 0.25))[2];
  RH100p50[i] <- quantile(gen100[[i+2]], probs = seq(0, 1, 0.25))[3];
  RH100p75[i] <- quantile(gen100[[i+2]], probs = seq(0, 1, 0.25))[4];
}
gen.abun <- data.frame(taxon = colnames(gen50[,-(1:2)]),                        # make df
                       RH50freq, RH50p25, RH50p50, RH50p75,
                       RH85freq, RH85p25, RH85p50, RH85p75,
                       RH100freq, RH100p25, RH100p50, RH100p75)
# clean environment
rm(list = c("gen50", "gen85", "gen100", "RH50freq", "RH50p25", "RH50p50",       
            "RH50p75", "RH85freq", "RH85p25", "RH85p50", "RH85p75", "RH100freq", 
            "RH100p25", "RH100p50", "RH100p75"))     


# summary statistics for species at humidity levels
RH50freq <- vector(); RH50p25 <- vector(); RH50p50 <- vector();                 # empty vectors
RH50p75 <- vector(); RH85freq <- vector(); RH85p25 <- vector();
RH85p50 <- vector(); RH85p75 <- vector(); RH100freq <- vector();
RH100p25 <- vector(); RH100p50 <- vector(); RH100p75 <- vector()
for (i in 1:ncol(spec50[,-(1:2)])) {                                            # freq & quantiles
  RH50freq[i] <- sum(spec50[[i+2]] > 0);
  RH50p25[i] <- quantile(spec50[[i+2]], probs = seq(0, 1, 0.25))[2];
  RH50p50[i] <- quantile(spec50[[i+2]], probs = seq(0, 1, 0.25))[3];
  RH50p75[i] <- quantile(spec50[[i+2]], probs = seq(0, 1, 0.25))[4];
}
for (i in 1:ncol(spec85[,-(1:2)])) {
  RH85freq[i] <- sum(spec85[[i+2]] > 0);
  RH85p25[i] <- quantile(spec85[[i+2]], probs = seq(0, 1, 0.25))[2];
  RH85p50[i] <- quantile(spec85[[i+2]], probs = seq(0, 1, 0.25))[3];
  RH85p75[i] <- quantile(spec85[[i+2]], probs = seq(0, 1, 0.25))[4];
}
for (i in 1:ncol(spec100[,-(1:2)])) {
  RH100freq[i] <- sum(spec100[[i+2]] > 0);
  RH100p25[i] <- quantile(spec100[[i+2]], probs = seq(0, 1, 0.25))[2];
  RH100p50[i] <- quantile(spec100[[i+2]], probs = seq(0, 1, 0.25))[3];
  RH100p75[i] <- quantile(spec100[[i+2]], probs = seq(0, 1, 0.25))[4];
}
spec.abun <- data.frame(taxon = colnames(spec50[,-(1:2)]),                      # make df
                        RH50freq, RH50p25, RH50p50, RH50p75,
                        RH85freq, RH85p25, RH85p50, RH85p75,
                        RH100freq, RH100p25, RH100p50, RH100p75)

# clean environment
rm(list = c("spec50", "spec85", "spec100", "RH50freq", "RH50p25", "RH50p50", 
            "RH50p75", "RH85freq", "RH85p25", "RH85p50", "RH85p75", "RH100freq", 
            "RH100p25", "RH100p50", "RH100p75"))    



# Kruskal-wallis test (i.e., humidity versus genus) followed by pairwise 
# Wilcoxon test if Kruskal-Wallis test significant (i.e., p < 0.05).
kwlist <- list()       
pwlist <- list()
kwtest <- vector()
pwtest85v50 <- vector()     
pwtest100v50 <- vector()
pwtest100v85 <- vector()
for (i in 1:ncol(gen[,-(1:2)])) {      
  kwlist[[i]] <- kruskal.test(gen[[i+2]] ~ RH, data = gen);
  kwtest[i] <- kwlist[[i]]$"p.value";
  if (kwlist[[i]]["p.value"] < 0.05) {
    pwlist[[i]] <- pairwise.wilcox.test(gen[[i+2]], gen[["RH"]], p.adjust = "none")
    pwtest85v50[i] <- pwlist[[i]]$p.value[1,1];
    pwtest100v50[i] <- pwlist[[i]]$p.value[2,1];
    pwtest100v85[i] <- pwlist[[i]]$p.value[2,2];
  } else {
    pwlist[[i]] <- NA;
    pwtest85v50[i] <- NA;
    pwtest100v50[i] <- NA;
    pwtest100v85[i] <- NA;
  }
}


# false discovery rate; combine results
calc.q <- data.frame(taxon = colnames(gen[,-(1:2)]), kwtest, pwtest85v50,
                     pwtest100v50, pwtest100v85)
p.vals <- c(calc.q$pwtest85v50, calc.q$pwtest100v50, calc.q$pwtest100v85)
q.vals <- p.adjust(p.vals, method = "BH") %>%
  matrix(ncol = 3) %>%
  data.frame()
colnames(q.vals) <- c("q85v50", "q100v50", "q100v85")
gen.result <- bind_cols(gen.abun, calc.q, q.vals)


# clean environment
rm(list = c("calc.q", "gen", "gen.abun", "kwlist", "kwtest", "p.vals", "pwlist",
            "pwtest100v50", "pwtest100v85", "pwtest85v50", "q.vals"))


# Kruskal-wallis test (i.e., humidity versus species) followed by pairwise 
# Wilcoxon test if Kruskal-Wallis test significant (i.e., p < 0.05).
kwlist <- list()       
pwlist <- list()
kwtest <- vector()
pwtest85v50 <- vector()     
pwtest100v50 <- vector()
pwtest100v85 <- vector()
for (i in 1:ncol(spec[,-(1:2)])) {      
  kwlist[[i]] <- kruskal.test(spec[[i+2]] ~ RH, data = spec);
  kwtest[i] <- kwlist[[i]]$"p.value";
  if (kwlist[[i]]["p.value"] < 0.05) {
    pwlist[[i]] <- pairwise.wilcox.test(spec[[i+2]], spec[["RH"]], p.adjust = "none")
    pwtest85v50[i] <- pwlist[[i]]$p.value[1,1];
    pwtest100v50[i] <- pwlist[[i]]$p.value[2,1];
    pwtest100v85[i] <- pwlist[[i]]$p.value[2,2];
  } else {
    pwlist[[i]] <- NA;
    pwtest85v50[i] <- NA;
    pwtest100v50[i] <- NA;
    pwtest100v85[i] <- NA;
  }
}


# false discovery rate; combine results
calc.q <- data.frame(taxon = colnames(spec[,-(1:2)]), kwtest, pwtest85v50,
                     pwtest100v50, pwtest100v85)
p.vals <- c(calc.q$pwtest85v50, calc.q$pwtest100v50, calc.q$pwtest100v85)
q.vals <- p.adjust(p.vals, method = "BH") %>%
  matrix(ncol = 3) %>%
  data.frame()
colnames(q.vals) <- c("q85v50", "q100v50", "q100v85")
spec.result <- bind_cols(spec.abun, calc.q, q.vals)


# clean environment
rm(list = c("calc.q", "kwlist", "kwtest", "p.vals", "pwlist", "pwtest100v50", 
            "pwtest100v85", "pwtest85v50", "q.vals", "spec", "spec.abun"))


# export results
# Results also in SI Tables S3 and S4
write_xlsx(gen.result, "diff.abun.genus.xlsx")
write_xlsx(spec.result, "diff.abun.species.xlsx")


















































