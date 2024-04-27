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

# Script name: MOLD_state_map_dust
# Created on: 7 June 2023
# Last updated: 21 April 2024
# Author: Neeraja Balasubrahmaniam
# Description: This script locates the homes that donated dust to 
# the study on a US state map. Inkscape [v.1.3] (https://inkscape.org/) was 
# used for further modifications of colors and markers on the .pdf map  
# that this script produced.

# ---- Setup -------------------------------------------------------------------

# install packages
install.packages("ggplot2")
install.packages("maps")
install.packages("ggmap")

# load packages
library(ggplot2)
library(maps)
library(ggmap)

# package versions
packageVersion("ggplot2")
packageVersion("maps")
packageVersion("ggmap")

# Reporting packages versions 
# > packageVersion("ggplot2")
# [1] ‘3.4.3’
# > packageVersion("maps")
# [1] ‘3.4.1’
# > packageVersion("ggmap")
# [1] ‘3.0.2’


# ---- Get data ----------------------------------------------------------------

# Get the US state map data for plotting
us_map <- map_data("state")

# Define the coordinates of the 7 locations included and the 1 location excluded 
# from the study due to poor RNA quality.
# The 8 (7 included + 1 excluded) location's latitude and longitude data are 
# removed due to IRB restrictions, however, the states are provided here along
# with the number of dust samples used in parenthesis:
# California(1), Washington(1), Colorado(1), Kansas(1), Michigan(1), Ohio (3), 
# Pennsylavia(1), Texas(excluded). 3-digit zip codes for the locations are 
# provided in Table S1.

cities <- data.frame(
  lat = c(00.0000, 00.0000, 00.0000, 00.0000, 00.0000, 00.0000, 00.0000, 00.0000),
  lon = c(00.0000, 00.0000, 00.0000, 00.0000, 00.0000, 00.0000, 00.0000, 00.0000),
  Dust_collected = c("Included", "Included","Included","Included","Included",
                     "Included","Included","Excluded")
)

# ---- Plotting ----------------------------------------------------------------
pdf(file = "Map_dust_Routput.pdf", width = 12, height = 8); #plots as a pdf

# Create a U.S. map with locations colored
ggplot() +
  geom_polygon(data = us_map, aes(x = long, y = lat, group = group), 
               fill = "lightblue", color = "white") +
  geom_point(data = cities, aes(x = lon, y = lat, color = Dust_collected), 
             size = 5) +
  scale_color_manual(values=c('blue','black')) +
  coord_fixed(ratio = 1.25) +
  theme_void()

dev.off()
