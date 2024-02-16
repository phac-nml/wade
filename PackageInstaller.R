#####################################
####   WADE Package Downloader   ####
####  Author: Shelley Peterson   ####
####     Date: 2023-12-21        ####
#####################################


if(!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("ggtree", "Biostrings"))

install.packages("plyr")
install.packages("tidyverse")
install.packages("tidyselect")
install.packages("tidysq")
install.packages("shiny")
install.packages("shinyWidgets")
install.packages("DT")
install.packages("readxl")
install.packages("beepr")

