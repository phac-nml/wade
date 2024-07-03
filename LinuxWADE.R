##################################################################
####  WGS Analysis and Detection of Molecular Markers (WADE)  ####
####       Authors: Walter Demczuk & Shelley Peterson         ####
####                    Date: 2024-04-24                      ####
##################################################################

library(plyr)
library(tidyverse)
library(tidyselect)
library(tidysq)
library(Biostrings)
library(shiny)
library(shinyWidgets)
library(DT)
library(readxl)
library(beepr)
library(WADE)  

####### <<<< Change this path to where the WADE directory is stored >>>> #######
curr_work_dir <- "~/Documents/WADE/"
#-------------------------------------------------------------------------------

removefiles(curr_work_dir)

MASTER_pipeline("GAS", "AMR", "list", "list", curr_work_dir)
labware_GAS_amr("GAS", curr_work_dir)
MLST_pipeline("GAS", "MLST", "list", "list", curr_work_dir)
EMM_pipeline("GAS", "list", curr_work_dir)
MASTER_pipeline("GAS", "rRNA16S", "list", "list", curr_work_dir)
rRNA16S_pipeline("GAS", curr_work_dir)
MASTER_pipeline("GAS", "TOXINS", "list", "list", curr_work_dir)
labware_gas_toxins("GAS", curr_work_dir)