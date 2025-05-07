#' M1UK typing pipeline from WGS SNVPhyl M1UK Analysis
#' February 5 2025, Walter Demczuk & Shelley Peterson
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param variant M1UK or M1DK - which SNP list should be referenced
#' @param locus Sample number associated with contig.fasta file
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details GAS emm1 M1UK typing:
#' Use the SNV_Table.csv generated from the SNVphyl Galaxy pipeline, with MGAS5005 (CP000017.2) as a mapping reference.
#' include a known M1UK as a reference so that all SNV positions will be included in Galaxy SNPhyl outputs.
#' This will make LabWareUpload_GAS_M1UK.csv in the local Output directory.
#' Output.csv has SampleNo, SNVprofile and Genotype
#' if SNV profile matches MGAS5005 reference = "Global" genotype
#' if SNV profile has all 27 different = "M1UK"
#' Otherwise = "Intermediate (Number of SNVs vs references)"
#' @return A table frame containing the results of the query

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#Org_id <- "GAS"
#variant <- "M1UK"                #M1UK or M1DK
#curr_work_dir <- "C:/WADE/"
#Blast_evalue <- "10e-50"         #sets sensitivity of Blast gene match 10e-50 to 10e-150; use 10e-5 for primers
#-------------------------------------------------------------------------------

################################################################################
################### Define Emm Variant-Specific Functions ######################
################################################################################
#' @export
# Build SNV Table for M1UK/M1DK
buildSNVtable <- function(SNVTable, SNPList, pattern_ref){
  df <- filter(SNVTable, Position %in% SNPList$Position) %>%
    select(-Status, -"#Chromosome")
  headers <- names(df)
  
  df_T <- as_tibble(t(df), .name_repair = "minimal")
  df_T$SampleNo <- headers
  df_T$SampleNo[df_T$SampleNo == "Position"] <- "SampleNo"
  colnames(df_T) <- df_T[1, ]
  
  df_T <- df_T[-1, ]
  df2 <- unite(df_T, SNVprofile, 1:(ncol(df_T)-1), sep = "", remove = TRUE)
  df3 <- tibble(SampleNo = df2$SampleNo, SNVprofile = df2$SNVprofile)
  df3$Genotype <- NA
  df3$patternref <- pattern_ref
  return(df3)
}
#' @export
# Count M1UK/M1DK SNPs
SNPcount <- function(df){
  mapply(function(x, y) {
    len <- min(length(x), length(y))
    sum(x[1:len] == y[1:len])
  }, strsplit(df$SNVprofile, ''), strsplit(df$patternref, ''))
}
#' @export
variantgenotype <- function(df, patternlength, variant){
  df$Genotype[df$SNV_total == 0] <- variant
  df$Genotype[df$SNV_total > 0] <- "Intermediate"
  df$Genotype[df$SNV_total == patternlength] <- "Global"  
  df$SNV_total[df$SNV_total == 0 | df$SNV_total == patternlength] <- ""
  df$Genotype <- paste0(df$Genotype, " (", (as.numeric(patternlength) - as.numeric(df$SNV_total)), ")") %>%
    str_remove_all(fixed(" (NA)"))
  df <- df %>% select(-patternref, -SNV_total)
  return(df)
}

################################################################################
########################### Emm Variant Pipeline ###############################
################################################################################
#' @export
EMM_V_pipeline <- function(Org_id, variant, curr_work_dir){

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  directorylist <- getdirectory(curr_work_dir, Org_id, "EMM")
  reflist <- refdirectory(directorylist, Org_id, "EMM")
  #-----------------------------------------------------------------------------

  ########################### Load Reference Files #############################
  if(variant == "M1UK")
  {
    SNP_File <- paste(reflist$Ref_Dir, "M1UK_SNP_Positions.csv", sep = "")
    SNPList.df <- as_tibble(read.csv(SNP_File, header = TRUE, sep = ",", stringsAsFactors = FALSE))
    pattern_ref <- "CACTGAGTGGCTGGGGCGGGCCCACGG"  #SNV pattern for MGAS5005 WT Ref Strain
  }

  if(variant == "M1DK")
  {
    SNP_File <- paste(reflist$Ref_Dir, "M1DK_SNP_Positions.csv", sep = "")
    SNPList.df <- as_tibble(read.csv(SNP_File, header = TRUE, sep = ",", stringsAsFactors = FALSE))
    pattern_ref <- "ATCGCGGCAAGGGTC" #SNV pattern for MGAs5005 WT Ref Strain
  }

  dna_file <- paste0(directorylist$temp_dir, "output_dna.fasta")
  unlink(dna_file) #this deletes the file!

  ##################### Load, Filter & Rearrange SNV Table #####################
  SNVTableFile <- file.choose() # *-snvTable.tsv from SNVPhyl on Galaxy
  SNVTable.df <- as_tibble(read.csv(SNVTableFile, header = TRUE, sep = "\t", check.names = "F", stringsAsFactors = FALSE))
  
  SNVTable <- buildSNVtable(SNVTable.df, SNPList.df, pattern_ref)

  ############################ Assign Genotypes ################################
  # count differences between strings
  SNVTable$SNV_total <- SNPcount(SNVTable)
  
  # assign genotype
  SNVTable <- variantgenotype(SNVTable, nchar(pattern_ref), variant)
  SNVTable <- filter(SNVTable, SampleNo != "Reference")
  
  # fix intermediates with 1-4 "-" so they're called WT
  SNVTable$Genotype[str_detect(SNVTable$Genotype, "Intermediate \\([1-4]\\)") &
                    str_detect(SNVTable$SNVprofile, "-")] <- "Global"
  
  ############################## Export Results ################################
  # Make multi-fasta file of the genome positions for each sample
  export <- paste0(">", SNVTable$SampleNo, "_", SNVTable$Genotype, "\n", 
              SNVTable$SNVprofile)
  writeLines(export, dna_file)

  # rename first column to fix labware upload bug maybe?
  SNVTable <- SNVTable %>% rename("SampleNo" = "SampleNo")
  # export results as csv file
  write.csv(SNVTable, paste0(directorylist$output_dir, "LabWareUpload_GAS_", variant, ".csv"), na = "", row.names = F)
  
  cat("\n\nDone\n\nLabWareUpload and output_dna.fasta created \n\n")

  return(SNVTable)

}

