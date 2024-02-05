#' Renames submitted contig files to NML LabWare numbers
#' February 5 2024, Walter Demczuk & Shelley Peterson
#' 
#' This script makes a new LabWareUpload_metrics_SUBMITTED.csv with submitted lab nos replaced with LabWare nos.
#'
#' @param Org_code Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details Copies Contig Assemblies and VCF files from submitted IRIDA sample no to LabWare No.
#' Reads in a list of new LabWare Nos. and submitted genome sample numbers from list.csv
#' Makes a new LabWareUpload_metrics_SUBMITTED.csv with subm.lab.nos replaced with LabWare Nos.
#'
#' Input:
#' SampleNo    Variable
#' LabWare id  Submitted id in IRIDA
#'
#' SampleNo is new LabWare No
#' Variable is the submitted IRIDA sample id.
#'
#' @return A table frame containing the results of the query
#' @export

#--------------------------------------------------------------------------------------------------------
#  For troubleshooting and debugging
# Org_id <- "GONO"                  #GAS, GBS, PNEUMO or GONO
# curr_work_dir <- "C:\\WADE\\"
#--------------------------------------------------------------------------------------------------------

rename_submitted <- function(Org_id, curr_work_dir) {
  cat("\nRenaming contigs, vcf and metrics...\n")
  
  # get directory structure
  directorylist <- getdirectory(curr_work_dir, Org_id, "RENAME")
  #------------------------------------------------------------------------------------------------------------

  # Setup sample list table 
  SampleList.df <- as_tibble(read.csv(directorylist$SampleList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  SampleList <- SampleList.df
  
  # rename (actually copy) assembly files
  SampleList$SourceFile <- paste0(directorylist$contigs_dir, SampleList$Variable, ".fasta")
  SampleList$DestFile <- paste0(directorylist$contigs_dir, SampleList$SampleNo, ".fasta") 
  SampleList$contigs_copy <- file.copy(SampleList$SourceFile, SampleList$DestFile, overwrite = T)
  
  # rename (actually copy) vcf files
  SampleList$Sourcevcf <- paste0(directorylist$vcf_dir, SampleList$Variable, ".vcf")
  SampleList$Destvcf <- paste0(directorylist$vcf_dir, SampleList$SampleNo, ".vcf")
  SampleList$vcf_copy <- file.copy(SampleList$Sourcevcf, SampleList$Destvcf, overwrite = T)
  
  # rename metrics
  metrics.df <- as_tibble(read.csv(paste0(directorylist$output_dir, "LabWareUpload_Metrics.csv"), header = TRUE, sep = ",", stringsAsFactors = FALSE))
  metrics2.df <- inner_join(metrics.df, SampleList.df, by = c("Sample" = "Variable"))
  metrics3.df <- metrics2.df[,c("SampleNo", "Num_Reads", "Coverage", "MeanContigLength", "N50ContigLength", "Comments")] %>%
    rename("SampleNo" = "Sample")
  write.csv(metrics3.df, paste0(directorylist$output_dir, "LabWareUpload_Metrics_SUBMITTED.csv"), na = "", row.names = F)

  cat("\n\nLabWareUpload_Metrics_SUBMITTED.csv created.\n")
  return(select(SampleList, SampleNo, contigs_copy, vcf_copy))
}
