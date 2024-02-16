#' 16SrRNA Pipeline
#' February 5 2024, Walter Demczuk & Shelley Peterson
#' Run this analysis to combine data output MasterBlastR profile to upload to LabWare.
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export
#'
#'

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#Org_id <- "GONO"                  #GAS, GBS, or PNEUMO
#curr_work_dir <- "C:\\WADE\\"
#Blast_evalue <- "10e-50"         #sets sensitivity of Blast gene match 10e-50 to 10e-150; use 10e-5 for primers
#-------------------------------------------------------------------------------

rRNA16S_pipeline <- function(Org_id, curr_work_dir) {

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  Test_id <- "rRNA16S"
  directorylist <- getdirectory(curr_work_dir, Org_id, Test_id)
  #-----------------------------------------------------------------------------

  Output.df <- as_tibble(read.csv(paste0(directorylist$output_dir, "output_profile_", Org_id, "_rRNA16S.csv"),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE))

  Labware.df <- tibble(SampleNo = Output.df$SampleNo,
                       Allele = Output.df$X16S_allele,
                       Identification = Output.df$X16S_comments,
                       Percent_Match = Output.df$X16S_PctWT)
  
  Labware.df$Percent <- as.numeric(str_extract(Labware.df$Percent_Match,  "(?<=\\().+?(?=\\%)"))
  Labware.df$Allele[Labware.df$Allele == ""| Labware.df$Percent < 97] <- "NF"
  Labware.df$Identification[Labware.df$Identification == ""| Labware.df$Percent < 97] <- "NF"
  
  Labware.df <- select(Labware.df, -Percent)
  
  write.csv(Labware.df, paste0(directorylist$output_dir, "LabWareUpload_", Org_id, "_rRNA16S.csv"), 
            quote = FALSE,  row.names = FALSE)

  cat("\n\nDone! ", directorylist$output_dir, "LabWareUpload_rRNA16S.csv is ready in output folder", "\n\n\n", sep = "")

return(Labware.df)

}
