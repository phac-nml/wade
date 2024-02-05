#' Labware Upload Formatter for PNEUMO VIRULENCE
#' February 5 2024, Walter Demczuk & Shelley Peterson
#' run MasterBlastR first, then this script outputs csv to LabWareUpload_GAS_TOXINS.csv in Output folder
#'
#' Then run this analysis to combine data the full amr profile to upload to LabWare.
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#Org_id <- "PNEUMO" 
#curr_work_dir <- "C:\\WADE\\"
#-------------------------------------------------------------------------------

labware_pneumo_virulence <- function(Org_id, curr_work_dir) {

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  directorylist <- getdirectory(curr_work_dir, Org_id, "VIRULENCE")
  #-----------------------------------------------------------------------------

  Output.df <- as_tibble(read.csv(paste0(directorylist$output_dir, "output_profile_PNEUMO_VIRULENCE.csv"),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE))

  Output.df$SampleProfile[Output.df$SampleProfile == ""] <- "Err"

  LabWare.df <- tibble(Output.df$SampleNo,
                       Output.df$PI1_result,
                       Output.df$PI2_result,
                       Output.df$lytA_result,
                       Output.df$ply_result,
                       Output.df$rpoB_allele,
                       Output.df$rpoB_mutations,
                       Output.df$rpoB_PctWT,
                       Output.df$SampleProfile)
  names(LabWare.df) <- list("SampleNo", "PI1", "PI2", "lytA", "ply", "rpoB", 
                            "rpoB_ID", "rpoB_match", "Profile")

  write.csv(LabWare.df, paste0(directorylist$output_dir, "LabWareUpload_PNEUMO_VIRULENCE.csv"), 
            quote = FALSE, row.names = FALSE)
  cat("\n\nDone! ", directorylist$output_dir, "LabWareUpload_PNEUMO_VIRULENCE.csv is ready in output folder", "\n\n\n", sep = "")

return(LabWare.df)

}