#' Labware Upload Formatter for GAS Toxins
#' February 5 2025, Walter Demczuk & Shelley Peterson
#' run MasterBlastR first, then this script outputs csv to LabWareUpload_GAS_TOXINS.csv in Output folder
#'
#' Then run this analysis to combine data the full amr profile to upload to LabWare.
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
#Org_id <- "GAS"                  #GAS, PNEUMO or GONO
#curr_work_dir <- "C:\\WADE\\"
#-------------------------------------------------------------------------------

labware_gas_toxins <- function(Org_id, curr_work_dir) {

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  directorylist <- getdirectory(curr_work_dir, Org_id, "TOXINS")
  #-----------------------------------------------------------------------------
  
  Output.df <- as_tibble(read.csv(paste0(directorylist$output_dir, "output_profile_GAS_TOXINS.csv"),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE))

  Output.df$SampleProfile[Output.df$SampleProfile == ""] <- "Error"
  Labware.df <- tibble(Output.df$SampleNo,
                       Output.df$speA_result,
                       Output.df$speC_result,
                       Output.df$speG_result,
                       Output.df$speH_result,
                       Output.df$speI_result,
                       Output.df$speJ_result,
                       Output.df$speK_result,
                       Output.df$speL_result,
                       Output.df$speM_result,
                       Output.df$smeZ_result,
                       Output.df$ssa_result,
                       Output.df$sagA_result,
                       Output.df$SampleProfile)
    
  names(Labware.df) <- c("SampleNo", "speA", "speC", "speG", "speH", "speI", 
                         "speJ", "speK", "speL", "speM","smeZ", "ssa", "sagA",
                         "TOXIN Profile")
  
  Labware.df$`TOXIN Profile` <- gsub("-sagA", "", Labware.df$`TOXIN Profile`)
  Labware.df$`TOXIN Profile` <- gsub("sagA", "No Toxins Found", Labware.df$`TOXIN Profile`)
  Labware.df$`TOXIN Profile` <- gsub("Susceptible", "No Toxins Found", Labware.df$`TOXIN Profile`)

  write.csv(Labware.df, paste(directorylist$output_dir, "LabWareUpload_GAS_TOXINS.csv", sep = ""),
            quote = FALSE, row.names = FALSE)

  cat("\n\nDone! ", directorylist$output_dir, "LabWareUpload_GAS_TOXINS.csv is ready in output folder", "\n\n\n", sep = "")

  return(Labware.df)
}