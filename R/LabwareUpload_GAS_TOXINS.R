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

labware_gas_toxins <- function(Org_id, curr_work_dir) {

  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  Directories.df <- as_tibble(read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir

  Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_GAS_TOXINS.csv", sep = ""),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE))

  Size.df <- dim(Output.df)
  NumSamples <- Size.df[1]
  NumLoci <- Size.df[2]  # if single locus queried num.columns = 10 : 2 (Sample id + Variable) + 8

  if (NumLoci > 10)
  {
    Output.df$SampleProfile[Output.df$SampleProfile == ""] <- "Error"
    LabWare.df <- tibble(Output.df$SampleNo,
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
                         Output.df$SampleProfile
                        )
  headers <- list("SampleNo", "speA", "speC", "speG", "speH", "speI", "speJ", "speK", "speL", "speM","smeZ", "ssa",
                  "sagA", "TOXIN Profile")

  names(LabWare.df) <- headers

  write.csv(LabWare.df, paste(local_output_dir, "LabWareUpload_GAS_TOXINS.csv", sep = ""), quote = FALSE, row.names = FALSE)

  cat("\n\nDone! ", local_output_dir, "LabWareUpload_GAS_TOXINS.csv is ready in output folder", "\n\n\n", sep = "")

  }else  # else if a single loci was run return the original output.csv
  {
    LabWare.df <- Output.df
  }

  return(LabWare.df)
}