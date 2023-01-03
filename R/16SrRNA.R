#' Then run this analysis to combine data output MasterBlastR profile to upload to LabWare.
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export
#'
#'
rRNA16S_pipeline <- function(Org_id, curr_work_dir) {

  #------------------------------------------------------------------------------------------------------------
  # get directory structure
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")

  Directories.df <- as_tibble(read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir

  Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_", Org_id, "_rRNA_16S.csv", sep = ""),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE))

  Size.df <- dim(Output.df)
  NumSamples <- Size.df[1]

  LabWare.df <- tibble(SampleNo = Output.df$SampleNo,
                       Allele = Output.df$X16S_allele,
                       Identification = Output.df$X16S_comments,
                       Percent_Match = Output.df$X16S_PctWT
                       )

  write.csv(LabWare.df, paste(local_output_dir, "LabWareUpload_16S.csv", sep = ""), quote = FALSE,  row.names = FALSE)

  cat("\n\nDone! ", local_output_dir, "LabWareUpload_16S.csv is ready in output folder", "\n\n\n", sep = "")

return(LabWare.df)

}
