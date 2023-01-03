#' Browse to raw downloaded .csv files to create the upload format table for Labware
#' Doesn't matter in what order files are selected
#' Will combine columns into Labware upload ready .csv
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export
#'
#'

metrics <- function(Org_id, curr_work_dir) {

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

Files12 <- choose.files(default = "", caption = "Select Matching Files", multi = TRUE)

File1 <- Files12[1]
File2 <- Files12[2]

File1.df <- as_tibble(read.csv(File1, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
File2.df <- as_tibble(read.csv(File2, header = TRUE, sep = "\t", stringsAsFactors = FALSE))

colnames(File1.df)[1] <- "StrainID"
colnames(File2.df)[1] <- "StrainID"

File_1_2.df <- full_join(File1.df, File2.df, by = "StrainID")
File_1_2.df <- File_1_2.df[order(File_1_2.df$StrainID),]

File_1_2.df$StrainID <- sub("\\__.*", "", File_1_2.df$StrainID)

File_1_2.df$N50.contig.length <- as.numeric(sub("\\,", "", File_1_2.df$N50.contig.length))

MetricsUpload.df <-  tibble(File_1_2.df$StrainID,
                                  File_1_2.df$X..of.Reads,
                                  File_1_2.df$Estimated.Coverage,
                                  File_1_2.df$Mean.contig.length,
                                  File_1_2.df$N50.contig.length,
                                  Comments = "Pass")

headers <- c("Sample", "Num_Reads", "Coverage", "MeanContigLength", "N50ContigLength", "Comments")

names(MetricsUpload.df) <- headers

#for some unknown reason need to write then read the table to make numeric values
write.csv(MetricsUpload.df, paste(local_output_dir, "output_METRICS.csv", sep = ""), quote = FALSE, row.names = F)

MetricsUpload.df$Comments[MetricsUpload.df$N50ContigLength < 30000] <- "Warning"

MetricsUpload.df$Comments[MetricsUpload.df$Coverage < 10 |
                            MetricsUpload.df$N50ContigLength < 10000 |
                            MetricsUpload.df$MeanContigLength < 1000] <- "Fail"

Metrics_bad.df <- filter(MetricsUpload.df, Comments == "Fail")
Metrics_good.df <- filter(MetricsUpload.df, Comments == "Pass")

#make list.csv from good metrics?
list.df <- tibble(SampleNo = MetricsUpload.df$Sample, Variable = "")
write.csv(list.df, paste(local_dir, "list.csv", sep = ""), quote = FALSE, row.names = F)

write.csv(MetricsUpload.df, paste(local_output_dir, "LabWareUpload_METRICS.csv", sep = ""), quote = FALSE, row.names = F)
write.csv(Metrics_good.df, paste(local_output_dir, "LabWareUpload_METRICS_good.csv", sep = ""), quote = FALSE, row.names = F)
write.csv(Metrics_bad.df, paste(local_output_dir, "LabWareUpload_METRICS_bad.csv", sep = ""), quote = FALSE, row.names = F)

cat("\n\nDONE   .. Output \ LabWareUpload_METRICS.csv has been created.\n")

return(MetricsUpload.df)
}
