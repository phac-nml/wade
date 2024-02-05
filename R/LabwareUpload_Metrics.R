#' Labware Upload Formatter for WGS Metrics
#' February 5 2024, Walter Demczuk & Shelley Peterson
#' Browse to raw downloaded .csv files to create the upload format table for Labware
#' Doesn't matter in what order files are selected
#' Will combine columns into Labware upload ready .csv
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#Org_id <- "GONO"                  #GAS, GBS, PNEUMO or GONO
#curr_work_dir <- "C:\\WADE\\"
#-------------------------------------------------------------------------------

metrics <- function(Org_id, curr_work_dir) {

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  Test_id <- "LW_METRICS"
  directorylist <- getdirectory(curr_work_dir, Org_id, Test_id)
  #-----------------------------------------------------------------------------

  # Upload metrics files (Assembly_stats.tsv and fastqc.tsv) and combine them into one dataframe
  Files12 <- choose.files(default = "", caption = "Select Matching Files", multi = TRUE)

  File1.df <- as_tibble(read.csv(Files12[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  File2.df <- as_tibble(read.csv(Files12[2], header = TRUE, sep = "\t", stringsAsFactors = FALSE))

  colnames(File1.df)[1] <- "StrainID"
  colnames(File2.df)[1] <- "StrainID"

  combined.df <- full_join(File1.df, File2.df, by = "StrainID")
  combined.df <- combined.df[order(combined.df$StrainID),]
  combined.df$StrainID <- sub("\\__.*", "", combined.df$StrainID)

  combined.df$N50.contig.length <- as.numeric(sub("\\,", "", combined.df$N50.contig.length))

  # Select the desired metrics to output for Labware
  MetricsUpload.df <-  tibble("Sample" = combined.df$StrainID,
                              "Num_Reads" = combined.df$X..of.Reads,
                              "Coverage" = combined.df$Estimated.Coverage,
                              "MeanContigLength" = combined.df$Mean.contig.length,
                              "N50ContigLength" = combined.df$N50.contig.length,
                              "Comments" = "Pass")

  # Add "Warning" and "Fail" to the samples with poor metrics
  MetricsUpload.df$Comments[MetricsUpload.df$N50ContigLength < 30000] <- "Warning"

  MetricsUpload.df$Comments[MetricsUpload.df$Coverage < 10 |
                            MetricsUpload.df$N50ContigLength < 10000 |
                            MetricsUpload.df$MeanContigLength < 1000] <- "Fail"

  # Generate files for pass/fail/all samples and make a list.csv for WADE from the good metrics
  Metrics_bad.df <- filter(MetricsUpload.df, Comments == "Fail"| Comments == "Warning")
  Metrics_good.df <- filter(MetricsUpload.df, Comments == "Pass" | Comments == "Warning")

  #make list.csv from good metrics
  list.df <- tibble(SampleNo = MetricsUpload.df$Sample, Variable = "")
  write.csv(list.df, paste(curr_work_dir, "list.csv", sep = ""), quote = FALSE, row.names = F)

  write.csv(MetricsUpload.df, paste(directorylist$output_dir, "LabWareUpload_METRICS.csv", sep = ""), quote = FALSE, row.names = F)
  write.csv(Metrics_good.df, paste(directorylist$output_dir, "LabWareUpload_METRICS_good.csv", sep = ""), quote = FALSE, row.names = F)
  write.csv(Metrics_bad.df, paste(directorylist$output_dir, "LabWareUpload_METRICS_bad.csv", sep = ""), quote = FALSE, row.names = F)

  cat("\n\nDONE   .. Output \ LabWareUpload_METRICS.csv has been created.\n")

return(MetricsUpload.df)
}
