#' Renames submitted WGS reads to NML LabWare numbers
#' July 3, 2024, Shelley Peterson
#' 
#' This script renames submitted PT reads to NML numbers for upload into IRIDA
#'
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details Copies submitted raw reads from submitted sample no to LabWare No.
#' Reads in a list of new LabWare Nos. and submitted read file numbers from metadata.csv
#'
#'
#' @return Copied files with the new NML lab numbers, outputs samples which failed to copy into the R console
#' @export

#--------------------------------------------------------------------------------------------------------
#  For troubleshooting and debugging
# curr_work_dir <- "C:\\WADE\\"
#--------------------------------------------------------------------------------------------------------

rename_fastqs <- function(x) {
  metadatafile <- file.choose()
  setwd(dirname(metadatafile))
  data <- read.csv(metadatafile, skip = 1, header = T)
  data$NewR1 <- gsub(".*_R1", "_R1", data$R1.File.Name)
  data$NewR1 <- paste0(data$NML.No., data$NewR1)
  data$NewR2 <- gsub(".*_R2", "_R2", data$R2.File.Name)
  data$NewR2 <- paste0(data$NML.No., data$NewR2)
  
  # rename (actually copy) files
  data$RenamedR1 <- file.copy(data$R1.File.Name, data$NewR1, overwrite = T)
  data$RenamedR2 <- file.copy(data$R2.File.Name, data$NewR2, overwrite = T) 
  
  #Make SampleList.csv
  samplelist <- data %>% filter(RenamedR1 == "TRUE" & RenamedR2 == "TRUE") %>%
    select(NML.No., NewR1, NewR2) %>%
    dplyr::rename("Sample_Name" = "NML.No.",
                  "File_Forward" = "NewR1",
                  "File_Reverse" = "NewR2") %>%
    dplyr::mutate(Project_ID = "2647", .after = "Sample_Name") %>%
    dplyr::add_row("Sample_Name" = "Sample_Name", "Project_ID" = "Project_ID", 
                   "File_Forward" = "File_Forward", "File_Reverse" = "File_Reverse",
                   .before = 1) %>%
    dplyr::add_row("Sample_Name" = "[Data]", "Project_ID" = "",
                   "File_Forward" = "", "File_Reverse" = "", .before = 1)
  write.table(samplelist, "SampleList.csv", sep = ",", 
              row.names = FALSE, col.names = FALSE)
  
  #Check for failed files
  failtest <- filter(data, (RenamedR1 == "FALSE" | RenamedR2 == "FALSE")) 
  if(nrow(failtest) > 0)
  {
    write.csv(failtest, "failed_rename.csv")
    failtest$NML.No. <- paste0(failtest$NML.No., "\n")
    cat("\n Done! \n The following samples failed to be renamed: \n", failtest$NML.No., "\n")
  } else 
  {
    cat("\n Done. All submissions successfully renamed \n")
  }
  
}


