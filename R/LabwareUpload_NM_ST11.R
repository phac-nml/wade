#' ST-11 Sublineage pipeline for N. meningitidis WGS assemblies
#' May 4 2026, Shelley Peterson
#'
#' Takes Sample Number, Locus, and a Variable at queries a contig.fasta file
#' @param SampleNo Sample number associated with contig.fasta file
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#  SampleNo <- "list"
#  curr_work_dir <- "C:\\WADE\\"
#-------------------------------------------------------------------------------

NMST11_pipeline <- function(curr_work_dir) {
  
  start_time <- Sys.time()
  
  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  Org_id <- "NM"
  Test_id <- "ST11"
  
  directorylist <- getdirectory(curr_work_dir, Org_id, Test_id)
  unlink(paste0(directorylist$temp_dir, "blastout.txt")) #this deletes the file!
  #-----------------------------------------------------------------------------

  MLSToutput <- as_tibble(read.csv(paste0(directorylist$output_dir, "output_profile_", Org_id, "_MLST.csv"),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  
  ST11 <- filter(MLSToutput, ST == "11")
  
  if(dim(ST11)[1] > 0)
  {
    list.df <- tibble(SampleNo = ST11$SampleNo, Variable = "")
    write.csv(list.df, paste0(curr_work_dir, "filteredlist.csv"), quote = FALSE, row.names = F)
    fumCResults <- MASTER_pipeline(Org_id, Test_id, "filteredlist", "fumC", curr_work_dir)
    fumCResults$Sublineage <- ifelse(fumCResults$fumC_motifs == "WT", "37", "15")
    Output <- select(fumCResults, c(SampleNo, Sublineage))
    
    write.csv(Output, paste0(directorylist$output_dir, "LabWareUpload_", Org_id, "_", Test_id, ".csv"), quote = FALSE, row.names = FALSE)
  
    unlink(paste0(directorylist$temp_dir, "filteredlist.csv"))
    elapsed <- format(Sys.time() - start_time)
    cat("\n\nDone!\n\n\n", "Elapsed time: ", elapsed, "\n\n output_profile_", Org_id,
        "_", Test_id, ".csv is ready in output folder\n", sep = "")
  } else
  {
    cat(paste("No ST-11 Found", "\n"))
  }

  return(Output) 
  
} 
