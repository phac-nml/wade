#' Detects potententially contaminated samples based on plate proximity
#' February 5 2024, Shelley Peterson
#' Run labware_GONO_AMR first, then run this analysis to evaluate sample proximity
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#Org_id <- "GONO"
#curr_work_dir <- "C:\\WADE\\"
#-------------------------------------------------------------------------------

proximity_test <- function(Org_id, curr_work_dir) {

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  directorylist <- getdirectory(curr_work_dir, Org_id, "AMR")
  #-----------------------------------------------------------------------------

  # Load datafiles
  AMR_data <- as_tibble(read.csv(paste0(directorylist$output_dir, "LabWareUpload_", Org_id, "_AMR.csv"),
                                 header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Locations_data <- read_excel(file.choose(), sheet = "Client Submission", skip = 16)
  
  #define variables
  myLetters <- toupper(letters[1:8])
  myLetters[9] <- "None"
  
  # Add samplelocation to AMR file
  AMR_data$samplelocation <- NA
  AMR_index <- match(AMR_data$lw_CurrSampleNo, Locations_data$`Sample Name`)
  AMR_data$samplelocation <- Locations_data$Well[AMR_index]
  
  # Calculate locations of all wells surrounding each sample
  AMR_data$rows <- as.numeric(gsub("[[:digit:]]", "", AMR_data$samplelocation) %>%
    match(., myLetters))
  AMR_data$cols <- as.numeric(gsub("[[:alpha:]]", "", AMR_data$samplelocation))
  AMR_data$minusrow <- AMR_data$rows - 1
  AMR_data$minusrow[AMR_data$minusrow == 0] <- 9
  AMR_data$plusrow <- AMR_data$rows + 1
  AMR_data$minuscol <- sprintf('%02d', as.numeric(AMR_data$cols - 1))
  AMR_data$minuscol[AMR_data$minuscol == 00] <- "None"
  AMR_data$pluscol <- sprintf('%02d', as.numeric(AMR_data$cols + 1))
  AMR_data$pluscol[AMR_data$pluscol == 13] <- "None"
  AMR_data$cols <- sprintf('%02d', as.numeric(AMR_data$cols))
  AMR_data$cell1 <- paste0(myLetters[AMR_data$minusrow], AMR_data$minuscol)
  AMR_data$cell2 <- paste0(myLetters[AMR_data$minusrow], AMR_data$cols)
  AMR_data$cell3 <- paste0(myLetters[AMR_data$minusrow], AMR_data$pluscol)
  AMR_data$cell4 <- paste0(myLetters[AMR_data$rows], AMR_data$minuscol)
  AMR_data$cell5 <- paste0(myLetters[AMR_data$rows], AMR_data$pluscol)
  AMR_data$cell6 <- paste0(myLetters[AMR_data$plusrow], AMR_data$minuscol)
  AMR_data$cell7 <- paste0(myLetters[AMR_data$plusrow], AMR_data$cols)
  AMR_data$cell8 <- paste0(myLetters[AMR_data$plusrow], AMR_data$pluscol)
  
  if(Org_id == "GONO")
  {
    genes <- c("lw_ermB", "lw_ermC", "lw_tetM", "lw_bla")
  } else if(Org_id == "GAS")
  {
    genes <- c("lw_ermB", "lw_ermT", "lw_mefAE", "lw_tetM", "lw_tetO", "lw_tetT",
               "lw_tetA", "lw_tetB", "lw_cat", "lw_catQ", "lw_dfrF", "lw_dfrG")
  }else if(Org_id == "GBS")
  {
    genes <- c()
  }else if(Org_id == "PNEUMO")
  {
    genes <- c()
  }
  
  proximity <- function(data, gene){
    data_filtered <- dplyr::filter(data, data[[gene]] == "POS")
    proximal_cells <- c(data_filtered$cell1, data_filtered$cell2, data_filtered$cell3,
                        data_filtered$cell4, data_filtered$cell5, data_filtered$cell6,
                        data_filtered$cell7, data_filtered$cell8)
    proximal_samples <- dplyr::filter(data_filtered, samplelocation %in% proximal_cells)
    sample_list <- proximal_samples$lw_CurrSampleNo
    return(sample_list)
  }
  
  samples <- lapply(genes, proximity, data=AMR_data)
  names(samples) <- genes
  
  samples.df <- map_df(samples, ~as_tibble(.x), .id = "Gene") %>%
    dplyr::rename("SampleNo" = "value") 
  samples.df$Gene <- gsub("lw_", "", samples.df$Gene) 
  location_index <- match(samples.df$SampleNo, AMR_data$lw_CurrSampleNo)
  samples.df$SampleLocation <- AMR_data$samplelocation[location_index]
  
  write.csv(samples.df, paste0(directorylist$output_dir, "ContaminationCheck.csv"), quote = FALSE, row.names = FALSE)  

  cat("\n\nDone! ", "ContaminationCheck.csv is ready in output folder", "\n\n\n", sep = "")
  
  return(samples.df)
}
