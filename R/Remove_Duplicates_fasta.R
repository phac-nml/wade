#' Fasta Duplicates Remover
#' February 5 2024, Shelley Peterson
#' Removes duplicate sequences from fasta file.
#'
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details Remove Duplicate Sequences
#'
#' Removes duplicate sequences from fasta file.
#' @return A fasta file containing the results of the query
#' @export

remove_duplicates_fasta <- function(curr_work_dir) {

  FileName <- file.choose()
  data <- read_fasta(FileName)
  data$sq <- as.character(data$sq)
  data2 <- distinct(data, sq, .keep_all = TRUE)
  data2 <- data2[,c("name", "sq")]
  data2$name <- paste0(">", data2$name)
  
  outfile <- paste0(curr_work_dir, "distinct_seqs.fasta")
  
  write.table(data2, file = outfile, row.names = FALSE, col.names = FALSE, 
              quote = FALSE, sep = "\n")

  cat("distinct_seqs.fasta created. \n\n")
  
  return(data2)
}
