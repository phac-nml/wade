#' Removes duplicate fasta sequences from output_dna_notfound.fasta when using master_blastR.
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details Remove Duplicate Sequences
#'
#' Removes duplicate sequences from dna_notfound.fasta and numbers them starting at the number
#' entered in locus entry area.  These can be copied and pasted directly into the allele_lookup fastas.
#' @return A table frame containing the results of the query
#' @export
#'
#'

remove_duplicate_fasta <- function(Org_id, allele_id, curr_work_dir) {

  #------------------------------------------------------------------------------------------------------------
  # get directory structure
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  cat("\n directory string passed to function: ", dir_file, "\n")
  Directories.df <- read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  Directories.df <- as_tibble(Directories.df)
  Directories_org.df <- filter(Directories.df, OrgID == Org_id) #this one passed = "Curator"
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  cat("\n local output directory string: ", local_output_dir, "\n")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir
  #------------------------------------------------------------------------------------------------------------

  Start_num <-  as.integer(allele_id)

  FileName <- paste(local_output_dir, "output_dna_notfound.fasta", sep = "")

  con <- file(FileName, open="r")
  linn <- readLines(con)
  close(con)

  seq_table <- tibble(FastaHeader = as.character(), DNASequence = as.character())

  sample <- 1L

  for (i in seq(1, length(linn), by = 2))
  {
    seq_table[sample, "FastaHeader"] <- linn[i]
    seq_table[sample, "DNASequence"] <- linn[i+1]
    sample <- sample + 1L
  }

  #remove duplicate sequences from table
  seq_table_2 <- filter(seq_table, !is.na(FastaHeader))
  seq_table_3 <- distinct(seq_table_2, DNASequence, .keep_all = TRUE)
  
  #make new fasta
  SizeTable <- dim(seq_table_3)
  NumRows <- SizeTable[1]

  #write new fasta file
  unlink(paste(local_output_dir, "output_dna_notfound_distinct.fasta", sep="")) #this deletes the file!
  sink(paste(local_output_dir, "output_dna_notfound_distinct.fasta", sep=""), split=FALSE, append = TRUE)
  for (j in 1:NumRows)
  {
    fasta_header <- as.character(seq_table_3[j, 1])
    fasta_header_parts <- unlist(strsplit(fasta_header, "_"))

    new_fasta_header <- paste(fasta_header_parts[1], "_", as.character(Start_num), "_", fasta_header_parts[2], "_", fasta_header_parts[3], sep = "")

    cat(new_fasta_header, "\n",  as.character(seq_table_3[j,2]), "\n", sep ="")
    Start_num <- Start_num + 1L

  }
  sink()

  cat("New output_dna_notfound_distinct.fasta created. \n\n")
  return(seq_table_3)
}
