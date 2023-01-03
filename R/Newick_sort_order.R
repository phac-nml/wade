#' Gets the sort order for samples from a newick tree file.
#' Browse to newick file and gets the sample id sort order.
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export
#'
#'

nwk_sort_order <- function(Org_id, curr_work_dir) {

#Code-----------------------------------------####
#library(dplyr)
#library(stringr)
#---------------------------------------------

  #------------------------------------------------------------------------------------------------------------
  # get directory structure
  curr_dir <- curr_work_dir
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  Directories.df <- read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir

FileName <- file.choose()
con <- file(FileName, open="r")
linn <- readLines(con)
close(con)

ListLine <- linn %>%
  str_match_all("[\\(\\,]?[\\(\\,](.*?)\\:")

Samples <- ListLine[[1]][,2]
Samples <- gsub("\\(", "", Samples)
NumSamples <- length(Samples)
TreeSort <-seq(1:NumSamples)
SampleList <- tibble(TreeSort, Samples)

write.csv(SampleList, paste(local_output_dir, "TreeOrder.csv", sep = ""), row.names = F)

cat("\n\nDONE! ... Output - TreeOrder.csv created.")

return(SampleList)

}
