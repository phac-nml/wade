#' Gets the sort order for samples from a newick tree file.
#' February 5 2024, Walter Demczuk & Shelley Peterson
#' 
#' Browse to newick file and gets the sample id sort order.
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

nwk_sort_order <- function(Org_id, curr_work_dir) {

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  directorylist <- getdirectory(curr_work_dir, Org_id, "TREES")
  #-----------------------------------------------------------------------------

  FileName <- file.choose()
  con <- file(FileName, open="r")
  linn <- readLines(con)
  close(con)

  ListLine <- linn %>% str_match_all("[\\(\\,]?[\\(\\,](.*?)\\:")

  Samples <- ListLine[[1]][,2]
  Samples <- gsub("\\(", "", Samples)
  TreeSort <-seq(1:length(Samples))
  SampleList <- tibble(TreeSort, Samples)

  write.csv(SampleList, paste0(directorylist$output_dir, "TreeOrder.csv"), row.names = F)

  cat("\n\nDONE! ... Output - TreeOrder.csv created.")

return(SampleList)

}
