#' Replace sequencing sample numbers in a newick tree file
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details Replace node labels in newick file:
#' Replace sequencing sample numbers in a newick tree file
#' with custom sample numbers in .csv files
#' with old number label in first column with name "Samples:
#' with new number label in second column with name "TreeLabel"
#' @return A table frame containing the results of the query
#' @export
#'
#'

nwk_rename_samples <- function(Org_id, curr_work_dir) {

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

  #------------------------------------------------------------------------------------------------------------

unlink(paste(local_output_dir, "TreeRenamed.newick", sep="")) #this deletes the file!

Files12 <- choose.files(default = "", caption = "Select Files", multi = TRUE)

if(str_detect(Files12[1], ".newick")) {treefile <- Files12[1]} else {trefile <- Files12[2]}
if(str_detect(Files12[2], ".csv")) {dataFile <- Files12[2]} else {dataFile <-Files12[1]}

con <- file(treefile, open="r")
linn <- readLines(con)
close(con)

TreeLabels.df <- read.csv(dataFile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
Size.df <- dim(TreeLabels.df)
NumSamples <- Size.df[1]

names(TreeLabels.df) <- c("Samples", "TreeLabel")

i<-1L
for (i in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Loop through sample table
{
  Old_Value <- TreeLabels.df$Samples[i]
  New_Value <- TreeLabels.df$TreeLabel[i]
  linn<-gsub(Old_Value, New_Value, linn)          # both of these work when excuted line by line, but not in for loop???
}

sink(paste(local_output_dir, "TreeRenamed.newick", sep = ""), split=FALSE, append = FALSE)
cat(linn, sep ="")
sink()

cat("\n\n DONE!  Output - TreeRenamed.newick created.")

return()

}
