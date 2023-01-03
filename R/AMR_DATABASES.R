#Generic Profiler : SampleNo_contigs.fasta vs. ARG_ANNOT database
#R Studio scripting on local version of R Studio
#2018-08-23
#Walter Demczuk

#----------------------------------------------------------------------------------------------
#' AMR_DATABASES.R
#' Generic AMR Profiler : SampleNo_contigs.fasta vs. ARG_ANNOT database
#'
#' Takes Organism, Sample Number and queries a contig.fasta file
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param SampleNo Sample number (or list of sample numbers) associated with contig.fasta file
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details GAS AMR loci:
#' GAS AMR loci
#' dfrF, dfrG, DHFR, ermA, ermB, ermT, folA, folP, gyrA, mefAE, msrD, parC, rpsJ, tetM, tetO
#' note that ermA == ermTR, ermC == ermT, mefA/E includes mefA and mefE; DHFR is the same as folA
#'
#'
#' @return A table frame containing the results of the query
#' @export

AMR_DATABASE_pipeline <- function(Org_id, SampleNo, curr_work_dir) {

#------------------------------------------------------------------------------------------------------####
Variable <- "NA"
Blast_evalue <- "10e-100"
#------------------------------------------------------------------------------------------------------------
# get directory structure
curr_dir <- curr_work_dir
dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
Directories.df <- read.csv(dir_file, header = TRUE, sep = ",")
Directories_org.df <- filter(Directories.df, OrgID == Org_id)
local_dir <- Directories_org.df$LocalDir
SampList <- paste(local_dir, "list.csv", sep = "")
local_output_dir <- paste(local_dir, "Output\\", sep = "")
local_temp_dir <- paste(local_dir, "temp\\", sep = "")
system_dir <- Directories_org.df$SystemDir
ContigsDir <- Directories_org.df$ContigsDir
#------------------------------------------------------------------------------------------------------------

AMR_Found.df <- tibble(SampleNo=character(), DataBase=character(), GeneID=character(), MatchID=character())
headers <-  c("SampleNo", "DataBase", "GeneID", "MatchID")

if(SampleNo == "list")
{
  SampleList.df <- read.csv(SampList, header = TRUE, sep = ",")
}else
{
  SampleList.df <- tibble(SampleNo, Variable)
}

Size.df <- dim(SampleList.df)
NumSamples <- Size.df[1]

m<-1L

for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
{
  CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
  CurrSampleVar <-as.character(SampleList.df[m, "Variable"])
  CurrSample.df <- filter(SampleList.df, SampleNo == CurrSampleNo)
  SampleProfile <- ""

  AlleleLine <- ""
  IDLine <- ""

  SampleQuery <- paste(ContigsDir, CurrSampleNo, ".fasta", sep = "")
  if(!file.exists(SampleQuery))
  {
    sample_error.df <- tibble(Loucs_ID = locus, Output = "Sample number not found.", stringsAsFactors = FALSE)
    return(sample_error.df)
  }else
  {
    database <- "ARG-ANNOT"
    cat("\n", CurrSampleNo, " -> ", database , sep = " ")
    
    BlastCommand2 <- paste("blastn -db ", system_dir, "Resistance_Genes\\arg_annot\\arg-annot.fasta -query ", ContigsDir, CurrSampleNo, ".fasta -out ", local_temp_dir, "blastout_arg.txt -evalue ", Blast_evalue,sep = "")

    try(system(BlastCommand2))

    FileName2 <- paste(local_temp_dir, "blastout_arg.txt", sep="")
    con <- file(FileName2, open="r")
    linn <- readLines(con)
    close(con)

    #parse blastout2

    for (i in 1:length(linn))
    {
      if (str_detect(linn[i], ">"))
      {
        AlleleLine <- unlist(linn[i])
        AlleleParts <- strsplit(AlleleLine, ":")
        AlleleLine2 <- unlist(AlleleParts)
        cat("\n", AlleleLine2[1], sep = " ")
      }
      if (str_detect(linn[i], "Identities"))
      {
        IDLine <- unlist(linn[i])
        IDLineParts <- strsplit(IDLine, ",")
        IDLine2 <- unlist(IDLineParts)
        cat("\t", IDLine2[1], sep = " ")
        AMR_sample.df <- tibble(CurrSampleNo, database, substr(AlleleLine2[1], 3, 100), IDLine2[1])
        names(AMR_sample.df) <- headers
        AMR_Found.df <- rbind(AMR_Found.df, AMR_sample.df)
      }
    }

    #...................................................
    database <- "ResFinder"

    cat("\n\n", CurrSampleNo, " -> ", database, sep = " ")

    BlastCommand2 <- paste("blastn -db ", system_dir, "Resistance_Genes\\resfinder2\\resfinder2.fasta -query ", ContigsDir, CurrSampleNo, ".fasta -out C:\\WGS_Typing\\temp\\blastout_resfinder.txt -evalue ", Blast_evalue, sep = "")

    try(system(BlastCommand2))

    FileName2 <- paste(local_temp_dir, "blastout_resfinder.txt", sep = "")
    con <- file(FileName2, open="r")
    linn <- readLines(con)
    close(con)

    #parse blastout2

    for (i in 1:length(linn))
    {
      if (str_detect(linn[i], ">"))
      {
        AlleleLine <- unlist(linn[i])
        AlleleParts <- strsplit(AlleleLine, ":")
        AlleleLine2 <- unlist(AlleleParts)
        cat("\n", AlleleLine2[1], sep = " ")
      }
      if (str_detect(linn[i], "Identities"))
      {
        IDLine <- unlist(linn[i])
        IDLineParts <- strsplit(IDLine, ",")
        IDLine2 <- unlist(IDLineParts)
        cat("\t", IDLine2[1], sep = " ")

        AMR_sample.df <- tibble(CurrSampleNo, database, substr(AlleleLine2[1], 3, 100), IDLine2[1])
        names(AMR_sample.df) <- headers
        AMR_Found.df <- rbind(AMR_Found.df, AMR_sample.df)
      }
    }

    #...................................................
    database <- "CARD"

    cat("\n\n", CurrSampleNo, " -> ", database, sep = " ")

    BlastCommand2 <- paste("blastn -db ", system_dir, "Resistance_Genes\\CARD\\CARD.fasta -query ", ContigsDir, CurrSampleNo, ".fasta -out C:\\WGS_Typing\\temp\\blastout_CARD.txt -evalue ", Blast_evalue, sep = "")

    try(system(BlastCommand2))

    FileName2 <- paste(local_temp_dir, "blastout_CARD.txt", sep = "")
    con <- file(FileName2, open="r")
    linn <- readLines(con)
    close(con)

    #parse blastout2

    for (i in 1:length(linn))
    {
      if (str_detect(linn[i], ">"))
      {
        AlleleLine <- unlist(linn[i])
        AlleleParts <- strsplit(AlleleLine, "|", fixed=TRUE)
        AlleleLine2 <- unlist(AlleleParts)
        cat("\n", AlleleLine2[6], sep = " ")
      }
      if (str_detect(linn[i], "Identities"))
      {
        IDLine <- unlist(linn[i])
        IDLineParts <- strsplit(IDLine, ",")
        IDLine2 <- unlist(IDLineParts)
        cat("\t", IDLine2[1], sep = " ")

        AMR_sample.df <- tibble(CurrSampleNo, database, AlleleLine2[6], IDLine2[1])
        names(AMR_sample.df) <- headers
        AMR_Found.df <- rbind(AMR_Found.df, AMR_sample.df)
      }
    }
  } #close bracket for sample exists
} #close brack for sample list loop<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

cat("\n\nDone!\n\n\n")

return(AMR_Found.df)
}
