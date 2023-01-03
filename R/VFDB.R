#Generic Profiler : SampleNo_contigs.fasta vs. VFDB (Virulence Factor DataBase)
#R Studio scripting on local version of R Studio
#2018-08-23
#Walter Demczuk

#----------------------------------------------------------------------------------------------
#' VFDB.R
#' Generic Virulence Factor Profiler : SampleNo_contigs.fasta vs. VFDB database
#'
#' Takes Organism, Sample Number and queries a contig.fasta file
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param SampleNo Sample number (or list of sample numbers) associated with contig.fasta file
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#'
#' @details searches for virulence factors in genomic data
#'
#' @return A table frame containing the results of the query
#' @export

VFDB_pipeline <- function(Org_id, SampleNo, curr_work_dir) {

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

Variable <- "NA"
Blast_evalue <- "10e-50"

Found.df <- tibble(SampleNo=character(), DataBase=character(), GeneID=character(), MatchID=character())

if(SampleNo == "list")
{
  SampleList.df <- read.csv(SampList, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  SampleList.df <- as_tibble((SampleList.df))
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

    Found.df$SampleNo <- CurrSampleNo
    Found.df$MatchID <- "Error - contig not found."
    return(Found.df)
  }else
  {
    database <- "VFDB"
    cat("\n", CurrSampleNo, " -> ", database , " Please wait...",  sep = " ")

    BlastCommand2 <- paste("blastn -db ", system_dir, "VFDB\\VFDB.fasta -query ", ContigsDir, CurrSampleNo, ".fasta -out ", system_dir, "VFDB\\temp\\blastout_VFDB.txt -evalue ", Blast_evalue,sep = "")

        try(system(BlastCommand2))

        FileName2 <- paste(system_dir, "VFDB\\temp\\blastout_VFDB.txt", sep = "")
        con <- file(FileName2, open="r")
        linn <- readLines(con)
        close(con)

        for (i in 1:length(linn))
        {

          if (str_detect(linn[i], ">"))
          {
            AlleleLine <- unlist(linn[i])
            AlleleParts <- strsplit(AlleleLine, " ")
            AlleleLine2 <- unlist(AlleleParts)
            cat("\n", AlleleLine, sep = " ")

          }
          if (str_detect(linn[i], "Identities"))
          {
            IDLine <- unlist(linn[i])
            IDLineParts <- strsplit(IDLine, ",")
            IDLine2 <- unlist(IDLineParts)
            #cat("\t", IDLine2[1], sep = " ")
            sample.df <- tibble(CurrSampleNo, database, AlleleLine, IDLine2[1])
            Found.df <- rbind(Found.df, sample.df)
          }
      } #close bracket for blast.txt loop
  } #close bracket for sample exists
} #close brack for sample list loop<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

header <- c("SampleNo", "Database", "GeneID", "MatchID")
names(Found.df) <- header

outfile <- paste(local_output_dir, "output_profile_VFDB.csv", sep = "")
write.csv(Found.df, outfile, quote = FALSE,  row.names = F)

return(Found.df)

}
