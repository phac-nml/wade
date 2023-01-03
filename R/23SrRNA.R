#' 23SrRNA.R
#' 23S rRNA pipeline for WGS assemblies to determine number of mutated alleles
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param SampleNo Sample number or "list" or "folder" of sample numbers associated with VCF file(s)
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @details
#' 23S rRNA pipeline for WGS assemblies to determine number of mutated alleles
#'
#' E.coli: A2059G and C2611T
#'
#' GONO:
#' Run SNP core pipeline with NCCP11945_23S4.fasta file: A2045G or C2597T

#' PNEUMO:
#' Run SNP core pipeline with 23S_R6.fasta file: A2061G or C2613T
#'
#' Takes Organism, Sample Number, Locus, and a Variable at queries a contig.fasta file
#' #Parses 23s rRNA mutations from VCF files
#'
#' ON GALAXY:
#'
#' Alternative_allele_proporition = 0.1
#' min_coverage = 15
#' min_mean_mapping = 30
#' run_name = 23S
#' in FreeBayes step of phylogeny set ploidy = 4
#'
#' Export VCF from Galaxy:
#' Under "Collection Tools" , use the tool "Export to Warehouse" on "Filter vcf on collection 210".
#' It will bundle them together and makes a history item.  They are in /Warehouse/Temporary/galaxy_exports/RANDOM_NAME.
#' Click the "eye" to view RANDOM_NAME. You can then copy them from there using the regular file explorer into:
#'             " folder
#'
#' @export


rRNA23S_pipeline <- function(Org_id, SampleNo, curr_work_dir) {

setwd(here())
work_dir <- getwd()
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
vcf_folder <- Directories_org.df$VCFDir
#------------------------------------------------------------------------------------------------------------

Variable <- NA

if (Org_id == "GONO")
{
  rRNA23S_position1 <- "23S4_NCCP11945	2045"
  rRNA23S_position2 <- "23S4_NCCP11945	2597"
}
if (Org_id == "PNEUMO")
{
  rRNA23S_position1 <- "23S_rRNA_R6_sprr02	2061"
  rRNA23S_position2 <- "23S_rRNA_R6_sprr02	2613"
}

OutputProfile.df <- tibble(SampleNo=character(), A2059G=integer(), C2611T=integer())
#------------------------------------------------------------------------------------------------
if(SampleNo == "list")
{
  SampleList.df <- read.csv(SampList, header = TRUE, sep = ",")
}else
{
  SampleList.df <- tibble(SampleNo, Variable)
}

Size.df <- dim(SampleList.df)
NumSamples <- Size.df[1]
#-----------------------------------------------------------------------------------------------

cat("\n\n23s rRNA counts from VCF\n")

m<-1
for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
{
  A2059G <- 99L
  C2611T <- 99L
  CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
  QueryFile <- paste(vcf_folder, "\\", CurrSampleNo, ".vcf", sep = "")

  if(!file.exists(QueryFile))
  {
    A2059G <- NA
    C2611T <- NA
    QueryFile <- "File not found"
    cat("Sample Number ", CurrSampleNo, " not found!\n", sep = "")
  }

  if (QueryFile != "File not found")
  {
    allele_fraction <- 0L
    con <- file(QueryFile, open="r")
    linn <- readLines(con)
    close(con)
    for (i in 1:length(linn))
    {
      if (str_detect(linn[i], rRNA23S_position1))  #2597 & 2045 for GONO; 2061 for pneumo
      {
        vcf_parts <- strsplit(linn[i], ";")
        vcf_parts <- unlist(vcf_parts)
        DP1 <- vcf_parts[8]
        DP<-as.integer(substr(DP1, 4, 10))
        AO1 <- vcf_parts[6]
        AO<-as.integer(substr(AO1, 4, 10))

        allele_fraction <- as.integer((AO/DP)*100)
        if (allele_fraction < 14){A2059G <- 0L}
        if ((allele_fraction >= 15) && (allele_fraction <= 34)) {A2059G <- 1L}
        if ((allele_fraction >= 35) && (allele_fraction <= 64)) {A2059G <- 2L}
        if ((allele_fraction >= 65) && (allele_fraction <= 84)) {A2059G <- 3L}
        if (allele_fraction >= 85) {A2059G <- 4L}
      }else
      {
        A2059G <- 0L
      }

    if (str_detect(linn[i], rRNA23S_position2))  #2597 & 2045 for GONO; 2061 for pneumo
    {
      vcf_parts2 <- strsplit(linn[i], ";")
      vcf_parts2 <- unlist(vcf_parts2)

      DP2_str <- vcf_parts2[8]
      DP2 <- as.integer(substr(DP2_str, 4, 10))   #total number of reads (depth of reads)
      AO2_str <- vcf_parts2[6]
      AO2 <- as.integer(substr(AO2_str, 4, 10))  #alternate observations from reference



      allele_fraction <- as.integer((AO2/DP2)*100)
      if (allele_fraction < 14){C2611T <- 0L}
      if ((allele_fraction >= 15) && (allele_fraction <= 34)) {C2611T <- 1L}
      if ((allele_fraction >= 35) && (allele_fraction <= 64)) {C2611T <- 2L}
      if ((allele_fraction >= 65) && (allele_fraction <= 84)) {C2611T <- 3L}
      if (allele_fraction >= 85) {C2611T <- 4L}

    }else
    {
      C2611T <- 0L
    }

  }
  }

cat(CurrSampleNo, "\tAC_2059: ", A2059G, "\tAC_2611: ", C2611T,"\n", sep = "")

SampleProfile.df <- tibble(CurrSampleNo, A2059G, C2611T)
OutputProfile.df <- rbind(OutputProfile.df, SampleProfile.df)

}

header <- list("SampleNo", "A2059G", "C2611T")
names(OutputProfile.df) <- header
outfile <- paste(local_output_dir, "output_profile_23S.csv", sep = "")
write.csv(OutputProfile.df, outfile, row.names = F)

cat("\n\nDone!\n\n\n")

return(OutputProfile.df)

}



