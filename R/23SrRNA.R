#' 23S rRNA pipeline for WGS assemblies to determine number of mutated alleles
#' July 3 2024, Walter Demczuk & Shelley Peterson
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

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#Org_id <- "GONO"                   #PNEUMO or GONO
#SampleNo <- "55555"                #Sample No or "list"  
#curr_work_dir <- "C:\\WADE\\"
#-------------------------------------------------------------------------------

rRNA23S_pipeline <- function(Org_id, SampleNo, curr_work_dir) {
  
  start_time <- Sys.time()
  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  Test_id <- "rRNA23S"
  directorylist <- getdirectory(curr_work_dir, Org_id, Test_id)

  Variable <- NA
  #-----------------------------------------------------------------------------

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
  #-----------------------------------------------------------------------------
  if(SampleNo == "list")
  {
    SampleList.df <- read.csv(directorylist$SampleList, header = TRUE, sep = ",")
  }else
  {
    SampleList.df <- tibble(SampleNo, Variable)
  }
  NumSamples <- dim(SampleList.df)[1]

  #-----------------------------------------------------------------------------

  cat("\n\n23s rRNA counts from VCF\n")

  #Progress Bar
  n_iter <- NumSamples
  init <- numeric(n_iter)
  end <- numeric(n_iter)
  
  m<-1
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Start Sample Loop
  {
    # for progress bar
    init[m] <- Sys.time()
    
    # actual code
    A2059G <- 99L
    C2611T <- 99L
    CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
    QueryFile <- paste(directorylist$vcf_dir, "/", CurrSampleNo, ".vcf", sep = "")

    if(!file.exists(QueryFile))
    {
      A2059G <- NA
      C2611T <- NA
      QueryFile <- "File not found"
      cat("Sample Number ", CurrSampleNo, " not found!\n", sep = "")
    }

    if(QueryFile != "File not found")
    {
      allele_fraction <- 0L
      con <- file(QueryFile, open="r")
      linn <- readLines(con)
      close(con)
      for (i in 1:length(linn))
      {
        if(str_detect(linn[i], rRNA23S_position1))  #2597 & 2045 for GONO; 2061 for pneumo
        {
          vcf_parts <- strsplit(linn[i], ";")
          vcf_parts <- unlist(vcf_parts)
          DP1 <- vcf_parts[8]
          DP<-as.integer(substr(DP1, 4, 10))
          AO1 <- vcf_parts[6]
          AO<-as.integer(substr(AO1, 4, 10))

          allele_fraction <- as.integer((AO/DP)*100)
          if(allele_fraction < 14){A2059G <- 0L}
          if((allele_fraction >= 15) && (allele_fraction <= 34)) {A2059G <- 1L}
          if((allele_fraction >= 35) && (allele_fraction <= 64)) {A2059G <- 2L}
          if((allele_fraction >= 65) && (allele_fraction <= 84)) {A2059G <- 3L}
          if(allele_fraction >= 85) {A2059G <- 4L}
        }else
        {
          A2059G <- 0L
        }

        if(str_detect(linn[i], rRNA23S_position2))  #2597 & 2045 for GONO; 2061 for pneumo
        {
          vcf_parts2 <- strsplit(linn[i], ";")
          vcf_parts2 <- unlist(vcf_parts2)
          DP2_str <- vcf_parts2[8]
          DP2 <- as.integer(substr(DP2_str, 4, 10))   #total number of reads (depth of reads)
          AO2_str <- vcf_parts2[6]
          AO2 <- as.integer(substr(AO2_str, 4, 10))  #alternate observations from reference

          allele_fraction <- as.integer((AO2/DP2)*100)
          if(allele_fraction < 14){C2611T <- 0L}
          if((allele_fraction >= 15) && (allele_fraction <= 34)) {C2611T <- 1L}
          if((allele_fraction >= 35) && (allele_fraction <= 64)) {C2611T <- 2L}
          if((allele_fraction >= 65) && (allele_fraction <= 84)) {C2611T <- 3L}
          if(allele_fraction >= 85) {C2611T <- 4L}
        }else
        {
          C2611T <- 0L
        }
      }
    }

    cat(CurrSampleNo, "\tAC_2059: ", A2059G, "\tAC_2611: ", C2611T,"\n", sep = "")

    SampleProfile.df <- tibble(CurrSampleNo, A2059G, C2611T)
    OutputProfile.df <- rbind(OutputProfile.df, SampleProfile.df)
  
    #Progress Bar
    end[m] <- Sys.time()
    time <- round(seconds_to_period(sum(end - init)), 0)
    
    # Estimated remaining time based on the
    # mean time that took to run the previous iterations
    est <- NumSamples * (mean(end[end != 0] - init[init != 0])) - time
    remaining <- round(seconds_to_period(est), 0)
    
    cat(paste(" // Estimated time remaining:", remaining), "\n")
  } #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> End Sample Loop

  header <- c("SampleNo", "A2059G", "C2611T")
  names(OutputProfile.df) <- header
  write.csv(OutputProfile.df, directorylist$outfile, row.names = F)

  elapsed <- format(Sys.time() - start_time)
  
  cat("\n\nDone!\n\n\n", "Elapsed time: ", elapsed)

return(OutputProfile.df)
}



