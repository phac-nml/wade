#' emm typing pipeline from WGS assemblies
#' February 5 2024, Walter Demczuk & Shelley Peterson
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param SampleNo Sample number associated with contig.fasta file
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details GAS emm typing:
#' Takes Organism, Sample Number and queries a contig.fasta file.
#' For the most current emm types go to:
#' https://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/
#' download the alltrimmed.fta file
#' rename to "emm_trimmed.fasta"
#' and place in the allele_lkup_dna folder in WADE GAS/emm
#' then be sure to index it using MakeBlastDB
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
# Org_id <- "GAS"
# SampleNo <- "list"        # Sample No or "list"
# curr_work_dir <- "C:\\WADE\\"
# Blast_evalue <- "10e-50"         #sets sensitivity of Blast gene match 10e-50 to 10e-150; use 10e-5 for primers
#-------------------------------------------------------------------------------

EMM_pipeline <- function(Org_id, SampleNo, curr_work_dir){

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  Test_id <- "EMM"
  directorylist <- getdirectory(curr_work_dir, Org_id, Test_id)
  reflist <- refdirectory(directorylist, Org_id, Test_id)
  #-----------------------------------------------------------------------------

  Variable <- NA
  emmType <- ""
  emmSubtype <- ""
  emmTypeRep <- ""
  BP_ids <-""
  emmComment <- ""

  LocusLkupDNA <- paste0(reflist$Lkup_Dir, "emm_trimmed.fasta")
  CurrLocusLen <- 180L

  ########################### Setup sample list table ##########################
  if(SampleNo == "list")
  {
    SampleList.df <-as_tibble(read.csv(directorylist$SampleList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  }else
  {
    SampleList.df <- tibble(SampleNo, Variable)
  }
  
  NumSamples <- dim(SampleList.df)[1]

  bad_emm_file <- paste(reflist$Ref_Dir, "bad_emm_types.csv", sep = "")
  df.bad_emm <- as_tibble(read.csv(bad_emm_file, header = TRUE, sep = ",", stringsAsFactors = FALSE))

  cat("\n\n", "Sample No.", "\n", sep = "")

  m<-1L
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Loop through sample table
  {
    emmType <- ""
    emmSubtype <- ""
    emmTypeFinal <- ""
    emmComment <- ""
    
    CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
    QueryFile <- paste(directorylist$contigs_dir, CurrSampleNo, ".fasta", sep = "")
    
    if(!file.copy(QueryFile, reflist$DestFile, overwrite = T))
    {
      emmType <- "Sample_Err"
      emmSubtype <- "Sample_Err"
      emmTypeRep <- "Sample_Err"
      BP_ids <- 0L
      emmComment <- "Sample not found"
    }
    
    if(emmType != "Sample_Err") # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Look for emm if no sample error
    {
      Blast_Out_File <- blastquery(directorylist, reflist, "emm_trimmed", 10e-50)
      
      info = file.info(Blast_Out_File)
      if(info$size == 0)
      {
        emmType <- "NF"
        emmSubtype <- "NF"
        emmTypeRep <- "NF"
        BP_ids <- 0L
        emmComment <- "No emm gene"
      }else #=================================================================== Found EMM
      {
        df.blastout <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
        names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
        df.blastout <- distinct(df.blastout, Allele, .keep_all = TRUE) # sometimes there are multiple alleles

        #first find bad alleles
        df.blastout_bad <- inner_join(df.blastout, df.bad_emm, by = "Allele")
        df.blastout_bad_100 <- filter(df.blastout_bad, Ident == 100 & Align == CurrLocusLen)
        emmComment <- tolower(paste0(df.blastout_bad_100$Allele, collapse = "/"))

        #now find a match
        df.blastout_2 <- anti_join(df.blastout, df.bad_emm, by = "Allele")  #remove bad emm's
        dfSize_blastout_2 <- nrow(df.blastout_2)
        
        # if there are no good emm types even with partial match (blastout_2 is empty)
        # but only bad ones, use the bad one as the alternative emm type as if it were
        # a partial match
        if(dfSize_blastout_2 == 0)
        {
          df.blastout_2 <- df.blastout_bad_100
        }
        
        df.blastout_2$emmType <- sub("\\.\\d+","", df.blastout_2$Allele)
        
        df.blastout100 <- filter(df.blastout_2, Ident == 100 & Align == CurrLocusLen)
        
        if(nrow(df.blastout100 >= 1))
        {
          emmdf <- df.blastout100
        } else
        {
          emmdf <- dplyr::slice(df.blastout_2,1)
        }
        
        emmType <- tolower(paste0(emmdf$emmType, collapse = "/"))
        emmSubtype <- tolower(paste0(emmdf$Allele, collapse = "/"))
        BP_ids <- emmdf$Align[1] - emmdf$Mismatches[1]

        # assess quality of match
        if(BP_ids <= 149 | is.na(BP_ids))
        {
          emmType <- "NF"
          emmSubtype <- "NF"
          emmTypeRep <- "NF"
          emmComment <- tolower(df.blastout_2$Allele[1])
          if(is.na(BP_ids) & nrow(df.blastout_bad > 0))
          {
            emmComment <- "partial match to bad allele"
          }
        }else
        {
          emmTypeRep[BP_ids >= 150] <- "NT"
          emmTypeRep[BP_ids == 180] <- emmSubtype
        }
      }#======================================================================== End found EMM
    }#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End look for EMM

    Output.df <- tibble(CurrSampleNo, emmType, emmTypeRep, emmSubtype, BP_ids, emmComment)
    names(Output.df) <- c("Sample", "Type", "Subtype_Rep", "Subtype", "bp_id", "Comments")

    #if this is the first sample, copy output, else rowbind to add next sample profiles to output
    if(m==1)
    {
      SampleOutput.df <- tibble(Output.df)
    }else
    {
      SampleOutput.df <- bind_rows(SampleOutput.df, Output.df)
    }

    cat(CurrSampleNo, "\t", emmType, "\t", emmTypeRep, "\t", emmSubtype, "\t", BP_ids, "\t", emmComment, "\n", sep = "")

  } #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> end sample loop

  ##############################################################################
  # Export Results to csv file
  ##############################################################################
  write.csv(SampleOutput.df, paste0(directorylist$output_dir, "LabWareUpload_GAS_emm.csv"), quote = FALSE,  row.names = FALSE)
  write.csv(SampleOutput.df, paste0(directorylist$output_dir, "output_profile_emm.csv"), quote = FALSE, row.names = FALSE)
  cat("\n\nDone!\n\n\n")

  return(SampleOutput.df)

}

