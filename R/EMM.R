#' emm typing pipeline from WGS assemblies - October 30, 2020
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param SampleNo Sample number associated with contig.fasta file
#' @param locus Sample number associated with contig.fasta file
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details GAS emm typing:
#' Takes Organism, Sample Number, Locus, and a Variable at queries a contig.fasta file.
#' For the most current emm types go to:
#' https://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/
#' download the alltrimmed.fta file
#' rename to "emm_trimmed.fasta"
#' and place in the allele_lkup_dna folder in WADE GAS/emm_R
#' then be sure to index it using MakeBlastDB
#' @return A table frame containing the results of the query
#' @export

EMM_pipeline <- function(Org_id, SampleNo, locus, curr_work_dir){

  # get directory structure
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  Directories.df <- read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir

  Lkup_Dir <- paste(system_dir, "GAS\\emm_R\\allele_lkup_dna\\", sep = "")
  Tmp_Dir <- paste(system_dir, "GAS\\emm_R\\temp\\", sep = "")
  Loci_List <- paste(Tmp_Dir, "loci.csv", sep = "")

  Variable <- NA
  emmType <- ""
  emmSubtype <- ""
  emmTypeRep <- ""
  BP_ids <-""
  emmComment <- ""

  LocusLkupDNA <- file.path(Lkup_Dir, "emm_trimmed.fasta")
  CurrLocusLen <- 180L

  #------------------------------------------------------------------------------------Set up Sample list table
  if(SampleNo == "list")
  {
    SampleList.df <- read.csv(SampList, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    SampleList.df <- as_tibble(SampleList.df)
  }else
  {
    SampleList.df <- tibble(SampleNo, Variable)
  }

  Size.df <- dim(SampleList.df)
  NumSamples <- Size.df[1]

  bad_emm_file <- paste(Tmp_Dir, "bad_emm_types.csv", sep = "")
  df.bad_emm <- read.csv(bad_emm_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  df.bad_emm <- as_tibble(df.bad_emm)

  cat("\n\n", "Sample No.", "\n", sep = "")

  m<-1L
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Loop through sample table
  {
    emmType <- ""
    emmSubtype <- ""
    emmTypeFinal <- ""
    emmComment <- ""

    CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
    QueryFile <- paste(ContigsDir, CurrSampleNo, ".fasta", sep = "")
    DestFile <- paste(local_temp_dir, "queryfile.fasta", sep = "")
    if (!file.copy(QueryFile, DestFile, overwrite = T))
    {
      emmType <- "Sample_Err"
      emmSubtype <- "Sample_Err"
      emmTypeRep <- "Sample_Err"
      BP_ids <- 0L
      emmComment <- "Sample not found"
    }

      if (emmType != "Sample_Err")
      {
        if (!(file.exists(LocusLkupDNA)))
        {
          cat("emm_trimmed.fasta file not found!")
          sample_error.df <- tibble(Loucs_ID = locus, Output = "File Not Found")

          return(sample_error.df)
        }

        Blast_Out_File <- paste(local_temp_dir, "blastout.txt", sep = "")
        BlastCommand <- paste("blastn -query ", DestFile, " -db ", LocusLkupDNA, " -out ", Blast_Out_File, " -num_alignments 10 -evalue 10e-50 -outfmt 6")

        shell(BlastCommand)

        info = file.info(Blast_Out_File)
        if(info$size == 0)
        {
          emmType <- "NF"
          emmSubtype <- "NF"
          emmTypeRep <- "NF"
          BP_ids <- 0L
          emmComment <- "No emm gene"
        }else
        {
          df.blastout <- read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
          df.blastout <- as_tibble(df.blastout)
          names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
          df.blastout <- distinct(df.blastout, Allele, .keep_all = TRUE) # sometime there are multiple alleles

          #first find bad alleles
          df.blastout_bad <- inner_join(df.blastout, df.bad_emm, by = "Allele")
          df.blastout_bad <- filter(df.blastout_bad, Ident == 100 & Align == CurrLocusLen)
          dfSize_bad <- nrow(df.blastout_bad)

          if (dfSize_bad > 0)
          {
            n<-1
            for (n in 1L:dfSize_bad)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Loop through bad emms for comment entry
            {
              if (n==1)
              {
                emmComment <- tolower(df.blastout_bad$Allele[n])
              }else
              {
                emmComment <- paste(emmComment, "/", tolower(df.blastout_bad$Allele[n]), sep = "")
              }
            }
          }

          #now find a match
          df.blastout_2 <- anti_join(df.blastout, df.bad_emm, by = "Allele")  #remove bad emm's
          #2023-02-09 WD: If there is no good emm types even with partial match (blastout_2 is empty),
          #               but only bad ones, use the bad one as
          #               the alternative emm type as if it were a partial match.
          dfSize_blastout_2 <- nrow(df.blastout_2)
          if (dfSize_blastout_2 == 0)
          {
            df.blastout_2 <- df.blastout_bad
          }
          df.blastout100 <- filter(df.blastout_2, Ident == 100 & Align == CurrLocusLen)

          dfSize <- nrow(df.blastout100)
          if (dfSize >= 1L)
          {
            p<-1L
            for (p in 1L:dfSize)
            {
              if (p==1)
              {
                emmSubtype <- tolower(df.blastout100$Allele[p])
                emmParts <- unlist(strsplit(emmSubtype, "\\."))
                emmType <- tolower(emmParts[1])
                BP_ids <- 180
              }else
              {
                emmSubtype2 <- tolower(df.blastout100$Allele[p])
                emmSubtype <- paste(emmSubtype, "/", emmSubtype2, sep = "")
                emmParts <- unlist(strsplit(emmSubtype2, "\\."))
                emmType2 <- tolower(emmParts[1])
                emmType <- paste(emmType, "/", emmType2, sep = "")
                BP_ids <- 180
              }
            }
          }else
          {
            emmSubtype <- tolower(df.blastout_2$Allele[1])
            emmParts <- unlist(strsplit(emmSubtype, "\\."))
            emmType <- tolower(emmParts[1])
            BP_ids <- (df.blastout_2$Align[1] - df.blastout_2$Mismatches[1])
          }

          #------------------
          if (BP_ids == 180)
          {
            emmTypeRep <- emmSubtype
          } else
          {
            if (BP_ids >= 150)
            {
            emmTypeRep <- "NT"
            }
            if (BP_ids <= 149)
            {
            emmType <- "NF"
            emmSubtype <- "NF"
            emmTypeRep <- "NF"
            emmComment <- tolower(df.blastout_2$Allele[1])
            }
          }
          #------------------
        }#end if no emm gene found
      }#end if not sample error

      Output.df <- tibble(CurrSampleNo, emmType, emmTypeRep, emmSubtype, BP_ids, emmComment)
      names(Output.df) <- c("Sample", "Type", "Subtype_Rep", "Subtype", "bp_id", "Comments")

    #-------------------- if this is the first sample, copy output of first, else rowbind to add next sample profiles to output
    if(m==1)
    {
      SampleOutput.df <- tibble(Output.df)
    }else
    {
      SampleOutput.df <- bind_rows(SampleOutput.df, Output.df)
    }

    cat(CurrSampleNo, "\t", emmType, "\t", emmTypeRep, "\t", emmSubtype, "\t", BP_ids, "\t", emmComment, "\n", sep = "")

  } #end sample loop


  #============================================================  Write output files
  if (NumSamples == 1L)
  {
    write.csv(SampleOutput.df, paste(local_output_dir, "LabWareUpload_", Org_id, "_emm.csv", sep = ""), quote = FALSE,  row.names = FALSE)
    if (emmType != "Sample_Err")
    {
      write.csv(df.blastout, paste(local_output_dir, "output_profile_emm.csv", sep = ""), row.names = F)
      return(df.blastout)
    } else return(SampleOutput.df)

    cat("\n\nDone!\n\n\n")
  }else
  {
    write.csv(SampleOutput.df, paste(local_output_dir, "LabWareUpload_GAS_emm.csv", sep = ""), quote = FALSE,  row.names = FALSE)
    write.csv(SampleOutput.df, paste(local_output_dir, "output_profile_emm.csv", sep = ""), quote = FALSE, row.names = FALSE)
    cat("\n\nDone!\n\n\n")

    return(SampleOutput.df)
  }

} # end function call

