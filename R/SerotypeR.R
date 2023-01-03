#' Molecular typing pipeline for WGS assemblies
#'
#' Takes Organism, Sample Number, Locus, to query a contig.fasta file
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param Test_id AMR, TOXINS, VIRULENCE, NGSTAR, use MASTER for 16S id
#' @param SampleNo Sample number associated with contig.fasta file
#' @param LocusID The locus to query, or enter "list" to use a list of alleles
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details How it works:
#'
#' Pneumococcus serotyping based on PneumoCaT and SeroBA libraries.
#' Copies contig assembly from Warehouse drive
#' BLAST's copied assembly vs. reference data fasta of whole CPS regions for each serotype
#' If CPS region found that requires snp based analysis, queries each reference gene for that serogroup vs. loci list
#' Interprets snps to serotype level
#' The algorithm first determines the serogroup by blasting CPS_reference
#' Each locus in the locus list (temp folder) associated with that serogroup will then be blasted
#' There are 4 stages of locus analysis:
#' 1)result = POS/NEG (presence/abscence)
#' 2)pseudo = pseudogene (disrupted/intact)
#' 3)mutations = serotype determining amino acid substitutions as listed in locus_mutations
#' 4)allele = entire gene sequence match of conserved serotype determining genes as found in the allele_lookup folders
#'
#' The relevant result type is listed in the locus list table (temp folder)
#' As each result is evaluated for each relevant locus result, the serotype lookup table locus_lookup (temp folder)
#' is filtered
#' If a single row is left in the lookup table by the end of the sample, that is the serotype
#' Otherwise a Fail<serogroup> answer will be returned
#' Supporting files like blast outputs, results, extracted fasta sequences can be found
#' in the user's local output or temp folders
#'
#' @return A table frame containing the results of the query
#' @export

###for testing
#curr_work_dir <- "C:\\WADE\\"
#Org_id <- "PNEUMO"
#SampleNo <- "list"
#LocusID <- "list"

SEROTYPE_pipeline <- function(Org_id, SampleNo, LocusID, Test_id, curr_work_dir){

  #--------------------------------------------------------------SET VARIABLES MANUALLY FOR DEBUGGING

  Variable <- NA
  Run_MakeBlastDB <- FALSE

  #------------------------------------------------------------------------------------------------------------
  # get directory structure
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  Directories.df <- read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "\\", "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "\\Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "\\temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  Contigs_Dir <- Directories_org.df$ContigsDir
  Main_Dir <- paste(system_dir, Org_id, "\\Serotype_R\\", sep = "")
  Lkup_Dir <- paste(Main_Dir, "allele_lkup_dna\\", sep = "")
  Tmp_Dir <- paste(Main_Dir, "temp\\", sep = "")
  Wild_Dir <- paste(Main_Dir, "wildgenes\\", sep = "")
  LocusLkupDNA_CPS <- paste(Lkup_Dir, "reference_CPS.fasta", sep = "")
  CPSlist <- paste(Tmp_Dir, "CPS_types_for_blast.csv", sep = "")
  CPS_types.df <- as_tibble(read.csv(CPSlist, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  CPS_types <- CPS_types.df$Serogroup
  #------------------------------------------------------------------------------------------

  #------------------------------------
  # if anything changed in reference_CPS.fasta, need to re-index using 2 lines below
  #FormatCommand <- paste("makeblastdb -in ", LocusLkupDNA_CPS, " -dbtype nucl", sep = "")
  #shell(FormatCommand, intern = TRUE)
  #------------------------------------

  evalueList <- paste(Tmp_Dir, "blast_evalues.csv", sep = "")
  blast_evalues.df <-  read.csv(evalueList, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  Blast_evalue <- as.character( blast_evalues.df$contig[1])
  evalue_CPS <- as.character(blast_evalues.df$cps[1])
  evalue_allele <- as.character( blast_evalues.df$allele[1])

  outfile_nf <- paste(local_output_dir, "output_dna_notfound.fasta", sep = "")
  unlink(outfile_nf) #this deletes the file!

  outfile <- paste(local_output_dir, "LabWareUpload_", Org_id, "_SEROTYPE.csv", sep = "")
  unlink(outfile)

  if(SampleNo == "list")
  {
    SampleList.df <- as_tibble(read.csv(SampList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  }else
  {
    SampleList.df <- tibble(SampleNo, Variable)
  }

  Size.df <- dim(SampleList.df)
  NumSamples <- Size.df[1]

  m <- 1L
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  {
    CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
    CurrSampleVar <-as.character(SampleList.df[m, "Variable"])

    SampleProfile <- ""
    serogroup <- ""
    serotype <- ""
    OutputLocus.df <- tibble(SampleNo = as.character(CurrSampleNo))  #all results

    cat("\n\n", CurrSampleNo, ": ", m, " of ", NumSamples, "\n", sep = "")

    #----------------------------  Get serogroup CPS locus for sample first
    ContigFile <- paste(Contigs_Dir, CurrSampleNo, ".fasta", sep = "")
    ContigFileLocal <- paste(local_temp_dir, "queryfile.fasta", sep="")
    if (!file.copy(ContigFile, ContigFileLocal, overwrite = T))
    {SampleFound <- FALSE
     serogroup <- "Sample_Err"
     serotype <- "Sample_Err"
     serotype1 <- "Sample_Err"
     SampleProfile <- "Sample_Err"
    }else
      {SampleFound <- TRUE
      #makeblastdb from contig file for later when blasting wildgenes vs. contig to extract genes
      FormatCommand <- paste("makeblastdb -in ", ContigFileLocal, " -dbtype nucl", sep = "")
      shell(FormatCommand, intern = TRUE)
      #use contig file to blast against the CPS data to see what the CPS group is
      Blast_Out_File <- paste(local_temp_dir, "blastout_serogroup.txt", sep = "")
      BlastCommand <- paste("blastn -query ", ContigFileLocal, " -db ", LocusLkupDNA_CPS, " -out ", Blast_Out_File, " -num_alignments 10 ", "-evalue ", evalue_CPS,  " -outfmt 6")
      shell(BlastCommand)
      info = file.info(Blast_Out_File)
      if(info$size == 0)  #no blast result from CPS locus lookup i.e. not pneumo, bad sequencing
      {
        serotype <- "Unknown"
        serotype1 <- "unknown"
        serogroup <- "unknown"
      }else
      {
        df.blastout <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
        names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
        df.blastout <- arrange(df.blastout, desc(bit))
        serotype1 <- df.blastout$Allele[1]

        SerotypeParts <- unlist(strsplit(serotype1, "_"))
        SerotypeParts <- SerotypeParts[1]
        serotype <- ifelse(substr(SerotypeParts, 1, 1) == "0", sub("^.", "", SerotypeParts), SerotypeParts)

        if (Org_id == "PNEUMO")
        {
          serogroup <- as.character(substr(SerotypeParts, 1, 2))
        }
        if (Org_id == "GBS")
        {
          serogroup <- "GBS"
        }
        outfile_sero <- paste(local_output_dir, "output_blastout_serogroup.csv", sep = "")
        write.csv(df.blastout, outfile_sero, quote = FALSE, row.names = F )
      }
    }#end Sample Found

    cat("CPS Blast Serogroup ", serogroup, " - Serotype ", serotype, " !\n")

    if (serogroup == "unknown" | serogroup == "Sample_Err")  # no CPS match or no contig
    {
      serotype <- paste("NT[", serogroup, "]", sep = "")
      SampleProfile <- serogroup
    }else  # valid serogroup found so setup loci and get locuslist screen for serotypes requiring no snp analysis
      if (serotype %in% CPS_types)
      {
       SampleProfile <- "CPS operon type match"
       locus <- ""
       head(df.blastout, n = 5L)

      }else
      { SampleProfile <- ""
        LocusList <- paste(Tmp_Dir, serogroup, "_loci.csv", sep = "")
        if (file.exists(LocusList))
          {LocusListPresent <- TRUE
          LocusList.df <- read.csv(LocusList, header = TRUE, sep = ",", stringsAsFactors = FALSE)

          if(LocusID != "list")
            {LocusList.df <- filter(LocusList.df, Locus_id == LocusID)}
          SizeList <- dim(LocusList.df)
          NumLoci <- SizeList[1]
        }else {
              LocusListPresent <- FALSE
              NumLoci <- 0
              }

        LocusMutationsFile <- paste(Tmp_Dir, serogroup, "_loci_mutations.csv", sep = "")
        if(file.exists(LocusMutationsFile))
        {Mutns <- "Yes"
        LocusMutations.df <- read.csv(LocusMutationsFile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
        } else {Mutns <- "No"}

        SerotypeLookupFile <- paste(Tmp_Dir, serogroup, "_loci_lookup.csv", sep = "")
        if(file.exists(SerotypeLookupFile))
        {SeroLookupPresent <- "Yes"
        SerotypeLookup.df <- read.csv(SerotypeLookupFile, header = TRUE, sep = ",", stringsAsFactors = FALSE)

        if (Org_id == "PNEUMO")
        {
        headers <- names(SerotypeLookup.df) #need to remove the first X from each header name
        headers2 <- substr(headers, 2, 50 )
        headers2[1]<-"Serotype"
        names(SerotypeLookup.df) <- headers2

        headers2 <- substr(headers, 2, 50 )
        headers2[1]<-"Serotype"
        }

        } else {SeroLookupPresent <- "No"}

        if ((Run_MakeBlastDB) && (LocusListPresent))
        {
          for(q in 1L:NumLoci)
          {
            locus <- as.character(LocusList.df[q,1])
            LocusLkupDNA <- paste(Lkup_Dir, locus, ".fasta", sep = "")
            if(file.exists(LocusLkupDNA))
            {
              BlastFormatCommand <- paste("makeblastdb -in ", LocusLkupDNA, " -dbtype nucl", sep = "")
              try(system(BlastFormatCommand))
            }
          }
        }

    #XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  Loci List

    OutputLocus.df <- tibble(SampleNo = as.character(CurrSampleNo))  #all results
    OutputLocusProfile.df <- tibble(SampleNo = as.character(CurrSampleNo)) #for serotype
    p <- 1L
    for(p in 1L:NumLoci)
    {
      ifelse(SampleProfile == "", ProfileSeparator <- "", ProfileSeparator<-";")
      locus <- as.character(LocusList.df$Locus_id[p])

      LocusFile <- paste(Wild_Dir, locus, ".fasta", sep = "")
      if(!file.exists(LocusFile))
      {
        sample_error.df <- tibble(Loucs_ID = LocusID, Output = "Reference wildtype file not found.")
        return(sample_error.df)
      }

      locus2 <- gsub("GBS_", "", locus) #makes profile smaller for GBS LabWare output
      locus_result_type <- as.character(LocusList.df$Result_type[p])

      LocusLkupDNA <- paste(Lkup_Dir, locus, ".fasta", sep = "")
      ifelse(file.exists(LocusLkupDNA), LocusLkupDNApresent <- TRUE, LocusLkupDNApresent <- FALSE)

      col2_name <- paste(locus, "_result", sep = "")
      col3_name <- paste(locus, "_allele", sep = "")
      col4_name <- paste(locus, "_mutations", sep = "")
      col5_name <- paste(locus, "_pseudo", sep = "")

      OutputLocus.df[1, col2_name] <- NA
      OutputLocus.df[1, col3_name] <- NA
      OutputLocus.df[1, col4_name] <- NA
      OutputLocus.df[1, col5_name] <- NA

      LocusFile <- paste(Wild_Dir, locus, ".fasta", sep = "")
      if(!file.exists(LocusFile))
      {
        OutputLocus.df[1, col2_name] <- "Locus_Err"  #result
        OutputLocus.df[1, col3_name] <- "Locus_Err"  #allele
        OutputLocus.df[1, col4_name] <- "Locus_Err"  #mutations
        OutputLocus.df[1, col5_name] <- "Locus_Err"  #pseudo
        WildTypeLocusPresent <- FALSE
      }else   #only do the rest of the code if there is a wild type gene
        {
        WildTypeLocusPresent <- TRUE

        BlastCommand <- paste("blastn -query ", LocusFile, " -db ", ContigFileLocal, " -out ", local_temp_dir, "blastout.txt -evalue ", Blast_evalue, sep = "")

        shell(BlastCommand, intern = TRUE)

        Blastout <- paste(local_temp_dir, "blastout.txt", sep = "")
        con <- file(Blastout, open="r")
        linn <- readLines(con)
        close(con)

        #check if gene was found in BLAST
        BlastResult <- NA
        for (x in 1:length(linn))
        {
          if (str_detect(linn[x], "No hits found"))
          {
            BlastResult <- "NEG"
            break()   #break out of blastout scan when line found with no hits
          } else
          {
            BlastResult <- "POS"
          }
        }
        OutputLocus.df[col2_name] <- BlastResult  #result

        if (locus_result_type == "result")
        {
          OutputLocusProfile.df[col2_name] <- BlastResult
          locus2 <- gsub("GBS_", "", locus)
          SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[", BlastResult, "]", sep="")
          cat(col2_name, BlastResult, "\n", sep = "\t\t")
          if (LocusID == "list"){SerotypeLookup.df <- filter(SerotypeLookup.df, (  ((!!sym(col2_name)) == BlastResult) | (is.na(!!sym(col2_name)) ))  )}
        }
        #============================================= add missing genes to profile for mutations and pseudogenes
        if ((locus_result_type == "pseudo") & (BlastResult == "NEG"))
        {
          SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[no pseudogene]", sep="")
          cat(col5_name, "no gene present", "\n", sep = "\t\t")
        }
        if ((locus_result_type == "mutations") & (BlastResult == "NEG"))
        {
          SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[no mutn gene]", sep="")
          cat(col4_name, "no gene present", "\n", sep = "\t\t")
        }
        if ((locus_result_type == "allele") & (BlastResult == "NEG"))
        {
          SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[no allele gene]", sep="")
          cat(col3_name, "no gene present", "\n", sep = "\t\t")
        }

        #==========================================================

        DNASeqLine_str <- ""
        WTDNASeqLine_str <- ""
        TimesThrough <- 0
        i<-1L
        if (BlastResult == "POS")
        {
          for (i in 1:length(linn))
          {
            if (str_detect(linn[i], "Score =")) #set counter to parse only first align block of blastout
            {
              TimesThrough <- TimesThrough + 1
            }

            if (TimesThrough == 1)
            {
              if (str_detect(linn[i], "Query "))
              {
                QueryLine <-  unlist(strsplit(linn[i], " "))
                QueryLine <- QueryLine[QueryLine != ""]
                WTDNASeqLine_str <- paste(WTDNASeqLine_str, QueryLine[3], sep = "")
              }

              if (str_detect(linn[i], "Sbjct "))
              {
                SbjctLine <-  unlist(strsplit(linn[i], " "))
                SbjctLine <- SbjctLine[SbjctLine != ""]
                DNASeqLine_str <- paste(DNASeqLine_str, SbjctLine[3], sep = "")
              }
            }
          }

          WTDNASeqLine <- DNAString(WTDNASeqLine_str)
          WTDNASeqLine_NoDash_str <- str_replace_all(WTDNASeqLine_str, "-", "")
          WTDNASeqLine_NoDash <- DNAString(WTDNASeqLine_NoDash_str)

          DNASeqLine <- DNAString(DNASeqLine_str)
          DNASeqLine_NoDash_str1 <- str_replace_all(DNASeqLine_str, "-", "")
          DNASeqLine_NoDash_str <- str_replace_all(DNASeqLine_NoDash_str1, "N", "")
          DNASeqLine_NoDash <- DNAString(DNASeqLine_NoDash_str)

          Seq_File <- paste(local_temp_dir, "querygene.fasta", sep="")
          sink(Seq_File, split=FALSE, append = FALSE)
          cat(">", CurrSampleNo, "_", locus , "\n", DNASeqLine_NoDash_str, sep ="")
          sink()

          if (LocusID != "list")
          {
            cat("\n\n>", locus , "(Wildtype)\n", WTDNASeqLine_str, "\n", sep ="")
            cat("\n\n>", CurrSampleNo, "_", locus , "\n", DNASeqLine_str, "\n", sep ="")
          }

          #-------------------------------------------------------------------------------make Protein sequence
          WTAASeqLine <- suppressWarnings(translate(WTDNASeqLine_NoDash))
          WTAASeqLine_str <- toString(WTAASeqLine)
          AASeqLine <- suppressWarnings(translate(DNASeqLine_NoDash))
          AASeqLine_str <- toString(AASeqLine)

          AASeqLine_length <- str_length(AASeqLine_str)
          AAseqLine_lastchr <- substring(AASeqLine_str, AASeqLine_length, AASeqLine_length)
          if (AAseqLine_lastchr == "*"){AASeqLine_str <- substring(AASeqLine_str, 1L, AASeqLine_length-1)}

          ifelse (str_detect(AASeqLine_str, "[*]"), nonsence_mutation <- "Disrupted", nonsence_mutation <- "Intact") #used for wcjE_11E

          #update locus ouput
          OutputLocus.df[col5_name] <- nonsence_mutation  #pseudo
          if (locus_result_type == "pseudo")
          {
            OutputLocusProfile.df[col5_name] <- nonsence_mutation
            SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[", nonsence_mutation, "]", sep="")
            cat(col5_name, nonsence_mutation, "\n", sep = "\t\t")
            if (LocusID == "list"){SerotypeLookup.df <- filter(SerotypeLookup.df, (  (!!sym(col5_name)) == nonsence_mutation | is.na(!!sym(col5_name))  ))}
          }

          if (LocusID != "list")
          {
            cat("\n\n>", locus , "(Wildtype)\n", WTAASeqLine_str, "\n", sep ="")
            cat("\n\n>", CurrSampleNo, "_", locus , "\n", AASeqLine_str, "\n", sep ="")
          }

          #-------------------------------------------------------------------------------------------mutations
          motifs <- ""
          if (Mutns == "Yes")
          {
            LocusMutationsCurr.df <- filter(LocusMutations.df, Locus_id == locus)
            SizeMutns.df <- dim(LocusMutationsCurr.df)
            NumMutations <- SizeMutns.df[1]
            aa_mut_name <- ""
            aa_start <- ""
            aa_end <- ""
            aa_mut <- ""
            aa_value <- ""
            motif<-""
            motifs <- ""

            if (NumMutations > 0 )
            {
              for(w in 1L:NumMutations)
              {
                aa_mut_name <- LocusMutationsCurr.df[w,"Name"]
                aa_start <- LocusMutationsCurr.df[w,"Posn_1"]
                aa_end <- LocusMutationsCurr.df[w,"Posn_2"]
                aa_mut <- LocusMutationsCurr.df[w,"Mutation"]
                aa_value<-substr(AASeqLine_str, aa_start, aa_end)
                if (motifs=="") {motif_sep <- ""}else{motif_sep <- "/"}

                 if (aa_mut == aa_value)
                 {motif<- paste(motif_sep, aa_start, aa_mut, sep = "")}else
                 {motif<-""}
                motifs <- paste(motifs, motif, sep = "")
              } #endfor
            }
            if (LocusID != "list"){cat("\nMotif: ", motifs, "\n\n", sep = "")}
          } #endif mutations found

          OutputLocus.df[col4_name] <- motifs  #mutations
          if (locus_result_type == "mutations")
          {
            OutputLocusProfile.df[col4_name] <- motifs  #mutations
            SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[", motifs, "]", sep="")
            cat(col4_name, motifs, "\n", sep = "\t\t")
            if (LocusID == "list")
            {
              SerotypeLookup.df <- filter(SerotypeLookup.df, ((!!sym(col4_name)) == motifs | is.na(!!sym(col4_name))))
            }
          }

          #----------------------------------------------------------------------------------Lookup Alleles DNA
          ExactMatchFound <- FALSE

          if(LocusLkupDNApresent)
          {
            #-----------------
            Blast_Out_File <- paste(local_temp_dir, "blastout2.txt", sep = "")
            BlastCommand <- paste("blastn -query ", Seq_File, " -db ", LocusLkupDNA, " -out ", Blast_Out_File,
                                  " -num_alignments 10 -evalue ", evalue_allele, " -outfmt 6")
            shell(BlastCommand)
            info = file.info(Blast_Out_File)
            if(info$size == 0)  #no blast result
            {
              Allele <- NA
            }else
            {
              df.blastout2 <- read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
              names(df.blastout2) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
              df.blastout <- arrange(df.blastout2, desc(bit))
              outfile_allele <- paste(local_temp_dir, "blastout_lkup_allele.csv", sep = "")
              write.csv(df.blastout2, outfile_allele, quote = FALSE, row.names = F )
              Allele <- df.blastout2[1, "Allele"]
              AlleleParts <- unlist(strsplit(Allele, "_"))


              Allele2 <- AlleleParts[2]

              if (df.blastout2$Ident[1] >= 97)
              {Allele <- Allele2}else
              {Allele<- "NF" #paste("NF(", Allele2, ")", sep = "")
                outfile_nf <- paste(local_output_dir, "output_dna_notfound.fasta", sep = "")
                sink(outfile_nf, split=FALSE, append = TRUE)
                cat(">", locus, "_", CurrSampleNo, "_", "\n", DNASeqLine_NoDash_str, "\n", sep ="")
                sink()
              }
            }

            OutputLocus.df[col3_name] <- Allele
            if (locus_result_type == "allele")
            {
              OutputLocusProfile.df[col3_name] <- Allele  #allele
              SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[", Allele, "]", sep="")
              cat(col3_name, Allele,"\n", sep = "\t\t")
              if (LocusID == "list") {SerotypeLookup.df <- filter(SerotypeLookup.df, ( ( ((!!sym(col3_name)) == Allele) |  (is.na(!!sym(col3_name))) ) )    ) }
            }

          }#close bracket for DNA lookup file exists check

        }# end BLAST positive

      } #close bracket for wild type gene file found.^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    } #end of locus list loop xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

    # Interpret molecular characterization to get serotype
    if (LocusID == "list")
    {
    SizeSerotypeResult <- dim(SerotypeLookup.df)
    NumResults <- SizeSerotypeResult[1]
    if (NumResults == 0)
    {
      serotype <- paste("NT[", serotype, "]", sep = "")
    }else
    {
      for(w in 1L:NumResults)
      {
        if (w == 1)
        {serotype <-  SerotypeLookup.df$Serotype[w]}else
        {serotype <- paste(serotype, "/", SerotypeLookup.df$Serotype[w], sep = "")}
      }
    }
    }

      }     #end not sample error

    if (m==1L)
    {
      SampleOutputProfile.df <- tibble(SampleNo = CurrSampleNo, Serotype = serotype, Profile = SampleProfile)

    }else
    {
      SampleOutputProfile2.df <- tibble(SampleNo = CurrSampleNo, Serotype = serotype, Profile = SampleProfile )
      SampleOutputProfile.df<- bind_rows(SampleOutputProfile.df, SampleOutputProfile2.df)
    }

    cat("\nProfile:  ", SampleProfile, sep = "")
    cat("\nSerotype: ", serotype, sep = "")

  } #close brack for sample list loop<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  if (LocusID == "list")
  {
  write.csv(SampleOutputProfile.df, outfile, quote = FALSE,  row.names = F)
  }

  if (NumSamples == 1L)
  {
    outfile_sero2 <- paste(local_output_dir, "output_profile_SEROTYPE.csv", sep = "")
    write.csv(OutputLocus.df, outfile_sero2, quote = FALSE, row.names = F )
  }

  cat("\n\nDone!\n\n\n")

  return(SampleOutputProfile.df)

}
