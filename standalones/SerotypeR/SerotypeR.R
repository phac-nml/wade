##################################################################
####                       SerotypeR                          ####
####  WGS Analysis and Detection of Molecular Markers (WADE)  ####
####                 Author: Walter Demczuk                   ####
####                    Date: 2022-11-08                      ####
##################################################################

# Load libraries####
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyselect)
library(stringr)
library(Biostrings)

# USER INPUT--------------------------------------------------------------------------------

Org_id <- "PNEUMO"              # GBS or PNEUMO
SampleNo <- "12345"             # Single sample no (without the ".fasta" extension) or "list" for list.csv  
Variable <- NA                  # Extra optional info to be displayed in outputs (Quellung, metadata)
curr_work_dir <- "C:\\SerotypeR\\"
Contigs_Dir <- "C:\\SerotypeR\\contigs\\"

#-------------------------------------------------------------------------------------------
# get directory structure + set location of directory for lookups, mapping tables
local_output_dir <- paste(curr_work_dir, "output\\", sep = "")
local_ref_dir <- paste(curr_work_dir, "reference\\", sep = "")
local_temp_dir <- paste(curr_work_dir, "temp\\", sep = "")
Lkup_Dir <- paste(curr_work_dir, "allele_lkup_dna\\", sep = "")
Wild_Dir <- paste( curr_work_dir, "wildgenes\\", sep = "")

# FILES-------------------------------------------------------------------------------------

LocusLkupDNA_CPS <- paste(Lkup_Dir, "reference_CPS.fasta", sep = "")
SampList <- paste(curr_work_dir, "list.csv", sep = "")
  
#-------------------------------------------------------------------------------------------
evalueList <- paste(local_ref_dir, "blast_evalues.csv", sep = "")
blast_evalues.df <-  as_tibble(read.csv(evalueList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
Blast_evalue <- as.character( blast_evalues.df$contig[1])
evalue_CPS <- as.character(blast_evalues.df$cps[1])
evalue_allele <- as.character( blast_evalues.df$allele[1])

CPS_types <- c("1", "2", "3", "4", "5", "8", "13", "14", "20", "21", "24", "27",
               "29", "31", "34", "36", "37", "38", "39", "40", "16F", "16A", "17A",
               "17F", "10D", "32A", "32F", "Swiss_NT", "Swiss(NT)", "Alternative-aliB(NT)")

# set output files - need to delete previous outputs or errors when left open by user.
outfile_nf <- paste(local_output_dir, "output_dna_notfound.fasta", sep = "")
unlink(outfile_nf) #this deletes the file!
outfile <- paste(local_output_dir, "output_", Org_id, "_SEROTYPE.csv", sep = "")
unlink(outfile)

#-----------------------------------------------------------------------------------------------

# Index the BLAST lookup database
LocusLkupDNA_CPS <- paste(curr_work_dir, "allele_lkup_dna\\reference_CPS.fasta", sep = "")
BlastFormatCommand <- paste("makeblastdb -in ", LocusLkupDNA_CPS, " -dbtype nucl", sep = "")
try(system(BlastFormatCommand))

#-------------------------------------------------------------------------------------------
if(SampleNo == "list")
{
  SampleList.df <- as_tibble(read.csv(SampList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
}else
{
  SampleList.df <- tibble(SampleNo, Variable)
}

names(SampleList.df) <- c("SampleNo", "Variable")
Size.df <- dim(SampleList.df)
NumSamples <- Size.df[1]

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Loop through samples
m <- 1L
for(m in 1L:NumSamples) 
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
  if(!file.copy(ContigFile, ContigFileLocal, overwrite = T))
  {
    SampleFound <- FALSE
    serogroup <- "Sample_Err"
    serotype <- "Sample_Err"
    serotype1 <- "Sample_Err"
    SampleProfile <- "Sample_Err"
  }else
  {
    SampleFound <- TRUE
    
    #makeblastdb from contig file for later when blasting wildgenes vs. contig to extract genes
    FormatCommand <- paste("makeblastdb -in ", ContigFileLocal, " -dbtype nucl", sep = "")
    shell(FormatCommand, intern = TRUE)
    
    #use contig file to blast against the CPS data to see what the CPS group is
    Blast_Out_File <- paste(local_temp_dir, "blastout_serogroup.txt", sep = "")
    BlastCommand <- paste("blastn -query ", ContigFileLocal, " -db ", LocusLkupDNA_CPS, " -out ", Blast_Out_File, " -num_alignments 10 ", "-evalue ", evalue_CPS,  " -outfmt 6")
    shell(BlastCommand)
    info = file.info(Blast_Out_File)
    if(info$size == 0)  #no blast result from CPS locus lookup... not pneumo or GBS, bad sequencing
    {
      serotype <- "Unknown"
      serotype1 <- "Unknown"
      serogroup <- "Unknown"
    }else
    {
      df.blastout <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
      names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
      df.blastout <- arrange(df.blastout, -bit)
      serotype1 <- df.blastout$Allele[1]

      SerotypeParts <- unlist(strsplit(serotype1, "_"))
      # remove zeros from start of serotype: ^ means first character, . means anything
      serotype <- ifelse(substr(SerotypeParts, 1, 1) == "0", sub("^.", "", SerotypeParts), SerotypeParts)

      if(Org_id == "PNEUMO")
      {
        serogroup <- as.character(substr(SerotypeParts, 1, 2))
      }
      if(Org_id == "GBS")
      {
        serogroup <- "GBS"
      }
      outfile_sero <- paste(local_output_dir, "output_blastout_serogroup.csv", sep = "")
      write.csv(df.blastout, outfile_sero, quote = FALSE, row.names = F )
    }
  }#end Sample Found

  cat("CPS Blast Serogroup ", serogroup, " - Serotype ", serotype, " !\n")
  
  if(serogroup == "unknown" | serogroup == "Sample_Err")  # no CPS match or no contig
  {
    serotype <- paste("NT[", serogroup, "]", sep = "")
    SampleProfile <- serogroup
  }else  # valid serogroup found so setup loci and get locuslist screen for serotypes requiring no SNP analysis
  if(serotype %in% CPS_types)
  {
    SampleProfile <- "CPS operon type match"
    locus <- ""
    head(df.blastout, n = 5L)
  }else
  { 
    SampleProfile <- ""
    LocusList <- paste(local_ref_dir, serogroup, "_loci.csv", sep = "")
    if(file.exists(LocusList))
    {
      LocusListPresent <- TRUE
      LocusList.df <- as_tibble(read.csv(LocusList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
      SizeList <- dim(LocusList.df)
      NumLoci <- SizeList[1]
    }else 
    {
      LocusListPresent <- FALSE
      NumLoci <- 0
    }

    LocusMutationsFile <- paste(local_ref_dir, serogroup, "_loci_mutations.csv", sep = "")
    if(file.exists(LocusMutationsFile))
    {
      Mutns <- "Yes"
      LocusMutations.df <- as_tibble(read.csv(LocusMutationsFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
    }else 
    {
      Mutns <- "No"
    }

    SerotypeLookupFile <- paste(local_ref_dir, serogroup, "_loci_lookup.csv", sep = "")
    if(file.exists(SerotypeLookupFile))
    {
      SeroLookupPresent <- "Yes"
      SerotypeLookup.df <- as_tibble(read.csv(SerotypeLookupFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))

      if(Org_id == "PNEUMO")
      {
        headers <- names(SerotypeLookup.df) #need to remove the first X from each header name
        headers2 <- substr(headers, 2, 50 )
        headers2[1]<-"Serotype"
        names(SerotypeLookup.df) <- headers2

        headers2 <- substr(headers, 2, 50 )
        headers2[1]<-"Serotype"
      }
    }else 
    {
      SeroLookupPresent <- "No"
    }

    #XXXXXXXXXXXXXXXXXXXXXXxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx  Loci List
    OutputLocus.df <- tibble(SampleNo = as.character(CurrSampleNo))  #all results
    OutputLocusProfile.df <- tibble(SampleNo = as.character(CurrSampleNo)) #for serotype
    p <- 1L
    for(p in 1L:NumLoci)
    {
      ifelse(SampleProfile == "", ProfileSeparator <- "", ProfileSeparator<-";")
      locus <- as.character(LocusList.df$Locus_id[p])
      locus2 <- gsub("GBS_", "", locus) #makes profile smaller for GBS LabWare output
      locus_result_type <- as.character(LocusList.df$Result_type[p])

      LocusLkupDNA <- paste(Lkup_Dir, locus, ".fasta", sep = "")
      ifelse(file.exists(LocusLkupDNA), LocusLkupDNApresent <- TRUE, LocusLkupDNApresent <- FALSE)

      col2_name <- paste(locus, "_result", sep = "")
      col3_name <- paste(locus, "_allele", sep = "")
      col4_name <- paste(locus, "_mutations", sep = "")
      col5_name <- paste(locus, "_pseudo", sep = "")

      OutputLocus.df[1, col2_name] <- NA_character_
      OutputLocus.df[1, col3_name] <- NA_character_
      OutputLocus.df[1, col4_name] <- NA_character_
      OutputLocus.df[1, col5_name] <- NA_character_

      LocusFile <- paste(Wild_Dir, locus, ".fasta", sep = "")
      if(!file.exists(LocusFile))
      {
        OutputLocus.df[1, col2_name] <- "Locus_Err"  #result
        OutputLocus.df[1, col3_name] <- "Locus_Err"  #allele
        OutputLocus.df[1, col4_name] <- "Locus_Err"  #mutations
        OutputLocus.df[1, col5_name] <- "Locus_Err"  #pseudo
        WildTypeLocusPresent <- FALSE
        cat(locus, "\tWildtype reference file not found.", "\n", sep = "\t\t")
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
        for(x in 1:length(linn))
        {
          if(str_detect(linn[x], "No hits found"))
          {
            BlastResult <- "NEG"
            break()   #break out of blastout scan when line found with no hits
          }else
          {
            BlastResult <- "POS"
          }
        }
        OutputLocus.df[col2_name] <- BlastResult  #result

        if(locus_result_type == "result")
        {
          OutputLocusProfile.df[col2_name] <- BlastResult
          locus2 <- gsub("GBS_", "", locus)
          SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[", BlastResult, "]", sep="")
          cat(col2_name, BlastResult, "\n", sep = "\t\t")
          SerotypeLookup.df <- filter(SerotypeLookup.df, (  ((!!sym(col2_name)) == BlastResult) | (is.na(!!sym(col2_name)) ))  )
        }
        #====================================================================== add missing genes to profile for mutations and pseudogenes
        if((locus_result_type == "pseudo") & (BlastResult == "NEG"))
        {
          SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[no pseudogene]", sep="")
          cat(col5_name, "no gene present", "\n", sep = "\t\t")
        }
        if((locus_result_type == "mutations") & (BlastResult == "NEG"))
        {
          SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[no mutn gene]", sep="")
          cat(col4_name, "no gene present", "\n", sep = "\t\t")
        }
        if((locus_result_type == "allele") & (BlastResult == "NEG"))
        {
          SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[no allele gene]", sep="")
          cat(col3_name, "no gene present", "\n", sep = "\t\t")
        }

        #======================================================================

        DNASeqLine_str <- ""
        WTDNASeqLine_str <- ""
        TimesThrough <- 0
        i<-1L
        if(BlastResult == "POS")
        {
          for(i in 1:length(linn))
          {
            if(str_detect(linn[i], "Score =")) #set counter to parse only first align block of blastout
            {
              TimesThrough <- TimesThrough + 1
            }

            if(TimesThrough == 1)
            {
              if(str_detect(linn[i], "Query "))
              {
                QueryLine <-  unlist(strsplit(linn[i], " "))
                QueryLine <- QueryLine[QueryLine != ""]
                WTDNASeqLine_str <- paste(WTDNASeqLine_str, QueryLine[3], sep = "")
              }

              if(str_detect(linn[i], "Sbjct "))
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

          #-------------------------------------------------------------------- Make Protein sequence
          WTAASeqLine <- suppressWarnings(translate(WTDNASeqLine_NoDash))
          WTAASeqLine_str <- toString(WTAASeqLine)
          AASeqLine <- suppressWarnings(translate(DNASeqLine_NoDash))
          AASeqLine_str <- toString(AASeqLine)

          AASeqLine_length <- str_length(AASeqLine_str)
          AAseqLine_lastchr <- substring(AASeqLine_str, AASeqLine_length, AASeqLine_length)
          if(AAseqLine_lastchr == "*")
          {
            AASeqLine_str <- substring(AASeqLine_str, 1L, AASeqLine_length-1)
          }

          ifelse (str_detect(AASeqLine_str, "[*]"), nonsence_mutation <- "Disrupted", nonsence_mutation <- "Intact") #used for wcjE_11E

          #update locus ouput
          OutputLocus.df[col5_name] <- nonsence_mutation  #pseudo
          if(locus_result_type == "pseudo")
          {
            OutputLocusProfile.df[col5_name] <- nonsence_mutation
            SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[", nonsence_mutation, "]", sep="")
            cat(col5_name, nonsence_mutation, "\n", sep = "\t\t")
            SerotypeLookup.df <- filter(SerotypeLookup.df, (  (!!sym(col5_name)) == nonsence_mutation | is.na(!!sym(col5_name))  ))
          }

          #-------------------------------------------------------------------- GET MUTATIONS
          motifs <- ""
          if(Mutns == "Yes")
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

            if(NumMutations > 0 )
            {
              for(w in 1L:NumMutations)
              {
                aa_mut_name <- LocusMutationsCurr.df[w,"Name"]
                aa_start <- LocusMutationsCurr.df[w,"Posn_1"]
                aa_end <- LocusMutationsCurr.df[w,"Posn_2"]
                aa_mut <- LocusMutationsCurr.df[w,"Mutation"]
                aa_value<-substr(AASeqLine_str, aa_start, aa_end)
                if(motifs=="") 
                {
                  motif_sep <- ""
                }else
                {
                  motif_sep <- "/"
                }

                if(aa_mut == aa_value)
                {
                  motif<- paste(motif_sep, aa_start, aa_mut, sep = "")
                }else
                {
                  motif<-""
                }
                motifs <- paste(motifs, motif, sep = "")
              } #endfor
            }
          } #endif mutations found

          #update locus ouput
          OutputLocus.df[col4_name] <- motifs  #mutations
          if(locus_result_type == "mutations")
          {
            OutputLocusProfile.df[col4_name] <- motifs  #mutations
            SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[", motifs, "]", sep="")
            cat(col4_name, motifs, "\n", sep = "\t\t")
            SerotypeLookup.df <- filter(SerotypeLookup.df, ((!!sym(col4_name)) == motifs | is.na(!!sym(col4_name))))
          }

          #-------------------------------------------------------------------- GET ALLELES
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
              df.blastout2 <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
              names(df.blastout2) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
              df.blastout <- arrange(df.blastout2, -bit)
              outfile_allele <- paste(local_temp_dir, "blastout_lkup_allele.csv", sep = "")
              write.csv(df.blastout2, outfile_allele, quote = FALSE, row.names = F )
              Allele <- as.character(df.blastout2[1, "Allele"])
              AlleleParts <- unlist(strsplit(Allele, "_"))

              Allele2 <- AlleleParts[2]

              if(df.blastout2$Ident[1] >= 97)
              {
                Allele <- Allele2
              }else
              {
                Allele <- "NF" 
                outfile_nf <- paste(local_output_dir, "output_dna_notfound.fasta", sep = "")
                sink(outfile_nf, split=FALSE, append = TRUE)
                cat(">", locus, "_", CurrSampleNo, "_", "\n", DNASeqLine_NoDash_str, "\n", sep ="")
                sink()
              }
            }

            #update locus ouput
            OutputLocus.df[col3_name] <- Allele
            if(locus_result_type == "allele")
            {
              OutputLocusProfile.df[col3_name] <- Allele  #allele
              SampleProfile <- paste(SampleProfile, ProfileSeparator, locus2, "[", Allele, "]", sep="")
              cat(col3_name, Allele,"\n", sep = "\t\t")
              SerotypeLookup.df <- filter(SerotypeLookup.df, ( ( ((!!sym(col3_name)) == Allele) |  (is.na(!!sym(col3_name))) ) )    )
            }

          }#close bracket for DNA lookup file exists check

        }# end BLAST positive

      } #close bracket for wild type gene file found

    } #end of locus list loop 
    
    # Interpret molecular characterization to get serotype
    SizeSerotypeResult <- dim(SerotypeLookup.df)
    NumResults <- SizeSerotypeResult[1]
    if(NumResults == 0)
    {
      serotype <- paste("NT[", serotype, "]", sep = "")
    }else
    {
      for(w in 1L:NumResults)
      {
        if(w == 1)
        {
          serotype <-  SerotypeLookup.df$Serotype[w]
        }else
        {
          serotype <- paste(serotype, "/", SerotypeLookup.df$Serotype[w], sep = "")
        }
      }
    }
  }     #end not sample error

  if(m==1L)
  {
    SampleOutputProfile.df <- tibble(SampleNo = CurrSampleNo, Serotype = serotype, Profile = SampleProfile)
  }else
  {
    SampleOutputProfile2.df <- tibble(SampleNo = CurrSampleNo, Serotype = serotype, Profile = SampleProfile )
    SampleOutputProfile.df<- bind_rows(SampleOutputProfile.df, SampleOutputProfile2.df)
  }

  cat("\nProfile:  ", SampleProfile, sep = "")
  cat("\nSerotype: ", serotype, sep = "")

} #close sample list loop

write.csv(SampleOutputProfile.df, outfile, quote = FALSE,  row.names = F)
  
cat("\n\nDone!\n\n\n")
