##################################################################
####  WGS Analysis and Detection of Molecular Markers (WADE)  ####
####                 Author: Walter Demczuk                   ####
####                    Date: 2022-10-13                      ####
##################################################################

# Load libraries####
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyselect)
library(stringr)
library(Biostrings)

# USER INPUT -----------------------------------------------------------------------------------

    SampleNo <- "12345"           # Single sample no or "list" for list.csv
    Variable <- NA                # Information to be displayed in outputs (MIC value, metadata)
    LocusID <- "parC"             # Single locus gene or "list" for loci.csv
    curr_work_dir <- "C:\\MasterBlastR\\"
    ContigsDir <- "C:\\MasterBlastR\\contigs\\"

# FILES ----------------------------------------------------------------------------------------
    SampleListName <- "list.csv"  # File name of sample num list if SampleNo="list", otherwise ignored
    LocusListName <- "loci.csv"   # File name of loci list if LociID="list", otherwise ignored
    LocusMutationsName <- "loci_mutations.csv"
# ----------------------------------------------------------------------------------------------

cat("\n\n", "SampleNo", "\tMasterBlastR\n", sep = "")

############################################################################### Organize directories for lookups + mapping tables
# Make sure the two slashes are on the end of the directory paths   
if(str_sub(ContigsDir, start = -1, end = -1) != "\\") 
  {ContigsDir <- paste(ContigsDir, "\\", sep = "")}
if(str_sub(curr_work_dir, start = -1, end = -1) != "\\")
  {curr_work_dir <- paste(curr_work_dir, "\\", sep = "")}

# make directory structure
local_output_dir <- paste(curr_work_dir, "output\\", sep = "")
local_temp_dir <- paste(curr_work_dir, "temp\\", sep = "")
LkupDir <- paste(curr_work_dir, "allele_lkup_dna\\", sep = "")
SampListFile <- paste(curr_work_dir, SampleListName, sep = "")
LocusListFile <- paste(curr_work_dir, LocusListName, sep = "")
LocusMutationsFile <- paste(curr_work_dir, LocusMutationsName, sep = "")
 
# Define evalue variables 
evalueList <- paste(curr_work_dir, "blast_evalues.csv", sep = "")
blast_evalues.df <-  read.csv(evalueList, header = TRUE, sep = ",", stringsAsFactors = FALSE)
blast_evalues.df <- as_tibble(blast_evalues.df)
Blast_evalue <- as.character(blast_evalues.df$contig[1])
evalue_allele <- as.character(blast_evalues.df$allele[1])
blast_wt_id_threshold <- as.numeric(blast_evalues.df$wt_id[1])

# remove old output files if present
unlink(paste(local_output_dir, "output_dna.fasta", sep = ""))
unlink(paste(local_output_dir, "output_dna_notfound.fasta", sep = ""))
unlink(paste(local_output_dir, "output_aa.fasta", sep = ""))
unlink(paste(local_output_dir, "output_profile.csv", sep = ""))
#-----------------------------------------------------------------------------------------------

############################################################################### Read List Files
if(SampleNo == "list")
{
  SampleList.df <- as_tibble(read.csv(SampListFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
}else
{
  SampleList.df <- tibble(SampleNo, Variable)
}
names(SampleList.df) <- c("SampleNo", "Variable")
Size.df <- dim(SampleList.df)
NumSamples <- Size.df[1]

if(LocusID == "list")
{
  LocusList.df <- as_tibble(read.csv(LocusListFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
}else
{
  LocusList.df <- tibble(LocusID)
}

names(LocusList.df) <- c("Locus_id")
SizeList <- dim(LocusList.df)
NumLoci <- SizeList[1]

if(file.exists(LocusMutationsFile))
{
  Mutns <- "Yes"
  LocusMutations.df <- as_tibble(read.csv(LocusMutationsFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
}else 
{
  Mutns <- "No"
}

#-----------------------------------------------------------------------------------------------
  
# If there is an allele lookup file (allele_lkup_dna folder) for this locus, index the BLAST lookup database
for(q in 1L:NumLoci)
{
  locus <- as.character(LocusList.df[q,1])
  LocusLkupDNA <- paste(curr_work_dir, "allele_lkup_dna\\",  locus, ".fasta", sep = "")
  if(file.exists(LocusLkupDNA))
  {
    BlastFormatCommand <- paste("makeblastdb -in ", LocusLkupDNA, " -dbtype nucl", sep = "")
    try(system(BlastFormatCommand))
  }
}

#-----------------------------------------------------------------------------------------------
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Loop through samples in list.csv
m <- 1L
for(m in 1L:NumSamples)  
{
  CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
  CurrSampleVar <- as.character(SampleList.df[m, "Variable"])
  CurrSample.df <- filter(SampleList.df, SampleNo == CurrSampleNo)
  SampleProfile <- ""
  cat(CurrSampleNo, "\t",  m, " of ", NumSamples, "\n")
    
  #============================================================================ Loop through loci in loci.csv for each sample
  p <- 1L
  for(p in 1L:NumLoci)     
  {
    locus <- as.character(LocusList.df[p,1])
    LocusLkupDNA <- paste(curr_work_dir, "allele_lkup_dna\\", locus, ".fasta", sep = "")
    if(file.exists(LocusLkupDNA))
    {
      LocusLkupDNApresent = TRUE
    }else
    {
      LocusLkupDNApresent = FALSE
    }

    AlleleLine <- ""
    AlleleInfo <- NA
    AlleleInfo[1] <- NA
    AlleleInfo[2] <- NA
    AlleleInfo[3] <- NA
    AlleleInfo[4] <- NA
    IdentLine <- ""
    motifs <- ""

    LocusFile <- paste(curr_work_dir, "wildgenes\\", locus, ".fasta", sep = "")
    if(!file.exists(LocusFile))
    {
      sample_error.df <- tibble(Loucs_ID = locus, Output = "File Not Found", stringsAsFactors = FALSE)
      return(sample_error.df)
    }

    #.......................................................................... BLAST and parse blastout.txt
    QueryFile <- paste(ContigsDir, CurrSampleNo, ".fasta", sep = "")
    DestFile <- paste(local_temp_dir, "queryfile.fasta", sep = "")
    if(!file.copy(QueryFile, DestFile, overwrite = T))
    {
      AlleleInfo[1] <- "Sample_Err"
      AlleleInfo[2] <- NA
      AlleleInfo[3] <- NA
      AlleleInfo[4] <- NA
      IdLine <- ""
      IDpercent2 <- ""
      IdLine_trimmed <- ""
    }else
    #.......................................................................... make blast database of contig file, then blast locus against it
    {
      FormatCommand <- paste("makeblastdb -in ", DestFile, " -dbtype nucl", sep = "")
      shell(FormatCommand, intern = TRUE)

      BlastCommand <- paste("blastn -query ", LocusFile, " -db ", DestFile, " -out ", local_temp_dir, 
                          "blastout.txt -evalue ", Blast_evalue, sep = "")
      shell(BlastCommand, intern = TRUE)

      Blastout <- paste(local_temp_dir, "blastout.txt", sep = "")
      con <- file(Blastout, open="r")
      linn <- readLines(con)
              close(con)

      #........................................................................ check if gene was found in BLAST
      BlastResult <- NA
      for(i in 1:length(linn))
      {
        if(str_detect(linn[i], "No hits found"))
        {
          BlastResult <- "NEG"
          AlleleInfo[1] <- "NEG"
          AlleleInfo[2] <- ""
          AlleleInfo[3] <- ""
          AlleleInfo[4] <- ""
          break()
        }else
        {
          BlastResult <- "POS"
          AlleleInfo[1] <- "POS"
          AlleleInfo[2] <- ""
          AlleleInfo[3] <- ""
          AlleleInfo[4] <- ""
        }
      }

      #define blast parsing variables
      DNASeqLine_str <- ""
      WTDNASeqLine_str <- ""
      TimesThrough <- 0
      WTlengthLine <- ""
      WTlength <- 0
      IdLine <- ""
      IDlength <- ""
      IDpercent <- ""
      IdNum <- ""
      IdLine_trimmed <- ""

      #........................................................................ Parse blast_out.txt
      if(BlastResult == "POS")
      {
        for(i in 1:length(linn))
        {
          if(str_detect(linn[i], "Score =")) #set counter to parse only first align block of blastout
          {
            TimesThrough <- TimesThrough + 1
          }

          if((str_detect(linn[i], "Length=" )) & (WTlength == 0)) #get length of WT gene
          {
            WTlengthLine <- unlist(linn[i])
            WTlength <- as.numeric(substr(WTlengthLine, 8, 50))
          }

          if(TimesThrough == 1)
          {
            if(str_detect(linn[i], "Identities"))
            {
              IdLine <- unlist(linn[i])
              IdLine <- substr(IdLine, 15, 50)
              IdLineParts <- strsplit(IdLine, "/")
              IDLineParts2 <- unlist(IdLineParts)
              IDlength <- as.numeric(IDLineParts2[1])
              IDcoverage <- ((IDlength / WTlength) * 100)
              IDpercent <- sprintf("%3.1f%%", IDcoverage)
              IDLineParts3 <- strsplit(IdLine, ",")
              IDLineParts3_2 <- unlist(IDLineParts3)
              IdLine_trimmed <- IDLineParts3_2[1]
            }
            
            if(str_detect(linn[i], "Query "))
            {
              QueryLine <- unlist(strsplit(linn[i], " "))
              QueryLine <- QueryLine[QueryLine != ""]
              WTDNASeqLine_str <- paste(WTDNASeqLine_str, QueryLine[3], sep = "")
            }

            if(str_detect(linn[i], "Sbjct "))
            {
              SbjctLine <- unlist(strsplit(linn[i], " "))
              SbjctLine <- SbjctLine[SbjctLine != ""]
              DNASeqLine_str <- paste(DNASeqLine_str, SbjctLine[3], sep = "")
            }
          }
        }

        WTDNASeqLine <- DNAString(WTDNASeqLine_str)
        WTDNASeqLine_NoDash_str <- str_replace_all(WTDNASeqLine_str, "-", "")
        WTDNASeqLine_NoDash <- DNAString(WTDNASeqLine_NoDash_str)
                
        if(LocusID != "list")
        {
          cat("\n\n>", locus , "(Wildtype)\n", WTDNASeqLine_NoDash_str, "\n", sep ="")
        }

        DNASeqLine <- DNAString(DNASeqLine_str)
        DNASeqLine_NoDash_str1 <- str_replace_all(DNASeqLine_str, "-", "")
        DNASeqLine_NoDash_str <- str_replace_all(DNASeqLine_NoDash_str1, "N", "")
        DNASeqLine_NoDash <- DNAString(DNASeqLine_NoDash_str)
                
        if(LocusID != "list")
        {
          cat("\n\n>", CurrSampleNo, "_", locus , "\n", DNASeqLine_NoDash_str, "\n", sep ="")
        }
                
        #...................................................................... MAKE PROTEIN SEQUENCE
        WTAASeqLine <- translate(WTDNASeqLine_NoDash)
        WTAASeqLine_str <- toString(WTAASeqLine)
                
        if(LocusID != "list")
        {
          cat("\n\n>", locus , "(Wildtype)\n", WTAASeqLine_str, "\n", sep ="")
        }
                
        AASeqLine <- translate(DNASeqLine_NoDash)
        AASeqLine_str <- toString(AASeqLine)

        if(LocusID != "list")
        {
          cat("\n\n>", CurrSampleNo, "_", locus , "\n", AASeqLine_str, "\n", sep ="")
        }

        #-----------make sure AASeqLine, WTAASeqLine are listed in correct order!
        globalAlign_AA <- pairwiseAlignment(AASeqLine, WTAASeqLine, substitutionMatrix = "BLOSUM50", 
                                            gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
        WTAASeqLine_aln <- subject(globalAlign_AA)
        WTAASeqLine_aln_str <- toString(WTAASeqLine_aln)
        AASeqLine_aln <- pattern(globalAlign_AA)
        AASeqLine_aln_str <- toString(AASeqLine_aln)
                
        DNASeqLine_aln <- ""
        for(j in 1:str_length(WTDNASeqLine_str))
        {
          if(str_sub(WTDNASeqLine_str, j, j) == str_sub(DNASeqLine_str, j, j))
          {
            DNASeqLine_aln <- paste(DNASeqLine_aln, ".", sep = "")
          }else
          {
            DNASeqLine_aln <- paste(DNASeqLine_aln, str_sub(DNASeqLine_str, j, j), sep = "")
          }
        }

        if(LocusID != "list")
        {
          cat("\n\nDNA Alignment:\n", DNASeqLine_aln, "\n", sep = "")
        }

        #---------------------------------------------------------------------- MUTATIONS AND MOTIFS
        motifs <- ""
        aa_mut_name <- ""
        aa_start <- ""
        aa_end <- ""
        aa_wt <- ""
        motif <- ""

        if(Mutns == "Yes")
        {
          LocusMutationsCurr.df <- filter(LocusMutations.df, Locus_id == locus)
          SizeMutns.df <- dim(LocusMutationsCurr.df)
          NumMutations <- SizeMutns.df[1]

          if(NumMutations > 0)
          {
            for(w in 1L:NumMutations)
            {
              aa_mut_name <- LocusMutationsCurr.df[w,"Name"]
              aa_start <- LocusMutationsCurr.df[w,"Posn_1"]
              aa_end <- LocusMutationsCurr.df[w,"Posn_2"]
              aa_wt <- LocusMutationsCurr.df[w,"WildType"]
              motif<-substr(AASeqLine_aln_str, aa_start, aa_end)
                             
              if(motif == aa_wt)
              {
                motif <- "WT"
                aa_mut_name <- ""
              }
              
              if(w == 1)
              {
                motifs <- paste(aa_mut_name, motif, sep = "")
              }else
              {
                motifs <- paste(motifs, "/", aa_mut_name, motif, sep = "")
              }
            }
          }
        }

        AASeqLine_aln_disp <- ""
        mutations <- ""

        for(k in 1:str_length(WTAASeqLine_aln_str))
        {
          if(str_sub(WTAASeqLine_aln_str, k, k) == str_sub(AASeqLine_aln_str, k, k))
          {
            AASeqLine_aln_disp <- paste(AASeqLine_aln_disp, ".", sep = "")
          }else
          {
            AASeqLine_aln_disp <- paste(AASeqLine_aln_disp, str_sub(AASeqLine_str, k, k), sep = "")
            mutations <- paste(mutations, str_sub(WTAASeqLine_aln_str, k, k), 
                               k, str_sub(AASeqLine_str, k, k), " ", sep = "")
          }
        }
        
        if(LocusID != "list")
        {
          cat("\n\nProtein Alignment:\n", AASeqLine_aln_disp, "\n", sep = "")
          cat("\nMutations: ", mutations, "\n\n", sep = "")
          cat("\nMotifs: ", motifs, "\n\n", sep = "")
        }
          
        #---------------------------------------------------------------------- ALLELE LOOKUP in allele_lkup_dna folder
          
        # write DNASeqLine_NoDash_str to a file named = querygene.fasta,
        # BLAST vs. lookup table,
        # parse out the allele numbers

        Seq_File <- paste(local_temp_dir, "querygene.fasta", sep="")
        sink(Seq_File, split=FALSE, append = FALSE)
        cat(">", CurrSampleNo, "_", locus , "\n", DNASeqLine_NoDash_str, sep ="")
        sink()
        ExactMatchFound <- FALSE

        if(LocusLkupDNApresent)
        {
          #BLAST lookup table
          BlastCommand2 <- paste("blastn -query ", Seq_File, " -db ", LocusLkupDNA, " -out ", local_temp_dir, "blastout2.txt -evalue ", evalue_allele , " -num_alignments 1", sep = "")
          shell(BlastCommand2, intern = TRUE)
          Blastout2 <- paste(local_temp_dir, "blastout2.txt", sep = "")
          con <- file(Blastout2, open="r")
          linn <- readLines(con)
          close(con)

          #parse blastout2
          for (i in 1:length(linn))
          {
            if(str_detect(linn[i], "Identities"))
            {
              IdentLine <-  unlist(linn[i])
              IdentLine <- substr(IdentLine, 15, 50)
              if(str_detect(IdentLine, "100%"))
              {
                ExactMatchFound <- TRUE
              }
            }

            if (str_detect(linn[i], ">"))
            {
              AlleleLine <- unlist(linn[i])
              AlleleLine <- substr(AlleleLine, 2, 50)
              AlleleParts <- strsplit(AlleleLine, "_")
              AlleleParts2 <- unlist(AlleleParts)
            }
          }

          if(ExactMatchFound)   #if only not found in fasta make files here.
          {
            #AlleleParts2[1] holds POS/NEG; [2] allele number; [3] mutations or WT; [4] extra info
            AlleleInfo[2] <- AlleleParts2[2]  #allele number
            AlleleInfo[3] <- AlleleParts2[3]  #mutations
            AlleleInfo[4] <- AlleleParts2[4]  #comments
          }else
          {
            AlleleInfo[2] <- "NF"
            AlleleInfo[3] <- "???"
            AlleleInfo[4] <- ""
          }
        }#close bracket for DNA lookup file exists check

        #---------------------------------------------------------------------- WRITE OUTPUT fasta files output_dna.fasta and output_aa.fasta
        dna_file <- paste(local_output_dir, "output_dna.fasta", sep = "")
        dna_nf_file <- paste(local_output_dir, "output_dna_notfound.fasta", sep = "")
        aa_file <- paste(local_output_dir, "output_aa.fasta", sep = "")

        sink(dna_file, split=FALSE, append = TRUE)
        cat(">", locus, "_", CurrSampleNo, "_", locus, AlleleInfo[2], "_", AlleleInfo[3], "_", CurrSampleVar, "_", motifs, "\n", DNASeqLine_NoDash_str, "\n", sep ="")
        sink()

        if(AlleleInfo[2] == "NF")
        {
          sink(dna_nf_file, split=FALSE, append = TRUE)
          cat(">", locus, "_", motifs, "_", CurrSampleNo, "_", CurrSampleVar, "_", "\n", DNASeqLine_NoDash_str, "\n", sep ="")
          sink()
        }

        sink(aa_file, split=FALSE, append = TRUE)
        cat(">", locus, "_", CurrSampleNo, "_", CurrSampleVar, "\n", AASeqLine_str, "\n", sep ="")
        sink()

        IDpercent2 <- paste(IDlength, "/", WTlength, " (", IDpercent, ")", sep = "") #made a blast ID line based on full WT gene

        if(IDcoverage <= blast_wt_id_threshold)
        {
          AlleleInfo[1] <- "NEG"  #blast result
        }
      }else  #close bracket for BLAST positive
      {
        IdLine <-""
        IDpercent2 <- ""
        IdLine_trimmed <-""
      }
    } #close bracket for contig file found.

    col1_name <- paste(locus, "_result", sep = "")
    col2_name <- paste(locus, "_allele", sep = "")
    col3_name <- paste(locus, "_mutations", sep = "")
    col4_name <- paste(locus, "_comments", sep = "")
    col5_name <- paste(locus, "_BlastID", sep = "")
    col6_name <- paste(locus, "_PctWT", sep = "")
    col7_name <- paste(locus, "_motifs", sep = "")

    headers <- c(col1_name, col2_name, col3_name, col4_name, col5_name, col6_name, col7_name)

    if(p==1)
    {
      OutputLocus.df <- tibble(AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], IdLine_trimmed, IDpercent2, motifs)
        names(OutputLocus.df) <- headers
    }else
    {
      OutputLocus2.df <- tibble(AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], IdLine_trimmed, IDpercent2, motifs)
      names(OutputLocus2.df) <- headers
      OutputLocus.df <- cbind(OutputLocus.df, OutputLocus2.df)
    }

    ProfileEntry <- ""
    if(AlleleInfo[1] == "POS")
    {
      ProfileEntry<-locus
      if(LocusLkupDNApresent)
      {
        if(AlleleInfo[3] == "" | is.na(AlleleInfo[3]) | AlleleInfo[3] == "???")
        {
          ProfileEntry<- paste(locus, AlleleInfo[3], sep = " ")
        }else
        if(AlleleInfo[3]=="WT" | AlleleInfo[3]=="WT/WT" | AlleleInfo[3]=="WT/WT/WT" | AlleleInfo[3]=="WT/WT/WT/WT")
        {
          ProfileEntry<-""
        }else
        {
          ProfileEntry <- paste(locus, AlleleInfo[3], sep = " ")
        }
      }

      if(SampleProfile == "")
      {
        SampleProfile <- ProfileEntry
      }else
      {
        if(ProfileEntry!=""){ProfileEntry<-paste("-", ProfileEntry, sep = "")}
        SampleProfile <- paste(SampleProfile, ProfileEntry, sep = "")
      }
    }

    if(AlleleInfo[1] == "Sample Err")
    {
      SampleProfile <- "Sample Err"
    }

    cat(CurrSampleNo, locus, AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], 
        IdLine_trimmed, IDpercent2, motifs, "\n", sep = "\t")
    
  } #end of locus list loop ================================================================================

  SampleProfile.df <- tibble(SampleProfile)
  cat("Molecular Profile:  ", SampleProfile, "\n\n", sep = "")
  OutputLocus1.df <- cbind(CurrSample.df, OutputLocus.df, SampleProfile.df)

  if(m==1)  #if first sample make one row profile table, otherwise add new row to table
  {
    OutputProfile.df <- tibble(OutputLocus1.df)
    headers2 <- names(OutputLocus1.df)
    names(OutputProfile.df) <- headers2
  }else
  {
    OutputProfile.df <- rbind(OutputProfile.df, OutputLocus1.df)
  }
} #close bracket for sample list loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

outfile <- paste(local_output_dir, "output_profile.csv", sep = "")

write.csv(OutputProfile.df, outfile, quote = FALSE,  row.names = F)

cat("\n\nDone!\n\n\n")
