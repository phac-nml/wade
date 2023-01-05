##################################################################
####                          WAMR                            ####
####  WGS Analysis and Detection of Molecular Markers (WADE)  ####
####                 Author: Walter Demczuk                   ####
####                    Date: 2022-11-18                      ####
##################################################################

# Load libraries####
library(plyr)
library(dplyr)
library(tidyverse)
library(tidyselect)
library(stringr)
library(Biostrings)

#  USER INPUT ------------------------------------------------------------------------------

SampleNo <- "SC21-0169-P"     # Single sample no or "list" for list.csv
Variable <- NA         # Information to be displayed in outputs (MIC value, metadata)
LocusID <- "list"      # Single locus gene or "list" for loci.csv
curr_work_dir <- "C:\\WamR-Pneumo\\"
ContigsDir <- "C:\\WamR-Pneumo\\contigs\\
vcf_folder <- "C:\\WamR-Pneumo\\vcf\\"
    
#-------------------------------------------------------------------------------------------
Org_id <- "PNEUMO"
    
# get directory structure + set location of directory for lookups, mapping tables
local_output_dir <- paste(curr_work_dir, "output\\", sep = "")
local_temp_dir <- paste(curr_work_dir, "temp\\", sep = "")
ref_dir <- paste(curr_work_dir, "reference\\", sep = "")
LkupDir <- paste(curr_work_dir, "allele_lkup_dna\\", sep = "")
SampListFile <- paste(curr_work_dir, "list.csv", sep = "")
LocusListFile <- paste(ref_dir, "loci.csv", sep = "")
LocusMutationsFile <- paste(ref_dir, "loci_mutations.csv", sep = "")
evalueList <- paste(ref_dir, "blast_evalues.csv", sep = "")
blast_evalues.df <-  read.csv(evalueList, header = TRUE, sep = ",", stringsAsFactors = FALSE)
blast_evalues.df <- as_tibble(blast_evalues.df)
Blast_evalue <- as.character( blast_evalues.df$contig[1])
evalue_allele <- as.character( blast_evalues.df$allele[1])
blast_wt_id_threshold <- as.numeric( blast_evalues.df$wt_id[1])

# set output files - need to delete previous outputs or errors when left open by user.
unlink(paste(local_output_dir, "output_dna.fasta", sep = "")) #this deletes the old file!
unlink(paste(local_output_dir, "output_dna_notfound.fasta", sep = ""))
unlink(paste(local_output_dir, "output_aa.fasta", sep = ""))
unlink(paste(local_output_dir, "output_profile_AMR_markers.csv", sep = ""))
unlink(paste(local_output_dir, "output_profile_23S.csv", sep = ""))
unlink(paste(local_output_dir, "output_profile_interpret_AMR.csv", sep = ""))

#-------------------------------------------------------------------------------------------

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

#-------------------------------------------------------------------------------------------
if(SampleNo == "list")
{SampleList.df <- as_tibble(read.csv(SampListFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
}else
{SampleList.df <- tibble(SampleNo, Variable)}

names(SampleList.df) <- c("SampleNo", "Variable")
Size.df <- dim(SampleList.df)
NumSamples <- Size.df[1]
  
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Part 1) MasterBlastR for the list of AMR determinants in loci.csv
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

cat("\n\n", "SampleNo", "\t Part 1: AMR Determinants.\n", sep = "")

if(LocusID == "list")
{LocusList.df <- as_tibble(read.csv(LocusListFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
}else
{LocusList.df <- tibble(LocusID)}

names(LocusList.df) <- c("Locus_id")
SizeList <- dim(LocusList.df)
NumLoci <- SizeList[1]

if(file.exists(LocusMutationsFile))
{
  Mutns <- "Yes"
  LocusMutations.df <- as_tibble(read.csv(LocusMutationsFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
}else 
{Mutns <- "No"}

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  Loop through samples in list.csv
m <- 1L
for(m in 1L:NumSamples)  
{
  CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
  CurrSampleVar <-as.character(SampleList.df[m, "Variable"])
  CurrSample.df <- filter(SampleList.df, SampleNo == CurrSampleNo)
  SampleProfile <- ""
  cat(m, " of ", NumSamples, "\n")
  
  #===========================================================================  Loop through loci in loci.csv for each sample
  p <- 1L
  for(p in 1L:NumLoci)     
  {
    locus <- as.character(LocusList.df[p,1])
    LocusLkupDNA <- paste(curr_work_dir, "allele_lkup_dna\\", locus, ".fasta", sep = "")
    if(file.exists(LocusLkupDNA))
    {LocusLkupDNApresent = TRUE
    }else
    {LocusLkupDNApresent = FALSE}

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
    { #make blast database of contig file, then blast locus against it
      FormatCommand <- paste("makeblastdb -in ", DestFile, " -dbtype nucl", sep = "")
      shell(FormatCommand, intern = TRUE)
      BlastCommand <- paste("blastn -query ", LocusFile, " -db ", DestFile, " -out ", local_temp_dir, "blastout.txt -evalue ", Blast_evalue, sep = "")
      shell(BlastCommand, intern = TRUE)
      Blastout <- paste(local_temp_dir, "blastout.txt", sep = "")
      con <- file(Blastout, open = "r")
      linn <- readLines(con)
      close(con)

      #check if gene was found in BLAST
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

      # Parse blast_out.txt ----------------------------------------------------
        if(BlastResult == "POS")
        {
          for(i in 1:length(linn))
          {
            if(str_detect(linn[i], "Score =")) #set counter to parse only first align block of blastout
            {TimesThrough <- TimesThrough + 1}
            
            if((str_detect(linn[i], "Length =" )) & (WTlength == 0)) #get length of WT gene
            {
              WTlengthLine <-  unlist(linn[i])
              WTlength <- as.numeric(substr(WTlengthLine, 8, 50))
            }

            if(TimesThrough == 1)
            {
              if(str_detect(linn[i], "Identities"))
              {
                IdLine <-  unlist(linn[i])
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
          
          if(LocusID != "list")
          {cat("\n\n>", locus , "(Wildtype)\n", WTDNASeqLine_NoDash_str, "\n", sep ="")}

          DNASeqLine <- DNAString(DNASeqLine_str)
          DNASeqLine_NoDash_str1 <- str_replace_all(DNASeqLine_str, "-", "")
          DNASeqLine_NoDash_str <- str_replace_all(DNASeqLine_NoDash_str1, "N", "")
          DNASeqLine_NoDash <- DNAString(DNASeqLine_NoDash_str)
          
          if(LocusID != "list")
          {cat("\n\n>", CurrSampleNo, "_", locus , "\n", DNASeqLine_NoDash_str, "\n", sep = "")}
          
          #-------------------------------------------------------------------- make protein sequence
          WTAASeqLine <- translate(WTDNASeqLine_NoDash)
          WTAASeqLine_str <- toString(WTAASeqLine)
          
          if(LocusID != "list")
          {cat("\n\n>", locus , "(Wildtype)\n", WTAASeqLine_str, "\n", sep ="")}
          
          AASeqLine <- translate(DNASeqLine_NoDash)
          AASeqLine_str <- toString(AASeqLine)
          
          if(LocusID != "list")
          {cat("\n\n>", CurrSampleNo, "_", locus , "\n", AASeqLine_str, "\n", sep = "")}

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
            {DNASeqLine_aln <- paste(DNASeqLine_aln, ".", sep = "")
            }else
            {DNASeqLine_aln <- paste(DNASeqLine_aln, str_sub(DNASeqLine_str, j, j), sep = "")}
          }

          if(LocusID != "list")
          {cat("\n\nDNA Alignment:\n", DNASeqLine_aln, "\n", sep = "")}

          #-------------------------------------------------------------------- MUTATIONS AND MOTIFS
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

            if(NumMutations > 0 )
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
                {motifs <- paste(aa_mut_name, motif, sep = "")
                }else
                {motifs <- paste(motifs, "/", aa_mut_name, motif, sep = "")}
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
              mutations <- paste(mutations, str_sub(WTAASeqLine_aln_str, k, k), k, str_sub(AASeqLine_str, k, k), " ", sep = "")
            }
          }
          
          if(LocusID != "list")
          {
            cat("\n\nProtein Alignment:\n", AASeqLine_aln_disp, "\n", sep = "")
            cat("\nMutations: ", mutations, "\n\n", sep = "")
            cat("\nMotifs: ", motifs, "\n\n", sep = "")
          }
          
          #-------------------------------------------------------------------- ALLELE LOOKUP in allele_lkup_dna folder
          
          # write DNASeqLine_NoDash_str to a file named = querygene.fasta,
          # BLAST against lookup table,
          # parse out the allele numbers

          Seq_File <- paste(local_temp_dir, "querygene.fasta", sep = "")
          sink(Seq_File, split = FALSE, append = FALSE)
          cat(">", CurrSampleNo, "_", locus , "\n", DNASeqLine_NoDash_str, sep = "")
          sink()
          
          ExactMatchFound <- FALSE

          if(LocusLkupDNApresent)
          {
            BlastCommand2 <- paste("blastn -query ", Seq_File, " -db ", LocusLkupDNA, " -out ", local_temp_dir, "blastout2.txt -evalue ", evalue_allele , " -num_alignments 1", sep = "")
            shell(BlastCommand2, intern = TRUE)
            Blastout2 <- paste(local_temp_dir, "blastout2.txt", sep = "")
            con <- file(Blastout2, open = "r")
            linn <- readLines(con)
            close(con)

            for(i in 1:length(linn))
            {
              if(str_detect(linn[i], "Identities"))
              {
                IdentLine <-  unlist(linn[i])
                IdentLine <- substr(IdentLine, 15, 50)
                
                if (str_detect(IdentLine, "100%"))
                {ExactMatchFound <- TRUE}
              }

              if(str_detect(linn[i], ">"))
              {
                AlleleLine <- unlist(linn[i])
                AlleleLine <- substr(AlleleLine, 2, 50)
                AlleleParts <- strsplit(AlleleLine, "_")
                AlleleParts2 <- unlist(AlleleParts)
              }
            }

            if(ExactMatchFound)   #if only not found in fasta make files here.
            #AlleleParts2[1] holds POS/NEG; [2] allele number; [3] mutations or WT; [4] extra info
            {
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

          #-------------------------------------------------------------------- WRITE OUTPUT fasta files output_dna.fasta and output_aa.fasta
          dna_file <- paste(local_output_dir, "output_dna.fasta", sep = "")
          dna_nf_file <- paste(local_output_dir, "output_dna_notfound.fasta", sep = "")
          aa_file <- paste(local_output_dir, "output_aa.fasta", sep = "")
          sink(dna_file, split = FALSE, append = TRUE)
          cat(">", locus, "_", CurrSampleNo, "_", locus, AlleleInfo[2], "_", AlleleInfo[3], "_", CurrSampleVar, "_", motifs, "\n", DNASeqLine_NoDash_str, "\n", sep ="")
          sink()

          if(AlleleInfo[2] == "NF")
          {
            sink(dna_nf_file, split = FALSE, append = TRUE)
            cat(">", locus, "_", motifs, "_", CurrSampleNo, "_", CurrSampleVar, "_", "\n", DNASeqLine_NoDash_str, "\n", sep = "")
            sink()
          }

          sink(aa_file, split = FALSE, append = TRUE)
          cat(">", locus, "_", CurrSampleNo, "_", CurrSampleVar, "\n", AASeqLine_str, "\n", sep = "")
          sink()
          
          IDpercent2 <- paste(IDlength, "/", WTlength, " (", IDpercent, ")", sep = "") #made a blast ID line based on full WT gene

          if(IDcoverage <= blast_wt_id_threshold)
          {AlleleInfo[1] <- "NEG"}  #blast result

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

      if(p == 1)
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
          {ProfileEntry<- paste(locus, AlleleInfo[3], sep = " ")
          }else if(AlleleInfo[3] == "WT" | AlleleInfo[3] == "WT/WT" | AlleleInfo[3] == "WT/WT/WT" | AlleleInfo[3] == "WT/WT/WT/WT")
          {ProfileEntry<-""
          }else
          {ProfileEntry <- paste(locus, AlleleInfo[3], sep = " ")}
        }

        if(SampleProfile == "")
        {
          SampleProfile <- ProfileEntry
        }else
        {
          if(ProfileEntry != "")
          {ProfileEntry<-paste("-", ProfileEntry, sep = "")}
          SampleProfile <- paste(SampleProfile, ProfileEntry, sep = "")
        }
      }

      if(AlleleInfo[1] == "Sample Err")
      {SampleProfile <- "Sample Err"}

      cat(CurrSampleNo, locus, AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], IdLine_trimmed, IDpercent2, motifs, "\n", sep = "\t")

  } #end of locus list loop ================================================================

  SampleProfile.df <- tibble(SampleProfile)
  cat("Molecular Profile:  ", SampleProfile, "\n\n", sep = "")
  OutputLocus1.df <- cbind(CurrSample.df, OutputLocus.df, SampleProfile.df)

    if(m == 1)  #if first sample make one row profile table, otherwise add new row to table
    {
      OutputProfile.df <- tibble(OutputLocus1.df)
      headers2 <- names(OutputLocus1.df)
      names(OutputProfile.df) <- headers2
    }else
    {
      OutputProfile.df <- rbind(OutputProfile.df, OutputLocus1.df)
    }
} #close bracket for sample list loop <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

outfile <- paste(local_output_dir, "output_profile_AMR_markers.csv", sep = "")
write.csv(OutputProfile.df, outfile, quote = FALSE,  row.names = F)
cat("\n\nAMR determinants done!\n\n\n")

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Part 2) 23S rRNA analysis
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# These 23S rRNA position texts are put into the VCF files from the fasta header in
#    the 23S rRNA reference sequences (in the wildgenes folder).  Use these files
#    as a mapping reference and then download the VCF files.

if(LocusID == "list")
{
  cat("\n\n", "SampleNo", "\t Part 2: 23S rRNA mutations.\n", sep = "")
  
  if (Org_id == "PNEUMO")
  {
    rRNA23S_position1 <- "23S_rRNA_R6_sprr02	2061"
    rRNA23S_position2 <- "23S_rRNA_R6_sprr02	2613"
  }
  
  OutputProfile.df <- tibble(SampleNo = character(), A2059G = integer(), C2611T = integer())

  m <- 1
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
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
    
    if(QueryFile != "File not found")
    {
      allele_fraction <- 0L
      con <- file(QueryFile, open = "r")
      linn <- readLines(con)
      close(con)
      
      for(i in 1:length(linn))
      {
        if(str_detect(linn[i], rRNA23S_position1))  #2061 for pneumo
        {
          vcf_parts <- strsplit(linn[i], ";")
          vcf_parts <- unlist(vcf_parts)
          DP1 <- vcf_parts[8]
          DP<-as.integer(substr(DP1, 4, 10))
          AO1 <- vcf_parts[6]
          AO<-as.integer(substr(AO1, 4, 10))
          
          allele_fraction <- as.integer((AO/DP)*100)
          if(allele_fraction < 14) 
          {A2059G <- 0L}
          if((allele_fraction >= 15) && (allele_fraction <= 34)) 
          {A2059G <- 1L}
          if((allele_fraction >= 35) && (allele_fraction <= 64)) 
          {A2059G <- 2L}
          if((allele_fraction >= 65) && (allele_fraction <= 84)) 
          {A2059G <- 3L}
          if(allele_fraction >= 85) 
          {A2059G <- 4L}
        }else
        {A2059G <- 0L}
        
        if(str_detect(linn[i], rRNA23S_position2))  #2061 for pneumo
        {
          vcf_parts2 <- strsplit(linn[i], ";")
          vcf_parts2 <- unlist(vcf_parts2)
          DP2_str <- vcf_parts2[8]
          DP2 <- as.integer(substr(DP2_str, 4, 10))  #total number of reads (depth of reads)
          AO2_str <- vcf_parts2[6]
          AO2 <- as.integer(substr(AO2_str, 4, 10))  #alternate observations from reference
          
          allele_fraction <- as.integer((AO2/DP2)*100)
          if(allele_fraction < 14)
          {C2611T <- 0L}
          if((allele_fraction >= 15) && (allele_fraction <= 34)) 
          {C2611T <- 1L}
          if((allele_fraction >= 35) && (allele_fraction <= 64)) 
          {C2611T <- 2L}
          if((allele_fraction >= 65) && (allele_fraction <= 84)) 
          {C2611T <- 3L}
          if(allele_fraction >= 85) 
          {C2611T <- 4L}
          
        }else
        {C2611T <- 0L}
      }
    }
    
    cat(CurrSampleNo, "\tAC_2059: ", A2059G, "\tAC_2611: ", C2611T,"\n", sep = "")
    SampleProfile.df <- tibble(CurrSampleNo, A2059G, C2611T)
    OutputProfile.df <- rbind(OutputProfile.df, SampleProfile.df)
  }
  
  names(OutputProfile.df) <- c("SampleNo", "A2059G", "C2611T")
  outfile <- paste(local_output_dir, "output_profile_23S.csv", sep = "")
  write.csv(OutputProfile.df, outfile, row.names = F)
  cat("\n\n23S rRNA allele analysis done!\n\n\n")
  
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Part 3) InterpretAMR combines previous 2 steps; calculates MICs, AMR profiles
  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  
  cat("\n\n Part 3: MIC Interpretations.\n", sep = "")
  
  AMR_Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_AMR_markers.csv", sep = ""),
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Size.df <- dim(AMR_Output.df)
  NumSamples <- Size.df[1]
  NumLoci <- ((Size.df[2]-3) / 7)
  
  rRNA23S_Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_23S.csv", sep = ""),
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Combined_Output.df <- full_join(AMR_Output.df, rRNA23S_Output.df, by = "SampleNo")
  
  Size.df <- dim(Combined_Output.df)
  NumSamples <- Size.df[1]
  
  pen_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_penicillin.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  cro_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_ceftriaxone.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  cfm_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_cefuroxime.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  ery_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_erythromycin.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  azi_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_azithromycin.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  cla_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_clarithromycin.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  cli_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_clindamycin.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  lev_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_levofloxacin.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  mox_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_moxifloxacin.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  tet_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_tetracycline.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  dox_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_doxycycline.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  sxt_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_trimeth_sulfa.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  chl_mic.df <- as_tibble(read.csv(paste(curr_work_dir, "reference\\inc_mic_chloramphenicol.csv", sep = ""),
                                   header = TRUE, sep = ",", stringsAsFactors = FALSE))
  sepr <- " - "
  
  m <- 1
  ErrorFound <- FALSE
  for(m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Build Molecular Profile
  {
    molec_profile <- NA
    lw_comments <- ""
    lw_CurrSampleNo <- as.character(Combined_Output.df[m, "SampleNo"])
    lw_ermB <- as.character(Combined_Output.df[m, "ermB_result"])
    
    if(lw_ermB != "Sample_Err")
    {
      lw_ermB_allele <- as.character(Combined_Output.df[m, "ermB_allele"])  #check for NF(Not Found found in allele lkup)
      if(lw_ermB_allele == "NF") #if ermB-Resistant then no allele req'd, if found then get allele
      {lw_ermB_allele <- "+"}
      if(lw_ermB == "NEG")
      {lw_ermB_allele <- "-"}
      
      lw_pbp1a_allele <- as.character(Combined_Output.df[m, "pbp1a_allele"])
      lw_pbp2b_allele <- as.character(Combined_Output.df[m, "pbp2b_allele"])
      lw_pbp2x_allele <- as.character(Combined_Output.df[m, "pbp2x_allele"])
      lw_folA_allele <- as.character(Combined_Output.df[m, "folA_allele"])
      lw_folP_allele <- as.character(Combined_Output.df[m, "folP_allele"])
      lw_gyrA_allele <- as.character(Combined_Output.df[m, "gyrA_allele"])
      lw_parC_allele <- as.character(Combined_Output.df[m, "parC_allele"])
      lw_allele_profile <- paste("pbp1a", lw_pbp1a_allele, ":",
                                 "pbp2b", lw_pbp2b_allele, ":",
                                 "pbp2x", lw_pbp2x_allele, ":",
                                 "ermB", lw_ermB_allele, ":",
                                 "folA", lw_folA_allele, ":",
                                 "folP", lw_folP_allele, ":",
                                 "gyrA", lw_gyrA_allele, ":",
                                 "parC", lw_parC_allele)
      
      #-------------------------------------------------------------------------------------
      lw_ermB_mut <- as.character(Combined_Output.df[m, "ermB_mutations"])  #check for ermB(S) -- Susceptible
      
      if(lw_ermB == "POS")
      {
        molec_profile <- "ermB"
        if(lw_ermB_mut == "S")
        {
          lw_ermB <- "NEG"
          molec_profile <- "ermB(S)"
        }
      }else
      {
        molec_profile <- NA
      }
      
      if(lw_ermB_mut == "S")
      {lw_ermB <- "NEG"}
      
      lw_ermTR <- as.character(Combined_Output.df[m, "ermTR_result"])
      if(lw_ermTR == "POS")
      {
        if(is.na(molec_profile))
        {molec_profile <- "ermTR"
        }else
        {molec_profile <- paste(molec_profile, sepr, "ermTR", sep = "")}
      }
      
      lw_mefAE <- as.character(Combined_Output.df[m, "mefAE_result"])
      if(lw_mefAE == "POS")
      {
        if(is.na(molec_profile))
        {molec_profile <- "mefAE"
        }else
        {molec_profile <- paste(molec_profile, sepr, "mefAE", sep = "")}
      }
      
      lw_23S_prof <- NA
      
      if(is.na(Combined_Output.df[m, "A2059G"]))
      {
        lw_23S_A2059G <- 0L
        lw_23S_prof <- "23S Err"
        lw_comments <- paste(lw_comments, "23S VCF file not found.", sep = "")
      }else
      {
        lw_23S_A2059G <- as.integer(Combined_Output.df[m, "A2059G"])
      }
      
      if(is.na(Combined_Output.df[m, "C2611T"]))
      {
        lw_23S_C2611T <- 0L
        lw_23S_prof <- "23S Err"
      }else
      {
        lw_23S_C2611T <- as.integer(Combined_Output.df[m, "C2611T"])
      }
      
      if(lw_23S_A2059G != 0L)
      {
        if(is.na(lw_23S_prof))
        {lw_23S_prof <- paste("23S rRNA A2059G ", lw_23S_A2059G, "/4", sep = "")
        }else
        {lw_23S_prof <- paste(ls_23S_prof, "/A2059G ", lw_23S_A2059G, "/4",  sep = "")}
      }
      
      if(lw_23S_C2611T != 0L)
      {
        if(is.na(lw_23S_prof))
        {lw_23S_prof <- paste("23S rRNA C2611T ", lw_23S_C2611T, "/4", sep = "")
        }else
        {lw_23S_prof <- paste(lw_23S_prof, "/C2611T ", lw_23S_C2611T, "/4",  sep = "")}
      }
      
      if(!is.na(lw_23S_prof))
      {
        if(is.na(molec_profile))
        {molec_profile <- lw_23S_prof
        }else
        {molec_profile <- paste(molec_profile, sepr, lw_23S_prof, sep = "")}
      }else
      {
        lw_23S_prof <- ""
      }
      
      lw_tetM <- as.character(Combined_Output.df[m, "tetM_result"])
      if(lw_tetM == "POS")
      {
        if(is.na(molec_profile))
        {molec_profile <- "tetM"
        }else
        {molec_profile <- paste(molec_profile, sepr, "tetM", sep = "")}
      }
      
      lw_tetO <- as.character(Combined_Output.df[m, "tetO_result"])
      if(lw_tetO == "POS")
      {
        if(is.na(molec_profile))
        {molec_profile <- "tetO"
        }else
        {molec_profile <- paste(molec_profile, sepr, "tetO", sep = "")}
      }
      
      lw_cat <- as.character(Combined_Output.df[m, "cat_result"])
      if(lw_cat == "POS")
      {
        if(is.na(molec_profile))
        {molec_profile <- "cat" 
        }else
        {molec_profile <- paste(molec_profile, sepr, "cat", sep = "")}
      }
      
      lw_folA <- as.character(Combined_Output.df[m, "folA_motifs"])
      lw_folA_result <- as.character(Combined_Output.df[m, "folA_result"])
      if(lw_folA != "WT")
      {
        if(is.na(molec_profile))
        {molec_profile <- paste("folA ", lw_folA, sep = "")
        }else
        {molec_profile <- paste(molec_profile, sepr, "folA ", lw_folA, sep = "")}
      }
      
      lw_folP <- as.character(Combined_Output.df[m, "folP_motifs"])
      lw_folP_result <- as.character(Combined_Output.df[m, "folP_result"])
      if(lw_folP_result == "NEG")
      {
        lw_folP <- "Err"
        lw_comments <- paste(lw_comments, "no folP gene", sep = "")
      }
      
      if(lw_folP != "WT")
      {
        if(is.na(molec_profile))
        {molec_profile <- paste("folP ", lw_folP, sep = "")
        }else
        {molec_profile <- paste(molec_profile, sepr, "folP ", lw_folP, sep = "")}
      }
       
      lw_pbp1a <- as.character(Combined_Output.df[m, "pbp1a_motifs"])
      lw_pbp1a_result <- as.character(Combined_Output.df[m, "pbp1a_result"])
      if(lw_pbp1a_result == "NEG")
      {
        lw_pbp1a <- "Err/Err/Err/Err"
        lw_comments <- paste(lw_comments, "pbp1a missing", sep = "")
      }
      
      if(lw_pbp1a != "WT/WT/WT/WT")
      {
        if(is.na(molec_profile))
        {molec_profile <- paste("pbp1a ", lw_pbp1a, sep = "")
        }else
        {molec_profile <- paste(molec_profile, sepr, "pbp1a ", lw_pbp1a, sep = "")}
      }
      
      lw_pbp2b <- as.character(Combined_Output.df[m, "pbp2b_motifs"])
      lw_pbp2b_result <- as.character(Combined_Output.df[m, "pbp2b_result"])
      if(lw_pbp2b_result == "NEG")
      {
        lw_pbp2b <- "Err/Err/Err/Err"
        lw_comments <- paste(lw_comments, "pbp2b missing", sep = "")
      }
      
      if(lw_pbp2b != "WT/WT/WT/WT")
      {
        if(is.na(molec_profile))
        {molec_profile <- paste("pbp2b ", lw_pbp2b, sep = "")
        }else
        {molec_profile <- paste(molec_profile, sepr, "pbp2b ", lw_pbp2b, sep = "")}
      }
      
      lw_pbp2x <- as.character(Combined_Output.df[m, "pbp2x_motifs"])
      lw_pbp2x_result <- as.character(Combined_Output.df[m, "pbp2x_result"])
      if(lw_pbp2x_result == "NEG")
      {
        lw_pbp2x <- "Err/Err/Err/Err"
        lw_comments <- paste(lw_comments, "pbp2x missing", sep = "")
      }
      
      if(lw_pbp2x != "WT/WT/WT/WT")
      {
        if(is.na(molec_profile))
        {molec_profile <- paste("pbp2x ", lw_pbp2x, sep = "")
        }else
        {molec_profile <- paste(molec_profile, sepr, "pbp2x ", lw_pbp2x, sep = "")}
      }
      
      lw_pbp_profile <- paste("PBP" ,as.character(Combined_Output.df[m, "pbp1a_allele"]),
                              as.character(Combined_Output.df[m, "pbp2b_allele"]),
                              as.character(Combined_Output.df[m, "pbp2x_allele"]), sep = ":")
      lw_gyrA <- as.character(Combined_Output.df[m, "gyrA_motifs"])
      lw_gyrA_result <- as.character(Combined_Output.df[m, "gyrA_result"])
      
      if(lw_gyrA_result == "NEG")
      {
        lw_gyrA <- "Err"
        lw_comments <- paste(lw_comments, "gyrA missing", sep = "")
      }
      
      if(lw_gyrA != "WT")
      {
        if(is.na(molec_profile))
        {molec_profile <- paste("gyrA ", lw_gyrA, sep = "")
        }else
        {molec_profile <- paste(molec_profile, sepr, "gyrA ", lw_gyrA, sep = "")}
      }
      
      lw_parC <-   as.character(Combined_Output.df[m, "parC_motifs"])
      lw_parC_result <- as.character(Combined_Output.df[m, "parC_result"])
      if(lw_parC_result == "NEG")
      {
        lw_parC <- "Err/Err/Err"
        lw_comments <- paste(lw_comments, "parC missing", sep = "")
      }
  
      lw_parC_parts <- unlist(strsplit(lw_parC, "/"))
      lw_parC_S79 <-  lw_parC_parts[1]
      lw_parC_D83 <- lw_parC_parts[2]
      lw_parC_N91 <- lw_parC_parts[3]
      lw_parC_prof <- NA
      
      if(lw_parC_result != "Sample_Err")
      {
        if(lw_parC_S79 != "WT")
        {
          if(is.na(lw_parC_prof))
          {lw_parC_prof <- paste("parC ", lw_parC_S79, sep = "")
          }else
          {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_S79, sep = "")}
        }
        
        if(lw_parC_D83 != "WT")
        {
          if(is.na(lw_parC_prof))
          {lw_parC_prof <- paste("parC ", lw_parC_D83, sep = "")
          }else
          {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_D83, sep = "")}
        }
          
        if(lw_parC_N91 != "WT")
        {
          if(is.na(lw_parC_prof))
          {lw_parC_prof <- paste("parC ", lw_parC_N91, sep = "")
          }else
          {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_N91, sep = "")}
        }
        
        if(!is.na(lw_parC_prof))
        {
          if(is.na(molec_profile))
          {molec_profile <- lw_parC_prof
          }else
          {molec_profile <- paste(molec_profile, sepr, lw_parC_prof, sep = "")}
        }
      }
      
      if((molec_profile == "") | is.na(molec_profile)) 
      {molec_profile <- "Wild Type"}

      #------------------------------------------------------------------------ Fluoroquinolone MICs and Interpretations
      
      lw_gyrA_value_S81F <- 0
      lw_gyrA_value_S81Y <- 0
      lw_gyrA_value_S81L <- 0
      lw_parC_value_S79_any <- 0
      lw_parC_value_D83_any <- 0
      
      if(lw_gyrA == "S81F")
      {lw_gyrA_value_S81F <- 1}
      if (lw_gyrA == "S81Y")
      {lw_gyrA_value_S81Y <- 1}
      if (lw_gyrA == "S81L")
      {lw_gyrA_value_S81L <- 1}
      if (lw_parC_S79 != "WT") #any mutation
      {lw_parC_value_S79_any <- 1}
      if (lw_parC_D83 != "WT") #any mutation
      {lw_parC_value_D83_any <- 1}

      lev_mic_inc <- round(-0.218 +
                          (lw_gyrA_value_S81F * 2.028)+
                          (lw_gyrA_value_S81Y * 1.564)+
                          (lw_gyrA_value_S81L * 3.564)+
                          (lw_parC_value_S79_any * 1.654)+
                          (lw_parC_value_D83_any * 0.834))
      lev_mic <- round( (2^lev_mic_inc), digits = 3)
      lev_mic2 <- paste(lev_mic, " ug/ml", sep = "")
      
      if(lev_mic < 2)
      {lev_mic2 <- paste("<= 1 ug/ml", sep = "")
      }else if(lev_mic > 16)
      {lev_mic2 <- paste(">= 32 ug/ml", sep = "")}
      
      lev_interp <- lev_mic.df$Interp[lev_mic.df$MIC == lev_mic]
      
      mox_mic_inc <-  round(-2.819 +
                           (lw_gyrA_value_S81F * 3.130)+
                           (lw_gyrA_value_S81Y * 3.907)+
                           (lw_gyrA_value_S81L * 4.907)+
                           (lw_parC_value_S79_any * 0.911))
      mox_mic <- round( (2^mox_mic_inc), digits = 3)
      mox_mic2 <- paste(mox_mic, " ug/ml", sep = "")
      
      if(mox_mic < 0.25)
      {mox_mic2 <- paste("<= 0.125 ug/ml", sep = "")
      }else if(mox_mic > 16)
      {mox_mic2 <- paste(">= 32 ug/ml", sep = "")}
      
      mox_interp <- mox_mic.df$Interp[mox_mic.df$MIC == mox_mic]
      
      if(lw_gyrA == "Err" | lw_parC == "Err/Err/Err" | lw_gyrA == "short" | lw_parC == "short")
      {
        lev_mic2 <- "Error"
        lev_interp <- "Error"
        mox_mic2 <- "Error"
        mox_interp <- "Error"
      }
      
      if(lw_gyrA_result=="NEG")
      {
        lw_gyrA <- "no gene"
        lw_comments <- paste(lw_comments, "gyrA missing.", sep="")
      }
      
      if (lw_parC_result == "NEG")
      {
        lw_parC <- "no gene"
        lw_parC_S79 <- "no gene"
        lw_parC_D83 <- "no gene"
        lw_parC_N91 <- "no gene"
        lw_comments <- paste(lw_comments, "parC missing.", sep = "")
      }
      
      #------------------------------------------------------------------------ Penicillin  MICs and Interpretations
      lw_pbp1a_parts <- unlist(strsplit(lw_pbp1a, "/")) #get motif regions
      lw_pbp1a_1 <- lw_pbp1a_parts[1]
      lw_pbp1a_2 <- lw_pbp1a_parts[2]
      lw_pbp1a_3 <- lw_pbp1a_parts[3]
      lw_pbp1a_4 <- lw_pbp1a_parts[4]
      if(is.na(lw_pbp1a_4)) #if truncated gene, last character is missing
      {lw_pbp1a_4 <- ""} 
      
      lw_pbp2b_parts <- unlist(strsplit(lw_pbp2b, "/")) #get motif regions
      lw_pbp2b_1 <- lw_pbp2b_parts[1]
      lw_pbp2b_2 <- lw_pbp2b_parts[2]
      lw_pbp2b_3 <- lw_pbp2b_parts[3]
      lw_pbp2b_4 <- lw_pbp2b_parts[4]
      if(is.na(lw_pbp2b_4)) #if truncated gene, last character is missing
      {lw_pbp2b_4 <- ""}
      
      lw_pbp2x_parts <- unlist(strsplit(lw_pbp2x, "/"))
      lw_pbp2x_1 <- lw_pbp2x_parts[1]
      lw_pbp2x_2 <- lw_pbp2x_parts[2]
      lw_pbp2x_3 <- lw_pbp2x_parts[3]
      lw_pbp2x_4 <- lw_pbp2x_parts[4]
      if(is.na(lw_pbp2x_4)) #if truncated gene, last character is missing
      {lw_pbp2x_4 <- ""} 
      
      lw_pbp1a_motif1_value_any <- 0
      lw_pbp1a_motif1_value_SAMK <- 0
      lw_pbp1a_motif1_value_SSMK <- 0
      lw_pbp1a_motif2_value_any <- 0
      lw_pbp1a_motif3_value_any <- 0
      lw_pbp1a_motif4_value_any <- 0
      
      lw_pbp2b_motif1_value_any <- 0
      lw_pbp2b_motif2_value_any <- 0
      lw_pbp2b_motif3_value_any <- 0
      lw_pbp2b_motif4_value_any <- 0
      
      lw_pbp2x_motif1_value_any <- 0
      lw_pbp2x_motif1_value_SAFK <- 0
      lw_pbp2x_motif1_value_other <- 0
      lw_pbp2x_motif2_value_any <- 0
      lw_pbp2x_motif3_value_any <- 0
      lw_pbp2x_motif3_value_EDT <- 0
      lw_pbp2x_motif3_value_KEA <- 0
      lw_pbp2x_motif4_value_VKSG <- 0
      lw_pbp2x_motif4_value_any <- 0
      
      if(lw_pbp1a_1 != "WT") #any mutation
      {lw_pbp1a_motif1_value_any <- 1}
      if(lw_pbp1a_1 == "SAMK")
      {lw_pbp1a_motif1_value_SAMK <- 1}
      if(lw_pbp1a_1 == "SSMK")
      {lw_pbp1a_motif1_value_SSMK <- 1}
      if(lw_pbp1a_4 != "WT") #any mutation
      {lw_pbp1a_motif4_value_any <- 1}
      if(lw_pbp2b_2 != "WT") #any mutation
      {lw_pbp2b_motif2_value_any <- 1}
      if(lw_pbp2b_3 != "WT") #any mutation
      {lw_pbp2b_motif3_value_any <- 1}
      
      if(lw_pbp2x_1 == "SAFK")
      {
        lw_pbp2x_motif1_value_SAFK <- 1
      }else
      {
        if(lw_pbp2x_1 != "WT")
        {lw_pbp2x_motif1_value_other <- 1}
      }
      
      if(lw_pbp2x_3 == "EDT")
      {lw_pbp2x_motif3_value_EDT <- 1}
      if (lw_pbp2x_3 == "KEA")
      {lw_pbp2x_motif3_value_EDT <- 1}
      if (lw_pbp2x_4 == "VKSG")
      {lw_pbp2x_motif4_value_VKSG <- 1}
      
      pen_mic_inc <-  round(-4.6099 +
                           (lw_pbp1a_motif1_value_any * 1.5466)+
                           (lw_pbp1a_motif4_value_any * 0.9491)+
                           (lw_pbp2b_motif2_value_any * 1.2025)+
                           (lw_pbp2b_motif3_value_any * 0.3563)+
                           (lw_pbp2x_motif1_value_SAFK * 1.6259)+
                           (lw_pbp2x_motif3_value_EDT * 1.5476)+
                           (lw_pbp2x_motif3_value_KEA * 0.6760)+
                           (lw_pbp2x_motif4_value_VKSG * 0.7536))
      pen_mic <- round( (2^pen_mic_inc), digits = 3)
      pen_mic2 <- paste(pen_mic, " ug/ml", sep = "")
      
      if(pen_mic < 0.032)
      {pen_mic2 <- paste("<= 0.03 ug/ml", sep = "")
      }else if(pen_mic > 3)
      {pen_mic2 <- paste(">= 4 ug/ml", sep = "")}
      
      pen_interp <- pen_mic.df$Interp[pen_mic.df$MIC == pen_mic]
      
      cro_mic_inc <- round(-2.709 +
                          (lw_pbp1a_motif1_value_any * 1.25)+
                          (lw_pbp2x_motif1_value_SAFK * 2.72)+
                          (lw_pbp2x_motif3_value_EDT * 0.76)+
                          (lw_pbp2x_motif4_value_VKSG * 0.989))
      
      cro_mic <- round((2^cro_mic_inc), digits = 3)
      cro_mic2 <- paste(cro_mic, " ug/ml", sep = "")
      
      if(cro_mic < 0.13)
      {cro_mic2 <- paste("<= 0.125 ug/ml", sep = "")
      }else if (cro_mic > 3)
      {cro_mic2 <- paste(">= 4 ug/ml", sep = "")}
      
      cro_interp <- cro_mic.df$Interp[cro_mic.df$MIC == cro_mic]
      
      cfm_mic_inc <- round(-1.018 +
                          (lw_pbp1a_motif1_value_SAMK * 1.509)+
                          (lw_pbp1a_motif1_value_SSMK * 2.170)+
                          (lw_pbp2x_motif1_value_SAFK * 2.322)+
                          (lw_pbp2x_motif1_value_other * 0.256)+
                          (lw_pbp2x_motif4_value_VKSG * 1.026))
      
      cfm_mic <- round( (2^cfm_mic_inc), digits = 3)
      cfm_mic2 <- paste(cfm_mic, " ug/ml", sep = "")
      
      if (cfm_mic < 0.51)
      {cfm_mic2 <- paste("<= 0.5 ug/ml", sep = "")
      }else if (cfm_mic > 15)
      {cfm_mic2 <- paste(">= 16 ug/ml", sep = "")}
      
      cfm_interp <- cfm_mic.df$Interp[cfm_mic.df$MIC == cfm_mic]
      
      if ((lw_pbp1a == "Err/Err/Err/Err") | (lw_pbp2b == "Err/Err/Err/Err") | (lw_pbp2x == "Err/Err/Err/Err"))
      {
        pen_mic2 <- "Error"
        cro_mic2 <- "Error"
        cfm_mic2 <- "Error"
        pen_interp <- "Error"
        cro_interp <- "Error"
        cfm_interp <- "Error"
      }
      
      if(lw_pbp1a_result == "NEG")
      {
        lw_pbp1a <- "no gene"
        lw_comments <- paste(lw_comments, "pbp1a missing.", sep = "")
      }
      
      if(lw_pbp2b_result == "NEG")
      {
        lw_pbp2b <- "no gene"
        lw_comments <- paste(lw_comments, "pbp2b missing.", sep = "")
      }
      
      if(lw_pbp2x_result == "NEG")
      {
        lw_pbp2x <- "no gene"
        lw_comments <- paste(lw_comments, "pbp2x missing.", sep = "")
      }
      
      #------------------------------------------------------------------------ Trimethoprim/Sulfamethoxazole
      lw_folA_value <- 0
      lw_folP_value <- 0
      
      if(lw_folA == "I100L")
      {lw_folA_value <- 1}
      if(lw_folP != "WT")
      {lw_folP_value <- 1}
      
      sxt_mic_inc <- round(-2.265 +
                          (lw_folA_value * 2.113)+
                          (lw_folP_value * 2.668))
      sxt_mic <- round( (2^sxt_mic_inc), digits = 3)
      sxt_mic2 <- paste(sxt_mic, " ug/ml", sep = "")
      
      if(sxt_mic < 0.51)
      {sxt_mic2 <- paste("<= 0.5/9.5 ug/ml", sep = "")
      }else if (sxt_mic > 3)
      {sxt_mic2 <- paste(">= 4/76 ug/ml", sep = "")}
      
      sxt_interp <- sxt_mic.df$Interp[sxt_mic.df$MIC == sxt_mic]
      
      if((lw_folA == "Err") | (lw_folP == "Err"))
      {
        sxt_mic2 <- "Error"
        sxt_interp <- "Error"
      }
      
      if(lw_folA_result == "NEG")
      {
        lw_folA <- "no gene"
        lw_comments <- paste(lw_comments, "folA missing.", sep = "")
      }
      if(lw_folP_result == "NEG")
      {
        lw_folP <- "no gene"
        lw_comments <- paste(lw_comments, "folP missing.", sep = "")
      }
      
      #------------------------------------------------------------------------ MACROLIDEs MICs and Interpretations
      
      # lw_23S_A2059G set above
      # lw_23s_C2611T set above
      lw_ermB_value <- 0
      lw_ermTR_value <- 0
      lw_mefAE_value <- 0
      
      if (lw_ermB == "POS")
      {lw_ermB_value <- 1} # if ermB(S), then lw_ermB is set to NEG above
      if (lw_ermTR == "POS")
      {lw_ermTR_value <- 1}
      if (lw_mefAE == "POS")
      {lw_mefAE_value <- 1}
      
      # if both ermB and mef are positive, only use ermB value to calculate
      if ((lw_ermB == "POS") && (lw_mefAE == "POS"))
      {lw_mefAE_value <- 0}
      
      ery_mic_inc <-  round(-2.975 +
                           (lw_23S_A2059G * 1.993)+
                           (lw_ermB_value * 7.680)+
                           (lw_mefAE_value * 4.808))
      
      ery_mic <- round((2^ery_mic_inc), digits = 5)
      ery_mic2 <- paste(ery_mic, " ug/ml", sep = "")
      
      if(ery_mic <= 0.25)
      {ery_mic2 <- paste("<= 0.125 ug/ml", sep = "")
      }else if(ery_mic > 16)
      {ery_mic2 <- paste(">= 32 ug/ml", sep = "")}

      ery_interp <- ery_mic.df$Interp[ery_mic.df$MIC == ery_mic]
      
      #--------
      azi_mic_inc <- round(-1.9722 +
                          (lw_23S_A2059G * 1.1204)+
                          (lw_23S_C2611T * 0.9917)+
                          (lw_ermB_value * 3.9722)+
                          (lw_mefAE_value * 3.8122))
      
      azi_mic <- round((2^azi_mic_inc), digits = 5)
      azi_mic2 <- paste(azi_mic, " ug/ml", sep = "")
      
      if(azi_mic <= 0.25)
      {azi_mic2 <- paste("<= 0.125 ug/ml", sep = "")
      }else if(azi_mic >= 2)
      {azi_mic2 <- paste(">= 2 ug/ml", sep = "")}

      azi_interp <- azi_mic.df$Interp[azi_mic.df$MIC == azi_mic]
      #--------
      
      cla_mic_inc <- round(-4.984 +
                          (lw_23S_A2059G * 1.819)+
                          (lw_23S_C2611T * 1.246)+
                          (lw_ermB_value * 10.820)+
                          (lw_mefAE_value * 6.064))
      
      cla_mic <- round((2^cla_mic_inc), digits = 5)
      cla_mic2 <- paste(cla_mic, " ug/ml", sep = "")
      
      if (cla_mic < 0.06)
      {cla_mic2 <- paste("<= 0.03 ug/ml", sep = "")
      }else if(cla_mic > 32)
      {cla_mic2 <- paste(">= 64 ug/ml", sep = "")}
      
      cla_interp <- cla_mic.df$Interp[cla_mic.df$MIC == cla_mic]
      
      cli_mic_inc <- round(-2.8145 +
                          (lw_23S_A2059G * 0.456)+
                          (lw_ermB_value * 9.048))
      
      cli_mic <- round((2^cli_mic_inc), digits = 5)
      
      cli_mic2 <- paste(cli_mic, " ug/ml", sep = "")
      
      if (cli_mic <= 0.25)
      {cli_mic2 <- paste("<= 0.125 ug/ml", sep = "")
      }else if(cli_mic > 32)
      {cli_mic2 <- paste(">= 64 ug/ml", sep = "")}

      cli_interp <- cli_mic.df$Interp[cli_mic.df$MIC == cli_mic]
      if (lw_ermTR == "POS")
      {cli_interp <- "Inducible"}

      if(lw_ermB == "Sample_Err")
      {
        ery_mic2 <- "Sample_Err"
        ery_interp <- "Sample_Err"
        cli_mic2 <- "Sample_Err"
        cli_interp <- "Sample_Err"
        cla_mic2 <- "Sample_Err"
        cla_interp <- "Sample_Err"
        azi_mic2 <- "Sample_Err"
        azi_interp <- "Sample_Err"
      }
      
      if(lw_ermB == "NF")
      {
        ery_mic2 <- "Error"
        ery_interp <- "Curation Error"
        cli_mic2 <- "Error"
        cli_interp <- "Curation Error"
        cla_mic2 <- "Error"
        cla_interp <- "Curation Error"
        azi_mic2 <- "Error"
        azi_interp <- "Curation Error"
      }
      
      
      if(lw_23S_prof == "23S Err")
      {
        lw_23S_A2059G <- 99L
        lw_23S_C2611T <- 99L
        ery_mic2 <- "Error"
        ery_interp <- "No VCF"
        cli_mic2 <- "Error"
        cli_interp <- "No VCF"
        cla_mic2 <- "Error"
        cla_interp <- "No VCF"
        azi_mic2 <- "Error"
        azi_interp <- "No VCF"
      }

      #-----------------------------------------------------------------------------------------
      chl_mic <- 4
      if(lw_cat == "POS"){chl_mic <- 8} #list in increasing mic order
      
      if(chl_mic == 4)
      {chl_mic2 <- paste("<= 4 ug/ml", sep = "")
      }else if (chl_mic == 8)
      {chl_mic2 <- paste(">= 8 ug/ml", sep = "")}
      
      chl_interp <- chl_mic.df$Interp[chl_mic.df$MIC == chl_mic]
      
      if(lw_cat == "Sample_Err")
      {
        chl_mic2 <- "Sample_Err"
        chl_interp <- "Sample_Err"
      }
      
      #-----------------------------------------------------------------------------------------
      
      tet_mic <- 1
      if(lw_tetM == "POS"){tet_mic <- 16} #list in increasing mic order
      
      if(tet_mic == 1)
      {tet_mic2 <- paste("<= 1 ug/ml", sep = "")
      }else if (tet_mic == 16)
      {tet_mic2 <- paste(">= 16 ug/ml", sep = "")}
      
      tet_interp <- tet_mic.df$Interp[tet_mic.df$MIC == tet_mic]
      dox_mic <- 0.25
      
      if(lw_tetM == "POS")
      {dox_mic <- 4} #list in increasing mic order
      
      if(dox_mic == 0.25)
      {dox_mic2 <- paste("<= 0.25 ug/ml", sep = "")
      }else if(dox_mic == 4)
      {dox_mic2 <- paste(">= 4 ug/ml", sep = "")}
      
      dox_interp <- dox_mic.df$Interp[dox_mic.df$MIC == dox_mic]
      
      if(lw_tetM == "Sample_Err")
      {
        tet_mic2 <- "Sample_Err"
        tet_interp <- "Sample_Err"
        dox_mic2 <- "Sample_Err"
        dox_interp <- "Sample_Err"
      }
      
      #--------------------------------------------------------------  MAKE AMR PROFILE
      amr_profile <- "Susceptible"
    
      if((str_detect(molec_profile, "folA Err")) | (str_detect(molec_profile, "folP Err")))
      {
        amr_profile <- "Error"
        sxt_interp <- "Error"
      }
      
      if((str_detect(molec_profile, "gyrA Err")) | (str_detect(molec_profile, "parC Err")))
      {
        lev_interp <- "Error"
        mox_interp <- "Error"
        amr_profile <- "Error"
      }
      
      if(str_detect(molec_profile, "Err"))
      {amr_profile <- "Error"}
      
      if(lw_ermB == "Sample_Err")
      {amr_profile <- "Sample Error"}

      if(amr_profile == "Susceptible")
      {
        sepr2 <- "/"
        
        if(pen_interp == "Intermediate")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "PEN-I"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "PEN-I", sep = "")}
        }
        
        if(pen_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "PEN-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "PEN-R", sep = "")}
        }
        
        if(cro_interp == "Intermediate")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "CRO-I"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "CRO-I", sep = "")}
        }
        
        if(cro_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "CRO-R"
          } else
          {amr_profile <- paste(amr_profile, sepr2, "CRO-R", sep = "")}
        }
        
        if(cfm_interp == "Intermediate")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "CFM-I"
          } else
          {amr_profile <- paste(amr_profile, sepr2, "CFM-I", sep = "")}
        }
        
        if(cfm_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "CFM-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "CFM-R", sep = "")}
        }
        
        if(ery_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "ERY-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "ERY-R", sep = "")}
        }
        
        if(azi_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "AZI-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "AZI-R", sep = "")}
        }
        
        if(cla_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "CLA-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "CLA-R", sep = "")}
        }
        
        if(cli_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "CLI-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "CLI-R", sep = "")}
        }
        
        if(cli_interp == "Inducible")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "CLI-Ind"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "CLI-Ind", sep = "")}
        }
        
        if(chl_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "CHL-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "CHL-R", sep = "")}
        }
        
        if(lev_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "LEV-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "LEV-R", sep = "")}
        }
        
        if(mox_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "MOX-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "MOX-R", sep = "")}
        }
        
        if(tet_interp == "Resistant")
        {
          if (amr_profile == "Susceptible")
          {amr_profile <- "TET-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "TET-R", sep = "")}
        }
        
        if(dox_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "DOX-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "DOX-R", sep = "")}
        }
        
        if(sxt_interp == "Resistant")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "SXT-R"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "SXT-R", sep = "")}
        }
        
        if(sxt_interp == "Intermediate")
        {
          if(amr_profile == "Susceptible")
          {amr_profile <- "SXT-I"
          }else
          {amr_profile <- paste(amr_profile, sepr2, "SXT-I", sep = "")}
        }
        
      }
      
    } else
    {
      lw_pbp1a <- "Sample Error"
      lw_pbp2b <- "Sample Error"
      lw_pbp2x <- "Sample Error"
      lw_23S_A2059G <- 99L
      lw_23S_C2611T <- 99L
      lw_ermB <- "Sample Error"
      lw_ermTR <- "Sample Error"
      lw_mefAE <- "Sample Error"
      lw_folA <- "Sample Error"
      lw_folP <- "Sample Error"
      lw_gyrA <- "Sample Error"
      lw_parC <- "Sample Error"
      lw_tetM <- "Sample Error"
      lw_tetO <- "Sample Error"
      lw_cat <- "Sample Error"
      molec_profile <- "Sample Error"
      pen_mic2 <- "Sample Error"
      pen_interp <- "Sample Error"
      cro_mic2 <- "Sample Error"
      cro_interp <- "Sample Error"
      cfm_mic2 <- "Sample Error"
      cfm_interp <- "Sample Error"
      azi_mic2 <- "Sample Error"
      azi_interp <- "Sample Error"
      ery_mic2 <- "Sample Error"
      ery_interp <- "Sample Error"
      cla_mic2 <- "Sample Error"
      cla_interp <- "Sample Error"
      cli_mic2 <- "Sample Error"
      cli_interp <- "Sample Error"
      chl_mic2 <- "Sample Error"
      chl_interp <- "Sample Error"
      lev_mic2 <- "Sample Error"
      lev_interp <- "Sample Error"
      mox_mic2 <- "Sample Error"
      mox_interp <- "Sample Error"
      tet_mic2 <- "Sample Error"
      tet_interp <- "Sample Error"
      dox_mic2 <- "Sample Error"
      dox_interp <- "Sample Error"
      sxt_mic2 <- "Sample Error"
      sxt_interp <- "Sample Error"
      amr_profile <- "Sample Error"
      lw_allele_profile <- "Contig file not found"
    }
    
    #-------------------------------------------------------------------------- New LabWare Uploader structure
    sample_data.df <- tibble(lw_CurrSampleNo, lw_pbp1a, lw_pbp2b, lw_pbp2x,
                             lw_23S_A2059G, lw_23S_C2611T,  lw_ermB, lw_ermTR, lw_mefAE,
                             lw_folA, lw_folP,
                             lw_gyrA, lw_parC,
                             lw_tetM, lw_tetO, lw_cat, molec_profile,
                             pen_mic2, pen_interp,
                             cro_mic2, cro_interp,
                             cfm_mic2, cfm_interp,
                             azi_mic2, azi_interp,
                             ery_mic2, ery_interp,
                             cla_mic2, cla_interp,
                             cli_mic2, cli_interp,
                             chl_mic2, chl_interp,
                             lev_mic2, lev_interp,
                             mox_mic2, mox_interp,
                             tet_mic2, tet_interp,
                             dox_mic2, dox_interp,
                             sxt_mic2, sxt_interp,
                             amr_profile,
                             lw_allele_profile)
    #---------------------------------------------------------------------------------------
    if(m == 1)
    {lw_Output.df <- tibble(sample_data.df)
    } else 
    {lw_Output.df <- bind_rows(lw_Output.df, sample_data.df)}
  }#end samples loop

  write.csv(lw_Output.df, paste(local_output_dir, "output_profile_interpret_AMR.csv", sep = ""), quote = FALSE, row.names = FALSE)
  cat("\n\n AMR interpretation analysis done! output_AMR_interpret is ready in the output folder. \n\n\n")
  view(lw_Output.df)
 
} # end if LocusID != list
