#' Molecular typing pipeline for WGS assemblies
#' February 5 2024, Walter Demczuk & Shelley Peterson
#'
#' Takes Organism, Sample Number, Locus, to query a contig.fasta file
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param Test_id AMR, TOXINS, VIRULENCE, NGSTAR, use MASTER for 16S id
#' @param SampleNo Sample number associated with contig.fasta file
#' @param LocusID The locus to query, or enter "list" to use a list of alleles
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details GAS AMR loci:
#'
#' dfrF, dfrG, DHFR, ermA, ermB, ermT, folA, folP, gyrA, mefAE, msrD, parC, rpsJ, tetM, tetO
#' note that ermA == ermTR, ermC == ermT, mefA/E includes mefA and mefE; DHFR is the same as folA
#' ermB + = ERY-R/CLI-R; ermT = ERY-R/CLI-?; mefA + = ERY-R/CLI-S; ermA/TR = ERY-R/CLI-Ind.
#'
#' PNEUMO AMR loci:
#' cat ermB, ermTR, folA, folP, gyrA, parC, mefA, mefE, msrD, tetM, tetO, pbp1a, pbp2b, pbp2x
#' note that ermA == ermTR
#'
#' GAS VIRULENCE FACTORS:
#' covR => V128A
#' covS => I332V; E226G (1228 deletion)
#' rocA => V47A; R259Q; V333A; D396N; I404T; P405S
#' rpoB/rgg => V169I
#' nga => promoter mutations: A8G (K3R) or T13G (F5V) or C17T (A6V) (protein): F5V+A6V=pnga3
#' nga => D349G is the D330G mutation
#' cpA, fctA, fctB, lepA and srtC1 all pilin proteins - use only fctA (major pilin protein)
#' fbp54, grab, ideS_mac, mf_spd, speB, scpA, ska always positive?
#' hasA,B,C are the same operon, use only hasA
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#   Org_id <- "PNEUMO"                #GAS, GBS, PNEUMO or GONO
#   Test_id <- "AMR"                  #AMR, TOXINS, VIRULENCE, NGSTAR, rRNA_16S
#   SampleNo <- "SC23-8992-A"
#   LocusID <- "list"
#   curr_work_dir <- "C:\\WADE\\"
#   Blast_evalue <- "10e-50"          #sets sensitivity of Blast gene match 10e-50 to 10e-150; use 10e-5 for primers
#-------------------------------------------------------------------------------

MASTER_pipeline <- function(Org_id, Test_id, SampleNo, LocusID, curr_work_dir){

  start_time <- Sys.time()
  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  Variable <- NA
  Test_id[Org_id == "GONO" & 
          (LocusID %in% c("penA", "mtrR", "porB", "ponA", "gyrA", "parC", "rRNA23S")) &
          Test_id == "AMR_ALL"] <- "NGSTAR"
  Test_id[Test_id == "AMR_ALL"] <- "AMR"
  Test_id[Org_id == "PNEUMO" &
          Test_id %in% c("AMR_ALL", "AMR_2")] <- "AMR"

  directorylist <- getdirectory(curr_work_dir, Org_id, Test_id)
  reflist <- refdirectory(directorylist, Org_id, Test_id)
  
  # create a fasta file for all sequences output_dna.fasta and output_aa.fasta
  dna_file <- paste0(directorylist$output_dir, "output_dna.fasta")
  dna_nf_file <- paste0(directorylist$output_dir, "output_dna_notfound.fasta")
  aa_file <- paste0(directorylist$output_dir, "output_aa.fasta")
  #-----------------------------------------------------------------------------
  # Blast Evalues
  blast_evalues.df <- as_tibble(read.csv(paste0(reflist$Ref_Dir, "blast_evalues.csv"),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Blast_evalue <- as.character(blast_evalues.df$contig[1])
  evalue_allele <- as.character(blast_evalues.df$allele[1])
  blast_wt_id_threshold <- as.numeric(blast_evalues.df$wt_id[1])
  #-----------------------------------------------------------------------------
  
  ########################### Setup sample list table ##########################
  if(SampleNo == "list")
  {
    SampleList.df <-as_tibble(read.csv(directorylist$SampleList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  }else
  {
    SampleList.df <- tibble(SampleNo, Variable)
  }
  
  NumSamples <- dim(SampleList.df)[1]
  
  ########################### Setup locus list table ###########################
  if(LocusID == "list")
  {
    LocusList.df <- read.csv(reflist$Loci_List, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    LocusList.df <- as_tibble(LocusList.df)
  }else
  {
    LocusList.df <- tibble(LocusID)
  }
  NumLoci <- dim(LocusList.df)[1]

  LocusMutationsFile <- paste0(reflist$Ref_Dir, "loci_mutations.csv")
  if(file.exists(LocusMutationsFile))
  {
    LocusMutations.df <- as_tibble(read.csv(LocusMutationsFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
    Mutns <- "Yes"
  } else 
  {
    Mutns <- "No"
  }

  if(Test_id == "MASTERBLASTR")
  {
    quickBlastIndex(NumLoci, LocusList.df, reflist)
  }

  #Progress Bar
  n_iter <- NumSamples
  init <- numeric(n_iter)
  end <- numeric(n_iter)
  
  m <- 1L
  for(m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Start Sample Loop
  {
    init[m] <- Sys.time()
    
    sampleinfo <- getsampleinfo(SampleList.df, directorylist, reflist, m)
    CurrSample.df <- filter(SampleList.df, SampleNo == sampleinfo$CurrSampleNo)
    SampleProfile <- ""
    cat(m, " of ", NumSamples, "\n")
    
    if(sampleinfo$Allele[1] != "Sample_Err"){
      #make blast database of contig file
      BlastFormatCommand <- paste0("makeblastdb -in ", reflist$DestFile, " -dbtype nucl")
      try(system(BlastFormatCommand))
    }

    p <- 1L
    for(p in 1L:NumLoci)   # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Start Locus Loop
    {
      AlleleLine <- ""
      AlleleInfo <- ""
      AlleleInfo[1:4] <- ""
      blast_parsed <- ""
      blastAA <- ""
      IDpercent2 <- ""
      BlastResult <- NA
      motifs <- ""

      ##########################################################################
      # BLAST and parse blastout.txt
      ##########################################################################
      CurrLocus <- as.character(LocusList.df[p,1])
      LocusLkupDNA <- paste0(reflist$Lkup_Dir, CurrLocus, ".fasta")
      LocusLkupDNApresent <- ifelse(file.exists(LocusLkupDNA), TRUE, FALSE)
      
      LocusFile <- paste0(reflist$WT_Dir, CurrLocus, ".fasta")
      if(!file.exists(LocusFile))
      {
        sample_error.df <- tibble(Locus_ID = CurrLocus, Output = "File Not Found", stringsAsFactors = FALSE)
        return(sample_error.df)
      }
      
      AlleleInfo[1][sampleinfo$Allele == "Sample_Err"] <- "Sample_Err"
      if(AlleleInfo[1] != "Sample_Err") #....................................... Not Sample Error
      { 
        BlastResult <- blasthits(reflist, directorylist, LocusFile, Blast_evalue)
        AlleleInfo[1] <- ifelse(BlastResult == "POS", "POS", "NEG")
        
        if(BlastResult == "POS") #============================================== Positive BLAST Result Loop
        {
          # Parse out to get WT and query DNA sequences, match % etc.
          blastoutput <- readLines(paste0(directorylist$temp_dir, "blastout.txt"))
          blast_parsed <- parseblast(blastoutput, LocusID)
          
          # Print both WT and query DNA sequences
          if(LocusID != "list")
          {
            cat("\n\n>", CurrLocus, "(Wildtype)\n", blast_parsed$WTDNASeqLine_NoDash_str, "\n", sep ="")
            cat("\n\n>", sampleinfo$CurrSampleNo, "_", CurrLocus , "\n", blast_parsed$DNASeqLine_NoDash_str, "\n", sep ="")
          }
          
          WTDNASeqLine <- DNAString(blast_parsed$WTDNASeqLine_str[1])
          WTDNASeqLine_NoDash <- DNAString(blast_parsed$WTDNASeqLine_NoDash_str[1])
          DNASeqLine <- DNAString(blast_parsed$DNASeqLine_str[1])
          DNASeqLine_NoDash <- DNAString(blast_parsed$DNASeqLine_NoDash_str[1])

          ######################################################################
          # Make Protein Sequence
          ######################################################################
          blastAA <- AAconvert(WTDNASeqLine_NoDash, DNASeqLine_NoDash)
          
          if(LocusID != "list")
          {
            cat("\n\n>", CurrLocus , "(Wildtype)\n", blastAA$WTAASeqLine_str, "\n", sep ="")
          }
          
          if (LocusID != "list")
          {
            cat("\n\n>", sampleinfo$CurrSampleNo, "_", CurrLocus , "\n", blastAA$AASeqLine_str, "\n", sep ="")
          }

          # Align protein sequences - make sure AASeqLine, WTAASeqLine are listed in correct order!
          WTAASeqLine <- AAString(blastAA$WTAASeqLine_str)
          AASeqLine <- AAString(blastAA$AASeqLine_str)
          globalAlign_AA <- pairwiseAlignment(AASeqLine, WTAASeqLine, substitutionMatrix = "BLOSUM50", 
                                              gapOpening = -2, gapExtension = -8, scoreOnly = FALSE)
          WTAASeqLine_aln <- subject(globalAlign_AA)
          WTAASeqLine_aln_str <- toString(WTAASeqLine_aln)
          AASeqLine_aln <- pattern(globalAlign_AA)
          AASeqLine_aln_str <- toString(AASeqLine_aln)

          # Create and print DNA "alignment" sequence that only shows mutations
          DNASeqLine_aln <- ""
          for(j in 1:str_length(blast_parsed$WTDNASeqLine_str))
          {
            if(str_sub(blast_parsed$WTDNASeqLine_str, j, j) == str_sub(blast_parsed$DNASeqLine_str, j, j))
            {
              DNASeqLine_aln <- paste0(DNASeqLine_aln, ".")
            } else
            {
              DNASeqLine_aln <- paste0(DNASeqLine_aln, str_sub(blast_parsed$DNASeqLine_str, j, j))
            }
          }

          if(LocusID != "list")
          {
            cat("\n\nDNA Alignment:\n", DNASeqLine_aln, "\n", sep = "")
          }
          
          ######################################################################
          # Mutations
          ######################################################################
          mutationspresent <- ""
          # Create and print AA "alignment" sequence that only shows mutations
          AASeqLine_aln_disp <- ""
          mutations <- ""

          for(k in 1:str_length(WTAASeqLine_aln_str))
          {
            if(str_sub(WTAASeqLine_aln_str, k, k) == str_sub(AASeqLine_aln_str, k, k))
            {
              AASeqLine_aln_disp <- paste0(AASeqLine_aln_disp, ".")
            } else
            {
              AASeqLine_aln_disp <- paste0(AASeqLine_aln_disp, str_sub(blastAA$AASeqLine_str, k, k))
              mutations <- paste0(mutations, str_sub(WTAASeqLine_aln_str, k, k), k, str_sub(blastAA$AASeqLine_str, k, k), " ")
            }
          }
          
          if (Mutns == "Yes")
          {
            motifs <- getmotifs(LocusMutations.df, Test_id, CurrLocus, LocusID, AASeqLine_aln_str, AASeqLine_aln_disp)
          }
          
          #Check for disrupted genes
          AlleleInfo[4] <- ifelse(str_detect(AASeqLine_aln_disp, "[*]"), "Disrupted", 
                                  ifelse(blast_parsed$IDcoverage < 80, "Disrupted", ""))
          
          if(LocusID != "list")
          {
            cat("\n\nProtein Alignment:\n", AASeqLine_aln_disp, "\n", sep = "")
            cat("\nMutations: ", mutations, "\n\n", sep = "")
            cat("\nMotifs: ", motifs, "\n\n", sep = "")
          }

          ######################################################################
          # Lookup DNA alleles
          ######################################################################
          # write DNASeqLine_NoDash_str to a file named = querygene.fasta,
          # BLAST vs. lookup table,
          # parse out the allele numbers

          Seq_File <- paste0(directorylist$temp_dir, "querygene.fasta")
          sink(Seq_File, split=FALSE, append = FALSE)
          cat(">", sampleinfo$CurrSampleNo, "_", CurrLocus , "\n", blast_parsed$DNASeqLine_NoDash_str, sep ="")
          sink()

          ExactMatchFound <- FALSE

          if(LocusLkupDNApresent) #********************************************* Locus Lookup DNA Loop
          {
            #BLAST lookup table
            BlastCommand <- paste0("blastn -query ", Seq_File, 
                                   " -db ", LocusLkupDNA, 
                                   " -out ", directorylist$temp_dir, "blastout2.txt", 
                                   " -evalue ", evalue_allele, 
                                   " -num_alignments 1")
            system(BlastCommand, intern = TRUE)
        
            blastoutput2 <- readLines(paste0(directorylist$temp_dir, "blastout2.txt"))

            querylength <- grep("Length=", blastoutput2, value = TRUE)
            querylength <- sub("Length=", "", querylength)
            lengthmatch <- ifelse(querylength[1] == querylength[2], TRUE, FALSE)
            IDLine2 <- grep("Identities", blastoutput2, value = TRUE)
            ExactMatchFound <- sum(str_detect(IDLine2, "(100%)"))
            
            if(lengthmatch == TRUE & ExactMatchFound > 0)
            {
              AlleleLine <- grep(">", blastoutput2, value = TRUE)
              AlleleParts <- unlist(strsplit(AlleleLine, "_"))
              AlleleInfo[2] <- AlleleParts[2]
              AlleleInfo[3] <- AlleleParts[3]
              AlleleInfo[4] <- AlleleParts[4]
            } else
            {
              AlleleInfo[2] <- "NF"
              AlleleInfo[3] <- "???"
              AlleleInfo[4] <- ""
            }
          } #******************************************************************* End Locus Lookup DNA Loop
      
          ######################################################################
          # Create Sample Output
          ######################################################################
          sink(dna_file, split=FALSE, append = TRUE)
          cat(">", CurrLocus, "_", sampleinfo$CurrSampleNo, "_", CurrLocus, AlleleInfo[2], "_", 
              AlleleInfo[3], "_", sampleinfo$CurrSampleVar, "_", motifs, "\n", blast_parsed$DNASeqLine_NoDash_str, "\n", sep ="")
          sink()

          if(AlleleInfo[2] == "NF")
          {
            sink(dna_nf_file, split=FALSE, append = TRUE)
            cat(">", CurrLocus, "_", motifs, "_", sampleinfo$CurrSampleNo, "_", 
                sampleinfo$CurrSampleVar, "_", "\n", blast_parsed$DNASeqLine_NoDash_str, "\n", sep ="")
            sink()
          }

          sink(aa_file, split=FALSE, append = TRUE)
          cat(">", CurrLocus, "_", sampleinfo$CurrSampleNo, "_", 
              sampleinfo$CurrSampleVar, "\n", blastAA$AASeqLine_str, "\n", sep ="")
          sink()
      
          #make a blast ID line based on full WT gene
          IDpercent2 <- paste0(blast_parsed$IDlength, "/", blast_parsed$WTlength, 
                              " (", blast_parsed$IDpercent, ")") 

          # special notes for porB or poor quality BLAST
          AlleleInfo[4][CurrLocus == "porB" & blast_parsed$IDcoverage < 90] <- "porB1a"
          AlleleInfo[4][CurrLocus == "porB" & blast_parsed$IDcoverage >= 90] <- "porB1b"
          AlleleInfo[1][blast_parsed$IDcoverage <= blast_wt_id_threshold] <- "NEG"
        } else #================================================================ End Positive BLAST Result Loop
        {
          blast_parsed <- tibble(IDLine = "",
                                 IDLine_trimmed = "")
          IDpercent2 <- ""
          motifs <- ""
        }
      } else #.................................................................. Close Not Sample Error Loop
      {
        AlleleInfo[2] <- "Sample_Err"
        AlleleInfo[3] <- "Sample_Err"
        AlleleInfo[4] <- "Sample_Err"
        blast_parsed <- tibble(IDLine = "Sample_Err",
                                IDLine_trimmed = "Sample_Err")
        IDpercent2 <- "Sample_Err"
        motifs <- "Sample_Err"
      }
      
      headers <- c("result", "allele", "mutations", "comments", "BlastID",
                   "PctWT", "motifs")
      headers <- paste(CurrLocus, headers, sep = "_")

      if(p==1)
      {
        OutputLocus.df <- tibble(AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], blast_parsed$IDLine_trimmed, IDpercent2, motifs)
        names(OutputLocus.df) <- headers
      }else
      {
        OutputLocus2.df <- tibble(AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], blast_parsed$IDLine_trimmed, IDpercent2, motifs)
        names(OutputLocus2.df) <- headers
        OutputLocus.df <- cbind(OutputLocus.df, OutputLocus2.df)
      }

      # Output Profile
      ProfileEntry <- ""
      if(AlleleInfo[1] == "POS")
      {
        ProfileEntry <- CurrLocus
        if(LocusLkupDNApresent == TRUE)
        {
          ProfileEntry <- paste(CurrLocus, AlleleInfo[3], sep = " ")
          ProfileEntry[AlleleInfo[3] %in% c("WT", "WT/WT", "WT/WT/WT", "WT/WT/WT/WT")] <- ""
        }
        SampleProfile[SampleProfile != "" & ProfileEntry != ""] <- paste(SampleProfile, ProfileEntry, sep = "-")
        SampleProfile[SampleProfile == ""] <- ProfileEntry
      }
      SampleProfile[AlleleInfo[1] == "Sample_Err"] <- "Sample_Err"

      cat(sampleinfo$CurrSampleNo, CurrLocus, AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], 
          AlleleInfo[4], blast_parsed$IDLine_trimmed, IDpercent2, motifs, "\n", sep = "\t")

    } #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End Locus List Loop

    SampleProfile <- ifelse(SampleProfile == "", "Susceptible", SampleProfile)
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
    
    #Progress Bar
    end[m] <- Sys.time()
    time <- round(seconds_to_period(sum(end - init)), 0)
    
    # Estimated remaining time based on the
    # mean time that took to run the previous iterations
    est <- NumSamples * (mean(end[end != 0] - init[init != 0])) - time
    remaining <- round(seconds_to_period(est), 0)
    
    cat(paste(" // Estimated time remaining:", remaining), "\n")
    
  } #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< End Sample Loop

  write.csv(OutputProfile.df, directorylist$outfile, quote = FALSE,  row.names = F)
  
  if(file.exists(dna_nf_file))
  {
    nf_fasta <- read_fasta(dna_nf_file)
    nf_fasta$sq <- as.character(nf_fasta$sq)
    unique_nf <- distinct(nf_fasta, sq, .keep_all = TRUE)
    unique_nf <- unique_nf[,c("name", "sq")]
    unique_nf$name <- paste0(">", unique_nf$name)
  
    unique_nf_file <- paste0(directorylist$output_dir, "output_dna_notfound_distinct.fasta")
    write.table(unique_nf, file = unique_nf_file, row.names = FALSE, col.names = FALSE, 
                quote = FALSE, sep = "")
  }
  
  elapsed <- format(Sys.time() - start_time)
  
  cat("\n\nDone!\n\n\n", "Elapsed time: ", elapsed, 
      "\n\n output_profile_", Org_id, "_",  Test_id, ".csv is ready in output folder", "\n", sep = "")

  return(OutputProfile.df)
}