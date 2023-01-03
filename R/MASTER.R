#' Molecular typing pipeline for WGS assemblies
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


MASTER_pipeline <- function(Org_id, Test_id, SampleNo, LocusID, curr_work_dir){
  #--------------------------------------------------------------------------------------------------------
  #  For troubleshooting and debugging
#   Org_id <- "GONO"               #GAS, PNEUMO or GONO
#   Test_id <- "AMR_ALL"              #AMR, TOXINS, VIRULENCE, NGSTAR use MASTER for 16S id
#   SampleNo <- "list"            
#   LocusID <- "list"             
#   curr_work_dir <- "C:\\WADE\\"
  #--------------------------------------------------------------------------------------------------------

  Variable <- NA
  #Blast_evalue <- "10e-50"            #sets sensitivity of Blast gene match 10e-50 to 10e-150; use 10e-5 for primers

  #------------------------------------------------------------------------------------------------------------
  # get directory structure
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  Directories.df <- read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  Directories.df <- as_tibble(Directories.df)
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir

  #------------------------------------------------------------------------------------------------------------

  cat("\n\n", "SampleNo", "\t", Test_id,"\n", sep = "")
  Test_id_out <- Test_id

  if ((Org_id == "GONO") & (LocusID %in% c("penA", "mtrR", "porB", "ponA", "gyrA", "parC", "rRNA23S")) &
      Test_id %in% c("AMR", "AMR_LW", "AMR_ALL"))
  {
    Test_id <- "NGSTAR"
    Test_id_out <- "AMR"
  }
  if ((Org_id == "GONO") &
      Test_id %in% c("AMR", "AMR_LW", "AMR_ALL"))
  {
    Test_id_out <- "AMR"
  }

  if ((Org_id == "PNEUMO") &
      Test_id %in% c("AMR_ALL", "AMR_2"))
  {
    Test_id <- "AMR"
    Test_id_out <- "AMR"
  }

  switch(Test_id,
         AMR={test_dir <- "Wamr_R\\"},
         AMR_2={test_dir <- "Wamr_R\\"},
         AMR_LW={test_dir <- "Wamr_R\\"},
         AMR_ALL={test_dir <- "Wamr_R\\"},
         TOXINS={test_dir <- "Toxins_R\\"},
         VIRULENCE={test_dir <- "Virulence_R\\"},
         MASTER={test_dir <- "Master_Blaster_R\\"},
         NGMAST={test_dir <- "NGMAST_R\\"},
         NGSTAR={test_dir <- "NGSTAR_R\\"},
         rRNA_16S={test_dir <- "16S_rRNA\\"}
        )
  LkupDir <- paste(system_dir, Org_id, "\\", test_dir, sep = "")

  #------------------------------------------------------------------------------------------------------####
  #  Code starts here ####
  unlink(paste0(local_output_dir,"output_dna.fasta")) #this deletes the file!
  unlink(paste0(local_output_dir, "output_dna_notfound.fasta"))
  unlink(paste0(local_output_dir, "output_aa.fasta"))

  evalueList <- paste(LkupDir, "temp\\blast_evalues.csv", sep = "")
  blast_evalues.df <-  read.csv(evalueList, header = TRUE, sep = ",", stringsAsFactors = FALSE)
  blast_evalues.df <- as_tibble(blast_evalues.df)
  Blast_evalue <- as.character( blast_evalues.df$contig[1])
  evalue_allele <- as.character( blast_evalues.df$allele[1])
  blast_wt_id_threshold <- as.numeric( blast_evalues.df$wt_id[1])

  if(SampleNo == "list")
  {
    SampleList.df <- as_tibble(read.csv(SampList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  }else
  {
    SampleList.df <- tibble(SampleNo, Variable)
  }

  Size.df <- dim(SampleList.df)
  NumSamples <- Size.df[1]

  if(LocusID == "list")
  {
    LocusList <- paste(LkupDir, "temp\\loci.csv", sep = "")
    LocusList.df <- read.csv(LocusList, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    LocusList.df <- as_tibble(LocusList.df)
  }else
  {
    LocusList.df <- tibble(LocusID)
  }

  SizeList <- dim(LocusList.df)
  NumLoci <- SizeList[1]

  LocusMutationsFile <- paste(LkupDir, "temp\\loci_mutations.csv", sep = "")
  if(file.exists(LocusMutationsFile))
  {Mutns <- "Yes"
   LocusMutations.df <- as_tibble(read.csv(LocusMutationsFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  } else {Mutns <- "No"}

  if (Test_id == "MASTER")
  {
  #wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww  Index BLAST lookup files
  for(q in 1L:NumLoci)
  {
    locus <- as.character(LocusList.df[q,1])
    LocusLkupDNA <- paste(LkupDir, "allele_lkup_dna\\", locus, ".fasta", sep = "")
    if(file.exists(LocusLkupDNA))
    {
      BlastFormatCommand <- paste("makeblastdb -in ", LocusLkupDNA, " -dbtype nucl", sep = "")

      try(system(BlastFormatCommand))
    }
  }
  #wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww
  }

  m <- 1L
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  {
    CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
    CurrSampleVar <-as.character(SampleList.df[m, "Variable"])

    CurrSample.df <- filter(SampleList.df, SampleNo == CurrSampleNo)
    SampleProfile <- ""
    cat(m, " of ", NumSamples, "\n")
    p <- 1L
    for(p in 1L:NumLoci)     #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    {
      locus <- as.character(LocusList.df[p,1])

      LocusLkupDNA <- paste(LkupDir, "allele_lkup_dna\\", locus, ".fasta", sep = "")
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

      LocusFile <- paste(LkupDir, "wildgenes\\", locus, ".fasta", sep = "")
      if(!file.exists(LocusFile))
      {
        #stop(paste(locus, " reference file not found!"))
        sample_error.df <- tibble(Loucs_ID = locus, Output = "File Not Found", stringsAsFactors = FALSE)
        return(sample_error.df)
      }

      #  .....................................................................BLAST and parse blastout.txt

      QueryFile <- paste(ContigsDir, CurrSampleNo, ".fasta", sep = "")
      DestFile <- paste(local_temp_dir, "queryfile.fasta", sep = "")
      if (!file.copy(QueryFile, DestFile, overwrite = T))
      {
        AlleleInfo[1] <- "Sample_Err"
        AlleleInfo[2] <- NA
        AlleleInfo[3] <- NA
        AlleleInfo[4] <- NA
        IdLine <- ""
        IDpercent2 <- ""
        IdLine_trimmed <- ""
      } else
      { #make blast database of contig file, then blast locus against it
        FormatCommand <- paste("makeblastdb -in ", DestFile, " -dbtype nucl", sep = "")
        shell(FormatCommand, intern = TRUE)

        BlastCommand <- paste("blastn -query ", LocusFile, " -db ", DestFile, " -out ", local_temp_dir, "blastout.txt -evalue ", Blast_evalue, sep = "")

        shell(BlastCommand, intern = TRUE)

        Blastout <- paste(local_temp_dir, "blastout.txt", sep = "")
        con <- file(Blastout, open="r")
        linn <- readLines(con)
        close(con)

        #check if gene was found in BLAST
        BlastResult <- NA
        for (i in 1:length(linn))
        {
          if (str_detect(linn[i], "No hits found"))
          {
            BlastResult <- "NEG"
            AlleleInfo[1] <- "NEG"
            AlleleInfo[2] <- ""
            AlleleInfo[3] <- ""
            AlleleInfo[4] <- ""
            break()
          } else
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

        if (BlastResult == "POS")
        {
          for (i in 1:length(linn))
          {
            if (str_detect(linn[i], "Score =")) #set counter to parse only first align block of blastout
            {
              TimesThrough <- TimesThrough + 1
            }

            if ((str_detect(linn[i], "Length=" )) & (WTlength == 0)) #get length of WT gene
            {
              WTlengthLine <-  unlist(linn[i])
              WTlength <- as.numeric(substr(WTlengthLine, 8, 50))
            }

            if (TimesThrough == 1)
            {
              if (str_detect(linn[i], "Identities"))
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
          if (LocusID != "list")
          {
            cat("\n\n>", locus , "(Wildtype)\n", WTDNASeqLine_NoDash_str, "\n", sep ="")
          }

          DNASeqLine <- DNAString(DNASeqLine_str)
          DNASeqLine_NoDash_str1 <- str_replace_all(DNASeqLine_str, "-", "")
          DNASeqLine_NoDash_str <- str_replace_all(DNASeqLine_NoDash_str1, "N", "")
          DNASeqLine_NoDash <- DNAString(DNASeqLine_NoDash_str)
          if (LocusID != "list")
          {
            cat("\n\n>", CurrSampleNo, "_", locus , "\n", DNASeqLine_NoDash_str, "\n", sep ="")
          }
          #-------------------------------------------------------------------------------make Protein sequence
          WTAASeqLine <- translate(WTDNASeqLine_NoDash)
          WTAASeqLine_str <- toString(WTAASeqLine)
          if (LocusID != "list")
          {
            cat("\n\n>", locus , "(Wildtype)\n", WTAASeqLine_str, "\n", sep ="")
          }
          AASeqLine <- translate(DNASeqLine_NoDash)
          AASeqLine_str <- toString(AASeqLine)

          AASeqLine_length <- str_length(AASeqLine_str)
          AAseqLine_lastchr <- substring(AASeqLine_str, AASeqLine_length, AASeqLine_length)
          if (AAseqLine_lastchr == "*"){AASeqLine_str <- substring(AASeqLine_str, 1L, AASeqLine_length-1)}

          if (str_detect(AASeqLine_str, "[*]"))
              {
                AlleleInfo[4] <- "Disrupted"
              }

          if (LocusID != "list")
          {
            cat("\n\n>", CurrSampleNo, "_", locus , "\n", AASeqLine_str, "\n", sep ="")
          }

          #-----------------------------------make sure AASeqLine, WTAASeqLine are listed in correct order!!!!!!!!!!!!
          globalAlign_AA <- pairwiseAlignment(AASeqLine, WTAASeqLine, substitutionMatrix = "BLOSUM50", gapOpening = -2,
                                              gapExtension = -8, scoreOnly = FALSE)

          WTAASeqLine_aln <- subject(globalAlign_AA)
          WTAASeqLine_aln_str <- toString(WTAASeqLine_aln)

          AASeqLine_aln <- pattern(globalAlign_AA)
          AASeqLine_aln_str <- toString(AASeqLine_aln)

          DNASeqLine_aln <- ""
          for (j in 1:str_length(WTDNASeqLine_str))
          {
            if (str_sub(WTDNASeqLine_str, j, j) == str_sub(DNASeqLine_str, j, j))
            {
              DNASeqLine_aln <- paste(DNASeqLine_aln, ".", sep = "")
            } else
            {
              DNASeqLine_aln <- paste(DNASeqLine_aln, str_sub(DNASeqLine_str, j, j), sep = "")
            }
          }
          if (LocusID != "list")
          {
            cat("\n\nDNA Alignment:\n", DNASeqLine_aln, "\n", sep = "")
          }
          #------------------------------------------------------------------------------mutations

          motifs <- ""
          aa_mut_name <- ""
          aa_start <- ""
          aa_end <- ""
          aa_wt <- ""
          motif<-""
          #======================================================

          if (Mutns == "Yes")
          {
          LocusMutationsCurr.df <- filter(LocusMutations.df, Locus_id == locus)
          SizeMutns.df <- dim(LocusMutationsCurr.df)
          NumMutations <- SizeMutns.df[1]

          if (NumMutations > 0 )
          {
          for(w in 1L:NumMutations)
          {
          aa_mut_name <- LocusMutationsCurr.df[w,"Name"]
          aa_start <- LocusMutationsCurr.df[w,"Posn_1"]
          aa_end <- LocusMutationsCurr.df[w,"Posn_2"]
          aa_wt <- LocusMutationsCurr.df[w,"WildType"]
          motif<-substr(AASeqLine_aln_str, aa_start, aa_end)
          if (motif == aa_wt)
          {
            motif <- "WT"
            aa_mut_name <- ""
          }
          if (w == 1)
          {motifs <- paste(aa_mut_name, motif, sep = "")}else
          {motifs <- paste(motifs, "/", aa_mut_name, motif, sep = "")}
          } #endfor
          }
          }

          AASeqLine_aln_disp <- ""
          mutations <- ""

          for (k in 1:str_length(WTAASeqLine_aln_str))
          {
            if (str_sub(WTAASeqLine_aln_str, k, k) == str_sub(AASeqLine_aln_str, k, k))
            {
              AASeqLine_aln_disp <- paste(AASeqLine_aln_disp, ".", sep = "")
            } else
            {
              AASeqLine_aln_disp <- paste(AASeqLine_aln_disp, str_sub(AASeqLine_str, k, k), sep = "")
              mutations <- paste(mutations, str_sub(WTAASeqLine_aln_str, k, k), k, str_sub(AASeqLine_str, k, k), " ", sep = "")
            }
          }
          if (LocusID != "list")
          {
            cat("\n\nProtein Alignment:\n", AASeqLine_aln_disp, "\n", sep = "")
            cat("\nMutations: ", mutations, "\n\n", sep = "")
            cat("\nMotifs: ", motifs, "\n\n", sep = "")
          }
          #----------------------------------------------------------------------------------Lookup Alleles DNA
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
              if (str_detect(linn[i], "Identities"))
              {
                IdentLine <-  unlist(linn[i])
                IdentLine <- substr(IdentLine, 15, 50)
                if (str_detect(IdentLine, "100%"))
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

          #-----------------------write a fasta file of all sequences output_dna.fasta and output_aa.fasta

          dna_file <- paste(local_output_dir, "output_dna.fasta", sep = "")
          dna_nf_file <- paste(local_output_dir, "output_dna_notfound.fasta", sep = "")
          aa_file <- paste(local_output_dir, "output_aa.fasta", sep = "")

          sink(dna_file, split=FALSE, append = TRUE)
          cat(">", locus, "_", CurrSampleNo, "_", locus, AlleleInfo[2], "_", AlleleInfo[3], "_", CurrSampleVar, "_", motifs, "\n", DNASeqLine_NoDash_str, "\n", sep ="")
          sink()

          if (AlleleInfo[2] == "NF")
          {
            sink(dna_nf_file, split=FALSE, append = TRUE)
            cat(">", locus, "_", motifs, "_", CurrSampleNo, "_", CurrSampleVar, "_", "\n", DNASeqLine_NoDash_str, "\n", sep ="")
            sink()
          }

          sink(aa_file, split=FALSE, append = TRUE)
          cat(">", locus, "_", CurrSampleNo, "_", CurrSampleVar, "\n", AASeqLine_str, "\n", sep ="")
          sink()

          IDpercent2 <- paste(IDlength, "/", WTlength, " (", IDpercent, ")", sep = "") #made a blast ID line based on full WT gene

          if ((locus == "porB") & (IDcoverage < 90))
          {
            AlleleInfo[4] <- "porB1a"
          }else if ((locus == "porB") & (IDcoverage >= 90))
          {
            AlleleInfo[4] <- "porB1b"
          }

          if (IDcoverage <= blast_wt_id_threshold)
          {
            AlleleInfo[1] <- "NEG"  #blast result
            #AlleleInfo[2] <- ""  #allele number
            #AlleleInfo[3] <- ""  #mutations
            #AlleleInfo[4] <- ""  #comments
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
      if (AlleleInfo[1] == "POS")
      {
        ProfileEntry<-locus

        if (LocusLkupDNApresent)
        {
          if (AlleleInfo[3] == "" | is.na(AlleleInfo[3]) | AlleleInfo[3] == "???")
          {
            ProfileEntry<- paste(locus, AlleleInfo[3], sep = " ")
          }else
          if (AlleleInfo[3]=="WT" | AlleleInfo[3]=="WT/WT" | AlleleInfo[3]=="WT/WT/WT" | AlleleInfo[3]=="WT/WT/WT/WT")
                {ProfileEntry<-""}
          else
                {ProfileEntry <- paste(locus, AlleleInfo[3], sep = " ")}
        }

        if (SampleProfile == "")
        {
          SampleProfile <- ProfileEntry
        }else
        {
          if (ProfileEntry!=""){ProfileEntry<-paste("-", ProfileEntry, sep = "")}
          SampleProfile <- paste(SampleProfile, ProfileEntry, sep = "")
        }
      }

      if (AlleleInfo[1] == "Sample Err")
      {
        SampleProfile <- "Sample Err"
      }

      cat(CurrSampleNo, locus, AlleleInfo[1], AlleleInfo[2], AlleleInfo[3], AlleleInfo[4], IdLine_trimmed, IDpercent2, motifs, "\n", sep = "\t")

    } #end of locus list loop xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

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
  } #close brack for sample list loop<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  if (Test_id_out == "AMR")  { Test_id <- "AMR"}

  outfile <- paste(local_output_dir, "output_profile_", Org_id, "_", Test_id, ".csv", sep = "")

  write.csv(OutputProfile.df, outfile, quote = FALSE,  row.names = F)

  cat("\n\nDone! output_profile_", Org_id, "_MASTER.csv is ready in output folder", "\n", sep = "")

  return(OutputProfile.df)

}