#' Serotyping pipeline for WGS assemblies
#' February 5 2024, Walter Demczuk & Shelley Peterson
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
#' 1)result = POS/NEG (presence/absence)
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

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#Org_id <- "PNEUMO"               #GBS or PNEUMO
#Test_id <- "SERO"
#SampleNo <- "list"               #Sample No or "list"
#LocusID <- "list"
#curr_work_dir <- "C:\\WADE\\"
#Blast_evalue <- "10e-50"         #sets sensitivity of Blast gene match 10e-50 to 10e-150; use 10e-5 for primers
#-------------------------------------------------------------------------------

SEROTYPE_pipeline <- function(Org_id, SampleNo, LocusID, Test_id, curr_work_dir){

  start_time <- Sys.time()
  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  directorylist <- getdirectory(curr_work_dir, Org_id, Test_id)
  reflist <- refdirectory(directorylist, Org_id, Test_id)
  unlink(paste0(directorylist$output_dir, "LabWareUpload_", Org_id, "_SEROTYPE.csv"))
  #-----------------------------------------------------------------------------
  # Load CPS types info
  LocusLkupDNA_CPS <- paste0(reflist$Lkup_Dir, "reference_CPS.fasta")
  CPSlist <- paste0(directorylist$temp_dir, "CPS_types_for_blast.csv")
  CPS_types <- read.csv(paste0(reflist$Ref_Dir, "CPS_types_for_blast.csv"), 
                               header = TRUE, sep = ",", stringsAsFactors = FALSE)
  CPS_types <- CPS_types$Serogroup
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
    SampleList.df <- tibble(SampleNo, Variable = NA)
  }
  
  NumSamples <- dim(SampleList.df)[1]

  #Progress Bar
  n_iter <- NumSamples
  init <- numeric(n_iter)
  end <- numeric(n_iter)
  
  m <- 1L
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Start Sample Loop
  {
    init[m] <- Sys.time()
    sampleinfo <- getsampleinfo(SampleList.df, directorylist, reflist, m)
    
    SampleProfile <- ""
    serogroup <- ""
    serotype <- ""
    OutputLocus.df <- tibble(SampleNo = as.character(sampleinfo$CurrSampleNo))  #all results

    cat("\n\n", sampleinfo$CurrSampleNo, ": ", m, " of ", NumSamples, "\n", sep = "")

    ############################################################################
    # Get serogroup CPS locus for sample first
    ############################################################################
    if (sampleinfo$Allele[1] == "Sample_Err")
    {
      serogroup <- "Sample_Err"
      serotype <- "Sample_Err"
      SampleProfile <- "Sample_Err"
    }else #..................................................................... Not Sample Error
    {
      #makeblastdb from contig file for later when blasting wildgenes vs. contig to extract genes
      BlastFormatCommand <- paste0("makeblastdb -in ", reflist$DestFile, " -dbtype nucl")
      try(system(BlastFormatCommand))
      
      #use contig file to blast against the CPS data to see what the CPS group is
      Blast_Out_File <- blastquery(directorylist, reflist, "reference_CPS", Blast_evalue)
      info = file.info(Blast_Out_File)
      if(info$size == 0)  #no blast result from CPS locus lookup i.e. not pneumo, bad sequencing
      {
        serotype <- "unknown"
        serogroup <- "unknown"
      }else
      {
        df.blastout <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))      
        names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
        df.blastout <- arrange(df.blastout, dplyr::desc(bit))
        
        serotype1 <- df.blastout$Allele[1]
        SerotypeParts <- unlist(strsplit(serotype1, "_"))
        SerotypeParts <- SerotypeParts[1]
        serotype <- ifelse(substr(SerotypeParts, 1, 1) == "0", sub("^.", "", SerotypeParts), SerotypeParts)
        serotype[serotype %in% c("15B", "15C")] <- "15BC"

        if (Org_id == "PNEUMO")
        {
          serogroup <- as.character(substr(SerotypeParts, 1, 2))
          if(serogroup == "40") {serogroup <- "07"}
          if(serogroup == "44" || serogroup == "46"){serogroup <- "12"}
          if(serogroup == "38") {serogroup <- "25"}
          if(serogroup == "37") {serogroup <- "33"}
          if(serogroup == "42") {serogroup <- "35"}
        }
        if (Org_id == "GBS")
        {
          serogroup <- "GBS"
        }
        outfile_sero <- paste0(directorylist$output_dir, "output_blastout_serogroup.csv")
        write.csv(df.blastout, outfile_sero, quote = FALSE, row.names = F )
      }
    }#.......................................................................... close Not Sample Error loop

    cat("CPS Blast Serogroup ", serogroup, " - Serotype ", serotype, " !\n")

    if (serogroup == "unknown" | serogroup == "Sample_Err")  # no CPS match or no contig
    {
      serotype <- paste("NT[", serogroup, "]", sep = "")
      SampleProfile <- serogroup
    }else # //////////////////////////////////////////////////////////////////// Serogroup Found 
    {
      ##########################################################################
      # setup loci and get locuslist screen for serotypes requiring no SNP analysis
      ##########################################################################
      if(serotype %in% CPS_types)
      {
        SampleProfile <- "CPS operon type match"
        locus <- ""
        head(df.blastout, n = 5L)
      }else #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ Not in CPS_types list
      {
        LocusList <- paste0(reflist$Ref_Dir, serogroup, "_loci.csv")
        if(file.exists(LocusList))
        {
          LocusList.df <- read.csv(LocusList, header = TRUE, sep = ",", stringsAsFactors = FALSE)
          SerotypeLookup.df <- read.csv(paste0(reflist$Ref_Dir, serogroup, "_loci_lookup.csv"), 
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE,
                                        check.names = FALSE)
          if(LocusID != "list")
          {
            LocusList.df <- filter(LocusList.df, Locus_id == LocusID)
          }
          NumLoci <- dim(LocusList.df)[1]
        }else 
        {
          NumLoci <- 0
        }

        LocusMutationsFile <- paste0(reflist$Ref_Dir, serogroup, "_loci_mutations.csv")
        if(file.exists(LocusMutationsFile))
        {
          Mutns <- "Yes"
          LocusMutations.df <- read.csv(LocusMutationsFile, header = TRUE, sep = ",", stringsAsFactors = FALSE)
        } else 
        {
          Mutns <- "No"
        }

        OutputLocus.df <- tibble(SampleNo = as.character(sampleinfo$CurrSampleNo))  #all results
        OutputLocusProfile.df <- OutputLocus.df #for serotype
    
        p <- 1L
        for(p in 1L:NumLoci) #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Locus Loop
        {
          ProfileSeparator <- ifelse(SampleProfile == "", "", ";")
          CurrLocus <- as.character(LocusList.df$Locus_id[p])
          locus <- gsub("GBS_", "", CurrLocus) #makes profile smaller for GBS LabWare output
          locus_result_type <- as.character(LocusList.df$Result_type[p])
        
          LocusLkupDNA <- paste0(reflist$Lkup_Dir, locus, ".fasta")
          LocusLkupDNApresent <- ifelse(file.exists(LocusLkupDNA), TRUE, FALSE)

          resultcol <- paste0(locus, "_result")
          allelecol <- paste0(locus, "_allele")
          pseudocol <- paste0(locus, "_pseudo")
          mutationscol <- paste0(locus, "_mutations")
          OutputLocus.df[1, resultcol] <- NA
          OutputLocus.df[1, allelecol] <- NA
          OutputLocus.df[1, mutationscol] <- NA
          OutputLocus.df[1, pseudocol] <- NA

          LocusFile <- paste0(reflist$WT_Dir, locus, ".fasta")
          if(!file.exists(LocusFile))
          {
            sample_error.df <- tibble(Locus_ID = LocusID, Output = "Reference wildtype file not found.")
            OutputLocus.df <- OutputLocus.df %>% replace(is.na(.), "Locus_Err")
            return(sample_error.df)
          }else #*************************************************************** WT gene found
          {
            BlastResult <- blasthits(reflist, directorylist, LocusFile, Blast_evalue)
            OutputLocus.df[1,2] <- BlastResult  #result
            blastoutput <- readLines(paste0(directorylist$temp_dir, "blastout.txt"))

            if(locus_result_type == "result")
            {
              SampleProfile <- paste0(SampleProfile, ProfileSeparator, locus, "[", BlastResult, "]")
              cat(paste0(locus, "_result"), BlastResult, "\n", sep = "\t\t")
              
              if(LocusID == "list")
              {
                SerotypeLookup.df <- SerotypeLookup.df %>% filter(.[[resultcol]] == BlastResult|
                                                                   is.na(.[[resultcol]]))
              }
            }
            #------------------------------------------------------------------- 
            #add missing genes to profile for mutations and pseudogenes
            if((locus_result_type == "pseudo") & (BlastResult == "NEG"))
            {
              SampleProfile <- paste(SampleProfile, ProfileSeparator, locus, "[no pseudogene]", sep="")
              cat(pseudocol, "no gene present", "\n", sep = "\t\t")
            }
            if((locus_result_type == "mutations") & (BlastResult == "NEG"))
            {
              SampleProfile <- paste(SampleProfile, ProfileSeparator, locus, "[no mutn gene]", sep="")
              cat(mutationscol, "no gene present", "\n", sep = "\t\t")
            }
            if((locus_result_type == "allele") & (BlastResult == "NEG"))
            {
              SampleProfile <- paste(SampleProfile, ProfileSeparator, locus, "[no allele gene]", sep="")
              cat(allelecol, "no gene present", "\n", sep = "\t\t")
            }
            #-------------------------------------------------------------------

            DNASeqLine_str <- ""
            WTDNASeqLine_str <- ""
            blastdetails <- ""
            if(BlastResult == "POS") #========================================== Positive BLAST Result Loop
            {
              # Parse out to get WT and query DNA sequences, match % etc.
              blast_parsed <- parseblast(blastoutput, LocusID)
              
              # Print both WT and query DNA sequences
              if(LocusID != "list")
              {
                cat("\n\n>", CurrLocus, "(Wildtype)\n", blast_parsed$WTDNASeqLine_NoDash_str, "\n", sep ="")
                cat("\n\n>", sampleinfo$CurrSampleNo, "_", CurrLocus , "\n", blast_parsed$DNASeqLine_NoDash_str, "\n", sep ="")
              }
              
              WTDNASeqLine <- DNAString(blast_parsed$WTDNASeqLine_str)
              WTDNASeqLine_NoDash <- DNAString(blast_parsed$WTDNASeqLine_NoDash_str)
              DNASeqLine <- DNAString(blast_parsed$DNASeqLine_str)
              DNASeqLine_NoDash <- DNAString(blast_parsed$DNASeqLine_NoDash_str)

              ##################################################################
              # Make Protein Sequence + Check if gene is intact
              ##################################################################
              blastAA <- AAconvert(WTDNASeqLine_NoDash, DNASeqLine_NoDash)

              if(locus_result_type == "pseudo")
              {
                #update locus output
                WTend <- head(unlist(gregexpr('[*]', blastAA$WTAASeqLine_str[1])), n = 1)
                WTend <- ifelse(is.na(WTend), "0", WTend)
                queryend <- head(unlist(gregexpr('[*]', blastAA$AASeqLine_str[1])), n = 1)
                queryend <- ifelse(is.na(queryend), 0, queryend)
                diff <- ifelse(queryend == -1, 0, (WTend - queryend))
                lengthdiff <- blast_parsed$WTlength - blast_parsed$IDlength

                nonsense_mutation <- ifelse(diff > 6, "Disrupted", 
                                            ifelse((blast_parsed$WTlength - blast_parsed$IDlength > 25),
                                            "Disrupted", "Intact"))
                
                OutputLocus.df[pseudocol] <- nonsense_mutation
                OutputLocusProfile.df[pseudocol] <- nonsense_mutation
                SampleProfile <- paste0(SampleProfile, ProfileSeparator, locus, "[", nonsense_mutation, "]")
                cat(pseudocol, nonsense_mutation, "\n", sep = "\t\t")
                if(LocusID == "list")
                {
                  SerotypeLookup.df <- SerotypeLookup.df %>% filter(.[[pseudocol]] == nonsense_mutation|
                                                                    is.na(.[[pseudocol]]))
                }
              }

              if(LocusID != "list")
              {
                cat("\n\n>", locus , "(Wildtype)\n", blastAA$WTAASeqLine_str, "\n", sep ="")
                cat("\n\n>", sampleinfo$CurrSampleNo, "_", locus , "\n", blastAA$AASeqLine_str, "\n", sep ="")
              }
            
              ##################################################################
              # Mutations
              ##################################################################
              motifs <- ""
              mutationspresent <- ""
            
              if(Mutns == "Yes")
              {
                motifs <- getmotifs(LocusMutations.df, Test_id, locus, LocusID, blastAA$AASeqLine_str)
              }
              
              #motifs <- ifelse(motifs == "?", "NA", motifs)
              OutputLocus.df[mutationscol] <- motifs
              if(locus_result_type == "mutations")
              {
                OutputLocusProfile.df[mutationscol] <- motifs  #mutations
                SampleProfile <- paste(SampleProfile, ProfileSeparator, locus, "[", motifs, "]", sep="")
                cat(mutationscol, motifs, "\n", sep = "\t\t")
                if(LocusID == "list")
                {
                  SerotypeLookup.df <- SerotypeLookup.df %>% filter(.[[mutationscol]] == motifs|
                                                                    is.na(.[[mutationscol]]))
                }
              }
            
              ##################################################################
              # Lookup DNA alleles
              ##################################################################
              # write DNASeqLine_NoDash_str to DestFile,
              # BLAST vs. lookup table,
              # parse out the allele numbers
              sink(reflist$DestFile, split=FALSE, append = FALSE)
              cat(">", sampleinfo$CurrSampleNo, "_", locus , "\n", blast_parsed$DNASeqLine_NoDash_str, sep ="")
              sink()

              if(locus_result_type == "allele")
              {
                #BLAST lookup table
                Allele <- locusblast(locus, LocusList.df, sampleinfo, reflist, Test_id, Blast_evalue, 
                                     directorylist = directorylist)
                #Allele <- ifelse(Allele == "?", NA, Allele)
                OutputLocus.df[allelecol] <- Allele
                OutputLocusProfile.df[allelecol] <- Allele  #allele
                SampleProfile <- paste0(SampleProfile, ProfileSeparator, locus, "[", Allele, "]")
                cat(allelecol, Allele,"\n", sep = "\t\t")
                if(LocusID == "list") 
                {
                  SerotypeLookup.df <- SerotypeLookup.df %>% filter(.[[allelecol]] == Allele|
                                                                    is.na(.[[allelecol]]))
                }
              }
            }#================================================================== End BLAST positive
          }#******************************************************************** End WT gene found 
        }#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End Locus Loop

        # Interpret molecular characterization to get serotype
        if(LocusID == "list")
        {
          NumResults <- dim(SerotypeLookup.df)[1]
          if(NumResults == 0)
          {
            serotype <- paste0("NT[", serotype, "]")
          }else
          {
            for(w in 1L:NumResults)
            {
              if(w == 1)
              {
                serotype <-  SerotypeLookup.df$Serotype[w]
              }else
              {
                serotype <- paste0(serotype, "/", SerotypeLookup.df$Serotype[w])
              }
            }
          }
        }
      } #^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ End Not in CPS Types List
    } #///////////////////////////////////////////////////////////////////////// End Serogroup Found

    if (m==1L)
    {
      SampleOutputProfile.df <- tibble(SampleNo = sampleinfo$CurrSampleNo, Serotype = serotype, Profile = SampleProfile)
    }else
    {
      SampleOutputProfile2.df <- tibble(SampleNo = sampleinfo$CurrSampleNo, Serotype = serotype, Profile = SampleProfile )
      SampleOutputProfile.df<- bind_rows(SampleOutputProfile.df, SampleOutputProfile2.df)
    }

    cat("\nProfile:  ", SampleProfile, sep = "")
    cat("\nSerotype: ", serotype, sep = "")
    
    #Progress Bar
    end[m] <- Sys.time()
    time <- round(seconds_to_period(sum(end - init)), 0)
    
    # Estimated remaining time based on the
    # mean time that took to run the previous iterations
    est <- NumSamples * (mean(end[end != 0] - init[init != 0])) - time
    remaining <- round(seconds_to_period(est), 0)
    
    cat(paste(" // Estimated time remaining:", remaining), "\n")
  } #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< End Sample Loop

  outfile <- paste0(directorylist$output_dir, "LabWareUpload_", Org_id, "_SEROTYPE.csv")
  if (LocusID == "list")
  {
    write.csv(SampleOutputProfile.df, outfile, quote = FALSE,  row.names = F)
  }

  if (NumSamples == 1L)
  {
    outfile_sero2 <- paste(directorylist$output_dir, "output_profile_SEROTYPE.csv", sep = "")
    write.csv(OutputLocus.df, outfile_sero2, quote = FALSE, row.names = F )
  }

  dna_nf_file <- paste0(directorylist$output_dir, "output_dna_notfound.fasta")
  if(file.exists(dna_nf_file))
  {
    nf_fasta <- read_fasta(dna_nf_file)
    nf_fasta$sq <- as.character(nf_fasta$sq)
    unique_nf <- distinct(nf_fasta, sq, .keep_all = TRUE)
    unique_nf <- unique_nf[,c("name", "sq")]
    unique_nf$name <- paste0(">", unique_nf$name)
    
    unique_nf_file <- paste0(directorylist$output_dir, "output_dna_notfound_distinct.fasta")
    write.table(unique_nf, file = unique_nf_file, row.names = FALSE, col.names = FALSE, 
                quote = FALSE, sep = "\n")
  }
  
  elapsed <- format(Sys.time() - start_time)
  
  cat("\n\nDone!\n\n\n", "Elapsed time: ", elapsed, 
      "\n\n Serotyping Results are ready in output folder\n", sep = "")

  return(SampleOutputProfile.df)
}
