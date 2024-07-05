#' MLST pipeline for WGS assemblies
#' July 3, 2024, Walter Demczuk & Shelley Peterson
#'
#' Takes Organism, Sample Number, Locus, and a Variable at queries a contig.fasta file
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param SampleNo Sample number associated with contig.fasta file
#' @param locus The locus to query, or enter list to use a list of alleles
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details MLST pipeline:
#' Takes Organism, Sample Number, Locus, and a Variable at queries a contig.fasta file vs downloaded pubMLST data.
#' @return A table frame containing the results of the query
#' @export
#'

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#  Org_id <- "GONO"                  #GAS, GBS, PNEUMO or GONO
#  SampleNo <- "list"
#  LocusID <- "list"
#  Test_id <- "NGSTAR"               #NGSTAR, NGMAST, MLST, NGSTAR_AMR
#  curr_work_dir <- "C:\\WADE\\"
 # Blast_evalue <- "10e-50"          #sets sensitivity of Blast gene match 10e-50 to 10e-150; use 10e-5 for primers
#-------------------------------------------------------------------------------

MLST_pipeline <- function(Org_id, Test_id, SampleNo, LocusID, curr_work_dir) {
  
  start_time <- Sys.time()

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  if(Test_id == "NGSTAR_AMR")
  {
    Test_id <- "NGSTAR"
    Test <- "NGSTAR_AMR"
  } else
  {
    Test <- "other"
  }
  directorylist <- getdirectory(curr_work_dir, Org_id, Test_id)
  reflist <- refdirectory(directorylist, Org_id, Test_id)
  unlink(paste0(directorylist$temp_dir, "blastout.txt")) #this deletes the file!
  #-----------------------------------------------------------------------------
  # Blast Evalues
  blast_evalues.df <- as_tibble(read.csv(paste0(reflist$Ref_Dir, "blast_evalues.csv"),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Blast_evalue <- as.character(blast_evalues.df$contig[1])
  #-----------------------------------------------------------------------------
  
  ########################### Setup locus list table ###########################
  LocusList.df <- as_tibble(read.csv(reflist$Loci_List, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  if(LocusID != "list")
  {
    LocusList.df <- filter(LocusList.df, Locus_id == LocusID)
  }
  locuslist <- LocusList.df$Locus_id
  Output.df <- tibble(Locus = c("SampleNo", locuslist))
  NumLoci <- dim(LocusList.df)[1]

  Variable <- NA

  ########################### Setup sample list table ##########################
  if(SampleNo == "list")
  {
    SampleList.df <- as_tibble(read.csv(directorylist$SampleList, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  }else
  {
    SampleList.df <- tibble(SampleNo, Variable)
  }

  NumSamples <- dim(SampleList.df)[1]

  ######################## Load sequence type profiles #########################
  if(Test == "other")
  {
    profiles.df <-as_tibble(read.csv(reflist$Profiles, header = TRUE, sep = "\t", stringsAsFactors = FALSE))
  }
  
  if(Test_id == "NGMAST")
  {
    names(profiles.df) <- c("ST", "porB", "tbpB")
    df.bad_alleles <- as_tibble(read.csv(paste0(reflist$Ref_Dir, "bad_alleles.csv"),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE))
  }

  #Progress Bar
  n_iter <- NumSamples
  init <- numeric(n_iter)
  end <- numeric(n_iter)
  
  ##############################################################################
  # BLAST loci + determine allele numbers
  ##############################################################################
  m <- 1L
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Start Sample Loop
  {
    init[m] <- Sys.time()
    
    cat(m, " of ", NumSamples, "\n")
    
    #actual Code
    sampleinfo <- getsampleinfo(SampleList.df, directorylist, reflist, m)
    if(Test == "NGSTAR_AMR")
    {
      Sample_Output <- as_tibble(sapply(locuslist, mutationblast, sampleinfo = sampleinfo, 
                                        directorylist = directorylist,
                                        reflist = reflist, Test_id = Test_id, 
                                        df.bad_alleles = df.bad_alleles))
    } else
    {
      Sample_Output <- as_tibble(sapply(LocusList.df$Locus_id, locusblast, 
                                        LocusList.df = LocusList.df,
                                        sampleinfo = sampleinfo, 
                                        reflist = reflist, Test_id = Test_id, 
                                        Blast_evalue = Blast_evalue, 
                                        df.bad_alleles = df.bad_alleles, 
                                        directorylist = directorylist))
    }
    names(Sample_Output) <- sampleinfo$CurrSampleNo
    Sample_Output <- rbind(sampleinfo$CurrSampleNo, Sample_Output)

    print(list(Sample_Output))

    Output.df <- bind_cols(Output.df, Sample_Output)

    #Progress Bar
    end[m] <- Sys.time()
    time <- round(seconds_to_period(sum(end - init)), 0)

    # Estimated remaining time based on the
    # mean time that took to run the previous iterations
    est <- NumSamples * (mean(end[end != 0] - init[init != 0])) - time
    remaining <- round(seconds_to_period(est), 0)

    cat(paste(" // Estimated time remaining:", remaining), "\n")
  }  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< End Sample Loop

  # Transpose columns so they're in the correct format
  Output.df <- data.table::transpose(Output.df)
  names(Output.df) <- c("SampleNo", locuslist)
  Output.df <- filter(Output.df, !SampleNo == "SampleNo")

  ##############################################################################
  # Lookup ST profiles
  ##############################################################################

  if(Test == "NGSTAR_AMR" | LocusID != "list")
  {
    Output <- Output.df
  }else
  {
    if(Test_id == "NGSTAR")
    {
      profiles.df <- profiles.df %>% dplyr::rename("rRNA23S" = "X23S")
    }
    profiles.df <- profiles.df %>% dplyr::mutate(across(everything(), as.character))
    Output <- left_join(Output.df, profiles.df)
    Output$clonal_complex <- NA
    Output <- select(Output, -clonal_complex)
  }

  ##############################################################################
  # Export Results to csv file
  ##############################################################################

  if(Test_id == "NGMAST")
  {
    Output$porB <- paste0("porB-", Output$porB)
    Output$tbpB <- paste0("tbpB-", Output$tbpB)
    Output$ST <- paste0("ST-", Output$ST)
  }

  if(Test == "NGSTAR_AMR" | LocusID != "list")
  {
    write.csv(Output, paste0(directorylist$output_dir, "output_profile_mut_GONO_NGSTAR.csv"), quote = FALSE,  row.names = FALSE)
  } else
  {
    if(Test_id == "NGMAST")
    {
      Output_good <- dplyr::filter(Output, ST != "ST-NA")
      Output_bad <- dplyr::filter(Output, ST == "ST-NA")
    } else
    {
      Output_good <- filter(Output, !is.na(ST))
      Output_bad <- filter(Output, is.na(ST)) 
    }
    write.csv(Output_good, paste0(directorylist$output_dir, "LabWareUpload_", Org_id, "_", Test_id, "_good.csv"), quote = FALSE, row.names = FALSE)
    write.csv(Output_bad, paste0(directorylist$output_dir, "LabWareUpload_", Org_id, "_", Test_id, "_bad.csv"), quote = FALSE, row.names = FALSE)
    write.csv(Output, paste0(directorylist$output_dir, "LabWareUpload_", Org_id, "_", Test_id, ".csv"), quote = FALSE, row.names = FALSE)
  }
  write.csv(Output, paste0(directorylist$output_dir, "output_profile_", Org_id, "_", Test_id, ".csv"), quote = FALSE,  row.names = FALSE)

  elapsed <- format(Sys.time() - start_time)
  
  cat("\n\nDone!\n\n\n", "Elapsed time: ", elapsed, "\n\n output_profile_", Org_id,
      "_", Test_id, ".csv is ready in output folder\n", sep = "")

  return(Output) 
  
} 
