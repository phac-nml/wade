#' Common helper functions for WADE analysis
#' February 5, 2024, Shelley Peterson
#' 
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param Test_id AMR, TOXINS, VIRULENCE, NGSTAR, use MASTER for 16S id
#' @param SampleNo Sample number associated with contig.fasta file
#' @param LocusID The locus to query, or enter "list" to use a list of alleles
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @param directorylist the list generated from the getdirectory function
#' @param SampleList.df the list of samples being used 
#' @param reflist the list generated from the refdirectory function
#' @param NumLoci the number of loci being analyzed
#' @param LocusList.df the list of loci being analyzed
#' @param Blast_evalue the evalue parameter for the Blast function
#' @param CurrLocus the current locus being analyzed
#' @param LocusFile the location for the reference locus file
#' @param blastoutput the BLAST output file to be parsed 
#' @details These are helper functions to streamline the other WADE functions
#' @return a dataframe or list to be fed into other WADE functions
#' 
#' @export
# get directory structure
getdirectory <- function(curr_work_dir, Org_id, Test_id){ 
  dir_file <- paste0(curr_work_dir, "DirectoryLocations.csv")
  Directories.df <- as_tibble(read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)) %>%
    filter(OrgID == Org_id)
  SampleList <- paste0(Directories.df$LocalDir, "list.csv")
  output_dir <- paste0(Directories.df$LocalDir, "Output/")
  temp_dir <- paste0(Directories.df$LocalDir, "temp/")
  ref_dir <- paste0(Directories.df$LocalDir, "reference/")
  system_dir <- Directories.df$SystemDir
  contigs_dir <- Directories.df$ContigsDir
  vcf_dir <- Directories.df$VCFDir
  outfile <- paste0(output_dir, "output_profile_", Org_id, "_", Test_id, ".csv")

  return(list("SampleList" = SampleList,
              "output_dir" = output_dir,
              "temp_dir" = temp_dir,
              "ref_dir" = ref_dir,
              "system_dir" = system_dir,
              "contigs_dir" = contigs_dir,
              "vcf_dir" = vcf_dir,
              "outfile" = outfile))
}

#' @export
# get reference files structure
refdirectory <- function(directorylist, Org_id, Test_id){
  Lkup_Dir <- paste0(directorylist$system_dir, Org_id, "/", Test_id, "/allele_lkup_dna/")
  Ref_Dir <- paste0(directorylist$system_dir, Org_id, "/", Test_id, "/reference/")
  WT_Dir <- paste0(directorylist$system_dir, Org_id, "/", Test_id, "/wildgenes/")
  Loci_List <- paste0(Ref_Dir, "loci.csv")
  Profiles <- paste0(Ref_Dir, "profiles.txt")
  DestFile <- paste0(directorylist$temp_dir, "queryfile.fasta")
  
  return(list("Lkup_Dir" = Lkup_Dir,
              "Ref_Dir" = Ref_Dir,
              "WT_Dir" = WT_Dir,
              "Loci_List" = Loci_List,
              "Profiles" = Profiles,
              "DestFile" = DestFile,
              "Variable" = NA))
}

#' @export
# remove old DNA output files
removefiles <- function(curr_work_dir){
  unlink(paste0(curr_work_dir, "Output/output_dna.fasta")) #this deletes the file!
  unlink(paste0(curr_work_dir, "Output/output_dna_notfound.fasta"))
  unlink(paste0(curr_work_dir, "Output/output_dna_notfound_distinct.fasta"))
  unlink(paste0(curr_work_dir, "Output/output_aa.fasta"))
}

#' @export
# summarize sample info
getsampleinfo <- function(SampleList.df, directorylist, reflist, m){
  Allele <- ""
  AlleleNo <- ""

  CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
  CurrSampleVar <- as.character(SampleList.df[m, "Variable"])
  QueryFile <- paste0(directorylist$contigs_dir, CurrSampleNo, ".fasta")
  
  if(!file.copy(QueryFile, reflist$DestFile, overwrite = T))
  {
    Allele <- "Sample_Err"
    AlleleNo <- "Sample_Err"
  }
  return(list("Allele" = Allele,
              "AlleleNo" = AlleleNo,
              "CurrSampleNo" = CurrSampleNo,
              "CurrSampleVar" = CurrSampleVar,
              "QueryFile" = QueryFile))
}

#' @export
# makeblastDB  
quickBlastIndex <- function(NumLoci, LocusList.df, reflist){
  q <- 1
  for(q in 1L:NumLoci)
  {
    CurrLocus <- as.character(LocusList.df[q,1])
    LocusLkupDNA <- paste0(reflist$Lkup_Dir, CurrLocus, ".fasta")
    if(file.exists(LocusLkupDNA))
    {
      BlastFormatCommand <- paste0("makeblastdb -in ", LocusLkupDNA, " -dbtype nucl")
      try(system(BlastFormatCommand))
    }
  }
}

#' @export
# run BLAST for Allele Lookups
blastquery <- function(directorylist, reflist, CurrLocus, Blast_evalue){
  LocusLkupDNA <- paste0(reflist$Lkup_Dir, CurrLocus, ".fasta")
  Blast_Out_File <- paste0(directorylist$temp_dir, "blastout.txt")
  BlastCommand <- paste0("blastn -query ", reflist$DestFile,
                        " -db ", LocusLkupDNA,
                        " -out ", Blast_Out_File,
                        " -num_alignments 100 -evalue ", Blast_evalue, " -outfmt 6")
  system(BlastCommand)
  return(Blast_Out_File)
}

#' @export
#run BLAST for gene detection
blasthits <- function(reflist, directorylist, LocusFile, Blast_evalue){
  BlastCommand <- paste0("blastn -query ", LocusFile, 
                         " -db ", reflist$DestFile, 
                         " -out ", directorylist$temp_dir, "blastout.txt ", 
                         "-evalue ", Blast_evalue)
  system(BlastCommand)
  blastoutput <- readLines(paste0(directorylist$temp_dir, "blastout.txt"))

  #check if gene was found in BLAST
  nohits <- any(grepl("No hits found", blastoutput))
  BlastResult <- ifelse(nohits == TRUE, "NEG", "POS")
  return(BlastResult)
}

#' @export
#parse blast to get info on best alignment
parseblast <- function(blastoutput, LocusID){
  WTlength <- as.numeric(sub("Length=", "", grep("Length=", blastoutput, value = TRUE)[1]))
  IDLine <- grep("Identities", blastoutput, value = TRUE)[1]
  IDLine_trimmed <- sub(" Identities = (.+),.+", "\\1", IDLine)
  IDlength <- as.numeric(sub(" Identities = (\\d+).+", "\\1", IDLine))
  IDpercent <- paste0(round((IDlength / WTlength * 100),1),"%")
  IDcoverage <- ((IDlength / WTlength) * 100)
  
  #trim from blast file so it only shows first match
  contigsline <- grep(">contig", blastoutput)
  if(length(contigsline) > 1)
  {
    blastoutput <- blastoutput[1:contigsline[2]]   
  }

  #get sequence from first match
  QueryLine <- grep("Query ", blastoutput, value = TRUE)
  WTDNASeqLine_str <- gsub("[^AaCcTtGgN-]", "", paste(QueryLine, collapse = ""))
  WTDNASeqLine_NoDash_str <- gsub("[-N]", "", WTDNASeqLine_str)
  SbjctLine <- grep("Sbjct ", blastoutput, value = TRUE)
  SbjctLine <- gsub("Sbjct","", SbjctLine)
  DNASeqLine_str <- gsub("[^AaCcTtGgN-]", "", paste(SbjctLine, collapse = ""))       
  DNASeqLine_NoDash_str <- gsub("[-N]","",DNASeqLine_str)
  
  result <- tibble(WTlength = WTlength, IDLine_trimmed = IDLine_trimmed, IDlength = IDlength,
                   IDpercent = IDpercent, IDcoverage = IDcoverage, WTDNASeqLine_str = WTDNASeqLine_str,
                   WTDNASeqLine_NoDash_str = WTDNASeqLine_NoDash_str, DNASeqLine_str = DNASeqLine_str,
                   DNASeqLine_NoDash_str = DNASeqLine_NoDash_str)
  return(result)
}

#' @export
#convert from DNA strings to AA strings
AAconvert <- function(WTDNASeqLine_NoDash, DNASeqLine_NoDash){
  WTAASeqLine <- suppressWarnings(Biostrings::translate(WTDNASeqLine_NoDash))
  WTAASeqLine_str <- toString(WTAASeqLine)
  AASeqLine <- suppressWarnings(Biostrings::translate(DNASeqLine_NoDash))
  AASeqLine_str <- toString(AASeqLine)
  AASeqLine_length <- str_length(AASeqLine_str)
  AAseqLine_lastchr <- substring(AASeqLine_str, AASeqLine_length, AASeqLine_length)
  
  if(AAseqLine_lastchr == "*")
  {
    AASeqLine_str <- substring(AASeqLine_str, 1L, AASeqLine_length-1)
  }
  result <- tibble(WTAASeqLine_str = WTAASeqLine_str, AASeqLine_str = AASeqLine_str,
                   AASeqLine_length = AASeqLine_length, AAseqLine_lastchr = AAseqLine_lastchr)
  return(result)
}

#' @export  
# determine motifs
getmotifs <- function(LocusMutations.df, Test_id, CurrLocus, LocusID, AAseq, AAcompare){
  mutationspresent <- filter(LocusMutations.df, Locus_id == CurrLocus)
  
  NumMutations <- dim(mutationspresent)[1]
  
  if(NumMutations > 0)
  {
    mutationspresent$motif <- AAseq
    mutationspresent$motif <- substr(mutationspresent$motif, mutationspresent$aa_start, mutationspresent$aa_end)
    if(Test_id == "SERO")
    {
      mutationspresent$motif <- ifelse(mutationspresent$motif == mutationspresent$Mutation, 
                                       paste0(mutationspresent$aa_start, mutationspresent$Mutation),
                                       "")
      motifs <- paste0(mutationspresent$motif, collapse = "/")
      motifs <- gsub("\\/+", "/", motifs)
      motifs <- gsub("\\/$", "", motifs)
      motifs <- gsub("^\\/", "", motifs)
    } else
    {
      mutationspresent$AAcompare <- AAcompare
      mutationspresent$AAcompare <- substr(mutationspresent$AAcompare, mutationspresent$aa_start, mutationspresent$aa_end)
      mutationspresent$motif <- ifelse(str_detect(mutationspresent$Name, "INS"), 
                                       ifelse(mutationspresent$WildType == mutationspresent$motif,
                                              "WT", mutationspresent$motif),
                                       ifelse(str_detect(mutationspresent$AAcompare, "[:alpha:]"),
                                              mutationspresent$motif, "WT"))
      mutationspresent$fullmotif <- ifelse(str_detect(mutationspresent$motif, "WT"), "WT", 
                                             paste0(mutationspresent$Name, mutationspresent$motif))

      motifs <- paste0(mutationspresent$fullmotif, collapse = "/")
    }
  } else
  {
    motifs <- ""
  }
  if (LocusID != "list"){cat("\nMotif: ", motifs, "\n\n", sep = "")}
  return(motifs)
}

#' @export
# Locus Info + run BLAST on each locus to get allele numbers
locusblast <- function(CurrLocus, LocusList.df, sampleinfo, reflist, Test_id, Blast_evalue, df.bad_alleles, directorylist){
  if(sampleinfo$AlleleNo[1] != "Sample_Err")
  {
    Blast_Out_File <- blastquery(directorylist, reflist, CurrLocus, Blast_evalue)
    
    info = file.info(Blast_Out_File)
    if(info$size == 0)
    {
      AlleleNo <- "0"
      if(Test_id == "SERO")
      {
        AlleleNo <- NA  
      }
    }else
    {
      df.blastout <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
      names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", 
                              "Gaps", "SampleStart", "SampleEnd", "AlleleStart", 
                              "AlleleEnd", "eValue", "bit")
      if(Test_id %in% c("MLST", "NGMAST", "NGSTAR"))
      {
        size <- dplyr::filter(LocusList.df, Locus_id == CurrLocus)
        size <- size$size
        df.blastout100 <- filter(df.blastout, Ident == 100 & Mismatches == 0 & Gaps == 0 & Align > (size - 10))
      }else if(Test_id == "NGMAST")
      {
        min <- size$min
        df.blastout100 <- filter(df.blastout, Ident == 100 & Mismatches == 0 & Gaps == 0 & Align >= min)
      }
      {
        df.blastout100 <- filter(df.blastout, Ident == 100 & Mismatches == 0 & Gaps == 0)        
      }

      if(Test_id == "NGMAST")
      {df.blastout100 <- anti_join(df.blastout100, df.bad_alleles, by = "Allele")}  
      dfSize <- nrow(df.blastout100)
      if (dfSize > 0)
      {
        Allele <- df.blastout100$Allele[1]
        AlleleParts <- unlist(strsplit(Allele, "_"))
        if(Test_id == "NGMAST")
        {AlleleNo <- AlleleParts[3]}
        else
        {AlleleNo <- AlleleParts[2]}
      }else
      {
        AlleleNo <- "?"
        if(Test_id == "SERO") 
        {
          if(df.blastout$Ident[1] > 97)
          {
            Allele <- df.blastout$Allele[1]
            AlleleParts <- unlist(strsplit(Allele, "_"))
            AlleleNo <- AlleleParts[2]
          }
          else
          {
            AlleleNo <- "NF" #paste("NF(", Allele2, ")", sep = "")
          }
        }
        LocusFile <- paste0(reflist$WT_Dir, CurrLocus, ".fasta")
        BlastFormatCommand <- paste0("makeblastdb -in ", reflist$DestFile, " -dbtype nucl")
        try(system(BlastFormatCommand))
        blasthits(reflist, directorylist, LocusFile, "1")
        blastoutput <- readLines(paste0(directorylist$temp_dir, "blastout.txt"))
        blast_parsed <- parseblast(blastoutput, CurrLocus)

        outfile_nf <- paste0(directorylist$output_dir, "output_dna_notfound.fasta")
        sink(outfile_nf, split=FALSE, append = TRUE)
        cat(">", CurrLocus, "_", sampleinfo$CurrSampleNo, "_", "\n", blast_parsed$DNASeqLine_NoDash_str, "\n", sep ="")
        sink()
      }  
    }
  }else
  {
    AlleleNo <- "Sample_Err"
  }
  return(AlleleNo)
}

#' @export
# Locus Info + run BLAST on each locus to get mutation list
mutationblast <- function(CurrLocus, sampleinfo, directorylist, reflist, Test_id, df.bad_alleles){
  if (sampleinfo$AlleleNo != "Sample_Err")
  {
    Blast_Out_File <- blastquery(directorylist, reflist, CurrLocus, 10e-1)

    info = file.info(Blast_Out_File)
    if(info$size == 0)
    {
      Mutations <- "?"
    }else
    {
      df.blastout <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
      names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", 
                              "Gaps", "SampleStart", "SampleEnd", "AlleleStart", 
                              "AlleleEnd", "eValue", "bit")
      df.blastout100 <- filter(df.blastout, Ident == 100 & Mismatches == 0 & Gaps == 0)
      dfSize <- nrow(df.blastout100)
      if (dfSize > 0)
      {
        Allele <- df.blastout100$Allele[1]
        AlleleParts <- unlist(strsplit(Allele, "_"))
        Mutations <- AlleleParts[3]
      }else
      {
        Mutations <- "?"
      }
    }
  }else
  {
    Mutations <- "Sample_Err"
  }
  return(Mutations)
}

#' @export
# For labware upload scripts - get values for genes with pos/neg results
posneg_gene <- function(x, gene, m){
  allelename <- paste0(gene, "_allele")
  allele <- as.character(x[m, allelename])
  result <- paste0(gene, "_result")
  lw_gene_result <- as.character(x[m, result])
  if(lw_gene_result == "POS")
  {
    value <- 1L
    molec_profile <- gene
    allelevalue <- ifelse(allele == "NA", NA, paste(gene, allele, sep = " "))
  } else {
    value <- 0L
    molec_profile <- NA
    allelevalue <- NA
  }
  
  df <- tibble(lw_gene_result, value, molec_profile, allelevalue)
  names(df) <- c(paste0("lw_", gene), paste0("v_", gene), paste0("molec_profile_", gene), paste0("allele_", gene))
  return(df)
}

#' @export
# for labware upload scripts - get values for genes with different SNPs in 1 allele
allele1SNP <- function(x, gene, mutationcol, m){
  allelename <- paste0(gene, "_allele")
  allele <- as.character(x[m, allelename])
  mutations <- paste(gene, mutationcol, sep = "_")
  lw_gene_result <- as.character(x[m, mutations])

  if(lw_gene_result == "" | is.na(lw_gene_result)|
     lw_gene_result == "Sample_Err")
  {lw_gene_result <- "Err"}

  if(lw_gene_result != "WT" & lw_gene_result != "Err")
  {
    value <- 1L
    molec_profile <- paste(gene, lw_gene_result, sep = " ")
  } else if(lw_gene_result == "Err")
  {
    value <- "Err"
    molec_profile <- paste0(gene, " Err")
  } else
  {
    value <- 0L
    molec_profile <- NA
  }
  
  allelevalue <- ifelse(allele == "NA", NA, paste(gene, allele, sep = " "))
  
  df <- tibble(lw_gene_result, value, molec_profile, allelevalue)
  names(df) <- c(paste0("lw_", gene), paste0("v_", gene), paste0("molec_profile_", gene), paste0("allele_", gene))
  return(df)
}

#' @export
# for labware upload scripts - get values for genes with different SNPs in 2 alleles
allele2SNPs <- function(x, gene, genename, mutationcol, allele1, allele2, m){
  allelename <- paste0(gene, "_allele")
  allele <- as.character(x[m, allelename])
  allelename <- paste0("allele_", gene)
  mutations <- paste(gene, mutationcol, sep = "_")
  lw_gene <- paste0("lw_", gene)
  lw1 <- paste("lw", gene, allele1, sep = "_")
  lw2 <- paste("lw", gene, allele2, sep = "_")
  molec_profile <- paste0("molec_profile_", gene)
  value1 <- paste("v", gene, allele1, sep = "_")
  value2 <- paste("v", gene, allele2, sep = "_")

  df <- tibble(lw_gene = as.character(x[m, mutations]))
  names(df) <- lw_gene
  df[[lw_gene]][df[[lw_gene]] == "???" | df[[lw_gene]] == "x" | df[[lw_gene]] == "?"|
                df[[lw_gene]] == "NA" | df[[lw_gene]] == ""| df[[lw_gene]] == "Err"|
                df[[lw_gene]] == "NEG"| df[[lw_gene]] == "Sample_Err"] <- "Err/Err"

  df <- separate(df, lw_gene, c(lw1, lw2), sep = "/", remove = FALSE)
  df[[molec_profile]][df[[lw1]] != "WT" & df[[lw2]] == "WT"] <- paste(genename, df[[lw1]], sep = " ")
  df[[molec_profile]][df[[lw1]] == "WT" & df[[lw2]] != "WT"] <- paste(genename, df[[lw2]], sep = " ")
  df[[molec_profile]][df[[lw1]] != "WT" & df[[lw2]] != "WT"] <- paste0(genename, " ", df[[lw1]], "/", df[[lw2]])

  allelevalue <- ifelse(allele == "NA", NA, paste(gene, allele, sep = " "))
  df[[allelename]] <- allelevalue
  
  if(df[[lw_gene]] != "Err/Err")
  {
    df[[value1]] <- ifelse(df[[lw1]] == "WT", 0, 1)
    df[[value2]] <- ifelse(df[[lw2]] == "WT", 0, 1)
  } else
  {
    df[[value1]] <- "Err"
    df[[value2]] <- "Err"
  }
    
  return(df)
}

#' @export
# for labware upload scripts - get values for genes with different SNPs in 3 alleles
allele3SNPs <- function(x, gene, genename, mutationcol, allele1, allele2, allele3, m){
  allelename <- paste0(gene, "_allele")
  allele <- as.character(x[m, allelename])
  allelename <- paste0("allele_", gene)
  mutations <- paste(gene, mutationcol, sep = "_")
  lw_gene <- paste0("lw_", gene)
  lw1 <- paste("lw", gene, allele1, sep = "_")
  lw2 <- paste("lw", gene, allele2, sep = "_")
  lw3 <- paste("lw", gene, allele3, sep = "_")
  molec_profile <- paste0("molec_profile_", gene)
  value1 <- paste("v", gene, allele1, sep = "_")
  value2 <- paste("v", gene, allele2, sep = "_")
  value3 <- paste("v", gene, allele3, sep = "_")

  df <- tibble(lw_gene = as.character(x[m, mutations]))
  names(df) <- lw_gene
  df[[lw_gene]][df[[lw_gene]] == "???" | df[[lw_gene]] == "x" | df[[lw_gene]] == "?"|
                  df[[lw_gene]] == "NA" | df[[lw_gene]] == ""|
                  df[[lw_gene]] == "NEG"| df[[lw_gene]] == "Sample_Err"] <- "Err/Err/Err"

  df <- separate(df, lw_gene, c(lw1, lw2, lw3), sep = "/", remove = FALSE)

  df[[molec_profile]][df[[lw1]] != "WT" & df[[lw2]] == "WT" & df[[lw3]] == "WT"] <- paste(genename, df[[lw1]], sep = " ")
  df[[molec_profile]][df[[lw1]] == "WT" & df[[lw2]] != "WT" & df[[lw3]] == "WT"] <- paste(genename, df[[lw2]], sep = " ")
  df[[molec_profile]][df[[lw1]] == "WT" & df[[lw2]] == "WT" & df[[lw3]] != "WT"] <- paste(genename, df[[lw3]], sep = " ")
  df[[molec_profile]][df[[lw1]] != "WT" & df[[lw2]] != "WT" & df[[lw3]] == "WT"] <- paste0(genename, " ", df[[lw1]], "/", df[[lw2]])
  df[[molec_profile]][df[[lw1]] != "WT" & df[[lw2]] == "WT" & df[[lw3]] != "WT"] <- paste0(genename, " ", df[[lw1]], "/", df[[lw3]])
  df[[molec_profile]][df[[lw1]] == "WT" & df[[lw2]] != "WT" & df[[lw3]] != "WT"] <- paste0(genename, " ", df[[lw2]], "/", df[[lw3]])
  df[[molec_profile]][df[[lw1]] != "WT" & df[[lw2]] != "WT" & df[[lw3]] != "WT"] <- paste0(genename, " ", df[[lw1]], "/", df[[lw2]], "/", df[[lw3]])

  allelevalue <- ifelse(allele == "NA", NA, paste(gene, allele, sep = " "))
  df[[allelename]] <- allelevalue
  
  if(df[[lw_gene]] != "Err/Err/Err")
  {
    df[[value1]] <- ifelse(df[[lw1]] == "WT", 0, 1)
    df[[value2]] <- ifelse(df[[lw2]] == "WT", 0, 1)
    df[[value3]] <- ifelse(df[[lw3]] == "WT", 0, 1)
  } else
  {
    df[[value1]] <- "Err"
    df[[value2]] <- "Err"
    df[[value3]] <- "Err"
  }

  return(df)
}

#' @export
# for labware upload scripts - get values for genes with different SNPs in 4 alleles
allele4SNPs <- function(x, gene, genename, mutationcol, allele1, allele2, allele3, allele4, m){
  allelename <- paste0(gene, "_allele")
  allele <- as.character(x[m, allelename])
  allelename <- paste0("allele_", gene)
  mutations <- paste(gene, mutationcol, sep = "_")
  lw_gene <- paste0("lw_", gene)
  lw1 <- paste("lw", gene, allele1, sep = "_")
  lw2 <- paste("lw", gene, allele2, sep = "_")
  lw3 <- paste("lw", gene, allele3, sep = "_")
  lw4 <- paste("lw", gene, allele4, sep = "_")
  molec_profile <- paste0("molec_profile_", gene)
  value1 <- paste("v", gene, allele1, sep = "_")
  value2 <- paste("v", gene, allele2, sep = "_")
  value3 <- paste("v", gene, allele3, sep = "_")
  value4 <- paste("v", gene, allele4, sep = "_")

  df <- tibble(lw_gene = as.character(x[m, mutations]))
  names(df) <- lw_gene
  df[[lw_gene]][df[[lw_gene]] == "???" | df[[lw_gene]] == "x" | df[[lw_gene]] == "?"|
                  df[[lw_gene]] == "NA" | df[[lw_gene]] == ""|
                  df[[lw_gene]] == "NEG"| df[[lw_gene]] == "Sample_Err"] <- "Err/Err/Err/Err"

  df <- separate(df, lw_gene, c(lw1, lw2, lw3, lw4), sep = "/", remove = FALSE)

  df[[molec_profile]] <- NA
  if(df[[lw1]] != "WT" | df[[lw2]] != "WT" | df[[lw3]] != "WT" | df[[lw4]] != "WT"){df[[molec_profile]] <- paste0(genename, " ")}
  if(df[[lw1]] != "WT"){df[[molec_profile]] <- paste(df[[molec_profile]], df[[lw1]], sep = "")}
  if(df[[lw2]] != "WT"){df[[molec_profile]] <- paste(df[[molec_profile]], df[[lw2]], sep = "/")}
  if(df[[lw3]] != "WT"){df[[molec_profile]] <- paste(df[[molec_profile]], df[[lw3]], sep = "/")}
  if(df[[lw4]] != "WT"){df[[molec_profile]] <- paste(df[[molec_profile]], df[[lw4]], sep = "/")}

  df[[molec_profile]] <- gsub("\\s/", " ", df[[molec_profile]])

  allelevalue <- ifelse(allele == "NA", NA, paste(gene, allele, sep = " "))
  df[[allelename]] <- allelevalue
  
  if(df[[lw_gene]] != "Err/Err/Err/Err")
  {
    df[[value1]] <- ifelse(df[[lw1]] == "WT", 0, 1)
    df[[value2]] <- ifelse(df[[lw2]] == "WT", 0, 1)
    df[[value3]] <- ifelse(df[[lw3]] == "WT", 0, 1)
    df[[value4]] <- ifelse(df[[lw4]] == "WT", 0, 1)
  }else
  {
    df[[value1]] <- "Err"
    df[[value2]] <- "Err"
    df[[value3]] <- "Err"
    df[[value4]] <- "Err"
  }

  return(df)
}

#' @export
# for labware upload scripts - get values for genes with different SNPs in 5 alleles
allele5SNPs <- function(x, gene, genename, mutationcol, allele1, allele2, allele3, allele4, allele5, m){
  allelename <- paste0(gene, "_allele")
  allele <- as.character(x[m, allelename])
  allelename <- paste0("allele_", gene)
  mutations <- paste(gene, mutationcol, sep = "_")
  lw_gene <- paste0("lw_", gene)
  lw1 <- paste("lw", gene, allele1, sep = "_")
  lw2 <- paste("lw", gene, allele2, sep = "_")
  lw3 <- paste("lw", gene, allele3, sep = "_")
  lw4 <- paste("lw", gene, allele4, sep = "_")
  lw5 <- paste("lw", gene, allele5, sep = "_")
  molec_profile <- paste0("molec_profile_", gene)
  value1 <- paste("v", gene, allele1, sep = "_")
  value2 <- paste("v", gene, allele2, sep = "_")
  value3 <- paste("v", gene, allele3, sep = "_")
  value4 <- paste("v", gene, allele4, sep = "_")
  value5 <- paste("v", gene, allele5, sep = "_")

  df <- tibble(lw_gene = as.character(x[m, mutations]))
  names(df) <- lw_gene
  df[[lw_gene]][df[[lw_gene]] == "???" | df[[lw_gene]] == "x" | df[[lw_gene]] == "?"|
                  df[[lw_gene]] == "NA" | df[[lw_gene]] == ""|
                  df[[lw_gene]] == "NEG"| df[[lw_gene]] == "Sample_Err"] <- "Err/Err/Err/Err/Err"

  df <- separate(df, lw_gene, c(lw1, lw2, lw3, lw4, lw5), sep = "/", remove = FALSE)

  df[[molec_profile]] <- NA
  if(df[[lw1]] != "WT" | df[[lw2]] != "WT" | df[[lw3]] != "WT" | df[[lw4]] != "WT" | df[[lw5]] != "WT"){df[[molec_profile]] <- paste0(genename, " ")}
  if(df[[lw1]] != "WT"){df[[molec_profile]] <- paste(df[[molec_profile]], df[[lw1]], sep = "")}
  if(df[[lw2]] != "WT"){df[[molec_profile]] <- paste(df[[molec_profile]], df[[lw2]], sep = "/")}
  if(df[[lw3]] != "WT"){df[[molec_profile]] <- paste(df[[molec_profile]], df[[lw3]], sep = "/")}
  if(df[[lw4]] != "WT"){df[[molec_profile]] <- paste(df[[molec_profile]], df[[lw4]], sep = "/")}
  if(df[[lw5]] != "WT"){df[[molec_profile]] <- paste(df[[molec_profile]], df[[lw5]], sep = "/")}

  df[[molec_profile]] <- gsub("\\s/", " ", df[[molec_profile]])

  allelevalue <- ifelse(allele == "NA", NA, paste(gene, allele, sep = " "))
  df[[allelename]] <- allelevalue
  if(df[[lw_gene]] != "Err/Err/Err/Err/Err")
  {
    df[[value1]] <- ifelse(df[[lw1]] == "WT", 0, 1)
    df[[value2]] <- ifelse(df[[lw2]] == "WT", 0, 1)
    df[[value3]] <- ifelse(df[[lw3]] == "WT", 0, 1)
    df[[value4]] <- ifelse(df[[lw4]] == "WT", 0, 1)
    df[[value5]] <- ifelse(df[[lw5]] == "WT", 0, 1)
  }else
  {
    df[[value1]] <- "Err"
    df[[value2]] <- "Err"
    df[[value3]] <- "Err"
    df[[value4]] <- "Err"
    df[[value5]] <- "Err"
  }

  return(df)
}

#' @export
# For labware upload scripts - MIC calculations
MICcalc <- function(MIC_table, antibiotic, calculation, lower, upper){
  if(calculation <= lower)
  {
    MIC_calc <- round(2 ^ lower, 3)
    MIC <- paste("<=", MIC_calc, "ug/ml")     
  } else if(calculation >= upper)
  {
    MIC_calc <- round(2 ^ upper, 3)
    MIC <- paste(">=", MIC_calc, "ug/ml")
  } else
  {
    MIC_calc <- round(2 ^ calculation, 3)
    MIC <- paste0(MIC_calc, " ug/ml")
  }

  MIC_column <- paste0(antibiotic, "_MIC")
  interp <- paste0(antibiotic, "_interp")
  interpretation <- MIC_table[[antibiotic]][MIC_table$MIC == MIC_calc]
  profile <- paste0(antibiotic, "_profile")
  if(interpretation == "Resistant")
  {
    result <- paste(toupper(antibiotic), "R", sep = "-")
  } else if(interpretation == "Decreased Susceptibility")
  {
    result <- paste(toupper(antibiotic), "DS", sep = "-")
  } else if(interpretation == "Intermediate")
  {
    result <- paste(toupper(antibiotic), "I", sep = "-")
  } else
  {result <- NA}
  MIC_result <- tibble(MIC, interpretation, result)
  names(MIC_result) <- c(MIC_column, interp, profile)
  MIC_result[[antibiotic]] <- MIC

  return(MIC_result)
}
