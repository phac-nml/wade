#MLST Blaster : Blast a WGS_contig.fasta vs. locus_library.fasta
#R Studio scripting on local version of R Studio
#2017-06-15; 2018-06-01
#Walter Demczuk

#-------------------------------------------------------------------------------------Get Inputs ####

#' NG-STAR MLST pipeline for WGS assemblies
#'
#' Takes Organism, Sample Number, Locus, and a Variable at queries a contig.fasta file
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param SampleNo Sample number (or list of sample numbers) associated with contig.fasta file
#' @param locus The locus to query, or enter list to use a list of alleles
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

# Org_id <- "GONO"
# SampleNo <- "51123"
# locus <- "list"
# curr_work_dir <- here()

NGSTAR_MLST_pipeline <- function(Org_id, SampleNo, locus, curr_work_dir) {

  #------------------------------------------------------------------------------------------------------------
  # get directory structure
  curr_dir <- curr_work_dir
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  Directories.df <- as_tibble(read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir

  #------------------------------------------------------------------------------------------------------------

Lkup_Dir <- paste(system_dir, "GONO\\NGSTAR_R\\allele_lkup_dna\\", sep = "")
Tmp_Dir <- paste(system_dir, "GONO\\NGSTAR_R\\temp\\", sep = "")
Loci_List <- paste(Tmp_Dir, "loci.csv", sep = "")
SampleListFile <- SampList
Profiles <- paste(Tmp_Dir, "profiles.csv", sep = "")
Variable <- NA

#--------------------------------------------------------------------------------------Setup locus list table ####
LocusList.df <- as_tibble(read.csv(Loci_List, header = TRUE, sep = ",", stringsAsFactors = FALSE))
if(locus != "list")
{
  LocusList.df <- filter(LocusList.df, Locus_id == locus)
}
Size.df <- dim(LocusList.df)
NumLoci <- Size.df[1]

#wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww  Index BLAST lookup files
# for(q in 1L:NumLoci)
# {
#   CurrLocus <- as.character(LocusList.df[q,1])
#   LocusLkupDNA <- paste(Lkup_Dir, CurrLocus, ".fasta", sep = "")
#   if(file.exists(LocusLkupDNA))
#   {
#     #BlastFormatCommand <- paste("formatdb -i ", LocusLkupDNA, " -p F", sep = "")
#     BlastFormatCommand <- paste("makeblastdb -in ", LocusLkupDNA, " -dbtype nucl", sep = "")
#     #try(system(BlastFormatCommand))
#     shell(BlastFormatCommand)
#   }
# }
#wwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwwww

cat("\n\n", "SampleNo NG-STAR", "\n", sep = "")

#------------------------------------------------------------------------------------Set up Sample list table
if(SampleNo == "list")
{
  SampleList.df <-as_tibble(read.csv(SampleListFile, header = TRUE, sep = ",", stringsAsFactors = FALSE))
}else
{
  SampleList.df <- tibble(SampleNo, Variable)
}

Size.df <- dim(SampleList.df)
NumSamples <- Size.df[1]

#----------------------------------------------------------------------------------------- Load sequence type profiles
profiles.df <-as_tibble(read.csv(Profiles, header = TRUE, sep = ",", stringsAsFactors = FALSE))
#------------------------------------------------------------------------------------------

m<- 1L

for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Loop through sample table
{
Allele <- ""
AlleleNo <- ""

CurrSampleNo <- as.character(SampleList.df[m, "SampleNo"])
QueryFile <- paste(ContigsDir, CurrSampleNo, ".fasta", sep = "")
DestFile <- paste0(Tmp_Dir, "queryfile.fasta")
if (!file.copy(QueryFile, DestFile, overwrite = T))
{
  Allele <- "Sample Number Error"
  AlleleNo <- "Sample_Err"

}

  n<-1L
  for (n in 1L:NumLoci)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< loop through loci table
  {

    CurrLocus <- as.character(LocusList.df[n, "Locus_id"])
    CurrLocusName <- paste("x", CurrLocus, sep = "")
    CurrLocusLen <- as.integer(LocusList.df[n, "size"])

    if (AlleleNo != "Sample_Err")
    {
    LocusLkupDNA <- paste(Lkup_Dir, CurrLocus, ".fasta", sep = "")
    Blast_Out_File <- paste0(Tmp_Dir, "blastout.txt")
    BlastCommand <- paste("blastn -query ", DestFile,
                          " -db ", LocusLkupDNA,
                          " -out ", Blast_Out_File,
                          " -num_alignments 10 -evalue 10e-1 -outfmt 6")

    shell(BlastCommand)

    info = file.info(Blast_Out_File)
    if(info$size == 0)
    {
      Allele <- "No gene present"
      AlleleNo <- "x"
      Mutations <- "?"
      df.blastout <- tibble(SampleNo, CurrLocus)
    }else
    {
      df.blastout <- as_tibble(read.csv(Blast_Out_File, header = FALSE, sep = "\t", stringsAsFactors = FALSE))
      names(df.blastout) <- c("SampleNo", "Allele", "Ident", "Align", "Mismatches", "Gaps", "SampleStart", "SampleEnd", "AlleleStart", "AlleleEnd", "eValue", "bit")
      df.blastout100 <- filter(df.blastout, Ident == 100 & Mismatches == 0 & Gaps == 0)

      dfSize <- nrow(df.blastout100)
      if (dfSize > 0)
      {
        Allele <- df.blastout100$Allele[1]
        AlleleParts <- unlist(strsplit(Allele, "_"))
        AlleleNo <- AlleleParts[2]
        Mutations <- AlleleParts[3]
      }else
      {
        Allele <- "Not Found"
        AlleleNo <- "?"
        Mutations <- "?"
      }
    }

    }#end if not sample error

    #--------------------------------------------if this is the first locus, add headers to output, else cbind next locus to output.
    if(n==1)
    {
      Output.df <- tibble(CurrSampleNo, AlleleNo)
      names(Output.df) <- c("SampleNo", CurrLocus)
      Output_mut.df <- tibble(CurrSampleNo, Mutations)
      names(Output_mut.df) <- c("SampleNo", CurrLocus)
    }else
    {
      LocusOutput.df <- tibble(AlleleNo)
      names(LocusOutput.df) <- c(CurrLocus)

      Output.df <- bind_cols(Output.df, LocusOutput.df)

      LocusOutput_mut.df <- tibble(Mutations)
      names(LocusOutput_mut.df) <- c(CurrLocus)
      Output_mut.df <- bind_cols(Output_mut.df, LocusOutput_mut.df)
    }
  }#end locus loop


if(locus == "list")
{
  #--------------------------------------------------- lookup profiles
  profile2.df <- tibble(profiles.df)
  p<-1L
  for (p in 1L:NumLoci)
  {
    if (!empty(profile2.df))
    {
      profile2.df <- filter(profile2.df, profile2.df[,p] == as.character(Output.df[1,p+1]))
    }
  }

  if (empty(profile2.df))
  {
    ST <- NA
    MLSTtype.df <- tibble(ST)
  }else
  {
    MLSTtype.df <- select(profile2.df, ST)
  }

  Output.df <- bind_cols(Output.df, MLSTtype.df)

}else
{
  MLSTtype.df <- tibble(Output.df[2])
}

#-------------------- if this is the first sample, copy output of first, else rowbind to add next sample profiles to output
if(m==1)
{
  SampleOutput.df <- tibble(Output.df)
  SampleOutput_mut.df <- tibble(Output_mut.df)
}else
{
  SampleOutput.df <- bind_rows(SampleOutput.df, Output.df)
  SampleOutput_mut.df <- bind_rows(SampleOutput_mut.df, Output_mut.df)
}

if (locus == "list")
{
  cat(CurrSampleNo, "\t", "ST-", MLSTtype.df$ST[1], "\n", sep = "")

}else
{
  cat(CurrSampleNo, "\t", locus, "\t", AlleleNo, "\t", Mutations, "\n", sep = "")
}

} #end sample loop

if (NumLoci == 1L)
{
  write.csv(df.blastout, paste(local_output_dir, "output_profile_NGSTAR.csv", sep = ""), row.names = F)
  return(df.blastout)
}else
{
  SampleOutput_good.df <- filter(SampleOutput.df, !is.na(ST))
  SampleOutput_bad.df <- filter(SampleOutput.df, is.na(ST))


  write.csv(SampleOutput.df, paste(local_output_dir, "LabWareUpload_GONO_NGSTAR.csv", sep = ""), quote = FALSE, row.names = FALSE)
  write.csv(SampleOutput_good.df, paste(local_output_dir, "LabWareUpload_GONO_NGSTAR_good.csv",sep = ""), quote = FALSE, row.names = FALSE)
  write.csv(SampleOutput_bad.df, paste(local_output_dir, "LabWareUpload_GONO_NGSTAR_bad.csv", sep = ""), quote = FALSE, row.names = FALSE)



  write.csv(SampleOutput.df, paste(local_output_dir, "output_profile_GONO_NGSTAR.csv", sep = ""), quote = FALSE,  row.names = FALSE)
  write.csv(SampleOutput_mut.df, paste(local_output_dir, "output_profile_mut.csv", sep = ""), quote = FALSE,  row.names = FALSE)
  return(SampleOutput.df)

}

cat("\n\nDone!\n\n\n")

}