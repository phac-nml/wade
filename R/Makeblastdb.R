#' Makeblastdb.R indexes Blast fasta databases
#' February 5 2024, Walter Demczuk & Shelley Peterson
#' 
#' Index BLAST fasta datases of lookup data
#'
#' Takes Organism, Test, and Locus and create BLAST databases
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param Test_id Sample number associated with contig.fasta file
#' @param LocusID Locus_name to build blast database
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#Org_id <- "PNEUMO"                #GAS, GBS, PNEUMO or GONO
#Test_id <- "AMR_ALL"              #AMR_R, TOXINS_R, VIRULENCE_R, NGSTAR_R use MASTER for 16S id
#LocusID <- "list"
#curr_work_dir <- "C:\\WADE\\"
#Blast_evalue <- "10e-50"         #sets sensitivity of Blast gene match 10e-50 to 10e-150; use 10e-5 for primers
#-------------------------------------------------------------------------------

MakeblastDB_pipeline <- function(Org_id, Test_id, LocusID, curr_work_dir){

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  switch(Test_id,
         AMR_ALL={Test_id <- "AMR"},
         rRNA23S={Test_id <- "NGSTAR"})

  directorylist <- getdirectory(curr_work_dir, Org_id, Test_id)
  reflist <- refdirectory(directorylist, Org_id, Test_id)
  
  #-----------------------------------------------------------------------------

  if(Test_id %in% c("rRNA23S", "TOXINS_R"))
  {
    cat("No indexing required for", Test_id)
    done_signal.df <- tibble(Output = "No indexing required.")
    return(done_signal.df)
  } else
  {
    if(LocusID == "list")
    {
      LocusList.df <- as_tibble(read.csv(reflist$Loci_List, header = TRUE, sep = ",", stringsAsFactors = FALSE))
    } else
    {
      LocusList.df <- tibble(LocusID)
    }

    NumLoci <- dim(LocusList.df)[1]
    quickBlastIndex(NumLoci, LocusList.df, reflist)

    cat("\n\nDONE.  BLAST databases indexed! \n")
    done_signal.df <- tibble(Output = "Blast databases indexed.")

    return(done_signal.df)
  }
}