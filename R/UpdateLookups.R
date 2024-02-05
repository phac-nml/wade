# Updates fasta lookup tables
#' February 5 2024, Shelley Peterson
#'
#' Downloads lookup table data from URLs
#'
#' Takes Organism, and test to update lookup tables
#' @param Test_id Test for which databases need to be updated
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#Org_Test <- "GONO_MLST"                  #GONO_MLST, GONO_NGSTAR, GONO_NGMAST, GAS_MLST, GAS_EMM, GBS_MLST, PNEUMO_MLST
#curr_work_dir <- "C:\\WADE\\"
#Blast_evalue <- "10e-50"         #sets sensitivity of Blast gene match 10e-50 to 10e-150; use 10e-5 for primers
#-------------------------------------------------------------------------------

update_lookups <- function(Org_Test, curr_work_dir){
  
  Org_id <- gsub("_.*", "", Org_Test)
  Test_id <- gsub(".*_", "", Org_Test)
  directorylist <- getdirectory(curr_work_dir, Org_id, Test_id)
  reflist <- refdirectory(directorylist, Org_id, Test_id)
  
  switch(Org_Test,
         GONO_MLST={download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=abcZ", 
                                  paste0(reflist$Lkup_Dir, "abcZ.fasta"), mode = "wb")
                    download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=adk", 
                                  paste0(reflist$Lkup_Dir, "adk.fasta"), mode = "wb")
                    download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=aroE", 
                                  paste0(reflist$Lkup_Dir, "aroE.fasta"), mode = "wb")
                    download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=fumC", 
                                  paste0(reflist$Lkup_Dir, "fumC.fasta"), mode = "wb")
                    download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=gdh", 
                                  paste0(reflist$Lkup_Dir, "gdh.fasta"), mode = "wb")
                    download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=pdhC", 
                                  paste0(reflist$Lkup_Dir, "pdhC.fasta"), mode = "wb")
                    download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=pgm", 
                                  paste0(reflist$Lkup_Dir, "pgm.fasta"), mode = "wb")
                    download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadProfiles&scheme_id=1", 
                                  paste0(reflist$Ref_Dir, "profiles.txt"), mode = "wb")},
         GONO_NGMAST={download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=NG-MAST_porB", 
                                    paste0(reflist$Lkup_Dir, "porB.fasta"), mode = "wb")
                      download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadAlleles&locus=NG-MAST_tbpB", 
                                    paste0(reflist$Lkup_Dir, "tbpB.fasta"), mode = "wb")
                      download.file("https://pubmlst.org/bigsdb?db=pubmlst_neisseria_seqdef&page=downloadProfiles&scheme_id=71", 
                                    paste0(reflist$Ref_Dir, "profiles.txt"), mode = "wb")},
         GAS_MLST={download.file("https://pubmlst.org/bigsdb?db=pubmlst_spyogenes_seqdef&page=downloadAlleles&locus=gki", 
                                 paste0(reflist$Lkup_Dir, "gki.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_spyogenes_seqdef&page=downloadAlleles&locus=gtr", 
                                 paste0(reflist$Lkup_Dir, "gtr.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_spyogenes_seqdef&page=downloadAlleles&locus=murI", 
                                 paste0(reflist$Lkup_Dir, "murI.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_spyogenes_seqdef&page=downloadAlleles&locus=mutS", 
                                 paste0(reflist$Lkup_Dir, "mutS.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_spyogenes_seqdef&page=downloadAlleles&locus=recP", 
                                 paste0(reflist$Lkup_Dir, "recP.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_spyogenes_seqdef&page=downloadAlleles&locus=xpt", 
                                 paste0(reflist$Lkup_Dir, "xpt.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_spyogenes_seqdef&page=downloadAlleles&locus=yqiL", 
                                 paste0(reflist$Lkup_Dir, "yqiL.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_spyogenes_seqdef&page=downloadProfiles&scheme_id=1", 
                                 paste0(reflist$Ref_Dir, "profiles.txt"), mode = "wb")},
         GAS_EMM={download.file("https://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/alltrimmed.tfa", 
                                paste0(reflist$Lkup_Dir, "emm_trimmed.fasta"), mode = "wb")},
         GBS_MLST={download.file("https://pubmlst.org/bigsdb?db=pubmlst_sagalactiae_seqdef&page=downloadAlleles&locus=adhP", 
                                 paste0(reflist$Lkup_Dir, "adhP.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_sagalactiae_seqdef&page=downloadAlleles&locus=pheS", 
                                 paste0(reflist$Lkup_Dir, "pheS.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_sagalactiae_seqdef&page=downloadAlleles&locus=atr", 
                                 paste0(reflist$Lkup_Dir, "atr.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_sagalactiae_seqdef&page=downloadAlleles&locus=glnA", 
                                 paste0(reflist$Lkup_Dir, "glnA.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_sagalactiae_seqdef&page=downloadAlleles&locus=sdhA", 
                                 paste0(reflist$Lkup_Dir, "sdhA.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_sagalactiae_seqdef&page=downloadAlleles&locus=glcK", 
                                 paste0(reflist$Lkup_Dir, "glcK.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_sagalactiae_seqdef&page=downloadAlleles&locus=tkt", 
                                 paste0(reflist$Lkup_Dir, "tkt.fasta"), mode = "wb")
                   download.file("https://pubmlst.org/bigsdb?db=pubmlst_sagalactiae_seqdef&page=downloadProfiles&scheme_id=1", 
                                 paste0(reflist$Ref_Dir, "profiles.txt"), mode = "wb")},
         PNEUMO_MLST={download.file("https://pubmlst.org/bigsdb?db=pubmlst_sagalactiae_seqdef&page=downloadAlleles&locus=adhP", 
                                    paste0(reflist$Lkup_Dir, "adhP.fasta"), mode = "wb")
                      download.file("https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=downloadAlleles&locus=aroE", 
                                    paste0(reflist$Lkup_Dir, "aroE.fasta"), mode = "wb")
                      download.file("https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=downloadAlleles&locus=gdh", 
                                    paste0(reflist$Lkup_Dir, "gdh.fasta"), mode = "wb")
                      download.file("https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=downloadAlleles&locus=gki", 
                                    paste0(reflist$Lkup_Dir, "gki.fasta"), mode = "wb")
                      download.file("https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=downloadAlleles&locus=recP", 
                                    paste0(reflist$Lkup_Dir, "recP.fasta"), mode = "wb")
                      download.file("https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=downloadAlleles&locus=spi", 
                                    paste0(reflist$Lkup_Dir, "spi.fasta"), mode = "wb")
                      download.file("https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=downloadAlleles&locus=xpt", 
                                    paste0(reflist$Lkup_Dir, "xpt.fasta"), mode = "wb")
                      download.file("https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=downloadAlleles&locus=ddl", 
                                    paste0(reflist$Lkup_Dir, "ddl.fasta"), mode = "wb")
                      download.file("https://pubmlst.org/bigsdb?db=pubmlst_spneumoniae_seqdef&page=downloadProfiles&scheme_id=1", 
                                    paste0(reflist$Ref_Dir, "profiles.txt"), mode = "wb")}
  )

  cat("\n\nDONE.  Lookup Files Downloaded! \n")
  done_signal.df <- tibble(Output = "Lookup Files Downloaded!")
  
  return(done_signal.df)
}

