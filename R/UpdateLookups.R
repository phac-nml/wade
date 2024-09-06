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
         GONO_MLST={download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/abcZ/alleles_fasta", 
                                  paste0(reflist$Lkup_Dir, "abcZ.fasta"), mode = "wb")
                    download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/adk/alleles_fasta", 
                                  paste0(reflist$Lkup_Dir, "adk.fasta"), mode = "wb")
                    download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/aroE/alleles_fasta", 
                                  paste0(reflist$Lkup_Dir, "aroE.fasta"), mode = "wb")
                    download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/fumC/alleles_fasta", 
                                  paste0(reflist$Lkup_Dir, "fumC.fasta"), mode = "wb")
                    download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/gdh/alleles_fasta", 
                                  paste0(reflist$Lkup_Dir, "gdh.fasta"), mode = "wb")
                    download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/pdhC/alleles_fasta", 
                                  paste0(reflist$Lkup_Dir, "pdhC.fasta"), mode = "wb")
                    download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/pgm/alleles_fasta", 
                                  paste0(reflist$Lkup_Dir, "pgm.fasta"), mode = "wb")
                    download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/1/profiles_csv", 
                                  paste0(reflist$Ref_Dir, "profiles.txt"), mode = "wb")},
         GONO_NGMAST={download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NG-MAST_porB/alleles_fasta", 
                                    paste0(reflist$Lkup_Dir, "porB.fasta"), mode = "wb")
                      download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/loci/NG-MAST_tbpB/alleles_fasta", 
                                    paste0(reflist$Lkup_Dir, "tbpB.fasta"), mode = "wb")
                      download.file("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/71/profiles_csv", 
                                    paste0(reflist$Ref_Dir, "profiles.txt"), mode = "wb")},
         GAS_MLST={download.file("https://rest.pubmlst.org/db/pubmlst_spyogenes_seqdef/loci/gki/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "gki.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_spyogenes_seqdef/loci/gtr/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "gtr.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_spyogenes_seqdef/loci/murI/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "murI.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_spyogenes_seqdef/loci/mutS/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "mutS.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_spyogenes_seqdef/loci/recP/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "recP.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_spyogenes_seqdef/loci/xpt/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "xpt.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_spyogenes_seqdef/loci/yqiL/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "yqiL.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_spyogenes_seqdef/schemes/1/profiles_csv", 
                                 paste0(reflist$Ref_Dir, "profiles.txt"), mode = "wb")},
         GAS_EMM={download.file("https://ftp.cdc.gov/pub/infectious_diseases/biotech/tsemm/alltrimmed.tfa", 
                                paste0(reflist$Lkup_Dir, "emm_trimmed.fasta"), mode = "wb")},
         GBS_MLST={download.file("https://rest.pubmlst.org/db/pubmlst_sagalactiae_seqdef/loci/adhP/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "adhP.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_sagalactiae_seqdef/loci/pheS/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "pheS.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_sagalactiae_seqdef/loci/atr/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "atr.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_sagalactiae_seqdef/loci/glnA/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "glnA.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_sagalactiae_seqdef/loci/sdhA/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "sdhA.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_sagalactiae_seqdef/loci/glcK/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "glcK.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_sagalactiae_seqdef/loci/tkt/alleles_fasta", 
                                 paste0(reflist$Lkup_Dir, "tkt.fasta"), mode = "wb")
                   download.file("https://rest.pubmlst.org/db/pubmlst_sagalactiae_seqdef/schemes/1/profiles_csv", 
                                 paste0(reflist$Ref_Dir, "profiles.txt"), mode = "wb")},
         PNEUMO_MLST={download.file("https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/loci/aroE/alleles_fasta", 
                                    paste0(reflist$Lkup_Dir, "aroE.fasta"), mode = "wb")
                      download.file("https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/loci/gdh/alleles_fasta", 
                                    paste0(reflist$Lkup_Dir, "gdh.fasta"), mode = "wb")
                      download.file("https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/loci/gki/alleles_fasta", 
                                    paste0(reflist$Lkup_Dir, "gki.fasta"), mode = "wb")
                      download.file("https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/loci/recP/alleles_fasta", 
                                    paste0(reflist$Lkup_Dir, "recP.fasta"), mode = "wb")
                      download.file("https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/loci/spi/alleles_fasta", 
                                    paste0(reflist$Lkup_Dir, "spi.fasta"), mode = "wb")
                      download.file("https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/loci/xpt/alleles_fasta", 
                                    paste0(reflist$Lkup_Dir, "xpt.fasta"), mode = "wb")
                      download.file("https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/loci/ddl/alleles_fasta", 
                                    paste0(reflist$Lkup_Dir, "ddl.fasta"), mode = "wb")
                      download.file("https://rest.pubmlst.org/db/pubmlst_spneumoniae_seqdef/schemes/1/profiles_csv", 
                                    paste0(reflist$Ref_Dir, "profiles.txt"), mode = "wb")}
  )

  cat("\n\nDONE.  Lookup Files Downloaded! \n")
  done_signal.df <- tibble(Output = "Lookup Files Downloaded!")
  
  return(done_signal.df)
}

