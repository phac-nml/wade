#' Labware Upload Formatter for GAS AMR
#' February 5 2024, Walter Demczuk & Shelley Peterson
#' Run AMR first
#' Then run this analysis to combine data the full amr profile to upload to LabWare.
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#  Org_id <- "GAS"                  #GAS, PNEUMO or GONO
#  curr_work_dir <- "C:\\WADE\\"
#-------------------------------------------------------------------------------

labware_gas_amr <- function(Org_id, curr_work_dir) {

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  directorylist <- getdirectory(curr_work_dir, Org_id, "AMR")
  #-----------------------------------------------------------------------------
  
  # Load datafile
  Output.df <- as_tibble(read.csv(paste0(directorylist$output_dir, "output_profile_GAS_AMR.csv"),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE))
  NumSamples <- dim(Output.df)[1]
  NumLoci <- ((dim(Output.df)[2]-3) / 7)

  m <- 1
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Start of Sample Loop
  {
    ############################################################################
    # Build Molecular Profile
    ############################################################################
    molec_profile <- NA
    allele_profile <- NA
    amr_profile <- NA
    lw_CurrSampleNo <- as.character(Output.df[m, "SampleNo"])
      
    if(Output.df[m,"SampleProfile"] == "Sample_Err")
    {
      sample_data.df <- tibble(lw_CurrSampleNo, lw_ermA = "Sample_Err", 
                               lw_ermB = "Sample_Err", lw_ermT = "Sample_Err",
                               lw_mefAE = "Sample_Err", lw_gyrA = "Sample_Err", 
                               lw_parC = "Sample_Err", lw_tetM = "Sample_Err", 
                               lw_tetO = "Sample_Err", lw_tetT = "Sample_Err", 
                               lw_tetA = "Sample_Err", lw_tetB = "Sample_Err", 
                               lw_cat = "Sample_Err", lw_dfrF = "Sample_Err", 
                               lw_dfrG = "Sample_Err", lw_folA = "Sample_Err",
                               lw_folP = "Sample_Err", lw_pbp2x = "Sample_Err", 
                               molec_profile = "Sample_Err", 
                               ery_MIC = "Sample_Err", ery = "Sample_Err", 
                               chl_MIC = "Sample_Err", chl = "Sample_Err", 
                               lev_MIC = "Sample_Err", lev = "Sample_Err", 
                               cli_MIC = "Sample_Err", cli = "Sample_Err", 
                               tet_MIC = "Sample_Err", tet = "Sample_Err", 
                               sxt_MIC = "Sample_Err", sxt = "Sample_Err", 
                               pen_MIC = "Sample_Err", pen = "Sample_Err", 
                               amr_profile = "Sample_Err", lw_allele_profile = "Sample_Err")
    }else
    {
      ##### ermA #####
      lw_ermA <- as.character(Output.df[m, "ermA_result"])
      lw_ermAp_allele <- NA
      molec_profile_ermA <- NA

      if (lw_ermA == "POS")
      {
        # ermA promoter analysis
        lw_ermAp_result <- as.character(Output.df[m, "ermAp_result"]) #POS/NEG
        lw_ermAp_allele <- as.character(Output.df[m, "ermAp_allele"]) #ID number
        lw_ermAp_mutation <- as.character(Output.df[m, "ermAp_mutations"]) # differences from wildtype
        lw_ermAp_comment <- as.character(Output.df[m, "ermAp_comments"]) #Susceptible(Inducible)/Resistant
        lw_ermAp_motif <- as.character(Output.df[m, "ermAp_motifs"]) #amino acid substitutions from MasterBlastR

        if (lw_ermAp_motif == "WT/WT") {lw_ermAp_motif <- ""}
        if (lw_ermAp_result == "NEG")
        {
          molec_profile_ermA <- "ermA R"
          lw_ermAp_allele <- "ermAp 0"
        } else # susceptible promoter found, might have mutations for CLI-R
        {
          if (lw_ermAp_motif == "N90H/D97G")
          {
            molec_profile_ermA <- "ermA S N90H/D97G"    # CLI-Inducible
            lw_ermAp_allele <- paste0("ermAp ", lw_ermAp_allele)
          }else
          {
            molec_profile_ermA <- "ermA R"
            lw_ermAp_allele <-"ermA 0"
          }
        }
      }

      ##### ermB #####
      ermB <- posneg_gene(Output.df, "ermB", m)
      ermBSNP <- allele1SNP(Output.df, "ermB", "motifs", m)
      if(ermB$lw_ermB == "POS" & ermBSNP$lw_ermB == "N100S")
      {
        ermB$v_ermB <- 0
        ermB$molec_profile_ermB <- "ermB N100S"
      }

      ##### ermT #####
      ermT <- posneg_gene(Output.df, "ermT", m)

      ##### mefAE #####
      mefAE <- posneg_gene(Output.df, "mefAE", m)

      ##### msrD #####
      msrD <- posneg_gene(Output.df, "msrD", m) ################################ Not included in labware export, but is included in molec_profile

      ##### gyrA #####
      gyrA <- allele1SNP(Output.df, "gyrA", "motifs",m)

      ##### parC #####
      parC <- allele3SNPs(Output.df, "parC", "parC", "motifs", "D78", "S79", "D83", m)

      ##### tetM #####
      tetM <- posneg_gene(Output.df, "tetM", m)

      ##### tetO #####
      tetO <- posneg_gene(Output.df, "tetO", m)

      ##### tetT #####
      tetT <- posneg_gene(Output.df, "tetT", m)

      ##### tetA #####
      tetA <- posneg_gene(Output.df, "tetA", m)

      ##### tetB #####
      tetB <- posneg_gene(Output.df, "tetB", m)

      ##### cat #####
      cat <- posneg_gene(Output.df, "cat", m)

      ##### catQ #####
      catQ <- posneg_gene(Output.df, "catQ", m)

      ##### dfrF #####
      dfrF <- posneg_gene(Output.df, "dfrF", m)

      ##### dfrG #####
      dfrG <- posneg_gene(Output.df, "dfrG", m)

      ##### folA #####
      folA <- allele1SNP(Output.df, "folA", "motifs", m)

      ##### folP #####
      folP <- allele1SNP(Output.df, "folP", "motifs", m)

      ##### pbp2x #####
      pbp2x <- allele1SNP(Output.df, "pbp2x", "motifs", m)

      # -------------------- Vancomycin - Not included in LW export, but is included in molec_profile
      ##### vanA #####
      vanA <- posneg_gene(Output.df, "vanA", m)

      ##### vanB #####
      vanB <- posneg_gene(Output.df, "vanB", m)

      ##### vanC #####
      vanC <- posneg_gene(Output.df, "vanC", m)

      ##### vanD #####
      vanD <- posneg_gene(Output.df, "vanD", m)

      ##### vanE #####
      vanE <- posneg_gene(Output.df, "vanE", m)

      ##### vanG #####
      vanG <- posneg_gene(Output.df, "vanG", m)

      ######################## Now put them all together #######################
      molec_profile <- paste(molec_profile_ermA, ermB$molec_profile_ermB,
                             ermT$molec_profile_ermT, mefAE$molec_profile_mefAE,
                             msrD$molec_profile_msrD, gyrA$molec_profile_gyrA,
                             parC$molec_profile_parC, tetM$molec_profile_tetM,
                             tetO$molec_profile_tetO, tetT$molec_profile_tetT,
                             tetA$molec_profile_tetA, tetB$molec_profile_tetB,
                             cat$molec_profile_cat, catQ$molec_profile_catQ,
                             dfrF$molec_profile_dfrF, dfrG$molec_profile_dfrG,
                             folA$molec_profile_folA, folP$molec_profile_folP,
                             pbp2x$molec_profile_pbp2x, vanA$molec_profile_vanA,
                             vanB$molec_profile_vanB, vanC$molec_profile_vanC,
                             vanD$molec_profile_vanD, vanE$molec_profile_vanE,
                             vanG$molec_profile_vanG, sep = "; ")
      molec_profile <- gsub("NA; ", "", molec_profile)
      molec_profile <- sub("; NA", "", molec_profile)
      if(molec_profile == "NA") {molec_profile <- "Wild Type"}

      lw_allele_profile <- paste(lw_ermAp_allele, ermB$allele_ermB, ermT$allele_ermT,
                                 mefAE$allele_mefAE, msrD$allele_msrD, gyrA$allele_gyrA,
                                 parC$allele_parC, tetM$allele_tetM, tetO$allele_tetO,
                                 tetT$allele_tetT, tetA$allele_tetA, tetB$allele_tetB,
                                 cat$allele_cat, catQ$allele_catQ, dfrF$allele_dfrF,
                                 dfrG$allele_dfrG, folA$allele_folA, folP$allele_folP,
                                 pbp2x$allele_pbp2x, vanA$allele_vanA, vanB$allele_vanB,
                                 vanC$allele_vanC, vanD$allele_vanD, vanE$allele_vanE,
                                 vanG$allele_vanG, sep = ":")
      lw_allele_profile <- gsub("NA:", "", lw_allele_profile)
      lw_allele_profile <- sub(":NA", "", lw_allele_profile)

      ##########################################################################
      # Calculate MICs
      ##########################################################################

      ##### Erythromycin (ERY) #####
      if(str_detect(molec_profile, paste(c("ermA", "ermB", "ermT", "mefAE"),collapse = '|')))
      {
        ery_MIC <- ">= 2 ug/ml"
        ery <- "Resistant"
        ery_AMR <- "ERY-R"
      }else
      {
        ery_MIC <- "<= 0.25 ug/ml"
        ery <- "Susceptible"
        ery_AMR <- NA
      }

      ##### Clindamycin (CLI) #####
      if(str_detect(molec_profile, paste(c("ermB", "ermA R"), collapse = '|')))
      {
        cli_MIC <- ">= 1 ug/ml"
        cli <- "Resistant"
        cli_AMR <- "CLI-R"
      }else
      {
        cli_MIC <- "<= 0.12 ug/ml"
        cli <-"Susceptible"
        cli_AMR <- NA
      }

      if(str_detect(molec_profile, "ermB N100S"))
      {
        cli_MIC <- "<= 0.12 ug/ml"
        cli <- "Susceptible"
        cli_AMR <- NA
      }

      if(str_detect(molec_profile, "ermA S"))
      {
        cli_MIC <- "<= 0.12 ug/ml"
        cli <- "Inducible"
        cli_AMR <- "CLI-Ind"
      }

      ##### Chloramphenicol (CHL) #####
      if(str_detect(molec_profile, paste(c("cat", "catQ"),collapse = '|')))
      {
        chl_MIC <- ">= 32 ug/ml"
        chl <- "Resistant"
        chl_AMR <- "CHL-R"
      }else
      {
        chl_MIC <- "<= 4 ug/ml"
        chl <- "Susceptible"
        chl_AMR <- NA
      }

      ##### Levofloxacin (LEV) #####
      if(gyrA$lw_gyrA == "Err" | parC$lw_parC == "Err")
      {
        lev_MIC <- "Error"
        lev <- "Error"
        lev_AMR <- "LEV-Err"
      }else if(str_detect(molec_profile, "gyrA S81"))
      {
        lev_MIC <- ">= 8 ug/ml"
        lev <- "Resistant"
        lev_AMR <- "LEV-R"
      }else if((str_detect(molec_profile, "D78")))
      {
        lev_MIC <- "1 ug/ml"
        lev <- "Susceptible"
          lev_AMR <- NA
      }else if((str_detect(molec_profile, paste(c("S79", "D83"), collapse = '|'))))
      {
        lev_MIC <- "2 ug/ml"
        lev <- "Susceptible"
        lev_AMR <- NA
      }else
      {
        lev_MIC <- "<= 0.5 ug/ml"
        lev <- "Susceptible"
        lev_AMR <- NA
      }

      ##### Tetracycline (TET) #####
      if(str_detect(molec_profile, paste(c("tetM", "tetO", "tetT", "tetA", "tetB"), collapse = '|')))
      {
        tet_MIC <- ">= 8 ug/ml"
        tet <- "Resistant"
        tet_AMR <- "TET-R"
      }else
      {
        tet_MIC <- "<= 1 ug/ml"
        tet <- "Susceptible"
        tet_AMR <- NA
      }

      if(str_detect(molec_profile, "tetM Disrupted"))
      {
        tet_MIC <- "<= 1 ug/ml"
        tet <- "Susceptible"
        tet_AMR <- NA
      }

      ##### Trimethoprim/Sulfamethoxazole (SXT) #####
      if(folA$lw_folA == "Err" | folP$lw_folP == "Err")
      {
        sxt_MIC <- "Error"
        sxt <- "Error"
        sxt_AMR <- "SXT-Err"
      }else if(str_detect(molec_profile, paste(c("dfrG", "dfrF", "folA", "folP"), collapse = '|')))
      {
        if(str_detect(molec_profile, paste(c("dfrG", "dfrF"), collapse = '|')))
        {
          sxt_MIC <- ">= 4 ug/ml" #=4/76
          sxt <- "Resistant"
          sxt_AMR <- "SXT-R"
        } else if(str_detect(molec_profile, paste(c("folA", "folP"), collapse = '|')))
        {
          sxt_MIC <- "<= 0.5 ug/ml" #=0.5/9.5
          sxt <- "Susceptible"
          sxt_AMR <- NA
          if((str_detect(molec_profile, "folA")) & (str_detect(molec_profile, "folP")))
          {
            sxt_MIC <- ">= 4 ug/ml"
            sxt <- "Resistant"
            sxt_AMR <- "SXT-R"
          }
        }
      }else
      {
        sxt_MIC <- "<= 0.5 ug/ml"
        sxt <- "Susceptible"
        sxt_AMR <- NA
      }

      ##### Penicillin (PEN) #####
      if(pbp2x$lw_pbp2x == "Err")
      {
        pen_MIC <- "Error"
        pen <- "Error"
        pen_AMR <- "PEN-Err"
      }else
      {
        if(str_detect(molec_profile, "pbp2x"))
        {
          pen_MIC <- "> 0.12 ug/ml"
          pen <- "Unknown"
          pen_AMR <- "pen-NS"
        }else
        {
          pen_MIC <- "<= 0.12 ug/ml"
          pen <- "Susceptible"
          pen_AMR <- NA
        }
      }

      ##### Vancomycin (VAN) #####
      if(str_detect(molec_profile, paste(c(
         "vanA", "vanB", "vanC", "vanD", "vanE", "vanG"), collapse = '|')))
      {
        van_MIC <- ">= 2 ug/ml"
        van <- "Resistant"
      van_AMR <- "VAN-R"
      }else
      {
        van_MIC <- "<= 1 ug/ml"
        van <- "Susceptible"
        van_AMR <- NA
      }

      ###################### Combine them for AMR profile ######################
      amr_profile <- paste(ery_AMR, cli_AMR, chl_AMR, lev_AMR, tet_AMR, sxt_AMR,
                           pen_AMR, van_AMR, sep = "/")
      amr_profile <- gsub("NA/", "", amr_profile)
      amr_profile <- sub("/NA", "", amr_profile)

      if(amr_profile == "NA") {amr_profile <- "Susceptible"}

      ##########################################################################
      # Put all data together
      ##########################################################################

      # combine cat and catQ results into single column for labware. Specified will be in molecular profile.
      lw_cat <- ifelse(str_detect(molec_profile, paste(c("cat", "catQ"),collapse = '|')), "POS", "NEG")

      sample_data.df <- tibble(lw_CurrSampleNo, lw_ermA, ermB$lw_ermB,
                               ermT$lw_ermT,mefAE$lw_mefAE, gyrA$lw_gyrA,
                               parC$lw_parC, tetM$lw_tetM, tetO$lw_tetO,
                               tetT$lw_tetT, tetA$lw_tetA, tetB$lw_tetB,
                               cat$lw_cat, dfrF$lw_dfrF, dfrG$lw_dfrG,
                               folA$lw_folA, folP$lw_folP, pbp2x$lw_pbp2x,
                               molec_profile, ery_MIC, ery, chl_MIC, chl,
                               lev_MIC, lev, cli_MIC, cli, tet_MIC, tet,
                               sxt_MIC, sxt, pen_MIC, pen, amr_profile,
                               lw_allele_profile)
      sample_data.df <- sample_data.df %>% rename_at(vars(contains("$")), ~sub(".*\\$","",.))
    }

    if(m==1)  #if first sample make one row profile table, otherwise add new row to table
    {
      lw_Output.df <- tibble(sample_data.df)
    }else
    {
      lw_Output.df <- rbind(lw_Output.df, sample_data.df)
    }
  } # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> End of sample loop

  lw_output_bad.df <- filter(lw_Output.df, amr_profile == "Error")
  write.csv(lw_Output.df, paste(directorylist$output_dir, "LabWareUpload_GAS_AMR.csv", sep = ""), quote = FALSE,  row.names = FALSE)
  write.csv(lw_output_bad.df, paste(directorylist$output_dir, "LabWareUpload_GAS_AMR_bad.csv", sep = ""), quote = FALSE,  row.names = FALSE)

  cat("\n\nDone! ", directorylist$output_dir, "LabWareUpload_GAS_AMR.csv is ready in output folder", "\n\n\n", sep = "")

  return(lw_Output.df)
}
