#' Labware Upload Formatter for PNEUMO AMR
#' August 28 2025, Walter Demczuk & Shelley Peterson
#' Run AMR first, then run the 23S allele counts,
#' Then run this analysis to combine data from AMR, 23S rRNA
#' to prepare full amr profile to upload to LabWare.
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#Org_id <- "PNEUMO"
#curr_work_dir <- "C:\\WADE\\"
#-------------------------------------------------------------------------------

labware_pneumo_amr <- function(Org_id, curr_work_dir) {

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  directorylist <- getdirectory(curr_work_dir, Org_id, "AMR")
  #-----------------------------------------------------------------------------

  # Load datafiles + combine them into one dataframe
  AMR_Output.df <- as_tibble(read.csv(paste0(directorylist$output_dir, "output_profile_PNEUMO_AMR.csv"),
                             header = TRUE, sep = ",", stringsAsFactors = FALSE))
  rRNA23S_Output.df <- as_tibble(read.csv(paste0(directorylist$output_dir, "output_profile_PNEUMO_rRNA23S.csv"),
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Combined_Output.df <- full_join(AMR_Output.df, rRNA23S_Output.df, by = "SampleNo")

  NumSamples <- dim(Combined_Output.df)[1]
  NumLoci <- ((dim(Combined_Output.df)[2]-3) / 7)

  # Load MIC chart for WGS MIC calculations
  MIC_table <- as_tibble(read.csv(paste0(directorylist$system_dir, "PNEUMO/AMR/reference/MIC_table.csv"),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE))

  m <- 1
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Start of Sample Loop
  {
    ############################################################################
    # Build Molecular Profile
    ############################################################################
    molec_profile <- NA
    allele_profile <- NA
    amr_profile <- NA
    lw_CurrSampleNo <- as.character(Combined_Output.df[m, "SampleNo"])

    cat("Sample", m, " of ", NumSamples, "\n")
    if(Combined_Output.df[m,"SampleProfile"] == "Sample_Err")
    {
      sample_data.df <- tibble(lw_CurrSampleNo, lw_pbp1a = "Sample_Err", 
                               lw_pbp2b = "Sample_Err", lw_pbp2x = "Sample_Err",
                               lw_23S_A2059G = "Sample_Err", lw_23S_C2611T <- "Sample_Err",
                               lw_ermB = "Sample_Err", lw_ermTR = "Sample_Err", 
                               lw_mefAE = "Sample_Err", lw_folA = "Sample_Err", 
                               lw_folP = "Sample_Err", lw_gyrA = "Sample_Err", 
                               lw_parC = "Sample_Err", lw_tetM = "Sample_Err",
                               lw_tetO = "Sample_Err", lw_cat = "Sample_Err",
                               molec_profile = "Sample_Err", 
                               pen_MIC = "Sample_Err", pen_interp = "Sample_Err",
                               cro_MIC = "Sample_Err", cro_interp = "Sample_Err", 
                               cxm_MIC = "Sample_Err", cxm_interp = "Sample_Err", 
                               azi_MIC = "Sample_Err", azi_interp = "Sample_Err", 
                               ery_MIC = "Sample_Err", ery_interp = "Sample_Err", 
                               cla_MIC = "Sample_Err", cla_interp = "Sample_Err",
                               cli_MIC = "Sample_Err", cli_interp = "Sample_Err", 
                               chl_MIC = "Sample_Err", chl_interp = "Sample_Err",
                               lev_MIC = "Sample_Err", lev_interp = "Sample_Err",
                               mox_MIC = "Sample_Err", mox_interp = "Sample_Err",
                               tet_MIC = "Sample_Err", tet_interp = "Sample_Err",
                               dox_MIC = "Sample_Err", dox_interp = "Sample_Err",
                               sxt_MIC = "Sample_Err", sxt_interp = "Sample_Err",
                               amr_profile = "Sample_Err", lw_allele_profile = "Sample_Err")
      
    } else #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ Sample Not Error
    {
      ##### pbp1a #####
      pbp1a <- allele4SNPs(Combined_Output.df, "pbp1a", "pbp1a", "motifs", "motif1", "motif2", "motif3", "motif4", m)
      pbp1a$v_pbp1a_motif1_SAMK <- ifelse(pbp1a$lw_pbp1a_motif1 == "SAMK", 1, 0)
      pbp1a$v_pbp1a_motif1_SSMK <- ifelse(pbp1a$lw_pbp1a_motif1 == "SSMK", 1, 0)
      pbp1a$molec_profile_pbp1a[!is.na(pbp1a$molec_profile_pbp1a)] <- paste("pbp1a", pbp1a$lw_pbp1a, sep = " ")

      ##### pbp2b #####
      pbp2b <- allele4SNPs(Combined_Output.df, "pbp2b", "pbp2b", "motifs", "motif1", "motif2", "motif3", "motif4", m)
      pbp2b$molec_profile_pbp2b[!is.na(pbp2b$molec_profile_pbp2b)] <- paste("pbp2b", pbp2b$lw_pbp2b, sep = " ")

      ##### pbp2x #####
      pbp2x <- allele4SNPs(Combined_Output.df, "pbp2x", "pbp2x", "motifs", "motif1", "motif2", "motif3", "motif4", m)
      pbp2x$v_pbp2x_motif1_SAFK <- ifelse(pbp2x$lw_pbp2x_motif1 == "SAFK", 1, 0)
      pbp2x$v_pbp2x_motif1_other <- ifelse(pbp2x$v_pbp2x_motif1_SAFK == 0 & pbp2x$v_pbp2x_motif1 == 1, 1, 0)
      pbp2x$v_pbp2x_motif3_EDT <- ifelse(pbp2x$lw_pbp2x_motif3 == "EDT", 1, 0)
      pbp2x$v_pbp2x_motif3_KEA <- ifelse(pbp2x$lw_pbp2x_motif3 == "KEA", 1, 0)
      pbp2x$v_pbp2x_motif4_VKSG <- ifelse(pbp2x$lw_pbp2x_motif4 == "VKSG", 1, 0)
      pbp2x$molec_profile_pbp2x[!is.na(pbp2x$molec_profile_pbp2x)] <- paste("pbp2x", pbp2x$lw_pbp2x, sep = " ")

      ##### rRNA 23S #####
      rRNA23S <- tibble(lw_23S_A2059G = as.integer(Combined_Output.df[m, "A2059G"]),
                        lw_23S_C2611T = as.integer(Combined_Output.df[m, "C2611T"]),
                        v_23S_A2059G = as.integer(Combined_Output.df[m, "A2059G"]),
                        v_23S_C2611T = as.integer(Combined_Output.df[m, "C2611T"]),
                        lw_23S_prof = NA,
                        molec_profile_23S = NA)
      rRNA23S$lw_23S_prof <- ifelse((is.na(rRNA23S$lw_23S_A2059G)|is.na(rRNA23S$lw_23S_C2611T)), "23S Err", NA)
      rRNA23S$molec_profile_23S <- ifelse((is.na(rRNA23S$lw_23S_prof) & rRNA23S$lw_23S_A2059G > 0),
                                          paste0("23S rRNA A2059G: ", rRNA23S$lw_23S_A2059G, "/4 "), NA)
      rRNA23S$molec_profile_23S <- ifelse((is.na(rRNA23S$lw_23S_prof) & rRNA23S$lw_23S_C2611T > 0),
                                           paste0(rRNA23S$molec_profile_23S, "23S rRNA C2611T: ",
                                                  rRNA23S$lw_23S_C2611T, "/4"), NA)
      rRNA23S$molec_profile_23S <- sub("NA", "", rRNA23S$molec_profile_23S)
      rRNA23S$lw_23S_prof[is.na(rRNA23S$lw_23S_prof)] <- rRNA23S$molec_profile_23S

      ##### ermB #####
      ermB <- posneg_gene(Combined_Output.df, "ermB", m)
      ermB$mutations_ermB <- Combined_Output.df[[m, "ermB_mutations"]]
      ermB$allele_ermB <- Combined_Output.df[[m, "ermB_allele"]]
      ermB$allele_ermB[ermB$allele_ermB == "NF"] <- "ermB +"
      ermB$allele_ermB[ermB$lw_ermB == "NEG"] <- "ermB -"
      ermB$lw_ermB[ermB$mutations_ermB == "S"] <- "NEG"
      ermB$v_ermB[ermB$mutations_ermB == "S"] <- 0
      ermB$molec_profile_ermB[ermB$mutations_ermB == "S"] <- "ermB(S)"

      ##### ermTR #####
      ermTR <- posneg_gene(Combined_Output.df, "ermTR", m)

      ##### mefAE #####
      mefAE <- posneg_gene(Combined_Output.df, "mefAE", m)

      ##### folA #####
      folA <- allele1SNP(Combined_Output.df, "folA", "motifs", m)

      ##### folP #####
      folP <- allele1SNP(Combined_Output.df, "folP", "motifs", m)

      ##### gyrA #####
      gyrA <- allele1SNP(Combined_Output.df, "gyrA", "motifs", m)
      gyrA$v_gyrA_S81F <- ifelse(gyrA$lw_gyrA == "S81F", 1, 0)
      gyrA$v_gyrA_S81Y <- ifelse(gyrA$lw_gyrA == "S81Y", 1, 0)
      gyrA$v_gyrA_S81L <- ifelse(gyrA$lw_gyrA == "S81L", 1, 0)

      ##### parC #####
      parC <- allele3SNPs(Combined_Output.df, "parC", "parC", "motifs", "S79", "D83", "N91", m)

      ##### tetM #####
      tetM <- posneg_gene(Combined_Output.df, "tetM", m)

      ##### tetO #####
      tetO <- posneg_gene(Combined_Output.df, "tetO", m)

      ##### cat #####
      cat <- posneg_gene(Combined_Output.df, "cat", m)

      # -------------------- Vancomycin - Not included in LW export, but is included in molec_profile
  
      ##### vanA #####
      vanA <- posneg_gene(Combined_Output.df, "vanA", m)

      ##### vanB #####
      vanB <- posneg_gene(Combined_Output.df, "vanB", m)

      ##### vanC #####
      vanC <- posneg_gene(Combined_Output.df, "vanC", m)

      ##### vanD #####
      vanD <- posneg_gene(Combined_Output.df, "vanD", m)
  
      ##### vanE #####
      vanE <- posneg_gene(Combined_Output.df, "vanE", m)

      ##### vanG #####
      vanG <- posneg_gene(Combined_Output.df, "vanG", m)

      ######################## Now put them all together #######################
      molec_profile <- paste(ermB$molec_profile_ermB, ermTR$molec_profile_ermTR,
                             mefAE$molec_profile_mefAE, rRNA23S$molec_profile_23S,
                             tetM$molec_profile_tetM, tetO$molec_profile_tetO,
                             cat$molec_profile_cat, folA$molec_profile_folA,
                             folP$molec_profile_folP, pbp1a$molec_profile_pbp1a,
                             pbp2b$molec_profile_pbp2b, pbp2x$molec_profile_pbp2x,
                             gyrA$molec_profile_gyrA, parC$molec_profile_parC,
                             vanA$molec_profile_vanA, vanB$molec_profile_vanB,
                             vanC$molec_profile_vanC, vanD$molec_profile_vanD,
                             vanE$molec_profile_vanE, vanG$molec_profile_vanG, sep = "; ")
      molec_profile <- gsub("NA; ", "", molec_profile)
      molec_profile <- sub("; NA", "", molec_profile)

      if(molec_profile == "NA") {molec_profile <- "Wild Type"}
  
      lw_allele_profile <- paste(pbp1a$allele_pbp1a, pbp2b$allele_pbp2b,
                                 pbp2x$allele_pbp2x, ermB$allele_ermB,
                                 folA$allele_folA, folP$allele_folP,
                                 gyrA$allele_gyrA, parC$allele_parC, sep = " : ")
      lw_allele_profile <- gsub("NA : ", "", lw_allele_profile)
      lw_allele_profile <- sub(" : NA", "", lw_allele_profile)

      ##########################################################################
      # Calculate MICs
      ##########################################################################

      ##### Levofloxacin (LEV) #####
      lev <- NA

      if(str_detect(gyrA$lw_gyrA, "Err") | str_detect(parC$lw_parC, "Err"))
      {
        lev <- tibble(lev_MIC = "Error",
                      lev_interp = "Error",
                      lev_profile = "Error")
      }else
      {
        lev_MIC_calc <-  round(-0.218 +
                              (gyrA$v_gyrA_S81F * 2.028)+
                              (gyrA$v_gyrA_S81Y * 1.564)+
                              (gyrA$v_gyrA_S81L * 3.564)+
                              (parC$v_parC_S79 * 1.654)+
                              (parC$v_parC_D83 * 0.834))
        lev <-  MICcalc(MIC_table, "lev", lev_MIC_calc, 0L, 5L)
      }

      ##### Moxifloxacin (MOX) #####
      mox <- NA

      if(str_detect(gyrA$lw_gyrA, "Err") | str_detect(parC$lw_parC, "Err"))
      {
        mox <- tibble(mox_MIC = "Error",
                      mox_interp = "Error",
                      mox_profile = "Error")
      }else
      {
        mox_MIC_calc <-  round(-2.819 +
                              (gyrA$v_gyrA_S81F * 3.130)+
                              (gyrA$v_gyrA_S81Y * 3.907)+
                              (gyrA$v_gyrA_S81L * 4.907)+
                              (parC$v_parC_S79 * 0.911))
        mox <-  MICcalc(MIC_table, "mox", mox_MIC_calc, 0, 5L)
      }

      ##### Penicillin (PEN) #####
      pen <- NA

      if(str_detect(pbp1a$lw_pbp1a, "Err") | str_detect(pbp2b$lw_pbp2b, "Err")|
         str_detect(pbp2x$lw_pbp2x, "Err"))
      {
        pen <- tibble(pen_MIC = "Error",
                      pen_interp = "Error",
                      pen_profile = "Error")
      }else
      {
        pen_MIC_calc <- round(-4.6099 +
                             (pbp1a$v_pbp1a_motif1 * 1.5466)+
                             (pbp1a$v_pbp1a_motif4 * 0.9491)+
                             (pbp2b$v_pbp2b_motif2 * 1.2025)+
                             (pbp2b$v_pbp2b_motif3 * 0.3563)+
                             (pbp2x$v_pbp2x_motif1_SAFK * 1.6259)+
                             (pbp2x$v_pbp2x_motif3_EDT * 1.5476)+
                             (pbp2x$v_pbp2x_motif3_KEA * 0.6760)+
                             (pbp2x$v_pbp2x_motif4_VKSG * 0.7536))
        pen <-  MICcalc(MIC_table, "pen", pen_MIC_calc, -5, 2L)
        pen$pen_MIC[pen$pen_MIC == "<= 0.031 ug/ml"] <- "<= 0.03 ug/ml"
      }

      ##### Cephalosporins: Ceftriaxone (CRO), Cefuroxime (CXM)
      cro <- NA
      cxm <- NA

      if(str_detect(pbp1a$lw_pbp1a, "Err") | str_detect(pbp2x$lw_pbp2x, "Err"))
      {
        cro <- tibble(cro_MIC = "Error",
                      cro_interp = "Error",
                      cro_profile = "Error")
        cxm <- tibble(cxm_MIC = "Error",
                      cxm_interp = "Error",
                      cxm_profile = "Error")
      }else
      {
        cro_MIC_calc <- round(-2.709 +
                             (pbp1a$v_pbp1a_motif1 * 1.25)+
                             (pbp2x$v_pbp2x_motif1_SAFK * 2.72)+
                             (pbp2x$v_pbp2x_motif3_EDT * 0.76)+
                             (pbp2x$v_pbp2x_motif4_VKSG * 0.989))
        cro <- MICcalc(MIC_table, "cro", cro_MIC_calc, -3, 2L)

        cxm_MIC_calc <- round(-1.018 +
                             (pbp1a$v_pbp1a_motif1_SAMK * 1.509)+
                             (pbp1a$v_pbp1a_motif1_SSMK * 2.170)+
                             (pbp2x$v_pbp2x_motif1_SAFK * 2.322)+
                             (pbp2x$v_pbp2x_motif1_other * 0.256)+
                             (pbp2x$v_pbp2x_motif4_VKSG * 1.026))
        cxm <- MICcalc(MIC_table, "cxm", cxm_MIC_calc, -1, 3L)
      }

      ##### Trimethoprim/Sulfamethoxazole (SXT) #####
      sxt <- NA

      if(str_detect(folA$lw_folA, "Err") | str_detect(folP$lw_folP, "Err"))
      {
        sxt <- tibble(sxt_MIC = "Error",
                      sxt_interp = "Error",
                      sxt_profile = "Error")
      }else
      {
        sxt_MIC_calc <- round(-2.265 +
                             (folA$v_folA * 2.113)+
                             (folP$v_folP * 2.668))
        sxt <- MICcalc(MIC_table, "sxt", sxt_MIC_calc, -1, 4L)
        sxt$sxt_MIC[sxt_MIC_calc < 0] <- "<= 0.5/9.5 ug/ml"
        sxt$sxt_MIC[sxt_MIC_calc > 1] <- ">= 4/76 ug/ml"
      }
      
      ##### Macrolides: Erythromycin (ERY), Azithromycin (AZI)  #####
      #####             Clarithromycin (CLA), Clindamycin (CLI) #####
      ery <- NA
      azi <- NA
      cla <- NA
      cli <- NA

      if(str_detect(rRNA23S$lw_23S_A2059G, "Err") | str_detect(rRNA23S$lw_23S_C2611T, "Err")|
         str_detect(ermB$lw_ermB, "Err") | str_detect(mefAE$lw_mefAE, "Err"))
      {
        ery <- tibble(ery_MIC = "Error",
                      ery_interp = "Error",
                      ery_profile = "Error")
        azi <- tibble(azi_MIC = "Error",
                      azi_interp = "Error",
                      azi_profile = "Error")
        cla <- tibble(cla_MIC = "Error",
                      cla_interp = "Error",
                      cla_profile = "Error")
        cli <- tibble(cli_MIC = "Error",
                      cli_interp = "Error",
                      cli_profile = "Error")
      }else
      {
        # if both ermB and mef are positive, only use ermB value to calculate
        mefAE$v_mefAE[1] <- ifelse(ermB$lw_ermB[1] == "POS", 0, mefAE$v_mefAE[1])

        ery_MIC_calc <- round(-2.975 +
                             (rRNA23S$v_23S_A2059G * 1.993)+
                             (ermB$v_ermB * 7.680)+
                             (mefAE$v_mefAE * 4.808))
        ery <- MICcalc(MIC_table, "ery", ery_MIC_calc, -3, 2L)

        azi_MIC_calc <- round(-1.9722 +
                             (rRNA23S$v_23S_A2059G * 1.1204)+
                             (rRNA23S$v_23S_C2611T * 0.9917)+
                             (ermB$v_ermB * 3.9722)+
                             (mefAE$v_mefAE * 3.8122))
        azi <- MICcalc(MIC_table, "azi", azi_MIC_calc, -2, 2L)

        cla_MIC_calc <- round(-4.984 +
                             (rRNA23S$v_23S_A2059G * 1.819)+
                             (rRNA23S$v_23S_C2611T * 1.246)+
                             (ermB$v_ermB * 10.820)+
                             (mefAE$v_mefAE * 6.064))
        cla <-  MICcalc(MIC_table, "cla", cla_MIC_calc, -5, 6L)
        cla$cla_MIC[cla$cla_MIC == "<= 0.031 ug/ml"] <- "<= 0.03 ug/ml"

        cli_MIC_calc <- round(-2.8145 +
                             (rRNA23S$v_23S_A2059G * 0.456)+
                             (ermB$v_ermB * 9.048))
        cli <-  MICcalc(MIC_table, "cli", cli_MIC_calc, -3, 1L)
        cli$cli_interp <- ifelse(ermTR$lw_ermTR == "POS", "Inducible", cli$cli_interp)
        if(cli$cli_interp == "Inducible")
        {
          cli$cli_profile <- "CLI-Ind"
        }
      }

      ##### Chloramphenicol (CHL) #####
      chl <- NA
  
      if(str_detect(cat$lw_cat, "Err"))
      {
        chl <- tibble(chl_MIC = "Error",
                      chl_interp = "Error",
                      chl_profile = "Error")
      }else if(cat$lw_cat == "POS")
      {
        chl <- tibble(chl_MIC = ">= 8 ug/ml",
                      chl_interp = "Resistant",
                      chl_profile = "CHL-R")
      }else
      {
        chl <- tibble(chl_MIC = "<= 4 ug/ml",
                      chl_interp = "Susceptible",
                      chl_profile = NA)
      }
  
      ##### Tetracycline (TET), Doxycycline (DOX) #####
      tet <- NA
      dox <- NA

      if(str_detect(tetM$lw_tetM, "Err"))
      {
        tet <- tibble(tet_MIC = "Error",
                      tet_interp = "Error",
                      tet_profile = "Error")
        dox <- tibble(dox_MIC = "Error",
                      dox_interp = "Error",
                      dox_profile = "Error")
      }else if(tetM$lw_tetM == "POS")
      {
        tet <- tibble(tet_MIC = ">= 16 ug/ml",
                      tet_interp = "Resistant",
                      tet_profile = "TET-R")
        dox <- tibble(dox_MIC = ">= 4 ug/ml",
                      dox_interp = "Resistant",
                      dox_profile = "DOX-R")
      }else
      {
        tet <- tibble(tet_MIC = "<= 1 ug/ml",
                      tet_interp = "Susceptible",
                      tet_profile = NA)
        dox <- tibble(dox_MIC = "<= 0.25 ug/ml",
                      dox_interp = "Susceptible",
                      dox_profile = NA)
      }

      ##### Vancomycin (VAN) #####
      van <- NA

      if(str_detect(molec_profile, paste(c("vanA", "vanB", "vanC", "vanD", "vanE", "vanG"), collapse = '|')))
      {
        van <- tibble(van_MIC = ">= 2 ug/ml",
                      van_interp = "Resistant",
                      van_profile = "VAN-R")
      }else
      {
        van <- tibble(van_MIC = "<= 1 ug/ml",
                      van_interp = "Susceptible",
                      van_profile = NA)
      }

      ###################### Combine them for AMR profile ######################
      amr_profile <- paste(pen$pen_profile, cro$cro_profile, cxm$cxm_profile,
                           ery$ery_profile, azi$azi_profile, cla$cla_profile,
                           cli$cli_profile, chl$chl_profile, lev$lev_profile,
                           mox$mox_profile, tet$tet_profile, dox$dox_profile,
                           sxt$sxt_profile, van$van_profile,  sep = "/")
      amr_profile <- gsub("NA/", "", amr_profile)
      amr_profile <- gsub("NA", "", amr_profile)
      amr_profile <- sub("\\/$", "", amr_profile)
      amr_profile <- ifelse(amr_profile == "", "Susceptible", amr_profile)

      ##########################################################################
      # Put all data together
      ##########################################################################
      sample_data.df <- tibble(lw_CurrSampleNo, pbp1a$lw_pbp1a, pbp2b$lw_pbp2b, pbp2x$lw_pbp2x,
                               as.character(rRNA23S$lw_23S_A2059G), 
                               as.character(rRNA23S$lw_23S_C2611T), ermB$lw_ermB,
                               ermTR$lw_ermTR, mefAE$lw_mefAE, folA$lw_folA, folP$lw_folP,
                               gyrA$lw_gyrA, parC$lw_parC, tetM$lw_tetM, tetO$lw_tetO,
                               cat$lw_cat, molec_profile, pen$pen_MIC, pen$pen_interp,
                               cro$cro_MIC, cro$cro_interp, cxm$cxm_MIC, cxm$cxm_interp,
                               azi$azi_MIC, azi$azi_interp, ery$ery_MIC, ery$ery_interp,
                               cla$cla_MIC, cla$cla_interp, cli$cli_MIC, cli$cli_interp,
                               chl$chl_MIC, chl$chl_interp, lev$lev_MIC, lev$lev_interp,
                               mox$mox_MIC, mox$mox_interp, tet$tet_MIC, tet$tet_interp,
                               dox$dox_MIC, dox$dox_interp, sxt$sxt_MIC, sxt$sxt_interp,
                               amr_profile, lw_allele_profile)
      sample_data.df <- sample_data.df %>% rename_at(vars(contains("$")), ~sub(".*\\$","",.))
    } #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ End of Sample Not Error loop
    
    if(m==1)
    {
      lw_output.df <- tibble(sample_data.df)
    }else
    {
      lw_output.df <- bind_rows(lw_output.df, sample_data.df)
    }
  } #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> End of sample loop

  lw_output_bad.df <- filter(lw_output.df, amr_profile == "Error")
  lw_output_good.df <- filter(lw_output.df, amr_profile != "Error")

  write.csv(lw_output.df, paste0(directorylist$output_dir, "LabWareUpload_PNEUMO_AMR.csv"), quote = FALSE, row.names = FALSE)
  write.csv(lw_output_good.df, paste0(directorylist$output_dir, "LabWareUpload_PNEUMO_AMR_good.csv"), quote = FALSE, row.names = FALSE)
  write.csv(lw_output_bad.df, paste0(directorylist$output_dir, "LabWareUpload_PNEUMO_AMR_bad.csv"), quote = FALSE, row.names = FALSE)

  cat("\n\nDone! ", directorylist$output_dir, "LabWareUpload_PNEUMO_AMR.csv is ready in output folder", "\n\n\n", sep = "")

  return(lw_output.df)

} 