#' Labware Upload Formatter for GONO AMR
#' January 22 2025, Walter Demczuk & Shelley Peterson
#' #' Run AMR first, then run the 23S allele counts, and then the NGSTAR-MLST analyses
#' Then run this analysis to combine data from AMR, 23S rRNA and NG-STAR
#' to prepare full amr profile to upload to LabWare.
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export
#'
#'

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#  Org_id <- "GONO"
#  curr_work_dir <- "C:/WADE/"
#-------------------------------------------------------------------------------

labware_gono_amr <- function(Org_id, curr_work_dir) {

  #-----------------------------------------------------------------------------
  # get directory structure and remove previous output files
  directorylist <- getdirectory(curr_work_dir, Org_id, "AMR")
  #-----------------------------------------------------------------------------

  # Load datafiles + combine them into one dataframe
  AMR_Output.df <- as_tibble(read.csv(paste(directorylist$output_dir, "output_profile_GONO_AMR.csv", sep = ""),
                                      header = TRUE, sep = ",", stringsAsFactors = FALSE))
  NGSTAR_Output.df <- as_tibble(read.csv(paste(directorylist$output_dir, "output_profile_mut_GONO_NGSTAR.csv", sep = ""),
                                         header = TRUE, sep = ",", stringsAsFactors = FALSE))
  rRNA23S_Output.df <- as_tibble(read.csv(paste(directorylist$output_dir, "output_profile_GONO_rRNA23S.csv", sep = ""),
                                          header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Combined_Output.df <- full_join(AMR_Output.df, NGSTAR_Output.df, by = "SampleNo")
  Combined_Output.df <- full_join(Combined_Output.df, rRNA23S_Output.df, by = "SampleNo")
  Combined_Output.df <- Combined_Output.df %>% dplyr::rename("penA_mutations" = "penA",
                                                             "mtrR_mutations" = "mtrR",
                                                             "porB_mutations" = "porB",
                                                             "ponA_mutations" = "ponA",
                                                             "gyrA_mutations" = "gyrA",
                                                             "parC_mutations" = "parC",
                                                             "rRNA23S_mutations" = "rRNA23S")

  NumSamples <- dim(Combined_Output.df)[1]

  # Load MIC chart for WGS MIC calculations
  MIC_table <- as_tibble(read.csv(paste0(directorylist$system_dir, "GONO/AMR/reference/MIC_table.csv"),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE))

  # ----------------------------------------------------------------------------
  m <- 1
  for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Start of Sample Loop
  {
    ############################################################################
    # Build Molecular Profile
    ############################################################################
    cat(m, " of ", NumSamples, "\n")
    molec_profile <- NA
    amr_profile <- NA
    lw_CurrSampleNo <- as.character(Combined_Output.df[m, "SampleNo"])

    if(Combined_Output.df[m, "SampleProfile"] == "Sample_Err")
    {
      sample_data.df <- tibble(lw_CurrSampleNo, lw_ermB = "Sample_Err", lw_ermC = "Sample_Err", 
                               lw_rpsJ = "Sample_Err", lw_tetM = "Sample_Err", 
                               lw_bla = "Sample_Err", lw_penA_prof = "Sample_Err", 
                               lw_mtrR_p = "Sample_Err", lw_mtrR_A39 = "Sample_Err", 
                               lw_mtrR_G45 = "Sample_Err", lw_porB_struct = "Sample_Err",
                               lw_porB_G120 = "Sample_Err", lw_porB_A121 = "Sample_Err", 
                               lw_ponA = "Sample_Err", lw_gyrA_S91 = "Sample_Err", 
                               lw_gyrA_D95 = "Sample_Err", lw_parC_D86 = "Sample_Err",
                               lw_parC_S87 = "Sample_Err", lw_parC_S88 = "Sample_Err", 
                               lw_parC_E91 = "Sample_Err", lw_23S_A2059G = "Sample_Err", 
                               lw_23S_C2611T = "Sample_Err", lw_16S_C1192T = "Sample_Err", 
                               lw_16S_T1458C = "Sample_Err", amr_profile = "Sample_Err", 
                               azi = "Sample_Err", azi_interp = "Sample_Err", 
                               cro = "Sample_Err", cro_interp = "Sample_Err", 
                               cfm = "Sample_Err", cfm_interp = "Sample_Err",
                               cip = "Sample_Err", cip_interp = "Sample_Err", 
                               tet = "Sample_Err", tet_interp = "Sample_Err",
                               pen = "Sample_Err", pen_interp = "Sample_Err", 
                               spe = "Sample_Err", spe_interp = "Sample_Err")
    } else
    {
      ##### ermB #####
      ermB <- posneg_gene(Combined_Output.df, "ermB", m)

      ##### ermC #####
      ermC <- posneg_gene(Combined_Output.df, "ermC", m)

      ##### rpsJ #####
      lw_rpsJ_blast <- as.character(Combined_Output.df[m, "rpsJ_result"])
      rpsJ <- allele1SNP(Combined_Output.df, "rpsJ", "motifs", m)

      ##### tetM #####
      tetM <- posneg_gene(Combined_Output.df, "tetM", m)

      ##### bla #####
      bla <- posneg_gene(Combined_Output.df, "bla", m)

      ##### rRNA 16S #####
      rRNA16S <- allele2SNPs(Combined_Output.df, "rRNA16S", "16S rRNA", "mutations", "C1192T", "T1458C", m)
      rRNA16S <- rRNA16S %>% dplyr::rename("lw_16S_C1192T" = "lw_rRNA16S_C1192T",
                                           "lw_16S_T1458C" = "lw_rRNA16S_T1458C")

      ##### penA #####
      Combined_Output.df$penA_allele <- NA
      penA <- allele5SNPs(Combined_Output.df, "penA", "penA", "mutations", "A311V", "A501", "N513Y", "A517G", "G543S", m)
      penA$lw_penA_prof <- ifelse(penA$lw_penA == "WT/WT/WT/WT/WT", "WT/WT/WT/WT/WT",
                                  sub("penA ", "", penA$molec_profile_penA))

      penA$v_penA_A501P <- ifelse(penA$lw_penA_A501 == "A501P", 1, 0)
      penA$v_penA_A501V <- ifelse(penA$lw_penA_A501 == "A501V", 1, 0)
      penA$v_penA_A501T <- ifelse(penA$lw_penA_A501 == "A501T", 1, 0)
  
      #####  mtrR #####
      Combined_Output.df$mtrR_allele <- NA
      mtrR <- allele3SNPs(Combined_Output.df, "mtrR", "mtrR", "mutations", "p", "A39", "G45", m)

      mtrR$v_mtrR_35Adel <- ifelse(mtrR$lw_mtrR_p == "-35Adel", 1, 0)
      mtrR$v_mtrR_MEN <- ifelse(mtrR$lw_mtrR_p == "MEN", 1, 0)
      mtrR$v_mtrR_Disrupted <- ifelse(mtrR$lw_mtrR_p == "Disrupted", 1, 0)
      mtrR$v_mtrR_MENDIS <- ifelse(mtrR$lw_mtrR_p == "MEN" | mtrR$lw_mtrR_p == "Disrupted", 1, 0)
      mtrR$v_mtrR_ANY <- ifelse(mtrR$lw_mtrR_p != "WT", 1, 0)

      ##### porB #####
      Combined_Output.df$porB_allele <- NA
      if(Combined_Output.df$porB_mutations[m] == "porB1a")
      {
        porB <- tibble(lw_porB = "porB1a",
                       lw_porB_struct = "porB1a",
                       lw_porB_G120 = "NA",
                       lw_porB_A121 = "NA",
                       molec_profile_porB = "NA",
                       v_porB_G120 = 0,
                       v_porB_A121 = 0,
                       v_porB_porB1b = 0,
                       v_porB_G120D = 0,
                       v_porB_G120K = 0,
                       v_porB_G120N = 0,
                       v_porB_A121D = 0,
                       v_porB_A121G = 0)
      }else
      {
        porB <- allele2SNPs(Combined_Output.df, "porB", "porB", "mutations", "G120", "A121", m)
        if(porB$lw_porB == "Err/Err")
        {
          porB$lw_porB_struct <- "Err"
          porB$v_porB_G120 <- 0
          porB$v_porB_A121 <- 0
          porB$v_porB_porB1b <- 0
          porB$v_porB_G120D <- 0
          porB$v_porB_G120K <- 0
          porB$v_porB_G120N <- 0
          porB$v_porB_A121D <- 0
          porB$v_porB_A121G <- 0
        }
        else
        {
          porB$lw_porB_struct <- "porB1b"
          porB <- relocate(porB, lw_porB, lw_porB_struct, everything())
          porB$v_porB_porB1b <- ifelse(porB$lw_porB_struct == "porB1b", 1, 0)
          porB$v_porB_G120D <- ifelse(porB$lw_porB_G120 == "G120D", 1, 0)
          porB$v_porB_G120K <- ifelse(porB$lw_porB_G120 == "G120K", 1, 0)
          porB$v_porB_G120N <- ifelse(porB$lw_porB_G120 == "G120N", 1, 0)
          porB$v_porB_A121D <- ifelse(porB$lw_porB_A121 == "A121D", 1, 0)
          porB$v_porB_A121G <- ifelse(porB$lw_porB_A121 == "A121G", 1, 0)
        }
      }
  
      ##### ponA #####
      Combined_Output.df$ponA_allele <- NA
      ponA <- allele1SNP(Combined_Output.df, "ponA", "mutations", m)

      ##### gyrA #####
      Combined_Output.df$gyrA_allele <- NA
      gyrA <- allele2SNPs(Combined_Output.df, "gyrA", "gyrA", "mutations", "S91", "D95", m)

      ##### parC #####
      Combined_Output.df$parC_allele <-NA
      parC <- allele4SNPs(Combined_Output.df, "parC", "parC", "mutations", "D86", "S87", "S88", "E91", m)

      parC$v_parC_S87R <- ifelse(parC$lw_parC_S87 == "S87R", 1, 0)
      parC$v_parC_S87I <- ifelse(parC$lw_parC_S87 == "S87I", 1, 0)
      parC$v_parC_S87C <- ifelse(parC$lw_parC_S87 == "S87C", 1, 0)
      parC$v_parC_S87N <- ifelse(parC$lw_parC_S87 == "S87N", 1, 0)

      ##### 23S #####
      rRNA23S <- tibble(lw_23S_A2059G = as.integer(Combined_Output.df[m, "A2059G"]),
                        lw_23S_C2611T = as.integer(Combined_Output.df[m, "C2611T"]),
                        v_23S_A2059G = as.integer(Combined_Output.df[m, "A2059G"]),
                        v_23S_C2611T = as.integer(Combined_Output.df[m, "C2611T"]),
                        molec_profile_23S = paste0("23S rRNA A2059G: ", as.integer(Combined_Output.df[m, "A2059G"]), "/4 ",
                                                   "C2611T: ", as.integer(Combined_Output.df[m, "C2611T"]), "/4"))

      ######################## Now put them all together #######################
      molec_profile <- paste(ermB$molec_profile_ermB, ermC$molec_profile_ermC,
                             rpsJ$molec_profile_rpsJ, tetM$molec_profile_tetM,
                             bla$molec_profile_bla, rRNA16S$molec_profile_rRNA16S,
                             penA$molec_profile_penA, mtrR$molec_profile_mtrR,
                             porB$molec_profile_porB, ponA$molec_profile_ponA,
                             gyrA$molec_profile_gyrA, parC$molec_profile_parC,
                             rRNA23S$molec_profile_23S, sep = ", ")
      molec_profile <- gsub("NA, ", "", molec_profile)
      molec_profile <- sub(" A2059G: 0/4", "", molec_profile)
      molec_profile <- sub(" C2611T: 0/4", "", molec_profile)
      molec_profile <- sub(", 23S rRNA$", "", molec_profile)
      molec_profile <- sub("23S rRNA$", "", molec_profile)

      ##########################################################################
      # Calculate MICs
      ##########################################################################

      ##### Azithromycin (AZI) #####
      azi <- NA

      if(str_detect(molec_profile, paste(c("mtrR Err", "A2059 Err", "C2611T Err", "NA/4"),collapse = '|')))
      {
        azi <- tibble(azi = "Err",
                      azi_interp = "Err",
                      azi_profile = "AZI-Err")
      }else
      {
        azi_MIC_calc <- round(-3.014+
                             (2.596*rRNA23S$lw_23S_A2059G)+
                             (1.313*rRNA23S$lw_23S_C2611T)+
                             (2.893*mtrR$v_mtrR_MEN)+
                             (5.014*mtrR$v_mtrR_Disrupted)+
                             (0.443*mtrR$v_mtrR_35Adel)+
                             (2.825*ermB$v_ermB)+
                             (4.038*ermC$v_ermC)+
                             (0.840*porB$v_porB_G120))
        azi <- MICcalc(MIC_table, "azi", azi_MIC_calc, -3L, 9L)
      }

      ##### Penicillin (PEN) #####
      pen <- NA

      if(str_detect(molec_profile, paste(c("mtrR Err", "porB Err", "ponA Err", "penA Err"),collapse = '|')))
      {
        pen <- tibble(pen = "Err",
                      pen_interp = "Err",
                      pen_profile = "PEN-Err")
      }else
      {
        pen_MIC_calc <- round(-3.44+
                             (7.11*bla$v_bla)+
                             (0.19*mtrR$v_mtrR_MEN)+
                             (0.32*mtrR$v_mtrR_G45)+
                             (0.23*penA$v_penA_A501T)+
                             (1.93*penA$v_penA_N513Y)+
                             (1.36*penA$v_penA_A517G)+
                             (0.74*penA$v_penA_G543S)+
                             (0.54*ponA$v_ponA)+
                             (0.68*porB$v_porB_G120D)+
                             (1.54*porB$v_porB_G120K)+
                             (0.86*porB$v_porB_G120N)+
                             (0.63*porB$v_porB_A121D)+
                             (0.13*porB$v_porB_A121G))
        pen <- MICcalc(MIC_table, "pen", pen_MIC_calc, -3L, 7L)
      }

      ##### Cephalosporins: Ceftriaxone (CRO/CX), Cefixime (CFM/CE) #####
      cro <- NA
      cfm <- NA

      if(str_detect(molec_profile, paste(c("mtrR Err", "porB Err", "ponA Err", "penA Err"),collapse = '|')))
      {
        cro <- tibble(cro = "Err",
                      cro_interp = "Err",
                      cro_profile = "CRO-Err")
        cfm <- tibble(cfm = "Err",
                      cfm_interp = "Err",
                      cfm_profile = "CFM-Err")
      }else
      {
        cro_MIC_calc <- round(-7.72+
                             (0.54*mtrR$v_mtrR_MENDIS)+
                             (1.38*porB$v_porB_G120)+
                             (0.67*ponA$v_ponA)+
                             (3.90*penA$v_penA_A311V)+
                             (5.15*penA$v_penA_A501P)+
                             (1.51*penA$v_penA_A501T)+
                             (1.92*penA$v_penA_A501V)+
                             (1.53*penA$v_penA_N513Y)+
                             (0.43*penA$v_penA_A517G)+
                             (0.48*penA$v_penA_G543S))
        cro <- MICcalc(MIC_table, "cro", cro_MIC_calc, -8L, 1L)

        cfm_MIC_calc <- round(-7.198+
                             (0.386*mtrR$v_mtrR_MENDIS)+
                             (0.932*mtrR$v_mtrR_G45)+
                             (0.407*porB$v_porB_G120)+
                             (5.422*penA$v_penA_A311V)+
                             (4.494*penA$v_penA_A501P)+
                             (1.463*penA$v_penA_A501T)+
                             (1.185*penA$v_penA_A501V)+
                             (4.297*penA$v_penA_N513Y)+
                             (0.497*penA$v_penA_A517G))
        cfm <- MICcalc(MIC_table, "cfm", cfm_MIC_calc, -7L, 2L)
      }

      ##### Ciprofloxacin (CIP) #####
      cip <- NA

      if(str_detect(gyrA$lw_gyrA, "Err") | str_detect(parC$lw_parC, "Err"))
      {
        cip <- tibble(cip = "Err",
                      cip_interp = "Err",
                      cip_profile = "CIP-Err")
      }else
      {
        cip_MIC_calc <- round(-7.83+
                             (7.603*gyrA$v_gyrA_S91)+
                             (2.996*parC$v_parC_D86)+
                             (3.112*parC$v_parC_S87R)+
                             (3.124*parC$v_parC_S87I)+
                             (4.223*parC$v_parC_S87C)+
                             (1.591*parC$v_parC_S88)+
                             (3.175*parC$v_parC_E91))
        cip <- MICcalc(MIC_table, "cip", cip_MIC_calc, -8L, 6L)
      }

      ##### Tetracycline (TET) #####
      tet <- NA
      tetfail <- FALSE

      if(rpsJ$lw_rpsJ == "Err")
      {
        tet <- tibble(tet = "Err",
                      tet_interp = "Err",
                      tet_profile = "TET-Err")
        tetfail <- TRUE
      }
      if(lw_rpsJ_blast == "NEG")
      {
        rpsJ$lw_rpsJ <- "NEG"

        if(tetM$lw_tetM == "NEG")
        {
          tet <- tibble(tet = "Err",
                        tet_interp = "Err",
                        tet_profile = "TET-Err")
          tetfail <- TRUE
        } else   # if TET-M is positive (big MIC), override the rpsJ Err
        {rpsJ$v_rpsJ <- 0L}
      }

      if(tetfail == FALSE)
      {
        if(str_detect(molec_profile, paste(c("mtrR Err", "porB Err"),collapse = '|')))
        {
          tet <- tibble(tet = "Err",
                        tet_interp = "Err",
                        tet_profile = "TET-Err")
        }else
        {
          tet_MIC_calc <- round(-1.64+
                               (0.33*mtrR$v_mtrR_ANY)+
                               (0.33*porB$v_porB_porB1b)+
                               (0.91*porB$v_porB_G120K)+
                               (0.12*porB$v_porB_A121)+
                               (1.73*rpsJ$v_rpsJ)+
                               (4.41*tetM$v_tetM))
          tet <- MICcalc(MIC_table, "tet", tet_MIC_calc, -3L, 7L)
        }
      }

      ##### Spectinomycin (SPE) #####
      spe <- NA
    
      spe <- tibble(spe = "<= 32 ug/ml",
                    spe_interp = "Susceptible",
                    spe_profile = "")

      if(str_detect(molec_profile, "16S rRNA Err/Err"))
      {
        spe <- tibble(spe = "Err",
                      spe_interp = "Err",
                      spe_profile = "SPE-Err")
      }else
      {
        if(rRNA16S$lw_16S_C1192T == "C1192T")
        {
          spe <- tibble(spe = ">= 128 ug/ml",
                        spe_interp = "Resistant",
                        spe_profile = "SPE-R")
        }
        if(rRNA16S$lw_16S_T1458C == "T1485C")
        {
          spe <- tibble(spe = "64 ug/ml",
                        spe_interp = "Resistant",
                        spe_profile = "SPE-R")
        }
      }

      ###################### Combine them for AMR profile ######################
      amr_profile <- paste(azi$azi_profile, pen$pen_profile, cro$cro_profile,
                             cfm$cfm_profile, cip$cip_profile, tet$tet_profile,
                             spe$spe_profile, sep = "/")
      amr_profile <- gsub("NA/", "", amr_profile)
      amr_profile <- sub("\\/$", "", amr_profile)

      ##########################################################################
      # Put all data together
      ##########################################################################
      rRNA23S$lw_23S_A2059G <- as.character(rRNA23S$lw_23S_A2059G)
      rRNA23S$lw_23S_C2611T <- as.character(rRNA23S$lw_23S_C2611T)      
      sample_data.df <- tibble(lw_CurrSampleNo, ermB$lw_ermB, ermC$lw_ermC, rpsJ$lw_rpsJ,
                               tetM$lw_tetM, bla$lw_bla, penA$lw_penA_prof, mtrR$lw_mtrR_p,
                               mtrR$lw_mtrR_A39, mtrR$lw_mtrR_G45, porB$lw_porB_struct,
                               porB$lw_porB_G120, porB$lw_porB_A121, ponA$lw_ponA,
                               gyrA$lw_gyrA_S91, gyrA$lw_gyrA_D95, parC$lw_parC_D86,
                               parC$lw_parC_S87, parC$lw_parC_S88, parC$lw_parC_E91,
                               rRNA23S$lw_23S_A2059G, rRNA23S$lw_23S_C2611T, 
                               rRNA16S$lw_16S_C1192T, rRNA16S$lw_16S_T1458C, 
                               amr_profile, azi$azi, azi$azi_interp,
                               cro$cro, cro$cro_interp, cfm$cfm, cfm$cfm_interp,
                               cip$cip, cip$cip_interp, tet$tet, tet$tet_interp,
                               pen$pen, pen$pen_interp, spe$spe, spe$spe_interp)
      sample_data.df <- sample_data.df %>% rename_at(vars(contains("$")), ~sub(".*\\$","",.))
    }
    
    if(m==1)
    {
      lw_Output.df <- tibble(sample_data.df)
    }else
    {
      lw_Output.df <- bind_rows(lw_Output.df, sample_data.df)
    }
  } # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> End of sample loop

  lw_Output.df$lw_23S_A2059G <- as.character(lw_Output.df$lw_23S_A2059G)
  lw_Output.df$lw_23S_C2611T <- as.character(lw_Output.df$lw_23S_C2611T) 
 
  # Export results to csv file
  lw_output_bad.df <- filter(lw_Output.df, str_detect(amr_profile, "Err"))
  lw_output_good.df <- filter(lw_Output.df, !str_detect(amr_profile, "Err"))

  write.csv(lw_Output.df, paste0(directorylist$output_dir, "LabWareUpload_GONO_AMR.csv"), quote = FALSE, row.names = FALSE)
  write.csv(lw_output_good.df, paste0(directorylist$output_dir, "LabWareUpload_GONO_AMR_good.csv"), quote = FALSE, row.names = FALSE)
  write.csv(lw_output_bad.df, paste0(directorylist$output_dir, "LabWareUpload_GONO_AMR_bad.csv"), quote = FALSE, row.names = FALSE)
  
  # filter and export only high MIC results
  highMICs <- lw_Output.df %>% 
    mutate_all(function(x) gsub(" ug/ml|>=|<=","",x))
  suppressWarnings(highMICs <- highMICs %>% mutate_at(c('azi', 'cro', 'cfm'), as.numeric))
  highMICs <- filter(highMICs, as.numeric(azi) > 0.5 | as.numeric(cro) > 0.05 | as.numeric(cfm) > 0.1)
  write.csv(highMICs, paste0(directorylist$output_dir, "Elevated_MICs.csv"), quote = FALSE, row.names = FALSE)
  
  cat("\n\nDone! ", directorylist$output_dir, "LabWareUpload_GONO_AMR.csv is ready in output folder", "\n\n\n", sep = "")
  
  return(lw_Output.df)
}
