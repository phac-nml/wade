#' Run AMR first, then run the 23S allele counts, and then the NGSTAR-MLST analyses
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

labware_gono_amr <- function(Org_id, curr_work_dir) {

  #------------------------------------------------------------------------------------------------------------
  # get directory structure from DirectoryLocations.csv file
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
  Directories.df <- as_tibble(read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir

AMR_Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_GONO_AMR.csv", sep = ""),
                                    header = TRUE, sep = ",", stringsAsFactors = FALSE))
NGSTAR_Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_mut.csv", sep = ""),
                                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
rRNA23S_Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_23S.csv", sep = ""),
                                        header = TRUE, sep = ",", stringsAsFactors = FALSE))
Combined_Output_first.df <- full_join(AMR_Output.df, NGSTAR_Output.df, by = "SampleNo")
Combined_Output.df <- full_join(Combined_Output_first.df, rRNA23S_Output.df, by = "SampleNo")

Size.df <- dim(Combined_Output.df)
NumSamples <- Size.df[1]

tet_mic.df <- as_tibble(read.csv(paste(system_dir, "GONO\\NGSTAR_R\\temp\\inc_mic_tetracycline.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))

azi_mic.df <- as_tibble(read.csv(paste(system_dir, "GONO\\NGSTAR_R\\temp\\inc_mic_azithromycin.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))

pen_mic.df <- as_tibble(read.csv(paste(system_dir, "GONO\\NGSTAR_R\\temp\\inc_mic_penicillin.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))

cfx_mic.df <- as_tibble(read.csv(paste(system_dir, "GONO\\NGSTAR_R\\temp\\inc_mic_ceftriaxone.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))

cfm_mic.df <- as_tibble(read.csv(paste(system_dir, "GONO\\NGSTAR_R\\temp\\inc_mic_cefixime.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))

cip_mic.df <- as_tibble(read.csv(paste(system_dir, "GONO\\NGSTAR_R\\temp\\inc_mic_ciprofloxacin.csv", sep = ""),
                       header = TRUE, sep = ",", stringsAsFactors = FALSE))
sepr <- " - "

m <- 1
ErrorFound <- FALSE

for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< Build Molecular Profile
{
  molec_profile <- NA

  lw_CurrSampleNo <- as.character(Combined_Output.df[m, "SampleNo"])

  lw_ermB <- as.character(Combined_Output.df[m, "ermB_result"])
    if (lw_ermB == "POS")
      {
      molec_profile <- "ermB"
      V_ermB <- 1L} else (v_ermB <- 0L)

  lw_ermC <- as.character(Combined_Output.df[m, "ermC_result"])
      if (lw_ermC == "POS")
      {
        if (is.na(molec_profile))
        {molec_profile <- "ermC"} else
        {molec_profile <- paste(molec_profile, sepr, "ermC", sep = "")}
        v_ermC <- 1L
      } else (v_ermC <- 0L)

  lw_rpsJ <- as.character(Combined_Output.df[m, "rpsJ_motifs"])

  lw_rpsJ_blast <- as.character(Combined_Output.df[m, "rpsJ_result"])

  if (lw_rpsJ == "" | is.na(lw_rpsJ)) {lw_rpsJ <- "Err"}
  if (lw_rpsJ != "WT" & lw_rpsJ != "Err")
    {
      v_rpsJ <- 1L
      if (is.na(molec_profile))
      {molec_profile <- paste("rpsJ ", lw_rpsJ, sep = "")} else
      {molec_profile <- paste(molec_profile, sepr, "rpsJ ", lw_rpsJ, sep = "")}
    }else {v_rpsJ <- 0L}

  lw_tetM <- as.character(Combined_Output.df[m, "tetM_result"])
  if (lw_tetM == "POS")
  {
    v_tetM <- 1L
    if (is.na(molec_profile))
    {molec_profile <- "tetM"} else
    {molec_profile <- paste(molec_profile, sepr, "tetM", sep = "")}
  }else (v_tetM <- 0L)

  #-------------------------------------------------------------------------------------- bla

  lw_bla <- as.character(Combined_Output.df[m, "bla_result"])
  if (lw_bla == "POS")
  {
    v_bla <- 1L

    if (is.na(molec_profile))
    {molec_profile <- "bla"} else
    {molec_profile <- paste(molec_profile, sepr, "bla", sep = "")}
  }else {v_bla <- 0L}
  #--------------------------------------------------------------------------------------16S
  lw_16S <-   as.character(Combined_Output.df[m, "X16S_mutations"])

  if (lw_16S == "?" | lw_16S == "x" | lw_16S == "???" | lw_16S == "NA" | lw_16S == "") {lw_16S <- "Err/Err"}
  lw_16S_parts <- unlist(strsplit(lw_16S, "/"))
  lw_16S_C1192T <-  lw_16S_parts[1]
  lw_16S_T1458C <- lw_16S_parts[2]

  lw_16S_prof <- NA

  if (lw_16S_C1192T != "WT")
  {
    if (is.na(lw_16S_prof))
    {lw_16S_prof <- paste("16S rRNA ", lw_16S_C1192T, sep = "")} else
    {lw_16S_prof <- paste(lw_16S_prof, "/", lw_16S_C1192T, sep = "")}
  }
  if (lw_16S_T1458C != "WT")
  {
    if (is.na(lw_16S_prof))
    {lw_16S_prof <- paste("23S rRNA ", lw_16S_T1458C, sep = "")} else
    {lw_16S_prof <- paste(lw_16S_prof, "/", lw_16S_T1458C, sep = "")}
  }

  if (!is.na(lw_16S_prof))
  {
    if (is.na(molec_profile))
    {molec_profile <- lw_16S_prof} else
    {molec_profile <- paste(molec_profile, sepr, lw_16S_prof, sep = "")}
  }

  #---------------------------------------------------------------------------  NG-STAR penA

  lw_penA <-   as.character(Combined_Output.df$penA[m])
  lw_penA_prof <- NA
  if (lw_penA == "?" | lw_penA == "x") {lw_penA <- "Err/Err/Err/Err/Err"}

  lw_penA_parts <- unlist(strsplit(lw_penA, "/"))
  lw_penA_A311V <-  lw_penA_parts[1]
  lw_penA_A501 <- lw_penA_parts[2]
  lw_penA_N513Y <- lw_penA_parts[3]
  lw_penA_A517G <- lw_penA_parts[4]
  lw_penA_G543S <- lw_penA_parts[5]

  if (lw_penA_A311V != "WT")
  {
    if (is.na(lw_penA_prof))
    {lw_penA_prof <- lw_penA_A311V} else
    {lw_penA_prof <- paste(lw_penA_prof, "/", lw_penA_A311V, sep = "")}
  }
  if (lw_penA_A501 != "WT")
  {
    if (is.na(lw_penA_prof))
    {lw_penA_prof <- lw_penA_A501} else
    {lw_penA_prof <- paste(lw_penA_prof, "/", lw_penA_A501, sep = "")}
  }
  if (lw_penA_N513Y != "WT")
  {
    if (is.na(lw_penA_prof))
    {lw_penA_prof <- lw_penA_N513Y} else
    {lw_penA_prof <- paste(lw_penA_prof, "/", lw_penA_N513Y, sep = "")}
  }
  if (lw_penA_A517G != "WT")
  {
    if (is.na(lw_penA_prof))
    {lw_penA_prof <- lw_penA_A517G} else
    {lw_penA_prof <- paste(lw_penA_prof, "/", lw_penA_A517G, sep = "")}
  }
  if (lw_penA_G543S != "WT")
  {
    if (is.na(lw_penA_prof))
    {lw_penA_prof <- lw_penA_G543S} else
    {lw_penA_prof <- paste(lw_penA_prof, "/", lw_penA_G543S, sep = "")}
  }

  if (!is.na(lw_penA_prof))
  {
    if (is.na(molec_profile))
    {molec_profile <- paste("penA ", lw_penA_prof, sep = "")} else
    {molec_profile <- paste(molec_profile, sepr, "penA ", lw_penA_prof, sep = "")}
  } else
  {
    lw_penA_prof <- lw_penA
  }

  ifelse(str_detect(lw_penA_prof, "A311V"), (v_penA_A311V <- 1L), (v_penA_A311V <- 0L))
  ifelse(str_detect(lw_penA_prof, "A501P"), (v_penA_A501P <- 1L), (v_penA_A501P <- 0L))
  ifelse(str_detect(lw_penA_prof, "A501V"), (v_penA_A501V <- 1L), (v_penA_A501V <- 0L))
  ifelse(str_detect(lw_penA_prof, "A501T"), (v_penA_A501T <- 1L), (v_penA_A501T <- 0L))
  ifelse(str_detect(lw_penA_prof, "N513Y"), (v_penA_N513Y <- 1L), (v_penA_N513Y <- 0L))
  ifelse(str_detect(lw_penA_prof, "A517G"), (v_penA_A517G <- 1L), (v_penA_A517G <- 0L))
  ifelse(str_detect(lw_penA_prof, "G543S"), (v_penA_G543S <- 1L), (v_penA_G543S <- 0L))

  #---------------------------------------------------------------------------  mtrR
  lw_mtrR <-   as.character(Combined_Output.df[m, "mtrR"])
  lw_mtrR_prof <- NA

    if (lw_mtrR == "?" | lw_mtrR == "x") {lw_mtrR <- "Err/Err/Err"}

    lw_mtrR_parts <- unlist(strsplit(lw_mtrR, "/"))
    lw_mtrR_p <-  lw_mtrR_parts[1]
    lw_mtrR_A39 <- lw_mtrR_parts[2]
    lw_mtrR_G45 <- lw_mtrR_parts[3]

    if (lw_mtrR_p != "WT")
    {
      if (is.na(lw_mtrR_prof))
      {lw_mtrR_prof <- paste("mtrR ", lw_mtrR_p, sep = "")} else
      {lw_mtrR_prof <- paste(lw_mtrR_prof, "mtrR ", lw_mtrR_p, sep = "")}
    }
    if (lw_mtrR_A39 == "A39T")
    {
      if (is.na(lw_mtrR_prof))
      {lw_mtrR_prof <- "mtrR A39T"} else
      {lw_mtrR_prof <- paste(lw_mtrR_prof, "/A39T", sep = "")}
    }
    if (lw_mtrR_G45 == "G45D")
    {
      if (is.na(lw_mtrR_prof))
      {lw_mtrR_prof <- "mtrR G45D"} else
      {lw_mtrR_prof <- paste(lw_mtrR_prof, "/G45D", sep = "")}
    }
    if (!is.na(lw_mtrR_prof))
    {
      if (is.na(molec_profile))
      {molec_profile <- lw_mtrR_prof} else
      {molec_profile <- paste(molec_profile, sepr, lw_mtrR_prof, sep = "")}
    }

  if (is.na(lw_mtrR_prof)) {lw_mtrR_prof <- ""}
  if (str_detect(lw_mtrR_prof, "-35Adel")){v_mtrR35Adel <- 1L}else {v_mtrR35Adel <- 0L}
  if (str_detect(lw_mtrR_prof, "MEN")){v_mtrRMEN <- 1L}else {v_mtrRMEN <- 0L}
  if (str_detect(lw_mtrR_prof, "Disrupted")){v_mtrRDisrupted <- 1L}else {v_mtrRDisrupted <- 0L}
  if (str_detect(lw_mtrR_prof, "A39T")){v_mtrR_A39T <- 1L}else {v_mtrR_A39T <- 0L}
  if (str_detect(lw_mtrR_prof, "G45D")){v_mtrR_G45D <- 1L}else {v_mtrR_G45D <- 0L}

  if (str_detect(lw_mtrR_prof, "MEN") || str_detect(lw_mtrR_prof, "Disrupted")
      ){v_mtrRMENDIS <- 1L}else {v_mtrRMENDIS <- 0L}

  if (str_detect(lw_mtrR_prof, "MEN") || str_detect(lw_mtrR_prof, "Disrupted") || str_detect(lw_mtrR_prof, "-35Adel")
    ){v_mtrRANY <- 1L}else {v_mtrRANY <- 0L}

  #--------------------------------------------------------------------------------------porB
  lw_porB <-   as.character(Combined_Output.df[m, "porB"])

  if (lw_porB == "?" | lw_porB == "x")
    {
    lw_porB <- "Err/Err"
    lw_porB_struct <- "Err"
    lw_porB_G120 <- "Err"
    lw_porB_A121 <- "Err"
    }else
    {
      lw_porB_parts <- unlist(strsplit(lw_porB, "/"))
      if (lw_porB_parts[1] == "porB1a")
      {
        lw_porB_struct <- lw_porB_parts[1]
        lw_porB_G120 <- NA
        lw_porB_A121 <- NA
      }else
      {
        lw_porB_struct <- "porB1b"
        lw_porB_G120 <- lw_porB_parts[1]
        lw_porB_A121 <- lw_porB_parts[2]
      }
    }

  lw_porB_prof <- NA

  if (lw_porB_struct == "porB1a")
  {
    lw_porB_prof <- "porB1a"
    v_porB_G120 <- 0L
    v_porB_A121 <- 0L
  }else
  {
    if (lw_porB_G120 != "WT")
    {
      if (is.na(lw_porB_prof))
      {lw_porB_prof <- paste("porB ", lw_porB_G120, sep = "")} else
      {lw_porB_prof <- paste(lw_porB_prof, "/", lw_porB_G120, sep = "")}

      v_porB_G120 <- 1L
    }else {v_porB_G120 <- 0L}
    if (lw_porB_A121 != "WT")
    {
      v_porB_A121 <- 1L
      if (is.na(lw_porB_prof))
      {lw_porB_prof <- paste("porB ", lw_porB_A121, sep = "")} else
      {lw_porB_prof <- paste(lw_porB_prof, "/", lw_porB_A121, sep = "")}
    }else {v_porB_A121 <- 0L}
  }
  if (!is.na(lw_porB_prof))
  {
    if (is.na(molec_profile))
    {molec_profile <- lw_porB_prof
    } else
    {molec_profile <- paste(molec_profile, sepr, lw_porB_prof, sep = "")
    }
  }

  #--------------------------------------------------------------------------------------ponA

  lw_ponA <-   as.character(Combined_Output.df[m, "ponA"])
  if (lw_ponA == "?" | lw_ponA == 'x') {lw_ponA <- "Err"}

  if (lw_ponA != "WT")
  {
    if (is.na(molec_profile))
    {molec_profile <- paste("ponA ", lw_ponA, sep = "")} else
    {molec_profile <- paste(molec_profile, sepr, "ponA ", lw_ponA, sep = "")}

    v_ponA_L421P <- 1L
  } else
  {
    v_ponA_L421P <- 0L
  }

  #--------------------------------------------------------------------------------------gyrA
  lw_gyrA <-   as.character(Combined_Output.df[m, "gyrA"])
  if (lw_gyrA == "?" | lw_gyrA == "x") {lw_gyrA <- "Err/Err"}
  lw_gyrA_parts <- unlist(strsplit(lw_gyrA, "/"))
  lw_gyrA_S91 <-  lw_gyrA_parts[1]
  lw_gyrA_D95 <- lw_gyrA_parts[2]

  lw_gyrA_prof <- NA

    if (lw_gyrA_S91 != "WT")
    {
      v_gyrA_S91 <- 1L
      if (is.na(lw_gyrA_prof))
      {lw_gyrA_prof <- paste("gyrA ", lw_gyrA_S91, sep = "")} else
      {lw_gyrA_prof <- paste(lw_gyrA_prof, "/", lw_gyrA_S91, sep = "")}
    }else {v_gyrA_S91 <- 0L}
  if (lw_gyrA_D95 != "WT")
  {
    v_gyrA_D95 <- 1L
    if (is.na(lw_gyrA_prof))
    {lw_gyrA_prof <- paste("gyrA ", lw_gyrA_D95, sep = "")} else
    {lw_gyrA_prof <- paste(lw_gyrA_prof, "/", lw_gyrA_D95, sep = "")}
  }else {v_gyrA_D95 <- 0L}

    if (!is.na(lw_gyrA_prof))
    {
      if (is.na(molec_profile))
      {molec_profile <- lw_gyrA_prof} else
      {molec_profile <- paste(molec_profile, sepr, lw_gyrA_prof, sep = "")}
    }


  #--------------------------------------------------------------------------------------parC
  lw_parC <-   as.character(Combined_Output.df[m, "parC"])
  if (lw_parC == "?" | lw_parC == "x") {lw_parC <- "Err/Err/Err/Err"}

  lw_parC_parts <- unlist(strsplit(lw_parC, "/"))
  lw_parC_D86 <-  lw_parC_parts[1]
  lw_parC_S87 <- lw_parC_parts[2]
  lw_parC_S88 <- lw_parC_parts[3]
  lw_parC_E91 <- lw_parC_parts[4]


  lw_parC_prof <- NA

  if (lw_parC_D86 != "WT")
  {
    v_parC_D86 <- 1L
    if (is.na(lw_parC_prof))
    {lw_parC_prof <- paste("parC ", lw_parC_D86)} else
    {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_D86, sep = "")}
  }else {v_parC_D86 <- 0L}


  if (lw_parC_S87 != "WT")
  {
    v_parC_S87 <- 1L
    if (is.na(lw_parC_prof))
    {lw_parC_prof <- paste("parC ", lw_parC_S87)} else
    {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_S87, sep = "")}
  }else {v_parC_S87 <- 0L}

  if (lw_parC_S87 == "S87R") {v_parC_S87R <- 1L} else {v_parC_S87R <- 0L}
  if (lw_parC_S87 == "S87I") {v_parC_S87I <- 1L} else {v_parC_S87I <- 0L}
  if (lw_parC_S87 == "S87C") {v_parC_S87C <- 1L} else {v_parC_S87C <- 0L}
  if (lw_parC_S87 == "S87N") {v_parC_S87N <- 1L} else {v_parC_S87N <- 0L}


  if (lw_parC_S88 != "WT")
  {
    V_parC_S88 <- 1L
    if (is.na(lw_parC_prof))
    {lw_parC_prof <- paste("parC ", lw_parC_S88)} else
    {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_S88, sep = "")}
  }else {V_parC_S88 <- 0L}
  if (lw_parC_E91 != "WT")
  {
    V_parC_E91 <- 1L
    if (is.na(lw_parC_prof))
    {lw_parC_prof <- paste("parC ", lw_parC_E91)} else
    {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_E91, sep = "")}
  }else {V_parC_E91 <- 0L}


  if (!is.na(lw_parC_prof))
  {
    if (is.na(molec_profile))
    {molec_profile <- lw_parC_prof} else
    {molec_profile <- paste(molec_profile, sepr, lw_parC_prof, sep = "")}
  }

  lw_23S_A2059G <- as.integer(Combined_Output.df[m, "A2059G"])
  lw_23S_C2611T <- as.integer(Combined_Output.df[m, "C2611T"])

  lw_23S_prof <- NA

  if (is.na(lw_23S_A2059G)) {lw_23S_A2059G <- "Err"}
  if (is.na(lw_23S_C2611T)) {lw_23S_C2611T <- "Err"}

  if (lw_23S_A2059G != 0L)
  {
    if (is.na(lw_23S_prof))
    {lw_23S_prof <- paste("23S rRNA A2059G ", lw_23S_A2059G, "/4", sep = "")} else
    {lw_23S_prof <- paste(ls_23S_prof, "/A20159G ", lw_23S_A2059G, "/4",  sep = "")}
  }

  if (lw_23S_C2611T != 0L)
  {
    if (is.na(lw_23S_prof))
    {lw_23S_prof <- paste("23S rRNA C2611T ", lw_23S_C2611T, "/4", sep = "")} else
    {lw_23S_prof <- paste(lw_23S_prof, "/C2611T ", lw_23S_C2611T, "/4",  sep = "")}
  }

  if (!is.na(lw_23S_prof))
  {
    if (is.na(molec_profile))
    {molec_profile <- lw_23S_prof} else
    {molec_profile <- paste(molec_profile, sepr, lw_23S_prof, sep = "")}
  }

  #----------------------------------------------------------------------------------------------------
  #                             Assign Resistance to AZI, CEPH, CIP and TET, later assign MIC equation
  #----------------------------------------------------------------------------------------------------
  if (is.na(molec_profile))
  {
    molec_profile <- ""
  }

  amr_profile <- "Susceptible"
  sepr2 <- "/"

  #---------------------------------------------------------------------------------------------- AZI
  azi <- "Susceptible"

  if (str_detect(molec_profile, paste(c("mtrR Err", "A2059 Err", "C2611T Err"),collapse = '|')))
  {
    azi <- "Error"
    azi_interp <- "Error"
    amr_profile <- "Error"
  }else

  {
    azi_MIC_inc <- round(-3.014+
                           (2.596*lw_23S_A2059G)+
                           (1.313*lw_23S_C2611T)+
                           (2.893*v_mtrRMEN)+
                           (5.014*v_mtrRDisrupted)+
                           (0.443*v_mtrR35Adel)+
                           (2.825*v_ermB)+
                           (4.038*v_ermC)+
                           (0.840*v_porB_G120)
    )

    if (azi_MIC_inc < -3L) {azi_MIC_inc <- -3L}

    if (azi_MIC_inc > 9L) {azi_MIC_inc <- 9L}

    azi_MIC <- round((2 ^ azi_MIC_inc), digits = 3)
    azi <- paste(azi_MIC, " ug/ml", sep = "")
    azi_interp <- azi_mic.df$Interp[azi_mic.df$MIC == azi_MIC]
    if (azi_interp != "S" )
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- paste("AZI-", azi_interp, sep = ""  )} else
      {amr_profile <- paste(amr_profile, sepr2, "AZI-", azi_interp, sep = "")}
    }
  }

  #---------------------------------------------------------------------------------------------- penicillin
  pen <- NA

  if (str_detect(lw_penA_prof, "Err"))
  {
   pen <- "Error"
   pen_interp <- "Error"
   amr_profile <- "Error"
  }else
  {
  pen_MIC_inc <- round(-3.21+
                         (6.42*v_bla)+
                         (0.56*v_mtrRMENDIS)+
                         (1.34*v_porB_G120)+
                         (1.55*v_ponA_L421P)+
                         (1.25*v_penA_N513Y)+
                         (1.28*v_penA_A517G)+
                         (0.42*v_penA_G543S)
                      )
  if (pen_MIC_inc < -3L) {cfm_MIC_inc <- -3L}
  if (pen_MIC_inc > 7L) {cfm_MIC_inc <- 7L}
  pen_MIC <- 2 ^ pen_MIC_inc

  pen <- paste(pen_MIC, " ug/ml", sep = "")
  pen_interp <- pen_mic.df$Interp[pen_mic.df$MIC == pen_MIC]

  if (pen_interp != "S" )
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- paste("PEN-", pen_interp, sep = ""  )} else
    {amr_profile <- paste(amr_profile, sepr2, "PEN-", pen_interp, sep = "")}
  }

  }

  #---------------------------------------------------------------------------------------------- cephalosporins

  cfx <- NA
  cfm <- NA

  if (str_detect(lw_penA_prof, "Err"))
  {
    cfx_interp <- "Error"
    cfx <- "Error"
    cfm_interp <- "Error"
    cfm <- "Error"
    #amr_profile <- "Error"
  }else
  {

  cfx_MIC_inc <- round(-7.72+
                         (0.54*v_mtrRMENDIS)+
                         (1.38*v_porB_G120)+
                         (0.67*v_ponA_L421P)+
                         (3.90*v_penA_A311V)+
                         (5.15*v_penA_A501P)+
                         (1.51*v_penA_A501T)+
                         (1.92*v_penA_A501V)+
                         (1.53*v_penA_N513Y)+
                         (0.43*v_penA_A517G)+
                         (0.48*v_penA_G543S)

                       )
  if (cfx_MIC_inc < -8L) {cfm_MIC_inc <- -8L}
  if (cfx_MIC_inc > 1L) {cfm_MIC_inc <- 1L}
  cfx_MIC <- round((2 ^ cfx_MIC_inc), digits = 3)

  cfx <- paste(cfx_MIC, " ug/ml", sep = "")
  cfx_interp <- cfx_mic.df$Interp[cfx_mic.df$MIC == cfx_MIC]

  cfm_MIC_inc <- round(-7.198+
                         (0.386*v_mtrRMENDIS)+
                         (0.932*v_mtrR_G45D)+
                         (0.407*v_porB_G120)+
                         (5.422*v_penA_A311V)+
                         (4.494*v_penA_A501P)+
                         (1.463*v_penA_A501T)+
                         (1.185*v_penA_A501V)+
                         (4.297*v_penA_N513Y)+
                         (0.497*v_penA_A517G)
                       )


  if (cfm_MIC_inc < -7L) {cfm_MIC_inc <- -7L}

  if (cfm_MIC_inc > 2L) {cfm_MIC_inc <- 2L}

  cfm_MIC <- round((2 ^ cfm_MIC_inc), digits = 3)

  cfm <- paste(cfm_MIC, " ug/ml", sep = "")
  cfm_interp <- cfm_mic.df$Interp[cfm_mic.df$MIC == cfm_MIC]
  }

  if (cfm_interp != "S" )
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- paste("CE-", cfm_interp, sep = ""  )} else
    {amr_profile <- paste(amr_profile, sepr2, "CE-", cfm_interp, sep = "")}
  }

  if (cfx_interp != "S" )
  {
    if (amr_profile == "Susceptible")
    {amr_profile <- paste("CX-", cfx_interp, sep = ""  )} else
    {amr_profile <- paste(amr_profile, sepr2, "CX-", cfx_interp, sep = "")}
  }

  #---------------------------------------------------------------------------------------------- fluorquilones
  if (is.na(lw_gyrA_prof)) {lw_gyrA_prof<-""}
  if (is.na(lw_parC_prof)) {lw_parC_prof<-""}

  cip <- NA

  if (str_detect(lw_gyrA_prof, "Err") | str_detect(lw_parC_prof, "Err"))
  {
    cip_interp <- "Error"
    cip <- "Error"
  }else
  {

  cip_MIC_inc <- round(-7.62+
                         (5.67*v_gyrA_S91)+
                         (5.05*v_parC_D86)+
                         (5.67*v_parC_S87R)+
                         (4.18*v_parC_S87I)+
                         (5.95*v_parC_S87C)+
                         (1.79*v_parC_S87N)+
                         (1.45*V_parC_S88)+
                         (5.43*V_parC_E91)
                      )
  if (cip_MIC_inc < -8L) {cip_MIC_inc <- -8L}  #MIC range 0.004 - 64 ug/ml
  if (cip_MIC_inc > 6L) {cip_MIC_inc <- 6L}

  cip_MIC <- round( (2^cip_MIC_inc), digits = 3)

  cip <- paste(cip_MIC, " ug/ml", sep = "")
  cip_interp <- cip_mic.df$Interp[cip_mic.df$MIC == cip_MIC]


if (cip_interp != "S" )
{
  if (amr_profile == "Susceptible")
  {amr_profile <- paste("CIP-", cip_interp, sep = ""  )} else
  {amr_profile <- paste(amr_profile, sepr2, "CIP-", cip_interp, sep = "")}
}

  }

  #---------------------------------------------------------------------------------------------- tetracycline

    tet <- "Susceptible"

    if (lw_rpsJ == "Err")
    {
      tet <- "Error"
    }
    if (lw_rpsJ_blast == "NEG")
    {
      lw_rpsJ <- "NEG"

      if (lw_tetM == "NEG")
      {
        tet <- "Error"
      } else                   # if TET-M is positive (big MIC), over-ride the rpsJ error
      {
        tet <- ""
        v_rpsJ <- 0L
      }
    }

    if (tet != "Error")
    {
    tet_MIC_inc <- round(-1.83+
                           (0.62*v_mtrRANY)+
                           (0.26*v_mtrR_A39T)+
                           (0.79*v_porB_G120)+
                           (0.22*v_porB_A121)+
                           (2.11*v_rpsJ)+
                           (4.15*v_tetM)
    )
    if (tet_MIC_inc < -3L) {tet_MIC_inc <- -3L}  #MIC range 0.125 - 128
    if (tet_MIC_inc > 7L) {tet_MIC_inc <- 7L}

    tet_MIC <- round( (2^tet_MIC_inc), digits = 3)

    tet <- paste(tet_MIC, " ug/ml", sep = "")
    tet_interp <- tet_mic.df$Interp[tet_mic.df$MIC == tet_MIC]


    if (tet_interp != "S" )
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- paste("TET-", tet_interp, sep = ""  )} else
      {amr_profile <- paste(amr_profile, sepr2, "TET-", tet_interp, sep = "")}
    }

  }

  #-----------------------------------------------------------------------------------  spectinomycin

  spe <- "<= 32 ug/ml"
  if (str_detect(molec_profile, "16S rRNA Err/Err"))
  {spe <- "Error"

  if (amr_profile == "Susceptible")
  {amr_profile <- "SPE-Err"} else
  {amr_profile <- paste(amr_profile, sepr2, "SPE-Err", sep = "")}

  }else
  {

    if (lw_16S_C1192T == "C1192T")
    {
      spe <- ">=128 ug/ml"
    }
    if (lw_16S_T1458C == "T1485C")
    {
      spe <- "64 ug/ml"
    }

    if (spe == ">=128 ug/ml" | spe == "64 ug/ml")
    {
      if (amr_profile == "Susceptible")
      {amr_profile <- "SPE-R"} else
      {amr_profile <- paste(amr_profile, sepr2, "SPE-R", sep = "")}
    }

  }

  #--------------------------------------------------------------------------------------

  sample_data.df <- tibble(lw_CurrSampleNo, lw_ermB, lw_ermC, lw_rpsJ, lw_tetM, lw_bla, lw_penA_prof,
                               lw_mtrR_p, lw_mtrR_A39, lw_mtrR_G45,
                               lw_porB_struct, lw_porB_G120, lw_porB_A121, lw_ponA,
                               lw_gyrA_S91, lw_gyrA_D95,
                               lw_parC_D86, lw_parC_S87, lw_parC_S88, lw_parC_E91,
                               lw_23S_A2059G, lw_23S_C2611T, lw_16S_C1192T, lw_16S_T1458C, amr_profile,
                               azi, cfx, cfm, cip, tet, pen, spe)

  if(m==1)
  {
    lw_Output.df <- tibble(sample_data.df)
  }else
  {
    lw_Output.df <- bind_rows(lw_Output.df, sample_data.df)
  }
}

lw_output_bad.df <- filter(lw_Output.df, amr_profile == "Error")
lw_output_good.df <- filter(lw_Output.df, amr_profile != "Error")

write.csv(lw_Output.df, paste(local_output_dir, "LabWareUpload_GONO_AMR.csv", sep = ""), quote = FALSE, row.names = FALSE)
write.csv(lw_output_good.df, paste(local_output_dir, "LabWareUpload_GONO_AMR_good.csv", sep = ""), quote = FALSE, row.names = FALSE)
write.csv(lw_output_bad.df, paste(local_output_dir, "LabWareUpload_GONO_AMR_bad.csv", sep = ""), quote = FALSE, row.names = FALSE)

return(lw_Output.df)

}
