#' Run AMR first
#' Then run this analysis to combine data the full amr profile to upload to LabWare.
#'
#'
#' @param Org_id Organism to query: GAS, PNEUMO or GONO
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export
#'
#'

labware_gas_amr <- function(Org_id, curr_work_dir) {

  #------------------------------------------------------------------------------------------------------------
  # get directory structure
  dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")

  Directories.df <- as_tibble(read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Directories_org.df <- filter(Directories.df, OrgID == Org_id)
  local_dir <- Directories_org.df$LocalDir
  SampList <- paste(local_dir, "list.csv", sep = "")
  local_output_dir <- paste(local_dir, "Output\\", sep = "")
  local_temp_dir <- paste(local_dir, "temp\\", sep = "")
  system_dir <- Directories_org.df$SystemDir
  ContigsDir <- Directories_org.df$ContigsDir

  Output.df <- as_tibble(read.csv(paste(local_output_dir, "output_profile_GAS_AMR.csv", sep = ""),
                                  header = TRUE, sep = ",", stringsAsFactors = FALSE))
  Size.df <- dim(Output.df)
  NumSamples <- Size.df[1]
  NumLoci <- ((Size.df[2]-3) / 7)

  if (NumLoci > 1)
  {
    sepr <- "; "

    m <- 1

    for (m in 1L:NumSamples)  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    {
      molec_profile <- ""
      lw_comments <- ""
      lw_allele_profile <- ""
      sepr3 <- ":" #for allele profile

      lw_CurrSampleNo <- as.character(Output.df[m, "SampleNo"])
      lw_ermA <- as.character(Output.df[m, "ermA_result"])

      if (lw_ermA != "Sample_Err")
      {

        if (lw_ermA == "POS")
        {
          # put ermA promoter analysis here
          lw_ermAp_result <- as.character(Output.df[m, "ermAp_result"]) #POS/NEG
          lw_ermAp_allele <- as.character(Output.df[m, "ermAp_allele"]) #ID number
          lw_ermAp_mutation <- as.character(Output.df[m, "ermAp_mutations"]) # differences from wildtype
          lw_ermAp_comment <- as.character(Output.df[m, "ermAp_comments"]) #Susceptible(Inducible)/Resistant
          lw_ermAp_motif <- as.character(Output.df[m, "ermAp_motifs"]) #amino acid substitutions from MasterBlastR

          if (lw_ermAp_motif == "WT/WT") {lw_ermAp_motif <- ""}

          if (lw_ermAp_result == "NEG")
          {
            molec_value <- "ermA R"
            lw_ermAp_allele <- " 0"
          } else # susceptible promoter found, might have mutations for CLI-R
            {
              if (lw_ermAp_motif == "N90H/D97G")
              {
                molec_value <- "ermA S"    # CLI-Inducible
                lw_ermAp_allele <- paste(" ", lw_ermAp_allele, sep = "")
              } else
              {
                molec_value <- "ermA R"
                lw_ermAp_allele <-" 0"
              }
            }

          if (molec_profile == "")
          {molec_profile <- paste(molec_value, " ", lw_ermAp_motif, sep = "")}else
          {molec_profile <- paste(molec_profile, sepr, molec_value, " ", lw_ermAp_motif, sep = "")}

          if (lw_allele_profile == "")
          {lw_allele_profile <- paste("ermAp", lw_ermAp_allele, sep = "")}else
          {lw_allele_profile <- paste(lw_allele_profile, sepr3, "ermAp", lw_ermAp_allele, sep = "")}
        }

        lw_ermB <- as.character(Output.df[m, "ermB_result"])
        if (lw_ermB == "POS")
        {
          lw_ermB_allele <- as.character(Output.df[m, "ermB_allele"]) #if gene broken will have ID, else "NF"
          if (lw_ermB_allele == "NF")
            {lw_ermB_mutations <- ""
             lw_ermB_allele <- "ermB"}else
            {lw_ermB_mutations <-  as.character(Output.df[m, "ermB_mutations"])
             lw_ermB_allele <- paste(as.character(Output.df[m, "ermB_allele"]),
                                     " ",
                                     as.character(Output.df[m, "ermB_comments"]),
                                     sep = "")
            }

          if (molec_profile == "")
          {molec_profile <- paste(molec_profile, "ermB ", lw_ermB_mutations, sep = "")}else
          {molec_profile <- paste(molec_profile, sepr, "ermB ", lw_ermB_mutations, sep = "")}

          if (lw_allele_profile == "")
          {lw_allele_profile <- paste(lw_allele_profile, "ermB ", lw_ermB_mutations, sep = "")}else
          {lw_allele_profile <- paste(lw_allele_profile, sepr3, "ermB ", lw_ermB_mutations, sep = "")}

        }

        lw_ermT <- as.character(Output.df[m, "ermT_result"])
        if (lw_ermT == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "ermT"}else
          {molec_profile <- paste(molec_profile, sepr, "ermT", sep = "")}
        }

        lw_mefAE <- as.character(Output.df[m, "mefAE_result"])
        if (lw_mefAE == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "mefAE"}else
          {molec_profile <- paste(molec_profile, sepr, "mefAE", sep = "")}
        }

        lw_msrD <- as.character(Output.df[m, "msrD_result"])
        if (lw_msrD == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "msrD"}else
          {molec_profile <- paste(molec_profile, sepr, "msrD", sep = "")}
        }

        lw_gyrA <-   as.character(Output.df[m, "gyrA_motifs"])  #use motifs to reduce curations errors
        lw_gyrA_result <-   as.character(Output.df[m, "gyrA_result"])
        lw_gyrA_allele <- as.character(Output.df[m, "gyrA_allele"])

        if (lw_gyrA_result == "NEG")
        {
          lw_gyrA <- "no gene"
          lw_comments <- paste(lw_comments, "gyrA missing.", sep = "")
        }

        if (!(lw_gyrA %in% c("WT", "no gene")))
        {
          if (molec_profile == "")
          {molec_profile <- paste("gyrA ", lw_gyrA, sep = "")}else
          {molec_profile <- paste(molec_profile, sepr, "gyrA ", lw_gyrA, sep = "")}
        }

        if (lw_allele_profile == "")
        {lw_allele_profile <- paste(lw_allele_profile, "gyrA ", lw_gyrA_allele, sep = "")}else
        {lw_allele_profile <- paste(lw_allele_profile, sepr3, "gyrA ", lw_gyrA_allele, sep = "")}

        lw_parC_result <-  as.character(Output.df[m, "parC_result"])
        lw_parC_allele <- as.character(Output.df[m, "parC_allele"])
        lw_parC_prof <- ""

        if (lw_parC_result == "NEG")
        {
          lw_parC <- "no gene"
          lw_comments <- paste(lw_comments, "parC missing", sep="")
        } else
        {
          lw_parC <-   as.character(Output.df[m, "parC_motifs"])
          lw_parC_parts <- unlist(strsplit(lw_parC, "/"))
          lw_parC_D78 <-  lw_parC_parts[1]
          lw_parC_S79 <- lw_parC_parts[2]
          lw_parC_D83 <- lw_parC_parts[3]

          if (lw_parC_D78 != "WT")
          {
            if (lw_parC_prof == "")
            {lw_parC_prof <- paste("parC ", lw_parC_D78, sep = "")} else
            {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_D78, sep = "")}
          }
          if (lw_parC_S79 != "WT")
          {
            if (lw_parC_prof == "")
            {lw_parC_prof <- paste("parC ", lw_parC_S79, sep = "")} else
            {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_S79, sep = "")}
          }
          if (lw_parC_D83 != "WT")
          {
            if (lw_parC_prof == "")
            {lw_parC_prof <- paste("parC ", lw_parC_D83, sep = "")} else
            {lw_parC_prof <- paste(lw_parC_prof, "/", lw_parC_D83, sep = "")}
          }
        }

        if (lw_parC_prof != "")
        {
          if (molec_profile == "")
          {molec_profile <- lw_parC_prof} else
          {molec_profile <- paste(molec_profile, sepr, lw_parC_prof, sep = "")}
        }

        if (lw_allele_profile == "")
        {lw_allele_profile <- paste(lw_allele_profile, "parC ", lw_parC_allele, sep = "")}else
        {lw_allele_profile <- paste(lw_allele_profile, sepr3, "parC ", lw_parC_allele, sep = "")}

        lw_tetM <- as.character(Output.df[m, "tetM_result"])
        if (lw_tetM == "POS")
        {
          lw_tetM_comments <- paste(as.character(Output.df[m, "tetM_comments"]), sep = "")

          if(lw_tetM_comments == "NA") {lw_tetM_comments <- ""}

          if (molec_profile == "")
          {molec_profile <- paste("tetM ", lw_tetM_comments, sep = "")}else
          {molec_profile <- paste(molec_profile, sepr, "tetM ", lw_tetM_comments, sep = "")}
        }

        lw_tetO <- as.character(Output.df[m, "tetO_result"])
        if (lw_tetO == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "tetO"}else
          {molec_profile <- paste(molec_profile, sepr, "tetO", sep = "")}
        }

        lw_tetT <- as.character(Output.df[m, "tetT_result"])
        if (lw_tetT == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "tetT"}else
          {molec_profile <- paste(molec_profile, sepr, "tetT", sep = "")}
        }

        lw_tetA <- as.character(Output.df[m, "tetA_result"])
        if (lw_tetA == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "tetA"}else
          {molec_profile <- paste(molec_profile, sepr, "tetA", sep = "")}
        }

        lw_tetB <- as.character(Output.df[m, "tetB_result"])
        if (lw_tetA == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "tetB"}else
          {molec_profile <- paste(molec_profile, sepr, "tetB", sep = "")}
        }

        lw_cat <- as.character(Output.df[m, "cat_result"])
        if (lw_cat == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "cat"}else
          {molec_profile <- paste(molec_profile, sepr, "cat", sep = "")}
        }

        lw_catQ <- as.character(Output.df[m, "catQ_result"])
        if (lw_catQ == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "catQ"}else
          {molec_profile <- paste(molec_profile, sepr, "catQ", sep = "")}
        }

        lw_dfrF <- as.character(Output.df[m, "dfrF_result"])
        if (lw_dfrF == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "dfrF"}else
          {molec_profile <- paste(molec_profile, sepr, "dfrF", sep = "")}
        }

        lw_dfrG <- as.character(Output.df[m, "dfrG_result"])
        if (lw_dfrG == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "dfrG"}else
          {molec_profile <- paste(molec_profile, sepr, "dfrG", sep = "")}
        }


        lw_folA_result <-   as.character(Output.df[m, "folA_result"])
        lw_folA_allele <- as.character(Output.df[m, "folA_allele"])

        if (lw_folA_result == "NEG")
        {
          lw_folA <- "no gene"
          lw_comments <- paste(lw_comments, "folA missing.", sep = "")
        } else
          {
          lw_folA <-   as.character(Output.df[m, "folA_motifs"])
          }

        if (!(lw_folA %in% c("WT", "no gene")))
        {
          if (molec_profile == "")
          {molec_profile <- paste("folA ", lw_folA, sep = "")}else
          {molec_profile <- paste(molec_profile, sepr, "folA ", lw_folA, sep = "")}
        }

        if (lw_allele_profile == "")
        {lw_allele_profile <- paste(lw_allele_profile, "folA ", lw_folA_allele, sep = "")}else
        {lw_allele_profile <- paste(lw_allele_profile, sepr3, "folA ", lw_folA_allele, sep = "")}

        lw_folP_result <-   as.character(Output.df[m, "folP_result"])
        lw_folP_allele <- as.character(Output.df[m, "folP_allele"])

        if (lw_folP_result == "NEG")
        {
          lw_folP <- "no gene"
          lw_comments <- paste(lw_comments, "folP missing.", sep = "")
        } else
        {
          lw_folP <-   as.character(Output.df[m, "folP_motifs"])
        }

        if (!(lw_folP %in% c("WT", "no gene")))
        {
          if (molec_profile == "")
          {molec_profile <- paste("folP ", lw_folP, sep = "")}else
          {molec_profile <- paste(molec_profile, sepr, "folP ", lw_folP, sep = "")}
        }

        if (lw_allele_profile == "")
        {lw_allele_profile <- paste(lw_allele_profile, "folP ", lw_folP_allele, sep = "")}else
        {lw_allele_profile <- paste(lw_allele_profile, sepr3, "folP ", lw_folP_allele, sep = "")}

        lw_pbp2x_result <-   as.character(Output.df[m, "pbp2x_result"])
        lw_pbp2x_allele <- as.character(Output.df[m, "pbp2x_allele"])

        if (lw_pbp2x_result == "NEG")
        {
          lw_pbp2x <- "no gene"
          lw_comments <- paste(lw_comments, "pbp2x missing.", sep = "")
        } else
        {
          lw_pbp2x <-   as.character(Output.df[m, "pbp2x_motifs"])
        }

        if (!(lw_pbp2x %in% c("WT", "no gene")))
        {
          if (molec_profile == "")
          {molec_profile <- paste("pbp2x ", lw_pbp2x, sep = "")}else
          {molec_profile <- paste(molec_profile, sepr, "pbp2x ", lw_pbp2x, sep = "")}
        }

        if (lw_allele_profile == "")
        {lw_allele_profile <- paste(lw_allele_profile, "pbp2x ", lw_pbp2x_allele, sep = "")}else
        {lw_allele_profile <- paste(lw_allele_profile, sepr3, "pbp2x ", lw_pbp2x_allele, sep = "")}

        #-----------------------------------------------------------------  Vancomycin

        lw_vanA <- as.character(Output.df[m, "vanA_result"])
        if (lw_vanA == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "vanA"}else
          {molec_profile <- paste(molec_profile, sepr, "vanA", sep = "")}
        }

        lw_vanB <- as.character(Output.df[m, "vanB_result"])
        if (lw_vanB == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "vanB"}else
          {molec_profile <- paste(molec_profile, sepr, "vanB", sep = "")}
        }

        lw_vanC <- as.character(Output.df[m, "vanC_result"])
        if (lw_vanC == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "vanC"}else
          {molec_profile <- paste(molec_profile, sepr, "vanC", sep = "")}
        }

        lw_vanD <- as.character(Output.df[m, "vanD_result"])
        if (lw_vanD == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "vanD"}else
          {molec_profile <- paste(molec_profile, sepr, "vanD", sep = "")}
        }

        lw_vanE <- as.character(Output.df[m, "vanE_result"])
        if (lw_vanE == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "vanE"}else
          {molec_profile <- paste(molec_profile, sepr, "vanE", sep = "")}
        }

        lw_vanG <- as.character(Output.df[m, "vanG_result"])
        if (lw_vanG == "POS")
        {
          if (molec_profile == "")
          {molec_profile <- "vanG"}else
          {molec_profile <- paste(molec_profile, sepr, "vanG", sep = "")}
        }

        #----------------------------------------------------------------------
        if (molec_profile == "") {molec_profile <- "Wild Type"}
        #----------------------------------------------------------------------  INTERPRETATIONS
        amr_profile <- "Susceptible"
        # ERY ---------------------------------
        if (str_detect(molec_profile, paste(c("ermA", "ermB", "ermT", "mefAE"),collapse = '|')))
        {
          ery_MIC <- ">= 2 ug/ml"
          ery <- "Resistant"
        }else{
          ery_MIC <- "<= 0.25 ug/ml"
          ery <- "Susceptible"}

        # CLI ---------------------------------
        if (str_detect(molec_profile, "ermB"))
        {
          cli_MIC <- ">= 1 ug/ml"
          cli <- "Resistant"
        }else {
          cli_MIC <- "<= 0.12 ug/ml"
          cli <-"Susceptible"}

        if (str_detect(molec_profile, "ermB N100S"))
        {
          cli_MIC <- "<= 0.12 ug/ml"
          cli <- "Susceptible"
        }

        if (str_detect(molec_profile, "ermA S"))
        {
          cli_MIC <- "<= 0.12 ug/ml"
          cli <- "Inducible"
        }
        if (str_detect(molec_profile, "ermA R"))
        {
          cli_MIC <- ">= 1 ug/ml"
          cli <- "Resistant"
        }

        if (str_detect(molec_profile, paste(c("cat", "catQ"),collapse = '|')))
        {
          chl_MIC <- ">= 32 ug/ml"
          chl <- "Resistant"
        }else{
          chl_MIC <- "<= 4 ug/ml"
          chl <- "Susceptible"}

        lev <- "Undetermined"
        lev_MIC <- "Undetermined"
        if (lw_gyrA == "no gene" | lw_parC == "no gene" )
        { lev_MIC <- "Error"
        lev <- "Error"
        }else
        {

          if (str_detect(molec_profile, "gyrA S81")) #gyrA mutations
          {
            lev_MIC <- ">= 8 ug/ml"
            lev <- "Resistant"
          } else
            if ((str_detect(molec_profile, "D78"))) # first parC mutation
            {
              lev_MIC <- "1 ug/ml"
              lev <- "Susceptible"
            } else
              if ((str_detect(molec_profile, paste(c("S79", "D83"),collapse = '|')))) # first parC mutation
              {
                lev_MIC <- "2 ug/ml"
                lev <- "Susceptible"
              } else
              {lev_MIC <- "<= 0.5 ug/ml"
              lev <- "Susceptible"}
        }


        if (str_detect(molec_profile, paste(c("tetM", "tetO", "tetT", "tetA", "tetB"),collapse = '|')))
        {
          tet_MIC <- ">= 8 ug/ml"
          tet <- "Resistant"
        }else
          {
          tet_MIC <- "<= 1 ug/ml"
          tet <- "Susceptible"
          }

        if (str_detect(molec_profile, "tetM Disrupted"))
        {
          tet_MIC <- "<= 1 ug/ml"
          tet <- "Susceptible"
        }

        if (lw_folA == "no gene" | lw_folP == "no gene")
        {
          sxt_MIC <- "Error"
          sxt <- "Error"
          #amr_profile <- "Error"
        }else

        {
          if (str_detect(molec_profile, paste(c("dfrG", "dfrF", "folA", "folP"),collapse = '|')))
          {

            if (str_detect(molec_profile, paste(c("dfrG", "dfrF"),collapse = '|')))
            {
              sxt_MIC <- ">= 4 ug/ml" #=4/76
              sxt <- "Resistant"
            } else if (str_detect(molec_profile, paste(c("folA", "folP"),collapse = '|')))
            {
              sxt_MIC <- "<= 0.5 ug/ml" #=0.5/9.5
              sxt <- "Susceptible"
              if ((str_detect(molec_profile, "folA")) & (str_detect(molec_profile, "folP")))
              {
                sxt_MIC <- ">= 4 ug/ml"
                sxt <- "Resistant"
              }
            }
          }else{
            sxt_MIC <- "<= 0.5 ug/ml"
            sxt <- "Susceptible"}
        }

        if (lw_pbp2x == "no gene")
        {
          pen_MIC <- "Error"
          pen <- "Error"
          #amr_profile <- "Error"
        }else
          {

          if (str_detect(molec_profile, "pbp2x"))
          {
            pen_MIC <- "> 0.12 ug/ml"
            pen <- "Unknown"
          }else
            {
            pen_MIC <- "<= 0.12 ug/ml"
            pen <- "Susceptible"
            }
          }

        # Vancomycin ---------------------------------
        if (str_detect(molec_profile, paste(c("vanA", "vanB", "vanC", "vanD", "vanE", "vanG"),collapse = '|')))
        {
          van_MIC <- ">= 2 ug/ml"
          van <- "Resistant"
        }else{
          van_MIC <- "<= 1 ug/ml"
          van <- "Susceptible"}

        #--------------------------------------------------------------  MAKE AMR PROFILE
        if (amr_profile != "Error")
        {

          amr_profile <- "Susceptible"
          sepr2 <- "/"

          if (ery == "Resistant")
          {
            if (amr_profile == "Susceptible")
            {amr_profile <- "ERY-R"} else
            {amr_profile <- paste(amr_profile, sepr2, "ERY-R", sep = "")}
          }
          if (cli == "Resistant")
          {
            if (amr_profile == "Susceptible")
            {amr_profile <- "CLI-R"} else
            {amr_profile <- paste(amr_profile, sepr2, "CLI-R", sep = "")}
          }
          if (cli == "Inducible")
          {
            if (amr_profile == "Susceptible")
            {amr_profile <- "CLI-Ind"} else
            {amr_profile <- paste(amr_profile, sepr2, "CLI-Ind", sep = "")}
          }

          if (chl == "Resistant")
          {
            if (amr_profile == "Susceptible")
            {amr_profile <- "CHL-R"} else
            {amr_profile <- paste(amr_profile, sepr2, "CHL-R", sep = "")}
          }
          if (lev == "Resistant")
          {
            if (amr_profile == "Susceptible")
            {amr_profile <- "LEV-R"} else
            {amr_profile <- paste(amr_profile, sepr2, "LEV-R", sep = "")}
          }
          if (tet == "Resistant")
          {
            if (amr_profile == "Susceptible")
            {amr_profile <- "TET-R"} else
            {amr_profile <- paste(amr_profile, sepr2, "TET-R", sep = "")}
          }

          if (sxt == "Resistant")
          {
            if (amr_profile == "Susceptible")
            {amr_profile <- "SXT-R"} else
            {amr_profile <- paste(amr_profile, sepr2, "SXT-R", sep = "")}
          }
          if (sxt == "Intermediate")
          {
            if (amr_profile == "Susceptible")
            {amr_profile <- "SXT-I"} else
            {amr_profile <- paste(amr_profile, sepr2, "SXT-I", sep = "")}
          }
          if (pen == "Decreased Susceptiblity")
          {
            if (amr_profile == "Susceptible")
            {amr_profile <- "PEN-DS"} else
            {amr_profile <- paste(amr_profile, sepr2, "PEN-DS", sep = "")}
          }

          if (van == "Resistant")
          {
            if (amr_profile == "Susceptible")
            {amr_profile <- "VAN-R"} else
            {amr_profile <- paste(amr_profile, sepr2, "VAN-R", sep = "")}
          }

        }
      }else #end if lw_ermA != "Sample_Err"
      {
        lw_ermA <- "Sample_Err"
        lw_ermB <- "Sample_Err"
        lw_ermT <- "Sample_Err"
        lw_mefAE <- "Sample_Err"
        lw_gyrA <- "Sample_Err"
        lw_parC <- "Sample_Err"
        lw_tetM <- "Sample_Err"
        lw_tetO <- "Sample_Err"
        lw_tetT <- "Sample_Err"
        lw_tetA <- "Sample_Err"
        lw_tetB <- "Sample_Err"
        lw_cat <- "Sample_Err"
        lw_catQ <- "Sample_Err"
        lw_dfrF <- "Sample_Err"
        lw_dfrG <- "Sample_Err"
        lw_folA <- "Sample_Err"
        lw_folP <- "Sample_Err"
        lw_pbp2x <- "Sample_Err"
        molec_profile <- "Sample_Err"
        ery_MIC <- "Sample_Err"
        ery <- "Sample_Err"
        chl_MIC <- "Sample_Err"
        chl <- "Sample_Err"
        lev_MIC <- "Sample_Err"
        lev <- "Sample_Err"
        cli_MIC <- "Sample_Err"
        cli <- "Sample_Err"
        tet_MIC <- "Sample_Err"
        tet <- "Sample_Err"
        sxt_MIC <- "Sample_Err"
        sxt <- "Sample_Err"
        pen_MIC <- "Sample_Err"
        pen <- "Sample_Err"
        amr_profile <- "Sample_Err"
        lw_allele_profile <- "Sample_Err"
      }

      #-------------------------------------------- New LabWare Uploader; remove #'s when labware ready

      # combine cat and catQ results into single column for labware. Specified will be in molecular profile.

      if (str_detect(molec_profile, paste(c("cat", "catQ"),collapse = '|')))
      {
        lw_cat <- "POS"
      }else{
        lw_cat <- "NEG"}

      #make lw_allele_profile here...
      #lw_ermB_allele <- as.character(Output.df[m, "ermB_allele"])
      # lw_pbp2x_allele <- as.character(Output.df[m, "pbp2x_allele"])
      # lw_folA_allele <- as.character(Output.df[m, "folA_allele"])
      # lw_folP_allele <- as.character(Output.df[m, "folP_allele"])
      # lw_gyrA_allele <- as.character(Output.df[m, "gyrA_allele"])
      # lw_parC_allele <- as.character(Output.df[m, "parC_allele"])
      # lw_allele_profile <- paste(
      #   "ermB", lw_ermB_allele, ":", #calculated above
      #   "gyrA", lw_gyrA_allele, ":",
      #   "parC", lw_parC_allele, ":",
      #   "folA", lw_folA_allele, ":",
      #   "folP", lw_folP_allele, ":",
      #   "pbp2x", lw_pbp2x_allele
      #)

      LabWare_Sample.df <- tibble(
        lw_CurrSampleNo,
        lw_ermA,
        lw_ermB,
        lw_ermT,
        lw_mefAE,
        lw_gyrA,
        lw_parC,
        lw_tetM,
        lw_tetO,
        lw_tetT,
        lw_tetA,
        lw_tetB,
        lw_cat,
        lw_dfrF,
        lw_dfrG,
        lw_folA,
        lw_folP,
        lw_pbp2x,
        molec_profile,
        ery_MIC,
        ery,
        chl_MIC,
        chl,
        lev_MIC,
        lev,
        cli_MIC,
        cli,
        tet_MIC,
        tet,
        sxt_MIC,
        sxt,
        pen_MIC,
        pen,
        amr_profile,
        lw_allele_profile
      )
      #-------------------------------------------------------------------------------------
      if(m==1)  #if first sample make one row profile table, otherwise add new row to table
      {
        LabWare.df <- tibble(LabWare_Sample.df)
      }else
      {
        LabWare.df <- rbind(LabWare.df, LabWare_Sample.df)
      }
    }

    lw_output_bad.df <- filter(LabWare.df, amr_profile == "Error")
    write.csv(LabWare.df, paste(local_output_dir, "LabWareUpload_GAS_AMR.csv", sep = ""), quote = FALSE,  row.names = FALSE)
    write.csv(lw_output_bad.df, paste(local_output_dir, "LabWareUpload_GAS_AMR_bad.csv", sep = ""), quote = FALSE,  row.names = FALSE)

    cat("\n\nDone! ", local_output_dir, "LabWareUpload_GAS_AMR.csv is ready in output folder", "\n\n\n", sep = "")

  }else  # else if a single loci was run return the original output.csv
  {
    LabWare.df <- Output.df
  }

return(LabWare.df)
}
