#' Duplicate Samples Finder
#' February 5 2024, Walter Demczuk & Shelley Peterson
#'
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @details Return duplicate samples
#' Will prompt for a Labware LineList Export .xslx file
#' Scans for duplicate samples from the same patient 
#' based on prov-age-sex-DOB-organism-results
#' Generates Duplicates.csv file which will show whether or not samples have
#' been assigned linked sample numbers.
#' @return A table frame containing the results of the query
#' @export

find_sample_duplicates <- function(Org_id, curr_work_dir) {

  library(readxl)
  library(writexl)
  library(zoo)
  library(lubridate)
  
  # load file
  FileName <- file.choose() #Labware Linelist Export
  data <- read_excel(FileName, guess_max = 10000)
  
  if(Org_id == "GONO")
  {
    # remove minor age differences
    data$age <- floor(data$Age_at_Isolation_Date)
    
    # Consolidate phenotypic and genotypic results
    data$NGMAST <- ifelse(is.na(data$Sequence_Type), data$`WGS_NG-MAST_Type`, data$Sequence_Type)
    
    # Make table with ID code to identify duplicates: Prov-Age-Gender-DOB-Org-Type
    samples <- tibble(NMLno = data$`NML_No`,
                      SenderNo = data$`Sender_No`,
                      LinkedNo = data$`Linked_to_Sample`,
                      SubmLab = data$`Center`,
                      Province = data$Province,
                      IsolnSite = data$`Isolation_Site`,
                      DateCollected = data$`Date_Isolated`,
                      ID = paste(data$Province,
                                 data$Center,
                                 as.character(data$age),
                                 data$Sex,
                                 as.character(data$`Date_of_Birth`),
                                 data$Organism,
                                 data$NGMAST,
                                 data$`WGS_NG-STAR_Type`,
                                 data$WGS_MLST_Type,
                                 sep = "-" ))
  }
  if(Org_id == "STREP")
  {
    # Remove minor age differences
    data$Age <- floor(data$Age)
    
    # Consolidate phenotypic and genotypic results
    data$emmType <- ifelse(is.na(data$'WGS emm Type') & 
                           data$Organism == "Streptococcus pyogenes (GAS)",
                           as.character(data$'emm Type'),
                           substr(as.character(data$'WGS emm Type'), 4, 5))
    
    data$seroSPN <- ifelse(is.na(data$'WGS SPN Serotype') & 
                           data$Organism == "Streptococcus pneumoniae",
                           as.character(data$Serotype),
                           as.character(data$'WGS SPN Serotype'))
    
    data$seroGBS <- ifelse(is.na(data$'WGS Predicted serotype') & 
                           data$Organism == "Streptococcus agalactiae (GBS)",
                           as.character(data$Serotype),
                           as.character(data$'WGS Predicted serotype'))
    
    # Get a single type either emm type, GBS serotpe or SPN serotype (or NA)
    data$Type <- paste0(data$emmType, data$seroSPN, data$seroGBS)
    
    # Make table with ID code to identify duplicates: Prov-Age-Gender-DOB-Org-Type
    samples <- tibble(NMLno = data$`NML Number`, 
                      SubmLabNo = data$`Submitting Lab Number`,
                      LinkedNo = data$`Linked Sample`,
                      SubmLab = data$`Submitting Lab`,
                      ClinSource = data$`Clinical Source`,
                      DateCollected = data$`Date Collected`,
                      ID = paste(data$Province,
                                 as.character(data$Age),
                                 data$Gender,
                                 as.character(data$`Date of Birth`),
                                 data$Organism,
                                 data$Type,
                                 sep = "-" ))
  }
  
  samples <- arrange(samples, ID)
  Duplicates <- samples[duplicated(samples$ID)|duplicated(samples$ID, fromLast=TRUE),]

  # Group by ID then only include duplicate sets that have a collect date <= 30 days of each other
  dates <- Duplicates %>% group_by(ID) %>%
    summarise("diff" = difftime(max(DateCollected), min(DateCollected), units = "days")) %>%
    filter(diff < 31)

  Duplicates <- filter(Duplicates, ID %in% dates$ID)

  write.csv(Duplicates, "Duplicates.csv", na = "", row.names = F)

  cat("\n\nDuplicates.csv table ready!\n\n")

  # open upload file
  shell.exec("Duplicates.csv")
}





