#' Detects potententially contaminated samples based on plate proximity
#' February 5, 2024, Shelley Peterson
#' Run this analysis to evaluate consistency of MIC results
#'
#'
#' @param curr_work_dir Start up directory from pipeline project to locate system file structure
#' @return A table frame containing the results of the query
#' @export

#-------------------------------------------------------------------------------
#  For troubleshooting and debugging
#curr_work_dir <- "C:\\WADE\\"
#-------------------------------------------------------------------------------

MIC_compare <- function(curr_work_dir) {

  library(naniar)
  
  #-----------------------------------------------------------------------------
  # get directory structure
  directorylist <- getdirectory(curr_work_dir, "GONO", "AMR")
  #-----------------------------------------------------------------------------

  # Load datafile
  file <- file.choose()
  LineList <- as_tibble(read.csv(file, header = TRUE, sep = ",", stringsAsFactors = FALSE))
  
  # Select required columns
  # This makes troubleshooting easier
  antimicrobial <- c("Penicillin", "Tetracycline", "Ceftriaxone", "Cefixime", "Ciprofloxacin", "Azithromycin")
  data <- select(LineList, NML_No, Penicillin, Submitted_Penicillin, `WGS_Predicted_Penicillin_MIC_.ug.ml.`,
                 #Spectinomycin, Submitted_Spectinomycin, `WGS_Predicted_Spectinomycin_MIC_.ug.ml.`,
                 Tetracycline, Submitted_Tetracycline, `WGS_Predicted_Tetracycline_MIC_.ug.ml.`,
                 Ceftriaxone, Submitted_Ceftriaxone, `WGS_Predicted_Ceftriaxone_MIC_.ug.ml.`,
                 Cefixime, Submitted_Cefixime, `WGS_Predicted_Cefixime_MIC_.ug.ml.`,
                 Ciprofloxacin, Submitted_Ciprofloxacin, `WGS_Predicted_Ciprofloxacin_MIC_.ug.ml.`,
                 Azithromycin, Submitted_Azithromycin, `WGS_Predicted_Azithromycin_MIC_.ug.ml.`)
  
  # Clean Data
  #baddata <- c("contaminated", "Inconclusive", "negative", "Negative", "no growth", 
  #             "N/A", "NA", "NAN/A", "NotDone", "not done", "Not Done", "UNKNOWN", 
  #             "No Result", "Noresult", "Error", "Err", "Cannot Test", "")
  
  data <- data.frame(apply(data,2, function(x) gsub(" ug/ml|<=|< =|>=|> =|<|>|<= |>= 
                                                    |Susceptible|Resistant|Intermediate|
                                                    Notavailable|Not available|N/A|R|S| |", "", x)))
  #data <- data %>% replace_with_na_all(condition = ~.x %in% baddata)
  data <- data %>% mutate_at(c(2:19), as.numeric)
  
  # Find differences between values
  MICdiffs <- function(data, antimicrobial){
    submitted <- paste0("Submitted_", antimicrobial)
    predicted <- paste("WGS_Predicted", antimicrobial, "MIC_.ug.ml.", sep = "_")
    NML_No <- data$NML_No
    cutoff <- switch(antimicrobial, "Penicillin" = 2,
                                    "Tetracycline" = 2,
                                    "Ceftriaxone" = 0.125, 
                                    "Cefixime" = 0.25, 
                                    "Ciprofloxacin" = 1, 
                                    "Azithromycin" = 2)
    df <- data %>% select(all_of(c(antimicrobial, submitted, predicted)))
    df$min <- do.call(pmin, c(df, na.rm = TRUE))
    df$max <- do.call(pmax, c(df, na.rm = TRUE))
    df$diff <- df$max / df$min
    df$breakpoint <- ifelse(df$max >= cutoff & df$min < cutoff, "NotOK", "OK")
    df <- cbind(NML_No, df)
    df_filtered <- df %>% filter(df$diff > 4 | df$breakpoint == "NotOK")
    df_filtered$retest[df_filtered$breakpoint == "NotOK"] <- paste0(antimicrobial, " breakpoint")
    df_filtered$retest[df_filtered$diff > 4] <- paste0(antimicrobial, " MIC diff")
    df_filtered <- select(df_filtered, c("NML_No", "retest"))
    df_filtered <- df_filtered %>% rename("retest" = antimicrobial)
    return(df_filtered)
  }
  
  samples <- lapply(antimicrobial, MICdiffs, data=data)
  names(samples) <- antimicrobial

  samplelist <- merge(samples$Ceftriaxone, samples$Cefixime, all = TRUE)
  samplelist <- merge(samplelist, samples$Azithromycin, all = TRUE)
  samplelist <- merge(samplelist, samples$Penicillin, all = TRUE)
  samplelist <- merge(samplelist, samples$Tetracycline, all = TRUE)
  samplelist <- merge(samplelist, samples$Ciprofloxacin, all = TRUE)

  write.csv(samplelist, paste0(directorylist$output_dir, "MIC_Check.csv"), quote = FALSE, row.names = FALSE, na = "")  

  cat("\n\nDone! ", "MIC_Check.csv is ready in output folder", "\n\n\n", sep = "")
  
  return(samplelist)
}