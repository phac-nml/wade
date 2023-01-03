##################################################################
####  WGS Analysis and Detection of Molecular Markers (WADE)  ####
####                 Author: Walter Demczuk                   ####
####                    Date: 2022-11-23                      ####
##################################################################

library(plyr)
library(dplyr)
library(tidyverse)
library(tidyselect)
library(stringr)
library(Biostrings)
library(shiny)
library(DT)
library(here)
library(WADE)  #installed from github

#USER INPUT: set location of working directory where system, lookup, mapping, init files
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
curr_work_dir <- "C:\\WADE\\"
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Define UI ---------------------------------------------------------------------------------------
ui <- fluidPage(

  img(src = "Logo.png"),

  sidebarLayout(position = "left",
    sidebarPanel(

      selectInput("Org",
                  label = h4("Choose an Organism"),
                  choices = list("GAS",
                                 "GBS",
                                 "PNEUMO",
                                 "GONO"),
                  selected = "GAS"),

      conditionalPanel(condition = "input.Org == 'GONO'",
                       radioButtons("test", h4("Choose an analysis:"),
                       choices = list("NG-STAR Only" = "NGSTAR",
                                      "AMR Profile" = "AMR_ALL",
                                      "NG-MAST" = "NGMAST",
                                      "MLST" = "MLST"),
                       selected = "AMR")),

      conditionalPanel(condition = "input.Org == 'GAS'",
                       radioButtons("test2", h4("Choose an analysis:"),
                       choices = list("AMR Profile" = "AMR",
                                      "Spe Toxin Profile" = "TOXINS",
                                      "MLST" = "MLST",
                                      "emm Typing" = "EMM"),
                       selected = "EMM")),

      conditionalPanel(condition = "input.Org == 'GBS'",
                       radioButtons("test4", h4("Choose an analysis:"),
                       choices = list("AMR profile" = "AMR",
                                      "SeroTypR" = "SERO2",
                                      "MLST Type" = "MLST"),
                       selected = "SERO2")),

      conditionalPanel(condition = "input.Org == 'PNEUMO'",
                       radioButtons("test3", h4("Choose an analysis:"),
                       choices = list("AMR Profile" = "AMR_ALL",
                                      "SeroTypR" = "SERO",
                                      "Virulence Factors" = "VIRULENCE"),
                       selected = "SERO")),

      conditionalPanel(condition = "input.Org != 'Curator'",
                       textInput("locus", h4("Enter a locus to query or \"list\" for default loci list"),
                       value = "list")),

      conditionalPanel(condition = "input.Org != 'Curator'",
                       textInput("sample", h4("Enter sample number or \"list\" for multiple samples"),
                       value = "list")),

      actionButton("action", "Go"),
      actionButton("openxls", "Output"),
      actionButton("index", "MakeBlastdb")

    ),

    mainPanel(
      titlePanel(strong("Strep/STI WGS Analysis and Detection of Molecular Markers (WADE)")),
                 DT::dataTableOutput("profile_table")
    )
  )
)

# Define server logic -------------------------------------------------------------------------------
server <- function(input, output) {

  observeEvent(input$action,
  {
    if(input$Org == "GONO")
    {
      switch(input$test,
             MLST={output.df <- MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
             NGSTAR={output.df <- NGSTAR_MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
             NGMAST={output.df <- NGMAST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
             rRNA23S={output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)},
             AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             AMR_LW={output.df <- labware_gono_amr(input$Org, curr_work_dir)},
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             AMR_ALL={output.df <- MASTER_pipeline(input$Org, input$test, input$sample, "list", curr_work_dir)
                      output.df <- NGSTAR_MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)
                      output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)
                      output.df <- labware_gono_amr(input$Org, curr_work_dir)},
             {output.df <- MASTER_pipeline(input$Org, input$test, input$sample, input$locus, curr_work_dir)})
    }

    if(input$Org == "GAS")
    {
      switch(input$test2,
             MLST={output.df <- MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
             AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             EMM=(output.df <- EMM_pipeline(input$Org, input$sample, input$locus, curr_work_dir)),
             VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             AMR={output.df <- MASTER_pipeline(input$Org, input$test2, input$sample, input$locus, curr_work_dir)
                  output.df <- labware_gas_amr(input$Org, curr_work_dir)},
             rRNA_16S={output.df <- MASTER_pipeline(input$Org, input$test2, input$sample, input$locus, curr_work_dir)
                       output.df <- rRNA16S_pipeline(input$Org, curr_work_dir)},
             TOXINS={output.df <- MASTER_pipeline(input$Org, input$test2, input$sample, input$locus, curr_work_dir)
                     output.df <- labware_gas_toxins(input$Org, curr_work_dir)},
             {output.df <- MASTER_pipeline(input$Org, input$test2, input$sample, input$locus, curr_work_dir)})
    }

    if (input$Org == "GBS")
    {
      switch(input$test4,
             MLST={output.df <- MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
             AMR={output.df <- MASTER_pipeline(input$Org, input$test4, input$sample, input$locus, curr_work_dir)
                  output.df <- labware_gbs_amr(input$Org, curr_work_dir)},
             rRNA_16S={output.df <- MASTER_pipeline(input$Org, input$test4, input$sample, input$locus, curr_work_dir)
                       output.df <- rRNA16S_pipeline(input$Org, curr_work_dir)},
             AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             SERO=(output.df <- SEROTYPE_pipeline(input$Org, input$sample, input$locus, input$test4, curr_work_dir)),
             SERO2=(output.df <- SEROTYPE_pipeline(input$Org, input$sample, input$locus, input$test4, curr_work_dir)),
             {output.df <- MASTER_pipeline(input$Org, input$test4, input$sample, input$locus, curr_work_dir)})
    }

    if (input$Org == "PNEUMO")
    {
      switch(input$test3,
             AMR_ALL={output.df <- MASTER_pipeline(input$Org, input$test3, input$sample, "list", curr_work_dir)
                      output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)
                      output.df <- labware_pneumo_amr(input$Org, curr_work_dir)},
             AMR={output.df <- MASTER_pipeline(input$Org, input$test3, input$sample, input$locus, curr_work_dir)},
             rRNA_16S={output.df <- MASTER_pipeline(input$Org, input$test3, input$sample, input$locus, curr_work_dir)
                       output.df <- rRNA16S_pipeline(input$Org, curr_work_dir)},
             MLST={output.df <- MLST_pipeline(input$Org, input$sample, input$locus, curr_work_dir)},
             rRNA23S={output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)},
             AMR_LW={output.df <- labware_pneumo_amr(input$Org, curr_work_dir)},
             AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             SERO=(output.df <- SEROTYPE_pipeline(input$Org, input$sample, input$locus, input$test3, curr_work_dir)),
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             VIRULENCE={output.df <- MASTER_pipeline(input$Org, input$test3, input$sample, input$locus, curr_work_dir)
                        output.df <- labware_pneumo_virulence(input$Org, curr_work_dir)},
             {output.df <- MASTER_pipeline(input$Org, input$test3, input$sample, input$locus, curr_work_dir)})
    }

    dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
    Directories.df <- read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    Directories_org.df <- filter(Directories.df, OrgID == input$Org)
    local_dir <- Directories_org.df$LocalDir

    unlink(paste(local_dir, "\\Output\\output_profile.csv", sep = ""))
    write.csv(output.df, paste(local_dir, "\\Output\\output_profile.csv", sep = ""), row.names = F)

    output$profile_table <- renderDataTable({output.df})

  })

  observeEvent(input$openxls,
  {
    dir_file <- paste(curr_work_dir, "DirectoryLocations.csv", sep="")
    Directories.df <- read.csv(dir_file, header = TRUE, sep = ",", stringsAsFactors = FALSE)
    Directories_org.df <- filter(Directories.df, OrgID == input$Org)
    local_dir <- Directories_org.df$LocalDir

    if(input$locus == "list")
    {
      shell.exec(paste(local_dir, "\\Output\\output_profile.csv", sep = ""))
    } else
    {
      if(input$Org == "Curator" & input$test5 == "REMOVE_DUPLICATES")
      {
        shell.exec(paste(local_dir, "\\Output\\output_dna_notfound_distinct.fasta", sep = ""))
      } else
      {
        shell.exec(paste(local_dir, "\\Output\\output_dna.fasta", sep = ""))
      }
    }
  })

  observeEvent(input$index,
  {
    switch(input$Org,
           GONO = {runblast <- Index_pipeline(input$Org, input$test, input$locus, curr_work_dir)},
           GAS = {runblast <- Index_pipeline(input$Org, input$test2, input$locus, curr_work_dir)},
           GBS = {runblast <- Index_pipeline(input$Org, input$test4, input$locus, curr_work_dir)},
           PNEUMO = {runblast <- Index_pipeline(input$Org, input$test3, input$locus, curr_work_dir)})
  })
}

# Run the app -------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)