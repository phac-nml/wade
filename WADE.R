##################################################################
####  WGS Analysis and Detection of Molecular Markers (WADE)  ####
####       Authors: Walter Demczuk & Shelley Peterson         ####
####                    Date: 2025-05-07                      ####
##################################################################

library(plyr)
library(tidyverse)
library(tidyselect)
library(tidysq)
library(Biostrings)
library(shiny)
library(shinyWidgets)
library(DT)
library(readxl)
library(beepr)
library(WADE)  #installed from github

#USER INPUT: set location of working directory where system, lookup, mapping, init files
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
curr_work_dir <- "C:/WADE/"
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Define UI ---------------------------------------------------------------------------------------
ui <- fluidPage(
  tags$style('.container-fluid {background-color: #FFFFFF;}'),
  tags$head(tags$style(
    HTML("#label > input[type='radio'] {
          #  opacity: 0;
          #  position: absolute;
          #}
          label > input[type='radio'] + *::before {
            content: '';
            margin: 4px 0 0;
            width: 13px;
            height: 13px;
            position: absolute;
            margin-left: -20px;
            border-radius: 50%;
            border-style: solid;
            border-width: 0.1rem;
          }
          label > input[type='radio']:checked + *::before {
                        background: radial-gradient(white 0%, white 30%, #8F00FF 30%, #8F00FF);
                                border-color: #8F00FF;
                    }",
    "input:checked + span {
            color: #8F00FF;
         }",
        '#sidebar {
            background-color: #FFFFFF;
            border-color: #696969;
        }
        body, label, input, button, select { 
          font-family: "Arial";
        }')
  )),
  img(src = "PHAC.png"),

  titlePanel(strong(h1("Strep/STI WGS Analysis and Detection of Molecular Markers (WADE)", 
                       align = "center"))),
  
  sidebarLayout(position = "left", 
    sidebarPanel(id = "sidebar",

      selectInput("Org",
                  label = h4("Choose an Organism"),
                  choices = list("GAS",
                                 "GBS",
                                 "PNEUMO",
                                 "GONO",
                                 "Curator"),
                  selected = "GAS"),

      conditionalPanel(condition = "input.Org == 'GONO'",
                       radioButtons("gonotest", h4("Choose an analysis:"),
                       choices = list("AMR Profile" = "AMR_ALL",
                                      "AMR Alleles" = "AMR",
                                      "23S rRNA Alleles" = "rRNA23S",
                                      "MLST" = "MLST",
                                      "NG-STAR" = "NGSTAR",
                                      "NG-MAST" = "NGMAST",
                                      "All Routine Analyses" = "ALL",
                                      #"MasterBlastR" = "MASTERBLASTR",
                                      "WGS Metrics" = "LW_METRICS"),
                       selected = "AMR_ALL")),

      conditionalPanel(condition = "input.Org == 'GAS'",
                       radioButtons("gastest", h4("Choose an analysis:"),
                       choices = list("AMR Profile" = "AMR",
                                      "emm Typing" = "EMM",
                                      "M1UK Typing" = "M1UK",
                                      "MLST" = "MLST",
                                      "16S rRNA" = "rRNA16S",
                                      "Toxin Profile" = "TOXINS",
                                      "All Routine Analyses" = "ALL",
                                      #"MasterBlastR" = "MASTERBLASTR",
                                      "WGS Metrics" = "LW_METRICS"),
                       selected = "AMR")),

      conditionalPanel(condition = "input.Org == 'GBS'",
                       radioButtons("gbstest", h4("Choose an analysis:"),
                       choices = list("AMR profile" = "AMR",
                                      "MLST" = "MLST",
                                      "16S rRNA" = "rRNA16S",
                                      "SeroTypR" = "SERO",
                                      "All Routine Analyses" = "ALL",
                                      #"MasterBlastR" = "MASTERBLASTR",
                                      "WGS Metrics" = "LW_METRICS"),
                       selected = "AMR")),

      conditionalPanel(condition = "input.Org == 'PNEUMO'",
                       radioButtons("spntest", h4("Choose an analysis:"),
                       choices = list("AMR Profile" = "AMR_ALL",
                                      "AMR Alleles" = "AMR",
                                      "23S rRNA Alleles" = "rRNA23S",
                                      "MLST" = "MLST",
                                      "16S rRNA" = "rRNA16S",
                                      "SeroTypR" = "SERO",
                                      "Virulence Factors" = "VIRULENCE",
                                      "All Routine Analyses" = "ALL",
                                      #"MasterBlastR" = "MASTERBLASTR",
                                      "WGS Metrics" = "LW_METRICS"),
                       selected = "AMR_ALL")),
      
      conditionalPanel(condition = "input.gastest == 'M1UK'",
                       radioButtons("variant", h4("Choose a variant:"),
                                    choices = list("M1UK",
                                                   "M1DK"),
                                    selected = "M1UK")),

      conditionalPanel(condition = "input.Org != 'Curator'",
                       textInput("locus", h5("Enter a locus to query or \"list\" for default loci list"),
                       value = "list")),

      conditionalPanel(condition = "input.Org != 'Curator'",
                       textInput("sample", h5("Enter sample number or \"list\" for multiple samples"),
                       value = "list")),

      conditionalPanel(condition = "input.Org == 'Curator'",
                       radioButtons("curatortest", h4("Choose an analysis:"),
                       choices = list("Update Lookup tables" = "UPDATE_LOOKUPS",
                                      "Contamination Check" = "CONTAMINATION_CHECK",
                                      "GC MIC Comparison" = "MIC_CHECK",
                                      "Remove Duplicate Sequences" = "REMOVE_DUPLICATES"),
                       selected = "CONTAMINATION_CHECK")),
      
      conditionalPanel(condition = "input.curatortest == 'UPDATE_LOOKUPS' && input.Org == 'Curator'",
                       radioButtons("lookuptables", h5("Choose an analysis:"),
                                    choices = list("GONO MLST" = "GONO_MLST",
                                                   "GONO NGMAST" = "GONO_NGMAST",
                                                   "GAS MLST" = "GAS_MLST",
                                                   "GAS EMM" = "GAS_EMM",
                                                   "GBS MLST" = "GBS_MLST",
                                                   "PNEUMO MLST" = "PNEUMO_MLST"),
                                    selected = "GONO_MLST")),
      
      conditionalPanel(condition = "(input.curatortest == 'FIND_DUPLICATES' | input.curatortest == 'IRIDA_METADATA') && input.Org == 'Curator'",
                       radioButtons("duplicates_orgs", h5("Choose an organism:"),
                                    choices = list("GONO" = "GONO",
                                                   "STREP" = "STREP"),
                                    selected = "GONO")),
      
      conditionalPanel(condition = "(input.curatortest == 'CONTAMINATION_CHECK'|input.curatortest == 'FILE_RENAME') && input.Org == 'Curator'",
                       radioButtons("curation_orgs", h5("Choose an organism:"),
                                    choices = list("GONO" = "GONO",
                                                   "GAS" = "GAS",
                                                   "GBS" = "GBS",
                                                   "PNEUMO" = "PNEUMO"),
                                    selected = "GONO")),
      
      conditionalPanel(condition = "input.curatortest == 'IRIDA_METADATA' && input.Org == 'Curator'",
                       textInput("batch", h5("Enter Batch Number"),
                                 value = "379")),
      
      actionBttn("action", 
                 label = "Go",
                 style = "gradient",
                 color = "primary"),
      actionBttn("openxls", 
                 label = "Output",
                 style = "gradient",
                 color = "warning"),
      actionBttn("index", 
                 label = "MakeBlastdb",
                 style = "gradient",
                 color = "royal"),
      
      conditionalPanel(condition = "input.curatortest == 'READS_RENAME' && input.Org == 'Curator'",
                       br(),
                       actionBttn("rename", 
                                  label = "Rename Reads",
                                  style = "gradient",
                                  color = "danger")),
    ),

    mainPanel(
      br(),
      textOutput("selected_Org"),
      textOutput("selected_test"),
      textOutput("selected_test2"),
      textOutput('selected_test3'),
      textOutput("selected_test4"),
      textOutput("entered_locus"),
      textOutput("entered_sample"),
      textOutput("button_value"),
      br(),
      textOutput("error_text"),
      textOutput("Status"),
      
      DT::dataTableOutput("profile_table")
    )
  )
)

# Define server logic -------------------------------------------------------------------------------
server <- function(input, output) {

  observeEvent(input$action,
  {
    #remove previous WADE output fasta files
    removefiles(curr_work_dir)
    if(input$Org == "GONO")
    {
      switch(input$gonotest,
             MLST={output.df <- MLST_pipeline(input$Org, "MLST", input$sample, input$locus, curr_work_dir)},
             NGSTAR={output.df <- MLST_pipeline(input$Org, "NGSTAR", input$sample, input$locus, curr_work_dir)},
             NGMAST={output.df <- MLST_pipeline(input$Org, "NGMAST", input$sample, input$locus, curr_work_dir)},
             rRNA23S={output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)},
             #AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             #VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             #AMR = {output.df <- MASTER_pipeline(input$Org, "AMR_ALL", input$sample, "list", curr_work_dir)
             #      output.df <- MLST_pipeline(input$Org, "NGSTAR_AMR", input$sample, input$locus, curr_work_dir)},
             AMR={output.df <- MASTER_pipeline(input$Org, input$gonotest, input$sample, input$locus, curr_work_dir)},
             AMR_ALL={output.df <- MASTER_pipeline(input$Org, "AMR_ALL", input$sample, "list", curr_work_dir)
                      output.df <- MLST_pipeline(input$Org, "NGSTAR_AMR", input$sample, input$locus, curr_work_dir)
                      output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)
                      output.df <- labware_gono_amr(input$Org, curr_work_dir)},
             ALL={output.df <- MASTER_pipeline(input$Org, "AMR_ALL", input$sample, "list", curr_work_dir)
                  output.df <- MLST_pipeline(input$Org, "NGSTAR_AMR", input$sample, input$locus, curr_work_dir)
                  output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)
                  output.df <- labware_gono_amr(input$Org, curr_work_dir)
                  output.df <- MLST_pipeline(input$Org, "MLST", input$sample, input$locus, curr_work_dir)
                  output.df <- MLST_pipeline(input$Org, "NGSTAR", input$sample, input$locus, curr_work_dir)
                  output.df <- MLST_pipeline(input$Org, "NGMAST", input$sample, input$locus, curr_work_dir)
                  beep(5)},
             MASTERBLASTR={output.df <- MASTER_pipeline(input$Org, "MASTERBLASTR", input$sample, input$locus, curr_work_dir)})
    }

    if(input$Org == "GAS")
    {
      switch(input$gastest,
             MLST={output.df <- MLST_pipeline(input$Org, "MLST", input$sample, input$locus, curr_work_dir)},
             EMM={output.df <- EMM_pipeline(input$Org, input$sample, curr_work_dir)},
             M1UK={output.df <- EMM_V_pipeline(input$Org, input$variant, curr_work_dir)},
             #AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             #VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             AMR={output.df <- MASTER_pipeline(input$Org, input$gastest, input$sample, input$locus, curr_work_dir)
                  output.df <- labware_gas_amr(input$Org, curr_work_dir)},
             rRNA16S={output.df <- MASTER_pipeline(input$Org, input$gastest, input$sample, input$locus, curr_work_dir)
                       output.df <- rRNA16S_pipeline(input$Org, curr_work_dir)},
             TOXINS={output.df <- MASTER_pipeline(input$Org, input$gastest, input$sample, input$locus, curr_work_dir)
                     output.df <- labware_gas_toxins(input$Org, curr_work_dir)},
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             ALL={output.df <- MASTER_pipeline(input$Org, "AMR", input$sample, input$locus, curr_work_dir)
                  output.df <- labware_gas_amr(input$Org, curr_work_dir)
                  output.df <- MLST_pipeline(input$Org, "MLST", input$sample, input$locus, curr_work_dir)
                  output.df <- EMM_pipeline(input$Org, input$sample, curr_work_dir)
                  output.df <- MASTER_pipeline(input$Org, "rRNA16S", input$sample, input$locus, curr_work_dir)
                  output.df <- rRNA16S_pipeline(input$Org, curr_work_dir)
                  output.df <- MASTER_pipeline(input$Org, "TOXINS", input$sample, input$locus, curr_work_dir)
                  output.df <- labware_gas_toxins(input$Org, curr_work_dir)
                  beep(5)},
             {output.df <- MASTER_pipeline(input$Org, input$gastest, input$sample, input$locus, curr_work_dir)})
    }

    if (input$Org == "GBS")
    {
      switch(input$gbstest,
             MLST={output.df <- MLST_pipeline(input$Org, "MLST", input$sample, input$locus, curr_work_dir)},
             AMR={output.df <- MASTER_pipeline(input$Org, input$gbstest, input$sample, input$locus, curr_work_dir)
                  output.df <- labware_gbs_amr(input$Org, curr_work_dir)},
             rRNA16S={output.df <- MASTER_pipeline(input$Org, input$gbstest, input$sample, input$locus, curr_work_dir)
                      output.df <- rRNA16S_pipeline(input$Org, curr_work_dir)},
             #AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             #VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             SERO={output.df <- SEROTYPE_pipeline(input$Org, input$sample, input$locus, input$gbstest, curr_work_dir)},
             ALL={output.df <- MASTER_pipeline(input$Org, "AMR", input$sample, input$locus, curr_work_dir)
                  output.df <- labware_gbs_amr(input$Org, curr_work_dir)
                  output.df <- MLST_pipeline(input$Org, "MLST", input$sample, input$locus, curr_work_dir)
                  output.df <- MASTER_pipeline(input$Org, "rRNA16S", input$sample, input$locus, curr_work_dir)
                  output.df <- rRNA16S_pipeline(input$Org, curr_work_dir)
                  output.df <- SEROTYPE_pipeline(input$Org, input$sample, input$locus, "SERO", curr_work_dir)
                  beep(5)},
             {output.df <- MASTER_pipeline(input$Org, input$gbstest, input$sample, input$locus, curr_work_dir)})
    }

    if (input$Org == "PNEUMO")
    {
      switch(input$spntest,
             AMR_ALL={output.df <- MASTER_pipeline(input$Org, input$spntest, input$sample, input$locus, curr_work_dir)
                      output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)
                      output.df <- labware_pneumo_amr(input$Org, curr_work_dir)},
             AMR={output.df <- MASTER_pipeline(input$Org, input$spntest, input$sample, input$locus, curr_work_dir)},
             rRNA16S={output.df <- MASTER_pipeline(input$Org, input$spntest, input$sample, input$locus, curr_work_dir)
                      output.df <- rRNA16S_pipeline(input$Org, curr_work_dir)},
             MLST={output.df <- MLST_pipeline(input$Org, "MLST", input$sample, input$locus, curr_work_dir)},
             rRNA23S={output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)},
             #AMR_DB=(output.df <- AMR_DATABASE_pipeline(input$Org, input$sample, curr_work_dir)),
             #VFDB=(output.df <- VFDB_pipeline(input$Org, input$sample, curr_work_dir)),
             VIRULENCE={output.df <- MASTER_pipeline(input$Org, input$spntest, input$sample, input$locus, curr_work_dir)
                        output.df <- labware_pneumo_virulence(input$Org, curr_work_dir)},
             SERO={output.df <- SEROTYPE_pipeline(input$Org, input$sample, input$locus, input$spntest, curr_work_dir)},
             LW_METRICS={output.df <- metrics(input$Org, curr_work_dir)},
             ALL={output.df <- MASTER_pipeline(input$Org, "AMR_ALL", input$sample, "list", curr_work_dir)
                  output.df <- rRNA23S_pipeline(input$Org, input$sample, curr_work_dir)
                  output.df <- labware_pneumo_amr(input$Org, curr_work_dir)
                  output.df <- MLST_pipeline(input$Org, "MLST", input$sample, input$locus, curr_work_dir)
                  output.df <- MASTER_pipeline(input$Org, "VIRULENCE", input$sample, input$locus, curr_work_dir)
                  output.df <- labware_pneumo_virulence(input$Org, curr_work_dir)
                  output.df <- SEROTYPE_pipeline(input$Org, input$sample, input$locus, "SERO", curr_work_dir)
                  beep(5)},
             {output.df <- MASTER_pipeline(input$Org, input$spntest, input$sample, input$locus, curr_work_dir)})
    }

    if (input$Org == "Curator")
    {
      switch(input$curatortest,
             TREE_ORDER={output.df <- nwk_sort_order(input$Org, curr_work_dir)},
             FILE_RENAME={output.df <- rename_submitted(input$curation_orgs, curr_work_dir)},
             REMOVE_DUPLICATES={output.df <- remove_duplicate_fasta(input$Org, input$allele, curr_work_dir)},
             FIND_DUPLICATES={output.df <- find_sample_duplicates(input$duplicates_orgs, curr_work_dir)},
             CONTAMINATION_CHECK={output.df <- proximity_test(input$curation_orgs, curr_work_dir)},
             MIC_CHECK={output.df <- MIC_compare(curr_work_dir)}
      )
    }

    unlink(paste0(curr_work_dir, "Output/output_profile.csv"))
    write.csv(output.df, paste0(curr_work_dir, "Output/output_profile.csv"), row.names = F)

    output$profile_table <- renderDataTable({output.df})
  })

  observeEvent(input$openxls,
  {
    if(input$locus == "list")
    {
      shell.exec(paste0(curr_work_dir, "Output/output_profile.csv"))
    }
  })

  observeEvent(input$index,
  {
    switch(input$Org,
           GONO = {runblast <- MakeblastDB_pipeline(input$Org, input$gonotest, input$locus, curr_work_dir)},
           GAS = {runblast <- MakeblastDB_pipeline(input$Org, input$gastest, input$locus, curr_work_dir)},
           GBS = {runblast <- MakeblastDB_pipeline(input$Org, input$gbstest, input$locus, curr_work_dir)},
           PNEUMO = {runblast <- MakeblastDB_pipeline(input$Org, input$spntest, input$locus, curr_work_dir)},
           MGEN = {runblast <- MakeblastDB_pipeline(input$Org, input$mgentest, input$locus, curr_work_dir)})
  })
  
  observeEvent(input$rename,
  {
    rename_fastqs()
  })
}

# Run the app -------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)