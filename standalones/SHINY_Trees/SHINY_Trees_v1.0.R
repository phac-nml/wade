##################################################################
####        Shiny Trees Phylogenetic Tree Annotator           ####
####       Authors: Walter Demczuk & Shelley Peterson         ####
####                    Date: 2025-05-13                      ####
##################################################################

# Inputs: Newick file (NOT NEXUS) and metadata.csv
# Metadata file must have:
# must have a "Date.Collected" column for the temporal plots to work
#                            "Year" column to plot discrete year values.
#                            ** DO NOT put Tree Sort column in metadata.csv.
# Custom colours and apply shapes only available using "custom" plots
# Must have ColourList.csv and ShapeList.csv files for custom colours and shapes.
# Custom colours are applied alphabetically in selected column as per order in local custom colour list.

library(tidyverse)
library(scales)
library(shiny)
library(shinyWidgets)
library(ggtree)
library(phytools)
library(gridExtra)

#library(plyr)
#library(tidyselect)
#library(stringr)
#library(Biostrings)
#library(DT)
#library(data.table)
#library(RColorBrewer)
#library(lubridate)

############################ USER INPUT ########################################
curr_work_dir <- "C:/SHINY_Trees/"
local_output_dir <- paste(curr_work_dir, "Output/", sep="")
################################################################################

# Define UI --------------------------------------------------------------------
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

  titlePanel(strong(h1("Strep/STI Phylogenetic Tree Annotation", 
                       align = "center"))),

  sidebarLayout(position = "left",
    sidebarPanel(id = "sidebar",
      selectInput("plottomake",
                  label = h5("Plot to make"),
                  choices = list("Tree" = "rectangular",
                                 "Radial (unrooted)" = "equal_angle",
                                 "Temporal",
                                 "Tree + Temporal"
                                 #"Temporal facets"
                  ),
                  selected = "Tree"),

      radioButtons("plotcolours", label = h5("Plot Colours"), inline = TRUE,
                   choices = list("Default" = "default",
                                  "Custom" = "custom"),
                   selected = "default"  ),

      fileInput("SNVFile", h5("Newick Tree file")),
      fileInput("DataFile", h5("Metadata file")),

      uiOutput("tipcolours"),
      uiOutput("highlightcolour"),
      uiOutput("colourother"),
      uiOutput("tipshapes"),
      uiOutput("highlightshape"),
      uiOutput("nodesize"),
      uiOutput("makeplot"),
      uiOutput("viewplot")
    ),

    mainPanel(
      textOutput("edgetable"),
      plotOutput("Trees")
    )
  )
)

# Define server logic ----------------------------------------------------------
server <- function(input, output) {

  # populates pick lists based on file selected
  filedata <- reactive({
    infile_data <-input$DataFile
    if(is.null(infile_data)){return(NULL)}
    temp <- read.csv(infile_data$datapath)
    temp[order(temp[,1]),]
  })

  ##### Colours #####
  output$tipcolours <- renderUI({
    df <- filedata()
    if(is.null(df)) {return(NULL)} else
    {df <- tibble(df, None = "")}
    items <- names(df)

    selectInput("colourtipsby",
                label = h5("Colour Tips By:"),
                choices = items,
                selected = "Province")
  })

  output$highlightcolour <- renderUI({
    df <- filedata()
    if(is.null(df)) {return(NULL)} else
    {df <- tibble(df, None = "")}
    items2 <- c(lapply(unique(df[input$colourtipsby]), sort, decreasing = FALSE), "None")

    selectInput("coloursingle",
                label = h5("Highlight in Red:"),
                choices = items2,
                selected = "None")
  })

  output$colourother <- renderUI({
    df <- filedata()
    if(is.null(df)) {return(NULL)} else
    {df <- tibble(df, None = "")}
  
    radioButtons("colourother", h5("Colour Other Tips:"),
                 choices = list("Colour" = "ColourOther",
                                "No Colour" = "NoColourOther"),
                 selected = "ColourOther")
  })

  ##### Shapes #####
  output$tipshapes <- renderUI({
    df <- filedata()
    if(is.null(df)) {return(NULL)} else
    {df <- tibble(df, None = "")}
    items3 <- names(df)
    
    selectInput("shapetipsby",
                label = h5("Apply Shape to:"),
                choices = items3,
                selected = "None")
  })

  output$highlightshape <- renderUI({
    df <- filedata()
    if(is.null(df)) {return(NULL)} else
    {df <- tibble(df, None = "")}

    items4 <- c(lapply(unique(df[input$shapetipsby]), sort, decreasing = FALSE), "None")

    selectInput("shapesingle",
                label = h5("Apply Triangles to:"),
                choices = items4,
                selected = "None")
  })
  
  ##### Size #####
  output$nodesize <- renderUI({
    df <- filedata()
    if (is.null(df)) {return(NULL)} else
      sliderInput("node_size", "Node Size:",
                  min = 1, max = 5,
                  value = 3.0, step = 0.5)
  })

  ##### Buttons to generate plot #####
  output$makeplot <- renderUI({
    df <- filedata()
    if(is.null(df)) {return(NULL)} else
      actionBttn("update_plot", 
                 label = "Draw Plot",
                 style = "gradient",
                 color = "warning")
    })

  output$viewplot <- renderUI({
    df <- filedata()
    if(is.null(df)) {return(NULL)} else
      actionBttn("openxls", 
                 label = "View Plot",
                 style = "gradient",
                 color = "royal")
  })

  #make plot when button clicked
  observeEvent(input$update_plot,
  {
    output$Trees <- renderPlot({

      #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Set up plot Info
      # Format Metadata
      df <- filedata()
      if(is.null(df)) {return(NULL)} else
      {
        df <- tibble(df, None = "")
        df$Date.Collected <- as.Date(df$Date.Collected)
        df$Year <- as.character(df$Year)
        df <- df %>% rename_at(vars(any_of(c("NML.Number", "NML Number", "NML_Number"))), ~"Strain_ID")
      }

      # Get Newick Tree
      # Remember to export from figtree as newick format
      infile_tree <- input$SNVFile
      if (is.null(infile_tree)) return(NULL) else
      {
        # Get plot type
        plot_type <- input$plottomake     #tree, radial, temporal, etc
        plot_colours <- input$plotcolours #default or custom colours/shapes
        tipscategory <- input$colourtipsby
        nodesize <- input$node_size
        
        #Get tree
        phylo <- read.newick(infile_tree$datapath)

        # parameters for different plot types
        if(plot_type == "rectangular"| plot_type == "Tree + Temporal")
        {
          treestyle <- "rectangular"
          plotx <- NULL
          ploty <- -9
          fontsize <- 3
          linesize <- 0.5
          offset <- NULL
          pdfwidth <- 8.5
          pdfheight <- 10.5
        } else if(plot_type == "equal_angle")
        {
          treestyle <- "equal_angle"
          plotx <- 0
          ploty <- -0.1
          fontsize <- 3
          linesize <- 0.05
          offset <- -0.00025
          pdfwidth <- 5.0
          pdfheight <- 5.0
        } 
        
        # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> CUSTOM SETTINGS
        my_categories <- unique(df[[tipscategory]])
        
        ##### Colours #####
        # Get custom colour list or default colours
        ColourFile <- paste0(curr_work_dir, "ColourList.csV")
        if(file.exists(ColourFile) & plot_colours == "custom")
        {
          my_colours_table <- as_tibble(read.csv(ColourFile, stringsAsFactors = FALSE))
          my_colours <- my_colours_table$Colours[1:length(my_categories)]
        } else
        {
          #cat("\n\n Custom colour list not found. (local temp\\ColourList.csV)\n")
          my_colours <- hue_pal()(length(my_categories))
        }
        
        # Make colour list based on selected field
        names(my_colours) <- str_sort(my_categories)

        # Highlight Single Value
        if(input$coloursingle != "None")
        {
          my_colours[input$coloursingle] <- "red"
        }
        
        # Grey Out Other Tips
        if(input$colourother == "NoColourOther")
        {
          df[[tipscategory]][df[[tipscategory]] != input$coloursingle] <- "Other"
          my_colours["Other"] <- "lightgrey"
        }
        
        ##### Shapes #####
        if(input$shapetipsby != "None")
        {
          shapecategory <- input$shapetipsby
          my_shape_categories <- unique(df[[shapecategory]])
          # Get custom shapes list or default shapes
          ShapeFile <- paste0(curr_work_dir, "ShapeList.csV")
          if(file.exists(ShapeFile))
          {
            my_shapes_table <- as_tibble(read.csv(ShapeFile, stringsAsFactors = FALSE))
            my_shapes <- my_shapes_table$Shapes[1:length(my_shape_categories)]
          } else
          {
            cat("\n\n Custom shape list not found. (ShapeList.csV)\n")
            my_shapes[1:length(my_shape_categories)] <- 21
          }
          # Make shape list based on selected field
          names(my_shapes) <- str_sort(my_shape_categories)
        
          # Highlight Single Value
          if(input$shapesingle != "None")
          {
            my_shapes[input$shapesingle] <- 24
            df[[shapecategory]][df[[shapecategory]] != input$shapesingle] <- "Other"
            my_shapes["Other"] <- 21
          }  
        } 
        
        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Sort Order
        treedata <- fortify(phylo)
        treetips <- subset(treedata, isTip)
        tree_order <- tibble(treetips) %>% arrange(dplyr::desc(y)) %>% select(label)
        tree_order$TreeSort <- 1:dim(tree_order)[1]
        ordereddf <- merge(df, tree_order, by.x = "Strain_ID", by.y = "label") %>%
          arrange(TreeSort)
        write.csv(tree_order, paste0(local_output_dir, "TreeOrder_R.csv"), na = "", row.names = F)

        #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Plots
        #======================================================================= Tree
        if(plot_type == "rectangular"|plot_type == "equal_angle"|plot_type == "Tree + Temporal")
        {
          if(input$shapetipsby != "None"| input$shapesingle != "None") #-------- With shapes
          {
            tree <- ggtree(phylo, layout = treestyle) %<+% df +
              geom_tippoint(aes(fill = .data[[tipscategory]], shape = .data[[shapecategory]]), 
                            size = nodesize, color = "black", stroke = 0.25, na.rm = FALSE)+
              geom_treescale(x = plotx, y = -0.1, fontsize = fontsize, 
                             linesize = linesize, offset = -0.00025)+
              scale_fill_manual(values = my_colours)+
              scale_shape_manual(values = my_shapes)+
              theme(legend.position = "right")
          } else #-------------------------------------------------------------- No Shapes
          {
            tree <- ggtree(phylo, layout = treestyle) %<+% df +
              geom_tippoint(aes(fill = .data[[tipscategory]]), shape = 21, size = nodesize,
                            color = "black", stroke = 0.25, na.rm = FALSE)+
              geom_treescale(x = plotx, y = -0.1, fontsize = fontsize, 
                             linesize = linesize, offset = -0.00025)+
              scale_fill_manual(values = my_colours)+
              theme(legend.position = "right")
          }
          if(plot_type != "Tree + Temporal")
          {
            plot(tree)
        
            pdf(paste(local_output_dir, "tree_plot.pdf", sep = ""), width = 8.5, height = 10.5)
            plot(tree)
            dev.off()
          }
        }
         #====================================================================== Temporal
        if(plot_type == "Temporal"|plot_type == "Tree + Temporal") 
        {
          if(input$shapetipsby != "None"| input$shapesingle != "None") #-------- With shapes
          {
            temp_plot <- ggplot(ordereddf, aes(x = Date.Collected, y = TreeSort))+
              ggtitle("Temporal Distribution of Isolates")+
              geom_point(aes(fill = .data[[tipscategory]], shape = .data[[shapecategory]]), 
                         size=nodesize, colour = "black", stroke = 0.5, alpha =1, na.rm = FALSE)+
              theme_classic()+
              scale_y_reverse()+
              scale_fill_manual(values = my_colours)+
              scale_shape_manual(values = my_shapes)+
              theme(axis.text.x = element_text(vjust = 0.5, angle = 90, hjust = 0),
                    panel.grid.major.x = element_line(colour = 'grey92', linewidth = 0.5),
                    panel.grid.major.y = element_line(colour = 'grey92', linewidth = 0.5),
                    legend.position = "right")+
              scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month", date_labels = "%Y")
          } else #-------------------------------------------------------------- No Shapes
          {
            temp_plot <- ggplot(ordereddf, aes(x = Date.Collected, y = TreeSort))+
              ggtitle("Temporal Distribution of Isolates")+
              geom_point(shape = 21, aes(fill = .data[[tipscategory]]), size=nodesize, 
                         colour = "black", stroke = 0.5, alpha =1, na.rm = FALSE)+
              theme_classic()+
              scale_y_reverse()+
              scale_fill_manual(values = my_colours)+
              theme(axis.text.x = element_text(vjust = 0.5, angle = 90, hjust = 0),
                    panel.grid.major.x = element_line(colour = 'grey92', linewidth = 0.5),
                    panel.grid.major.y = element_line(colour = 'grey92', linewidth = 0.5),
                    legend.position = "right")+
              scale_x_date(date_breaks = "1 year", date_minor_breaks = "1 month", date_labels = "%Y")
          }
          if(plot_type != "Tree + Temporal")
          {
            plot(temp_plot)
          
            pdf(paste(local_output_dir, "tree_plot.pdf", sep = ""), width = 8.5, height = 11.0)
            plot(temp_plot)
            dev.off()
          }
        }
        #======================================================================= Tree + Temporal
        if(plot_type == "Tree + Temporal")
        {
          tree <- tree+
            theme(plot.margin = unit(c(0.5, 0, 1, 0), "cm"))+ #top, right, bottom, left
            theme(legend.position = "none")
          
          temp_plot <- temp_plot+
            ggtitle(NULL)+
            theme(axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.line.y = element_blank(),
                  axis.title.x = element_blank())+
            theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) #top, right, bottom, left
          
          grid.arrange(tree, temp_plot,ncol=2)
          
          pdf(paste(local_output_dir, "tree_plot.pdf", sep = ""), width = 11.0, height = 8.5)
          grid.arrange(tree, temp_plot,ncol=2)
          dev.off()
        }  
      }
    })  #end of renderplot
  }) # end of observe event

  observeEvent(input$openxls,
  {
    shell.exec(paste(local_output_dir, "tree_plot.pdf", sep = ""))
  })
} # end server

# Run the app -------------------------------------------------------------------------------------
shinyApp(ui = ui, server = server)



