#' The application server-side
#' 
#' @param input,output,session Internal parameters for {shiny}. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
#' 

# change maximum data upload size
options(shiny.maxRequestSize=300*1024^2) 

app_server <- function( input, output, session ) {
  # List the first level callModules here
  
  rna_data <- reactive({
    
    inFile <- input$dane_rna
    if (is.null(inFile))
      return(NULL)
    d <- readRDS(file = inFile$datapath)
    return(d)
  })
  
  table_genes <- reactive({
    
    data <- rna_data()
    
    table <- data$data_all
    
    return(table)
    
  })
  
  output$data_all_table <- renderDataTable(table_genes())
  
  volcano_plot <- reactive({
    
    data <- rna_data()
    
    p <- plot_volcano(fit = data$fit,
                      contrast_comp =  input$contrast,
                      design = data$design,
                      genes = data$genes,
                      fold_change = 1.5,
                      top_tag_plot = 10)
    
    return(p)
    
  })
  
  output$plot_vulcano <- renderPlot(volcano_plot())
  
  output$choose_contrast <- renderUI({
    if (is.null(input$dane_rna))
      return(NULL)

    dane <- rna_data()
    
    grupy <- dane$contrasts
    
    selectInput("contrast", "Choose contrast for vulcano plot",
                choices = grupy, selected = grupy[1])
    
  })
  
}
