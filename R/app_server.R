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
                      fold_change = input$logFC,
                      top_tag_plot = 10)
    
    return(p)
    
  })
  
  output$plot_vulcano <- renderPlot(volcano_plot()[[1]])
  
  plot_heatmap_contrast <- reactive({
    
    data <- rna_data()
    
    # wynik <- make_comparison(fit = data$fit, 
    #                          contrast_comp = input$contrast, 
    #                          design = data$design,
    #                          genes = data$genes, 
    #                          complete_list = FALSE, 
    #                          toptags_print = 10, 
    #                          plot_volcano = FALSE, 
    #                          plot_md = FALSE,
    #                          fold_change = 1.5)
    
    wynik <- volcano_plot()
    
    group <- data$strain
    
    p <- draw_heatmap(data_edger = data$data_edger,
                      comparison = wynik[[2]],
                      groups = group,
                      n = input$heatmap_contrast_n)
    
    return(p)
    
  })
  
  output$heatmap_contrast <- renderPlot(plot_heatmap_contrast())
  
  output$choose_contrast <- renderUI({
    if (is.null(input$dane_rna))
      return(NULL)
    
    dane <- rna_data()
    
    grupy <- dane$contrasts
    
    selectInput("contrast", "Choose contrast for the plot",
                choices = grupy, selected = grupy[1])
    
  })
  
  
  output$choose_strains <- renderUI({
    if (is.null(input$dane_rna))
      return(NULL)
    
    dane <- rna_data()
    
    grupy <- dane$strain
    
    # selectInput("strains", 
    #             "Choose strains for the plot",
    #             choices = grupy, 
    #             selected = grupy[1:4], 
    #             multiple = TRUE)
    
    checkboxGroupInput("strains", 
                       "Choose strains for the plot",
                       choices = grupy, 
                       selected = grupy[1:4])
    
  })
  
  
  plot_rpkm_cluster <- reactive({
    
    data <- rna_data()
    
    p <- plot_rpkm_genes(fitted = data$data_rpkm, 
                    genes_positions = data$genes_pos, 
                    strains = input$strains,
                    gene_start = input$gene_start, 
                    gene_end = input$gene_end,
                    plot_type = 'gene_position',
                    flank = 0,
                    log = FALSE)
    
    return(p)
  })
  
  output$plot_rpkm_cluster <- renderPlot(plot_rpkm_cluster())
  
  plot_heatmap_cluster <- reactive({
    
    data <- rna_data()
    
    p <- draw_heatmap_cluster(fit = data$fit,
                              gene_start = input$gene_start,
                              gene_end = input$gene_end,
                              strains = input$strains,
                              genes_positions = data$genes_pos)
    
    return(p)
    
  })
  
  output$heatmap_cluster <- renderPlot(plot_heatmap_cluster())
  
}
