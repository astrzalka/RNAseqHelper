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
  library(magrittr)
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
  
  vulcano_plot <- reactive({
    
    data <- rna_data()
    
    p <- plot_volcano(fit = data$fit,
                      contrast_comp =  input$contrast,
                      design = data$design,
                      genes = data$genes,
                      fold_change = input$logFC,
                      top_tag_plot = 10)
    
    return(p)
    
  })
  
  output$plot_vulcano <- renderPlot(vulcano_plot()[[1]])
  
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
    
    wynik <- vulcano_plot()
    
    group <- data$strain
    
    p <- draw_heatmap(data_edger = data$data_edger,
                      comparison = wynik[[2]],
                      groups = group,
                      n = input$heatmap_contrast_n)
    
    return(p)
    
  })
  
  output$heatmap_contrast <- renderPlot(plot_heatmap_contrast())
  
  
  filter_table_upset <- reactive({
    
    data <- rna_data()
    
    table <- data$data_all
    
    strains <- input$strains_upset
    
    stays  <- which(colnames(table) %in% c("genes", "name", "product", strains))
    
    table_zostaje <- table[,stays]
    
    table_signif <- table_zostaje %>% 
      tidyr::gather('comparison', 'logFC', -genes, -name, -product) %>%
      dplyr::group_by(genes) %>% 
      dplyr::mutate(gene_NA = anyNA(logFC)) %>% 
      dplyr::filter(!gene_NA) %>% 
      tidyr::spread(key = 'comparison', value = 'logFC') %>%
      dplyr::select(-gene_NA)
    
    table_upset <- table_zostaje %>% 
      tidyr::gather('comparison', 'logFC', -genes, -name, -product) %>%
      dplyr::mutate(logFC = ifelse(abs(logFC) < input$logFC_upset, 0, 1),
                    logFC = ifelse(is.na(logFC), 0, logFC)) %>%
      tidyr::spread(key = 'comparison', value = 'logFC') %>%
      dplyr::select(-product, -name)
    
    
    
    return(list(table_signif, table_upset))
    
  })
  
  output$table_upset <- renderDataTable(filter_table_upset()[[1]])
  
  plot_upset <- reactive({
    
    table <- filter_table_upset()[[2]]
    
    
    p <- UpSetR::upset(data = table,
                       order.by = 'freq',
                       nsets = 10,
                       text.scale = 2)
    
    return(p)
    
  })
  
  output$upset_plot <- renderPlot(plot_upset())
  
  
  ############################ create UI parts from the data ############################
  
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
  
  output$choose_contrasts_logfc <- renderUI({
    if (is.null(input$dane_rna))
      return(NULL)
    
    dane <- rna_data()
    
    grupy <- dane$contrasts
    
    checkboxGroupInput("contrasts", 
                       "Choose contrasts for the plot",
                       choices = grupy, selected = grupy[1:2])
    
  })
  
  output$choose_strains_upset <- renderUI({
    if (is.null(input$dane_rna))
      return(NULL)
    
    dane <- rna_data()
    
    grupy <- colnames(dane$data_all)
    
    grupy <- grupy[!(grupy %in% c("genes", "name", "product"))]
    
    # selectInput("strains", 
    #             "Choose strains for the plot",
    #             choices = grupy, 
    #             selected = grupy[1:4], 
    #             multiple = TRUE)
    
    checkboxGroupInput("strains_upset", 
                       "Choose columns for the plot",
                       choices = grupy, 
                       selected = grupy[1:2])
    
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
                         log = input$log_rpkm,
                         format = input$TPM,
                         fit = data$fit)
    
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
  
  plot_logFC <- reactive({
    
    data <- rna_data()
    
    p <- plot_logfc_genes(logFC = data$data_all_notsig, 
                          genes_positions = data$genes_pos, 
                          contrasts = input$contrasts,
                          gene_start = input$gene_start, 
                          gene_end = input$gene_end,
                          flank = 0)
    
    return(p)
  })
  
  output$plot_logFC <- renderPlot(plot_logFC())
  
  ########################### Download plots #######################################
  
  prepare_plot_1 <- reactive({
    
    if(input$plot_type_comp == 'vulcano'){
      return(vulcano_plot()[[1]])
    }
    
    if(input$plot_type_comp == 'heatmap'){
      return(plot_heatmap_contrast())
    }
    
    
  }) 
  
  output$download_1 <- downloadHandler(
    filename = function() { paste(input$dataset, '.png', sep='') },
    content = function(file) {
      png(file, res = input$res_1, width = input$width_1, input$height_1, unit = 'cm')
      print(prepare_plot_1())
      dev.off()
    })
  
  
  prepare_plot_2 <- reactive({
    
    if(input$plot_type == 'genes_cluster'){
      return(plot_rpkm_cluster())
    }
    
    if(input$plot_type == 'heatmap'){
      return(plot_heatmap_cluster())
    }
    
    
  }) 
  
  output$download_2 <- downloadHandler(
    filename = function() { paste(input$dataset, '.png', sep='') },
    content = function(file) {
      png(file, res = input$res_2, width = input$width_2, input$height_2, unit = 'cm')
      print(prepare_plot_2())
      dev.off()
    })
  
  prepare_plot_3 <- reactive({
    
    return(plot_upset())
    
  }) 
  
  output$download_3 <- downloadHandler(
    filename = function() { paste(input$dataset, '.png', sep='') },
    content = function(file) {
      png(file, res = input$res_3, width = input$width_3, input$height_3, unit = 'cm')
      print(prepare_plot_3())
      dev.off()
    })
  
}
