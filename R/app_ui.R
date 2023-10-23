#' The application User-Interface
#' 
#' @param request Internal parameter for `{shiny}`. 
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    # Leave this function for adding external resources
    golem_add_external_resources(),
    # List the first level UI elements here 
    navbarPage(
      "RNAseqHelper",
      theme = shinythemes::shinytheme("united"),
      tabPanel('Filter data',
               sidebarLayout(
                 sidebarPanel(
                   fileInput("dane_rna", 'Load RNA-seq data file in .rds format',
                             accept=c('.rds')),
                   'Filter comparisons',
                   radioButtons('plot_type_comp', 'Choose plot type',
                                choices = c('Volcano plot' = 'vulcano',
                                            'Heatmap' = 'heatmap')),
                   #conditionalPanel('input.plot_type_comp == "vulcano"',
                   uiOutput('choose_contrast'),
                   numericInput('logFC', 'Choose logFC for the comparison', 
                                value = 1.5, min = 1, step = 0.5),
                   conditionalPanel('input.plot_type_comp == "heatmap"',
                                    numericInput('heatmap_contrast_n', 'Choose number of genes for the heatmap',
                                                 value = 100, min = 1, step = 10)
                   ),
                   downloadButton('download_1', 'Download png plot'),
                   numericInput('width_1', 'Plot width [cm]', 25, min = 5, max = 75),
                   numericInput('height_1', 'Plot height [cm]', 20, min = 5, max = 75),
                   numericInput('res_1', 'Resolution', 200, min = 100, max = 500)
                   #)
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel('Table',
                              dataTableOutput('data_all_table')
                     ),
                     tabPanel('Plots',
                              conditionalPanel('input.plot_type_comp == "vulcano"',
                                               plotOutput('plot_vulcano', height = 700)
                              ),
                              conditionalPanel('input.plot_type_comp == "heatmap"',
                                               plotOutput('heatmap_contrast', height = 700)
                              )
                     )
                   )
                 )
               )),
      tabPanel('Compare genes expression',
               sidebarLayout(
                 sidebarPanel(width = 2,
                              radioButtons('plot_type', 'Choose plot type', 
                                           choices = c(
                                             #'Compare genes RPKM value' = 'only_genes',
                                             'RPKM values of a gene' = 'genes_cluster',
                                             'LogFC values' = 'logfc',
                                             'Heatmap' = 'heatmap')),
                              textInput('gene_start', 'Type in first gene name', value = 'SCO0596'),
                              textInput('gene_end', 'Type in last gene name', value = 'SCO0602'),
                              conditionalPanel('input.plot_type == "genes_cluster"',
                                               checkboxInput('log_rpkm', 'Use log10 scale for RPKM?', value = FALSE),
                                               radioButtons('TPM', 'Which data format should be used?', 
                                                            choices = c('RPKM', 'TPM'), selected = 'RPKM', inline = TRUE),
                                               uiOutput('choose_strains')),
                              conditionalPanel('input.plot_type == "logfc"',
                                               uiOutput('choose_contrasts_logfc')),
                              downloadButton('download_2', 'Download png plot'),
                              numericInput('width_2', 'Plot width [cm]', 35, min = 5, max = 50),
                              numericInput('height_2', 'Plot height [cm]', 18, min = 5, max = 75),
                              numericInput('res_2', 'Resolution', 200, min = 100, max = 500)
                 ),
                 mainPanel(width = 10,
                           conditionalPanel('input.plot_type == "genes_cluster"',
                                            plotOutput('plot_rpkm_cluster', height = 700)
                           ),
                           conditionalPanel('input.plot_type == "heatmap"',              
                                            plotOutput('heatmap_cluster', height = 700, width = 400)
                           ),
                           conditionalPanel('input.plot_type == "logfc"',              
                                            plotOutput('plot_logFC', height = 700)
                           )
                           
                 )
               )),
      tabPanel('UpSet Plots',
               sidebarLayout(
                 sidebarPanel(
                   uiOutput('choose_strains_upset'),
                   numericInput('logFC_upset', 'Choose logFC threshold value for the upset plot',
                                value = 1.5, min = 0.5, step = 0.1),
                   downloadButton('download_3', 'Download png plot'),
                   numericInput('width_3', 'Plot width [cm]', 30, min = 5, max = 50),
                   numericInput('height_3', 'Plot height [cm]', 18, min = 5, max = 75),
                   numericInput('res_3', 'Resolution', 200, min = 100, max = 500)
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel('Table',
                              dataTableOutput('table_upset')
                     ),
                     tabPanel('Upset Plot',
                              plotOutput('upset_plot', height = 700)
                     )
                   )
                 )
               ))
    )
  )
}

#' Add external Resources to the Application
#' 
#' This function is internally used to add external 
#' resources inside the Shiny application. 
#' 
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function(){
  
  add_resource_path(
    'www', app_sys('app/www')
  )
  
  tags$head(
    favicon(),
    bundle_resources(
      path = app_sys('app/www'),
      app_title = 'RNAseqHelper'
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert() 
  )
}

