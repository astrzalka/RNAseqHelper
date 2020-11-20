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
                   conditionalPanel('input.plot_type_comp == "vulcano"',
                   uiOutput('choose_contrast')
                   )
                 ),
                 mainPanel(
                   tabsetPanel(
                     tabPanel('Table',
                              dataTableOutput('data_all_table')
                     ),
                     tabPanel('Plots',
                              plotOutput('plot_vulcano', height = 800))
                   )
                 )
               )),
      tabPanel('Compare genes expression',
               sidebarLayout(
                 sidebarPanel(
                   radioButtons('plot_type', 'Choose plot type', 
                                choices = c('Compare genes RPKM value' = 'only_genes',
                                            'RPKM values of a gene cluster' = 'genes_cluster',
                                            'Heatmap' = 'heatmap')),
                   textInput('gene_start', 'Type in first gene name'),
                   textInput('gene_end', 'Type in last gene name')
                 ),
                 mainPanel()
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

