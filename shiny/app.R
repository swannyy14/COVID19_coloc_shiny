library(shiny)
library(tidyverse)
library(heatmaply)
library(plotly)
library(DT)

# additional UI script saved in R/ui_elements.R
ui <- fluidPage(
  titlePanel(
    title = "COVID-19 Coloc Results Explorer"
  ),
  top_control_ui,
  conditionalPanel(
    "input.query_select == 'Gene'",
    gene_summary_ui
  ),
  conditionalPanel(
    "input.query_select == 'Clump'",
    clump_summary_ui
  )
)

server <- function(input, output, session) {
  # subset coloc summary data frames by QTL type and PP4
  coloc_summary_server <- reactive({
    if (input$qtl_type_select != "All") {
      coloc_summary_tmp <- coloc_summary %>% filter(QTL_type == input$qtl_type_select)
    } else {
      coloc_summary_tmp <- coloc_summary
    }
    coloc_summary_tmp %>%
      filter(PP4 >= input$pp4_min)
  })
  clumps_df_server <- reactive({
    clumps_df %>%
      filter(Clump_name %in% coloc_summary_server()$GWAS_clump)
  })
  gene_summary_level <- reactive({
    if (input$gene_summary_level == "Gene") {
      "Gene"
    } else {
      "MP"
    }
  })
  genes_server <- reactive({
    sort(coloc_summary_server()[[gene_summary_level()]])
  })
  
  # update gene and clump selection list
  observe({
    updateSelectInput(
      session = session,
      inputId = "gene_select",
      choices = genes_server()
    )
  })
  observe({
    updateSelectInput(
      session = session,
      inputId = "clump_select",
      choices = clumps_df_server()$Clump_name
    )
  })
  
  # Update gene selection
  mygene <- reactive({
    if (input$query_select == "Gene") {
      input$gene_select
    }
  })
  
  # Update Clump selection (If summary level is gene, then clump including the gene is automatically selected)
  myclumps <- reactive({
    if (input$query_select == "Gene") {
      unique(coloc_summary_server()$GWAS_clump[coloc_summary_server()[[gene_summary_level()]] == mygene()])
    } else if (input$query_select == "Clump") {
      input$clump_select
    }
  })
  
  # plot significant eqtl map when gene is selected
  output$plot_signif_map <- renderPlotly({
    validate(need(mygene(), label = "Gene"))
    plotly_signif_qtl_map_for_gene(coloc_summary_server(), mygene(), gene_summary_level(), covid_phenotypes, input$qtl_type_select)
  })
  
  # plot heatmap of posterior probabilities for Gene vs QTL Map for specifed clump
  output$plot_gene_eqtl_heatmap <- renderPlotly({
    validate(need(myclumps(), label = "Clump"))
    plot_heatmaply_for_clump(coloc_summary_server(), myclumps(), clumps_df, gene_pos_df, input$height_scale_gene, choose_max_pp4_per_gene = input$max_pp4_gene)
  })
  
  # Print QTL summary stats table
  output$table_gene_eqtl_heatmap <- renderDataTable({
    validate(need(mygene(), label = "Gene"))
    if(gene_summary_level() == 'Gene') printtable <- unique(coloc_summary[coloc_summary$Gene == mygene(),])
    if(gene_summary_level() == 'MP') printtable <- unique(coloc_summary[coloc_summary$MP == mygene(),])
    if(input$qtl_type_select != "All") printtable <- printtable[printtable$QTL_type == input$qtl_type_select,]
    printtable <- printtable[printtable$PP4 >= input$pp4_min,]
    printtable
      }, selection='single')
  
  output$table_clump_eqtl_heatmap <- renderDataTable({
    validate(need(myclumps(), label = "Clump"))
    printtable <- unique(coloc_summary[coloc_summary$GWAS_clump == myclumps(),])
    if(input$qtl_type_select != "All") printtable <- printtable[printtable$QTL_type == input$qtl_type_select,]
    printtable <- printtable[printtable$PP4 >= input$pp4_min,]
    printtable
  }, selection='single')
  
  output$LocusZoomPlotGene <- renderImage({
    validate(need(mygene(), label = "Gene"))
    if(gene_summary_level() == 'Gene') printtable <- unique(coloc_summary[coloc_summary$Gene == mygene(),])
    if(gene_summary_level() == 'MP') printtable <- unique(coloc_summary[coloc_summary$MP == mygene(),])
    if(input$qtl_type_select != "All") printtable <- printtable[printtable$QTL_type == input$qtl_type_select,]
    printtable <- printtable[printtable$PP4 >= input$pp4_min,]
    path <- paste0(printtable[input$table_gene_eqtl_heatmap_rows_selected,'LZplot'])
    print(path)
    list(src=path, width = 500)
      }, deleteFile = F)

  output$LocusZoomPlotClump <- renderImage({
    validate(need(myclumps(), label = "Clump"))
    printtable <- unique(coloc_summary[coloc_summary$GWAS_clump == myclumps(),])
    if(input$qtl_type_select != "All") printtable <- printtable[printtable$QTL_type == input$qtl_type_select,]
    printtable <- printtable[printtable$PP4 >= input$pp4_min,]
    path <- paste0(printtable[input$table_clump_eqtl_heatmap_rows_selected,'LZplot'])
    print(path)
    list(src=path, width = 500)
  }, deleteFile = F)  
  
  # plot heatmap of posterior prbability for genes in a clump
  output$plot_clump_eqtl_heatmap <- renderPlotly({
    validate(need(myclumps(), label = "Clump"))
    plot_heatmaply_for_clump(coloc_summary_server(), myclumps(), clumps_df, gene_pos_df, input$height_scale_clump, choose_max_pp4_per_gene = input$max_pp4_clump)
  })
}

shinyApp(ui, server)