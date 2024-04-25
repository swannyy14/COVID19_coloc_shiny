library(plotly)
library(DT)
# ui for top control panel
top_control_ui <- tagList(
  selectInput(
    "query_select",
    "Summary Level:",
    choices = c(
      "Gene" = "Gene",
      "Clump" = "Clump"
    )
  ),
  radioButtons(
    "qtl_type_select",
    "QTL Type:",
    choices = c(
      "All", "eQTL", "mQTL", "pQTL", "sQTL"
    ),
    selected = c(
      "eQTL"
    )
  ),
  sliderInput(
    "pp4_min",
    "PP4 Minimum Threshold:",
    min = 0, max = 1, value = 0.75
  ),
  sliderInput(
    "height_scale",
    "Heatmap Height scale:",
    min = 0.9, max = 2, value = 1
  ),
  checkboxInput(
    "max_pp4",
    "Heatmap - Display only Max PP4 Per Gene-QTL Map",
    value = FALSE
  )
)

# UI for gene level summary
gene_summary_ui <- tagList(
  radioButtons(
    "gene_summary_level",
    "Search by Gene or Molecular Phenotype:",
    choices = c(
      "Gene" = "Gene",
      "Molecular Phenotype" = "MP"
    ), inline = FALSE
  ),
  selectInput(
    "gene_select",
    NULL,
    choices = NULL
  )
)

gene_summary_main <- tagList(
  tabsetPanel(
    tabPanel(
      "Gene-level",
      plotlyOutput("plot_signif_map", height = 600)
    ),
    tabPanel(
      "Clump-level",
      plotlyOutput("plot_gene_eqtl_heatmap")
    ),
    tabPanel(
      "Locus Zoom",
      br(),
      dataTableOutput("table_gene_eqtl_heatmap"),
      hr(),
      column(6,
        imageOutput("LocusZoomPlotGene")
      )
    )
  )
)

# UI for clump level summary
clump_summary_ui <- tagList(
  selectInput(
    "clump_select",
    "Select Clump:",
    choices = NULL
  )
)

clump_summary_main <- tagList(
  tabsetPanel(
    tabPanel(
      "Clump-level",
      plotlyOutput("plot_clump_eqtl_heatmap")
    ),
    tabPanel(
      "Locus Zoom",
      br(),
      dataTableOutput("table_clump_eqtl_heatmap"),
      hr(),
      column(6,
        imageOutput("LocusZoomPlotClump") 
      )
    )
  )
)
