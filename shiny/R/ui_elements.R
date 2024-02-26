library(plotly)
library(DT)
# ui for top control panel
top_control_ui <- tagList(
  wellPanel(
    fluidRow(
      # column for choosing summary level
      column(
        4, 
        selectInput(
          "query_select",
          "Summary Level:",
          choices = c(
            "Gene" = "Gene",
            "Clump" = "Clump"
          )
        )
      ),
      # column for choosing qtl type
      column(
        3,
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
        offset = 1
      ),
      # column for choosing pp4
      column(
        4,
        sliderInput(
          "pp4_min",
          "PP4 Minimum Threshold:",
          min = 0, max = 1, value = 0.75
        )
      )
    )
  )
)

# UI for gene level summary
gene_summary_ui <- tagList(
  wellPanel(
    fluidRow(
      # column for choosing gene
      column(
        6,
        column(
          11,
          fluidRow(
            radioButtons(
              "gene_summary_level",
              "Search by Gene or Molecular Phenotype:",
              choices = c(
                "Gene" = "Gene",
                "Molecular Phenotype" = "MP"
              ), inline = FALSE
            )
          ),
          fluidRow(
            selectInput(
              "gene_select",
              NULL,
              choices = NULL
            )
          ), offset = 1
        )
      ),
      # column for height slider and heatmap option
      column(
        6,
        column(
          10,
          fluidRow(
            sliderInput(
              "height_scale_gene",
              "Height scale:",
              min = 0.9, max = 2, value = 1
            )
          ), 
          fluidRow(
            checkboxInput(
              "max_pp4_gene",
              "Heatmap - Max PP4 Per Gene-Tissue Association",
              value = FALSE
            )
          ),
          offset = 1
        )
      )
    )
  ),
  # plot output
  fluidRow(
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
        "Select Locus Zoom",
        dataTableOutput("table_gene_eqtl_heatmap")
      ),
      tabPanel(
        "Locus Zoom",
          imageOutput("LocusZoomPlotGene")
      )
    )
  )
)

# UI for clump level summary
clump_summary_ui <- tagList(
  wellPanel(
    fluidRow(
      # column for selecting summary level
      column(4, selectInput(
        "clump_select",
        "Select Clump:",
        choices = NULL
      )),
      # column for height slider
      column(4, sliderInput(
        "height_scale_clump",
        "Height scale:",
        min = 0.9, max = 2, value = 1
      )),
      # slider for heatmap option
      column(4, checkboxInput(
        "max_pp4_clump",
        "Heatmap - Max PP4 Per Gene-QTL Map",
        value = FALSE
      ))
    )
  ),
  # plot output
  fluidRow(
    tabsetPanel(
      tabPanel(
        "Clump-level",
        plotlyOutput("plot_clump_eqtl_heatmap")
      ),
      tabPanel(
        "Select Locus Zoom",
        dataTableOutput("table_clump_eqtl_heatmap")
      ),
      tabPanel(
        "Locus Zoom",
        imageOutput("LocusZoomPlotClump")
      )
    )
  )
)
