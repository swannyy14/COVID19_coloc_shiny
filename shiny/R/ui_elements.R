library(plotly)
library(DT)
library(bslib)
library(bsicons)
library(shinyWidgets)

tab <- function(...) {
  shiny::tabPanel(..., class = "p-2 border border-top-0 rounded-bottom")
}

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

# UI for clump level summary
clump_summary_ui <- tagList(
  selectInput(
    "clump_select",
    "Select Clump:",
    choices = NULL
  )
)

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
  conditionalPanel(
    "input.query_select == 'Gene'",
    gene_summary_ui
  ),
  conditionalPanel(
    "input.query_select == 'Clump'",
    clump_summary_ui
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
  )#,
  # sliderInput(
  #   "height_scale",
  #   "Heatmap Height scale:",
  #   min = 0.9, max = 2, value = 1
  # ),
  # checkboxInput(
  #   "max_pp4",
  #   "Heatmap - Display only Max PP4 Per Gene-QTL Map",
  #   value = FALSE
  # )
)

gene_summary_main <- tagList(
  tabsetPanel(
    tab(
      "Gene-level",
      card(
        plotlyOutput("plot_signif_map", height = 600)
      )
    ),
    tab(
      "Clump-level",
      layout_columns(
        dropdownButton(
          sliderInput(
            "height_scale_gene_summ",
            "Heatmap Height scale:",
            min = 0.9, max = 2, value = 1
          ),
          checkboxInput(
            "max_pp4_gene_summ",
            "Heatmap - Display only Max PP4 Per Gene-QTL Map",
            value = FALSE
          ),
          circle = FALSE,
          status = "info",
          size = "sm",
          icon = bs_icon("gear-fill")
        ),
        card(
          plotlyOutput("plot_gene_eqtl_heatmap"),
          full_screen = TRUE
        ), 
        col_widths = c(1,-11,12)
      )
    ),
    tab(
      "Locus Zoom",
      br(),
      card(
        div(
          dataTableOutput("table_gene_eqtl_heatmap")
        )
      ),
      hr(),
      card(
        card_header("Locus Zoom Plot"),
        imageOutput("LocusZoomPlotGene"),
        height = 1000
      )
    )
  )
)

clump_summary_main <- tagList(
  tabsetPanel(
    tab(
      "Clump-level",
      layout_columns(
        dropdownButton(
          sliderInput(
            "height_scale_clump_summ",
            "Heatmap Height scale:",
            min = 0.9, max = 2, value = 1
          ),
          checkboxInput(
            "max_pp4_clump_summ",
            "Heatmap - Display only Max PP4 Per Gene-QTL Map",
            value = FALSE
          ),
          circle = FALSE,
          status = "info",
          size = "sm",
          icon = bs_icon("gear-fill")
        ),
        card(
          plotlyOutput("plot_clump_eqtl_heatmap"),
          full_screen = TRUE
        ), 
        col_widths = c(1,-11,12)
      )
    ),
    tab(
      "Locus Zoom",
      br(),
      card(
        div(
          dataTableOutput("table_clump_eqtl_heatmap")
        )
      ),
      hr(),
      card(
        card_header("Locus Zoom Plot"),
        imageOutput("LocusZoomPlotClump"),
        height = 1000
      )
    )
  )
)
