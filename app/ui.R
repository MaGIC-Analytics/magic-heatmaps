library(shiny)
require(shinyjs)
library(shinythemes)
require(shinycssloaders)
library(shinyWidgets)

library(DT)
library(tidyverse)
library(data.table)
library(colourpicker)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(msigdbr)

tagList(
    tags$head(
        includeHTML(("www/GA.html")),
        tags$style(type = 'text/css','.navbar-brand{display:none;}'),
        tags$style(HTML("
            .control-group-panel {
                border: 1px solid #ddd;
                border-radius: 6px;
                padding: 10px 12px;
                margin-bottom: 10px;
                background-color: #f9f9f9;
            }
            .control-group-title {
                font-weight: bold;
                font-size: 14px;
                color: #0F344C;
                margin-bottom: 8px;
            }
            #show_help_float {
                position: fixed;
                bottom: 28px;
                right: 28px;
                z-index: 9999;
                border-radius: 50%;
                width: 46px;
                height: 46px;
                font-size: 20px;
                padding: 0;
                box-shadow: 0 3px 8px rgba(0,0,0,0.25);
            }
        "))
    ),
    ## Global always-visible help button (fixed bottom-right)
    actionButton("show_help_float", label=NULL,
        icon=icon("circle-question"),
        title="Help & documentation",
        class="btn btn-info"
    ),
    fluidPage(theme = shinytheme('yeti'),
            windowTitle = "MaGIC RNA-seq Heatmap Tool",
            useShinyjs(),
            titlePanel(
                fluidRow(
                column(2, tags$a(href='http://www.bioinformagic.io/', tags$img(height=75, src="MaGIC_Icon_0f344c.svg")), align='center'),
                column(10, fluidRow(
                    column(10, h1(strong('MaGIC RNA-seq Heatmap Tool'), align='center', style="color:#0F344C;"))
                ))
                ),
                windowTitle = "MaGIC RNA-seq Heatmap Tool"),
                tags$style(type='text/css', '.navbar{font-size:20px;}'),
                tags$style(type='text/css', '.nav-tabs{padding-bottom:20px;}'),
                tags$style(type='text/css', '.navbar-default{background-color:#0F344C;}'),
                tags$style(type='text/css', HTML('.navbar { background-color: #0F344C;}
                          .tab-panel{ background-color: #0F344C;}
                          .navbar-default .navbar-nav > .active > a,
                           .navbar-default .navbar-nav > .active > a:focus,
                           .navbar-default .navbar-nav > .active > a:hover {
                                color: white;
                                background-color: #008cba;
                            }')
                          ),
                tags$head(tags$style(".modal-dialog{ width:1300px}")),

        navbarPage(title="", id='NAVTABS',

        ## Intro Page
##########################################################################################################################################################
            tabPanel('Introduction',
                fluidRow(
                    column(2),
                    column(8,
                        column(12, align="center", style="margin-bottom:25px;",
                            h3(markdown("Welcome to the RNA-seq Heatmap Tool by the
                            [Molecular and Genomics Informatics Core (MaGIC)](http://www.bioinformagic.io)."))),
                        hr(),
                        h4("How to Use This Tool", style="color:#0F344C;"),
                        tags$ol(
                            tags$li(strong("Navigate to the Data Input tab."),
                                " Upload your Expression Matrix and Sample Metadata files, or click 'Load Demo Data' to explore with a built-in example."),
                            tags$li(strong("Submit your data."),
                                " Click the Submit button. The Heatmap tab will become visible once data is successfully loaded."),
                            tags$li(strong("Select genes of interest."),
                                " On the Heatmap tab, display all genes, enter a custom gene list, or query public gene set databases (GO, KEGG, MSigDB Hallmark)."),
                            tags$li(strong("Customize your heatmap."),
                                " Use the control panels in the left sidebar to adjust colors, annotations, fonts, clustering, row/column splitting, and legend settings."),
                            tags$li(strong("Resize and download."),
                                " Use the Resize panel to fine-tune plot dimensions, then download in your preferred format (PNG, PDF, SVG, TIFF, etc.).")
                        ),
                        hr(),
                        h4("Required Input Data Formats", style="color:#0F344C;"),
                        fluidRow(
                            column(6,
                                div(class="control-group-panel",
                                    h5(strong("Expression Matrix"), style="color:#0F344C;"),
                                    tags$ul(
                                        tags$li("File format: CSV or TSV"),
                                        tags$li("Rows: Genes (one gene per row)"),
                                        tags$li("Columns: Samples (one sample per column)"),
                                        tags$li("First column: Gene identifiers (symbols or IDs)"),
                                        tags$li("All remaining columns: Numeric expression values (e.g., log2 normalized counts, VST, TPM)")
                                    ),
                                    tags$pre("Gene,   Sample1, Sample2\nBRCA1,  6.5,     7.1\nTP53,   8.2,     7.9")
                                )
                            ),
                            column(6,
                                div(class="control-group-panel",
                                    h5(strong("Sample Metadata"), style="color:#0F344C;"),
                                    tags$ul(
                                        tags$li("File format: CSV or TSV"),
                                        tags$li("Rows: Samples (one sample per row)"),
                                        tags$li("First column: Sample names — must match matrix column names exactly"),
                                        tags$li("Additional columns: Categorical or numeric metadata variables (e.g., Group, Gender, Tissue)")
                                    ),
                                    tags$pre("Sample,  Group,   Gender\nSample1, Control, Male\nSample2, Treated, Female")
                                )
                            )
                        ),
                        hr()
                    ),
                    column(2)
                )
            ),


        ## Data Input Page
##########################################################################################################################################################
            tabPanel('Data Input',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('Input Data', align='center'),
                            hr(),
                            materialSwitch("DemoData", label="Upload custom data", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.DemoData",
                                h4("Expression Matrix", style="color:#0F344C;"),
                                fileInput('matrix_file', 'Upload Matrix File (CSV/TSV)',
                                    accept=c('text/csv', 'text/comma-separated-values, text/plain', '.csv',
                                             'text/tsv', 'text/tab-separated-values, text/plain', '.tsv'),
                                    multiple=FALSE
                                ),
                                h4("Sample Metadata", style="color:#0F344C;"),
                                fileInput('metadata_file', 'Upload Metadata File (CSV/TSV)',
                                    accept=c('text/csv', 'text/comma-separated-values, text/plain', '.csv',
                                             'text/tsv', 'text/tab-separated-values, text/plain', '.tsv'),
                                    multiple=FALSE
                                ),
                                hr(),
                                uiOutput('gene_col_selector'),
                                actionButton('submit', "Submit Data", class='btn btn-info btn-block')
                            ),
                            conditionalPanel("input.DemoData==false",
                                p("Use the pre-loaded RNA-seq demo data to explore the tool's features."),
                                p(em("Demo dataset: 30 cancer-related genes across 9 samples in 3 treatment groups.")),
                                hr(),
                                actionButton('demo_submit', "Load Demo Data", class='btn btn-success btn-block')
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='InputTables',
                            tabPanel(title='Expression Matrix', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('matrix_table')
                                )
                            ),
                            tabPanel(title='Sample Metadata', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('metadata_table')
                                )
                            )
                        )
                    )
                )
            ),


        ## Heatmap Page (hidden until data submitted)
##########################################################################################################################################################
            tabPanel('Heatmap',
                fluidRow(
                    column(3,
                        wellPanel(

                            ## Gene Selection
                            h5(strong("Gene Selection"), style="color:#0F344C; margin-top:4px;"),
                            hr(),
                            selectInput("gene_source", label=NULL,
                                choices=c("All Genes"="all", "Custom Gene List"="custom", "Database Gene Set"="database"),
                                selected="all", width="100%"
                            ),
                            conditionalPanel("input.gene_source == 'custom'",
                                textAreaInput("custom_genes",
                                    "Enter gene names (one per line or comma-separated):",
                                    rows=5, placeholder="TP53\nBRCA1\nMYC\n...")
                            ),
                            conditionalPanel("input.gene_source == 'database'",
                                selectInput("db_source", "Database:",
                                    choices=c(
                                        "MSigDB Hallmark (H)"="H",
                                        "MSigDB C2 KEGG"="C2_KEGG",
                                        "MSigDB C5 GO Biological Process"="C5_BP",
                                        "MSigDB C5 GO Molecular Function"="C5_MF"
                                    )
                                ),
                                selectizeInput("db_geneset", "Gene Set:", choices=NULL,
                                    options=list(placeholder='Type to search gene sets...', maxOptions=5000)),
                                selectInput("db_species", "Species:",
                                    choices=c("Homo sapiens"="Homo sapiens", "Mus musculus"="Mus musculus"),
                                    selected="Homo sapiens")
                            ),
                            p(em("All other plot options update automatically. Click below to apply a new gene selection."), style="font-size:11px; color:#888; margin-top:6px;"),
                            actionButton("apply_genes", "Apply Gene Selection",
                                class="btn btn-info btn-block", style="margin-top:4px;"),
                            hr(),

                            ## Color
                            materialSwitch("show_color", label="Color Options", value=TRUE, right=TRUE, status='info'),
                            conditionalPanel("input.show_color",
                                hr(),
                                radioButtons("color_mode", label=NULL, inline=TRUE,
                                    choices=c("Two-color balance"="two", "Color palette"="palette"),
                                    selected="two"
                                ),
                                conditionalPanel("input.color_mode == 'two'",
                                    column(6, colourInput("color_high", "High color", "#d73027")),
                                    column(6, colourInput("color_low",  "Low color",  "#4575b4")),
                                    div(style="clear:both;")
                                ),
                                conditionalPanel("input.color_mode == 'palette'",
                                    selectInput("color_palette", "Color Palette:",
                                        choices=c(
                                            "RdBu (diverging)"="RdBu",
                                            "RdYlBu"="RdYlBu",
                                            "Spectral"="Spectral",
                                            "Blues"="Blues",
                                            "Reds"="Reds",
                                            "Greens"="Greens",
                                            "YlOrRd"="YlOrRd",
                                            "PuOr"="PuOr",
                                            "PRGn"="PRGn"
                                        ),
                                        selected="RdBu"
                                    ),
                                    materialSwitch("color_reverse", label="Reverse palette",
                                        value=FALSE, right=TRUE, status='warning')
                                ),
                                materialSwitch("scale_rows", label="Scale rows (z-score)",
                                    value=TRUE, right=TRUE, status='info')
                            ),

                            ## Annotation
                            materialSwitch("show_anno", label="Annotation", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_anno",
                                hr(),
                                uiOutput('annotation_col_selector'),
                                selectInput("anno_position", "Annotation Position:",
                                    choices=c("Top"="top", "Bottom"="bottom"),
                                    selected="top"),
                                sliderInput("anno_bar_size", "Bar Size (mm):", min=1, max=30, step=1, value=5),
                                sliderInput("anno_font_size", "Label Font Size (pt):", min=4, max=20, step=1, value=10),
                                uiOutput('anno_color_ui')
                            ),

                            ## Fonts & Labels
                            materialSwitch("show_fonts", label="Fonts & Labels", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_fonts",
                                hr(),
                                materialSwitch("show_row_names", label="Show row names (genes)",
                                    value=TRUE, right=TRUE, status='info'),
                                conditionalPanel("input.show_row_names",
                                    sliderInput("row_font_size", "Row label size:", min=4, max=20, step=1, value=8),
                                    sliderInput("row_font_angle", "Row label angle:", min=0, max=360, step=5, value=0)
                                ),
                                materialSwitch("show_col_names", label="Show column names (samples)",
                                    value=TRUE, right=TRUE, status='info'),
                                conditionalPanel("input.show_col_names",
                                    sliderInput("col_font_size", "Column label size:", min=4, max=20, step=1, value=8),
                                    sliderInput("col_font_angle", "Column label angle:", min=0, max=360, step=5, value=45)
                                )
                            ),

                            ## Clustering
                            materialSwitch("show_clustering", label="Clustering", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_clustering",
                                hr(),
                                materialSwitch("cluster_rows", label="Cluster rows",
                                    value=TRUE, right=TRUE, status='info'),
                                conditionalPanel("input.cluster_rows",
                                    sliderInput("row_dend_size", "Row dendrogram width (mm):", min=1, max=30, step=1, value=10),
                                    materialSwitch("reorder_rows", label="Reorder rows by dendrogram",
                                        value=TRUE, right=TRUE, status='warning')
                                ),
                                materialSwitch("cluster_cols", label="Cluster columns",
                                    value=TRUE, right=TRUE, status='info'),
                                conditionalPanel("input.cluster_cols",
                                    sliderInput("col_dend_size", "Column dendrogram height (mm):", min=1, max=30, step=1, value=10),
                                    materialSwitch("reorder_cols", label="Reorder columns by dendrogram",
                                        value=TRUE, right=TRUE, status='warning')
                                )
                            ),

                            ## Splitting
                            materialSwitch("show_splitting", label="Splitting", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_splitting",
                                hr(),
                                radioButtons("split_mode", label="Split rows/columns by:", inline=TRUE,
                                    choices=c("None"="none", "K-means"="kmeans", "Metadata Column"="metadata"),
                                    selected="none"
                                ),
                                conditionalPanel("input.split_mode == 'kmeans'",
                                    sliderInput("row_km", "Row k (clusters):", min=1, max=10, step=1, value=2),
                                    sliderInput("col_km", "Column k (clusters):", min=1, max=10, step=1, value=1)
                                ),
                                conditionalPanel("input.split_mode == 'metadata'",
                                    uiOutput('split_col_selector')
                                )
                            ),

                            ## Legend
                            materialSwitch("show_legend", label="Legend", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_legend",
                                hr(),
                                sliderInput("legend_title_size", "Title font size:", min=6, max=24, step=1, value=12),
                                sliderInput("legend_font_size", "Label font size:", min=6, max=24, step=1, value=10),
                                sliderInput("legend_bar_height", "Bar height (mm):", min=5, max=80, step=5, value=40),
                                sliderInput("legend_grid_height", "Grid cell height (mm):", min=1, max=20, step=1, value=4),
                                sliderInput("legend_grid_width", "Grid cell width (mm):", min=1, max=20, step=1, value=4)
                            ),

                            ## Resize
                            materialSwitch("show_resize", label="Resize Plot", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.show_resize",
                                hr(),
                                sliderInput("hm_height", "Plot height (px):", min=200, max=2000, step=50, value=700),
                                sliderInput("hm_width",  "Plot width (px):",  min=200, max=2000, step=50, value=900)
                            )

                        )# end wellPanel sidebar
                    ),
                    column(9,
                        tabsetPanel(id='HeatmapTabs',
                            tabPanel(title='Heatmap', hr(),
                                fluidRow(style="margin: 0 8px 4px 0;",
                                    column(12, align="right",
                                        actionButton("show_code_modal", label=NULL,
                                            icon=icon("file-code"),
                                            title="View R code to reproduce this plot",
                                            class="btn btn-default btn-sm",
                                            style="border-radius:6px; font-size:16px; padding:4px 8px;"
                                        )
                                    )
                                ),
                                hr(),
                                div(style="overflow-x:auto; width:100%;",
                                    withSpinner(type=6, color='#5bc0de',
                                        plotOutput("heatmap_out", height='100%')
                                    )
                                ),
                                div(style="margin-top:30px; text-align:center; padding-bottom:50px;",
                                    div(style="display:inline-block; width:250px; margin-bottom:10px;",
                                        selectInput("hm_download_format", "Download format:",
                                            choices=c('png','pdf','svg','tiff','jpeg','eps'))
                                    ),
                                    br(),
                                    downloadButton('download_heatmap', 'Download Heatmap')
                                )
                            ),
                            tabPanel(title='Active Gene List', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('gene_list_table')
                                )
                            )
                        )
                    )
                )
            ),


        ## Footer
##########################################################################################################################################################
            tags$footer(
                wellPanel(
                    fluidRow(
                        column(4, align='center',
                        tags$a(href="https://github.com/MaGIC-Analytics/magic-heatmap", icon("github", "fa-3x")),
                        tags$h4('GitHub to submit issues/requests')
                        ),
                        column(4, align='center',
                        tags$a(href="http://www.bioinformagic.io/", icon("magic", "fa-3x")),
                        tags$h4('MaGIC Home Page')
                        ),
                        column(4, align='center',
                        tags$a(href="https://github.com/MaGIC-Analytics", icon("address-card", "fa-3x")),
                        tags$h4("Developer's Page")
                        )
                    ),
                    fluidRow(
                        column(12, align='center',
                            HTML('<a href="https://www.youtube.com/watch?v=dQw4w9WgXcQ">
                            <p>&copy;
                                <script language="javascript" type="text/javascript">
                                var today = new Date()
                                var year = today.getFullYear()
                                document.write(year)
                                </script>
                            </p>
                            </a>
                            ')
                        )
                    )
                )
            )
        )# Ends navbarPage
    )# Ends fluidPage
)# Ends tagList
