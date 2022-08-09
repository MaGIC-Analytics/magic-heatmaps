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

library(pheatmap)

tagList(
    tags$head(
        includeHTML(("www/GA.html")),
        tags$style(type = 'text/css','.navbar-brand{display:none;}')
    ),
    fluidPage(theme = shinytheme('yeti'),
            windowTitle = "MaGIC Heatmap Tool",
            titlePanel(
                fluidRow(
                column(2, tags$a(href='http://www.bioinformagic.io/', tags$img(height =75 , src = "MaGIC_Icon_0f344c.svg")), align = 'center'), 
                column(10, fluidRow(
                  column(10, h1(strong('MaGIC Heatmap Tool'), align = 'center',style="color:#0F344C;")),
                  ))
                ),
                windowTitle = "MaGIC Heatmap Tool" ),
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
        navbarPage(title ="", id='NAVTABS',

        ## Intro Page
##########################################################################################################################################################
            tabPanel('Introduction',
                fluidRow(
                    column(2,
                    ),
                    column(8,
                        column(12, align = "center", 
                            style="margin-bottom:25px;",
                            h3(markdown("Welcome to the Heatmap Tool by the
                            [Molecular and Genomics Informatics Core (MaGIC)](http://www.bioinformagic.io)."))),
                        hr(),
                    ),
                    column(2,
                    )
                ),
                fluidRow(
                    column(2,
                    ),
                    column(8, align='center', style="margin-bottom:30px;",
                        img(src="Fire-animated.gif")
                    ),
                    column(2,
                    )
                ),
                fluidRow(
                    column(3,
                    ),
                    column(6,
                        column(12, align = "left", 
                            style="margin-bottom:25px;",
                            h4("About Heatmaps"),
                            markdown("Heatmaps are a great way of visualizing your count data across multiple conditions. This can be represented as counts per sample, average counts per group, or log2FC values vs control per group. This can be a great way to demonstrate both inter- and intra-group differences in your data. 
                            "),
                            hr(),
                            h4("Data Sources"),
                            markdown("The data for heatmaps can come from any type of count data. A few examples are RNA-seq normalized hit counts, Luminex assay counts, plant petal size/color, average expression values per sample, log2FC values for multiple conditions, etc.
                            Generally we recommend scaling your data by row (each gene/value), though if you are visualizing something like log2FC values we would recommend no scaling. 
                            "),
                            hr(),
                            h4("Minimum requirements"),
                            markdown("To use this tool, at minimum you must have a tsv/csv table containing with the first row containing your row identifiers (for example Gene IDs), followed by a column with count data per sample (for example VST normalized hit counts).
                            Additionally you need a metadata table. The first column should contain the sample names matching the columns in the count data. Each subsequent column can include any non-measured data, such as grouping variables.
                            
                            ")),
                        hr(),
                    ),
                    column(3,
                    )
                ),

            ),


        ## Input Page
##########################################################################################################################################################
            tabPanel('Input Data',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('Input Data', align='center'),
                            hr(),
                            materialSwitch("DemoData", label="Upload your own data", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.DemoData",
                                fileInput('counts_file','Select your counts file',
                                accept=c(
                                    'text/csv',
                                    'text/comma-separated-values, text/plain',
                                    '.csv',
                                    'text/tsv',
                                    'text/tab-separated-values, text/plain',
                                    '.tsv'), 
                                multiple=FALSE
                                ),
                                fileInput('metadata_file','Select your metadata file',
                                accept=c(
                                    'text/csv',
                                    'text/comma-separated-values, text/plain',
                                    '.csv',
                                    'text/tsv',
                                    'text/tab-separated-values, text/plain',
                                    '.tsv'), 
                                multiple=FALSE
                                ),
                                uiOutput('data_selector')
                            ),
                            conditionalPanel("input.DemoData==false",
                                uiOutput('data_demo')
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='InputTables',
                            tabPanel(title='Data table', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('input_table')
                                )
                            ),
                            tabPanel(title='Metadata table', hr(),
                                withSpinner(type=6, color='#5bc0de',
                                    dataTableOutput('meta_table')
                                )
                            )
                        )
                    )
                )
            ),

        ## Heatmaps Page
##########################################################################################################################################################
            tabPanel('Heatmaps',
                fluidRow(
                    column(3,
                        wellPanel(
                            h2('Plot Options', align='center'),
                            materialSwitch("ClustOpts", label='Show clustering options', value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.ClustOpts",
                                materialSwitch("CClust", label='Cluster the columns', value=TRUE, right=TRUE, status='info'),
                                materialSwitch("RClust", label='Cluster the rows', value=TRUE, right=TRUE, status='info'),
                                radioButtons('Scale',label='Scale by: ', inline=TRUE, choices=c('Row'='row','Column'='column','None'='none'), selected='row')
                            ),
                            materialSwitch("LabOpts", label='Label options', value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.LabOpts",
                                materialSwitch("CName", label='Show column names', value=FALSE, right=TRUE, status='info'),
                                materialSwitch("RName", label='Show row names', value=FALSE, right=TRUE, status='info'),
                                sliderInput("HXsize", 'Label Size', min=1, max=30, step=1, value=10),
                                radioButtons("Hang", label='X-Axis Angle', inline=TRUE,
                                    choices=c('0'=0,'45'=45,'90'=90,'270'=270,'315'=315), selected=45)
                            ),
                            materialSwitch("BarOpts", label="Show bar options", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.BarOpts", 
                                selectInput("HMGroup", label="Group by: ", choices=NULL)
                            ),
                            materialSwitch("ColorOpts", label="Color options", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.ColorOpts",
                                selectInput("CPal", label="Color palette: ", choices=NULL)
                            ),
                            materialSwitch("Resize", label="Resize Image", value=FALSE, right=TRUE, status='info'),
                            conditionalPanel("input.Resize",
                                sliderInput('HHeight', label='Plot Heights: ', min=50, max=2000, step=10, value=800),
                                sliderInput('HWidth', label='Plot Widths: ',  min=50, max=2000, step=10, value=800)
                            )
                        )
                    ),
                    column(9,
                        tabsetPanel(id='Plot',
                            tabPanel(title='Heatmap', hr(),
                                withSpinner(type=6, color='#5bc0de',  
                                    plotOutput("heatmap_out", height='100%')
                                ),
                                fluidRow(align='center',style="margin-top:25px;",
                                    column(12, selectInput("DownHFormat", label='Choose download format', choices=c('jpeg','png','tiff'))),
                                    column(12, downloadButton('DownloadH', 'Download the Heatmap'),style="margin-bottom:50px;")
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
                        tags$a(href="https://github.com/MaGIC-Analytics/{$TOOL}", icon("github", "fa-3x")),
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
        )#Ends navbarPage,
    )#Ends fluidpage
)#Ends tagList