
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#
# 
library(shiny)
library(plotly)
library(knitr)
source("helper.R")

shinyUI(navbarPage(inverse=TRUE,
  title="Snowman - Structural variation detection by genome-wide local assembly",
  tabPanel(title="Home", icon = icon("fa fa-home"),
           sidebarLayout(
              sidebarPanel(width=3,

                 p("",a("Snowman", href="https://github.com/broadinstitute/SnowmanSV", target="_blank" ),
                 "was created by:"
               ),
               #div(
                 img(src="https://avatars1.githubusercontent.com/u/6922120?v=3&s=460", height = 100, width=100, style='border-radius:50px', align='center'),
                 br(),
                 a("Jeremiah Wala", href="https://github.com/jwalabroad", target="_blank"),
                 br(),br(),
                # style = 'float:left; text-align:center'
               #),
               #div(
                 img(src="http://www.nygenome.org/wp-content/uploads/2015/10/Imielinski-headshot.jpeg", height = 100, width=100 , style='border-radius:50px'),
                 br(),
                 a("Marcin Imielinski", href="http://www.nygenome.org/lab-groups-overview/imielinski-lab/", target="_blank"),
               br(),br(),
                # style = 'float:left; text-align:center'
               #),
               #div(
                 img(src="http://springerlab.tch.harvard.edu/springer/uploads/Alumni/cheng-zhongzhang.jpg", height = 100, width=100 , style='border-radius:50px'),
                 br(),
                 a("Cheng-Zhong Zhang", href="http://www.ncbi.nlm.nih.gov/pubmed/?term=Zhang%20CZ%5Bauth%5D", target="_blank"),
          br(),br()
              #   style = 'float:left; text-align:center'
              # )#,
               #div(style = 'clear: left;')
             ),
             mainPanel(
               img(src='schematic_snowman.png', align = "right", width="100%")
             )
          )
  ),
  navbarMenu(title="Documentation",
    tabPanel(title="README",
             uiOutput('markdown')
    )
  ),
  # navbarMenu(title="De novo assemblies",
  #    tabPanel(title="HCC1143",
  #             h3("Whole-genome de novo assembly with DISCOVAR of HCC1143 with 250-bp PCR-free Illumina Paired-End Reads"),
  #             div(dataTableOutput("discovar_hcc1143_events"), style="font-size:70%")
  #    )
  # ),
  navbarMenu("NA12878",
  tabPanel(title="Data tables",icon = icon("fa fa-th-list"),
           h3("NA12878 deletions from 1000 Genomes (Mills et al)"),
           div(dataTableOutput("truth_NA12878_table"), style="font-size:70%"),
           h3("NA12878 deletions from 1000 Genomes + Validated LUMPY+DELLY+GASVPro+Pindel (Layer et al)"),
           div(dataTableOutput("truth_NA12878_table2"), style="font-size:70%"),
           h3("NA12878 deletions from Snowman"),
           div(dataTableOutput("snowman_NA12878_table"), style="font-size:70%"),
           h3("NA12878 deletions from Pindel"),
           div(dataTableOutput("pindel_NA12878_table"), style="font-size:70%")
           ),
  tabPanel(title="Comparison figures",
           h3("Events detected"),
           plotlyOutput("na12878_set")
           )
  ),
  navbarMenu("Simulated tumor",
             tabPanel(title="Events table", icon = icon("fa fa-th-list"),
                      h3("Simulated SVs"),
                      div(dataTableOutput("simulated_events_table"), style="font-size:70%")
                      ),
            tabPanel(title="Caller plots",
                     selectInput("sim_select", label="Select data to plot", choices=c("Snowman 10X", "Snowman 5X", "Snowman 2X", "Pindel 10X", "Pindel 5X", "Pindel 2X"), selected="Snowman 10X"),
                     selectInput("tum_support_sim", label="Select minimum T_ALT", choices=c("PASS-ONLY",2:10), selected="PASS-ONLY"),
                     selectInput("sim_plot_evdnc", label="Select EVDNC to plot", choices=c("ALL", "INDEL", "ALL-SV", "ASSMB","ASDIS","DSCRD","COMPL"), selected="ALL"),
                     
                     plotOutput("sim_plot"),
                     plotlyOutput("sim_type_plot", width="80%")
                     )
  ),
  navbarMenu("HCC1143",
             tabPanel(title="Snowman", icon = icon("fa fa-th-list"),
                      h3("Snowman somatic events (101-bp Illumina)"),
                      #div(dataTableOutput("pindel_hcc1143_101"), style="font-size:70%"),
                      h3("Snowman somatic events (250-bp PCR-free Illumina)")
                      #div(dataTableOutput("pindel_hcc1143_250"), style="font-size:70%")
             ),
             tabPanel(title="LUMPY", icon = icon("fa fa-th-list"),
                      h3("LUMPY somatic events (101-bp Illumina)"),
                      #div(dataTableOutput("pindel_hcc1143_101"), style="font-size:70%"),
                      h3("LUMPY somatic events (250-bp PCR-free Illumina)")
                      #div(dataTableOutput("pindel_hcc1143_250"), style="font-size:70%")
             ),
             tabPanel(title="Pindel", icon = icon("fa fa-th-list"),
                      h3("Pindel somatic events (101-bp Illumina)"),
                      #div(dataTableOutput("pindel_hcc1143_101"), style="font-size:70%"),
                      h3("Pindel somatic events (250-bp PCR-free Illumina)")
                      #div(dataTableOutput("pindel_hcc1143_250"), style="font-size:70%")
             ),
             tabPanel(title="DELLY", icon = icon("fa fa-th-list"),
                      h3("Delly somatic events (101-bp Illumina)"),
                      #div(dataTableOutput("delly_hcc1143_101"), style="font-size:70%"),
                      h3("Delly somatic events (250-bp PCR-free Illumina)")
                      #div(dataTableOutput("delly_hcc1143_250"), style="font-size:70%")
             ),
             tabPanel(title="Strelka", icon = icon("fa fa-th-list"),
                      h3("Strelka somatic events (101-bp Illumina)"),
                      #div(dataTableOutput("strelka_hcc1143_101"), style="font-size:70%"),
                      h3("Strelka somatic events (250-bp PCR-free Illumina)")
                      #div(dataTableOutput("strelka_hcc1143_250"), style="font-size:70%")
             ),
             tabPanel(title="Platypus", icon = icon("fa fa-th-list"),
                      h3("Platypus somatic events (101-bp Illumina)"),
                      #div(dataTableOutput("platypus_hcc1143_101"), style="font-size:70%"),
                      h3("Platypus somatic events (250-bp PCR-free Illumina)")
                      #div(dataTableOutput("platypus_hcc1143_250"), style="font-size:70%")
             )
  )
  ))


