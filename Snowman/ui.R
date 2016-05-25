
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

##source("helper.R")

shinyUI(navbarPage(
  title="Snowman - Structural variation detection by genome-wide local assembly",
  tabPanel(title="Description",
           #sidebarLayout(
          # sidebarPanel(width=10,
               
                 p("",a("Snowman", href="https://github.com/broadinstitute/SnowmanSV", target="_blank" ),
                 "was created by:"
               ),
               div(
                 img(src="https://avatars1.githubusercontent.com/u/6922120?v=3&s=460", height = 100, width=100, style='border-radius:50px'), 
                 br(),
                 a("Jeremiah Wala", href="https://github.com/jwalabroad", target="_blank"),
                 style = 'float:left; text-align:center'
               ),
               div(
                 img(src="http://www.nygenome.org/wp-content/uploads/2015/10/Imielinski-headshot.jpeg", height = 100, width=100 , style='border-radius:50px'),
                 br(),
                 a("Marcin Imielinski", href="http://www.nygenome.org/lab-groups-overview/imielinski-lab/", target="_blank"),
                 style = 'float:left; text-align:center'
               ),
               div(
                 img(src="http://springerlab.tch.harvard.edu/springer/uploads/Alumni/cheng-zhongzhang.jpg", height = 100, width=100 , style='border-radius:50px'),
                 br(),
                 a("Cheng-Zhong Zhang", href="http://www.ncbi.nlm.nih.gov/pubmed/?term=Zhang%20CZ%5Bauth%5D", target="_blank"),
                 style = 'float:left; text-align:center'
               )#,
               #div(style = 'clear: left;')
             #),
             #mainPanel(p("sdf"))
          # )
  ),
  tabPanel(title="Analysis of NA12878",
           div(dataTableOutput("snowman_NA12878_table"), style="font-size:70%")
           )
))