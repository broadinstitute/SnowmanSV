
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

source("helper.R")

shinyServer(function(input, output) {

  output$snowman_NA12878_table <- renderDataTable(snow.dels.NA12878, options=list(pageLength = 10))
  output$truth_NA12878_table <- renderDataTable(truth.NA12878$dt, options=list(pageLength = 10))
  output$truth_NA12878_table2 <- renderDataTable(truth2.NA12878$dt, options=list(pageLength = 10))
  
  output$simulated_events_table <- renderDataTable(events.d1, options=list(pageLength=10, autoWidth=TRUE))
  
  output$distPlot <- renderPlot({

    # generate bins based on input$bins from ui.R
    x    <- faithful[, 2]
    bins <- seq(min(x), max(x), length.out = input$bins + 1)

    # draw the histogram with the specified number of bins
    hist(x, breaks = bins, col = 'darkgray', border = 'white')

  })

})
