MAPD_GUI <- function(pics){
  library(shiny)
  library(plotly)
  
  ui <- fluidPage(
    titlePanel("MADP"),
    sidebarLayout(
      sidebarPanel(
        sliderInput("inds", "index of pic:", min = 1, max = length(pics$pics), value = 1, step=1)
      ),
      mainPanel(
        plotlyOutput("Plot"),
        plotlyOutput("Plot1")
      )
    )
  )
  
  server <- function(input, output) {
    output$Plot <- renderPlotly({
      pic.i <- pics$pics[[input$inds]]
      peak.i <- pics$peaks[[input$inds]]
      plot_ly(x=pics$scantime[pic.i[,1]], y=pic.i[,2], type = 'scatter', mode = 'lines', showlegend= FALSE) %>%
        layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
               yaxis = list(title = 'Intensity')) %>%
        add_markers(x=pics$scantime[pic.i[peak.i$peakIndex,1]], y=pic.i[peak.i$peakIndex,2])
    })
    output$Plot1 <- renderPlotly({
      pic.i <- pics$pics[[input$inds]]
      pic.i <- pic.i[!is.na(pic.i[,3]),]
      plot_ly(x=pics$scantime[pic.i[,1]], y=pic.i[,3], color=pic.i[,2] ,type = 'scatter') %>%
        layout(xaxis = list(tick0 = 0, title = 'Retention Time (s)'),
               yaxis = list(title = 'M/Z'))
    })

  }
  shinyApp(ui = ui, server = server)
}