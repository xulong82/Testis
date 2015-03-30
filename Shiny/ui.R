library(shiny)
library(ggvis)

shinyUI(
  fluidPage(
    titlePanel("A Fun Project"),
    p(strong("Figure 1:"), "Expression profile of 12333 genes across 18 samples.", strong(a(href="norm.html", "Click here for the protocol"))),
    hr(),
    sidebarLayout(
      sidebarPanel(
	      radioButtons("norm", label = h3("ERCC Norm"), choices = list("raw" = 1, "norm" = 2), selected = 2),
        textInput("gene", label = h3("Input gene:"), value = "Hspa2"),
	      hr(),
	      p("EXCEL file"),
	      downloadButton('downloadData', 'Download')
      ),
      mainPanel(
        plotOutput("boxplot"),
	tableOutput("summary")
      )
    )
))
