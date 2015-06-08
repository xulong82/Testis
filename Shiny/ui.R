library(shiny)
library(ggvis)
library(Gviz)

load("./gvizList.rdt")
load("./geList.rdt")

sample <- names(gvizList$sample)
transcript <- seqnames(gvizList$seqs)

shinyUI(
  navbarPage("Eif4g3", 
             tabPanel("Expression", 
                      fluidPage(
                        p(strong("Expression"), "profiles of 15937 genes in the 18 samples."),
                        hr(),
                        fluidRow(
                          column(3,
                                 radioButtons("norm", label = "Norm. type", choices = names(geList), selected = "aligned"),
                                 textInput("gene", label = h3("Input gene:"), value = "Actb") 
                          ),
                          column(9,
                                 plotOutput("boxplot"),
                                 tableOutput("summary")
                          ))
                      )),
             tabPanel("Alignment", 
                      fluidPage(
                        p(strong("RNA-seq"), "reads that were aligned with Eif4g3."),
                        hr(),
                        fluidRow(
                          column(3,
                                 selectInput("sample", "Choose a sample", choices = sample, selected = "M1IN"),
                                 selectInput("transcript", "Choose a transcript", choices = transcript, selected = transcript[2]),
                                 helpText("Loading takes 5 sec to 1 mins."),
                                 hr(),
                                 h4("Genetic position in bp"),
                                 numericInput("from", label = "from", value = -1),
                                 numericInput("to", label = "to", value = -1)
                          ),
                          column(9,
                                 plotOutput("gviz")
                          ))
                      )),
             tabPanel("Documents",
                      fluidPage(
                        h4("Reports in html."),
                        p(strong(a(href="report.html", "Uploaded: April 15, 2015."))),
                        helpText("Press <f> to view in fullscreen; press <Win, +/-> to zoom in and out."),
                        hr(),
	h4("Joint lab meeting on May 8, 2015"),
	downloadButton('downloadData1', 'Download'),
	h4("Joint lab meeting on June 4, 2015"),
	downloadButton('downloadData2', 'Download')
                      ))
            
))
