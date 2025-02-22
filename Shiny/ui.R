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
	downloadButton('downloadData2', 'Download'),
	h4("Joint lab meeting on August 7, 2015"),
	downloadButton('downloadData3', 'Download'),
	h5("TableS1.", strong(a(href="TableS1(nonp_up).xlsx", "Up-regulated genes between non-polysomic samples (Hspa2)."))), 
	h5("TableS2.", strong(a(href="TableS2(nonp_down).xlsx", "Down-regulated genes between non-polysomic samples."))), 
	h5("TableS3.", strong(a(href="TableS3(IN).xlsx", "Differential genes between input samples."))), 
	h5("TableS4.", strong(a(href="TableS4(poly).xlsx", "Differential genes between polysomic samples."))), 
	h4("Joint lab meeting on September 4, 2015"),
	downloadButton('downloadData4', 'Download'),
                        h4("Linear regression method in html."),
                        p(strong(a(href="methods.html", "Uploaded: Oct 23, 2015."))),
	h5("Table(20151023).", strong(a(href="Table(20151023).xlsx", "Gene lists by empirical bayes and linear regression methods."))),
	h4("Joint lab meeting on Oct 23, 2015"),
	downloadButton('downloadData5', 'Download')
                      ))
            
))
