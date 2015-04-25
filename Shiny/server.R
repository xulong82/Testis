library(shiny)
library(ggplot2)
library(Gviz)
library(Rsamtools)

load("./geList.rdt")
load("./summary.rdt")
load("./alignList.rdt")

uid <- gsub("[123]", "", colnames(geList$raw))
group <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
geno <- gsub("^(W|M).*", "\\1", uid)

seqs <- alignList$seqs
tx_map <- alignList$tx_map
width <- seqs@sequence@ranges@width
names(width) <- seqnames(seqs)

options(ucscChromosomeNames=FALSE) 

shinyServer(function(input, output) {
  
  output$gviz <- renderPlot({
    tx1 <- tx_map$TX[tx_map$ID == input$transcript]
    axis <- GenomeAxisTrack()
    read1 <- alignList$sample[[input$sample]]
    from1 <- 1
    to1 <- width[tx1]
    if(input$from > 0) from1 <- input$from
    if(input$to > input$from & input$to < to1) to1 <- input$to
    plotTracks(list(seqs, axis, read1), chromosome = tx1, from = from1, to = to1, add53 = T)
  })
  
  output$boxplot <- renderPlot({
    
    geneId <- input$gene
    ge <- geList[[input$norm]]
    
    expr <- c(as.matrix(ge[geneId, ]))
    
    graph.dt <- data.frame(value = expr, 
      group = factor(uid, levels = group),
      geno = factor(geno, levels = c("W", "M")))
    rownames(graph.dt) <- NULL
    colnames(graph.dt) <- c("value", "group", "geno")
    
    ggplot(graph.dt, aes(x = group, y = value, fill = geno)) + 
    geom_boxplot() + theme_bw() +
    xlab("") + ylab("") + ggtitle(input$gene) +
    scale_fill_manual(values = c("white", "firebrick1")) +
    theme(plot.title = element_text(size=20, face="bold", vjust=2))
  })

  output$summary <- renderTable({
    geneId <- input$gene
    table[table$query == geneId, ]
  }, include.rownames = FALSE)

 output$downloadData <- downloadHandler(
   filename = "data.xlsx", content = function (file) file.copy("./data.xlsx", file)
 )
  
})
