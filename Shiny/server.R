library(shiny)
library(ggplot2)

load("./norm.rdt")
load("./summary.rdt")

uid <- gsub("[123]", "", colnames(genes))
cond <- c("WIN", "MIN", "WNONP", "MNONP", "WPLM", "MPLM") 
geno <- gsub("^(W|M).*", "\\1", uid)

shinyServer(function(input, output) {
  
  output$boxplot <- renderPlot({
    norm <- input$norm
    geneId <- input$gene
    
    if (norm == 1) load("./raw.rdt")

    if (geneId %in% rownames(genes))
      expr <- c(as.matrix(genes[geneId, ]))
    else
      expr <- rep(0, length(uid))
    
    graph.dt <- data.frame(value = expr, 
      cond = factor(uid, levels = cond),
      geno = factor(geno, levels = c("W", "M")))
    rownames(graph.dt) <- NULL
    colnames(graph.dt) <- c("value", "cond", "geno")
    
    ggplot(graph.dt, aes(x = cond, y = value, fill = geno)) + 
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
