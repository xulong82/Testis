library(xlsx)

rm(list = ls())
setwd("~/Dropbox/GitHub/Testis/")

# LOAD DATA
load("data/newList.rdt")
# for(obj in names(newList)) assign(obj, newList[[obj]])
for(obj in names(newList$more)) assign(obj, newList$more[[obj]])

cluster <- read.xlsx("Cate2015/nature14267-s1.xls", stringsAsFactors = F, sheetIndex = 1)

sub_g <- filter(cluster, Subunit == "g")
geneId <- sub_g$Gene.Name %>% unique

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

interMs <- getLDS("hgnc_symbol", "hgnc_symbol", geneId, human, attributesL = "external_gene_name", martL = mouse)
(interMs <- interMs[, 2])

interMs <- getLDS("hgnc_symbol", "hgnc_symbol", interHg$V1, human, attributesL = "external_gene_name", martL = mouse)

(interPoly <- lapply(poly, function(x) intersect(interMs, names(x))))
(interNonpoly <- lapply(nonpoly, function(x) intersect(interMs, names(x))))
